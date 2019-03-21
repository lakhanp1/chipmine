
#' Annotate narrowPeak using TxDB
#'
#' @param peakFile A narroPeak file
#' @param txdb TxDB object which will be used for annotation
#' @param includeFractionCut Number between [0, 1]. If a peak covers more than this
#' fraction of feature/gene, it will be marked as include_tx/include_CDS. Default: 0.7
#' @param bindingInGene Logical: whether the ChIPseq TF binds in gene body. This is
#' useful for polII ChIPseq data. Default: FALSE
#' @param promoterLength Promoter length in number of nucleotides. Default: 500
#' @param insideSkewToEndCut A floating point number in range [0, 1]. If a peak is
#' present inside feature/gene and the relative summit position is > insideSkewToEndCut,
#' it is closer to the end of the feature. Default: 0.7
#' @param bidirectionCut Distance cutoff to decide a peak as bidirectional for two
#' target genes. Both the target genes should be withing these bp distance from peak.
#' Default: 500
#' @param excludeType Types of transcripts to exclude from annotation. Should be a
#' character vector. Default: \code{c("tRNA", "rRNA", "snRNA", "snoRNA", "ncRNA")}
#'
#' @return A GenomicRanges object with peak annotation
#' @export
#'
#' @examples NA
narrowPeak_annotate <- function(peakFile, txdb, includeFractionCut = 0.7,
                                bindingInGene = FALSE, promoterLength = 500,
                                insideSkewToEndCut = 0.7, bidirectionCut = 500,
                                excludeType = c("tRNA", "rRNA", "snRNA", "snoRNA", "ncRNA")){



  ## calculate peak related features
  peaks <- rtracklayer::import(con = peakFile, format = "narrowPeak")
  if(is.null(mcols(peaks)$peak)){
    mcols(peaks)$peak <- as.integer(width(peaks) / 2)
  }

  mcols(peaks)$peakChr <- as.character(seqnames(peaks))
  mcols(peaks)$peakStart <- start(peaks)
  mcols(peaks)$peakEnd <- end(peaks)
  mcols(peaks)$peakSummit <- GenomicRanges::start(peaks) + mcols(peaks)$peak

  mcols(peaks)$relativeSummitPos <- as.numeric(
    sprintf("%.3f", (mcols(peaks)$peakSummit - start(peaks)) / width(peaks))
  )



  ## 5' UTR annotation
  fiveUtrGrl <- GenomicFeatures::fiveUTRsByTranscript(txdb)
  fiveUtrTargets <- UTR_annotate(queryGr = peaks, subjectGrl = fiveUtrGrl, utrType = "5UTR", txdb = txdb)

  ## 3' UTR region annotations
  threeUtrGrl <- GenomicFeatures::threeUTRsByTranscript(txdb)
  threeUtrTargets <- UTR_annotate(queryGr = peaks, subjectGrl = threeUtrGrl, utrType = "3UTR", txdb = txdb)


  ## CDS region annotations
  cdsGrl <- GenomicFeatures::cdsBy(x = txdb, by = "tx")
  cdsGr <- unlist(range(cdsGrl))
  mcols(cdsGr)$tx_id <- names(cdsGr)
  cdsTargets <- region_overlap_annotate(queryGr = peaks,
                                        subjectGr = cdsGr,
                                        includeFractionCut = includeFractionCut,
                                        name = "CDS")

  # Transcript region annotations
  transcriptsGr <- GenomicFeatures::transcripts(txdb)
  transcriptTargets <- region_overlap_annotate(queryGr = peaks,
                                               subjectGr = transcriptsGr,
                                               includeFractionCut = includeFractionCut,
                                               name = "tx")

  ## annotate upstream targets
  upstreamTargets <- upstream_annotate(peaksGr = peaks, featuresGr = transcriptsGr, txdb = txdb)


  ## prepare target preference list and peak category list
  ## this order is IMP for: select 3UTR between 3UTR and inside_tx
  peakTypes <- data.frame(
    peakType = c("include_tx", "include_CDS", "5UTR", "CDS_start", "tx_start", "3UTR", "inside_tx", "inside_CDS", "upstream", "pseudo_upstream", "tx_end", "CDS_end"),
    peakPosition = c("TSS", "TSS", "TSS", "TSS", "TSS", "TES", "TSS", "TSS", "TSS", "TSS", "TES", "TES"),
    preference = 1:12,
    stringsAsFactors = FALSE)

  peakCategories <- list(
    upstreamTss = c("upstream", "pseudo_upstream"),
    nearStart = c("5UTR", "CDS_start", "tx_start"),
    nearEnd = c("3UTR", "tx_end", "CDS_end"),
    peakInFeature = c("inside_tx", "inside_CDS"),
    featureInPeak = c("include_tx", "include_CDS")
  )

  peakCategoryDf <- map_dfr(.x = peakCategories,
                            .f = function(x){data.frame(peakType = x, stringsAsFactors = F)},
                            .id = "peakCategory")

  txToGene <- suppressMessages(
    AnnotationDbi::select(x = txdb, keys = AnnotationDbi::keys(x = txdb, keytype = "TXID"),
                                    columns = c("GENEID", "TXNAME", "TXTYPE"), keytype = "TXID")) %>%
    dplyr::mutate(TXID = as.character(TXID)) %>%
    dplyr::rename(geneId = GENEID, txName = TXNAME, txType = TXTYPE)

  ## combine annotations
  peakAnnotations <- c(fiveUtrTargets, threeUtrTargets, cdsTargets, transcriptTargets, upstreamTargets)

  ## remove "tRNA", "rRNA", "snRNA", "snoRNA", "ncRNA"
  allTargetsDf <- as.data.frame(peakAnnotations) %>%
    dplyr::left_join(y = txToGene, by = c("tx_id" = "TXID")) %>%
    dplyr::left_join(y = peakTypes, by = c("peakType" = "peakType")) %>%
    dplyr::left_join(y = peakCategoryDf, by = c("peakType" = "peakType")) %>%
    dplyr::filter(! txType %in% excludeType)


  ## extract one best TSS and TES region peak for each peak-gene combination using the preference
  bestPeakGeneTargets <- allTargetsDf %>%
    dplyr::mutate(bidirectional = 0) %>%
    dplyr::group_by(name, tx_id) %>%
    dplyr::arrange(preference, .by_group = TRUE) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup()


  bestPeakGeneTargetsGr <- makeGRangesFromDataFrame(df = bestPeakGeneTargets, keep.extra.columns = T)

  # ## for summary and debugging
  # tmpDf <- dplyr::group_by(bestPeakGeneTargets, name) %>%
  #   dplyr::summarise(n = n(),
  #                    peakType = paste(peakType, collapse = ","),
  #                    peakCat = paste(peakCategory, collapse = ","),
  #                    gene = paste(GENEID, collapse = ",")) %>%
  #   dplyr::filter(n > 1) %>%
  #   dplyr::distinct(peakType, .keep_all = T) %>%
  #   as.data.frame()

  # ## for testing select_optimal_targets()
  # tempTargetGrl <- GenomicRanges::split(x = bestPeakGeneTargetsGr, f = mcols(bestPeakGeneTargetsGr)$name)
  # select_optimal_targets(peakGr = tempTargetGrl[[1]])

  ## for each peak, find optimum target/s
  peakTargetGrl <- endoapply(
    X = GenomicRanges::split(x = bestPeakGeneTargetsGr, f = mcols(bestPeakGeneTargetsGr)$name),
    FUN = select_optimal_targets,
    insideSkewToEndCut = insideSkewToEndCut,
    promoterLength = promoterLength)

  peakTargetsGr <- unlist(peakTargetGrl, use.names = FALSE)

  ## rename columns: "peakId", "peakEnrichment", "peakPval", "peakQval"
  mcols(peakTargetsGr)$peakId <- mcols(peakTargetsGr)$name
  mcols(peakTargetsGr)$peakEnrichment <- mcols(peakTargetsGr)$signalValue
  mcols(peakTargetsGr)$peakPval <- mcols(peakTargetsGr)$pValue
  mcols(peakTargetsGr)$peakQval <- mcols(peakTargetsGr)$qValue

  ## remove unnecessary columns
  mcols(peakTargetsGr)$name <- NULL
  mcols(peakTargetsGr)$signalValue <- NULL
  mcols(peakTargetsGr)$pValue <- NULL
  mcols(peakTargetsGr)$qValue <- NULL
  mcols(peakTargetsGr)$peak <- NULL
  mcols(peakTargetsGr)$score <- NULL
  mcols(peakTargetsGr)$targetStart <- NULL
  mcols(peakTargetsGr)$targetEnd <- NULL
  mcols(peakTargetsGr)$targetStrand <- NULL
  mcols(peakTargetsGr)$peakCategory <- NULL
  mcols(peakTargetsGr)$tx_id <- NULL
  mcols(peakTargetsGr)$txType <- NULL
  names(peakTargetsGr) <- NULL

  return(peakTargetsGr)
}




##################################################################################


#' Annotate peak targets in UTR regions
#'
#' @param queryGr GRanges object generated from narrowPeak or broadPeak file
#' @param subjectGrl GRangesList object for CDS generated by \code{GenomicFeatures::fiveUTRsByTranscript()}
#' or \code{GenomicFeatures::threeUTRsByTranscript()}
#' @param utrType A character string either of \code{"5UTR"} or \code{"3UTR"}. This string is added
#' as peakType column value for the peaks which overlap the respective UTR region.
#' @param txdb TxDB object
#'
#' @return A modified peak GRanges object with additional columns: \code{ tx_id,
#' peakType, featureCovFrac, peakDist, summitDist}
#' @export
#'
#' @examples NA
UTR_annotate <- function(queryGr, subjectGrl, utrType, txdb){

  if(class(subjectGrl) != "CompressedGRangesList"){
    stop("subjectGrl should be a GRangesList object")
  }

  if(! any(utrType %in% c("5UTR", "3UTR"))){
    stop("Wrong argument utrType: should be one of \"5UTR\" or \"3UTR\")")
  }

  ## remember to combine the multi-exon UTRs from UTR GRangesList
  utrGr <- unlist(range(subjectGrl))
  mcols(utrGr)$tx_id <- names(utrGr)

  utrOvlp <- GenomicRanges::findOverlaps(query = queryGr, subject = utrGr)

  queryTargets <- queryGr[utrOvlp@from]
  mcols(queryTargets)$tx_id <- mcols(utrGr)$tx_id[utrOvlp@to]
  mcols(queryTargets)$peakType <- utrType
  mcols(queryTargets)$peakDist <- 0

  ## get respective transcripts for each UTR
  ## featureCovFrac is at transcript level
  transcriptsGrl <- GenomicFeatures::mapIdsToRanges(x = txdb,
                                                    keys = list(tx_id = mcols(queryTargets)$tx_id),
                                                    type = "tx")
  transcriptGr <- unlist(transcriptsGrl)



  mcols(queryTargets)$targetStart = start(transcriptGr)
  mcols(queryTargets)$targetEnd = end(transcriptGr)
  mcols(queryTargets)$targetStrand = strand(transcriptGr)
  mcols(queryTargets)$featureCovFrac <- as.numeric(
    sprintf(fmt = "%.3f", width(pintersect(x = queryTargets, y = transcriptGr)) / width(transcriptGr))
  )

  ## calculate summit distance and change relativeSummitPos based on target gene
  ## summitDist > (targetEnd - targetStart) : peak overlap at end and summit is after end
  ## summitDist < 0 : peak overlap at start and summit before start
  utrTargetsDf <- as.data.frame(queryTargets) %>%
    dplyr::mutate(
      summitDist = dplyr::case_when(
        targetStrand == "+" ~ peakSummit - targetStart,
        targetStrand == "-" ~ targetEnd - peakSummit
      )
    ) %>%
    dplyr::mutate(
      relativeSummitPos = dplyr::case_when(
        summitDist > (targetEnd - targetStart) ~ 1,
        summitDist < 0 ~ 0,
        summitDist > 0 ~ as.numeric(sprintf("%.3f", (peakSummit - targetStart) / (targetEnd - targetStart))),
        TRUE ~ relativeSummitPos
      )
    ) %>%
    dplyr::mutate(
      relativeSummitPos = dplyr::if_else(
        condition = targetStrand == "-", true = 1 - relativeSummitPos, false = relativeSummitPos)
    )

  utrTargetsGr <- makeGRangesFromDataFrame(df = utrTargetsDf, keep.extra.columns = TRUE)
  return(utrTargetsGr)
}


##################################################################################


#' Map peaks to given GRanges regions
#'
#' This function annotates the peaks onto a regions into e.g. \code{CDS_start, CDS_end, include_CDS,
#' inside_CDS} categories.
#' In addition, relativeSummitPos value is updated w.r.t. region for the peaks which are
#' annotated as e.g. \code{inside_CDS}
#'
#' @param queryGr GRanges object generated from narrowPeak or broadPeak file
#' @param subjectGr GRanges object for regions on which peaks needs to be mapped. E.g:
#'  \code{GenomicFeatures::cdsBy()} or \code{GenomicFeatures::transcripts()}
#' @param includeFractionCut A cutoff on peak coverage value of CDS. If peak covers more
#' than this proportion of CDS, it is annotated as \code{inside_CDS}. Default: 0.7
#' @param name Feature type to be used as suffix in peak type annotation. Eg. CDS, gene etc.
#'
#' @return A modified peak GRanges object with additional columns: \code{ tx_id,
#' peakType, featureCovFrac, peakDist, summitDist}
#' @export
#'
#' @examples NA
region_overlap_annotate <- function(queryGr, subjectGr, includeFractionCut = 0.7, name = "CDS"){

  if(includeFractionCut < 0 || includeFractionCut > 1){
    stop("Invalid includeFractionCut value. Should be in range [0, 1]")
  }

  ## define feature types
  insideFeature <- paste("inside_", name, sep = "")
  includeFeature <- paste("include_", name, sep = "")
  overlapStart <- paste(name, "_start", sep = "")
  overlapEnd <- paste(name, "_end", sep = "")


  ovlpHits <- GenomicRanges::findOverlaps(query = queryGr, subject = subjectGr)

  queryTargets <- queryGr[ovlpHits@from]
  mcols(queryTargets)$tx_id <- mcols(subjectGr)$tx_id[ovlpHits@to]
  mcols(queryTargets)$peakType <- insideFeature
  mcols(queryTargets)$peakDist <- 0

  ##
  mcols(queryTargets)$targetStart = start(subjectGr[ovlpHits@to])
  mcols(queryTargets)$targetEnd = end(subjectGr[ovlpHits@to])
  mcols(queryTargets)$targetStrand = strand(subjectGr[ovlpHits@to])
  mcols(queryTargets)$featureCovFrac <- as.numeric(
    sprintf(fmt = "%.3f", width(pintersect(x = queryTargets, y = subjectGr[ovlpHits@to])) / width(subjectGr[ovlpHits@to])))

  ## assign appropriate target location
  ## if the peak is inside_CDS, update the relativeSummitPos w.r.t. CDS
  ##
  ##     |>=====>=====>======>======>======>|            |<=====<=====<======<======<======<|
  ##   -------                                                                           --------
  ##                                    -------       -------
  ##                 -------                                          --------
  ##  -------------------------------------------     -------------------------------------------
  ##
  ## calculate summit distance and change relativeSummitPos based on target gene
  targetsDf <- as.data.frame(queryTargets) %>%
    dplyr::mutate(
      peakType = dplyr::case_when(
        peakStart <= targetStart & peakEnd >= targetEnd ~ includeFeature,
        featureCovFrac >= includeFractionCut ~ includeFeature,
        targetStrand == "+" & peakStart <= targetStart & peakEnd < targetEnd ~ overlapStart,
        targetStrand == "+" & targetStart < peakStart & targetEnd <= peakEnd ~ overlapEnd,
        targetStrand == "-" & peakStart <= targetStart & peakEnd < targetEnd ~ overlapEnd,
        targetStrand == "-" & targetStart < peakStart & targetEnd <= peakEnd ~ overlapStart,
        TRUE ~ insideFeature
      )
    ) %>%
    dplyr::mutate(
      summitDist = dplyr::case_when(
        peakType != overlapEnd & targetStrand == "+" ~ peakSummit - targetStart,
        peakType != overlapEnd & targetStrand == "-" ~ targetEnd - peakSummit,
        peakType == overlapEnd & targetStrand == "+" ~ peakSummit - targetEnd,
        peakType == overlapEnd & targetStrand == "-" ~ targetStart - peakSummit
      ),
      relativeSummitPos = dplyr::if_else(
        condition = peakType == insideFeature,
        true = as.numeric(sprintf("%.3f", (peakSummit - targetStart) / (targetEnd - targetStart))),
        false = relativeSummitPos
      )
    ) %>%
    dplyr::mutate(
      relativeSummitPos = dplyr::if_else(
        condition = targetStrand == "-", true = 1 - relativeSummitPos, false = relativeSummitPos)
    )

  ## convert back to GRanges
  queryTargets <- makeGRangesFromDataFrame(df = targetsDf, keep.extra.columns = TRUE)

  return(queryTargets)
}

##################################################################################

#' Set peakType to pseudo
#'
#' @param target A dataframe or GRanges object which has peakType column
#'
#' @return Same object with \code{pseudo} prefix to peakType column values
#' @export
#'
#' @examples
set_peakTarget_to_pseudo <- function(target){
  if(any(class(target) %in% "GRanges")){
    mcols(target)$peakType <- paste("pseudo_", mcols(target)$peakType, sep = "")
  } else if(any(class(target) %in% "data.frame")){
    target$peakType <- paste("pseudo_", target$peakType, sep = "")
  }

  return(target)
}

##################################################################################


#' Annotate upstream peaks on transcripts
#'
#' @param peaksGr GRanges object for peak data
#' @param featuresGr A transcript GRanges object
#' @param txdb Optional TxDB object. A UTR less region is created using TxDB and used
#' to find the overlap between peak and its downstream target. A max fractional
#' overlap of 0.2 with geneA is allowed in a case when peak overlaps with a geneA and
#' is upstream of geneB. This is useful for the peaks which are near TES of a geneA.
#' If TxDB object is not provided, featuresGr is used. Default: featuresGr is used.
#' @param ... Other arguments for \code{nearest_upstream_bidirectional()} function
#'
#' @return A modified peak GRanges object with additional columns: \code{ tx_id,
#' peakType, featureCovFrac, peakDist, summitDist}
#' @export
#'
#' @examples NA
upstream_annotate <- function(peaksGr, featuresGr, txdb = NULL, ...){

  ## select immediate downstream feature to the peak
  peakDownFeatures <- GenomicRanges::precede(x = peaksGr, subject = featuresGr,
                                             select = "first", ignore.strand = TRUE)

  ## select the upstreamHits which have strand == "+"
  ## because peak is downstream to the feature with strand == "-"
  peakDownHits <- data.frame(
    from = 1:length(peaksGr), peakId = mcols(peaksGr)$name,
    to = peakDownFeatures, featureStrand = strand(featuresGr[peakDownFeatures]),
    txName = featuresGr$tx_name[peakDownFeatures],
    stringsAsFactors = FALSE) %>%
    dplyr::filter(featureStrand == "+")


  ## select immediate upstream feature to the peak
  peakUpFeatures <- GenomicRanges::follow(x = peaksGr, subject = featuresGr,
                                          select = "last", ignore.strand = TRUE)

  ## select the upstreamHits which have strand == "-"
  ## because peak is downstream to the feature with strand == "+"
  peakUpHits <- data.frame(
    from = 1:length(peaksGr), peakId = mcols(peaksGr)$name,
    to = peakUpFeatures, featureStrand = strand(featuresGr[peakUpFeatures]),
    txName = featuresGr$tx_name[peakUpFeatures],
    stringsAsFactors = FALSE) %>%
    dplyr::filter(featureStrand == "-")


  ## merge the putative upstream hits
  upstreamHits <- dplyr::bind_rows(peakDownHits, peakUpHits) %>%
    dplyr::arrange(from) %>%
    dplyr::mutate(id = row_number())


  ## find the number of genes between peak and its target.
  ## only those targets are true where there is no other feature inbetween
  ## build a GRanges object of the gap region between peak and target gene

  ## generate peak-target gap GRanges
  peakTargetGapsGr <- GenomicRanges::pgap(x = peaksGr[upstreamHits$from],
                                          y = featuresGr[upstreamHits$to])

  peakTargetGapsGr <- unstrand(peakTargetGapsGr)
  names(peakTargetGapsGr) <- upstreamHits$id

  ## build a subject GRanges for tx - (5UTR + 3UTR)
  ## Such custom regions are used because genes have 3' UTR. A peak in UTR region of
  ## a gene can be upstream of another gene
  txMinusUtrs <- featuresGr
  if(!is.null(txdb)){
    fiveUtrGr <- unlist(range(GenomicFeatures::fiveUTRsByTranscript(txdb)))
    threeUtrGr <- unlist(range(GenomicFeatures::threeUTRsByTranscript(txdb)))
    txMinusFiveUtr <- GenomicRanges::setdiff(x = GenomicFeatures::transcripts(txdb),
                                             y = fiveUtrGr,
                                             ignore.strand= TRUE)

    txMinusUtrs <- GenomicRanges::setdiff(x = txMinusFiveUtr,
                                          y = threeUtrGr,
                                          ignore.strand= TRUE)
  }

  ## find first overlapping gene/feature in gap GRanges
  ## no need to find all.
  featureInGap <- GenomicRanges::findOverlaps(query = peakTargetGapsGr,
                                              subject = txMinusUtrs,
                                              select = "first",
                                              ignore.strand = TRUE)

  isFeatureInBetweenDf <- data.frame(
    id = as.numeric(names(peakTargetGapsGr)),
    gapGrRow = 1:length(featureInGap),
    firstOverlapFeature = featureInGap,
    # fractionOvlp = 1,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(firstOverlapFeature))

  ## calculate fraction overlap: needed for the peak which is inside a gene
  isFeatureInBetweenDf$intersectWd <- width(
    GenomicRanges::pintersect(
      x = peakTargetGapsGr[isFeatureInBetweenDf$gapGrRow],
      y = txMinusUtrs[isFeatureInBetweenDf$firstOverlapFeature],
      ignore.strand = TRUE)
  )

  isFeatureInBetweenDf$ovlpFeatureWd <- width(txMinusUtrs[isFeatureInBetweenDf$firstOverlapFeature])
  isFeatureInBetweenDf$fractionOvlp <- isFeatureInBetweenDf$intersectWd / isFeatureInBetweenDf$ovlpFeatureWd

  upstreamHitsFiltered <- dplyr::left_join(x = upstreamHits, y = isFeatureInBetweenDf, by = "id") %>%
    tidyr::replace_na(list(intersectWd = 0, ovlpFeatureWd = 0, fractionOvlp = 0)) %>%
    dplyr::filter(fractionOvlp <= 0.2)


  ## build upstream peaks data and filter unnecessary peaks where there is/are genes between peak and target
  upstreamPeaks <- peaksGr[upstreamHitsFiltered$from]

  mcols(upstreamPeaks)$tx_id <- mcols(featuresGr)$tx_id[upstreamHitsFiltered$to]
  mcols(upstreamPeaks)$peakType <- "upstream"
  mcols(upstreamPeaks)$peakDist <- GenomicRanges::distance(x = upstreamPeaks,
                                                           y = featuresGr[upstreamHitsFiltered$to])
  mcols(upstreamPeaks)$peakDist <- mcols(upstreamPeaks)$peakDist * -1

  mcols(upstreamPeaks)$targetStart = start(featuresGr[upstreamHitsFiltered$to])
  mcols(upstreamPeaks)$targetEnd = end(featuresGr[upstreamHitsFiltered$to])
  mcols(upstreamPeaks)$targetStrand = strand(featuresGr[upstreamHitsFiltered$to])

  ## calculate summit distance and change relativeSummitPos based on target gene
  upstreamPeaksDf <- as.data.frame(upstreamPeaks, stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      featureCovFrac = 0,
      summitDist = dplyr::case_when(
        targetStrand == "+" ~ peakSummit - targetStart,
        targetStrand == "-" ~ targetEnd - peakSummit
      )
    ) %>%
    dplyr::mutate(
      relativeSummitPos = dplyr::if_else(
        condition = targetStrand == "-", true = 1 - relativeSummitPos, false = relativeSummitPos)
    )


  ## find pseudo_upstream targets
  bidirectionalPairs <- dplyr::group_by(upstreamPeaksDf, name) %>%
    dplyr::arrange(desc(peakDist), .by_group = T) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::do(nearest_upstream_bidirectional(bdirTargets = .)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(seqnames, start) %>%
    as.data.frame() %>%
    dplyr::select(-n)

  upstreamPeaksAn <- makeGRangesFromDataFrame(bidirectionalPairs, keep.extra.columns = T)

  return(upstreamPeaksAn)

}


##################################################################################



#' True target for bidirectional peak
#'
#' This function
#'
#'     target1              peak                         target2
#'
#' ==<=====<=====<===     -------             ===>=====>=====>=====>==
#'
#'                                |
#'                    center between two targets
#'
#' True target: target1: more than 80% of the peak lies on target1 side
#'
#' @param bdirTargets A dataframe with two rows for bidirectional targets
#' @param skewFraction Minimum fraction of peak region allowed on the side of false target
#' from the midpoint of two target genes. Default: 0.2
#' @param minTSS_gapForPseudo Valid distance between two target genes TSS to mark one as psuedo
#'
#' @return Same dataframe in which peakType for one of the target is marked as pseudo
#' @export
#'
#' @examples NULL
nearest_upstream_bidirectional <- function(bdirTargets, skewFraction = 0.2, minTSS_gapForPseudo = 500){


  ##     target1              peak                         target2
  ## ==<=====<=====<===     -------             ===>=====>=====>=====>==
  ##                                |
  ##                    center between two targets
  ## True target: target1: more than 80% of the peak lies on target1 side

  if(nrow(bdirTargets) == 1){
    return(bdirTargets)
  } else if(nrow(bdirTargets) > 2){
    stop("Unhandled behaviour: 3 targets for bidirectional peaks")
  }

  if(all(bdirTargets$targetStrand %in% c("+", "-")) && bdirTargets$targetStrand[1] != bdirTargets$targetStrand[2]){

    posTg <- bdirTargets[ which(bdirTargets$targetStrand == "+"), ]
    negTg <- bdirTargets[ which(bdirTargets$targetStrand == "-"), ]

    ## gap center
    gapStart <- negTg$targetEnd
    gapEnd <- posTg$targetStart
    gapWidth <- gapEnd - gapStart
    gapCenter <- negTg$targetEnd + (gapWidth / 2)

    peakWidth <- posTg$peakEnd - posTg$peakStart
    peakFraction <- peakWidth * skewFraction


    if(gapWidth > minTSS_gapForPseudo){

      if((posTg$peakEnd - peakFraction) <= gapCenter){
        ## negative strand target is true; set positive strand target to pseudo
        posTg <- set_peakTarget_to_pseudo(target = posTg)
      } else if((posTg$peakStart + peakFraction) >= gapCenter){
        ## positive strand target is true; set negative strand target to pseudo
        negTg <- set_peakTarget_to_pseudo(target = negTg)
      }
    } else{
      ## if the distance between the TSS of two bidirectional targets < minTSS_gapForPseudo:
      ## CANNOT decide the pseudo target confidently
      if(posTg$peakEnd <= gapCenter){
        ## negative strand target is true; set positive strand target to pseudo
        posTg <- set_peakTarget_to_pseudo(target = posTg)
      } else if(posTg$peakStart >= gapCenter){
        ## positive strand target is true; set negative strand target to pseudo
        negTg <- set_peakTarget_to_pseudo(target = negTg)
      }
    }

    return(dplyr::bind_rows(posTg, negTg))

  }

  return(bdirTargets)
}


##################################################################################


#' Assign best target/s for each peak
#'
#' This function checks the different targets assigned for a peak and returns the
#' optimum target gene/s
#'
#' @param peakGr peak annotation in form of GRanges object. This is generated by
#' combining \code{UTR_annotate(), region_overlap_annotate(), upstream_annotate()}
#' functions.
#' @param insideSkewToEndCut A floating point number in range [0, 1]. If a peak is
#' present inside feature/gene and the relative summit position is > insideSkewToEndCut,
#' it is closer to the end of the feature. Default: 0.7
#' @param promoterLength Promoter length in number of nucleotides
#' @param bindingInGene Logical: whether the ChIPseq TF binds in gene body. This is
#' useful for polII ChIPseq data. Default: FALSE
#'
#' @return Same GRanges object as peakGr but with modified peakType column or excluding
#' some targets which were not optimum.
#' @export
#'
#' @examples NA
select_optimal_targets <- function(peakGr, insideSkewToEndCut = 0.7, promoterLength = 500,
                                bindingInGene = FALSE){
  ## if only one target, return as it is
  if(length(peakGr) == 1){
    return(peakGr)
  }

  ## apply various filters to select best target/s

  ## make GRangesList based on peakCategory
  typesGrl <- GenomicRanges::split(x = peakGr, f = mcols(peakGr)$peakCategory)

  peakFound <- list(
    upstreamTss = FALSE,
    nearStart = FALSE,
    nearEnd = FALSE,
    peakInFeature = FALSE,
    featureInPeak = FALSE
  )

  peakFound[names(typesGrl)] <- TRUE

  ## for nearStart peak
  if(peakFound$nearStart){

    ## 1) one target is nearStart and second target is upstreamTss
    if(peakFound$upstreamTss){
      ## decide pseudo_upstream based on nearest_upstream_bidirectional()
      bidirectDicision <- nearest_upstream_bidirectional(
        bdirTargets = as.data.frame(c(typesGrl$nearStart[1], typesGrl$upstreamTss[1])),
        minTSS_gapForPseudo = promoterLength
      )

      ## update upstreamTss peak. if it is marked as pseudo, it will be overwritten
      typesGrl$upstreamTss <- makeGRangesFromDataFrame(
        bidirectDicision[which(bidirectDicision$peakCategory == "upstreamTss"), ],
        keep.extra.columns = T
      )

      mcols(typesGrl$upstreamTss)$bidirectional <- 1

      ## upstreamTss peaks with peakDist > promoterLength: set to NULL
      if(abs(mcols(typesGrl$upstreamTss)$peakDist) > promoterLength){
        typesGrl$upstreamTss <- NULL
        peakFound$upstreamTss <- FALSE
      }
    }

    ## 2) one target is nearStart and second target is near TES: set TES target to NULL
    if(peakFound$nearEnd){
      typesGrl$nearEnd <- NULL
      peakFound$nearEnd <- FALSE
    }
  }


  ## for featureInPeak type peak
  if(peakFound$featureInPeak){
    ## 3) remove nearEnd target
    typesGrl$nearEnd <- NULL
    peakFound$nearEnd <- FALSE

    ## 4) upstreamTss peaks with peakDist < promoterLength: select; else reject
    if(peakFound$upstreamTss){
      mcols(typesGrl$upstreamTss)$bidirectional <- 2
      typesGrl$upstreamTss <- typesGrl$upstreamTss[which(abs(mcols(typesGrl$upstreamTss)$peakDist) < promoterLength)]
    }
  }


  ## for upstreamTss peak
  ## 8) two upstreamTss peaks: remove pseudo_upstream if peakDist > promoterLength
  if(length(typesGrl$upstreamTss) == 2){
    upstreamGrl <- GenomicRanges::split(x = typesGrl$upstreamTss, f = mcols(typesGrl$upstreamTss)$peakType)

    ## set pseudo_upstream peak to NULL if it far than promoterLength
    if(!is.null(upstreamGrl$pseudo_upstream)){
      if(abs(mcols(upstreamGrl$pseudo_upstream)$peakDist) > promoterLength){
        upstreamGrl$pseudo_upstream <- NULL
      }
    }

    ## rebuild the original upstreamTss GRanges
    typesGrl$upstreamTss <- unlist(upstreamGrl, use.names = FALSE)
  }


  ##
  if(bindingInGene){
    ## 5) for the TF which has known binding over gene body (E.g. polII ChIP)
    ## preference is for featureInPeak. all other targets are pseudo
    if(peakFound$featureInPeak){
      typesGrl$upstreamTss <- NULL
      peakFound$upstreamTss <- FALSE
    }
  } else{

    ## for upstreamTss peak
    if(peakFound$upstreamTss){

      ## 6) if there is upstreamTss peak and also nearEnd peak,
      ## set the upstreamTss to NULL if it is far than promoterLength
      ## ELSE just set the nearEnd peak to pseudo
      if(peakFound$nearEnd){

        if(abs(mcols(typesGrl$upstreamTss)$peakDist[1]) > promoterLength){
          typesGrl$upstreamTss <- NULL
          peakFound$upstreamTss <- FALSE
        } else{
          typesGrl$nearEnd <- set_peakTarget_to_pseudo(target = typesGrl$nearEnd)
        }
      }

      ## 7) if there is upstreamTss peak and also peakInFeature type peak
      if(peakFound$peakInFeature){

        ## if upstreamTss is within promoter range
        if(abs(mcols(typesGrl$upstreamTss)$peakDist[1]) <= promoterLength){
          if(mcols(typesGrl$peakInFeature)$relativeSummitPos[1] > insideSkewToEndCut){
            ## peak lies near end for peakInFeature. set peakInFeature to pseudo
            typesGrl$peakInFeature <- set_peakTarget_to_pseudo(target = typesGrl$peakInFeature)
          } else{
            ## set upstreamTss to pseudo
            typesGrl$upstreamTss <- set_peakTarget_to_pseudo(target = typesGrl$upstreamTss)
          }
        } else{
          typesGrl$upstreamTss <- NULL
          peakFound$upstreamTss <- FALSE
        }
      }
    }
  }

  peakGr <- unlist(typesGrl, use.names = FALSE)

  return(peakGr)

}


##################################################################################
