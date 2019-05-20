
#' Annotate narrowPeak using TxDB
#'
#' This function annotate the MACS2 called peaks with appropriate target transcript and
#' gene from TxDB object. Please refer to the \strong{Details} section for some
#' suggestions while doing peak annotation. \cr
#' Peaks are annnotated with following \strong{broad categories} and \emph{specific
#' types} (listed in decreasing order of preference):
#' \enumerate{
#' \item \strong{featureInPeak:} \emph{"include_tx", "include_CDS"}
#' \item \strong{nearStart:} \emph{"5UTR", "CDS_start", "tx_start"}
#' \item \strong{nearEnd:} \emph{"3UTR", "tx_end", "CDS_end"}
#' \item \strong{peakInFeature:} \emph{"inside_tx", "inside_CDS"}
#' \item \strong{upstreamTss:} \emph{"promoter", "upstream"}
#' }
#' Additionally, a \emph{pseudo} prefix is added to the peakType where a peak is
#' annotated to two target genes/features and one of it is more optimum than other.
#' The less optimum target type is prefixed with \emph{pseudo}. Please refer to the
#' \strong{Details} section for specific information on this.
#'
#' Some important observations to do before annotating ChIPseq data:
#' \enumerate{
#' \item Whether the signal is like TF/polII i.e. factor binds across whole gene or not.
#' Also see if binding is throughout the genome like CTCF factor.
#' See \code{bindingInGene, promoterLength} arguments for the details.
#' \item For the genes which are within peak region, what is the gene size (are
#' genes shorter in length than normal) and how far is the next downstream gene.
#' See \code{includeFractionCut} argument for the details.
#' \item Are there any TES or 3' UTR peaks and how confident are they?
#' \item Check the TXTYPE in TxDB object and see which type of features are of
#' interest to you. Usually tRNA, rRNA are not needed. See \code{excludeType}
#' argument for the details.
#' }
#' These observations will help to decide appropriate parameters while annotating
#' the peaks using TxDB object.
#'
#'
#' @param peakFile A narroPeak or broadPeak file. If a broadPeak file, peak center is
#' used as summit as broadPeak file does not report summit
#' @param fileFormat Format of the peak file. One of "narrowPeak" (Default) or "broadPeak".
#' @param txdb TxDB object which will be used for annotation
#' @param includeFractionCut Number between [0, 1]. If a peak covers more than this
#' fraction of feature/gene, it will be marked as include_tx/include_CDS. Default: 0.7
#' @param bindingInGene Logical: whether the ChIPseq TF binds in gene body. This is
#' useful for polII ChIPseq data. Default: FALSE
#' @param promoterLength Promoter length in number of nucleotides. Default: 500
#' @param insideSkewToEndCut A floating point number in range [0, 1]. If a peak is
#' present inside feature/gene and the relative summit position is > insideSkewToEndCut,
#' it is closer to the end of the feature. Default: 0.7
#' @param excludeType Types of transcripts to exclude from annotation. Should be a
#' character vector. Default: \code{c("tRNA", "rRNA", "snRNA", "snoRNA", "ncRNA")}
#' @param output Optionally store the annotation output to a file
#' @param reportPseudo Logical. Whether to report peak targets which are marked as
#' pseudo. Default: TRUE
#'
#' @return A GenomicRanges object with peak annotation
#' @export
#'
#' @examples NA
narrowPeak_annotate <- function(peakFile, fileFormat = "narrowPeak", txdb, includeFractionCut = 0.7,
                                bindingInGene = FALSE, promoterLength = 500,
                                insideSkewToEndCut = 0.7, reportPseudo = TRUE,
                                excludeType = c("tRNA", "rRNA", "snRNA", "snoRNA", "ncRNA"),
                                output = NULL){

  stopifnot(is(object = txdb, class2 = "TxDb"))

  ## started working for peak_annotation on larger genomes
  fileFormat <- match.arg(arg = fileFormat, choices = c("narrowPeak", "broadPeak"))

  ## new environment for global variables which takes time to generate
  ## create once and use multiple times
  # txdbEnv <- new.env(parent = emptyenv())

  transcriptsGr <- get_txdb_transcripts_gr(txdb = txdb, excludeType = excludeType)

  ## calculate peak related features
  peaks <- rtracklayer::import(con = peakFile, format = fileFormat)

  if(length(peaks) == 0){
    warning("no peak found in peak file")
    return(NULL)
  }

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
  fiveUtrGrl <- get_txdb_fiveUtr_grl(txdb = txdb)
  fiveUtrTargets <- UTR_annotate(queryGr = peaks, subjectGrl = fiveUtrGrl, utrType = "5UTR", txdb = txdb)

  ## 3' UTR region annotations
  threeUtrGrl <- get_txdb_threeUtr_grl(txdb = txdb)
  threeUtrTargets <- UTR_annotate(queryGr = peaks, subjectGrl = threeUtrGrl, utrType = "3UTR", txdb = txdb)

  ## exons
  exonsGr <- get_txdb_exons_gr(txdb = txdb)

  ## introns
  intronsGr <- get_txdb_introns_gr(txdb)


  # ## CDS region annotations
  # cdsGrl <- GenomicFeatures::cdsBy(x = txdb, by = "tx")
  # cdsGr <- unlist(range(cdsGrl))
  # mcols(cdsGr)$tx_id <- names(cdsGr)
  # cdsTargets <- region_overlap_annotate(queryGr = peaks,
  #                                       subjectGr = cdsGr,
  #                                       includeFractionCut = includeFractionCut,
  #                                       name = "CDS")

  ## Transcript region annotations
  transcriptTargets <- region_overlap_annotate(queryGr = peaks,
                                               subjectGr = transcriptsGr,
                                               includeFractionCut = includeFractionCut,
                                               name = "tx")

  ## annotate upstream targets: IMP to give excludeType so that rRNA, tRNA, snRNAs will be removed
  upstreamTargets <- upstream_annotate(peaksGr = peaks, featuresGr = transcriptsGr,
                                       txdb = txdb, excludeType = excludeType,
                                       promoterLength = promoterLength)


  ## prepare target preference list and peak category list
  ## this is internal preference list
  ## this order is IMP for: select 3UTR between 3UTR and inside_tx as it is more specific
  peakTypes <- data.frame(
    peakType = c("include_tx", "include_CDS", "5UTR", "CDS_start", "tx_start", "3UTR", "tx_end",
                 "CDS_end", "inside_tx", "inside_CDS", "promoter", "upstream", "pseudo_promoter", "pseudo_upstream"),
    peakPosition = c("TSS", "TSS", "TSS", "TSS", "TSS", "TES", "TES",
                     "TES", "TSS", "TSS", "TSS", "TSS", "TSS", "TSS"),
    preference = 1:14,
    stringsAsFactors = FALSE)

  ## later these peak categories will be used to decide which type of target to prefer
  peakCategories <- list(
    featureInPeak = c("include_tx", "include_CDS"),
    nearStart = c("5UTR", "CDS_start", "tx_start"),
    nearEnd = c("3UTR", "tx_end", "CDS_end"),
    peakInFeature = c("inside_tx", "inside_CDS"),
    upstreamTss = c("promoter", "upstream", "pseudo_promoter", "pseudo_upstream")
  )

  peakCategoryDf <- map_dfr(.x = peakCategories,
                            .f = function(x){data.frame(peakType = x, stringsAsFactors = F)},
                            .id = "peakCategory")

  ## combine annotations
  peakAnnotations <- NULL
  if(!is.null(fiveUtrTargets)){ peakAnnotations <- append(peakAnnotations, fiveUtrTargets) }
  if(!is.null(threeUtrTargets)){ peakAnnotations <- append(peakAnnotations, threeUtrTargets) }
  if(!is.null(cdsTargets)){ peakAnnotations <- append(peakAnnotations, cdsTargets) }
  if(!is.null(transcriptTargets)){ peakAnnotations <- append(peakAnnotations, transcriptTargets) }
  if(!is.null(upstreamTargets)){ peakAnnotations <- append(peakAnnotations, upstreamTargets) }

  if(!is.null(peakAnnotations)){

    ## remove "tRNA", "rRNA", "snRNA", "snoRNA", "ncRNA"
    allTargetsDf <- as.data.frame(peakAnnotations) %>%
      dplyr::left_join(y = txToGene, by = c("tx_id" = "TXID")) %>%
      dplyr::left_join(y = peakTypes, by = c("peakType" = "peakType")) %>%
      dplyr::left_join(y = peakCategoryDf, by = c("peakType" = "peakType")) %>%
      dplyr::filter(! txType %in% excludeType)


    ## extract best transcript region annotation for each peak-gene combination using the preference
    ## using geneId instead of tx_id: for gene which have multiple tx, 5UTR/3UTR/Intron can get
    ## annotated with any tx
    bestPeakGeneTargets <- allTargetsDf %>%
      dplyr::mutate(bidirectional = 0) %>%
      dplyr::group_by(name, geneId) %>%
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
    # select_optimal_targets(peakGr = tempTargetGrl$An_kdmB_20h_HA_1_withCtrl_peak_851)
    # endoapply(X = tempTargetGrl[1:100], FUN = select_optimal_targets)

    ## for each peak, find optimum target/s
    peakTargetGrl <- endoapply(
      X = GenomicRanges::split(x = bestPeakGeneTargetsGr, f = mcols(bestPeakGeneTargetsGr)$name),
      FUN = select_optimal_targets,
      insideSkewToEndCut = insideSkewToEndCut,
      promoterLength = promoterLength, bindingInGene = bindingInGene
    )

    peakTargetsGr <- unlist(peakTargetGrl, use.names = FALSE)

    ## add the unannotated peaks
    peakTargetsGr <- c(peakTargetsGr,
                       peaks[which(!peaks$name %in% peakTargetsGr$name)],
                       ignore.mcols=FALSE)

    ## optionally filter peak targets which are marked as pseudo
    if(!reportPseudo){
      peakTargetsGr <- peakTargetsGr[grepl(pattern = "pseudo_", x = mcols(peakTargetsGr)$peakType)]
    }

  } else{
    ## use the original peakset
    peakTargetsGr <- peaks

    mcols(peakTargetsGr)$peakType <- NA
    mcols(peakTargetsGr)$peakDist <- NA
    mcols(peakTargetsGr)$featureCovFrac <- NA
    mcols(peakTargetsGr)$summitDist <- NA
    mcols(peakTargetsGr)$geneId <- NA
    mcols(peakTargetsGr)$txName <- NA
    mcols(peakTargetsGr)$peakPosition <- NA
    mcols(peakTargetsGr)$peakCategory <- NA
    mcols(peakTargetsGr)$bidirectional <- NA

  }



  peakTargetsGr <- sort(peakTargetsGr)

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
  mcols(peakTargetsGr)$tx_id <- NULL
  mcols(peakTargetsGr)$txType <- NULL
  mcols(peakTargetsGr)$preference <- NULL
  names(peakTargetsGr) <- NULL

  ## optionally store the data
  if(!is.null(output)){
    readr::write_tsv(x = as.data.frame(mcols(peakTargetsGr)), path = output)
  }

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

  utrType <- match.arg(arg = toupper(utrType), choices = c("5UTR", "3UTR"))
  stopifnot(is(object = queryGr, class2 = "GRanges"))
  stopifnot(is(object = subjectGrl, class2 = "CompressedGRangesList"))
  stopifnot(is(object = txdb, class2 = "TxDb"))

  ## remember to combine the multi-exon UTRs from UTR GRangesList
  utrGr <- unlist(range(subjectGrl))
  mcols(utrGr)$tx_id <- names(utrGr)

  utrOvlp <- GenomicRanges::findOverlaps(query = queryGr, subject = utrGr)

  if(length(utrOvlp) == 0){
    return(NULL)
  }

  queryTargets <- queryGr[utrOvlp@from]
  mcols(queryTargets)$tx_id <- mcols(utrGr)$tx_id[utrOvlp@to]
  mcols(queryTargets)$peakType <- utrType
  mcols(queryTargets)$peakDist <- 0

  ## get respective transcripts for each UTR
  ## featureCovFrac is at transcript level
  txSubsetGrl <- GenomicFeatures::mapIdsToRanges(x = txdb,
                                                    keys = list(tx_id = mcols(queryTargets)$tx_id),
                                                    type = "tx", columns = c("gene_id"))

  txSubsetGr <- unlist(txSubsetGrl)

  mcols(queryTargets)$targetStart = start(txSubsetGr)
  mcols(queryTargets)$targetEnd = end(txSubsetGr)
  mcols(queryTargets)$targetStrand = strand(txSubsetGr)
  mcols(queryTargets)$txWidth = width(txSubsetGr)
  mcols(queryTargets)$gene_id = unlist(mcols(txSubsetGr)$gene_id)
  mcols(queryTargets)$featureCovFrac <- as.numeric(
    sprintf(fmt = "%.3f", width(pintersect(x = queryTargets, y = txSubsetGr)) / width(txSubsetGr))
  )

  ## calculate summit distance and change relativeSummitPos based on target gene
  ## summitDist > (targetEnd - targetStart) : peak overlap at end and summit is after end
  ## summitDist < 0 : peak overlap at start and summit before start
  utrTargetsDf <- as.data.frame(queryTargets) %>%
    dplyr::group_by(name, gene_id) %>%
    dplyr::arrange(desc(width), .by_group = TRUE) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::select(-gene_id, -txWidth) %>%
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
  return(sort(utrTargetsGr))
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

  name <- match.arg(arg = name, choices = c("gene", "tx", "CDS", "region"))
  stopifnot(is(object = queryGr, class2 = "GRanges"))
  stopifnot(is(object = subjectGr, class2 = "GRanges"))

  if(includeFractionCut < 0 || includeFractionCut > 1){
    stop("Invalid includeFractionCut value. Should be in range [0, 1]")
  }

  ## define feature types
  insideFeature <- paste("inside_", name, sep = "")
  includeFeature <- paste("include_", name, sep = "")
  overlapStart <- paste(name, "_start", sep = "")
  overlapEnd <- paste(name, "_end", sep = "")


  ovlpHits <- GenomicRanges::findOverlaps(query = queryGr, subject = subjectGr)

  if(length(ovlpHits) == 0){
    return(NULL)
  }

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
#' @examples NA
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
#' This function annotates the peaks with nearest downstream target. See Details.
#'
#' There will be cases when a peak is inside a gene and it is upstream of other gene
#' Use of \code{upstreamOverlappingFraction} (default: 0.2)
#'  #                                                                         #
#'	#        target1                     target2                              #
#'	#      =====<=======<===       =====<=======<========<=======             #
#'	#                                ---            ----                      #
#'	#                              peak1           peak2                      #
#'	#                      |<------>|                                         #
#'	#                      |<------------------------->|                      #
#'	#                                                                         #
#' in above cases, peak1 can be annotated as Upstream of target1. However not peak2
#' because target2 has bigger fraction in-between [target1, peak2] range
#'
#' @param peaksGr GRanges object for peak data
#' @param featuresGr A transcript GRanges object
#' @param txdb Optional TxDB object. A UTR less region is created using TxDB and used
#' to find the overlap between peak and its downstream target. A max fractional
#' overlap of 0.2 with geneA is allowed in a case when peak overlaps with a geneA and
#' is upstream of geneB. This is useful for the peaks which are near TES of a geneA.
#' If TxDB object is not provided, featuresGr is used. Default: featuresGr is used.
#' @param excludeType tx types to exclude from TxDB
#' @param upstreamOverlappingFraction Default: 0.2
#' @param promoterLength Promoter region length. Upstream peaks within \code{promoterLength}
#' distance of feature start are annotated as \code{promoter} region peaks.
#' @param ... Other arguments for \code{nearest_upstream_bidirectional()} function
#'
#' @return A modified peak GRanges object with additional columns: \code{ tx_id,
#' peakType, featureCovFrac, peakDist, summitDist}
#' @export
#'
#' @examples NA
upstream_annotate <- function(peaksGr, featuresGr, txdb = NULL, excludeType = NULL,
                              upstreamOverlappingFraction = 0.2,
                              promoterLength, ...){

  stopifnot(is(object = peaksGr, class2 = "GRanges"))
  stopifnot(is(object = featuresGr, class2 = "GRanges"))

  ## precede(ignore.strand = FALSE) is sufficient to find a peak that preceds a subject
  ## but for bidirectional peaks, precede() will return only nearest feature which is
  ## preceded by peak. However, both the features (+ and - strand) should be returned
  ## Hence, ignore.strand = TRUE is used with both precede() and follow() to extract
  ## features upstream and downstream of the peak

  ## select immediate downstream feature to the peak
  peakDownFeatures <- GenomicRanges::precede(x = peaksGr, subject = featuresGr,
                                             select = "first", ignore.strand = TRUE)

  ## select the upstreamHits which have strand == "+"
  ## because peak is downstream to the feature with strand == "-"
  peakDownHits <- data.frame(
    from = 1:length(peaksGr), peakId = mcols(peaksGr)$name,
    to = peakDownFeatures,
    featureStrand = as.vector(strand(featuresGr))[peakDownFeatures],
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
    to = peakUpFeatures,
    featureStrand = as.vector(strand(featuresGr))[peakUpFeatures],
    txName = featuresGr$tx_name[peakUpFeatures],
    stringsAsFactors = FALSE) %>%
    dplyr::filter(featureStrand == "-")

  ## this has to be on dataframe and not vectors above because the dataframes
  ## are filtering for strand
  if(nrow(peakUpHits) == 0 && nrow(peakDownHits) == 0){
    return(NULL)
  }

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
  ## Such custom regions are used because genes have very long 3' UTR.
  ## A peak in UTR region of a gene can be upstream of another gene
  txMinusUtrs <- featuresGr
  if(!is.null(txdb)){
    stopifnot(is(object = txdb, class2 = "TxDb"))

    fiveUtrGr <- unlist(range(get_txdb_fiveUtr_grl(txdb = txdb)))
    threeUtrGr <- unlist(range(get_txdb_threeUtr_grl(txdb = txdb)))
    transcriptsGr <- get_txdb_transcripts_gr(txdb = txdb, excludeType = excludeType)

    txMinusFiveUtr <- GenomicRanges::setdiff(x = transcriptsGr,
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

  ## upstreamOverlappingFraction based filtering
  ## 0.2 is still very big for the large genomes such as human, mouse as genes are very long
  upstreamHitsFiltered <- dplyr::left_join(x = upstreamHits, y = isFeatureInBetweenDf, by = "id") %>%
    tidyr::replace_na(list(intersectWd = 0, ovlpFeatureWd = 0, fractionOvlp = 0)) %>%
    dplyr::filter(fractionOvlp <= upstreamOverlappingFraction)

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
        condition = targetStrand == "-", true = 1 - relativeSummitPos, false = relativeSummitPos),
      peakType = if_else(abs(peakDist) < promoterLength, "promoter", peakType)
    )


  ## find pseudo_upstream targets
  bidirectionalPairs <- dplyr::group_by(upstreamPeaksDf, name) %>%
    dplyr::arrange(desc(peakDist), .by_group = T) %>%
    # dplyr::mutate(n = n()) %>%
    dplyr::do(nearest_upstream_bidirectional(bdirTargets = .)) %>%
    dplyr::ungroup() %>%
    # dplyr::select(-n) %>%
    dplyr::arrange(seqnames, start) %>%
    as.data.frame()

  upstreamPeaksAn <- makeGRangesFromDataFrame(bidirectionalPairs, keep.extra.columns = T)

  return(upstreamPeaksAn)

}


##################################################################################



#' True target for bidirectional peak
#'
#' This function uses following logic to select or reject target from bidirectional peak.
#' \preformatted{
#' use of minTSS_gapForPseudo (500)* and skewFraction (0.2)
#' #                      *                                                   #
#' #                     |<---more than 500bp--->|                            #
#' #       target1                    |                    target2            #
#' #   ==<=====<=====<===      peak1  |           ===>=====>=====>=====>==    #
#' #                          ------- | peak2                                 #
#' #                               ---|----                                   #
#' #                                  |                                       #
#' #                      midpoint between two targets                        #
#' #                                                                          #
#' peak1 => target1: more than 80% of the peak lies on target1 side
#' peak2 => target1, targe2: peak lies on the center
#' If gap between two bidirectional genes < minTSS_gapForPseudo (default: 500bp),
#' peak is assigned to both the targets
#' }
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
        ## negative strand target is true (peak1); set positive strand target to pseudo
        posTg <- set_peakTarget_to_pseudo(target = posTg)
      } else if((posTg$peakStart + peakFraction) >= gapCenter){
        ## positive strand target is true (peak2); set negative strand target to pseudo
        negTg <- set_peakTarget_to_pseudo(target = negTg)
      }
    } else{
      ## if the distance between the TSS of two bidirectional targets < minTSS_gapForPseudo:
      ## CANNOT decide the pseudo target confidently
      if(posTg$peakEnd <= gapCenter){
        ## negative strand target is true (peak1); set positive strand target to pseudo
        posTg <- set_peakTarget_to_pseudo(target = posTg)
      } else if(posTg$peakStart >= gapCenter){
        ## positive strand target is true (peak2); set negative strand target to pseudo
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
#' \strong{Use of \code{insideSkewToEndCut} (0.7)* and \code{promoterLength} (500)**: }
#' \preformatted{
#' #                                                                #
#' #           target1                *            target2          #
#' #   0    0.25     0.5     0.75     |<--500-->|                   #
#' #   =======>=======>=======>=======         =====>=====>===      #
#' #                            ----                                #
#' #                            peak1                               #
#' #   |<--------0.7------->|                                       #
#' #   **                                                           #
#' #                                                                #
#' In above example, peak1 is inside target1 but it is near the end
#' Relative position of the peak is >0.7 in target1. peak1 is also
#' upstream of target2 and within 500bp. So new annotation is
#' target1: pasudo_inside
#' target2: upstream
#' }
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
select_optimal_targets <- function(peakGr, insideSkewToEndCut = 0.7,
                                   promoterLength,
                                   bindingInGene = FALSE){

  # stopifnot(is(object = peakGr, class2 = "GRanges"))

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

      if(any(abs(mcols(typesGrl$upstreamTss)$peakDist) < promoterLength)){
        mcols(typesGrl$upstreamTss)$bidirectional <- 2
        typesGrl$upstreamTss <- typesGrl$upstreamTss[which(
          abs(mcols(typesGrl$upstreamTss)$peakDist) < promoterLength)]
      } else{
        typesGrl$upstreamTss <- NULL
        peakFound$upstreamTss <- FALSE
      }
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
    ## 5) for the TF which has known binding over gene body (E.g. polII ChIP or histone marks)
    ## preference is for featureInPeak. all other targets are pseudo
    if(peakFound$featureInPeak){
      typesGrl$upstreamTss <- NULL
      peakFound$upstreamTss <- FALSE
    }
  } else{

    ## for upstreamTss peaks: make sure the keep 6) and 7) independent
    ## if they are kept under same block of if(peakFound$upstreamTss)
    ## and both peakFound$nearEnd and  peakFound$peakInFeature are TRUE,
    ## 7) will give error if 6) has already set typesGrl$upstreamTss <- NULL

    ## 6) if there is upstreamTss peak and also nearEnd peak,
    ## set the upstreamTss to NULL if it is far than promoterLength
    ## ELSE just set the nearEnd peak to pseudo
    if(peakFound$upstreamTss && peakFound$nearEnd){

      if(abs(mcols(typesGrl$upstreamTss)$peakDist[1]) > promoterLength){
        typesGrl$upstreamTss <- NULL
        peakFound$upstreamTss <- FALSE
      } else{
        typesGrl$nearEnd <- set_peakTarget_to_pseudo(target = typesGrl$nearEnd)
      }
    }


    ## 7) if there is upstreamTss peak and also peakInFeature type peak
    if(peakFound$upstreamTss && peakFound$peakInFeature){

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

  peakGr <- unlist(typesGrl, use.names = FALSE)

  return(peakGr)
}


##################################################################################
