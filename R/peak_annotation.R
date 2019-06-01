
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
#' \item \strong{peakInFeature:} \emph{"exon", "intron", "inside_tx", "inside_CDS"}
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
#' @param txIds A vector of transcript IDs to be used specifically in the annotation
#' process instead of full transcript set. These should be internal tx_ids from TxDB
#' object. Default: NULL
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
narrowPeak_annotate <- function(peakFile, fileFormat = "narrowPeak",
                                txdb, txIds = NULL,
                                excludeType = c("tRNA", "rRNA", "snRNA", "snoRNA", "ncRNA"),
                                includeFractionCut = 0.7, bindingInGene = FALSE,
                                promoterLength = 500, insideSkewToEndCut = 0.7,
                                reportPseudo = TRUE,
                                output = NULL){

  stopifnot(is(object = txdb, class2 = "TxDb"))

  if(! all(txIds %in% keys(txdb, keytype = "TXID"))){
    stop("unknown TXIDs in txIds: ",
         paste(c(head(txIds[which(!txIds %in% keys(txdb, keytype = "TXID"))]), "..."),
               collapse = " ")
         )
  }

  ## started working for peak_annotation on larger genomes
  fileFormat <- match.arg(arg = fileFormat, choices = c("narrowPeak", "broadPeak"))

  ## new environment for global variables which takes time to generate
  ## create once and use multiple times
  # txdbEnv <- new.env(parent = emptyenv())

  transcriptsGr <- get_txdb_transcripts_gr(txdb = txdb, excludeType = excludeType,
                                           tx = txIds)
  txToGene <- get(x = "txToGene", envir = txdbEnv)

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
  fiveUtrGr <- get_txdb_fiveUtr_gr(txdb = txdb, tx = txIds)
  fiveUtrTargets <- splicing_unit_annotate(peaksGr = peaks, featuresGr = fiveUtrGr,
                                           featureType = "5UTR", txdb = txdb)

  ## 3' UTR region annotations
  threeUtrGr <- get_txdb_threeUtr_gr(txdb = txdb, tx = txIds)
  threeUtrTargets <- splicing_unit_annotate(peaksGr = peaks, featuresGr = threeUtrGr,
                                            featureType = "3UTR", txdb = txdb)

  ## exons annotations
  exonsGr <- get_txdb_exons_gr(txdb = txdb, tx = txIds)
  exonTargets <- splicing_unit_annotate(peaksGr = peaks, featuresGr = exonsGr,
                                        featureType = "exon", txdb = txdb)

  ## introns annotations
  intronsGr <- get_txdb_introns_gr(txdb = txdb, tx = txIds)
  intronTargets <- splicing_unit_annotate(peaksGr = peaks, featuresGr = intronsGr,
                                          featureType = "intron", txdb = txdb)

  ## Transcript region annotations
  transcriptTargets <- region_annotate(peaksGr = peaks,
                                       featuresGr = transcriptsGr,
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
    peakType = c("include_tx", "include_CDS", "5UTR", "CDS_start", "tx_start", "3UTR",
                 "tx_end", "CDS_end", "EXON", "INTRON", "inside_tx", "inside_CDS",
                 "promoter", "upstream", "pseudo_promoter", "pseudo_upstream"),
    peakPosition = c("TSS", "TSS", "TSS", "TSS", "TSS", "TES",
                     "TES", "TES", "TSS", "TSS", "TSS", "TSS",
                     "TSS", "TSS", "TSS", "TSS"),
    preference = 1:16,
    stringsAsFactors = FALSE)

  ## later these peak categories will be used to decide which type of target to prefer
  peakCategories <- list(
    featureInPeak = c("include_tx", "include_CDS"),
    nearStart = c("5UTR", "CDS_start", "tx_start"),
    nearEnd = c("3UTR", "tx_end", "CDS_end"),
    peakInFeature = c("EXON", "INTRON", "inside_tx", "inside_CDS"),
    upstreamTss = c("promoter", "upstream", "pseudo_promoter", "pseudo_upstream")
  )

  peakCategoryDf <- map_dfr(.x = peakCategories,
                            .f = function(x){data.frame(peakType = x, stringsAsFactors = F)},
                            .id = "peakCategory")

  peakTypes <- dplyr::left_join(x = peakTypes, y = peakCategoryDf, by = "peakType")

  ## combine annotations
  peakAnnotations <- NULL
  if(!is.null(fiveUtrTargets)){ peakAnnotations <- append(peakAnnotations, fiveUtrTargets) }
  if(!is.null(threeUtrTargets)){ peakAnnotations <- append(peakAnnotations, threeUtrTargets) }
  if(!is.null(exonTargets)){ peakAnnotations <- append(peakAnnotations, exonTargets) }
  if(!is.null(intronTargets)){ peakAnnotations <- append(peakAnnotations, intronTargets) }
  if(!is.null(transcriptTargets)){ peakAnnotations <- append(peakAnnotations, transcriptTargets) }
  if(!is.null(upstreamTargets)){ peakAnnotations <- append(peakAnnotations, upstreamTargets) }

  if(!is.null(peakAnnotations)){

    ## remove "tRNA", "rRNA", "snRNA", "snoRNA", "ncRNA"
    allTargetsDf <- as.data.frame(peakAnnotations, row.names = NULL) %>%
      dplyr::left_join(y = txToGene, by = c("tx_id" = "TXID")) %>%
      dplyr::left_join(y = peakTypes, by = c("peakType" = "peakType"))

    ## extract best transcript region annotation for each peak-gene combination using the preference
    ## using geneId instead of tx_id: for gene which have multiple tx, 5UTR/3UTR/Intron can get
    ## annotated with any tx
    bestPeakGeneTargets <- allTargetsDf %>%
      dplyr::mutate(bidirectional = 0) %>%
      dplyr::group_by(name, geneId) %>%
      dplyr::arrange(preference, .by_group = TRUE) %>%
      dplyr::slice(1L) %>%
      dplyr::ungroup()

    bestPeakGeneTargetsGr <- GenomicRanges::makeGRangesFromDataFrame(
      df = bestPeakGeneTargets,
      keep.extra.columns = T)

    #######################
    # ## for summary and debugging
    # tmpDf <- dplyr::group_by(bestPeakGeneTargets, name) %>%
    #   dplyr::summarise(n = n(),
    #                    peakType = paste(sort(peakType), collapse = ","),
    #                    peakCat = paste(peakCategory, collapse = ","),
    #                    uniqCat = paste(sort(unique(peakCategory)), collapse = ","),
    #                    gene = paste(geneId, collapse = ",")) %>%
    #   dplyr::filter(n > 1) %>%
    #   dplyr::distinct(peakType, .keep_all = T) %>%
    #   as.data.frame()
    #######################

    ## for testing select_optimal_targets()
    tempTargetGrl <- GenomicRanges::split(x = bestPeakGeneTargetsGr, f = mcols(bestPeakGeneTargetsGr)$name)
    # table(elementNROWS(tempTargetGrl))
    singleHitPeaks <- unlist(tempTargetGrl[which(elementNROWS(tempTargetGrl) == 1)])
    multipleHitPeaks <- tempTargetGrl[which(elementNROWS(tempTargetGrl) != 1)]

    ## for each peak, find optimum target/s
    multipleHitPeaks <- select_optimal_targets(
      targetGr = sort(unlist(multipleHitPeaks, use.names = F)),
      promoterLength = promoterLength,
      bindingInGene = bindingInGene,
      insideSkewToEndCut = insideSkewToEndCut)


    peakTargetsGr <- c(singleHitPeaks, multipleHitPeaks, ignore.mcols=FALSE)

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


#' Annotate peaks with 5UTR/exon/intron/3UTR overlap
#'
#' @param peaksGr GRanges object for peak data
#' @param featuresGr A transcript GRanges object
#' @param featureType Feature type. One of \code{c("5UTR", "exon", "intron", "3UTR")}
#' @param txdb TxDB object. Transcript information for the exon/intron overlapping with peak
#' is extracted from TxDB object and used to calculate other features. In case a peak
#' overlap with exons/introns from multiple transcripts of same gene, longest transcript is
#' selected.
#'
#' @return A modified peak GRanges object with additional columns: \code{tx_id,
#' peakType, featureCovFrac, peakDist, summitDist}
#' @export
#'
#' @examples NA
splicing_unit_annotate <- function(peaksGr, featuresGr, featureType, txdb){
  featureType <- match.arg(arg = toupper(featureType),
                           choices = c("5UTR", "EXON", "INTRON", "3UTR"))

  stopifnot(is(object = peaksGr, class2 = "GRanges"))
  stopifnot(is(object = featuresGr, class2 = "GRanges"))
  stopifnot(is(object = txdb, class2 = "TxDb"))

  ovlpHits <- GenomicRanges::findOverlaps(query = peaksGr, subject = featuresGr)

  if(length(ovlpHits) == 0){
    return(NULL)
  }

  peakTargets <- peaksGr[ovlpHits@from]
  mcols(peakTargets)$tx_id <- mcols(featuresGr)$tx_id[ovlpHits@to]
  mcols(peakTargets)$peakType <- featureType
  mcols(peakTargets)$peakDist <- 0

  ## featureCovFrac is at transcript level
  txSubsetGrl <- GenomicFeatures::mapIdsToRanges(x = txdb,
                                                 keys = list(tx_id = mcols(peakTargets)$tx_id),
                                                 type = "tx", columns = c("gene_id"))

  txSubsetGr <- unlist(txSubsetGrl)

  mcols(peakTargets)$targetStart = start(txSubsetGr)
  mcols(peakTargets)$targetEnd = end(txSubsetGr)
  mcols(peakTargets)$targetStrand = strand(txSubsetGr)
  mcols(peakTargets)$txWidth = width(txSubsetGr)
  mcols(peakTargets)$gene_id = unlist(mcols(txSubsetGr)$gene_id)
  mcols(peakTargets)$featureCovFrac <- as.numeric(
    sprintf(fmt = "%.3f", width(pintersect(x = peakTargets, y = txSubsetGr)) / width(txSubsetGr))
  )

  ## calculate summit distance and change relativeSummitPos based on target gene
  ## summitDist > (targetEnd - targetStart) : peak overlap at end and summit is after end
  ## summitDist < 0 : peak overlap at start and summit before start
  targetsDf <- as.data.frame(peakTargets) %>%
    dplyr::group_by(name, gene_id) %>%
    dplyr::arrange(desc(txWidth), .by_group = TRUE) %>%
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

  targetsGr <- sort(makeGRangesFromDataFrame(df = targetsDf, keep.extra.columns = TRUE))
  return(targetsGr)
}


##################################################################################

#' Map peaks to given GRanges regions
#'
#' This function annotates the peaks onto a regions into e.g. \code{CDS_start, CDS_end, include_CDS,
#' inside_CDS} categories.
#' In addition, relativeSummitPos value is updated w.r.t. region for the peaks which are
#' annotated as e.g. \code{inside_CDS}
#'
#' @param peaksGr GRanges object generated from narrowPeak or broadPeak file
#' @param featuresGr GRanges object for regions on which peaks needs to be mapped. E.g:
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
region_annotate <- function(peaksGr, featuresGr, includeFractionCut = 0.7, name = "CDS"){

  name <- match.arg(arg = name, choices = c("gene", "tx", "CDS", "region"))
  stopifnot(is(object = peaksGr, class2 = "GRanges"))
  stopifnot(is(object = featuresGr, class2 = "GRanges"))

  if(includeFractionCut < 0 || includeFractionCut > 1){
    stop("Invalid includeFractionCut value. Should be in range [0, 1]")
  }

  ## define feature types
  insideFeature <- paste("inside_", name, sep = "")
  includeFeature <- paste("include_", name, sep = "")
  overlapStart <- paste(name, "_start", sep = "")
  overlapEnd <- paste(name, "_end", sep = "")


  ovlpHits <- GenomicRanges::findOverlaps(query = peaksGr, subject = featuresGr)

  if(length(ovlpHits) == 0){
    return(NULL)
  }

  queryTargets <- peaksGr[ovlpHits@from]
  mcols(queryTargets)$tx_id <- mcols(featuresGr)$tx_id[ovlpHits@to]
  mcols(queryTargets)$peakType <- insideFeature
  mcols(queryTargets)$peakDist <- 0

  ##
  mcols(queryTargets)$targetStart = start(featuresGr[ovlpHits@to])
  mcols(queryTargets)$targetEnd = end(featuresGr[ovlpHits@to])
  mcols(queryTargets)$targetStrand = strand(featuresGr[ovlpHits@to])
  mcols(queryTargets)$featureCovFrac <- as.numeric(
    sprintf(fmt = "%.3f", width(pintersect(x = queryTargets, y = featuresGr[ovlpHits@to])) / width(featuresGr[ovlpHits@to])))

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
#' Target gene inside gene case:
#'  #                                                                         #
#'	#              target1          target2 (longer)                          #
#'	#      =====<=======<=======<=======<========<=======                     #
#'  #               ==<==                                                     #
#'	#                                ---            ----                      #
#'	#                              peak1           peak2                      #
#'	#                    |<-------->|                                         #
#'	#                                                                         #
#' In above case, peak1 is inside targe2 and upstream of target1. These targets are
#' selected in \code{select_optimal_targets()} if peak lies within promoter range.
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

  ## build a subject GRanges for tx - 3UTR
  ## Such custom regions are used because genes have very long 3' UTR.
  ## A peak in UTR region of a gene can be upstream of another gene
  txMinusUtrs <- featuresGr
  if(!is.null(txdb)){
    stopifnot(is(object = txdb, class2 = "TxDb"))

    fiveUtrGr <- get_txdb_fiveUtr_gr(txdb = txdb)
    threeUtrGr <- get_txdb_threeUtr_gr(txdb = txdb)
    transcriptsGr <- get_txdb_transcripts_gr(txdb = txdb, excludeType = excludeType)

    txMinusUtrs <- GenomicRanges::setdiff(x = txMinusUtrs,
                                          y = threeUtrGr,
                                          ignore.strand= TRUE)

    ## 5UTR is not needed
    # txMinusUtrs <- GenomicRanges::setdiff(x = txMinusUtrs,
    #                                       y = fiveUtrGr,
    #                                       ignore.strand= TRUE)

  }

  ## find first overlapping gene/feature in gap GRanges
  ## no need to find all.
  featureInGap <- GenomicRanges::findOverlaps(query = peakTargetGapsGr,
                                              subject = txMinusUtrs,
                                              select = "first",
                                              ignore.strand = TRUE)

  isFeatureInBetweenDf <- data.frame(
    id = as.numeric(names(peakTargetGapsGr)),
    gapWidth = width(peakTargetGapsGr),
    gapGrRow = 1:length(featureInGap),
    firstOverlapFeature = featureInGap,
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
  isFeatureInBetweenDf$ovlpFeatureFrac <- isFeatureInBetweenDf$intersectWd / isFeatureInBetweenDf$ovlpFeatureWd
  isFeatureInBetweenDf$ovlpGapFrac <- isFeatureInBetweenDf$gapWidth / isFeatureInBetweenDf$intersectWd

  ## upstreamOverlappingFraction based filtering
  ## 0.2 is still very big for the large genomes such as human, mouse as genes are very long
  upstreamHitsFiltered <- dplyr::left_join(x = upstreamHits, y = isFeatureInBetweenDf, by = "id") %>%
    tidyr::replace_na(list(intersectWd = 0, ovlpFeatureWd = 0, ovlpFeatureFrac = 0, ovlpGapFrac = 0)) %>%
    dplyr::filter(ovlpFeatureFrac <= upstreamOverlappingFraction)

  ## if no upstream peak after filtering
  if(nrow(upstreamHitsFiltered) == 0){
    return(NULL)
  }

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
  upstreamPeaks <- as.data.frame(upstreamPeaks, stringsAsFactors = FALSE) %>%
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
    ) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)


  ######
  ## select correct target/s in case of bi-directional peak
  tempTargetGrl <- GenomicRanges::split(x = upstreamPeaks, f = mcols(upstreamPeaks)$name)
  singalTargetPeaks <- tempTargetGrl[which(elementNROWS(tempTargetGrl) == 1)]
  dualTargetPeaks <- tempTargetGrl[which(elementNROWS(tempTargetGrl) > 1)]

  if(any(elementNROWS(dualTargetPeaks) > 2)){
    stop("More than two upstream targets targets found for peaks")
  }

  dualTargetPeaksDf <- sort(unlist(dualTargetPeaks, use.names = FALSE)) %>%
    as.data.frame(row.names = NULL, stringsAsFactors = FALSE) %>%
    dplyr::mutate(rowIdx = 1:n(),
                  target = TRUE)

  ## row index table for bidirectional peak
  pairTable <- dualTargetPeaksDf %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(group = 1:n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(group = paste("t", group, sep = "")) %>%
    dplyr::select(name, group, rowIdx) %>%
    tidyr::spread(key = group, value = rowIdx)

  pseudoUpIdx <- nearest_upstream_bidirectional(targetDf = dualTargetPeaksDf,
                                                t1Idx = pairTable$t1,
                                                t2Idx = pairTable$t2)

  ## remove the targets which are too far based on nearest_upstream_bidirectional()
  dualTargetFiltered <- dualTargetPeaksDf[-pseudoUpIdx,] %>%
    dplyr::select(-rowIdx, -target) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  ######


  upstreamPeaksAn <- sort(c(unlist(singalTargetPeaks, use.names = FALSE),
                            dualTargetFiltered))
  names(upstreamPeaksAn) <- NULL

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
#' @param targetDf A dataframe which has bidirectional targets for each peak.
#' @param t1Idx target1 index vector
#' @param t2Idx target2 index vector. IMP: \code{length(t1Idx)} should be equal to
#' \code{length(t2Idx)}
#' @param skewFraction Minimum fraction of peak region allowed on the side of false target
#' from the midpoint of two target genes. Default: 0.2
#' @param minTSS_gapForPseudo Valid distance between two target genes TSS to mark one as psuedo
#'
#' @return A vector of row index for pseudo targets.
#' @export
#'
#' @examples NA
nearest_upstream_bidirectional <- function(targetDf, t1Idx, t2Idx,
                                           skewFraction = 0.2, minTSS_gapForPseudo = 500){

  targetPairDf <- tibble::tibble(t1Idx = t1Idx, t2Idx = t2Idx,
                                 t1Select = TRUE, t2Select = TRUE,
                                 dir = "opposite")

  targetA <- targetDf[targetPairDf$t1Idx, ]
  targetB <- targetDf[targetPairDf$t2Idx, ]

  if(!all(targetA$name == targetB$name)){
    stop("peak name mismatch in bidirectional target pair")
  }

  targetPairDf$peakId <- targetA$name
  targetPairDf$t1PeakDist <- targetA$peakDist
  targetPairDf$t2PeakDist <- targetB$peakDist
  peakGr <- GenomicRanges::makeGRangesFromDataFrame(df = targetA)
  targetPairDf$peakWidth <- width(peakGr)
  targetPairDf$peakFraction <- round(targetPairDf$peakWidth * skewFraction)

  sameDirTargets <- which(targetA$targetStrand == targetB$targetStrand)
  targetPairDf$dir[sameDirTargets] <- "same"

  targetPairDf$t1Select[sameDirTargets] <- targetA$preference[sameDirTargets] < targetB$preference[sameDirTargets]
  targetPairDf$t2Select[sameDirTargets] <- targetA$preference[sameDirTargets] > targetB$preference[sameDirTargets]


  targetAGr <- dplyr::select(targetA, seqnames, targetStart, targetEnd, targetStrand, name, rowIdx) %>%
    makeGRangesFromDataFrame(start.field = "targetStart", end.field = "targetEnd",
                             strand.field = "targetStrand",
                             keep.extra.columns = TRUE)

  targetBGr <- dplyr::select(targetB, seqnames, targetStart, targetEnd, targetStrand, name, rowIdx) %>%
    makeGRangesFromDataFrame(start.field = "targetStart", end.field = "targetEnd",
                             strand.field = "targetStrand",
                             keep.extra.columns = TRUE)


  if(any(distance(x = targetAGr, y = targetBGr, ignore.strand = TRUE) == 0)){
    stop("Distance between bidirectional target genes cannot be 0")
  }

  ## gap between two bidirectional targets
  targetGapGr <- unstrand(pgap(x = targetAGr, y = targetBGr, ignore.strand = TRUE))
  targetPairDf$gapWidth <- width(targetGapGr)

  ## midpoint of gap between bidirectional targets
  midpoint <- resize(x = targetGapGr, width = 1, fix = "center", ignore.strand = TRUE)

  targetPairDf$midpointDist <- distance(x = midpoint, y = peakGr, ignore.strand = TRUE)

  targetPairDf <- dplyr::mutate(
    targetPairDf,
    t1Select = dplyr::case_when(
      gapWidth <= minTSS_gapForPseudo ~ t1Select,
      dir == "opposite" & (abs(t1PeakDist) + peakFraction) > (gapWidth/2) ~ FALSE,
      TRUE ~ t1Select
    ),
    t2Select = dplyr::case_when(
      gapWidth <= minTSS_gapForPseudo ~ t2Select,
      dir == "opposite" & (abs(t2PeakDist) + peakFraction) > (gapWidth/2) ~ FALSE,
      TRUE ~ t2Select
    )
  )

  falseUpIdx <- c(targetPairDf$t1Idx[which(!targetPairDf$t1Select)],
                  targetPairDf$t2Idx[which(!targetPairDf$t2Select)])

  return(falseUpIdx)
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
#'
#' @param targetGr peak annotation in form of GRanges object. This is generated by
#' combining \code{UTR_annotate(), region_annotate(), upstream_annotate()}
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
select_optimal_targets <- function(targetGr, promoterLength, bindingInGene,
                                   insideSkewToEndCut){

  ## important to sort the targets based on peakDist and relativeSummitPos
  ## rowIdx will be used later to extract correct targets
  targetDf <- as.data.frame(targetGr) %>%
    dplyr::arrange(seqnames, start, preference, peakDist, relativeSummitPos) %>%
    dplyr::mutate(rowIdx = 1:n(),
                  target = TRUE)

  subData <- targetDf %>%
    dplyr::select(seqnames, start, end, name, peakCategory, peakDist,
                  relativeSummitPos, rowIdx, target)

  masterIndexDf <- tibble::tibble(name = unique(subData$name))
  peakFound <- masterIndexDf
  peakDistDf <- masterIndexDf
  peakPosDf <- masterIndexDf

  for (ctg in c("featureInPeak", "nearStart", "nearEnd", "peakInFeature", "upstreamTss")) {
    groupedDf <- dplyr::filter(subData, peakCategory == ctg) %>%
      dplyr::select(name, rowIdx, target, peakDist, relativeSummitPos) %>%
      dplyr::group_by(name) %>%
      dplyr::arrange(rowIdx, .by_group = T)

    ## get the peak target GRanges row index for each peak category
    indexDf <- dplyr::summarise(groupedDf, !!ctg := list(rowIdx))
    masterIndexDf <- dplyr::left_join(x = masterIndexDf, y = indexDf, by = "name")

    ## create logical vector for foundPeak data for each peak category
    isTargetDf <- dplyr::summarise(groupedDf, !!ctg := any(target))
    peakFound <- dplyr::left_join(x = peakFound, y = isTargetDf, by = "name")

    ## create peak distance df for each peak category
    distDf <- dplyr::summarise(groupedDf, !!ctg := list(peakDist))
    peakDistDf <- dplyr::left_join(x = peakDistDf, y = distDf, by = "name")

    ## create peak position df for each peak category
    posDf <- dplyr::summarise(groupedDf, !!ctg := list(relativeSummitPos))
    peakPosDf <- dplyr::left_join(x = peakPosDf, y = posDf, by = "name")
  }

  markPseudoIdx <- c()

  ## nearStart peak
  ###########
  ## A1) nearStart & upstreamTss:
  ruleA1 <- which(peakFound$nearStart & peakFound$upstreamTss)
  # targetDf[unlist(masterIndexDf$nearStart[ruleA1]), ]

  ruleA1_far <- ruleA1[purrr::map_lgl(.x = peakDistDf$upstreamTss[ruleA1],
                                      .f = function(x){abs(x[1]) > promoterLength})]

  ## set the upstreamTss to NULL if it is far than promoterLength
  masterIndexDf$upstreamTss[ruleA1_far] <- purrr::map(.x = masterIndexDf$upstreamTss[ruleA1_far],
                                                      .f = ~ NULL)
  peakFound$upstreamTss[ruleA1_far] <- FALSE

  ## mark other upstreamTss within promoterLength as pseudo
  ruleA1_pro <- setdiff(x = ruleA1, ruleA1_far)
  markPseudoIdx <- append(markPseudoIdx, ruleA1_pro)


  ###########
  ## A2) nearStart & nearEnd:
  ruleA2 <- which(peakFound$nearStart & peakFound$nearEnd)
  # targetDf[unlist(masterIndexDf$nearStart[ruleA2]), ]

  ## set TES target to NULL
  masterIndexDf$nearEnd[ruleA2] <- purrr::map(.x = masterIndexDf$nearEnd[ruleA2],
                                              .f = ~ NULL)
  peakFound$nearEnd[ruleA2] <- FALSE


  ###########
  ## A3) nearStart & peakInFeature:
  ruleA3 <- which(peakFound$nearStart & peakFound$peakInFeature)
  # targetDf[unlist(masterIndexDf$nearStart[ruleA3]), ]


  ###########
  ## A4) nearStart & featureInPeak:
  ruleA4 <- which(peakFound$nearStart & peakFound$featureInPeak)
  # targetDf[unlist(masterIndexDf$nearStart[ruleA4]), ]
  ## no need to filter

  ## featureInPeak peak
  ###########
  ## B1) featureInPeak & nearEnd:
  ruleB1 <- which(peakFound$featureInPeak & peakFound$nearEnd)
  # targetDf[unlist(masterIndexDf$featureInPeak[ruleB1]), ]

  ## remove nearEnd target
  masterIndexDf$nearEnd[ruleB1] <- purrr::map(.x = masterIndexDf$nearEnd[ruleB1],
                                              .f = ~ NULL)
  peakFound$nearEnd[ruleB1] <- FALSE


  ###########
  ## B2) featureInPeak & upstreamTss:
  ruleB2 <- which(peakFound$featureInPeak & peakFound$upstreamTss)
  # targetDf[unlist(masterIndexDf$featureInPeak[ruleB2]), ]

  ruleB2_far <- ruleB2[purrr::map_lgl(.x = peakDistDf$upstreamTss[ruleB2],
                                      .f = function(x){abs(x[1]) > promoterLength})]

  ## set the upstreamTss to NULL if it is far than promoterLength
  masterIndexDf$upstreamTss[ruleB2_far] <- purrr::map(.x = masterIndexDf$upstreamTss[ruleB2_far],
                                                      .f = ~ NULL)
  peakFound$upstreamTss[ruleB2_far] <- FALSE
  ## set bidirectional = 2

  ## for the TF which has known binding over gene body (E.g. polII ChIP or histone marks)
  ## preference is for featureInPeak. all other targets are pseudo
  if(bindingInGene){
    masterIndexDf$upstreamTss[ruleB2] <- purrr::map(.x = masterIndexDf$upstreamTss[ruleB2],
                                                    .f = ~ NULL)
    peakFound$upstreamTss[ruleB2] <- FALSE

  }


  ###########
  ## B3) featureInPeak & peakInFeature:
  ruleB3 <- which(peakFound$featureInPeak & peakFound$peakInFeature)
  # targetDf[unlist(masterIndexDf$featureInPeak[ruleB3]), ]

  ## for the TF which has known binding over gene body (E.g. polII ChIP or histone marks)
  ## set peakInFeature to NULL
  if(bindingInGene){
    masterIndexDf$peakInFeature[ruleB3] <- purrr::map(.x = masterIndexDf$peakInFeature[ruleB3],
                                                      .f = ~ NULL)
    peakFound$peakInFeature[ruleB3] <- FALSE
  } else{
    masterIndexDf$featureInPeak[ruleB3] <- purrr::map(.x = masterIndexDf$featureInPeak[ruleB3],
                                                      .f = ~ NULL)
    peakFound$featureInPeak[ruleB3] <- FALSE
  }

  ###########
  ## C1) upstreamTss & nearEnd:
  ruleC1 <- which(peakFound$upstreamTss & peakFound$nearEnd)
  # targetDf[unlist(masterIndexDf$upstreamTss[ruleC1]), ]

  ## set nearEnd to NULL if upstreamTss within promoter
  ruleC1_far <- ruleC1[purrr::map_lgl(.x = peakDistDf$upstreamTss[ruleC1],
                                      .f = function(x){abs(x[1]) > promoterLength})]

  ## set the upstreamTss to NULL if it is far than promoterLength
  masterIndexDf$upstreamTss[ruleC1_far] <- purrr::map(.x = masterIndexDf$upstreamTss[ruleC1_far],
                                                      .f = ~ NULL)
  peakFound$upstreamTss[ruleC1_far] <- FALSE

  ## ELSE just set the nearEnd peak to pseudo
  ruleC1_pro <- setdiff(x = ruleC1, ruleC1_far)


  ###########
  ## C2) upstreamTss & peakInFeature:
  ruleC2 <- which(peakFound$upstreamTss & peakFound$peakInFeature)
  # targetDf[unlist(masterIndexDf$upstreamTss[ruleC2]), ]

  ## set the upstreamTss to NULL if it is far than promoterLength
  ruleC2_far <- ruleC2[purrr::map_lgl(.x = peakDistDf$upstreamTss[ruleC2],
                                      .f = function(x){abs(x[1]) > promoterLength})]

  masterIndexDf$upstreamTss[ruleC2_far] <- purrr::map(.x = masterIndexDf$upstreamTss[ruleC2_far],
                                                      .f = ~ NULL)
  peakFound$upstreamTss[ruleC2_far] <- FALSE

  ## upstreamTss peaks which are within promoter region
  ruleC2_pro <- setdiff(x = ruleC2, ruleC2_far)

  ## if peakInFeature peak lies near start: set upstreamTss to pseudo
  ruleC2_pro_ovStart <- ruleC2_pro[purrr::map_lgl(
    .x = peakPosDf$peakInFeature[ruleC2_pro],
    .f = function(x){all(x < (1 - insideSkewToEndCut))})]

  markPseudoIdx <- append(markPseudoIdx, unlist(masterIndexDf$upstreamTss[ruleC2_pro_ovStart]))

  ## if peakInFeature peak lies near end: set peakInFeature to pseudo


  ###########
  ## D1) peakInFeature & nearEnd:
  ruleD1 <- which(peakFound$peakInFeature & peakFound$nearEnd)
  # targetDf[unlist(masterIndexDf$peakInFeature[ruleD1]), ]


  ###########
  ## for upstreamTss peak
  ## 8) two upstreamTss peaks: remove pseudo_upstream if peakDist > promoterLength


  targetDf$peakType[markPseudoIdx] <- paste("pseudo_", targetDf$peakType[markPseudoIdx], sep = "")

  ## get the true targets
  trueTargetIdx <- c()
  for (ctg in c("featureInPeak", "nearStart", "nearEnd", "peakInFeature", "upstreamTss")) {
    trueTargetIdx <- append(x = trueTargetIdx,
                            values = unlist(masterIndexDf[[ctg]]))
  }
  trueTargetIdx <- unique(trueTargetIdx)

  peakTargetsGr <- targetDf[trueTargetIdx, ] %>%
    dplyr::select(-rowIdx, -target) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

  return(peakTargetsGr)

}


##################################################################################
