
##################################################################################
# Flow of function calls:
# annotate_peaks()
# |- annotate_ranges()
#    |- splicing_unit_annotations()             ## 5' UTR annotation
#    |- splicing_unit_annotations()             ## 3' UTR region annotations
#    |- splicing_unit_annotations()             ## exons annotations
#    |- splicing_unit_annotations()             ## introns annotations
#    |- region_annotations()                    ## Transcript region annotations
#    |- upstream_annotations()                  ## upstream annotations
#       |- nearest_upstream_bidirectional()
#    |- select_optimal_targets()                ## select correct genes from multiple annotations
#       |- nearest_upstream_bidirectional()
#
#
##################################################################################



#' Annotate narrowPeak using \code{TxDB} object
#'
#' This function annotate the MACS2 called peaks with appropriate target transcript and
#' gene from \code{TxDB} object. Internally uses \code{annotate_ranges()} for Granges
#' annotation. Peaks are annnotated with the following \strong{broad categories} and
#' \emph{specific types} (listed in decreasing order of preference):
#' \enumerate{
#' \item \strong{featureInPeak:} \emph{"include_tx", "include_CDS"}
#' \item \strong{nearStart:} \emph{"5UTR", "CDS_start", "tx_start"}
#' \item \strong{nearEnd:} \emph{"3UTR", "tx_end", "CDS_end"}
#' \item \strong{peakInFeature:} \emph{"exon", "intron", "inside_tx", "inside_CDS"}
#' \item \strong{upstreamTss:} \emph{"promoter", "upstream"}
#' \item \strong{intergenic:} \emph{"intergenic"}
#' }
#' Additionally, a \emph{pseudo} prefix is added to the peakAnnotation where a peak is
#' annotated to two target genes/features and one of it is more optimum than the other.
#' The less optimum target type is prefixed with \emph{pseudo}. Please refer to the
#' \strong{Guidelines} section for specific information on this.
#'
#' @section Guidelines:
#' Some important observations to do before annotating ChIPseq data:
#' \enumerate{
#' \item Whether the signal is sharp peak (normal TF peaks) or over broader region
#' (polII signal over gene body). Also check if binding is throughout the genome like
#' CTCF factor. See \code{bindingInGene, promoterLength} arguments for the details.
#' \item For the genes which are within peak region, what is the gene size (are
#' genes shorter in length than normal) and how far is the next downstream gene.
#' See \code{includeFractionCut} argument for the details.
#' \item Are there any TES or 3' UTR peaks and how confident are they?
#' \item Check the \code{TXTYPE} column in \code{TxDB} object and see which type of
#' features are of interest to you. Usually tRNA, rRNA are not needed.See
#' \code{excludeType} argument for the details.
#' }
#' These observations will help to decide appropriate parameters while annotating
#' the peaks using TxDB object.
#'
#'
#' @param peakFile A narroPeak or broadPeak file. If a broadPeak file, peak center is
#' used as summit as broadPeak file does not report summit
#' @param fileFormat Format of the peak file. One of "narrowPeak" (Default) or "broadPeak".
#' @param output Optionally store the annotation output to a file
#' @param summitRegion Region width around peak summit to use for annotation purpose. This
#' allows peaks with uniform peak width centered around summit. If 0, whole peak region
#' is used. If > 0, it indicates how many basepairs to include upstream and downstream
#' of the peak summit.
#'
#' @inheritParams annotate_ranges
#'
#' @inheritSection annotate_ranges Use of arguments
#'
#' @return A GenomicRanges object with peak annotation
#' @export
#'
#' @examples NA
annotate_peaks <- function(peakFile, fileFormat = "narrowPeak",
                           summitRegion = 0,
                           txdb,
                           promoterLength, upstreamLimit,
                           bidirectionalDistance = 1000, bidirectionalSkew = 0.2,
                           includeFractionCut = 0.7, insideSkewToEndCut = 0.7,
                           txIds = NULL, blacklistRegions = NULL,
                           excludeType = c("tRNA", "rRNA", "snRNA", "snoRNA", "ncRNA"),
                           bindingInGene = FALSE,
                           removePseudo = FALSE,
                           output = NULL){

  ## started working for peak_annotation on larger genomes
  fileFormat <- match.arg(arg = fileFormat, choices = c("narrowPeak", "broadPeak"))

  ## calculate peak related features
  peaks <- rtracklayer::import(con = peakFile, format = fileFormat)

  if(length(peaks) == 0){
    warning("no peak found in peak file ", basename(peakFile))
    return(NULL)
  }

  if(is.null(mcols(peaks)$peak)){
    mcols(peaks)$peak <- as.integer(width(peaks) / 2)
  }

  mcols(peaks)$peakRegion <- paste(
    as.character(seqnames(peaks)), ":", start(peaks), "-", end(peaks), sep = ""
  )

  ## update the peak region used for the annotation.
  if(summitRegion > 0){
    ## start(peaks) + peaks$peak
    peaks <- GenomicRanges::resize(
      x = GenomicRanges::shift(x = peaks, shift = peaks$peak - summitRegion),
      width = summitRegion*2, fix = "start"
    )

    ## update the summit
    mcols(peaks)$peak <- as.integer(width(peaks) / 2)
  }

  ## use main function annotate_ranges() for annotation
  peakTargetsGr <- annotate_ranges(
    peaks = peaks, txdb = txdb,
    promoterLength = promoterLength, upstreamLimit = upstreamLimit,
    txIds = txIds, blacklistRegions = blacklistRegions, excludeType = excludeType,
    bidirectionalDistance = bidirectionalDistance,
    includeFractionCut = includeFractionCut, bindingInGene = bindingInGene,
    insideSkewToEndCut = insideSkewToEndCut, removePseudo = removePseudo
  )

  ## rename narrowPeak/broadPeak specific columns and remove unnecessary columns
  peakTargetsGr <- as.data.frame(peakTargetsGr) %>%
    dplyr::rename(
      peakEnrichment = signalValue,
      peakPval = pValue,
      peakQval = qValue,
      peakId = name
    ) %>%
    dplyr::select(-score) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  ## optionally store the data
  if(!is.null(output)){
    readr::write_tsv(x = as.data.frame(mcols(peakTargetsGr)), path = output)
  }

  return(peakTargetsGr)
}


##################################################################################

#' Annotate GRanges using \code{TxDB} object
#'
#' This is a generic function to annotate GRanges object with TxDB annotations. \cr
#' Regions are annnotated with following \strong{broad categories} and \emph{specific
#' types} (listed in decreasing order of preference):
#' \enumerate{
#' \item \strong{featureInPeak:} \emph{"include_tx", "include_CDS"}
#' \item \strong{nearStart:} \emph{"5UTR", "CDS_start", "tx_start"}
#' \item \strong{nearEnd:} \emph{"3UTR", "tx_end", "CDS_end"}
#' \item \strong{peakInFeature:} \emph{"exon", "intron", "inside_tx", "inside_CDS"}
#' \item \strong{upstreamTss:} \emph{"promoter", "upstream"}
#' \item \strong{intergenic:} \emph{"intergenic"}
#' }
#' Additionally, a \emph{pseudo} prefix is added to the peakAnnotation where a peak is
#' annotated to two target genes/features and one of it is more optimum than other.
#' The less optimum target type is prefixed with \emph{pseudo}.\cr
#' See \strong{Use of arguments} section for more details.
#'
#'
#' @param peaks A GRanges object with name column.
#' @param blacklistRegions A BED file or GRanges object with ChIPseq blacklist regions.
#' Peaks overlapping with these regions are not used for annotation.
#' @param removePseudo Logical: whether to remove peak targets which are marked as pseudo.
#' Default: FALSE
#'
#' @inheritParams get_txdb_transcripts_gr
#' @inheritParams region_annotations
#' @inheritParams upstream_annotations
#' @inheritParams select_optimal_targets
#'
#' @inheritSection upstream_annotations Use of arguments
#' @inheritSection select_optimal_targets Use of arguments
#'
#' @return A GenomicRanges object with peak annotation
#' @export
#'
#' @examples NA
annotate_ranges <- function(peaks, txdb, promoterLength, upstreamLimit,
                            bidirectionalDistance = 1000, bidirectionalSkew = 0.2,
                            includeFractionCut = 0.7, insideSkewToEndCut = 0.7,
                            txIds = NULL, blacklistRegions = NULL,
                            excludeType = c("tRNA", "rRNA", "snRNA", "snoRNA", "ncRNA"),
                            bindingInGene = FALSE,
                            removePseudo = FALSE){


  stopifnot(is(object = peaks, class2 = "GRanges"))
  stopifnot(is(object = txdb, class2 = "TxDb"))

  if(! all(txIds %in% keys(txdb, keytype = "TXID"))){
    stop("unknown TXIDs in txIds: ",
         paste(c(head(txIds[which(!txIds %in% keys(txdb, keytype = "TXID"))]), "..."),
               collapse = " ")
    )
  }

  if(is.null(mcols(peaks)$name)){
    warning("name attribute not found in regions. Using region center...")
    mcols(peaks)$name <- paste("region",1:length(peaks), sep = "_")
  } else if(any(duplicated(mcols(peaks)$name)) ){
    stop("Duplicate values in name column")
  }

  ## rename columns for readability
  mcols(peaks)$peakChr <- as.character(seqnames(peaks))
  mcols(peaks)$peakStart <- start(peaks)
  mcols(peaks)$peakEnd <- end(peaks)


  if(is.null(mcols(peaks)$peak)){
    mcols(peaks)$peak <- as.integer(width(peaks) / 2)
  }

  mcols(peaks)$peakSummit <- GenomicRanges::start(peaks) + mcols(peaks)$peak

  ## ensure that the peak is given as offset from peakStart
  if(any(mcols(peaks)$peakSummit > end(peaks))){
    stop("Peak summit cannot be beyond peakEnd. ",
         "Make sure that peak is a 0 based offset from peakStart (eg. narrowPeak file 10th column)")
  }

  mcols(peaks)$relativeSummitPos <- round((mcols(peaks)$peakSummit - start(peaks)) / width(peaks), 3)

  ## new environment for global variables which takes time to generate
  ## create once and use multiple times
  # txdbEnv <- new.env(parent = emptyenv())

  assign(x = "pointBasedAnnotation", value = FALSE, envir = txdbEnv)
  if(all(width(peaks) < 10)){
    assign(x = "pointBasedAnnotation", value = TRUE, envir = txdbEnv)
  }

  pointBasedAnnotation <- get(x = "pointBasedAnnotation", envir = txdbEnv)

  transcriptsGr <- get_txdb_transcripts_gr(txdb = txdb, excludeType = excludeType,
                                           tx = txIds)
  txToGene <- get(x = "txToGene", envir = txdbEnv)

  ## 5' UTR annotation
  fiveUtrGr <- get_txdb_fiveUtr_gr(txdb = txdb, tx = transcriptsGr$tx_id)
  fiveUtrTargets <- splicing_unit_annotations(peaksGr = peaks, featuresGr = fiveUtrGr,
                                              featureType = "5UTR", txdb = txdb)

  ## 3' UTR region annotations
  threeUtrGr <- get_txdb_threeUtr_gr(txdb = txdb, tx = transcriptsGr$tx_id)
  threeUtrTargets <- splicing_unit_annotations(peaksGr = peaks, featuresGr = threeUtrGr,
                                               featureType = "3UTR", txdb = txdb)

  ## exons annotations
  exonsGr <- get_txdb_exons_gr(txdb = txdb, tx = transcriptsGr$tx_id)
  exonTargets <- splicing_unit_annotations(peaksGr = peaks, featuresGr = exonsGr,
                                           featureType = "exon", txdb = txdb)

  ## introns annotations
  intronsGr <- get_txdb_introns_gr(txdb = txdb, tx = transcriptsGr$tx_id)
  intronTargets <- splicing_unit_annotations(peaksGr = peaks, featuresGr = intronsGr,
                                             featureType = "intron", txdb = txdb)

  ## Transcript region annotations
  transcriptTargets <- region_annotations(peaksGr = peaks,
                                          featuresGr = transcriptsGr,
                                          includeFractionCut = includeFractionCut,
                                          name = "tx")

  ## annotate upstream targets: IMP to give excludeType so that rRNA, tRNA, snRNAs will be removed
  upstreamTargets <- upstream_annotations(peaksGr = peaks, featuresGr = transcriptsGr,
                                          txdb = txdb, promoterLength = promoterLength,
                                          upstreamLimit = upstreamLimit,
                                          bidirectionalDistance = bidirectionalDistance,
                                          bidirectionalSkew = bidirectionalSkew)

  ## prepare target preference list and peak category list
  ## this is internal preference list
  ## this order is IMP for: select 3UTR between 3UTR and inside_tx as it is more specific
  annotationTypes <- data.frame(
    peakAnnotation = c("include_tx", "include_CDS", "5UTR", "CDS_start", "tx_start", "3UTR",
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
                            .f = function(x){data.frame(peakAnnotation = x, stringsAsFactors = F)},
                            .id = "peakCategory")

  annotationTypes <- dplyr::left_join(x = annotationTypes, y = peakCategoryDf, by = "peakAnnotation")

  ## combine annotations
  combinedAnnotations <- NULL
  if(!is.null(fiveUtrTargets)){ combinedAnnotations <- append(combinedAnnotations, fiveUtrTargets) }
  if(!is.null(threeUtrTargets)){ combinedAnnotations <- append(combinedAnnotations, threeUtrTargets) }
  if(!is.null(exonTargets)){ combinedAnnotations <- append(combinedAnnotations, exonTargets) }
  if(!is.null(intronTargets)){ combinedAnnotations <- append(combinedAnnotations, intronTargets) }
  if(!is.null(transcriptTargets)){ combinedAnnotations <- append(combinedAnnotations, transcriptTargets) }
  if(!is.null(upstreamTargets)){ combinedAnnotations <- append(combinedAnnotations, upstreamTargets) }

  if(!is.null(combinedAnnotations)){

    allTargetsDf <- as.data.frame(combinedAnnotations, row.names = NULL) %>%
      dplyr::left_join(y = txToGene, by = c("tx_id" = "TXID")) %>%
      dplyr::left_join(y = annotationTypes, by = c("peakAnnotation" = "peakAnnotation"))

    ## extract best transcript region annotation for each peak-gene combination using the preference
    ## using geneId instead of tx_id: for gene which have multiple tx, 5UTR/3UTR/Intron can get
    ## annotated with any tx
    bestPeakGeneTargets <- allTargetsDf %>%
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
    #                    peakAnnotation = paste(sort(peakAnnotation), collapse = ","),
    #                    peakCat = paste(peakCategory, collapse = ","),
    #                    uniqCat = paste(sort(unique(peakCategory)), collapse = ","),
    #                    gene = paste(geneId, collapse = ",")) %>%
    #   dplyr::filter(n > 1) %>%
    #   dplyr::distinct(peakAnnotation, .keep_all = T) %>%
    #   as.data.frame()
    #######################

    tempTargetGrl <- GenomicRanges::split(x = bestPeakGeneTargetsGr, f = mcols(bestPeakGeneTargetsGr)$name)
    # table(elementNROWS(tempTargetGrl))
    singleHitPeaks <- unlist(tempTargetGrl[which(elementNROWS(tempTargetGrl) == 1)], use.names = F)
    multipleHitPeaksGrl <- tempTargetGrl[which(elementNROWS(tempTargetGrl) != 1)]

    # #######################
    # optGr <- select_optimal_targets(
    #   targetGr = tempTargetGrl$An_kdmB_20h_HA_1_withCtrl_peak_1288,
    #   promoterLength = promoterLength,
    #   bindingInGene = bindingInGene,
    #   insideSkewToEndCut = insideSkewToEndCut)
    # #######################

    if(length(multipleHitPeaksGrl) > 0){
      ## for each peak, find optimum target/s
      multipleHitPeaks <- select_optimal_targets(
        targetGr = sort(unlist(multipleHitPeaksGrl, use.names = F)),
        promoterLength = promoterLength,
        upstreamLimit = upstreamLimit,
        bindingInGene = bindingInGene,
        insideSkewToEndCut = insideSkewToEndCut,
        bidirectionalSkew = bidirectionalSkew,
        bidirectionalDistance = bidirectionalDistance)

    } else{
      multipleHitPeaks <- NULL
    }


    peakTargetsGr <- c(singleHitPeaks, multipleHitPeaks, ignore.mcols=FALSE)

    ## add the unannotated peaks
    unannotatedPeaks <- peaks[which(!peaks$name %in% peakTargetsGr$name)]
    if(length(unannotatedPeaks) > 0){
      mcols(unannotatedPeaks)$peakAnnotation <- "intergenic"
      mcols(unannotatedPeaks)$peakCategory <- "intergenic"
    }


    peakTargetsGr <- c(peakTargetsGr, unannotatedPeaks, ignore.mcols=FALSE)

    ## set upstream peaks with distance > upstreamLimit to upstream_intergenic
    farUpstreamPeaks <- which(
      peakTargetsGr$peakAnnotation == "upstream" & abs(peakTargetsGr$peakDist) > upstreamLimit
    )

    peakTargetsGr$peakAnnotation[farUpstreamPeaks] <- "upstream_intergenic"
    peakTargetsGr$peakCategory[farUpstreamPeaks] <- "intergenic"


    ## optionally filter peak targets which are marked as pseudo
    if(removePseudo){
      peakTargetsGr <- peakTargetsGr[!grepl(pattern = "pseudo_", x = mcols(peakTargetsGr)$peakAnnotation)]
    }

  } else{
    ## use the original peakset
    peakTargetsGr <- peaks

    mcols(peakTargetsGr)$geneId <- NA
    mcols(peakTargetsGr)$txName <- NA
    mcols(peakTargetsGr)$peakAnnotation <- NA
    mcols(peakTargetsGr)$peakCategory <- NA
    mcols(peakTargetsGr)$peakPosition <- NA
    mcols(peakTargetsGr)$peakDist <- NA
    mcols(peakTargetsGr)$summitDist <- NA
    mcols(peakTargetsGr)$bidirectional <- NA
    mcols(peakTargetsGr)$relativePeakPos <- NA
    mcols(peakTargetsGr)$targetOverlap <- NA
    mcols(peakTargetsGr)$peakOverlap <- NA

  }


  peakTargetsGr <- sort(peakTargetsGr)
  names(x = peakTargetsGr) <- NULL

  ## make sure the peakChr, peakStart and peakEnd columns are selected first
  peakFileCols <- c("peakChr", "peakStart", "peakEnd")
  peakFileCols <- union(peakFileCols, names(mcols(peaks)))

  ## select output colums and rename columns to standard column names:
  annotationGr <- as.data.frame(peakTargetsGr) %>%
    dplyr::select(
      seqnames, start, end, !!! peakFileCols,
      geneId, txName, peakAnnotation, peakCategory, peakPosition, peakDist, summitDist,
      bidirectional, relativePeakPos, targetOverlap, peakOverlap
    ) %>%
    dplyr::rename(
      txId = txName
    ) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  ## remove unnecessary columns from query GRanges
  mcols(annotationGr)$peak <- NULL

  return(annotationGr)
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
#' peakAnnotation, targetOverlap, peakOverlap, peakDist, summitDist, relativePeakPos}
#' @export
#'
#' @examples NA
splicing_unit_annotations <- function(peaksGr, featuresGr, featureType, txdb){
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
  mcols(peakTargets)$peakAnnotation <- featureType
  mcols(peakTargets)$peakDist <- 0

  ## targetOverlap is at transcript level
  txSubsetGrl <- GenomicFeatures::mapIdsToRanges(
    x = txdb,
    keys = list(tx_id = mcols(peakTargets)$tx_id),
    type = "tx", columns = c("gene_id")
  )

  txSubsetGr <- unlist(txSubsetGrl)

  mcols(peakTargets)$targetStart = start(txSubsetGr)
  mcols(peakTargets)$targetEnd = end(txSubsetGr)
  mcols(peakTargets)$targetStrand = strand(txSubsetGr)
  mcols(peakTargets)$txWidth = width(txSubsetGr)
  mcols(peakTargets)$gene_id = unlist(mcols(txSubsetGr)$gene_id)
  mcols(peakTargets)$targetOverlap <- round(
    width(pintersect(x = peakTargets, y = txSubsetGr)) / width(txSubsetGr), 3
  )
  mcols(peakTargets)$peakOverlap <- round(
    width(pintersect(x = peakTargets, y = txSubsetGr)) / width(peakTargets), 3
  )

  ## calculate relativePeakPos using peak summit position
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
      relativePeakPos = dplyr::case_when(
        summitDist >= (targetEnd - targetStart) ~ 1,
        summitDist <= 0 ~ 0,
        summitDist > 0 ~ round((summitDist) / (targetEnd - targetStart), 3),
        TRUE ~ 0
      )
    ) %>%
    dplyr::mutate(
      bidirectional = NA,
      relativeSummitPos = dplyr::if_else(
        condition = targetStrand == "-",
        true = round(1 - relativeSummitPos, 3), false = relativeSummitPos)
    )

  if(featureType == "5UTR"){
    ## ensure that for 5UTR region peak, peak does not go beyond target end
    targetsDf <- dplyr::filter(
      targetsDf,
      (targetStrand == "+" & peakEnd < targetEnd) |
        (targetStrand == "-" & peakStart > targetStart)
    )
  } else if(featureType == "3UTR"){
    ## ensure that for 3UTR region peak, peak does not start before target start
    targetsDf <- dplyr::filter(
      targetsDf,
      (targetStrand == "+" & peakStart > targetStart) |
        (targetStrand == "-" & peakEnd < targetEnd)
    )
  }

  targetsGr <- sort(makeGRangesFromDataFrame(df = targetsDf, keep.extra.columns = TRUE))
  return(targetsGr)
}


##################################################################################

#' Map peaks to given GRanges regions
#'
#' This function annotates the peaks onto a regions into e.g. \code{tx_start, tx_end,
#' include_tx, inside_tx} categories.
#' In addition, relativeSummitPos value is updated w.r.t. region for the peaks which are
#' annotated as e.g. \code{inside_tx}
#'
#' @param peaksGr GRanges object generated from narrowPeak or broadPeak file
#' @param featuresGr GRanges object for regions on which peaks needs to be mapped. E.g:
#'  \code{GenomicFeatures::cdsBy()} or \code{GenomicFeatures::transcripts()}
#' @param includeFractionCut A floating point number between [0, 1]. If a peak covers more
#' than this proportion of feature, it is annotated as, eg. \code{include_tx}. Default: 0.7
#' @param name Feature type to be used as suffix in peak type annotation. Eg. CDS, gene etc.
#'
#' @return A modified peak GRanges object with additional columns: \code{tx_id,
#' peakAnnotation, targetOverlap, peakOverlap, peakDist, summitDist, relativePeakPos}
#' @export
#'
#' @examples NA
region_annotations <- function(peaksGr, featuresGr, includeFractionCut = 0.7, name = "tx"){

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
  mcols(queryTargets)$peakAnnotation <- insideFeature
  mcols(queryTargets)$peakDist <- 0

  ##
  mcols(queryTargets)$targetStart = start(featuresGr[ovlpHits@to])
  mcols(queryTargets)$targetEnd = end(featuresGr[ovlpHits@to])
  mcols(queryTargets)$targetStrand = strand(featuresGr[ovlpHits@to])
  mcols(queryTargets)$targetOverlap <- round(
    width(pintersect(x = queryTargets, y = featuresGr[ovlpHits@to])) / width(featuresGr[ovlpHits@to]), 3
  )
  mcols(queryTargets)$peakOverlap <- round(
    width(pintersect(x = queryTargets, y = featuresGr[ovlpHits@to])) / width(queryTargets), 3
  )

  ## assign appropriate relativePeakPos using peak summit (shown as *)
  ##
  ##     |>=====>=====>======>======>======>|           |<=====<=====<======<======<======<|
  ##   -*----- 0                                                                     0 -----*--
  ##                                    -----*- 1     ----*--0.97
  ##                 ---*---0.45                                    ----*--- 0.55
  ##  -------------------------------------------    -------------------------------------------
  ##
  ## calculate relativePeakPos using peak summit position
  ## summitDist > (targetEnd - targetStart) : peak overlap at end and summit is after end
  ## summitDist < 0 : peak overlap at start and summit before start
  targetsDf <- as.data.frame(queryTargets) %>%
    dplyr::mutate(
      peakAnnotation = dplyr::case_when(
        peakStart <= targetStart & peakEnd >= targetEnd ~ includeFeature,
        targetOverlap >= includeFractionCut ~ includeFeature,
        targetStrand == "+" & peakStart <= targetStart & peakEnd < targetEnd ~ overlapStart,
        targetStrand == "+" & targetStart < peakStart & targetEnd <= peakEnd ~ overlapEnd,
        targetStrand == "-" & peakStart <= targetStart & peakEnd < targetEnd ~ overlapEnd,
        targetStrand == "-" & targetStart < peakStart & targetEnd <= peakEnd ~ overlapStart,
        TRUE ~ insideFeature
      )
    ) %>%
    dplyr::mutate(
      summitDist = dplyr::case_when(
        targetStrand == "+" ~ peakSummit - targetStart,
        targetStrand == "-" ~ targetEnd - peakSummit
      )
    ) %>%
    dplyr::mutate(
      relativePeakPos = dplyr::case_when(
        summitDist >= (targetEnd - targetStart) ~ 1,
        summitDist <= 0 ~ 0,
        summitDist > 0 ~ round((summitDist) / (targetEnd - targetStart), 3),
        TRUE ~ 0
      )
    ) %>%
    dplyr::mutate(
      bidirectional = NA,
      relativeSummitPos = dplyr::if_else(
        condition = targetStrand == "-",
        true = round(1 - relativeSummitPos, 3), false = relativeSummitPos)
    )

  ## convert back to GRanges
  queryTargets <- makeGRangesFromDataFrame(df = targetsDf, keep.extra.columns = TRUE)

  return(queryTargets)
}

##################################################################################

#' Set peakAnnotation to pseudo
#'
#' @param target A dataframe or GRanges object which has peakAnnotation column
#'
#' @return Same object with \code{pseudo} prefix to peakAnnotation column values
#' @export
#'
#' @examples NA
set_peakTarget_to_pseudo <- function(target){
  if(any(class(target) %in% "GRanges")){
    mcols(target)$peakAnnotation <- paste("pseudo_", mcols(target)$peakAnnotation, sep = "")
  } else if(any(class(target) %in% "data.frame")){
    target$peakAnnotation <- paste("pseudo_", target$peakAnnotation, sep = "")
  }

  return(target)
}

##################################################################################


#' Annotate upstream peaks on transcripts
#'
#' This function annotates the peaks with nearest downstream target.
#' See \strong{Use of arguments} section for more details.
#'
#' @section Use of arguments:
#' \subsection{upstreamOverlappingFraction}{
#' There will be cases when a peak is inside a gene and it is upstream of other gene.
#' \preformatted{
#' #                                                                            #
#' #           target1                     target2                              #
#' #         =====<=======<===       =====<=======<========<=======             #
#' #                                   -^--           --^-                      #
#' #                                 peak1           peak2                      #
#' #                         |<------>|                                         #
#' #                         |<------------------------->|                      #
#' #                                                                            #
#' In above cases, peak1 can be annotated as Upstream of target1. However not peak2
#' because target2 has bigger fraction in-between [target1, peak2] range
#'
#' Target gene inside gene case:
#' #                                                                            #
#' #         =====<=======<=======<=======<========<======= target2             #
#' #          target1 ==<==                                                     #
#' #                                  -^--           --^-                      #
#' #                                 peak1           peak2                      #
#' #                       |<-------->|                                         #
#' #                                                                            #
#' }
#' In above case, peak1 is inside targe2 and upstream of target1. These targets are
#' selected in \code{select_optimal_targets()} if peak lies within promoter range.
#' }
#'
#'
#' @param peaksGr GRanges object for peak data
#' @param featuresGr A transcript GRanges object
#' @param txdb Optional TxDB object. A 3UTR less region is created using TxDB and used
#' to find the overlap between peak and its downstream target. A max fractional
#' overlap of 0.2 with geneA is allowed in a case when peak overlaps with a geneA and
#' is upstream of geneB. This is useful for the peaks which are near TES of a geneA.
#' If TxDB object is not provided, featuresGr is used. Default: featuresGr is used.
#' @param upstreamOverlappingFraction See details. Default: 0.2
#'
#' @inheritParams nearest_upstream_bidirectional
#'
#' @inheritSection nearest_upstream_bidirectional Use of arguments
#'
#' @return A modified peak GRanges object with additional columns: \code{ tx_id,
#' peakAnnotation, targetOverlap, peakOverlap, peakDist, summitDist, relativePeakPos}
#' @export
#'
#' @examples NA
upstream_annotations <- function(peaksGr, featuresGr, txdb = NULL,
                                 upstreamOverlappingFraction = 0.2,
                                 promoterLength, upstreamLimit,
                                 bidirectionalDistance, bidirectionalSkew){

  stopifnot(is(object = peaksGr, class2 = "GRanges"))
  stopifnot(is(object = featuresGr, class2 = "GRanges"))

  if(exists("pointBasedAnnotation", envir=txdbEnv, inherits=FALSE)) {
    pointBasedAnnotation <- get(x = "pointBasedAnnotation", envir = txdbEnv)
  } else{
    stop("pointBasedAnnotation not found in txdbEnv")
  }

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
    expectedFeatureStrand = "+",
    featureStrand = as.vector(strand(featuresGr))[peakDownFeatures],
    txName = featuresGr$tx_name[peakDownFeatures],
    stringsAsFactors = FALSE)

  ## select immediate upstream feature to the peak
  peakUpFeatures <- GenomicRanges::follow(x = peaksGr, subject = featuresGr,
                                          select = "last", ignore.strand = TRUE)

  ## select the upstreamHits which have strand == "-"
  ## because peak is downstream to the feature with strand == "+"
  peakUpHits <- data.frame(
    from = 1:length(peaksGr), peakId = mcols(peaksGr)$name,
    to = peakUpFeatures,
    expectedFeatureStrand = "-",
    featureStrand = as.vector(strand(featuresGr))[peakUpFeatures],
    txName = featuresGr$tx_name[peakUpFeatures],
    stringsAsFactors = FALSE)


  targetsAroundPeak <- dplyr::bind_rows(peakDownHits, peakUpHits) %>%
    dplyr::filter(!is.na(to))

  ## remove false downstream targets
  upstreamHits <- dplyr::filter(targetsAroundPeak, expectedFeatureStrand == featureStrand)

  ####################
  ## if only one gene is found i.e. either upstream of peak or downstream of peak,
  ## additionally select one immediate downstream gene in opposite direction.
  ## In following case, gene0 is found by follow() method which is upstream of peak1
  ## and is possible target. Gene2 can also be a possible target but precede() will
  ## select gene1 (ignore.strand = TRUE) and it will be filtered ultimately.
  ## So we select on gene in opposite direction with ignore.strand = TRUE. This way
  ## we will have two targets on each side of peak. if its case like gene1 and gene3,
  ## gene3 will be filtered using peak-gap overlap filter code below
  ##                           ====<=======<======<====== gene1
  ##                             ====>=====>======>======= gene2
  ##  ====<======<=                                    ====>======> gene3
  ##                   -----
  ##     gene0            peak1
  ##

  ## find if any other target overlap with false downstream target (gene1) which can
  ## be possible true downstream target
  falseTargets <- dplyr::filter(targetsAroundPeak, expectedFeatureStrand != featureStrand)

  if (nrow(falseTargets) > 0) {

    falseTargetsGr <- featuresGr[falseTargets$to]
    mcols(falseTargetsGr)$from  <- falseTargets$from
    mcols(falseTargetsGr)$to <- falseTargets$to
    mcols(falseTargetsGr)$peakId <- falseTargets$peakId
    mcols(falseTargetsGr)$expectedFeatureStrand <- falseTargets$expectedFeatureStrand

    falseTargetOverlaps <- GenomicRanges::findOverlaps(query = falseTargetsGr,
                                                       subject = featuresGr,
                                                       ignore.strand = TRUE)

    otherUpstream <- tibble::tibble(
      from = mcols(falseTargetsGr)$from[falseTargetOverlaps@from],
      peakId =  mcols(falseTargetsGr)$peakId[falseTargetOverlaps@from],
      to = falseTargetOverlaps@to,
      expectedFeatureStrand = mcols(falseTargetsGr)$expectedFeatureStrand[falseTargetOverlaps@from],
      featureStrand = as.vector(strand(featuresGr))[falseTargetOverlaps@to],
      txName = featuresGr$tx_name[falseTargetOverlaps@to],
      stringsAsFactors = FALSE) %>%
      dplyr::filter(expectedFeatureStrand == featureStrand)

    ## some of above otherUpstream target can overlap with peak. remove such targets
    ## as it is handled by splicing_unit_annotations() for 5' UTR overlap
    notOvlpWithPeak <- distance(x = peaksGr[otherUpstream$from], y = featuresGr[otherUpstream$to],
                                ignore.strand = TRUE) != 0

    otherUpstream <- otherUpstream[notOvlpWithPeak, ]

  } else{
    otherUpstream <- NULL
  }
  ####################


  ## combine putative upstream hits
  upstreamHits <- dplyr::bind_rows(upstreamHits, otherUpstream) %>%
    dplyr::arrange(from) %>%
    dplyr::mutate(hitId = row_number())

  ## return null if nothing found
  if(nrow(upstreamHits) == 0){
    return(NULL)
  }

  ## find the number of genes between peak and its target.
  ## only those targets are true where there is no other feature inbetween
  ## build a GRanges object of the gap region between peak and target gene

  ## generate peak-target gap GRanges
  peakTargetGapsGr <- GenomicRanges::pgap(x = peaksGr[upstreamHits$from],
                                          y = featuresGr[upstreamHits$to])

  peakTargetGapsGr <- unstrand(peakTargetGapsGr)
  names(peakTargetGapsGr) <- upstreamHits$hitId
  mcols(peakTargetGapsGr)$hitId <- upstreamHits$hitId

  ## build a subject GRanges for tx - 3UTR
  ## Such custom regions are used because genes have very long 3' UTR.
  ## A peak in UTR region of a gene can be upstream of another gene
  featureIdx <- tibble::tibble(
    tx_id = mcols(featuresGr)$tx_id, featureIndex = 1:length(featuresGr)
  )

  if(!is.null(txdb)){
    stopifnot(is(object = txdb, class2 = "TxDb"))

    fiveUtrGr <- get_txdb_fiveUtr_gr(txdb = txdb)
    threeUtrGr <- get_txdb_threeUtr_gr(txdb = txdb)

    threeUtrIdx <- tibble::tibble(
      tx_id = mcols(threeUtrGr)$tx_id,
      threeUtrIndex = 1:length(threeUtrGr)
    ) %>%
      dplyr::left_join(y = featureIdx, by = "tx_id")

    txMinus3Utrs <- GenomicRanges::psetdiff(
      x = featuresGr[threeUtrIdx$featureIndex],
      y = threeUtrGr[threeUtrIdx$threeUtrIndex],
      ignore.strand= TRUE)

    txMinusUtrs <- sort(
      c(txMinus3Utrs,
        featuresGr[setdiff(featureIdx$featureIndex, threeUtrIdx$featureIndex)],
        ignore.mcols = TRUE)
    )

    ## no need to remove 5UTR
  }

  ## find overlapping gene/features in gap GRanges.
  featuresInGap <- GenomicRanges::findOverlaps(query = peakTargetGapsGr,
                                               subject = txMinusUtrs,
                                               ignore.strand = TRUE)

  isFeatureInBetweenDf <- tibble::tibble(
    hitId = mcols(peakTargetGapsGr)$hitId[featuresInGap@from],
    gapWidth = width(peakTargetGapsGr)[featuresInGap@from],
    gapGrRow = featuresInGap@from,
    firstOverlapFeature = featuresInGap@to
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
  isFeatureInBetweenDf$ovlpGapFrac <- isFeatureInBetweenDf$intersectWd / isFeatureInBetweenDf$gapWidth

  isFeatureInBetweenDf <- dplyr::group_by(isFeatureInBetweenDf, hitId) %>%
    dplyr::arrange(desc(ovlpFeatureFrac), .by_group = TRUE) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup()

  ## upstreamOverlappingFraction based filtering OR
  ## select if target is within promoterLength distance
  ## 0.2 is still very big for the large genomes such as human, mouse as genes are very long
  upstreamHitsFiltered <- dplyr::left_join(x = upstreamHits, y = isFeatureInBetweenDf, by = "hitId") %>%
    tidyr::replace_na(list(intersectWd = 0, ovlpFeatureWd = 0, ovlpFeatureFrac = 0, ovlpGapFrac = 0)) %>%
    dplyr::filter(ovlpFeatureFrac <= upstreamOverlappingFraction | gapWidth < promoterLength)

  ## if no upstream peak after filtering
  if(nrow(upstreamHitsFiltered) == 0){
    return(NULL)
  }

  ## build upstream peaks data and filter unnecessary peaks where there is/are genes between peak and target
  upstreamPeaks <- peaksGr[upstreamHitsFiltered$from]

  mcols(upstreamPeaks)$tx_id <- mcols(featuresGr)$tx_id[upstreamHitsFiltered$to]
  mcols(upstreamPeaks)$peakAnnotation <- "upstream"
  mcols(upstreamPeaks)$peakDist <- GenomicRanges::distance(
    x = upstreamPeaks,
    y = featuresGr[upstreamHitsFiltered$to]
  )

  mcols(upstreamPeaks)$peakDist <- mcols(upstreamPeaks)$peakDist * -1
  mcols(upstreamPeaks)$targetStart = start(featuresGr[upstreamHitsFiltered$to])
  mcols(upstreamPeaks)$targetEnd = end(featuresGr[upstreamHitsFiltered$to])
  mcols(upstreamPeaks)$targetStrand = strand(featuresGr[upstreamHitsFiltered$to])

  ## calculate summit distance and change relativeSummitPos based on target gene
  upstreamPeaks <- as.data.frame(upstreamPeaks, stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      targetOverlap = 0,
      peakOverlap = 0,
      summitDist = dplyr::case_when(
        targetStrand == "+" ~ peakSummit - targetStart,
        targetStrand == "-" ~ targetEnd - peakSummit
      )
    ) %>%
    dplyr::mutate(
      relativePeakPos = 0,
      bidirectional = NA,
      relativeSummitPos = dplyr::if_else(
        condition = targetStrand == "-",
        true = round(1 - relativeSummitPos, 3), false = relativeSummitPos),
      peakAnnotation = if_else(abs(peakDist) < promoterLength, "promoter", peakAnnotation)
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

  dualTargetFiltered <- NULL

  if(length(dualTargetPeaks) > 0){

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

    pseudoUpIdx <- nearest_upstream_bidirectional(
      targetDf = dualTargetPeaksDf,
      t1Idx = pairTable$t1,
      t2Idx = pairTable$t2,
      promoterLength = promoterLength,
      upstreamLimit = upstreamLimit,
      bidirectionalDistance = bidirectionalDistance,
      pointBasedAnnotation = pointBasedAnnotation)

    correctUpstream <- setdiff(x = c(pairTable$t1, pairTable$t2), y = pseudoUpIdx)

    bidirectAIdx <- purrr::map2(
      .x = pairTable$t1, .y = pairTable$t2,
      .f = function(x, y){
        if(any(c(x, y) %in% pseudoUpIdx)){
          return(NULL)
        } else{
          return(c(x, y))
        }
      }) %>%
      purrr::flatten_int()

    dualTargetPeaksDf$bidirectional[bidirectAIdx] <- "A"

    ## remove the targets which are too far based on nearest_upstream_bidirectional()
    dualTargetFiltered <- dualTargetPeaksDf[correctUpstream, ] %>%
      dplyr::select(-rowIdx, -target) %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  }


  upstreamPeaksAn <- sort(c(unlist(singalTargetPeaks, use.names = FALSE),
                            dualTargetFiltered))
  names(upstreamPeaksAn) <- NULL

  return(upstreamPeaksAn)

}


##################################################################################



#' True target for bidirectional peak
#'
#' In case of peak at bidirectional promoter, this function assigns one or both genes
#' as peak target. See \strong{Use of arguments} section for more details.
#'
#' @section Use of arguments:
#' \subsection{bidirectionalDistance and bidirectionalSkew}{
#' For peak at bidirectional promoter, its skewness from midpoint and summit poistion
#' w.r.t. midpoint is used to decide correct target/s. \cr
#' \code{nearest_upstream_bidirectional(..., bidirectionalDistance = 1000,
#' bidirectionalSkew = 0.2)}
#' \preformatted{
#' #                        *                                                   #
#' #                    |<-------more than 500bp------>|                        #
#' #         target1                    |                    target2            #
#' #  ==<=====<=====<===                |               ===>===>====>====>==    #
#' #                              --^---|- peak1                                #
#' #                                 ---|^--- peak2                             #
#' #                                  --|^------ peak3                          #
#' #                peak4  --^--        |                                       #
#' #                                    |                                       #
#' #                                  { | } central 10\% region                 #
#' #                                    |                                       #
#' #                        midpoint between two targets                        #
#' #                                                                            #
#' peak1 => target1: more than 80% of the peak lies on target1 side from TSS midpoint
#'          and peak summit is not in central 20\% region.
#' peak2 => target1, targe2: peak lies on the center
#' peak3 => Additionally, if peak summit (^) is within central 10\% region
#'          (= bidirectionalSkew/2), peak is assigned to both the target genes
#' peak4 => target1
#' If gap between two bidirectional genes < \code{bidirectionalDistance}
#' (default: 1000bp), peak is assigned to both the target genes
#' }
#' }
#'
#'
#' @param targetDf A dataframe which has bidirectional targets for each peak.
#' @param t1Idx target1 index vector
#' @param t2Idx target2 index vector. IMP: \code{length(t1Idx)} should be equal to
#' \code{length(t2Idx)}
#' @param bidirectionalSkew Maximum fraction of peak region allowed on the side of
#' false target from the midpoint of two target genes. Default: 0.2
#' @param bidirectionalDistance When a peak is present at bidirectional promoter where
#' distance between two TSS is < \code{bidirectionalDistance}, any gene within
#' \code{promoterLength} distance of peak is assigned to the peak as annotation. Default: 1000
#' @param promoterLength Promoter region length. Upstream peaks within \code{promoterLength}
#' distance of feature start are annotated as \code{promoter} region peaks.
#' @param upstreamLimit Maximum distance of peak for upstream annotation. Peak beyond
#' this distance can be considered as intergenic instead.
#' @param pointBasedAnnotation Logical: whether peak annotation is based on just
#' the summit or whole peak region.
#'
#' @return A vector of row index for pseudo targets.
#' @export
#'
#' @examples NA
nearest_upstream_bidirectional <- function(targetDf, t1Idx, t2Idx,
                                           promoterLength, upstreamLimit,
                                           bidirectionalSkew = 0.2,
                                           bidirectionalDistance = 1000,
                                           pointBasedAnnotation = FALSE){

  targetPairDf <- tibble::tibble(t1Idx = t1Idx, t2Idx = t2Idx,
                                 t1Select = TRUE, t2Select = TRUE,
                                 dir = "opposite")

  targetA <- targetDf[targetPairDf$t1Idx, ]
  targetB <- targetDf[targetPairDf$t2Idx, ]

  if(!all(targetA$name == targetB$name)){
    stop("peak name mismatch in bidirectional target pair")
  }

  targetPairDf$peakId <- targetA$name
  targetPairDf$peakSummit <- targetA$peakSummit
  targetPairDf$t1PeakDist <- targetA$peakDist
  targetPairDf$t2PeakDist <- targetB$peakDist
  targetPairDf$t1SummitDist <- targetA$summitDist
  targetPairDf$t2SummitDist <- targetB$summitDist

  peakGr <- GenomicRanges::makeGRangesFromDataFrame(df = targetA)
  targetPairDf$peakWidth <- width(peakGr)
  targetPairDf$peakFraction <- round(targetPairDf$peakWidth * bidirectionalSkew)

  sameDirTargets <- which(targetA$targetStrand == targetB$targetStrand)
  targetPairDf$dir[sameDirTargets] <- "same"


  targetAGr <- dplyr::select(targetA, seqnames, targetStart, targetEnd, targetStrand, name, rowIdx) %>%
    makeGRangesFromDataFrame(start.field = "targetStart", end.field = "targetEnd",
                             strand.field = "targetStrand",
                             keep.extra.columns = TRUE)

  targetBGr <- dplyr::select(targetB, seqnames, targetStart, targetEnd, targetStrand, name, rowIdx) %>%
    makeGRangesFromDataFrame(start.field = "targetStart", end.field = "targetEnd",
                             strand.field = "targetStrand",
                             keep.extra.columns = TRUE)

  zeroGapTargets <- which(distance(x = targetAGr, y = targetBGr, ignore.strand = TRUE) == 0)
  # if(any(distance(x = targetAGr, y = targetBGr, ignore.strand = TRUE) == 0)){
  #   stop("Distance between bidirectional target genes cannot be 0")
  # }

  ## gap between two bidirectional targets
  targetGapGr <- unstrand(pgap(x = targetAGr, y = targetBGr, ignore.strand = TRUE))
  targetPairDf$gapWidth <- width(targetGapGr)
  targetPairDf$midpoint <- start(targetGapGr) + round(width(targetGapGr)/2)

  ## midpoint of gap between bidirectional targets
  # midpoint <- resize(x = targetGapGr, width = 1, fix = "center", ignore.strand = TRUE)
  # targetPairDf$midpointDist <- distance(x = midpoint, y = peakGr, ignore.strand = TRUE)

  centralNoSummitZone <- bidirectionalSkew
  targetPairDf$summitPosLim <- (targetPairDf$gapWidth + targetPairDf$gapWidth*centralNoSummitZone)/2

  ## gapWidth condition has to be 1st always to ensure that the targets which are
  ## very close are not marked pseudo
  targetPairDf <- targetPairDf %>% dplyr::mutate(
    t1Select = dplyr::case_when(
      gapWidth <= bidirectionalDistance & abs(t1PeakDist) < promoterLength ~ TRUE,
      # pointBasedAnnotation == FALSE & abs(t1PeakDist) < promoterLength ~ TRUE,
      abs(t1PeakDist) > upstreamLimit & abs(t2PeakDist) <= promoterLength ~ FALSE,
      dir == "opposite" & pointBasedAnnotation == FALSE &
        (abs(t1PeakDist) + peakFraction) > (gapWidth/2) ~ FALSE,
      dir == "opposite" & pointBasedAnnotation &
        abs(t1SummitDist) > summitPosLim ~ FALSE,
      TRUE ~ t1Select
    ),
    t2Select = dplyr::case_when(
      gapWidth <= bidirectionalDistance & abs(t2PeakDist) < promoterLength ~ TRUE,
      # pointBasedAnnotation == FALSE & abs(t2PeakDist) < promoterLength ~ TRUE,
      abs(t2PeakDist) > upstreamLimit & abs(t1PeakDist) <= promoterLength ~ FALSE,
      dir == "opposite" & pointBasedAnnotation == FALSE &
        (abs(t2PeakDist) + peakFraction) > (gapWidth/2) ~ FALSE,
      dir == "opposite" & pointBasedAnnotation &
        abs(t2SummitDist) > summitPosLim ~ FALSE,
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
#' optimum target gene/s. See \strong{Use of arguments} section for more details.
#'
#' @section Use of arguments:
#' \subsection{insideSkewToEndCut}{
#' Use of \code{select_optimal_targets(..., insideSkewToEndCut = 0.7, promoterLength = 500)}
#' \preformatted{
#' #                                                                            #
#' #                  d<500bp                                                   #
#' #      target1    |<---->|    target2           |<--500-->|      target3     #
#' #                        0  0.25  0.5  0.75  1                               #
#' #    ====<====<===       =====>=====>=====>=~~~~           ====>====>===     #
#' #                         --^--            --^--                             #
#' #                         peak1            peak2                             #
#' #                        |<----0.7---->|                                     #
#' #                        **                                                  #
#' #                                                                            #
#' In first example, peak1 is inside target2 and it is near the start of target2.
#' Even though peak1-target1 distance is < 500bp, target1 is marked as pseudo here.
#' In another example, peak2 is inside target2 but it is near the end at 3'UTR.
#' Relative position of the peak2 is >0.7 in target2. peak2 is also upstream of
#' target3 and within (upstreamLimit)bp. So correct annotations are:
#' peak1: target2
#' peak2: target2, target3
#' }
#' }
#'
#'
#' @param targetGr peak annotation in form of GRanges object. This is generated by
#' combining \code{UTR_annotate(), region_annotations(), upstream_annotations()}
#' functions.
#' @param insideSkewToEndCut A floating point number in range [0, 1]. If a peak is
#' present inside feature/gene and the relative summit position is >
#' \code{insideSkewToEndCut}, it is closer to the end of the feature. Default: 0.7
#' @param bindingInGene Logical: whether the ChIPseq TF binds in gene body. This is
#' useful for polII ChIPseq data. Default: FALSE
#'
#' @inheritParams nearest_upstream_bidirectional
#'
#' @return Same GRanges object as peakGr but with modified peakAnnotation column or excluding
#' some targets which were not optimum.
#' @export
#'
#' @examples NA
select_optimal_targets <- function(targetGr, promoterLength, upstreamLimit,
                                   bindingInGene, insideSkewToEndCut,
                                   bidirectionalSkew, bidirectionalDistance){

  pointBasedAnnotation <- get(x = "pointBasedAnnotation", envir = txdbEnv)

  ## important to sort the targets based on peakDist and relativePeakPos
  ## rowIdx will be used later to extract correct targets
  targetDf <- as.data.frame(targetGr) %>%
    dplyr::arrange(seqnames, start, preference, peakDist, relativePeakPos) %>%
    dplyr::mutate(rowIdx = 1:n(),
                  target = TRUE)

  subData <- targetDf %>%
    dplyr::select(seqnames, start, end, name, peakCategory, peakDist, summitDist,
                  relativePeakPos, targetOverlap, peakOverlap, rowIdx, target)

  masterIndexDf <- tibble::tibble(name = unique(subData$name))
  peakFound <- masterIndexDf
  peakDistDf <- masterIndexDf
  summitDistDf <- masterIndexDf
  summitPosDf <- masterIndexDf
  targetOvlpDf <- masterIndexDf
  peakOvlpDf <- masterIndexDf

  ## create dataframes with columns for each peak category
  ## each row is a peak and elements are stores as list of values of interest
  for (ctg in c("featureInPeak", "nearStart", "nearEnd", "peakInFeature", "upstreamTss")) {
    groupedDf <- dplyr::filter(subData, peakCategory == ctg) %>%
      dplyr::select(name, rowIdx, target, peakDist, summitDist,
                    relativePeakPos, targetOverlap, peakOverlap) %>%
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

    ## create summit distance df for each peak category
    summitDf <- dplyr::summarise(groupedDf, !!ctg := list(summitDist))
    summitDistDf <- dplyr::left_join(x = summitDistDf, y = summitDf, by = "name")

    ## create peak position df for each peak category
    posDf <- dplyr::summarise(groupedDf, !!ctg := list(relativePeakPos))
    summitPosDf <- dplyr::left_join(x = summitPosDf, y = posDf, by = "name")

    ## create target overlap df for each cateogry
    targetCov <- dplyr::summarise(groupedDf, !!ctg := list(targetOverlap))
    targetOvlpDf <- dplyr::left_join(x = targetOvlpDf, y = targetCov, by = "name")

    ## create peak overlap df for each category
    peakCov <- dplyr::summarise(groupedDf, !!ctg := list(peakOverlap))
    peakOvlpDf <- dplyr::left_join(x = peakOvlpDf, y = peakCov, by = "name")
  }

  ## Future development: use priority list of peak category to select peak targets
  markPseudoIdx <- c()
  bidirectBIdx <- c()
  bidirectCIdx <- c()

  ## "featureInPeak", "nearStart", "nearEnd", "peakInFeature", "upstreamTss"
  ## nearStart peak
  ###########
  ## A1) nearStart & upstreamTss:
  ruleA1 <- which(peakFound$nearStart & peakFound$upstreamTss)
  # targetDf[unlist(masterIndexDf$nearStart[ruleA1]), ]

  ruleA1_bidirectIndex <- tibble::tibble(
    masterIndex = ruleA1,
    nearStart = masterIndexDf$nearStart[ruleA1],
    upstreamTss = masterIndexDf$upstreamTss[ruleA1]
  ) %>%
    tidyr::unnest(cols = c(nearStart, upstreamTss)) %>%
    dplyr::mutate(
      nearStartStrand = targetDf$targetStrand[nearStart],
      upstreamTssStrand = targetDf$targetStrand[upstreamTss]
    ) %>%
    dplyr::mutate(
      strandMatch = if_else(condition = nearStartStrand == upstreamTssStrand,
                            true = "same", false = "opposite")
    )

  ruleA1_sameDirUp <- ruleA1_bidirectIndex$masterIndex[
    which(ruleA1_bidirectIndex$strandMatch == "same")]

  ruleA1_bidirectPseudo <- nearest_upstream_bidirectional(
    targetDf = targetDf,
    t1Idx = ruleA1_bidirectIndex$nearStart,
    t2Idx = ruleA1_bidirectIndex$upstreamTss,
    promoterLength = promoterLength,
    upstreamLimit = upstreamLimit,
    bidirectionalDistance = bidirectionalDistance,
    pointBasedAnnotation = pointBasedAnnotation
  )

  ruleA1_bidirect <- purrr::map2(
    .x = ruleA1_bidirectIndex$nearStart,
    .y = ruleA1_bidirectIndex$upstreamTss,
    .f = function(x, y){
      if(any(c(x, y) %in% ruleA1_bidirectPseudo)){
        return(NULL)
      } else{
        return(c(x, y))
      }
    }) %>%
    purrr::flatten_int()

  ## ACTION: mark upstreamTss as pseudo based on nearest_upstream_bidirectional()
  ## no need to use masterIndexDf as ruleA1_bidirectPseudo has indices from targetDf
  markPseudoIdx <- append(markPseudoIdx, ruleA1_bidirectPseudo)
  bidirectBIdx <- union(bidirectBIdx, ruleA1_bidirect)

  ruleA1_farUp <- ruleA1[purrr::map_lgl(.x = peakDistDf$upstreamTss[ruleA1],
                                        .f = function(x){abs(x[1]) > upstreamLimit})]

  ## additionally add targets which are in same direction to mark them as NULL
  ## use union to remove duplicates
  ruleA1_farUp <- union(x = ruleA1_farUp, ruleA1_sameDirUp)

  ## ACTION: set the upstreamTss to NULL if it is far than promoterLength
  masterIndexDf$upstreamTss[ruleA1_farUp] <- purrr::map(
    .x = masterIndexDf$upstreamTss[ruleA1_farUp],
    .f = ~ NULL)

  peakFound$upstreamTss[ruleA1_farUp] <- FALSE


  ###########
  ## A2) nearStart & nearEnd:
  ruleA2 <- which(peakFound$nearStart & peakFound$nearEnd)
  # targetDf[unlist(masterIndexDf$nearStart[ruleA2]), ]

  ## ACTION: set TES target to NULL
  masterIndexDf$nearEnd[ruleA2] <- purrr::map(.x = masterIndexDf$nearEnd[ruleA2],
                                              .f = ~ NULL)
  peakFound$nearEnd[ruleA2] <- FALSE


  ###########
  ## A3) nearStart & peakInFeature:
  ruleA3 <- which(peakFound$nearStart & peakFound$peakInFeature)
  # targetDf[unlist(masterIndexDf$nearStart[ruleA3]), ]
  ## ACTION: nothing

  ###########
  ## A4) nearStart & featureInPeak:
  ruleA4 <- which(peakFound$nearStart & peakFound$featureInPeak)
  # targetDf[unlist(masterIndexDf$nearStart[ruleA4]), ]
  ## ACTION: nothing

  ## featureInPeak peak
  ###########
  ## B1) featureInPeak & nearEnd:
  ruleB1 <- which(peakFound$featureInPeak & peakFound$nearEnd)
  # targetDf[unlist(masterIndexDf$featureInPeak[ruleB1]), ]

  ## ACTION: remove nearEnd target
  masterIndexDf$nearEnd[ruleB1] <- purrr::map(.x = masterIndexDf$nearEnd[ruleB1],
                                              .f = ~ NULL)
  peakFound$nearEnd[ruleB1] <- FALSE


  ###########
  ## B2) featureInPeak & upstreamTss:
  ruleB2 <- which(peakFound$featureInPeak & peakFound$upstreamTss)
  # targetDf[unlist(masterIndexDf$featureInPeak[ruleB2]), ]

  ## ********************* :future development: *******************************
  ## work only on masterIndexDf and peakFound. use index directly to comapre/extract
  ## various values from targetDf
  # ruleB2_index <- masterIndexDf[ruleB2, ] %>%
  #   dplyr::select(name, featureInPeak, upstreamTss) %>%
  #   tidyr::unnest(cols = c(featureInPeak, upstreamTss))
  ## ********************* :future development: *******************************

  ruleB2_farUp <- ruleB2[purrr::map_lgl(.x = peakDistDf$upstreamTss[ruleB2],
                                        .f = function(x){abs(x[1]) > upstreamLimit})]

  ruleB2_upLimit <- setdiff(ruleB2, ruleB2_farUp)

  ## ACTION: set the upstreamTss to NULL if it is far than promoterLength
  masterIndexDf$upstreamTss[ruleB2_farUp] <- purrr::map(
    .x = masterIndexDf$upstreamTss[ruleB2_farUp],
    .f = ~ NULL)

  peakFound$upstreamTss[ruleB2_farUp] <- FALSE

  ## ACTION: for the TF which has known binding over gene body (E.g. polII ChIP or histone marks)
  ## preference is for featureInPeak. all other targets are pseudo
  if(bindingInGene){
    masterIndexDf$upstreamTss[ruleB2] <- purrr::map(.x = masterIndexDf$upstreamTss[ruleB2],
                                                    .f = ~ NULL)
    peakFound$upstreamTss[ruleB2] <- FALSE

  } else{
    ## set bidirectional = 2
    bidirectCIdx <- c(unlist(masterIndexDf$upstreamTss[ruleB2_upLimit]),
                      unlist(masterIndexDf$featureInPeak[ruleB2_upLimit]))

  }


  ###########
  ## B3) featureInPeak & peakInFeature:
  ruleB3 <- which(peakFound$featureInPeak & peakFound$peakInFeature)
  # targetDf[unlist(masterIndexDf$featureInPeak[ruleB3]), ]

  if(bindingInGene){
    ## for the TF which has known binding over gene body (E.g. polII ChIP or histone marks)
    ## ACTION: set peakInFeature to NULL
    masterIndexDf$peakInFeature[ruleB3] <- purrr::map(.x = masterIndexDf$peakInFeature[ruleB3],
                                                      .f = ~ NULL)
    peakFound$peakInFeature[ruleB3] <- FALSE
  } else{
    ## ACTION: set featureInPeak to NULL
    masterIndexDf$featureInPeak[ruleB3] <- purrr::map(.x = masterIndexDf$featureInPeak[ruleB3],
                                                      .f = ~ NULL)
    peakFound$featureInPeak[ruleB3] <- FALSE
  }

  ###########
  ## C1) upstreamTss & nearEnd:
  ruleC1 <- which(peakFound$upstreamTss & peakFound$nearEnd)
  # targetDf[unlist(masterIndexDf$upstreamTss[ruleC1]), ]

  ruleC1_pro <- ruleC1[purrr::map_lgl(.x = peakDistDf$upstreamTss[ruleC1],
                                      .f = function(x){abs(x[1]) <= promoterLength})]

  ## peaks which also have upstreamTss but its far than upstreamLimit
  ruleC1_farUp <- ruleC1[purrr::map_lgl(.x = peakDistDf$upstreamTss[ruleC1],
                                        .f = function(x){abs(x[1]) > upstreamLimit})]

  ## peaks which upstream but within upstreamLimit range
  ruleC1_upLimit <- setdiff(x = ruleC1, y = ruleC1_farUp)
  ruleC1_nearUp <- setdiff(x = ruleC1, y = c(ruleC1_pro, ruleC1_farUp))

  ## among these ruleC1_nearUp peaks, check: nearEnd has peak summit after insideSkewToEndCut OR
  ## the peak overlap with nearEnd target is very small. i.e. find nearEnd targets where peaks
  ## are at very end of the target. upstreamTss targets should be TRUE for these
  ruleC1_veryEnd <- purrr::map2_lgl(
    .x = summitPosDf$nearEnd[ruleC1_nearUp],
    .y = peakOvlpDf$nearEnd[ruleC1_nearUp],
    .f = function(x, y){x[1] >= insideSkewToEndCut || y[1] < 0.1})

  ## #********************************************************************************##
  ## using 0.1 as cutoff for peakOverlap. later can be converted to function argument ##
  ## #********************************************************************************##

  ## ACTION: set upstreamTss = NULL if: upstreamTss is far & nearEnd is not at very end
  masterIndexDf$upstreamTss[ruleC1_farUp] <- purrr::map(
    .x = masterIndexDf$upstreamTss[ruleC1_farUp],
    .f = ~ NULL)

  peakFound$upstreamTss[ruleC1_farUp] <- FALSE

  ## add ruleC1_nearUp peaks which are also ruleC1_veryEnd to ruleC1_pro
  ruleC1_pro <- append(ruleC1_pro, ruleC1_nearUp[ruleC1_veryEnd])

  ## ACTION: set the nearEnd peaks to pseudo because:
  ## it is upstreamTss within promoter OR (nearUp range and it is at very end for nearEnd)
  markPseudoIdx <- append(markPseudoIdx, unlist(masterIndexDf$nearEnd[ruleC1_pro]))


  ###########
  ## C2) upstreamTss & peakInFeature:
  ruleC2 <- which(peakFound$upstreamTss & peakFound$peakInFeature)
  # targetDf[unlist(masterIndexDf$upstreamTss[ruleC2]), ]

  ruleC2_pro <- ruleC2[purrr::map_lgl(.x = peakDistDf$upstreamTss[ruleC2],
                                      .f = function(x){abs(x[1]) <= promoterLength})]

  ruleC2_farUp <- ruleC2[purrr::map_lgl(.x = peakDistDf$upstreamTss[ruleC2],
                                        .f = function(x){abs(x[1]) > upstreamLimit})]

  ## accept all upstreamTss < promoterLength i.e. ruleC2_pro and mark bidirectional
  ## reject all upstreamTss > upstreamLimit i.e. ruleC2_farUp

  ## promoterLength < upstreamTss < upstreamLimit: decide using nearest_upstream_bidirectional()
  ruleC2_nearUp <- setdiff(x = ruleC2, y = c(ruleC2_pro, ruleC2_farUp))


  ## identify ruleC2_nearUp peaks which are near start for peakInFeature
  ruleC2_nearUp_ovStart <- ruleC2_nearUp[purrr::map_lgl(
    .x = summitPosDf$peakInFeature[ruleC2_nearUp],
    .f = function(x){all(x < (1 - insideSkewToEndCut))})]

  ## solve this using nearest_upstream_bidirectional()
  ruleC2_bidirectIndex <- tibble::tibble(
    masterIndex = ruleC2_nearUp_ovStart,
    peakInFeature = masterIndexDf$peakInFeature[ruleC2_nearUp_ovStart],
    upstreamTss = masterIndexDf$upstreamTss[ruleC2_nearUp_ovStart]
  ) %>%
    tidyr::unnest(cols = c(peakInFeature, upstreamTss)) %>%
    dplyr::mutate(
      peakInFeatureStrand = targetDf$targetStrand[peakInFeature],
      upstreamTssStrand = targetDf$targetStrand[upstreamTss]
    ) %>%
    dplyr::mutate(
      strandMatch = if_else(condition = peakInFeatureStrand == upstreamTssStrand,
                            true = "same", false = "opposite")
    )

  # ruleC2_targetDf <- targetDf
  # ruleC2_targetDf$targetStart[ruleC2_bidirectIndex$peakInFeature] <- dplyr::if_else(
  #   condition = ruleC2_targetDf$targetStrand[ruleC2_bidirectIndex$peakInFeature] == "+",
  #   true = ruleC2_targetDf$peakSummit[ruleC2_bidirectIndex$peakInFeature],
  #   false = ruleC2_targetDf$targetStart[ruleC2_bidirectIndex$peakInFeature])
  #
  # ruleC2_targetDf$targetEnd[ruleC2_bidirectIndex$peakInFeature] <- dplyr::if_else(
  #   condition = ruleC2_targetDf$targetStrand[ruleC2_bidirectIndex$peakInFeature] == "-",
  #   true = ruleC2_targetDf$peakSummit[ruleC2_bidirectIndex$peakInFeature],
  #   false = ruleC2_targetDf$targetEnd[ruleC2_bidirectIndex$peakInFeature])

  ruleC2_bidirectPseudo <- nearest_upstream_bidirectional(
    targetDf = targetDf,
    t1Idx = ruleC2_bidirectIndex$peakInFeature,
    t2Idx = ruleC2_bidirectIndex$upstreamTss,
    promoterLength = promoterLength,
    upstreamLimit = upstreamLimit,
    bidirectionalDistance = bidirectionalDistance,
    pointBasedAnnotation = pointBasedAnnotation
  )

  ruleC2_bidirect <- purrr::map2(
    .x = ruleC2_bidirectIndex$peakInFeature,
    .y = ruleC2_bidirectIndex$upstreamTss,
    .f = function(x, y){
      if(any(c(x, y) %in% ruleC2_bidirectPseudo)){
        return(NULL)
      } else{
        return(c(x, y))
      }
    }) %>%
    purrr::flatten_int()

  ## no need to use masterIndexDf as ruleC2_bidirectPseudo has indices from targetDf
  markPseudoIdx <- append(markPseudoIdx, ruleC2_bidirectPseudo)
  bidirectBIdx <- union(bidirectBIdx, unlist(masterIndexDf$peakInFeature[ruleC2_pro]))
  bidirectBIdx <- union(bidirectBIdx, unlist(masterIndexDf$upstreamTss[ruleC2_pro]))
  bidirectBIdx <- union(bidirectBIdx, ruleC2_bidirect)

  ## remaining ruleC2_nearUp: set as NULL
  ruleC2_nearUp_remaining <- setdiff(x = ruleC2_nearUp, y = ruleC2_nearUp_ovStart)
  ruleC2_farUp <- append(ruleC2_farUp, ruleC2_nearUp_remaining)

  ## ACTION: set upstreamTss = NULL if: upstreamTss is beyond upstreamLimit
  masterIndexDf$upstreamTss[ruleC2_farUp] <- purrr::map(
    .x = masterIndexDf$upstreamTss[ruleC2_farUp],
    .f = ~ NULL)

  peakFound$upstreamTss[ruleC2_farUp] <- FALSE


  ###########
  ## D1) peakInFeature & nearEnd:
  ruleD1 <- which(peakFound$peakInFeature & peakFound$nearEnd)
  # targetDf[unlist(masterIndexDf$peakInFeature[ruleD1]), ]
  ## ACTION: nothing

  ###########
  ## for upstreamTss peak
  ## 8) two upstreamTss peaks: remove pseudo_upstream if peakDist > promoterLength
  ## ACTION: already taken care of in upstream_annotations()


  ## mark selected targets to pseudo
  targetDf$peakAnnotation[markPseudoIdx] <- paste("pseudo_", targetDf$peakAnnotation[markPseudoIdx], sep = "")

  ## add bidirectional information
  targetDf$bidirectional[bidirectBIdx] <- "B"
  targetDf$bidirectional[bidirectCIdx] <- "C"

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
