
#' Generate the combinatorial peak association matrix with merged master peak list
#'
#' This method reads all the narrowPeak files for given TFs and create a merged peak region list.
#' All the original peaks are then searched for overlap against this master list and a combinatorial
#' dataframe showing the presence of peak in the master peak list is returned.
#' Additionally, it also extracts the sequence around summit position
#'
#' @param sampleInfo Sample information dataframe
#' @param peakFormat Format of the peak file. One of \code{"narrowPeak", "broadPeak", "bed"}
#' @param peakCols Column to extract from peak file. Default: \code{c("peakId", "peakEnrichment", "peakPval")}
#' @param genome Optionally BSgenome object for extracting summit sequence
#' @param summitSeqLen Length of sequence to extract at summit position. Default: 200
#'
#' @return A dataframe with a masterlist of peak regions generated after merging all peak regions from all samples. For each sample, its association with regions in the masterlist is reported.
#' @export
#'
#' @examples NA
combinatorial_binding_matrix <- function(sampleInfo, peakFormat = "narrowPeak",
                                         peakCols = c("peakId", "peakEnrichment", "peakPval"),
                                         genome = NULL, summitSeqLen = 200){

  peakList <- GenomicRanges::GRangesList(lapply(X = sampleInfo$narrowpeakFile,
                                                FUN = rtracklayer::import, format = peakFormat))

  names(peakList) <- sampleInfo$sampleId


  ## create a GRangesList masterlist of peak regions
  peakRegions <- GenomicRanges::reduce(unlist(peakList, recursive = TRUE, use.names = T))
  mcols(peakRegions) = data.frame(name = paste("peak_region", 1:length(peakRegions), sep = "_"), stringsAsFactors = F)

  ## master list as dataframe
  masterList <- data.frame(
    chr = seqnames(peakRegions),
    start = start(peakRegions),
    end = end(peakRegions),
    name = peakRegions$name,
    stringsAsFactors = F
  )


  ## find the overlap of individual peaklist with the masterlist
  for (i in 1:nrow(sampleInfo)) {

    sampleName <- sampleInfo$sampleId[i]
    peakIdCol <- paste("peakId.", sampleName, sep = "")
    overlapPeakCol <- paste("overlap.", sampleName, sep = "")

    cat("Reading peak information for sample: ", sampleName, "\n")

    ## findOverlaps
    hits <- as.data.frame(GenomicRanges::findOverlaps(query = peakRegions, subject = peakList[[sampleName]]))
    hits$peakRegionName <- peakRegions$name[hits$queryHits]
    hits[[peakIdCol]] <- peakList[[sampleName]]$name[hits$subjectHits]

    hits$queryHits <- NULL
    hits$subjectHits <- NULL

    ## get the additional columns from peak file
    dt <- read_peak_annotation_file(title = sampleName,
                                    file = sampleInfo$narrowpeakAnno[i],
                                    cols = peakCols)

    dt[[overlapPeakCol]] <- TRUE

    hits <- dplyr::left_join(x = hits, y = dt, by = setNames(peakIdCol, peakIdCol)) %>%
      dplyr::select("peakRegionName", !!peakIdCol, !! overlapPeakCol, dplyr::everything())

    if(! is.null(genome)){
      ## get the sequence around summit
      summitSeq <- get_narrowpeak_summit_seq(npFile = sampleInfo$narrowpeakFile[i],
                                             id = sampleName,
                                             genome = genome,
                                             length = summitSeqLen)

      hits <- dplyr::left_join(x = hits, y = summitSeq, by = setNames("name", peakIdCol))
    }

    ## join with master data
    masterList <- dplyr::left_join(x = masterList, y = hits, by = c("name" = "peakRegionName"))

    cat("Done...\n")

  }

  return(masterList)

}



