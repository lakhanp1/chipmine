
#' Generate the combinatorial peak association matrix with merged master peak list
#'
#' This method reads all the narrowPeak files for given TFs and create a merged peak region list.
#' All the original peaks are then searched for overlap against this master list and a combinatorial
#' dataframe showing the presence of peak in the master peak list is returned. If a region has
#' more than one peak overlapping, best peak is selected using \code{peakPval} column.
#' This is done to avoid the exponential growth in number of rows with increasing samples.
#' Additionally, it also extracts the sequence around summit position
#'
#' @param sampleInfo Sample information dataframe
#' @param peakRegions Optional GRanges object which has master peak regions. If not provided,
#' peaks from narrowPeak files are merged to create new master peakset.
#' @param peakFormat Format of the peak file. One of \code{"narrowPeak", "broadPeak", "bed"}
#' @param peakCols Column to extract from peak file. Default: \code{c("peakId", "peakEnrichment", "peakPval")}
#' @param genome Optionally BSgenome object for extracting summit sequence
#' @param summitSeqLen Length of sequence to extract at summit position. Default: 200
#'
#' @return A dataframe with a masterDf of peak regions generated after merging all peak regions from all samples. For each sample, its association with regions in the masterDf is reported.
#' @export
#'
#' @examples NA
combinatorial_binding_matrix <- function(sampleInfo, peakRegions = NULL, peakFormat = "narrowPeak",
                                         peakCols = c("peakId", "peakEnrichment", "peakPval"),
                                         genome = NULL, summitSeqLen = 200){

  peakList <- GenomicRanges::GRangesList(lapply(X = sampleInfo$narrowpeakFile,
                                                FUN = rtracklayer::import, format = peakFormat))

  names(peakList) <- sampleInfo$sampleId


  ## create a GRangesList masterlist of peak regions
  if(is.null(peakRegions)){
    peakRegions <- GenomicRanges::reduce(unlist(peakList, recursive = TRUE, use.names = T))
  }

  ## region name column
  if(is.null(mcols(peakRegions)$name)){
    mcols(peakRegions)$name <- paste("peak_region", 1:length(peakRegions), sep = "_")
  }

  ## master list as dataframe
  masterDf <- as.data.frame(peakRegions) %>%
    dplyr::select(-width, -strand) %>%
    tibble::as_tibble()

  masterDf <- dplyr::mutate_if(.tbl = masterDf, .predicate = is.factor, .funs = as.character)

  ## find the overlap of individual peaklist with the masterDf
  for (i in 1:nrow(sampleInfo)) {

    sampleName <- sampleInfo$sampleId[i]
    peakIdCol <- paste("peakId.", sampleName, sep = "")
    overlapPeakCol <- paste("overlap.", sampleName, sep = "")

    ## findOverlaps
    hits <- as.data.frame(GenomicRanges::findOverlaps(query = peakRegions, subject = peakList[[sampleName]]))
    hits$peakRegionName <- peakRegions$name[hits$queryHits]
    hits[[peakIdCol]] <- peakList[[sampleName]]$name[hits$subjectHits]

    hits$queryHits <- NULL
    hits$subjectHits <- NULL

    dt <- import_peaks_as_df(file = sampleInfo$narrowpeakFile[i],
                             sampleId = sampleName,
                             peakCols = peakCols) %>%
      dplyr::distinct()


    dt[[overlapPeakCol]] <- TRUE

    ## select best peak within region
    hits <- dplyr::left_join(x = hits, y = dt, by = setNames(peakIdCol, peakIdCol)) %>%
      dplyr::select("peakRegionName", !!peakIdCol, !! overlapPeakCol, dplyr::everything()) %>%
      dplyr::group_by(peakRegionName) %>%
      dplyr::arrange_at(.vars = vars(starts_with("peakPval.")),
                        .funs = list(~dplyr::desc(.)), .by_group = TRUE) %>%
      dplyr::slice(1L) %>%
      dplyr::ungroup()

    # hits <- dplyr::group_by(hits, peakRegionName) %>%
    #   dplyr::summarise(!!sampleName := n()) %>%
    #   dplyr::ungroup()

    if(! is.null(genome)){
      ## get the sequence around summit
      summitSeq <- get_narrowpeak_summit_seq(npFile = sampleInfo$narrowpeakFile[i],
                                             id = sampleName,
                                             genome = genome,
                                             length = summitSeqLen)

      hits <- dplyr::left_join(x = hits, y = summitSeq, by = setNames("name", peakIdCol))
    }

    ## join with master data
    # masterDf <- dplyr::left_join(x = masterDf, y = hits, by = c("name" = "peakRegionName"))
    masterDf <- dplyr::left_join(x = masterDf, y = hits, by = c("name" = "peakRegionName")) %>%
      tidyr::replace_na(replace = purrr::set_names(list(FALSE), nm = overlapPeakCol))

  }


  return(masterDf)

}



