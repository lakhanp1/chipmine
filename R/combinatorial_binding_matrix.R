
#' Generate the combinatorial peak association matrix with merged master peak list
#'
#' This method reads all the narrowPeak files for given TFs and create a merged peak region list.
#' All the original peaks are then searched for overlap against this master list and a combinatorial
#' dataframe showing the presence of peak in the master peak list is returned.
#' Additionally, it also extracts the sequence around summit position
#'
#' @param sampleInfo Sample information dataframe
#' @param genome BSgenome object for extracting summit sequence
#' @param summitSeqLen Length of sequence to extract at summit position. Default: 200
#'
#' @return A dataframe with a masterlist of peak regions generated after merging all peak regions from all samples. For each sample, its association with regions in the masterlist is reported.
#' @export
#'
#' @examples NA
combinatorial_binding_matrix = function(sampleInfo, genome, summitSeqLen = 200){

  ## read the narrowpeak files in the GRanges object
  narrowpeakCols = c(signalValue = "numeric", pValue = "numeric",
                     qValue = "numeric", peak = "integer")

  peakList = GenomicRanges::GRangesList(lapply(X = sampleInfo$narrowpeakFile,
                                               FUN = rtracklayer::import, format = "bed", extraCols = narrowpeakCols))

  names(peakList) = sampleInfo$sampleId


  ## create a GRangesList masterlist of peak regions
  peakRegions = GenomicRanges::reduce(unlist(peakList, recursive = TRUE, use.names = T))
  mcols(peakRegions) = data.frame(name = paste("peak_region", 1:length(peakRegions), sep = "_"), stringsAsFactors = F)

  ## master list as dataframe
  masterList = data.frame(
    chr = seqnames(peakRegions),
    start = start(peakRegions),
    end = end(peakRegions),
    name = peakRegions$name,
    stringsAsFactors = F
  )

  # np = names(peakList)[1]

  ## find the overlap of individual peaklist with the masterlist
  # for (np in names(peakList)) {
  for (i in 1:nrow(sampleInfo)) {

    np = sampleInfo$sampleId[i]

    cat("Reading peak information for sample: ", np, "\n")

    ## findOverlaps
    hits = as.data.frame(GenomicRanges::findOverlaps(query = peakRegions, subject = peakList[[np]]))
    hits$peakRegionName = peakRegions$name[hits$queryHits]
    hits[[np]] = peakList[[np]]$name[hits$subjectHits]

    hits$queryHits = NULL
    hits$subjectHits = NULL


    peakIdCol = paste("peakId.", np, sep = "")
    overlapPeakCol = paste("overlap.", np, sep = "")

    dt = read_peak_annotation_file(title = np,
                                   file = sampleInfo$narrowpeakAnno[i],
                                   cols = c("peakId", "peakEnrichment", "peakPval", "peakSummit"))

    dt[[overlapPeakCol]] = TRUE

    ## get the sequence around summit
    summitSeq = get_narrowpeak_summit_seq(npFile = sampleInfo$narrowpeakFile[i],
                                          id = np,
                                          genome = genome,
                                          length = summitSeqLen)


    hits = dplyr::left_join(x = hits, y = dt, by = setNames(peakIdCol, np)) %>%
      dplyr::left_join(y = summitSeq, by = setNames("name", np))


    masterList = dplyr::left_join(x = masterList, y = hits, by = c("name" = "peakRegionName"))

    cat("Done...\n")

  }

  return(masterList)

}



