
#' Get good replicate's peakset from two or more ChIPseq replicates
#'
#' This function internally uses \code{chipmine::combinatorial_binding_matrix}
#' function to generate a peak overlap matrix and selects the peaks which are
#' present in all the replicates. Then it uses the best replicate to extract the
#' target gene infomation.
#'
#' @param sampleInfo Sample information dataframe
#' @param cdsFile BED file with gene information
#' @param ... Other arguments for \code{preProcess_macs2_results()} function
#'
#' @return A dataframe with target gene information for confident peaks
#' @export
#'
#' @examples
best_replicate_peakset <- function(sampleInfo, cdsFile, ...){

  ## generate the combinatorial binding matrix using peak overlap
  mat <- combinatorial_binding_matrix(sampleInfo = sampleInfo)

  bestRepIndex <- which(sampleInfo$bestRep == 1)
  bestRep <- sampleInfo$sampleId[bestRepIndex]
  bestRepPeakIdCol <- paste("peakId.", bestRep, sep = "")

  # dplyr::group_by_at(mat, vars(starts_with("overlap."))) %>%
  #   dplyr::summarise(n = n())

  ## identify common peaks between replicates
  commonPeaks <- dplyr::filter_at(.tbl = mat, .vars = vars(starts_with("overlap.")),
                                  .vars_predicate = all_vars(. == TRUE)) %>%
    dplyr::select(ends_with(bestRep)) %>%
    dplyr::select(starts_with("peakId.")) %>%
    dplyr::distinct()

  ## get the target gene information for best ChIPseq replicate
  peakAnn <- import_peak_annotation(
    sampleId = sampleInfo$sampleId[ bestRepIndex ],
    peakAnnoFile = sampleInfo$narrowpeakAnno[ bestRepIndex ],
    peakFile = sampleInfo$narrowpeakFile[ bestRepIndex ],
    bwFile = sampleInfo$bwFile[ bestRepIndex ],
    fcCutoff = 1, pvalCutoff = 1)

  bestRepAnn <- dplyr::left_join(x = commonPeaks, y = peakAnn,
                                 by = structure(bestRepPeakIdCol, names = bestRepPeakIdCol))

  bestRepTargets <- preProcess_macs2_results(
    sampleId = sampleInfo$sampleId[ bestRepIndex ],
    peakAnnotation = bestRepAnn,
    cdsFile = cdsFile,
    peakFile = sampleInfo$narrowpeakFile[ bestRepIndex ],
    bwFile = sampleInfo$bwFile[ bestRepIndex ],
    ...) %>%
    dplyr::select(gene, starts_with("hasPeak"), starts_with("peakId"), starts_with("peakPosition"),
                  starts_with("peakType"), starts_with("peakDist"), starts_with("peakPval"),
                  starts_with("peakEnrichment"), starts_with("relativeSummitPos"), starts_with("bidirectional"))

  # "hasPeak", "peakPosition", "peakType", "peakId", "peakEnrichment", "peakPval", "peakQval",
  # "peakSummit", "peakDist", "summitDist", "bidirectional", "featureCovFrac", "relativeSummitPos",
  # "peakRegion", "peakCoverage"

  ## generate confident peak and target gene list data
  confidentPeaks <- dplyr::left_join(x = commonPeaks, y = bestRepTargets,
                                     by = structure(bestRepPeakIdCol, names = bestRepPeakIdCol)) %>%
    dplyr::select(gene, starts_with("hasPeak."), everything())

  return(confidentPeaks)
}



