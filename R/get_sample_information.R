
##################################################################################
## generate experiment data
#' Create sample information dataframe
#'
#' @param exptInfoFile A tabular file with details for each sample. "sampleId" and "IP_tag" columns are must
#' @param dataPath Path where the data is stored
#' @param matrixSource Source of profile matrix. One of "deeptools", "miao", "normalizedmatrix",
#' "normalizedmatrix_5kb", "TSS_4kb_2kb_normalized", "TES_2kb_4kb_normalized", "TSS_3kb_3kb_normalized",
#' "TES_3kb_3kb_normalized". Default: deeptools
#' @param samples A vector of Sample IDs which are to be processed
#' @param profileType Type of profile. This will be added as suffix to the profile name
#'
#' @return Data frame with details for each sample
#' @export
#'
#' @examples NA
get_sample_information = function(exptInfoFile, samples, dataPath, matrixSource = "deeptools",
                                  profileType = "profile"){

  ## read the experiment sample details and select only those which are to be plotted
  exptData = data.table::fread(input = exptInfoFile,
                               sep = "\t",
                               stringsAsFactors = F,
                               header = T,
                               data.table = F) %>%
    dplyr::filter(sampleId %in% samples)

  exptData$sampleId = factor(exptData$sampleId, levels = samples)

  exptData = exptData[order(exptData$sampleId), ] %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::mutate(
      profileName = paste(sampleId, profileType, sep = "_"),
      bwFile = paste(dataPath, "/", sampleId, "/", sampleId, "_normalized.bw", sep = ""),
      matFile = case_when(
        matrixSource == "deeptools" ~
          paste(dataPath, "/", sampleId, "/", sampleId, "_normalized_profile.tab.gz", sep = ""),
        matrixSource == "normalizedmatrix" ~
          paste(dataPath, "/", sampleId, "/", sampleId, "_normalizedMatrix.tab.gz", sep = ""),
        matrixSource == "normalizedmatrix_5kb" ~
          paste(dataPath, "/", sampleId, "/", sampleId, "_normalizedMatrix_5kb.tab.gz", sep = ""),
        matrixSource == "TSS_4kb_2kb_normalized" ~
          paste(dataPath, "/", sampleId, "/", sampleId, "_normalizedMatrix_4kbTSS2kb.tab.gz", sep = ""),
        matrixSource == "TES_2kb_4kb_normalized" ~
          paste(dataPath, "/", sampleId, "/", sampleId, "_normalizedMatrix_2kbTES4kb.tab.gz", sep = ""),
        matrixSource == "TSS_3kb_3kb_normalized" ~
          paste(dataPath, "/", sampleId, "/", sampleId, "_normalizedMatrix_3kbTSS3kb.tab.gz", sep = ""),
        matrixSource == "TES_3kb_3kb_normalized" ~
          paste(dataPath, "/", sampleId, "/", sampleId, "_normalizedMatrix_3kbTES3kb.tab.gz", sep = "")
      ),
      polIIExpFile = dplyr::if_else(
        IP_tag == "polII",
        paste(dataPath, "/", sampleId, "/", sampleId, "_normalizedExpression.tab", sep = ""),
        "NA"
      ),
      polIIExpMat = dplyr::if_else(
        IP_tag == "polII",
        paste(dataPath, "/", sampleId, "/", sampleId, "_polii_expr.tab.rel.mat", sep = ""),
        "NA"
      ),
      clusterFile = dplyr::if_else(
        IP_tag == "polII",
        "NA",
        paste(dataPath, "/", sampleId, "/", sampleId, ".kmeans.clusters.txt", sep = "")
      ),
      tfPeakFile = dplyr::if_else(
        IP_tag == "polII",
        "NA",
        paste(dataPath, "/", sampleId, "/", sampleId, "_peaks.annotated.tab", sep = "")
      ),
      mergedDataFile = dplyr::if_else(
        IP_tag == "polII",
        "NA",
        paste(dataPath, "/", sampleId, "/", sampleId, "_allGenes_clusters.tab", sep = "")
      ),
      tfRegionsBed = dplyr::if_else(
        IP_tag == "polII",
        "NA",
        paste(dataPath, "/", sampleId, "/", sampleId, "_peakRegions.bed", sep = "")
      ),
      tfExpAtUpstreamFile = dplyr::if_else(
        IP_tag == "polII",
        "NA",
        paste(dataPath, "/", sampleId, "/", sampleId, "_peakExp_cdsUpstream.rel.mat", sep = "")
      ),
      tfExpAtPeakFile = dplyr::if_else(
        IP_tag == "polII",
        "NA",
        paste(dataPath, "/", sampleId, "/", sampleId, "_narrowpeakExp_peakBed.rel.mat", sep = "")
      ),
      narrowpeakFile = dplyr::if_else(
        IP_tag == "polII",
        "NA",
        paste(dataPath, "/", sampleId, "/", sampleId, "_narrow_peaks.narrowPeak", sep = "")
      ),
      narrowpeakAnno = dplyr::if_else(
        IP_tag == "polII",
        "NA",
        paste(dataPath, "/", sampleId, "/", sampleId, ".narrowPeak.nearestCDS.tab", sep = "")
      ),
      broadpeakFile = dplyr::if_else(
        IP_tag == "polII",
        "NA",
        paste(dataPath, "/", sampleId, "/", sampleId, "_broad_peaks.broadPeak", sep = "")
      ),
      broadpeakAnno = dplyr::if_else(
        IP_tag == "polII",
        "NA",
        paste(dataPath, "/", sampleId, "/", sampleId, ".broadPeak.nearestCDS.tab", sep = "")
      )
    ) %>%
    dplyr::mutate_at(c("polIIExpFile", "polIIExpMat", "clusterFile", "tfPeakFile", "mergedDataFile", "narrowpeakFile", "narrowpeakAnno", "broadpeakFile", "broadpeakAnno"),
                     .funs = dplyr::funs(dplyr::na_if(., "NA")))


  return(exptData)
}

##################################################################################

