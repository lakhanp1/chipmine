
##################################################################################
## generate experiment data
#' Create sample information dataframe
#'
#' @param exptInfoFile A tabular file with details for each sample. "sampleId" and "IP_tag" columns are must
#' @param dataPath Path where the data is stored
#' @param matrixSource Source of profile matrix. One of "deeptools", "miao", "normalizedmatrix",
#' "normalizedmatrix_5kb", "TSS_4kb_2kb_normalized", "TES_2kb_4kb_normalized", "TSS_3kb_3kb_normalized",
#' "TES_3kb_3kb_normalized". Default: normalizedmatrix
#' @param samples A vector of Sample IDs which are to be processed. Default: All samples are used
#' @param profileType Type of profile. This will be added as suffix to the profile name
#' @param macs2Control Whether control was used for macs2 peak calling. Default: TRUE
#'
#' @return Data frame with details for each sample
#' @export
#'
#' @examples NA
get_sample_information <- function(exptInfoFile, samples = NULL, dataPath, matrixSource = "normalizedmatrix",
                                   profileType = "profile", macs2Control = TRUE){

  ## read the experiment sample details and select only those which are to be plotted
  exptData <- suppressMessages(readr::read_tsv(file = exptInfoFile))

  if(!is.null(samples)){
    exptData <- dplyr::filter(exptData, sampleId %in% samples)
  }

  exptData$sampleId <- factor(exptData$sampleId, levels = unique(exptData$sampleId))

  exptData <- exptData[order(exptData$sampleId), ] %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::mutate(
      profileName = paste(sampleId, profileType, sep = "_"),
      bwFile = paste(dataPath, "/", sampleId, "/", sampleId, "_normalized.bw", sep = ""),
      matFile = dplyr::case_when(
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
      mergedDataFile = dplyr::if_else(
        IP_tag == "polII",
        "NA",
        paste(dataPath, "/", sampleId, "/", sampleId, "_allGenes_clusters.tab", sep = "")
      ),
      narrowpeakFile = dplyr::case_when(
        IP_tag != "polII" & isTRUE(macs2Control) ~
          paste(dataPath, "/", sampleId, "/", sampleId, "_withCtrl_peaks.narrowPeak", sep = ""),
        IP_tag != "polII" & isFALSE(macs2Control) ~
          paste(dataPath, "/", sampleId, "/", sampleId, "_withoutCtrl_peaks.narrowPeak", sep = ""),
        TRUE ~ "NA"
      ),
      narrowpeakAnno = dplyr::case_when(
        IP_tag != "polII" & isTRUE(macs2Control) ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withCtrl.narrowPeak.nearestCDS.tab", sep = ""),
        IP_tag != "polII" & isFALSE(macs2Control) ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withoutCtrl.narrowPeak.nearestCDS.tab", sep = ""),
        TRUE ~ "NA"
      ),
      broadpeakFile = dplyr::case_when(
        IP_tag != "polII" & isTRUE(macs2Control) ~
          paste(dataPath, "/", sampleId, "/", sampleId, "_withCtrl_peaks.broadPeak", sep = ""),
        IP_tag != "polII" & isFALSE(macs2Control) ~
          paste(dataPath, "/", sampleId, "/", sampleId, "_withoutCtrl_peaks.broadPeak", sep = ""),
        TRUE ~ "NA"
      ),
      broadpeakAnno = dplyr::case_when(
        IP_tag != "polII" & isTRUE(macs2Control) ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withCtrl.broadPeak.nearestCDS.tab", sep = ""),
        IP_tag != "polII" & isFALSE(macs2Control) ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withoutCtrl.broadPeak.nearestCDS.tab", sep = ""),
        TRUE ~ "NA"
      ),
      tfPeakFile = dplyr::case_when(
        IP_tag != "polII" & isTRUE(macs2Control) ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withCtrl_peaks.annotated.tab", sep = ""),
        IP_tag != "polII" & isFALSE(macs2Control) ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withoutCtrl_peaks.annotated.tab", sep = ""),
        TRUE ~ "NA"
      )
    ) %>%
    dplyr::mutate_at(c("polIIExpFile", "polIIExpMat", "clusterFile", "tfPeakFile", "mergedDataFile", "narrowpeakFile", "narrowpeakAnno", "broadpeakFile", "broadpeakAnno"),
                     .funs = dplyr::funs(dplyr::na_if(., "NA")))


  return(exptData)
}

##################################################################################

