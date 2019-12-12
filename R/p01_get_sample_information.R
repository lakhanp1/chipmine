
##################################################################################
## generate experiment data
#' Create sample information dataframe
#'
#' This is a utility function to prepare the sample information dataframe for sample
#' of interest using experiment information file. Experiment information file should
#' have columns: \code{sampleId, IP_tag, control, peakType}. \cr\cr
#' \strong{Possible values for compulsory columns and how they are handled:}
#' \itemize{
#' \item \strong{IP_tag:} One of \code{c("HA", "FLAG", "HIS", "MYC", "TAP", "TF", "polII")}.
#' This column value is used to determine whether the ChIPseq data is for polII ChIPseq
#' or transcription factor ChIPseq. If polII ChIPseq, \code{polIIExpFile, polIIExpMat}
#' columns are populated with appropriate file names for polII data.
#' \item \strong{peakType:} One of \code{c("narrow", "broad")}. This column value is used
#' to decide the macs2 peak output file name (i.e. narrowPeak or broadPeak).
#' \item \strong{control:} This column value is used to decide the suffix of macs2 output
#' files. If the value is blank or ".", \emph{withoutCtrl} is used in macs2 peak realated
#' file names and \emph{withCtrl} otherwise.
#' }
#'
#' @param exptInfoFile A tabular file with details for each sample. "sampleId" and "IP_tag" columns are must
#' @param dataPath Path where the data is stored
#' @param profileMatrixSuffix A character string to be used in the file name of profile
#' matrix file. Generating profile matrix takes long time and hence it is efficient to save
#' the profile matrix. If profile matrix is saved, this suffix is used to design the name of
#' profile matrix file. This will generate a file name "sample_ID.normalizedMatrix.tab.gz".
#' Different intutuve suffixes can be used to store profile matrix like
#' \emph{TSS_4kb_2kb_normalized, TSS_3kb_3kb_normalized, normalizedmatrix_5kb}, etc.
#' Default: \emph{normalizedmatrix}.
#' @param samples A vector of Sample IDs which are to be processed. Default: All samples are used
#' @param profileType Type of profile. This will be added as suffix to the profile name
#'
#' @return Data frame with columns \code{"sampleId", "IP_tag", "control", "peakType",
#' "sampleName", "profileName", "bwFile", "matFile", "clusterFile", "mergedDataFile",
#' "polIIExpFile", "polIIExpMat", "peakFile", "peakAnno", "peakTargetFile", "FE_bwFile"}.
#' @export
#'
#' @examples NA
get_sample_information <- function(exptInfoFile, samples = NULL, dataPath,
                                   profileMatrixSuffix = "normalizedmatrix",
                                   profileType = "profile"){

  tfChipTags <- c("HA", "FLAG", "HIS", "MYC", "TAP", "TF")
  ## read the experiment sample details and select only those which are to be plotted
  exptData <- suppressMessages(readr::read_tsv(file = exptInfoFile))

  if(!is.null(samples)){
    exptData <- dplyr::left_join(x = data.frame(sampleId = samples, stringsAsFactors = F),
                                 y = exptData, by = "sampleId")
  }

  if (nrow(exptData) == 0) {
      return(NULL)
  }

  exptData$sampleId <- factor(exptData$sampleId, levels = unique(exptData$sampleId))


  if( !any("sampleName" %in% colnames(exptData)) ){
    exptData$sampleName <- exptData$sampleId
  }

  if( !any("peakType" %in% colnames(exptData)) ){
    warning("peakType column not found in sample information. Using \"narrow\" for all samples")
    exptData$peakType <- "narrow"
  }

  if( !any("control" %in% colnames(exptData)) ){
    exptData$control <- "."
  }

  exptData <- exptData[order(exptData$sampleId), ] %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::mutate(
      IP_tag = toupper(IP_tag),
      hasControl = if_else(condition = is.na(control) | control == ".",
                           true = FALSE, false = TRUE, missing = FALSE)
    ) %>%
    ## generic data
    dplyr::mutate(
      profileName = paste(sampleId, profileType, sep = "_"),
      bwFile = paste(dataPath, "/", sampleId, "/", sampleId, "_normalized.bw", sep = ""),
      matFile = paste(dataPath, "/", sampleId, "/", sampleId, ".",
                      profileMatrixSuffix,".tab.gz", sep = ""),
      clusterFile = dplyr::if_else(
        toupper(IP_tag) == "POLII",
        "NA",
        paste(dataPath, "/", sampleId, "/", sampleId, ".kmeans.clusters.txt", sep = "")
      ),
      mergedDataFile = dplyr::if_else(
        toupper(IP_tag) == "POLII",
        "NA",
        paste(dataPath, "/", sampleId, "/", sampleId, ".allGenes_clusters.tab", sep = "")
      )
    ) %>%
    ## polII ChIPseq related data
    dplyr::mutate(
      polIIExpFile = dplyr::if_else(
        toupper(IP_tag) == "POLII",
        paste(dataPath, "/", sampleId, "/", sampleId, ".normalizedExpression.tab", sep = ""),
        "NA"
      ),
      polIIExpMat = dplyr::if_else(
        toupper(IP_tag) == "POLII",
        paste(dataPath, "/", sampleId, "/", sampleId, "_polii_expr.tab.rel.mat", sep = ""),
        "NA"
      )
    ) %>%
    ## TF ChIP related data
    dplyr::mutate(
      peakFile = dplyr::case_when(
        IP_tag %in% !!tfChipTags & hasControl & peakType == "narrow" ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withCtrl_peaks.narrowPeak", sep = ""),
        IP_tag %in% !!tfChipTags & !hasControl & peakType == "narrow" ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withoutCtrl_peaks.narrowPeak", sep = ""),
        IP_tag %in% !!tfChipTags & hasControl & peakType == "broad" ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withCtrl_peaks.broadPeak", sep = ""),
        IP_tag %in% !!tfChipTags & !hasControl & peakType == "broad" ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withoutCtrl_peaks.broadPeak", sep = "")
      ),
      peakAnno = dplyr::case_when(
        IP_tag %in% !!tfChipTags & hasControl & peakType == "narrow" ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withCtrl.narrowPeak.annotation.tab", sep = ""),
        IP_tag %in% !!tfChipTags & !hasControl & peakType == "narrow" ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withoutCtrl.narrowPeak.annotation.tab", sep = ""),
        IP_tag %in% !!tfChipTags & hasControl & peakType == "broad" ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withCtrl.broadPeak.annotation.tab", sep = ""),
        IP_tag %in% !!tfChipTags & !hasControl & peakType == "broad" ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withoutCtrl.broadPeak.annotation.tab", sep = "")
      ),
      peakTargetFile = dplyr::case_when(
        IP_tag %in% !!tfChipTags & hasControl ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withCtrl_peaks.targets.tab", sep = ""),
        IP_tag %in% !!tfChipTags & !hasControl ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withoutCtrl_peaks.targets.tab", sep = ""),
        TRUE ~ "NA"
      ),
      FE_bwFile = dplyr::case_when(
        IP_tag %in% !!tfChipTags & hasControl ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withCtrl.FE.bw", sep = ""),
        IP_tag %in% !!tfChipTags & !hasControl ~
          paste(dataPath, "/", sampleId, "/", sampleId, ".withoutCtrl.FE.bw", sep = ""),
        TRUE ~ "NA"
      )
    ) %>%
    dplyr::mutate_at(
      .vars = c("polIIExpFile", "polIIExpMat", "clusterFile", "peakTargetFile", "mergedDataFile",
                "peakFile", "peakAnno"),
      .funs = list(~ dplyr::na_if(., "NA"))
    ) %>%
    dplyr::select(-hasControl)


  return(exptData)
}

##################################################################################

