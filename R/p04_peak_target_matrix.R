#' Peak target gene matrix for multiple samples
#'
#' This function generates gene based peak target matrix for multiple samples. If
#' \code{sampleInfo} has just one row for a sample, this function will return
#' per gene peak target for just one sample. If there are multiple samples,
#' a full join is performed to extract peaks in all samples for each gene.
#'
#' A gene can have two types of peaks, one at TSS and another at TES. \code{position}
#' option allows to select the peak target in following manner.
#' \itemize{
#' \item \strong{TSS:} select genes with peak signal upstream or near TSS
#' \item \strong{TES:} select genes with peak signal near end of the gene
#' \item \strong{best:} If a gene have both TSS and TES region peak, select best based
#' on log10Pvalue of peak caller
#' \item \strong{both:} If a gene have both TSS and TES region peak, select both peaks
#' }
#'
#' @param sampleInfo sample information dataframe
#' @param position One of \code{c("TSS", "TES", "best", "both")}. This is used
#' for filtering the target genes for each sample. Default: "TSS". Refer to Details
#' section for more information.
#'
#' @return A gene wise peak target dataframe for multiple samples
#' @export
#'
#' @examples NA
peak_target_matrix <- function(sampleInfo, position = "TSS"){

  position <- match.arg(arg = position, choices = c("TSS", "TES", "best", "both"))

  mergedDf <- NULL
  peakPositionCol <- NULL
  colNames <- c("hasPeak", "peakId", "peakType", "peakPosition", "peakPval", "peakDist")

  for(i in 1:nrow(sampleInfo)){

    selectCols <- paste(colNames, sampleInfo$sampleId[i], sep = ".")
    peakPositionCol <- paste("peakPosition", sampleInfo$sampleId[i], sep = ".")

    df <- suppressMessages(readr::read_tsv(file = sampleInfo$peakTargetFile[i])) %>%
      dplyr::select(geneId, !!!selectCols)

    joinBy <- c("geneId" = "geneId")

    if(position == "best"){
      ## choose one from TSS and TES based on log10Pvalue
      df <- dplyr::group_by(df, geneId) %>%
        dplyr::arrange_at(.vars = vars(starts_with("peakPval")), .by_group = TRUE) %>%
        dplyr::slice(1L) %>%
        dplyr::ungroup()

    } else if(position == "both"){
      ## use geneId + peakPosition columns for full_join
      joinBy <- structure(c("geneId", peakPositionCol), names = c("geneId", "peakPosition"))

    } else{
      ## only TSS/TES
      df <- dplyr::filter(df, !!sym(peakPositionCol) == position)
    }


    if(is.null(mergedDf)){
      ## first sample
      mergedDf <- df
      if(position == "both"){
        mergedDf$peakPosition <- mergedDf[[peakPositionCol]]
      }

    } else{
      ## very IMP to use full join instead of left join
      ## left_join will just add rows which exist in 1st sample
      ## full_join will ensure that all the unique rows from subsequent samples are present
      mergedDf <- dplyr::full_join(x = mergedDf, y = df, by = joinBy)

    }

  }

  ## set hasPeak NA values to FALSE
  mergedDf <- dplyr::mutate_at(
    .tbl = mergedDf,
    .vars = vars(starts_with("hasPeak")),
    .funs = funs(if_else(is.na(.), true = FALSE, false = .))
  )

  return(mergedDf)

}


