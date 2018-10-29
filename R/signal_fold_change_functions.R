


#' Raw log2 fold change value between two samples as log2(s2/s1)
#'
#' @param df data frame
#' @param nmt column name for sample in numerator
#' @param dmt column name for sample in denominator
#' @param newCol new column name for log2(nmt/dmt) values. Default: lfc
#' @param isExpressedCols A named vector with column names for isExpressed data for
#' both polII samples
#'
#' @return dataframe with additional fold change column
#' @export
#'
#' @examples NA
get_fold_change <- function(df, nmt, dmt, newCol = "lfc", isExpressedCols = NULL){

  df <- dplyr::mutate(df,
                      !! newCol := log2(!!as.name(nmt) + 1) - log2(!!as.name(dmt) + 1))

  ## LFC correction: set lfc = 0 if both the signal values are not in top 10%
  if(!is.null(isExpressedCols)){

    if(! all(c(dmt, nmt) %in% names(isExpressedCols))){
      stop("Missing is_expressed.* columns for given samples in isExpressedCols vector.")
    }

    df <- dplyr::mutate(
      df,
      !! newCol := if_else(
        condition = !! as.name(isExpressedCols[[dmt]]) == TRUE | !! as.name(isExpressedCols[[nmt]]) == TRUE,
        true = !! as.name(newCol),
        false = 0))
  }

  return(df)
}


