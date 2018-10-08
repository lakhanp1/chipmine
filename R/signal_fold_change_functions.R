


#' Raw log2 fold change value between two samples
#'
#' @param df data frame
#' @param s1 column name for sample1
#' @param s2 column name for sample2
#' @param newCol new column name for log2(fold_change) values. Default: lfc
#' @param isExpressedCols A named vector with column names for isExpressed data for
#' both polII samples
#'
#' @return dataframe with additional fold change column
#' @export
#'
#' @examples NA
get_fold_change <- function(df, s1, s2, newCol = "lfc", isExpressedCols = NULL){

  df <- dplyr::mutate(df,
                      !! newCol := log2(!!as.name(s1) + 1) - log2(!!as.name(s2) + 1))

  ## LFC correction: set lfc = 0 if both the signal values are not in top 10%
  if(!is.null(isExpressedCols)){

    if(! all(c(s1, s2) %in% names(isExpressedCols))){
      stop("Missing is_expressed.* columns for given samples in isExpressedCols vector.")
    }

    df <- dplyr::mutate(
      df,
      !! newCol := if_else(
        condition = !! as.name(isExpressedCols[[s1]]) == TRUE | !! as.name(isExpressedCols[[s2]]) == TRUE,
        true = !! as.name(newCol),
        false = 0))
  }

  return(df)
}


