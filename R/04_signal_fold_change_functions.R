


#' Raw log2 fold change value between two samples as log2(s2/s1)
#'
#' For each sample, any value less than 1 is set to 1 to avoid negative
#' log2 value.
#'
#' @param df data frame
#' @param nmt column name for sample in numerator
#' @param dmt column name for sample in denominator
#' @param newCol new column name for \code{log2(nmt/dmt)} values. Default: lfc
#' @param isExpressedCols A named vector with column names for isExpressed data for
#' both polII samples
#' @param lfcLimit All the LFC values in range \code{[-lfcLimit, lfcLimit]} are set
#' to 0. Default: 0 i.e. no limit
#'
#' @return dataframe with additional fold change column
#' @export
#'
#' @examples NA
get_fold_change <- function(df, nmt, dmt, newCol = "lfc",
                            isExpressedCols = NULL, lfcLimit = 0){

  pseudoCount <- 1

  df <- dplyr::mutate(
    df,
    !! newCol := log2(pmax(!!as.name(nmt), pseudoCount)) - log2(pmax(!!as.name(dmt), pseudoCount)))

  ## set LFC values within range (-lfcLimit < LFC < +lfcLimit) to zero
  df[[newCol]][
    dplyr::between(x = df[[newCol]], left = -1 * lfcLimit, right = lfcLimit)
    ] <- 0

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


