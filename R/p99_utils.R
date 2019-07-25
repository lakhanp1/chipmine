

#' Extract the TES position as GenomicRanges
#'
#' This function extract the TES regions from the GenomicRanges. Internally, it invert the strand
#' information of GenomicRanges and then use the promoters() method to extract the promoter
#' region for this flipped ranges. Finally the strand information is switched back again.
#'
#' @param gr A GenomicRange object
#' @param up upstream region length. This is 'upstream' argument of  GenomicRanges::promoters
#' method
#' @param down downstream region length. This is 'downstream' argument of  GenomicRanges::promoters
#' method
#' @param ... Other arguments to GenomicRanges::promoters function
#'
#' @return A GenomicRange object with TES position
#' @export
#'
#' @examples NULL
get_TES <- function(gr, up = 0, down = 1, ...){

  ## flip the strand of original GRanges
  rg <- GenomicRanges::invertStrand(gr)

  ## extract TSS position of this flipped GRanges
  tesGr <- GenomicRanges::promoters(x = rg, upstream = up, downstream = down, ...)

  ## reset the strand
  tesGr <- GenomicRanges::invertStrand(tesGr)

  return(tesGr)
}
