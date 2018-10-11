
##################################################################################

#' Calculate mean coverage for regions.
#'
#' This function calculate the coverage for the regions provided as bed/narrowPeak/broadPeak
#' format file from a bigWig file. For each region, it adds the depth and divide it by the
#' region length. NOTE: It assumes that the bigWig file is already scalled (e.g. 1 million)
#' and hence does perform scalling or normalization.
#'
#' @param regions Either GenomicRanges object or a file in bed, narrowPeak or broadPeak format
#' @param format format of the peak file. Should be one of "bed", "narrowPeak" or "broadPeak".
#' This option does not matter if the 'regions' argument if a GRanges object.
#' @param bwFile bigWig file from which the coverage will be calculated
#' @param name Name of the coverage column. Default: coverage
#'
#' @return If the regions is a GRanges object, then same object is returned with additional
#' coverage column. Otherwise a dataframe is returned.
#' @export
#'
#' @examples
region_coverage <- function(regions, format = "narrowPeak", bwFile, name = "coverage"){

  regionGr <- NULL
  returnType <- "df"

  if(any(class(regions) %in% c("GRanges"))){
    regionGr <- regions
    returnType <- "GRanges"
  } else{
    regionGr <- rtracklayer::import(con = regions, format = format)
  }

  bwGr <- rtracklayer::import.bw(con = bwFile, as = "RleList")

  mcols(regionGr)[[ name ]] <- 0

  ## get coverage for each chromosome
  for (chr in seqinfo(regionGr)@seqnames) {
    chrScore <- sum(IRanges::Views(
      subject = bwGr[[chr]],
      ranges(regionGr[seqnames(regionGr) == chr]))) / width(regionGr[seqnames(regionGr) == chr])

    mcols(regionGr[seqnames(regionGr) == chr])[[ name ]] <- chrScore

  }

  if(returnType == "df"){
    return(as.data.frame(mcols(regionGr)))
  } else{
    return(regionGr)
  }

}


##################################################################################

