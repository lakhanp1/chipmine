
##################################################################################

#' Calculate mean coverage for regions.
#'
#' This function calculate the coverage for the regions provided as bed/narrowPeak/broadPeak
#' format file from a bigWig file. For each region, it adds the depth and divide it by the
#' region length. NOTE: It assumes that the bigWig file is already scalled (e.g. 1 million)
#' and hence does perform scalling or normalization.
#'
#' @param regions Either GenomicRanges object or a file in bed, narrowPeak or broadPeak format
#' @param format format of the peak file. Should be one of \emph{bed}, \emph{narrowPeak} or
#' \emph{broadPeak}. This options must be provided if \code{regions} is a file.
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

    mcols(regionGr[seqnames(regionGr) == chr])[[ name ]] <- as.numeric(sprintf("%.3f", chrScore))

  }

  if(returnType == "df"){
    return(as.data.frame(mcols(regionGr)))
  } else{
    return(regionGr)
  }

}


##################################################################################

#' Coverage matrix for regions of interest
#'
#' This function calls \code{region_coverage()} interanlly to calculate coverage for
#' each sample's bigWig file. Remember that the bigWig files should be normalized to
#' make the coverage scores comparable. Refer to \code{region_coverage()} documentation
#' for the details on coverage calculation.
#'
#' @param regions A GRanges object
#' @param exptInfo A dataframe with sample information. sampleId and bwFile columns
#' are used for coverage calculation.
#'
#' @return A coverage matrix.
#' @export
#'
#' @examples NA
region_coverage_matrix <- function(regions, exptInfo){

  if(!any(class(regions) %in% c("GRanges"))){
    stop("regions should be a GRanges  object")
  }

  if(is.null(mcols(regions)$name)){
    mcols(regions)$name <- names(regions)
  }

  ## calculate coverage
  for(i in 1:nrow(exptInfo)){
    regions <- region_coverage(regions = regions, bwFile = exptInfo$bwFile[i], name = exptInfo$sampleId[i])
  }

  covMat <- as.data.frame(mcols(regions)) %>%
    dplyr::select(name, !!!exptInfo$sampleId)

  return(covMat)
}

##################################################################################

