% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p09_coverage_functions.R
\name{region_coverage_matrix}
\alias{region_coverage_matrix}
\title{Coverage matrix for regions of interest}
\usage{
region_coverage_matrix(regions, exptInfo)
}
\arguments{
\item{regions}{A GRanges object}

\item{exptInfo}{A dataframe with sample information. sampleId and bwFile columns
are used for coverage calculation.}
}
\value{
A coverage matrix.
}
\description{
This function calls \code{region_coverage()} interanlly to calculate coverage for
each sample's bigWig file. Remember that the bigWig files should be normalized to
make the coverage scores comparable. Refer to \code{region_coverage()} documentation
for the details on coverage calculation.
}
\examples{
NA
}
