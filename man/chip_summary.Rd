% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p03_QC_ChIPseq_summary.R
\name{chip_summary}
\alias{chip_summary}
\title{TF ChIPseq data summary plots}
\usage{
chip_summary(
  sampleId,
  peakAnnotation,
  peakFile,
  peakFormat,
  pval_cutoff = 20,
  enrichment_cutoff = 3,
  markTargets = NULL,
  pointColor = NULL,
  pointAlpha = NULL
)
}
\arguments{
\item{sampleId}{sampleId}

\item{peakAnnotation}{macs2 peak annotation file generated by \code{narrowPeak_annotate()}}

\item{peakFile}{narrowPeak or broadPeak file generated by macs2}

\item{peakFormat}{One of narrowPeak or broadPeak}

\item{pval_cutoff}{A cutoff for -log10(p-value) for drawing horizontal line in beeswarm
plot. Default: 20}

\item{enrichment_cutoff}{A cutoff for peak enrichment for drawing horizontal line in
beeswarm plot. Default: 3}

\item{markTargets}{A named list of target genes which needs to be annotated in the
scatter plot. Default: NULL}

\item{pointColor}{A named vector with colors for each target category in \code{markTargets}
list. Default: NULL}

\item{pointAlpha}{A named vector of alpha vlues for each target category in
\code{markTargets} list. Default: NULL}
}
\value{
A list with following elements is returned.
\itemize{
\item \strong{data:} plot data used for plotting
\item \strong{figure:} A combined figure generated by \code{ggpubr::ggarrange}
\item \strong{plots:} Individual ggplot list for each plot.
\itemize{
\item \strong{table:} a \code{ggtable} object for quantile stats
\item \strong{distribution:} A list with two beeswarm plot objects: enrichment, pval.
\item \strong{annoPie:} A pie chart for peak annotation distribution
}}
}
\description{
This method generates various summary plots for TF ChIPseq data. \cr
\itemize{
\item A table summarizing quantils of peak enrichment, p-value and q-value
\item A distribution beeswarm plot for peak enrichment values.
\item A distribution beeswarm plot for peak p-value.
\item A distribution beeswarm plot for peak width.
\item A pie chart showing peak annotation distribution.
}
If there is no peak information, it returns blank plot with "No Data" text
}
\examples{
NA
}
