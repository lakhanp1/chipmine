% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p08_average_profile_functions.R
\name{draw_avg_profile_plot}
\alias{draw_avg_profile_plot}
\title{Average signal plot for multiple samples}
\usage{
draw_avg_profile_plot(
  exptInfo,
  profileMats,
  genes,
  lineColors = NULL,
  lineShape = NULL
)
}
\arguments{
\item{exptInfo}{Experiment information dataframe}

\item{profileMats}{a named list of profile matrix}

\item{genes}{a vector of gene IDs. Profile matrix will be subset based on this
gene list to calculate average profile}

\item{lineColors}{color for each sample}

\item{lineShape}{linetype for each sample}
}
\value{
a ggplot object with line plot
}
\description{
Average signal plot for multiple samples
}
\examples{
NA
}
