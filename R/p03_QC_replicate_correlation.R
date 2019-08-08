

#' TF ChIPseq replicate correlation plots
#'
#' This function generate the peak overlap statistics table for two TF ChIPseq
#' replicates. Using common peaks, it plots the distribution and scatter plot of
#' either -log10(p-value) or fold enrichment from peak file. Macs2 peak caller
#' output files in narrowPeak/broadPeak format are used to extract the values.
#'
#' @param sampleInfo sample information dataframe
#' @param compare One of \emph{pvalue} or \emph{enrichment}
#' @param title Comparison title
#' @param yintercept yintercept for horizontal line on beeswarm plot
#'
#' @return A list with following ggplot2 elements:
#' \itemize{
#' \item \strong{table:} peak overlap summary table for replicates as a \code{ggtable} object
#' \item \strong{distribution:} Beeswarm plots showing comparative distribution between
#' replicates for a metric of choice
#' \item \strong{valueScatter:} X-Y scatter plot using values a metric of choice
#' \item \strong{rankScatter:} X-Y scatter plot using \code{rank(values)} a metric of choice
#' }
#'
#' @export
#'
#' @examples NA
tf_replicate_plots <- function(sampleInfo, compare = "pvalue", title, yintercept = -Inf){

  compare <- match.arg(tolower(compare), choices = c("pvalue", "enrichment", "qvalue"))

  compareConfig <- list(
    pvalue = list(narrowPeak = "pValue",
                  chipmine = "peakPval",
                  plotLabel = "-log10(p-value)",
                  color = "#a2b67c"),
    enrichment = list(narrowPeak = "signalValue",
                      chipmine = "peakEnrichment",
                      plotLabel = "fold enrichment",
                      color = "#bf9a2d"),
    pvalue = list(narrowPeak = "qValue",
                  chipmine = "peakQval",
                  plotLabel = "-log10(q-value)",
                  color = "#a2b67c")
  )

  if(nrow(sampleInfo) != 2){
    stop("sampleInfo should have 2 rows for 2 replicates to be compared")
  }

  rep1Id <- sampleInfo$sampleId[1]
  rep2Id <- sampleInfo$sampleId[2]

  sampleInfoList <- purrr::transpose(sampleInfo)  %>%
    purrr::set_names(nm = purrr::map(., "sampleId"))

  rep1Data <- sampleInfoList[[ rep1Id ]]
  rep2Data <- sampleInfoList[[ rep2Id ]]

  tfCols <- sapply(
    X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
          "peakType", "bidirectional", "targetOverlap", "peakOverlap", "relativeSummitPos", "peakRegion",
          "peakPosition", "peakCoverage", "pvalFiltered", "summitSeq", "relativePeakPos"),
    FUN = function(x){ structure(paste(x, ".", sampleInfo$sampleId, sep = ""), names = sampleInfo$sampleId) },
    simplify = F, USE.NAMES = T)


  ## find the common and unique peaks between replicates
  peakFormat <- dplyr::case_when(
    rep1Data$peakType == "narrow" ~ "narrowPeak",
    rep1Data$peakType == "broad" ~ "broadPeak"
  )

  r1Peaks <- rtracklayer::import(con = rep1Data$peakFile, format = peakFormat)
  r2Peaks <- rtracklayer::import(con = rep2Data$peakFile, format = peakFormat)

  mcols(r1Peaks)$overlap <- factor("unique", levels = c("common", "unique"), ordered = TRUE)
  mcols(r2Peaks)$overlap <- factor("unique", levels = c("common", "unique"), ordered = TRUE)

  peakOvlp <- findOverlaps(r1Peaks, r2Peaks)

  mcols(r1Peaks)$overlap[peakOvlp@from] <- "common"
  mcols(r2Peaks)$overlap[peakOvlp@to] <- "common"

  #######################################
  ovlpSummary <- bind_rows(table(r1Peaks$overlap), table(r2Peaks$overlap)) %>%
    dplyr::mutate(total = common + unique,
                  sampleId = c(rep1Data$sampleId, rep2Data$sampleId)) %>%
    dplyr::mutate(common_per = sprintf(fmt = "%d (%.2f%%)", common, (100 * common / total)),
                  unique_per = sprintf(fmt = "%d (%.2f%%)", unique, (100 * unique / total))) %>%
    dplyr::select(sampleId, common = common_per, unique = unique_per, total)


  ## summary table
  gg_stable <- ggpubr::ggtexttable(ovlpSummary, rows = NULL,
                                   theme = ttheme("mOrange"))


  r1PeaksDf <- as.data.frame(r1Peaks) %>% dplyr::mutate(sampleId = rep1Data$sampleId)
  r2PeaksDf <- as.data.frame(r2Peaks) %>% dplyr::mutate(sampleId = rep2Data$sampleId)

  #######################################
  ## set outliers to 99.5 quantile
  if(nrow(r1PeaksDf) > 20){
    r1PeaksDf$signalValue <- pmin(r1PeaksDf$signalValue, quantile(r1PeaksDf$signalValue, 0.995))
    r1PeaksDf$pValue <- pmin(r1PeaksDf$pValue, quantile(r1PeaksDf$pValue, 0.99))
  }

  if(nrow(r2PeaksDf) > 20){
    r2PeaksDf$signalValue <- pmin(r2PeaksDf$signalValue, quantile(r2PeaksDf$signalValue, 0.995))
    r2PeaksDf$pValue <- pmin(r2PeaksDf$pValue, quantile(r2PeaksDf$pValue, 0.99))
  }

  plotData <- dplyr::bind_rows(r1PeaksDf, r2PeaksDf)

  pointColor <- structure(c("green", "red"), names = c("common", "unique"))
  pointAlpha <- structure(c(0.25, 1), names = c("common", "unique"))


  ## density scatter plot
  gg_dot <- ggplot(
    data = plotData,
    mapping = aes(x = sampleId, y = !!sym(compareConfig[[compare]]$narrowPeak)) ) +
    geom_hline(yintercept = yintercept, color = "black", linetype = "dashed") +
    ggbeeswarm::geom_quasirandom(mapping = aes(color = overlap, alpha = overlap)) +
    geom_boxplot(width=0.1, fill = NA, outlier.colour = NA, color = alpha("green", 1)) +
    scale_color_manual(name = "Replicate overlap",
                       values = pointColor) +
    scale_alpha_manual(values = pointAlpha) +
    labs(title = paste(compareConfig[[compare]]$plotLabel, "distribution"),
         y = compareConfig[[compare]]$plotLabel) +
    guides(alpha = FALSE,
           color = guide_legend(override.aes = list(size = 4))) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.title = element_text(size = 14),
      axis.text.y = element_text(size = 16),
      axis.text.x = element_text(size = 8),
      panel.grid = element_blank(),
      title = element_text(hjust = 0.5, size = 14),
      legend.position = "bottom"
    )

  #######################################
  ## plot XY scatter plots
  theme_scatter <- theme_bw() +
    theme(
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 14),
      panel.grid = element_blank()
    )

  peakOverlapMat <- combinatorial_binding_matrix(sampleInfo = sampleInfo,
                                                 peakFormat = peakFormat)


  commonPeaks <- dplyr::filter_at(.tbl = peakOverlapMat,
                                  .vars = vars(starts_with("overlap.")),
                                  .vars_predicate = all_vars(. == TRUE))

  rep1Col <- tfCols[[compareConfig[[compare]]$chipmine]][rep1Id]
  rep2Col <- tfCols[[compareConfig[[compare]]$chipmine]][rep2Id]

  densColor <- densCols(x = commonPeaks[[ rep1Col ]],
                        y = commonPeaks[[ rep2Col ]],
                        colramp = colorRampPalette(c("black", "white"))
  )

  commonPeaks$density <- col2rgb(densColor)[1,] + 1L

  ## value XY scatter plot
  gg_scatter_val <- ggplot(
    data = commonPeaks,
    mapping = aes(x = !!sym(rep1Col), y = !!sym(rep2Col), color = density)
  ) +
    geom_point() +
    geom_smooth(method=lm, se = FALSE, formula = y ~ x) +
    ggpubr::stat_cor(method = "spearman", size = 10) +
    scale_color_gradientn(
      name = "Density",
      colours = RColorBrewer::brewer.pal(7, "RdYlBu")) +
    labs(title = paste(compareConfig[[compare]]$plotLabel, "scatter plot")) +
    # guides(color = FALSE) +
    theme_bw() +
    theme_scatter

  ## rank(-value) XY scatter plot
  gg_scatter_rank <- ggplot(
    data = commonPeaks,
    mapping = aes(x = rank(- !!sym(rep1Col)), y = rank(- !!sym(rep2Col)), color = density)
  ) +
    geom_point() +
    scale_color_gradientn(
      name = "Density",
      colours = RColorBrewer::brewer.pal(7, "RdYlBu")) +
    labs(title = paste(compareConfig[[compare]]$plotLabel, "rank scatter plot")) +
    # guides(color = FALSE) +
    theme_bw() +
    theme_scatter

  #######################################
  # ## IDR analysis
  # idrData <- cbind(commonPeaks[[rep1Col]], commonPeaks[[rep2Col]])
  #
  # idrRes <- idr::est.IDR(idrData, mu=3, sigma=1, rho=.9, p=.5)
  #
  # commonPeaks$idr <- idrRes$idr
  #
  # ## value scatter plot
  # gg_scatter_val_idr <- ggplot(data = commonPeaks,
  #                              mapping = aes(x = !!sym(rep1Col),
  #                                            y = !!sym(rep2Col),
  #                                            color = idr)
  # ) +
  #   geom_point() +
  #   geom_smooth(method=lm, se = FALSE, formula = y ~ x) +
  #   ggpubr::stat_cor(method = "pearson") +
  #   scale_color_gradientn(
  #     name = "IDR",
  #     colours = rev(RColorBrewer::brewer.pal(7, "YlOrBr"))) +
  #   labs(title = paste(compareConfig[[compare]]$plotLabel, "scatter plot")) +
  #   theme_bw() +
  #   theme_scatter
  #
  # ## rank(-value) scatter plot
  # gg_scatter_rank_idr <- ggplot(data = commonPeaks,
  #                               mapping = aes(x = rank(- !!sym(rep1Col)),
  #                                             y = rank(- !!sym(rep2Col)),
  #                                             color = idr)
  # ) +
  #   geom_point() +
  #   scale_color_gradientn(
  #     name = "IDR",
  #     colours = rev(RColorBrewer::brewer.pal(7, "YlOrBr"))) +
  #   labs(title = paste(compareConfig[[compare]]$plotLabel, "rank scatter plot")) +
  #   theme_bw() +
  #   theme_scatter

  #######################################

  return(list(
    table = gg_stable,
    distribution = gg_dot,
    valueScatter = gg_scatter_val,
    rankScatter = gg_scatter_rank
  ))

}
















