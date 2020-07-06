

#' TF ChIPseq replicate correlation plots
#'
#' This function generate the peak overlap statistics table for two TF ChIPseq
#' replicates. Using common peaks, it plots the distribution and scatter plot of
#' either -log10(p-value) or fold enrichment from peak file. Macs2 peak caller
#' output files in narrowPeak/broadPeak format are used to extract the values.
#'
#' @param sampleInfo sample information dataframe
#' @param compare One of \emph{pvalue} or \emph{enrichment}
#' @param yintercept yintercept for horizontal line on beeswarm plot
#' @inheritParams combinatorial_binding_matrix
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
compare_ChIPseq_replicates <- function(sampleInfo, compare = "pvalue", yintercept = -Inf,
                                       summitRegion = 0, peakFormat){

  compare <- match.arg(tolower(compare), choices = c("pvalue", "enrichment", "qvalue"))
  peakFormat <- match.arg(peakFormat, choices = c("narrowPeak", "broadPeak"))

  compareConfig <- list(
    pvalue = list(narrowPeak = "pValue",
                  chipmine = "peakPval",
                  plotLabel = "-log10(p-value)",
                  color = "#a2b67c"),
    enrichment = list(narrowPeak = "signalValue",
                      chipmine = "peakEnrichment",
                      plotLabel = "fold enrichment",
                      color = "#bf9a2d"),
    qValue = list(narrowPeak = "qValue",
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
          "peakAnnotation", "bidirectional", "targetOverlap", "peakOverlap", "relativeSummitPos", "peakRegion",
          "peakPosition", "peakCoverage", "pvalFiltered", "summitSeq", "relativePeakPos"),
    FUN = function(x){ structure(paste(x, ".", sampleInfo$sampleId, sep = ""), names = sampleInfo$sampleId) },
    simplify = F, USE.NAMES = T)

  ## find the common and unique peaks between replicates
  r1Peaks <- rtracklayer::import(con = rep1Data$peakFile, format = peakFormat)
  r2Peaks <- rtracklayer::import(con = rep2Data$peakFile, format = peakFormat)

  combinedSeqLevels <- union(seqlevels(r1Peaks), seqlevels(r2Peaks))
  seqlevels(r1Peaks) <- combinedSeqLevels
  seqlevels(r2Peaks) <- combinedSeqLevels

  mcols(r1Peaks)$overlap <- factor(
    x = rep("unique", length(r1Peaks)), levels = c("common", "unique"), ordered = TRUE
  )
  mcols(r2Peaks)$overlap <- factor(
    x = rep("unique", length(r2Peaks)), levels = c("common", "unique"), ordered = TRUE
  )

  peakOvlp <- findOverlaps(r1Peaks, r2Peaks)

  mcols(r1Peaks)$overlap[peakOvlp@from] <- "common"
  mcols(r2Peaks)$overlap[peakOvlp@to] <- "common"

  #######################################
  ovlpSummary <- bind_rows(table(r1Peaks$overlap), table(r2Peaks$overlap)) %>%
    tidyr::replace_na(replace = list(common = 0, unique = 0)) %>%
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
    scale_x_discrete(limits = sampleInfo$sampleId) +
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
                                                 peakFormat = peakFormat,
                                                 summitRegion = summitRegion)


  commonPeaks <- dplyr::filter_at(.tbl = peakOverlapMat,
                                  .vars = vars(starts_with("overlap.")),
                                  .vars_predicate = all_vars(. == TRUE))

  rep1Col <- tfCols[[compareConfig[[compare]]$chipmine]][rep1Id]
  rep2Col <- tfCols[[compareConfig[[compare]]$chipmine]][rep2Id]

  if(nrow(commonPeaks) >= 10){
    densColor <- densCols(x = commonPeaks[[ rep1Col ]],
                          y = commonPeaks[[ rep2Col ]],
                          colramp = colorRampPalette(c("black", "white"))
    )

    commonPeaks$density <- col2rgb(densColor)[1,] + 1L
  } else{
    commonPeaks$density <- rep(1, nrow(commonPeaks))
  }

  colorScale <- "viridis"

  if(nrow(commonPeaks) > 0){

    ## value XY scatter plot
    gg_scatter_val <- ggplot(
      data = commonPeaks,
      mapping = aes(x = !!sym(rep1Col), y = !!sym(rep2Col), color = density)
    ) +
      geom_point() +
      geom_smooth(method=lm, se = FALSE, formula = y ~ x, color = "red") +
      ggpubr::stat_cor(method = "spearman", size = 6) +
      viridis::scale_color_viridis(name = "Density", option = colorScale) +
      labs(title = paste(compareConfig[[compare]]$plotLabel, "scatter plot")) +
      theme_bw() +
      theme_scatter

    ## rank(-value) XY scatter plot
    gg_scatter_rank <- ggplot(
      data = commonPeaks,
      mapping = aes(x = rank(!!sym(rep1Col)), y = rank(!!sym(rep2Col)), color = density)
    ) +
      geom_point() +
      viridis::scale_color_viridis(name = "Density", option = colorScale) +
      labs(title = paste(compareConfig[[compare]]$plotLabel, "rank scatter plot")) +
      theme_bw() +
      theme_scatter

  } else{
    gg_scatter_val <- ggplot() +
      geom_text(mapping = aes(x = 0.5, y = 0.5), label = "No\ncommon\npeaks", size = 12) +
      theme_void()
    gg_scatter_rank <- ggplot() +
      geom_text(mapping = aes(x = 0.5, y = 0.5), label = "No\ncommon\npeaks", size = 12) +
      theme_void()
  }


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



##################################################################################

#' Correlation plots for numeric data
#'
#' @param data A dataframe
#' @param rep1Col column name for replicate 1
#' @param rep2Col column name for replicate 1
#' @param value A character name for the value. Default: FPKM
#' @param trans trans argument from ggplot::scale_x_continuous() function.
#' @param pseudoCount A pseudo count to avoid log(0) = Inf error. All values less than
#' this number will be set to this number in XY scatter plot. Default: 0
#'
#' @return A list with following elements is returned.
#' \itemize{
#' \item \strong{data:} plot data used for plotting
#' \item \strong{figure:} A combined figure generated by \code{ggpubr::ggarrange}
#' \item \strong{plots:} Individual ggplot list for each plot.
#' \itemize{
#' \item \strong{table:} a \code{ggtable} object for quantile stats
#' \item \strong{density:} a FPKM density plot for both replicates
#' \item \strong{scatter:} A list with two scatter plot objects: value, rank
#' \item \strong{corrVariation:} A line plot showing changing correlation with quantile
#' subset of data.
#' }
#' }
#'
#' @export
#'
#' @examples NA
compare_replicates <- function(data, rep1Col, rep2Col, value = "FPKM", trans = "identity",
                               pseudoCount = 0){

  transformer <- scales::as.trans(trans)
  colorScale <- "viridis"

  data <- data %>%
    dplyr::mutate(
      avgSignal = purrr::pmap_dbl(
        .l = list(!!!syms(c(rep1Col, rep2Col))),
        .f = purrr::lift_vd(mean))
    )

  ## calculate density
  if(nrow(data) >= 10){
    densColor <- densCols(x = transformer$transform(data[[rep1Col]]),
                          y = transformer$transform(data[[rep2Col]]),
                          colramp = colorRampPalette(c("black", "white"))
    )

    data$density <- col2rgb(densColor)[1,] + 1L
  } else{
    data$density <- rep(1, nrow(data))
  }


  theme_scatter <- theme_bw() +
    theme(
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 14),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(1, "cm")
    )

  ## FPKM density plot
  gg_fpkm_density <- dplyr::select(.data = data, geneId, rep1Col, rep2Col) %>%
    tidyr::pivot_longer(cols = c(rep1Col, rep2Col), names_to = "sampleId", values_to = "FPKM") %>%
    dplyr::mutate(sampleId = forcats::as_factor(sampleId)) %>%
    ggplot() +
    geom_density(mapping = aes(x = FPKM, color = sampleId), size = 1) +
    geom_density(mapping = aes(x = pmax(FPKM, 1), color = sampleId), linetype = "twodash") +
    scale_x_continuous(trans = transformer) +
    scale_color_brewer(name = NULL, palette="Dark2") +
    labs(x = paste(transformer$name, "(", value, ")", sep = ""), y = "Density", title = "Density plot") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14),
      panel.grid = element_blank()
    )


  ## value XY scatter plot
  ## there is some problem when geom_smooth() and coord_trans() used together
  ## using coord_trans() as it does not change the original values but just the scale
  gg_scatter_val <- ggplot(
    data = data,
    mapping = aes(
      x = pmax(!!pseudoCount, !!sym(rep1Col)),
      y = pmax(!!pseudoCount, !!sym(rep2Col)),
      color = density)
  ) +
    geom_point() +
    coord_trans(x = trans, y = trans) +
    # geom_smooth(method="lm", se = FALSE, formula = y ~ x, color = "red") +
    ggpubr::stat_cor(
      method = "pearson", size = 6, color = "red", label.x.npc = 0, label.y.npc = 1
    ) +
    viridis::scale_color_viridis(name = "Density", option = colorScale) +
    scale_x_continuous(
      breaks = trans_breaks(trans = transformer$transform, inv = transformer$inverse, n = 4)
    ) +
    scale_y_continuous(
      breaks = trans_breaks(trans = transformer$transform, inv = transformer$inverse, n = 4)
    ) +
    labs(
      x = rep1Col, y = rep2Col,
      title = paste("scatter plot: ", transformer$name, "(", value, "+", pseudoCount, ")", sep = "")
    ) +
    theme_bw() +
    theme_scatter



  ## rank(-value) XY scatter plot
  gg_scatter_rank <- ggplot(
    data = data,
    mapping = aes(x = rank(!!sym(rep1Col)), y = rank(!!sym(rep2Col)), color = density)
  ) +
    geom_point() +
    ggpubr::stat_cor(method = "spearman", size = 6, color = "red",
                     label.x.npc = 0, label.y.npc = 1) +
    viridis::scale_color_viridis(name = "Density", option = colorScale) +
    labs(title = paste("rank scatter plot")) +
    theme_bw() +
    theme_scatter


  ## calculate correlation coefficients on cumulative decreasing data: 100% to 10%
  corDf <- purrr::map_dfr(
    .x = c(seq(0, 0.9, by = 0.1), 0.95, 0.99),
    .f = function(x){
      qtVal <- quantile(data$avgSignal, x)
      tmpExprDf <- data[which(data$avgSignal >= qtVal), ]
      corP <- cor(x = tmpExprDf[[rep1Col]], y = tmpExprDf[[rep2Col]], method = "pearson")
      corS <- cor(x = tmpExprDf[[rep1Col]], y = tmpExprDf[[rep2Col]], method = "spearman")

      return(
        tibble::tibble(
          quant = x, quantPer = scales::percent(x),
          fraction = 1-x, fractionPer = scales::percent(x = fraction, accuracy = 1),
          fpkm = round(qtVal, 2), geneCount = nrow(tmpExprDf),
          pearson = round(x = corP, digits = 3),
          spearman = round(x = corS, digits = 3)
        )
      )
    }
  )


  corPlotDf <- dplyr::mutate(.data = corDf, quantPer = forcats::as_factor(quantPer)) %>%
    tidyr::pivot_longer(cols = c(pearson, spearman),
                        names_to = "method", values_to = "cor")

  corRange <- c(0, 1)
  if(any(corPlotDf$cor < 0)){
    corRange <- c(-1, 1)
  }


  gg_line_cor <- ggplot(data = corPlotDf, mapping = aes(x = quantPer, y = cor, group = method)) +
    geom_line(mapping = aes(linetype = method), size = 1) +
    geom_point(mapping = aes(color = cor), size = 3) +
    # scale_y_continuous(limits = corRange) +
    scale_x_discrete(breaks = corDf$quantPer, labels = corDf$fractionPer) +
    scale_color_gradientn(
      name = "Cor", limits = c(-1, 1), colours = c("red", "yellow", "blue")) +
    labs(x = "top N% genes", y = "Correlation",
         title = "Correlation with subdata") +
    theme_bw() +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14),
      panel.grid = element_blank()
    )

  meanValName = paste("mean(", value, ")", sep = "")
  summaryTable <- purrr::map_dfr(
    .x = structure(c("quantPer", "fpkm", "geneCount"), names = c("Quantile", meanValName, "# of genes")),
    .f = function(x){
      as.list(as.character(corDf[[x]])) %>% purrr::set_names(nm = corDf$fractionPer)
    },
    .id = "top N% genes"
  )


  gg_stable <- ggpubr::ggtexttable(summaryTable, rows = NULL,
                                   theme = ttheme("mOrange"))

  summaryFig <- ggpubr::ggarrange(
    ggpubr::ggarrange(
      gg_fpkm_density,
      ggpubr::ggarrange(gg_scatter_val, gg_scatter_rank,
                        nrow = 1, ncol = 2, legend = "bottom", common.legend = TRUE),
      nrow = 1, widths = c(1, 2)
    ),

    ggpubr::ggarrange(gg_line_cor, gg_stable,
                      nrow = 2, heights = c(5, 2)),
    nrow = 2,
    heights = c(5, 5)
  ) +
    theme(
      plot.margin = unit(rep(0.5, 4), "cm")
    )


  summaryFig <- ggpubr::annotate_figure(
    p = summaryFig,
    top = text_grob(label = paste("Replicates", value, "correlation:", rep1Col, "vs", rep2Col),
                    size = 15, face = "bold")
  )


  plotList <- list(table = gg_stable,
                   density = gg_fpkm_density,
                   scatter = list(value = gg_scatter_val,
                                  rank = gg_scatter_rank),
                   corrVariation = gg_line_cor)

  return(list(data = corDf, figure = summaryFig, plots = plotList, summaryTable = summaryTable))

}




##################################################################################








