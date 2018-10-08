
## functions for TF - polII exploratory analysis

##################################################################################
## plot polII expression distribution violin plot with respect to TF binding status

#' Title
#'
#' @param df
#' @param tf
#' @param polII
#' @param isExpCol
#' @param isBoundCol
#'
#' @return
#' @export
#'
#' @examples
polIIExp_vs_TFbinding_violin = function(df, tf, polII, isExpCol, isBoundCol){
  # df = clusterData dataframe
  # tf = TF sample name
  # polII = polII sample name
  # isExpCol = column name for top 10% polII expression fraction OR expressed status for each gene
  # isBoundCol = column name for TF binding status information

  ## update the columns to factors for ordering in ggplot
  df[[isExpCol]] = factor(df[[isExpCol]])
  df[[isBoundCol]] = factor(df[[isBoundCol]])

  ## statistics to show as annotation on plot
  expCount = structure(paste(table(df[[isExpCol]]), "genes", sep = " "), names = levels(df[[isExpCol]]))
  expCount["All"] = paste(nrow(df), "genes", sep = " ")
  bindCount = structure(paste(table(df[[isBoundCol]]), "genes", sep = " "), names = levels(df[[isBoundCol]]))


  ## new transformation object to transform polII expression values to log2(val + 0.05)
  ## default log2 transformation gives warnings for the expression values == 0
  ## so created this transformation object
  trans_log2Plus = scales::trans_new(name = "trans_log2Plus",
                             transform = function(x){log2(x+0.01)},
                             inverse = function(x){(2**x)-0.01})


  polIIExpLim = min(df[polII][which(df[isExpCol] == "TRUE"), ])
  cat("PolII expression cutoff for", polII, "=", polIIExpLim, "\n")

  allGenes = rep("All", nrow(df))


  ############
  ## distribution of polII expression vs ispolII_expression(top 10%)
  title1 = paste("Distribution of", polII, "expression \n in quantiles 0 - 0.9 and 0.9 - 1", sep = " ")

  p1 = ggplot2::ggplot(data = df, mapping = aes_string(x = as.name(isExpCol), y = polII)) +
    geom_violin(mapping = aes_string(fill = as.name(isExpCol)),
                draw_quantiles = c(0.1, 0.9)) +
    geom_boxplot(width=0.1) +
    geom_violin(mapping = aes_string(x = "allGenes", y = polII),
                draw_quantiles = c(0.1, 0.9)) +
    geom_boxplot(mapping = aes_string(x = "allGenes", y = polII),
                 width = 0.1) +
    geom_hline(yintercept = polIIExpLim,
               color = "blue", linetype="dashed", size = 1) +
    annotate(geom = "text", x = 0.4, y = polIIExpLim+2, label = "polII expressed genes", angle = 90, hjust = 0, vjust = 1) +
    annotate(geom = "text", x = 1, y = 0, label = expCount["All"], hjust = 0.5, vjust = 1.5) +
    annotate(geom = "text", x = 2, y = 0, label = expCount["FALSE"], hjust = 0.5, vjust = 1.5) +
    annotate(geom = "text", x = 3, y = 0, label = expCount["TRUE"], hjust = 0.5, vjust = 1.5) +
    scale_fill_manual(name = paste(polII, "expressed genes", sep = "\n"),
                      values = c("TRUE" = "green4", "FALSE" = "red4"),
                      breaks = c("TRUE", "FALSE"),
                      labels = c("Yes", "No")
    ) +
    scale_y_continuous(name = "polII expression",
                       trans = trans_log2Plus,
                       breaks = c(0, 2 ** c(0:7)),
                       limits = c(0.01, 130),
                       oob = squish
    ) +
    scale_x_discrete(name = "Expression quantile",
                     breaks = c("All","TRUE", "FALSE"),
                     labels = c("All genes", "90 to 100%", "0 to 90%")) +
    ggtitle(title1) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text = element_text(size = 12),
          axis.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12, face = "bold"),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )



  ##############
  ## distribution of polII expression vs TF peak at TSS
  title2 = paste("Distribution of", polII, "expression \n with respect to", tf, "binding near TSS", sep = " ")

  p2 = ggplot2::ggplot(data = df, mapping = aes_string(x = as.name(isBoundCol), y = polII)) +
    geom_violin(mapping = aes_string(fill = as.name(isBoundCol)),
                draw_quantiles = c(0.1, 0.9)) +
    geom_boxplot(width=0.1) +
    geom_hline(yintercept = polIIExpLim,
               color = "blue", linetype="dashed", size = 1) +
    annotate(geom = "text", x = 0.4, y = polIIExpLim+2, label = "polII expressed genes", angle = 90, hjust = 0, vjust = 1) +
    annotate(geom = "text", x = 1, y = 0, label = bindCount["FALSE"], hjust = 0.5, vjust = 1.5) +
    annotate(geom = "text", x = 2, y = 0, label = bindCount["TRUE"], hjust = 0.5, vjust = 1.5) +
    scale_fill_manual(name = paste(tf, "binding near TSS", sep = "\n"),
                      values = c("TRUE" = "green4", "FALSE" = "red4"),
                      breaks = c("TRUE", "FALSE"),
                      labels = c("binding", "No binding")
    ) +
    scale_y_continuous(name = "polII expression",
                       trans = trans_log2Plus,
                       breaks = c(0, 2 ** c(0:7)),
                       limits = c(0.01, 130),
                       oob = squish) +
    scale_x_discrete(name = "TF binding near TSS",
                     breaks = c("TRUE", "FALSE"),
                     labels = c("binding", "no binding")) +
    ggtitle(title2) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text = element_text(size = 12),
          axis.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12, face = "bold"),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )


  #################

  ## get the summary table to show the counts
  summaryTab = df %>%
    dplyr::group_by(!!as.name(isExpCol), !!as.name(isBoundCol)) %>%
    dplyr::summarise(geneCount = n(), median = median(!!as.name(polII)), mean = mean(!!as.name(polII)), stdDev = sd(!!as.name(polII))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate_if(.predicate = is.double,
              .funs = funs(as.numeric(sprintf("%.3f", .)))
    ) %>%
    dplyr::mutate_at(.vars = c(isExpCol, isBoundCol),
              .funs = funs(as.character(.)))


  stable <- ggtexttable(summaryTab, rows = NULL,
                        theme = ttheme(base_style = "mOrangeWhite", base_size = 9))

  #################


  text = paste("\n\n\n\n", isExpCol,
               ": Genes with pol-II expression values in quantil 90-100% \n",
               isBoundCol,
               ": Genes with", tf, "peak near TSS (-500 to TSS)",
               sep = " ")

  text.p <- ggparagraph(text = text, size = 9, color = "black")

  #################

  r1 = ggarrange(p1, p2,
                 ncol = 2, nrow = 1,
                 labels = c("A", "B"),
                 align = "v")

  r2 = ggarrange(stable, text.p,
                 nrow = 1, ncol = 2,
                 labels = "C",
                 widths = c(2.5, 1), align = "v")

  plt = ggarrange(r1, r2,
                  ncol = 1, nrow = 2,
                  heights = c(3,1)
  )

  return(plt)

}

##################################################################################




##################################################################################
## scatter plot for polII expression and TF binding level at promoter/peak region
#' Title
#'
#' @param df
#' @param tf
#' @param polII
#' @param tfExp
#' @param polIIExp
#' @param isExpCol
#' @param isBoundCol
#'
#' @return
#' @export
#'
#' @examples
polII_tf_expression_scatter = function(df, tf, polII, tfExp, polIIExp, isExpCol, isBoundCol){

  df2 = dplyr::mutate(
    .data = df,
    tfPolII = dplyr::if_else(
      UQ(as.name(isExpCol)) == "TRUE",
      dplyr::if_else(UQ(as.name(isBoundCol)) == "TRUE", "Expressed_Bound", "Expressed_notBound"),
      dplyr::if_else(UQ(as.name(isBoundCol)) == "TRUE", "notExpressed_Bound", "notExpressed_notBound")
    )
  )


  ## new transformation object to transform polII expression values to log2(val + 0.05)
  ## default log2 transformation gives warnings for the expression values == 0
  ## so created this transformation object
  trans_log2Plus = scales::trans_new(name = "trans_log2Plus",
                             transform = function(x){log2(x+0.01)},
                             inverse = function(x){(2**x)-0.01})


  title1 = paste("Scatter plot of polII expression distribution vs TF level at promoter/peak region\n",
                 "polII: ", polII, "; TF: ", tf, sep = "")

  plt = ggplot(data = df2) +
    geom_point(mapping = aes_string(x = as.name(tfExp),
                                    y = as.name(polIIExp),
                                    color = "tfPolII"
    )
    ) +
    scale_x_continuous(name = "TF expression at peak/promoter region",
                       trans = trans_log2Plus,
                       breaks = c(0, 2 ** c(0:7)),
                       limits = c(0.1, 130),
                       oob = squish) +
    scale_y_continuous(name = "pol-II expression",
                       trans = trans_log2Plus,
                       breaks = c(0, 2 ** c(0:7)),
                       limits = c(0.1, 130),
                       oob = squish) +
    scale_color_manual(name = "PolII expression and TF binding",
                       values = brewer.pal(n = 4, name = "Set1"),
                       breaks = c("notExpressed_notBound", "notExpressed_Bound", "Expressed_Bound", "Expressed_notBound"),
                       labels = c("not Expressed and no TF binding",
                                  "not expressed and TF binding",
                                  "Expressed and TF binding",
                                  "Expressed and no TF binding")
    ) +
    ggtitle(title1) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text = element_text(size = 12),
          axis.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12, face = "bold"),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    ) +
    guides(color = guide_legend(nrow=2, byrow=TRUE))

}

##################################################################################




##################################################################################
## finds the overlap between the target genes of two TFs, draw Venn diagram and
## run topGO enrichment for each partition of the venn diagram
#' Title
#'
#' @param df
#' @param tf1
#' @param tf2
#' @param outPrefix
#' @param topGoMapFile
#'
#' @return
#' @export
#'
#' @examples
tf_binding_overlap = function(df, tf1, tf2, outPrefix, topGoMapFile){
  hasPeakTf1 = paste("hasPeak.", tf1, sep = "")
  hasPeakTf2 = paste("hasPeak.", tf2, sep = "")

  ## select the genes which have binding in both the samples
  hasPeakDf = dplyr::filter_at(.tbl = df, .vars = c(hasPeakTf1, hasPeakTf2), .vars_predicate = any_vars(. == "TRUE")) %>%
    dplyr::mutate(
      overlap = dplyr::if_else(
        UQ(as.name(hasPeakTf1)) == "TRUE" & UQ(as.name(hasPeakTf2)) == "TRUE",
        "common",
        dplyr::if_else(UQ(as.name(hasPeakTf1)) == "TRUE", tf1, tf2)))

  geneList = list()

  geneList[[tf1]] = hasPeakDf$gene[ hasPeakDf[[hasPeakTf1]] ]
  geneList[[tf2]] = hasPeakDf$gene[ hasPeakDf[[hasPeakTf2]] ]

  ## plot Venn diagram for variants
  vennOut = paste(outPrefix, "_venn.png",  sep = "")
  vennTitle = paste("TF binding comparison between", tf1, "and", tf2, sep = " ")

  vennAll = venn.diagram(x = geneList,
                         category.names = names(geneList),
                         main = vennTitle,
                         imagetype = "png",
                         filename = vennOut,
                         height = 3000,
                         width = 3000,
                         resolution = 500,
                         euler.d = F,
                         scaled = F,
                         fill = c("blue", "green"),
                         alpha = 0.5,
                         fontface = "bold",
                         main.cex = 1,
                         cex = 1,
                         cat.fontface = "bold",
                         cat.cex = 1,
                         # cat.dist = c(0.05, 0.05),
                         # cat.pos = c(0, 0),
                         cat.default.pos = "outer",
                         # cat.just = list(c(0,1), c(0,1)),
                         margin = 0.2
  )

  unlink(paste(vennOut, ".*.log", sep = ""))


  ## tf1 specific genes
  tf1Specific = hasPeakDf$gene[ which(hasPeakDf$overlap == tf1) ]

  goTitle = paste("GO enrichment for genes with only", tf1, "peak:", length(tf1Specific), "genes", sep = " ")
  tf1GOOut = paste(outPrefix, "_", tf1, "specific_GO.png", sep = "")

  tf1Go = go_and_scatterPlot(goToGeneFile = topGoMapFile,
                             genes = tf1Specific,
                             goTitle = goTitle,
                             plotOut = tf1GOOut)

  tf1Go$vennPartition = tf1


  ## tf2 specific genes
  tf2Specific = hasPeakDf$gene[ which(hasPeakDf$overlap == tf2) ]

  goTitle = paste("GO enrichment for genes with only", tf2, "peak:", length(tf2Specific), "genes", sep = " ")
  tf2GoOut = paste(outPrefix, "_", tf2, "specific_GO.png", sep = "")

  tf2Go = go_and_scatterPlot(goToGeneFile = topGoMapFile,
                             genes = tf2Specific,
                             goTitle = goTitle,
                             plotOut = tf2GoOut)

  tf2Go$vennPartition = tf2


  ## gene targets common between TF1 and TF2
  tfCommon = hasPeakDf$gene[ which(hasPeakDf$overlap == "common") ]

  goTitle = paste("GO enrichment for common target genes of", tf1, "and", tf2, ":", length(tfCommon), "genes", sep = " ")
  commonGoOut = paste(outPrefix, "_common_GO.png", sep = "")

  commonGo = go_and_scatterPlot(goToGeneFile = topGoMapFile,
                                genes = tfCommon,
                                goTitle = goTitle,
                                plotOut = commonGoOut)

  commonGo$vennPartition = "common"


  goTable = bind_rows(tf1Go, tf2Go, commonGo)

  fwrite(x = goTable, file = paste(outPrefix, "_binding_GO.tab", sep = ""), sep = "\t", col.names = T, quote = F)

  return(goTable)

}

##################################################################################





