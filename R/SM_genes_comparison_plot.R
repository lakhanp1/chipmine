

##################################################################################
## generate the multi sample comparison plot for SM genes
#' Title
#'
#' @param data
#' @param orderFile
#' @param exptInfo
#' @param outFile
#' @param name
#' @param res
#'
#' @return
#' @export
#'
#' @examples NA
SM_genes_comparison_plot = function(data, orderFile, exptInfo, outFile, name, res){

  plotTitle = paste("Secondary metabolite genes for comparison", name, "\n")

  ## prepare SM dataframe
  smData = data %>% filter(is_SM_gene == TRUE)

  clusterOrder = fread(input = orderFile, header = T, stringsAsFactors = F, data.table = F)
  smData$SM_cluster = factor(smData$SM_cluster, levels = clusterOrder$clusters)
  smData = smData[order(smData$SM_cluster), ]

  row.names(smData) = smData$Gene


  htSplit = smData["SM_cluster"]


  ## clusterAnnotation
  clusterIds = as.character(unique(smData$SM_cluster))
  cluster_colors = structure(rep(x = c("lightgrey", "darkgrey"), length.out = length(clusterIds)), names = clusterIds)

  clustAn = rowAnnotation(df = smData["SM_cluster"],
                          col = list(SM_cluster = cluster_colors),
                          width = unit(0.5, "cm"),
                          show_legend = FALSE
  )


  ## TF annotation
  tfAn = rowAnnotation(df = smData["is_TF"],
                       col = list(is_TF = c("TRUE" = "blue", "FALSE" = "grey95")),
                       width = unit(1, "cm"),
                       show_legend = FALSE
  )


  ## plot polII data
  polII_ids = exptInfo$sampleId[which(exptInfo$IP_tag == "polII")]

  polII_mat = as.matrix(smData[polII_ids])

  polII_log2_mat = log2(polII_mat + 1)

  quantile(polII_log2_mat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 1), na.rm = T)


  polII_color = colorRamp2(breaks = seq(quantile(polII_log2_mat, 0.5), quantile(polII_log2_mat, 0.995), length.out = 9), colors = brewer.pal(9, "Reds"))


  polII_ht = Heatmap(matrix = polII_log2_mat,
                     name = "log2(polII_expression)",
                     col = polII_color,
                     split = htSplit,
                     cluster_rows = TRUE,
                     cluster_columns = FALSE,
                     show_row_names = TRUE,
                     column_names_gp = gpar(fontsize = 12),
                     row_names_side = "right",
                     row_names_gp = gpar(fontsize = 3),
                     row_title_rot = 0,
                     row_title_gp = gpar(fontsize = 10),
                     show_heatmap_legend = FALSE,
                     heatmap_legend_param = list(title = "log2(polII_expression)\n",
                                                 title_gp = gpar(fontsize = 14, fontface = "bold"),
                                                 labels_gp = gpar(fontsize = 12),
                                                 legend_height = unit(4, "cm")
                                                 # legend_width = unit(2, "cm")
                     )
  )


  htLgd = Legend(title = "log2(polII_expression)",
                 at = as.numeric(sprintf("%.0f", attributes(polII_color)$breaks)),
                 # labels = sprintf("%.0f", attributes(polII_color)$breaks),
                 col_fun = polII_color,
                 type = "grid",
                 direction = "horizontal",
                 title_gp = gpar(fontsize = 12, fontface = "bold"),
                 labels_gp = gpar(fontsize = 8),
                 legend_width = unit(4, "cm"),
                 grid_height = unit(5, "mm")
                 # grid_width = unit(5, "mm")
  )


  ## is_expressed(polII_sample) annotations
  isExpressedDf = smData[paste("is_expressed.", polII_ids, sep = "")]

  expressColor = sapply(X = names(isExpressedDf),
                        FUN = function(x){structure(c("green4", "grey95"), names = c("TRUE", "FALSE"))},
                        simplify = F,
                        USE.NAMES = T
  )

  expressedAnn = HeatmapAnnotation(df = isExpressedDf,
                                   which = "row",
                                   col = expressColor,
                                   show_legend = F,
                                   annotation_width = unit(rep(1, times = length(names(isExpressedDf))), "cm"),
                                   gap = unit(3, "mm")
  )


  ## hasPeak.TF_sample annotations
  tfSampleIds = exptInfo$sampleId[which(exptInfo$IP_tag != "polII")]

  hasPeakDf = smData[paste("hasPeak.", tfSampleIds, sep = "")]

  tfColor = sapply(X = names(hasPeakDf),
                   FUN = function(x){structure(c("sienna", "grey95"), names = c("TRUE", "FALSE"))},
                   simplify = F,
                   USE.NAMES = T
  )

  hasPeakAnn = HeatmapAnnotation(df = hasPeakDf,
                                 which = "row",
                                 col = tfColor,
                                 show_legend = F,
                                 annotation_width = unit(rep(1, times = length(names(hasPeakDf))), "cm"),
                                 gap = unit(3, "mm")
  )


  annLgd = Legend(title = "\nGene categories",
                  at = c("Transcription Factor", "Expressed gene", "TF binding peak at TSS"),
                  type = "points",
                  pch = 15,
                  size = unit(5, "mm"),
                  grid_height = unit(5, "mm"),
                  grid_width = unit(5, "mm"),
                  legend_gp = gpar(col = c("blue", "green4", "sienna")),
                  title_gp = gpar(fontsize = 12, fontface = "bold"),
                  labels_gp = gpar(fontsize = 8)
  )



  ## heatmap drawing
  htList = clustAn + tfAn + polII_ht + expressedAnn + hasPeakAnn

  imageWidth = (length(polII_ids) * 1200) + (length(tfSampleIds) * 300) + 800
  imageHeight = 10000
  imageRes = imageWidth/10

  if(!missing(res)){
    imageRes = res
  }

  cat("Image width =", imageWidth, "\n")
  cat("Image height =", imageHeight, "\n")
  cat("Image resulution =", imageRes, "\n")

  png(filename = outFile, width = imageWidth, height = imageHeight, res = imageRes)

  ## draw Heatmap list
  draw(htList,
       column_title = plotTitle,
       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
       padding = unit(c(3, 1, 1, 1), "cm"),
       row_sub_title_side = "left",
       annotation_legend_list = list(htLgd, annLgd)
  )


  ## SM_cluster annotation name
  decorate_annotation(annotation = "SM_cluster",
                      code = {
                        grid.text(label = "SM clusters", x = unit(0.5, "npc"), y = unit(0, "npc") - unit(2, "mm") ,
                                  default.units = "npc", just = "right", rot = 90,
                                  gp = gpar(fontsize = 10))
                      },
                      slice = length(clusterIds)
  )


  ## is_TF annotation name
  decorate_annotation(annotation = "is_TF",
                      code = {
                        grid.text(label = "Transcription Factor",
                                  x = unit(0.5, "npc"), y = unit(0, "npc") - unit(2, "mm") ,
                                  default.units = "npc", just = "right", rot = 90,
                                  gp = gpar(fontsize = 10))
                      },
                      slice = length(clusterIds)
  )


  ## isExpressed annotation names
  for(an in colnames(isExpressedDf)){
    txt = gsub("is_expressed", "is expressed\n", an) %>% gsub("\\(|\\)", "", .)

    decorate_annotation(annotation = an,
                        code = {
                          grid.text(label = txt, x = unit(0.5, "npc"), y = unit(0, "npc") - unit(2, "mm"),
                                    default.units = "npc", just = "right", rot = 90,
                                    gp = gpar(fontsize = 10))
                        },
                        slice = length(clusterIds)
    )
  }


  ## hasPeak annotation names
  for(an in colnames(hasPeakDf)){
    txt = gsub("hasPeak.", "has TSS peak\n", an)

    decorate_annotation(annotation = an,
                        code = {
                          grid.text(label = txt, x = unit(0.5, "npc"), y = unit(0, "npc") - unit(2, "mm"),
                                    default.units = "npc", just = "right", rot = 90,
                                    gp = gpar(fontsize = 10))
                        },
                        slice = length(clusterIds)
    )
  }


  dev.off()


}


##################################################################################
