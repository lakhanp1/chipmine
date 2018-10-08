
##################################################################################

## for a set of genes, generate a multi TF heatmap
#' Title
#'
#' @param data
#' @param exptInfo
#' @param outFile
#' @param plotTitle
#' @param name
#' @param res
#' @param numClust
#'
#' @return
#' @export
#'
#' @examples NA
geneset_multiTF_plot = function(data, exptInfo, outFile, plotTitle, name, res, numClust){

  ## empty heatmap list
  htList = NULL
  lgdList = list()


  polII_ids = exptInfo$sampleId[which(exptInfo$IP_tag == "polII")]
  polII_expIds = paste("is_expressed.", polII_ids, sep = "")

  ## is_expressed(polII_sample) annotations
  expressColor = sapply(X = polII_expIds,
                        FUN = function(x){structure(c("green4", "grey95"), names = c("TRUE", "FALSE"))},
                        simplify = F,
                        USE.NAMES = T
  )

  expressedAn = HeatmapAnnotation(df = data[polII_expIds],
                                  which = "row",
                                  col = expressColor,
                                  show_legend = FALSE,
                                  annotation_width = unit(rep(1, times = length(polII_expIds)), "cm"),
                                  gap = unit(3, "mm")
  )



  ## SM gene annotation
  smGeneAn = HeatmapAnnotation(df = data["is_SM_gene"],
                               which = "row",
                               col = list(is_SM_gene = c("TRUE" = "orange", "FALSE" = "grey95")),
                               annotation_width = unit(1, "cm"),
                               show_legend = FALSE
  )


  ## TF annotation
  tfAn = rowAnnotation(df = data["is_TF"],
                       col = list(is_TF = c("TRUE" = "blue", "FALSE" = "grey95")),
                       annotation_width = unit(1, "cm"),
                       show_legend = FALSE
  )



  ## clustering polII expression data
  polII_mat = as.matrix(data[polII_ids])

  polII_log2_mat = log2(polII_mat + 1)
  dend = hclust(dist(polII_log2_mat), method = "ward.D")


  ## polII Heatmap
  quantile(polII_log2_mat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 1), na.rm = T)


  polII_color = colorRamp2(breaks = seq(quantile(polII_log2_mat, 0), quantile(polII_log2_mat, 0.99), length.out = 9), colors = brewer.pal(9, "Reds"))
  # colorRamp2(breaks = quantile(polII_log2_mat, 0, 0.99),colors =  c("white", "red"))

  polII_ht = Heatmap(matrix = polII_log2_mat,
                     name = name,
                     col = polII_color,
                     cluster_rows = dend,
                     row_dend_width = unit(3, "cm"),
                     row_dend_reorder = TRUE,
                     row_dend_side = "left",
                     cluster_columns = FALSE,
                     show_row_names = FALSE,
                     column_names_gp = gpar(fontsize = 10),
                     row_title_rot = 0,
                     row_title_gp = gpar(fontsize = 12),
                     show_heatmap_legend = FALSE,
                     heatmap_legend_param = list(title = "log2(polII_expression)",
                                                 title_gp = gpar(fontsize = 12, fontface = "bold"),
                                                 labels_gp = gpar(fontsize = 8),
                                                 color_bar = "continuous",
                                                 legend_direction = "horizontal",
                                                 legend_width = unit(4, "cm"),
                                                 # legend_height = unit(4, "cm"),
                                                 grid_height = unit(5, "mm")
                                                 # grid_width = unit(5, "mm")
                     ),
                     width = 1 * ncol(polII_log2_mat)
  )

  ## heatmap legend
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


  ## hasPeak.TF_sample annotations
  tfSampleIds = exptInfo$sampleId[which(exptInfo$IP_tag != "polII")]

  hasPeakDf = data[paste("hasPeak.", tfSampleIds, sep = "")]

  tfColor = sapply(X = names(hasPeakDf),
                   FUN = function(x){structure(c("sienna", "grey95"), names = c("TRUE", "FALSE"))},
                   simplify = F,
                   USE.NAMES = T
  )

  hasPeakAn = HeatmapAnnotation(df = hasPeakDf,
                                which = "row",
                                col = tfColor,
                                show_legend = F,
                                annotation_width = unit(rep(1, times = length(names(hasPeakDf))), "cm"),
                                gap = unit(3, "mm")
  )


  ## legend for all the annotations
  annLgd = Legend(title = "\nGene categories",
                  at = c("SM gene", "Transcription Factor", "Expressed gene", "TF binding peak at TSS"),
                  type = "points",
                  pch = 15,
                  size = unit(5, "mm"),
                  grid_height = unit(5, "mm"),
                  grid_width = unit(5, "mm"),
                  legend_gp = gpar(col = c("orange", "blue", "green4", "sienna")),
                  title_gp = gpar(fontsize = 12, fontface = "bold"),
                  labels_gp = gpar(fontsize = 8)
  )


  ## cut dendrogram and draw cluster annotation only if number of cuts given
  clustLgd = NULL
  clusterData = NULL

  if(!missing(numClust)){

    clusterCut = cutree(dend, numClust)
    clusterData = data.frame(gene = names(clusterCut),
                             cluster = paste("cluster", clusterCut, sep = "_"),
                             stringsAsFactors = F)

    clusterColor <- structure(brewer.pal(9, "Set1")[1:numClust], names = unique(clusterData$cluster))

    # cluster annotation
    clustAn <- HeatmapAnnotation(df = clusterData["cluster"],
                                 which = "row",
                                 width = unit(2, "mm"),
                                 col = list(cluster = clusterColor),
                                 show_legend = FALSE
    )

    htList = htList + clustAn

    clustLgd = Legend(title = "Gene clusters",
                      at = names(clusterColor),
                      type = "points",
                      pch = 15,
                      size = unit(5, "mm"),
                      grid_height = unit(5, "mm"),
                      grid_width = unit(5, "mm"),
                      legend_gp = gpar(col = clusterColor),
                      title_gp = gpar(fontsize = 12, fontface = "bold"),
                      labels_gp = gpar(fontsize = 8),
                      ncol = 2
    )

    data = left_join(x = data, y = clusterData, by = c("gene" = "gene"))
  }

  ## all plotting Heatmap and HeatmapAnnotation objects
  htList = htList + smGeneAn + tfAn + polII_ht + expressedAn + hasPeakAn

  lgdList = list(htLgd, annLgd)

  if(!missing(numClust)){
    lgdList = list(clustLgd, htLgd, annLgd)
  }


  ## generate plot

  imageWidth = (length(polII_ids) * 750) + (length(tfSampleIds) * 500) + 1500
  imageHeight = nrow(polII_log2_mat) * 5
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
       column_title_gp = gpar(fontsize = 16, fontface = "bold"),
       padding = unit(c(2, 1, 1, 1), "cm"),
       row_dend_side = "left",
       row_sub_title_side = "left",
       annotation_legend_list = lgdList
  )


  ## SM_cluster annotation name
  decorate_annotation(annotation = "is_SM_gene",
                      code = {
                        grid.text(label = "SM gene", x = unit(0.5, "npc"), y = unit(0, "npc") - unit(2, "mm") ,
                                  default.units = "npc", just = "right", rot = 90,
                                  gp = gpar(fontsize = 10)
                        )
                      }
  )


  ## is_TF annotation name
  decorate_annotation(annotation = "is_TF",
                      code = {
                        grid.text(label = "Transcription Factor",
                                  x = unit(0.5, "npc"), y = unit(0, "npc") - unit(2, "mm") ,
                                  default.units = "npc", just = "right", rot = 90,
                                  gp = gpar(fontsize = 10)
                        )
                      }
  )


  ## isExpressed annotation names
  for(an in polII_expIds){
    txt = gsub("is_expressed", "is expressed\n", an) %>% gsub("\\(|\\)", "", .)

    decorate_annotation(annotation = an,
                        code = {
                          grid.text(label = txt, x = unit(0.5, "npc"), y = unit(0, "npc") - unit(2, "mm"),
                                    default.units = "npc", just = "right", rot = 90,
                                    gp = gpar(fontsize = 10)
                          )
                        }
    )
  }


  ## hasPeak annotation names
  for(an in colnames(hasPeakDf)){
    txt = gsub("hasPeak.", "has TSS peak\n", an)

    decorate_annotation(annotation = an,
                        code = {
                          grid.text(label = txt, x = unit(0.5, "npc"), y = unit(0, "npc") - unit(2, "mm"),
                                    default.units = "npc", just = "right", rot = 90,
                                    gp = gpar(fontsize = 10)
                          )
                        }
    )
  }


  dev.off()

  return(list(
    "df" = data,
    "dend" = dend
  ))

}

##################################################################################
