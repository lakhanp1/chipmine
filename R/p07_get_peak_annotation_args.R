
##################################################################################
## peak type color and heatmap legend
#' Get peak heatmap annotation
#'
#' @param anName Annotation name
#'
#' @return A list of two elements: peak type color and peak type legend
#' @export
#'
#' @examples NA
get_peak_annotation_args = function(anName){
  ## peak types
  peakTypeCol = c("include_tx" = "#e41a1c", "include_CDS" = "#e41a1c",
                "5UTR" = "#377eb8", "CDS_start" = "#377eb8", "tx_start" = "#377eb8",
                "3UTR" = "#4daf4a", "tx_end" = "#4daf4a", "CDS_end" = "#4daf4a",
                "inside_tx" = "#984ea3", "inside_CDS" = "#984ea3",
                "promoter" = "#ff7f00", "upstream" = "#ff7f00",
                "pseudo_upstream" = "black", "pseudo_upstream" = "black")

  ## legend for all the annotations
  peakTypeLgd = Legend(title = "\nPeak type",
                       at = names(peakTypeCol),
                       type = "points",
                       pch = 15,
                       size = unit(5, "mm"),
                       grid_height = unit(5, "mm"),
                       grid_width = unit(5, "mm"),
                       legend_gp = gpar(col = peakTypeCol),
                       title_gp = gpar(fontsize = 12, fontface = "bold"),
                       labels_gp = gpar(fontsize = 8)
  )

  ## draw the legend
  # pushViewport(viewport(width = 0.9, height = 0.9))
  # grid.rect()
  # draw(peakTypeLgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"))
  # popViewport()

  return(list(
    col = list(peakTypeCol) %>% purrr::set_names(nm = anName),
    lgd = peakTypeLgd
  ))

}

##################################################################################


