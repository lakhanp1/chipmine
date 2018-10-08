
##################################################################################
## peak type color and heatmap legend
#' Get peak heatmap annotation
#'
#' @param tssCol column name of the TSS region peaks
#' @param tesCol column name of the TES region peaks
#'
#' @return A list of two elements: peak type color and peak type legend
#' @export
#'
#' @examples NA
get_peak_annotation_args = function(tssCol, tesCol){
  ## peak types
  peakTypes = c("upstream", "overlapStart", "insideOverlapStart", "includeFeature", "insideOverlapStartOverlapEnd", "inside", "insideOverlapEnd", "overlapEnd", "pseudo_upstream", "pseudo_overlapEnd", "pseudo_insideOverlapEnd", "pseudo_inside")


  peakColors = c(RColorBrewer::brewer.pal(8, "Paired"), "black", "black", "black", "black")

  peakTypeCol = sapply(
    X = c(tssCol, tesCol),
    FUN = function(x){
      structure(.Data = peakColors,
                names = peakTypes)},
    USE.NAMES = T,
    simplify = F
  )

  ## legend for all the annotations
  peakTypeLgd = Legend(title = "\nPeak type",
                       at = peakTypes,
                       type = "points",
                       pch = 15,
                       size = unit(5, "mm"),
                       grid_height = unit(5, "mm"),
                       grid_width = unit(5, "mm"),
                       legend_gp = gpar(col = peakColors),
                       title_gp = gpar(fontsize = 12, fontface = "bold"),
                       labels_gp = gpar(fontsize = 8)
  )

  return(list(
    col = peakTypeCol,
    lgd = peakTypeLgd
  ))

}

##################################################################################


