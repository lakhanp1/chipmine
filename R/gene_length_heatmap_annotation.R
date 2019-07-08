
#' Gene length annotation as points
#'
#' This function generates a row annotation of gene length using
#' \code{ComplexHeatmap::row_anno_points()} function
#'
#' @param bedFile a bed file with 6 columns. feature name will be used to select the genes
#' @param genes a vector with genes for which annotation needs to be generated. IMP: gene
#' order in this vector should be same as the gene/row order of the matrix used for the main
#' heatmap.
#' @param axis_param \code{axis_param} argument of \code{ComplexHeatmap::anno_points()}
#' function.
#'
#' @return An anno_points heatmap annotation object
#' @export
#'
#' @examples NA
gene_length_heatmap_annotation = function(bedFile, genes, axis_param = default_axis_param("row")){

  ## read the bed file
  geneSet = data.table::fread(file = bedFile, header = F,
                              col.names = c("chr", "start", "end", "geneId", "score", "strand")) %>%
    dplyr::mutate(length = end - start) %>%
    dplyr::select(geneId, length)

  allGl = structure(geneSet$length, names = geneSet$geneId)

  lenLim = round(quantile(allGl, 0.95))
  ## set max gene length for genes with length > 95% quantile = 95% quantile
  allGl[allGl > lenLim] = lenLim

  selectGl = allGl[genes]

  # ## doing left join on the genes dataframe as gene order needs to be maintained
  # df = data.frame(gene = genes, stringsAsFactors = F) %>%
  #   dplyr::left_join(y = geneSet, by = c("geneId" = "geneId"))

  an = rowAnnotation(
    name = "gene_length",
    gene_length = row_anno_points(
      x = selectGl,
      size = unit(1, "mm"),
      axis_param = axis_param,
      gp = gpar(col = "#00000040")
    ),
    annotation_width = unit(2, "cm"),
    show_annotation_name = TRUE,
    annotation_name_side = "top",
    annotation_name_rot = 90
  )

  return(list(
    an = an,
    maxLen = lenLim
  ))

}


##################################################################################

