
#' Gene length annotation as points
#'
#' This function generates a row annotation of gene length using row_anno_points() of
#' ComplexHeatmap.
#'
#' @param bedFile a bed file with 6 columns. feature name will be used to select the genes
#' @param genes a vector with genes for which annotation needs to be generated. IMP: gene
#' order in this vector should be same as the gene/row order of the matrix used for the main
#' heatmap.
#'
#' @return An anno_points heatmap annotation object
#' @export
#'
#' @examples NA
gene_length_heatmap_annotation = function(bedFile, genes){

  ## read the bed file
  geneSet = data.table::fread(file = bedFile, header = F,
                              col.names = c("chr", "start", "end", "name", "score", "strand")) %>%
    dplyr::mutate(length = end - start) %>%
    dplyr::select(name, length)

  allGl = structure(geneSet$length, names = geneSet$name)

  lenLim = round(quantile(allGl, 0.95))
  ## set max gene length for genes with length > 95% quantile = 95% quantile
  allGl[allGl > lenLim] = lenLim

  selectGl = allGl[genes]

  # ## doing left join on the genes dataframe as gene order needs to be maintained
  # df = data.frame(gene = genes, stringsAsFactors = F) %>%
  #   dplyr::left_join(y = geneSet, by = c("gene" = "name"))


  an = rowAnnotation(gene_length = row_anno_points(x = selectGl,
                                                   size = unit(1, "mm"),
                                                   gp = gpar(col = "#00000040")
  ),
  annotation_width = unit(2, "cm")
  )

  return(list(
    an = an,
    maxLen = lenLim
  ))

}


##################################################################################

