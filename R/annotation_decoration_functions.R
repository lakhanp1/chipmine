

#' Add annotation column titles on the heatmap
#'
#' This function use decorate_annotation function internally to add the title at the top
#' for each annotation
#'
#' @param annotations a character vector of annotation IDs for which lables need to be added
#' @param anTitle a named list with titles for each annotation
#' @param fontSize fontsize for the annotation text
#'
#' @return The function returns no value.
#' @export
#'
#' @examples NA
add_annotation_titles = function(annotations, anTitle, fontSize = 7){
  ## decorate the annotations
  for(an in annotations){

    cat("Adding heatmap annotation title for annotation ", an, "\n")

    decorate_annotation(annotation = an,
                        slice = 1,
                        code = {
                          grid.text(label = anTitle[[an]], x = unit(0.5, "npc"), y = unit(1, "npc") + unit(2, "mm"),
                                    default.units = "npc", just = "left", rot = 90,
                                    gp = gpar(fontsize = fontSize)
                          )
                        }
    )
  }

}


##################################################################################

#' Draw scale for heatmap annotation
#'
#' @param an name of the row annotation which needs to be modified
#' @param at a numeric vector x coordinates where the ticks and labels will be drawn. Eg.: c(0, 2000, 4000)
#' @param labels a vector of lables. Eg.: c("0kb", "2kb", ">4kb")
#' @param slice row slice
#'
#' @return The function returns no value.
#' @export
#'
#' @examples row_annotation_axis(c(0, 2000, 4000), c("0kb", "2kb", ">4kb"), "gene_length", 7)
row_annotation_axis = function(an, at, labels, slice){

  decorate_annotation(an, slice = slice, {
    grid.segments(x0 = at, y0 = unit(0, "npc"), x1 = at, y1 = unit(-1, "mm"),
                  default.units = "native")

    for (i in 1:length(labels)) {
      grid.text(label = labels[i], x = unit(at[i], "native") - unit(2, "mm"), y = unit(-2, "mm"),
                gp = gpar(fontsize = 8), just = c("left", "top"))

    }

  })

}


##################################################################################



