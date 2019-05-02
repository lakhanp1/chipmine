

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

    decorate_annotation(
      annotation = an,
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


