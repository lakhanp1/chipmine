
##################################################################################
##
#' Calculate average profile for a gene list
#'
#' This method calculates the average profile for gene of interest.
#'
#' @param exptInfo Experiment information dataframe
#' @param profileMats a named list of profile matrix.
#' @param genes a vector of gene IDs. Profile matrix will be subset based on this
#' gene list to calculate average profile
#' @param cluster A group name for current profile. Default: group_A
#'
#' @return A dataframe which can be used for plotting average profile
#' @export
#'
#' @examples NULL
geneset_average_profile <- function(exptInfo, profileMats, genes, cluster = "group_A"){

  meanProfile <- list()
  checkList <- c(attributes(profileMats[[1]])$dim, attributes(profileMats[[1]])$extend)

  # i <- 1
  ## calculate mean profile for each matrix
  for (i in 1:nrow(exptInfo)) {
    mat <- profileMats[[ exptInfo$sampleId[i] ]]

    ## make sure that all the matrices have same dimensions and extension
    if(! all(checkList == c(attributes(mat)$dim, attributes(mat)$extend))){
      stop("Error: matrix are of not same dimentions")
    }

    mat <- mat[genes, ]

    sdCol <- paste(exptInfo$sampleId[i], "_sd", sep = "")

    meanProfile[[ exptInfo$sampleId[i] ]] <- apply(mat, MARGIN = 2, mean)
    meanProfile[[ sdCol ]] <- apply(mat, MARGIN = 2, sd)

    names(meanProfile[[ exptInfo$sampleId[i] ]]) <- c(attributes(mat)$upstream_index,
                                                      attributes(mat)$target_index,
                                                      attributes(mat)$downstream_index)

    names(meanProfile[[ sdCol ]]) <- c(attributes(mat)$upstream_index,
                                       attributes(mat)$target_index,
                                       attributes(mat)$downstream_index)

    cat("Calculated mean profile for sample", exptInfo$sampleId[i], "\n")
  }

  meanProfileDf <- data.table::as.data.table(meanProfile)

  meanProfileDf$bin <- c(attributes(profileMats[[1]])$upstream_index,
                         attributes(profileMats[[1]])$target_index,
                         attributes(profileMats[[1]])$downstream_index)


  grp1Cols <- exptInfo$sampleId
  grp2Cols <- paste(exptInfo$sampleId, "_sd", sep = "")

  plotDf <- data.table::melt(data = meanProfileDf,
                            id.vars = c("bin"),
                            measure.vars = list(grp1Cols, grp2Cols),
                            value.name = c("mean", "sd"),
                            variable.name = "sample"
  ) %>%
    as.data.frame()

  levels(plotDf$sample) <- exptInfo$sampleId
  # plotDf$sample <- as.character(plotDf$sample)

  plotDf$cluster <- cluster
  plotDf$groupSize <- length(genes)
  plotDf$groupLabel <- paste(cluster, ":", length(genes), "genes")

  return(plotDf)
}

##################################################################################



##################################################################################
##
#' Average signal plot for multiple samples
#'
#' @param exptInfo Experiment information dataframe
#' @param profileMats a named list of profile matrix
#' @param genes a vector of gene IDs. Profile matrix will be subset based on this
#' gene list to calculate average profile
#' @param lineColors color for each sample
#' @param lineShape linetype for each sample
#'
#' @return a ggplot object with line plot
#' @export
#'
#' @examples NA
draw_avg_profile_plot <- function(exptInfo,
                                  profileMats,
                                  genes,
                                  lineColors = NULL,
                                  lineShape = NULL){


  plotDf <- geneset_average_profile(exptInfo = exptInfo, profileMats = profileMats, genes = genes)

  ## decide the axis labels
  axisBrk <- NULL
  axisLab <- NULL
  if(attributes(profileMats[[1]])$target_is_single_point){
    axisBrk <- c(
      attributes(profileMats[[1]])$upstream_index[1],
      attributes(profileMats[[1]])$target_index[1],
      tail(attributes(profileMats[[1]])$downstream_index, 1)
    )

    axisLab <- c(
      -attributes(profileMats[[1]])$extend[1],
      attributes(profileMats[[1]])$target_name,
      attributes(profileMats[[1]])$extend[2]
    )

    targetEnd <- axisBrk[2] + 1

  } else if(! attributes(profileMats[[1]])$target_is_single_point){
    axisBrk <- c(
      attributes(profileMats[[1]])$upstream_index[1],
      attributes(profileMats[[1]])$target_index[1],
      tail(attributes(profileMats[[1]])$target_index, 1),
      tail(attributes(profileMats[[1]])$downstream_index, 1)
    )

    axisLab <- c(
      -attributes(profileMats[[1]])$extend[1],
      "START", "END",
      attributes(profileMats[[1]])$extend[2]
    )

    targetEnd <- axisBrk[3]

  }

  p <- ggplot2::ggplot(data = plotDf) +
    geom_line(mapping = aes(x = bin, y = mean, group = sample, color = sample, linetype = sample),
              size = 0.8, alpha = 0.8) +
    geom_hline(yintercept = 0, size = 2, color = "grey70") +
    geom_segment(mapping = aes(x = axisBrk[2], y = 0, xend = targetEnd, yend = 0),
                 size = 10, lineend = "butt", color = "grey50") +
    # annotate(geom = "text", x = 230, y = 0, label = "Scaled gene body", hjust = 0, vjust = 0.3) +
    scale_color_manual(values = lineColors) +
    scale_linetype_manual(values = lineShape) +
    scale_x_continuous(breaks = axisBrk, labels = axisLab) +
    ylab("Read coverage") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 13, angle = 90),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold"),
          legend.position = c(1, 1),
          legend.justification = c(1.1, 1.1),
          legend.text = element_text(size = 13),
          legend.title = element_blank()
    )

  # p

  return(p)

}




##################################################################################


