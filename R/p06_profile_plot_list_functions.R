

##################################################################################

## function to get a heatmap list of provided TF, polII samples to be used directly for plotting
#' \code{EnrichedHeatmap} plots and \code{ComplexHeatmaps} for multiple samples
#'
#' This method generates profile heatmaps for multiple samples. Optionally, if the sample list
#' has any polII samples, it also generate a heatmap for the polII expression values.
#'
#' @param exptInfo experiment info as data frame with information like sampleID, type, path etc
#' @param genesToPlot A vector of gene IDs which are to be plotted
#' @param clusters Dataframe with cluster information. Two columns must be present in this
#' dataframe: 'cluster', 'geneId'. Default: NULL
#' @param drawClusterAn Logical: Whether to add row annotation heatmap to heatmap list. Default: TRUE
#' @param clusterColor cluster annotation color information. Default: NULL
#' @param clustOrd A character vector of cluster order in the plot. If not provided,
#' clusters are arranged as per character sort order
#' @param profileColors a named list of color objects for each of the sample. Color for profile
#'  heatmap should be generated by colorRamp2
#' @param matSource Method by which profile matrix was generated. One of \emph{deeptools, miao,
#' normalizedmatrix}. Default: normalizedmatrix
#' @param targetType One of "region", "point". If targetType is "region", target is used to decide
#' the number of bins over region. Otherwise, for "point", a signle point matrix is expected where
#' the region is extracted around single point. Default: region
#' @param targetName Name to be used for target. Eg: "gene", "TSS", "TES", "summit" etc. default: gene
#' @param matBins A vector with four elements representing the column counts for up, body,
#'  down regions and binSize for the profile matrix. Default: \code{c(200, 200, 100, 10)}
#' @param ylimFraction A named list with intensity scores to use as ylimit for top annotation.
#' If the value is single number, it has to be floating point number to extract the quantile and use
#' limit \code{[0, quantile(x)]}. If the value is numeric vector of length 2, the two numbers are used as lower
#' and upper limits for ylimit of top annotation. Default: NULL.
#' @param plotExpression Binry: whether to plot polII expression heatmap or not. Default: FALSE
#' @param expressionData an dataframe which has info: clustering, polII expression, TF binding
#'  status for each gene etc. Default: NULL
#' @param keep Same as \code{EnrichedHeatmap::normalizeToMatrix}. First value is used as lower quantile
#' and any value in profile matrix less than lower quantile is set to lower quantile. Second value is
#' used as upper quantile and any value greater than upper quantile is set to upper quantile.
#' Default: \code{c(0, 1)}
#' @param expressionColor color object for polII expression heatmap
#' @param rasterize Binary: whether to rasterize the profile heatmap or not. By default this is set as TRUE
#' for multiple profile plots.
#' @param rasterQual Raster quality. This is the \code{raster_quality} argument \code{Heatmap} function. in Default: 5
#' @param ... Other arguments for EnrichedHeatmap function
#'
#' @return A list object with following elements:
#' \itemize{
#' \item heatmapList: a heatmap list object which has profile heatmaps and optionally expression heatmaps
#' \item profileHeatmaps: a named list of profile heatmaps generated by \code{EnrichedHeatmap}
#' \item expressionHeatmaps: a named list of polII signal heatmaps
#' \item plotGaps: a vector with required gaps between heatmaps while generating multiprofile plots
#' \item profileColors: a named list object with color scale used for each profile heatmap
#' \item expressionColor: a named list object with color scale used for each expression heatmap
#' }
#'
#' @export
#'
#' @examples NA
multi_profile_plots <- function(exptInfo,
                                genesToPlot,
                                clusters = NULL,
                                drawClusterAn = TRUE,
                                clusterColor = NULL,
                                clustOrd = NULL,
                                profileColors = NULL,
                                matSource = "normalizedmatrix",
                                targetType = "region",
                                targetName = "gene",
                                matBins = c(200, 200, 100, 10),
                                ylimFraction = NULL,
                                plotExpression = FALSE,
                                expressionData = NULL,
                                expressionColor = NULL,
                                keep = c(0, 1),
                                rasterize = TRUE,
                                rasterQual = 5,
                                ...){


  cat("Generating profile heatmaps...\n")

  ## profile heatmap list for all samples except first
  profileList <- get_profile_plot_list(exptInfo = exptInfo,
                                       cluster = clusters,
                                       clusterColor = clusterColor,
                                       clustOrd = clustOrd,
                                       geneList = genesToPlot,
                                       colorList = profileColors,
                                       matrixSource = matSource,
                                       matrixBins = matBins,
                                       targetType = targetType,
                                       targetName = targetName,
                                       ylimFraction = ylimFraction,
                                       keep = keep,
                                       rasterPar = list(use = rasterize, qual = rasterQual),
                                       ...)

  # cat("Done!!!\n")

  expHtList <- list()
  expHtList[["polII_color"]] <- expressionColor

  if(isTRUE(plotExpression)){

    cat("Generating expression heatmaps...\n")

    ## expression heatmap list
    expHtList <- get_expression_heatmap_list(expDf = expressionData,
                                             exptInfo = exptInfo,
                                             genes = genesToPlot,
                                             htColor = expressionColor)

    # cat("Done!!!\n")
  }

  cat("Generating heatmap list object...\n")

  ## build the heatmap list for plotting
  plotGaps <- numeric()
  htList <- NULL

  if(isTRUE(drawClusterAn)){
    htList <- profileList[[ exptInfo$sampleId[1] ]][[ "rowGroupHt" ]]
    plotGaps <- c(2)
  }

  profileColorList <- list()

  for(i in 1:nrow(exptInfo)){
    sampleID <- exptInfo$sampleId[i]

    ## add plot elements to the heatmap list
    if(exptInfo$IP_tag[i] == "polII"){

      ## polII expression + profile
      if(isTRUE(plotExpression)){
        htList <- htList + expHtList[[sampleID]] + profileList[[sampleID]][["heatmap"]]
        plotGaps <- append(plotGaps, c(2, 10))
      } else{
        ## only polII profile
        htList <- htList + profileList[[sampleID]][["heatmap"]]
        plotGaps <- append(plotGaps, 10)
      }

    } else{
      ## TF profile
      htList <- htList + profileList[[sampleID]][["heatmap"]]
      plotGaps <- append(plotGaps, 10)
    }

    ## add profile colors to profileColorList
    profileColorList[[sampleID]] <- profileList[[sampleID]][["profileColor"]]

  }

  plotGaps <- head(plotGaps, -1)

  returnList <- list(
    "heatmapList" = htList,
    "profileHeatmaps" = profileList,
    "plotGaps" = plotGaps,
    "profileColors" = profileColorList,
    "expressionColor" = expHtList$polII_color
  )

  if(isTRUE(plotExpression)){
    returnList[["expressionHeatmaps"]] <- expHtList
  }

  return(returnList)

}

##################################################################################




##################################################################################
## function to get the profile plots in list() using EnrichedHeatmap package
#' Multiple Enriched Heatmaps as list
#'
#' @param exptInfo experiment info as data frame with information like sampleID, type,
#'  path etc
#' @param cluster Dataframe with cluster information. Two columns must be present
#'  in this dataframe: 'cluster', 'geneId'. Default: NULL
#' @param clusterColor cluster annotation color information. Default: NULL
#' @param clustOrd A character vector of cluster order in the plot. If not provided,
#' clusters are arranged as per character sort order
#' @param geneList A vector of gene IDs which are to be plotted
#' @param colorList a named list of color objects for each of the sample. Color for
#'  profile heatmap should be generated by colorRamp2
#' @param matrixSource Method by which profile matrix was generated. One of \emph{deeptools, miao,
#' normalizedmatrix}. Default: normalizedmatrix
#' @param matrixBins A numeric vector with four elements representing the column
#' counts for up, body, down regions and binSize for the profile matrix.
#' Default: \code{c(200, 200, 100, 10)}
#' @param targetType One of "region", "point". If targetType is "region", target is used to decide
#' the number of bins over region. Otherwise, for "point", a signle point matrix is expected where
#' the region is extracted around single point. Default: region
#' @param targetName Name to be used for target. Eg: "gene", "TSS", "TES", "summit" etc. default: gene
#' @param ylimFraction A named list with intensity scores to use as ylimit for top annotation.
#' If the value is single number, it has to be floating point number to extract the quantile and use
#' limit \code{[0, quantile(x)]}. If the value is numeric vector of length 2, the two numbers are used as lower
#' and upper limits for ylimit of top annotation. Default: NULL.
#' @param keep Same as \code{EnrichedHeatmap::normalizeToMatrix}. First value is used as lower quantile
#' and any value in profile matrix less than lower quantile is set to lower quantile. Second value is
#' used as upper quantile and any value greater than upper quantile is set to upper quantile.
#' Default: \code{c(0, 1)}
#' @param rasterPar rasterization parameters as list() with two elements: use, qual. Default:
#' \code{list(use = TRUE, qual = 5)}. \code{rasterPar$use} is used for \code{use_raster} argument and
#' \code{rasterPar$qual} is used for \code{raster_quality} argument of \code{Heatmap} function.
#' @param ... Other arguments for \code{EnrichedHeatmap} function
#'
#' @return A named list of Enriched Heatmap objects
#' @export
#'
#' @examples NA
get_profile_plot_list <- function(exptInfo,
                                  cluster = NULL,
                                  clusterColor = NULL,
                                  clustOrd = NULL,
                                  geneList,
                                  colorList = NULL,
                                  matrixSource = "normalizedmatrix",
                                  matrixBins = c(200, 200, 100, 10),
                                  targetType = "region",
                                  targetName = "gene",
                                  ylimFraction = NULL,
                                  keep = c(0, 1),
                                  rasterPar = list(use = TRUE, qual = 5),
                                  ...){
  allPlots <- list()

  ## generate plot for each sample
  for(i in 1:nrow(exptInfo)){
    sampleName <- exptInfo$sampleId[i]

    cat("Generating profile heatmap for:", sampleName, "...\n")

    mat <- import_profile_from_file(file = exptInfo$matFile[i],
                                    source = matrixSource,
                                    signalName = sampleName,
                                    selectGenes = geneList,
                                    up = matrixBins[1],
                                    target = matrixBins[2],
                                    down = matrixBins[3],
                                    binSize = matrixBins[4],
                                    keep = keep,
                                    targetType = targetType,
                                    targetName = targetName)


    qtDistr <- quantile(mat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 1), na.rm = T)
    print(qtDistr)

    ## color: generate if colorList = NULL
    tfCol <- colorRamp2(quantile(mat, c(0.50, 0.99), na.rm = T), c("white", "red"))
    polIICol <- colorRamp2(quantile(mat, c(0.01, 0.5, 0.995), na.rm = T), c("blue", "white", "red"))

    colFun <- tfCol
    if(exptInfo$IP_tag[i] == "polII"){
      colFun <- polIICol
    }


    ## if the colors are provided with colorList argument, use those colors instead
    if(!is.null(colorList[[sampleName]])){
      colFun <- colorList[[sampleName]]
    } else{
      cat("Profile color for sample", sampleName, "not provided. Using internal default quantiles to generate color\n")
    }


    ## build the profile heatmap
    profile2 <- profile_heatmap(profileMat = mat,
                                signalName = exptInfo$profileName[i],
                                columnTitle = exptInfo$sampleName[i],
                                geneGroups = cluster,
                                profileColor = colFun,
                                clusterColor = clusterColor,
                                clusterOrder = clustOrd,
                                ylimFraction = ylimFraction[[sampleName]],
                                rasterPar = rasterPar,
                                ...)

    ## append to the plot list
    allPlots[[ sampleName ]] <- profile2

  }

  return(allPlots)
}

##################################################################################




##################################################################################
## get the pol-II expression heatmap in list() of class ComplexHeatmap for many samples
#' Multiple PolII expression heatmaps as list
#'
#' @param expDf a dataframe with PolII expression columns
#' @param exptInfo experiment info as data frame with information like sampleID, type, path etc
#' @param htColor heatmap color generated by \code{ColorRamp2}
#' @param genes A vector of gene IDs which are to be plotted
#'
#' @return A named list object with multiple Heatmaps
#' @export
#'
#' @examples NA
get_expression_heatmap_list <- function(expDf,
                                        exptInfo,
                                        genes,
                                        htColor = NULL){

  polII_ids <- exptInfo$sampleId[which(exptInfo$IP_tag == "polII")]
  polII_expIds <- paste("is_expressed.", polII_ids, sep = "")

  ## select the gene of interst for plotting
  geneSet <- data.frame(geneId = unique(genes), stringsAsFactors = F) %>%
    dplyr::left_join(y = expDf, by = c("geneId" = "geneId"))

  ## clustering polII expression data
  rownames(geneSet) <- geneSet$geneId
  polII_mat <- as.matrix(geneSet[polII_ids])

  polII_log2_mat <- log2(polII_mat + 1)

  ## polII Heatmap
  quantile(polII_log2_mat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 1), na.rm = T)


  polII_color <- colorRamp2(
    breaks = c(0, quantile(polII_log2_mat, c(0.5, 0.8, 0.9, 0.92, 0.95, 0.97, 0.99, 0.995, 0.999))),
    colors = c("white", RColorBrewer::brewer.pal(n = 9, name = "RdPu"))
  )

  # colorRamp2(breaks = quantile(polII_log2_mat, 0, 0.99),colors =  c("white", "red"))

  if(!is.null(htColor)){
    polII_color <- htColor
  }

  allPlots <- list()
  firstPlot <- TRUE

  ## generate heatmaps
  for(id in polII_ids){
    ht <- signal_heatmap(log2_matrix = polII_log2_mat[, id, drop = FALSE],
                         col_title = id,
                         column_title_gp = gpar(fontsize = 12),
                         legend_title = "log2(polII_FPKM + 1)",
                         color = polII_color,
                         showLegend = firstPlot)

    firstPlot <- FALSE

    allPlots[[id]] <- ht
  }

  allPlots[["polII_color"]] <- polII_color

  return(allPlots)
}


##################################################################################

