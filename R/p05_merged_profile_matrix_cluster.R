## Merged profile matrix clustering
#' Merged profile matrix clustering
#'
#' This function join the multiple profile matrix and performs the k-means clustering
#' on the merged profile matrix
#'
#'
#' @param name name of the newly mearged profile matrix
#' @param exptInfo experiment info as data frame with information like sampleID, type, path etc
#' @param genes A vector of gene IDs which are to be plotted
#' @param clusterStorePath Path for storing/accessing the cluster information file.
#' @param matrixBins A numeric vector with four elements representing the column
#' counts for up, body, down regions and binSize for the profile matrix.
#' Default: c(200, 200, 100, 10)
#' @param source Method by which profile matrix was generated.
#'  One of "deeptools", "miao" or "normalizedmatrix".
#' @param matrixBins A numeric vector with four elements
#' @param k number of cluster centers
#' @param kmIter iter.max argument of stats::kmeans() function. Default: 200
#' @param kmStarts nstart argument of stats::kmeans() function. Default: 200
#' @param ... other arguments for stats::kmeans() function
#'
#' @return A list with two elements:
#' \itemize{
#' \item profileKm: kmeans clustering object returned by function profile_matrix_kmeans()
#' \item mat: a merged profile matrix as dataframe
#' }
#' @export
#'
#' @examples NA
merged_profile_matrix_cluster <- function(name, exptInfo, genes, clusterStorePath,
                                          matrixBins = c(200, 200, 100, 10), source = "deeptools",
                                          k, kmIter = 200, kmStarts = 200, ...){

  ## read the profile matrix as dataframe
  matDf1 <- import_profile_from_file(file = exptInfo$matFile[1],
                                     source = source,
                                     signalName = exptInfo$sampleId[1],
                                     selectGenes = genes,
                                     up = matrixBins[1],
                                     target = matrixBins[2],
                                     down = matrixBins[3],
                                     binSize = matrixBins[4],
                                     returnDf = T)


  matDf <- matDf1

  ## apend the remaining profile matrices
  for (i in 2:length(exptInfo$sampleId)) {
    ## read profile matrix as dataframe
    matDf2 <- import_profile_from_file(file = exptData$matFile[i],
                                       source = source,
                                       signalName = exptData$sampleId[i],
                                       selectGenes = genes,
                                       up = matrixBins[1],
                                       target = matrixBins[2],
                                       down = matrixBins[3],
                                       binSize = matrixBins[4],
                                       returnDf = T)

    ## merge profile matrices
    matDf <- dplyr::left_join(x = matDf, y = matDf2, by = c("geneId" = "geneId"))

  }

  matDf <- tibble::column_to_rownames(df = matDf, var = "geneId")

  profileMat <- data.matrix(matDf)


  ## remove the rows with all NA values
  allNaRows <- which(apply(profileMat, 1, function(x) all(is.na(x))))

  if (length(allNaRows) > 0) {
    warning("Removing ", length(allNaRows), " genes found with all NA values while reading profile matrix")
    profileMat <- profileMat[-allNaRows, ]
  }


  ## find rows with NA values and ipmute the NA values: replace NA with mean(row)
  naRows <- which(apply(profileMat, 1, function(x) any(is.na(x))))

  profileMat[naRows, ] <- do.call(rbind, lapply(naRows, na_impute, mat = profileMat))

  ## k-means clustering on the merged profile matrix
  profileKm <- profile_matrix_kmeans(mat = profileMat,
                                     km = k,
                                     clustFile = clusterStorePath,
                                     name = comparisonName,
                                     kmIter = kmIter,
                                     kmStarts = kmStarts,
                                     ...)

  return(list(
    "profileKm" = profileKm,
    "mat" = matDf)
  )

}











