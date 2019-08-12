

#' Profile Matrix k-means Clustering
#'
#' This function perform k-means clustering on the profile matrix using stats::kmeans() function.
#' It tries to cluster the matrix by default "Hartigan-Wong" algorithm. If this algorithm fails
#' at Quick-TRANSfer stage, clustering is done with "MacQueen"
#'
#' @param mat profile matrix
#' @param km number of cluster centers
#' @param clustFile File to which cluster data will be stored
#' @param name Sample name
#' @param kmIter iter.max	argument of stats::kmeans() function. Default: 100
#' @param kmStarts nstart argument of stats::kmeans() function. Default: 150
#' @param ... other arguments for stats::kmeans() function
#'
#' @return A list with two elements:
#' \itemize{
#' \item km: kmeans clustering object returned by stats::kmeans
#' \item geneClusters: a dataframe with columns: geneId, cluster
#' }
#' @export
#'
#' @examples NA
profile_matrix_kmeans = function(mat, km, clustFile, name, kmIter = 100, kmStarts = 150, ...){

  ##
  cat("Running k-means clustering for sample", name, "\nNumber of clusters:", km, "\n")

  clusterData = data.frame(geneId = rownames(mat), stringsAsFactors = F, row.names = rownames(mat))

  ## perform k-means clustering
  if(km > 1){

    km = kmeans(x = mat, centers = km, iter.max = kmIter, nstart = kmStarts)

    cat("ifault =", km$ifault, "\n")

    ## if "Hartigan-Wong" algorithm fails in "Quick-TRANSfer", do clustering using "MacQueen"
    if(km$ifault == 4){

      cat("Hartigan-Wong algorithm failed in Quick-TRANSfer stage. Clustering using MacQueen's algorithm\n")
      km = kmeans(x = mat, centers = km, iter.max = kmIter*2, nstart = kmStarts, algorithm = "MacQueen")

      cat("ifault =", km$ifault, "\n")

    }

    kmClust = data.frame(geneId = names(km$cluster), kmClusters = km$cluster, stringsAsFactors = F)

    ## find the max of each cluster center to decide the cluster order
    mx = apply(km$centers, MARGIN = 1, max)
    clustRank = rank(-mx)

    ## change the order of clusters based on the enrichment.
    # cluster with highest enrichment signal is cluster 1 and so on
    kmClust$cluster = clustRank[kmClust$kmClusters]

    # clusterData[names(km$cluster), "cluster"] = km$cluster
    clusterData = dplyr::left_join(x = clusterData, y = kmClust, by = c("geneId" = "geneId")) %>%
      dplyr::select(geneId, cluster)


  } else {
    clusterData$cluster = 1
  }

  clusterData = dplyr::mutate(clusterData, cluster = paste("Cluster_", cluster, sep = ""))

  ## store for future use
  write.table(x = clusterData, file = clustFile, sep = "\t", col.names = T, quote = F, row.names = F)

  return(
    list(
      "km" = km,
      "geneClusters" = clusterData
    )
  )

}
