
#' column wise z-score scaling of matrix
#'
#' This function is similar to the \code{base::scale()} function. This function is
#' copied from the tutorial posted by user \strong{strictlystat} on R-Bloggers.
#' Reference: \url{https://www.r-bloggers.com/a-faster-scale-function/}
#'
#' IMP: Remember that default \code{base::scale()} function has different behaviour
#' when \code{center == FALSE}. If \code{center} is \code{TRUE}, the scaling is done
#' by dividing centered columns of x by standard deviation. If \code{center} is
#' \code{FALSE}, the root mean square of column. Refer to \code{\link[base]{scale}}
#' for details.
#'
#'
#' @param x A matrix
#' @param center Logical: whether to center each column by column mean or not.
#' Default: TRUE
#' @param scale Logical: scale each column by standard deviation or not. Default: TRUE
#' @param add_attr Logical: whether to add attributes to scalled matrix. Default: FALSE
#' @param rows Row indices if only a subset of rows need to be used. Default: NULL
#' i.e. all rows are used
#' @param cols Column indices if only a subset of columns needs to be used. Default: NULL
#' i.e. all columns are used
#'
#' @return A column centered matrix
#' @export
#'
#' @examples NA
colScale = function(x,
                    center = TRUE,
                    scale = TRUE,
                    add_attr = TRUE,
                    rows = NULL,
                    cols = NULL) {

  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }


  # Get the column means
  cm = colMeans(x, na.rm = TRUE)

  # Get the column sd
  if (scale) {
    csd = matrixStats::colSds(x, center = cm)
  } else {
    # just divide by 1 if not
    csd = rep(1, length = length(cm))
  }

  if (!center) {
    # just subtract 0
    cm = rep(0, length = length(cm))
  }

  x = t( (t(x) - cm) / csd )

  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
  }

  return(x)
}

##################################################################################

#' row wise z-score scaling of matrix
#'
#' This function is similar to the \code{base::scale()}  function. This function is
#' copied from the tutorial posted by user \strong{strictlystat} on R-Bloggers.
#' Reference: \url{https://www.r-bloggers.com/a-faster-scale-function/}
#'
#' IMP: Remember that default \code{base::scale()} function has different behaviour
#' when \code{center == FALSE}. If \code{center} is \code{TRUE}, the scaling is done
#' by dividing centered columns of x by standard deviation. If \code{center} is
#' \code{FALSE}, the root mean square of column. Refer to \code{\link[base]{scale}}
#' for details.
#'
#' @param x A matrix
#' @param center Logical: whether to center each row by row mean or not.
#' Default: TRUE
#' @param scale Logical: scale each row by standard deviation or not. Default: TRUE
#' @param add_attr Logical: whether to add attributes to scalled matrix. Default: FALSE
#' @param rows Row indices if only a subset of rows need to be used. Default: NULL
#' i.e. all rows are used
#' @param cols Column indices if only a subset of columns needs to be used. Default: NULL
#' i.e. all columns are used
#'
#' @return A row scaled matrix
#' @export
#'
#' @examples NA
rowScale = function(x,
                    center = TRUE,
                    scale = TRUE,
                    add_attr = TRUE,
                    rows = NULL,
                    cols = NULL) {

  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }

  # Get the row means
  cm = rowMeans(x, na.rm = TRUE)

  # Get the row sd
  if (scale) {
    csd = matrixStats::rowSds(x, center = cm)
  } else {
    # just divide by 1 if not
    csd = rep(1, length = length(cm))
  }

  if (!center) {
    # just subtract 0
    cm = rep(0, length = length(cm))
  }

  x = (x - cm) / csd
  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
  }

  return(x)
}


##################################################################################
