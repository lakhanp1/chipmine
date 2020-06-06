
#' Generate MA plot
#'
#' This method generates MA plot using two numeric columns from dataframe
#' It can also color the genes based on the colorCol variable
#'
#' @param df dataframe with values two columns for which MA plot is to be generated
#' @param s1 sample 1 column name
#' @param s2 sample 2 column name
#' @param title title to be showed on plot
#' @param colorCol column name which should be used for coloring the points
#' @param pseudoCount A small number to be added to the values to avoid log(0) error.
#' Default: 0.1
#' @param ylim A two element numeric vector for Y-limits of MA plot. Default: c(-4, 4)
#'
#' @return a ggplot2 object which has MA plot
#' @export
#'
#' @examples NA
plot_MA_gg = function(df, s1, s2, title = NA, colorCol, pseudoCount = 0.1, ylim = NULL){

  ## check for valid ylim
  if(!is.null(ylim) & length(ylim) != 2 & !is.numeric(ylim)){
    stop("wrong ylim option. Please provide ylim as numeric vector of length 2")
  } else{
    ylim = sort(ylim)
  }

  df2 <- df %>%
    dplyr::mutate(
      M = log2(!!sym(s1) + pseudoCount) - log2(!!sym(s2) + pseudoCount),
      A = (log2(!!sym(s1) + pseudoCount) + log2(!!sym(s2) + pseudoCount)) / 2
    ) %>%
    {
      if(!is.null(ylim)){
        dplyr::mutate(.data = .,
                      M = scales::squish(M, range = ylim))
      } else{
        .
      }
    }


  ## if the column for color does not exists, create new column with value "genes"
  ## this also means no group wise coloring
  if(!any(colorCol %in% names(df2))){
    df2[[colorCol]] <- "genes"
  }

  ## if the color column is not character, convert it to character.
  ## for numeric and interger type column, the color is continuous which we dont want
  if(class(df2[[colorCol]]) != "character"){
    df2[[colorCol]] <- as.character(df2[[colorCol]])
  }

  if(is.na(title)){
    title <- paste("MA plot: ", s1, "/", s2)
  }

  maPlot <- ggplot() +
    geom_point(data = df2, mapping = aes(x = A, y = M, color = !!sym(colorCol))) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1, linetype = 2, alpha = 0.5) +
    geom_hline(yintercept = -1, linetype = 2, alpha = 0.5) +
    ggtitle(title) +
    xlab("Mean expression (log2)") +
    ylab("log2(fold change)") +
    theme_bw() +
    theme_classic()



  return(maPlot)
}


##################################################################################

#' Generate XY scatter plot for two sets
#'
#' @param df dataframe with values two columns for which scatter plot is to be generated
#' @param sx sample 1 column name
#' @param sy sample 2 column name
#' @param title title to be showed on plot
#' @param colorCol column name which should be used for coloring the points
#' @param trans trans argument from ggplot::scale_x_continuous() function.
#' Default: identity
#' @param pseudoCount A small number to be added to the values to avoid log(0) error.
#' Default: NULL
#'
#' @return a ggplot object
#' @export
#'
#' @examples NA
plot_scatter <- function(df, sx, sy, title = "Scatter plot", colorCol = "Density",
                         trans = "identity", pseudoCount = NULL){

  trans <- scales::as.trans(trans)

  if(is.numeric(pseudoCount)){
    df <- df %>%
      dplyr::mutate(
        !!sx := pmax(!!sym(sx), pseudoCount),
        !!sy := pmax(!!sym(sy), pseudoCount)
      )
  }

  ## if the column for color does not exists, use density color by default
  if(!any(colorCol %in% names(df))){
    colorCol = "Density"
    if(nrow(df) >= 10){
      densColor <- densCols(x = trans$transform(df[[ sx ]]),
                            y = trans$transform(df[[ sy ]]),
                            colramp = colorRampPalette(c("black", "white"))
      )

      df[[colorCol]] <- col2rgb(densColor)[1,] + 1L
    } else{
      df[[colorCol]] <- rep(1, nrow(df))
    }
  }

  discreteColor <- FALSE
  if(is.character(df[[colorCol]])){
    discreteColor <- TRUE
  }

  if(trans$name != "identity"){
    xname <- paste(trans$name, "(",sx, ")", sep = "")
    yname <- paste(trans$name, "(",sy, ")", sep = "")
  } else{
    xname <- sx
    yname <- sy
  }

  sp <- ggplot(data = df,
               mapping = aes(x = !!sym(sx), y = !!sym(sy), color = !!sym(colorCol))) +
    geom_point() +
    scale_x_continuous(trans = trans) +
    scale_y_continuous(trans = trans) +
    viridis::scale_color_viridis(option = "plasma", discrete = discreteColor) +
    labs(title = title, x = xname, y = yname) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20, face = "bold"),
          legend.box.background = element_rect(colour = "black"),
          legend.box.margin = margin(1, 1, 1, 1)
    )

  return(sp)
}

##################################################################################



