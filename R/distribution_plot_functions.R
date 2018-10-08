
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
#'
#' @return a ggplot2 object which has MA plot
#' @export
#'
#' @examples NA
plot_MA_gg = function(df, s1, s2, title = "MA plot", colorCol, pseudoCount = 0.1){


  df2 <- df %>%
    dplyr::mutate(
      M = log2(!!as.name(s1) + pseudoCount) - log2(!!as.name(s2) + pseudoCount),
      A = (log2(!!as.name(s1) + pseudoCount) + log2(!!as.name(s2) + pseudoCount)) / 2
    )

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

  maPlot <- ggplot() +
    geom_point(data = df2, mapping = aes(x = A, y = M, color = !!as.name(colorCol))) +
    ggtitle(title) +
    xlab("Mean expression (log2)") +
    ylab("log2(fold change)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20, face = "bold"),
          legend.box.background = element_rect(colour = "black"),
          legend.box.margin = margin(1, 1, 1, 1)
    )


  return(maPlot)
}


##################################################################################

#' Generate XY scatter plot for two sets
#'
#' @param df dataframe with values two columns for which scatter plot is to be generated
#' @param s1 sample 1 column name
#' @param s2 sample 2 column name
#' @param title title to be showed on plot
#' @param colorCol column name which should be used for coloring the points
#' @param transformation trans argument from ggplot::scale_x_continuous() function.
#' Default: log2
#' @param pseudoCount A small number to be added to the values to avoid log(0) error.
#' Default: 0.1
#'
#' @return a ggplot object
#' @export
#'
#' @examples NA
plot_scatter <- function(df, s1, s2, title = "Scatter plot", colorCol, transformation = "log2", pseudoCount = 0.1){

  df2 <- df %>%
    dplyr::mutate(
      s1 = !!as.name(s1) + pseudoCount,
      s2 = !!as.name(s2) + pseudoCount
    )


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


  sp <- ggplot() +
    geom_point(data = df2, mapping = aes(x = s2, y = s1, color = !!as.name(colorCol))) +
    ggtitle(title) +
    scale_x_continuous(name = paste(transformation, "(",s2, ")", sep = ""),
                       trans = transformation) +
    scale_y_continuous(name = paste(transformation, "(",s1, ")", sep = ""),
                       trans = transformation) +
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



