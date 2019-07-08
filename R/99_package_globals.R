## chipmineEnv
chipmineEnv = new.env(parent = emptyenv())

## custom colors
customColors = c("blue", "red", "yellowgreen", "brown", "green", "purple", "deeppink", "dodgerblue", "forestgreen", "goldenrod", "hotpink", "navyblue", "indianred", "khaki3", "magenta", "maroon", "orange", "coral", "chocolate", "cyan", "olivedrab", "burlywood", "orchid", "peachpuff", "peru", "plum", "rosybrown", "aquamarine", "saddlebrown", "sienna", "slateblue", "cornflowerblue", "tan", "thistle", "tomato", "turquoise", "violet", "violetred", "wheat4", "yellow")

assign(x = "customColors", value = customColors, envir = chipmineEnv)

clusterColorSet = structure(
  customColors,
  names = paste("Cluster", 1:(length(customColors)), sep = "_")
)

assign(x = "clusterColorSet", value = clusterColorSet, envir = chipmineEnv)

##################################################################################

## txdbEnv
txdbEnv <- new.env(parent = emptyenv())



