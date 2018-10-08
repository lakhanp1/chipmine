
pkg.env = new.env(parent = emptyenv())

## custom colors
pkg.env$customColors = c("blue", "red", "yellowgreen", "brown", "green", "purple", "deeppink", "dodgerblue", "forestgreen", "goldenrod", "hotpink", "navyblue", "indianred", "khaki3", "magenta", "maroon", "orange", "coral", "chocolate", "cyan", "olivedrab", "burlywood", "orchid", "peachpuff", "peru", "plum", "rosybrown", "aquamarine", "saddlebrown", "sienna", "slateblue", "cornflowerblue", "tan", "thistle", "tomato", "turquoise", "violet", "violetred", "wheat4", "yellow")


pkg.env$clusterColorSet = structure(pkg.env$customColors,
                                    names = paste("Cluster", 1:(length(pkg.env$customColors)), sep = "_")
)

