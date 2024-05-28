## -----------------------------------------------------------------------------
mediumcontrast <- color("medium contrast")
colors_mc <- as.character(c(mediumcontrast(6),"#a08c8e", "#608da2"))
reds=colorRampPalette(c("white", colors_mc[5]))
blues=colorRampPalette(c(colors_mc[6], "white"))
browns=colorRampPalette(c(colors_mc[4], "white"))

vibrant <- color("vibrant")
colors_vibrant <- as.character(vibrant(6))

muted <- color("muted")
colors_mut <- as.character(muted(9))
smooth_rainbow <- color("smooth rainbow")


## -----------------------------------------------------------------------------
#plot_scheme(colors_mc, colours = TRUE, names = TRUE, size = 0.9) 