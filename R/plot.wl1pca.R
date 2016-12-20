plot.wl1pca <- function(x, ...) {
  if(!inherits(x,"wl1pca"))
    stop("Not a wl1pca object")
  
  if(ncol(x$scores) == 1)
    stop("Need scores in at least two dimensions")

  plot(x$scores[,1:2])
}
