plot.awl1pca <- function(x, ...) {
  if(!inherits(x,"awl1pca"))
    stop("Not an awl1pca object")
  
  if(ncol(x$scores) == 1)
    stop("Need scores in at least two dimensions")

  plot(x$scores[,1:2])
}
