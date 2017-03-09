plot.sharpel1pca <- function(x, ...) {
  if(!inherits(x,"sharpel1pca"))
    stop("Not a sharpel1pca object")
  
  if(ncol(x$scores) == 1)
    stop("Need scores in at least two dimensions")

  plot(x$scores[,1:2])
}
