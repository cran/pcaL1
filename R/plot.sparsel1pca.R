plot.sparsel1pca <- function(x, ...) {
  if(!inherits(x,"sparsel1pca"))
    stop("Not a sharpel1rs object")
  
  if(ncol(x$scores) == 1)
    stop("Need scores in at least two dimensions")

  plot(x$scores[,1:2])
}
