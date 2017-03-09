plot.pcalp <- function(x, ...) {
  if(!inherits(x,"pcalp"))
    stop("Not an pcalp object")
  
  if(ncol(x$scores) == 1)
    stop("Need scores in at least two dimensions")

  plot(x$scores[,1:2])
}
