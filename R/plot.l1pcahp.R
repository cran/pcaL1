plot.l1pcahp <- function(x, ...) {
  if(!inherits(x,"l1pcahp"))
    stop("Not an l1pcahp object")
  
  if(ncol(x$scores) == 1)
    stop("Need scores in at least two dimensions")

  plot(x$scores[,1:2])
}
