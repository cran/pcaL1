plot.sharpel1rs <- function(x, ...) {
  if(!inherits(x,"sharpel1rs"))
    stop("Not a sharpel1rs object")
  
  if(ncol(x$scores) == 1)
    stop("Need scores in at least two dimensions")

  plot(x$scores[,1:2])
}
