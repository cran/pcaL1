l1projection <- function (X, loadings)
{
  if (!inherits(X, "matrix")) {
    if (inherits(X, "data.frame"))
      X <- as.matrix(X)
    else
      X <- matrix(X, ncol = 1)
  }
  if (!inherits(loadings, "matrix")) {
    if (inherits(loadings, "data.frame"))
      loadings <- as.matrix(loadings)
    else
      loadings <- matrix(loadings, ncol = 1)
  }

  projDim <- ncol(loadings)

  X <- t(X)

  projLength <- nrow(X) * ncol(X) 

  sol <- .C (C_l1projection, as.double(X), as.integer(dim(X)), as.integer(projDim), as.double(loadings), projPoints=double(projLength), alphas=double(ncol(X)*projDim), PACKAGE="pcaL1")

  solution <- new.env()

  solution$projPoints <- matrix(sol[["projPoints"]], ncol=nrow(X), byrow=TRUE)
  solution$scores <- matrix(sol[["alphas"]], ncol=projDim, byrow=FALSE)
  row.names(solution$projPoints) <- colnames(X)
  
  solution <- as.list(solution)
  class(solution) <- "l1projection"
  solution
}
