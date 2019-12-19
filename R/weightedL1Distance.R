weightedL1Distance <- function(X, loadings, weights) # take a subspace and data as input, get distances to subspace for each point
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

  # get distance to projection onto each basis vector.  Add distances.

  wBasisDist <- function(j) {
    sol <- .C (C_l1projection, as.double(X), as.integer(dim(X)), as.integer(1), as.double(loadings[,j]), projPoints=double(projLength), alphas=double(ncol(X)*1), PACKAGE="pcaL1")

    myProjPoints <- matrix(sol[["projPoints"]], ncol=nrow(X), byrow=TRUE)
    wDistance <- weights[j] * apply(abs(t(X) - myProjPoints), 1,sum)
    wDistance
  }
  wDistances <- lapply(1:ncol(loadings), wBasisDist)

  output <- new.env()
  output$wDistances <- wDistances
  output$totalDistance <- Reduce("+", wDistances)
  output <- as.list(output)
  class(output) <- "weightedL1Distance"
  output
}
