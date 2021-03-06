l1pca <- function (X, projDim=1, center=TRUE, projections="l1", initialize="l2pca", tolerance=0.0001, iterations=10)
{
  if (!inherits(X, "matrix")) {
    if (inherits(X, "data.frame"))
      X <- as.matrix(X)
    else
      X <- matrix(X, ncol = 1)
  }
  if(center){
    X <- apply(X,2,function(y) y - median(y));
  }

  if (is.matrix(initialize)) {
    initV <- initialize
  }
  else if (initialize == "l2pca") {
    mypca <- prcomp(X, center=center)
    initV <- mypca$rotation[,1:projDim]
  }
  else if (initialize == "random") {
    initV <- matrix(runif(ncol(X)*projDim), ncol=projDim)
  }
  else {
    # return an error
  }

  X <- t(X)

  pcLength    <- projDim * nrow(X)
  scoreLength <- projDim * ncol(X)
  
  sol <- .C (C_l1pca, as.double(X), as.integer(dim(X)), as.integer(projDim), as.double(tolerance), as.integer(iterations), as.double(initV), loadings=double(pcLength), scores=double(scoreLength), PACKAGE="pcaL1")

  solution <- new.env()
  solution$loadings <- matrix(sol[["loadings"]], ncol=projDim, byrow=FALSE)
  
  #solution$scores <- matrix(sol[["scores"]], ncol=projDim, byrow=FALSE)
  #row.names(solution$scores) <- colnames(X)
  
  #solution$projPoints <- as.matrix(t((solution$loadings) %*% t(solution$loadings) %*% X))
  if (projections == "l1") {
    myl1projection <- l1projection(t(X), as.matrix(solution$loadings[,1:projDim]))
    solution$projPoints <- myl1projection$projPoints
    row.names(solution$projPoints) <- colnames(X)
    solution$scores <- myl1projection$scores

    totalDisp <- sum(abs(X))
    scoreDisp <- apply(abs(solution$scores), 2, sum)
  } else {
    myl2projection <- l2projection(t(X), as.matrix(solution$loadings[,1:projDim]))
    solution$projPoints <- myl2projection$projPoints
    row.names(solution$projPoints) <- colnames(X)
    solution$scores <- myl2projection$scores

    totalDisp <- sum(apply(t(X),2,var))
    scoreDisp <- apply(solution$scores, 2, var)
  }

  solution$dispExp <- scoreDisp/totalDisp
  
  mysort            <- sort(solution$dispExp, decreasing=TRUE, index.return=TRUE)
  
  solution$scores <- as.matrix(solution$scores[,mysort$ix])
  solution$loadings <- matrix(solution$loadings[,mysort$ix], ncol=projDim, byrow=FALSE)
  solution$dispExp  <- solution$dispExp[mysort$ix]

  solution <- as.list(solution)
  class(solution) <- "l1pca"
  solution
  
}
