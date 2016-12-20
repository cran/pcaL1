l1pcahp <- function (X, projDim=1, center=TRUE, projPoints=FALSE, scores=TRUE, initialize="l2pca", threshold=0.0001)
{
  
  
  if (class(X) != "matrix") {
    if (class(X) == "data.frame")
      X <- as.matrix(X)
    else
      X <- matrix(X, ncol = 1)
  }
  if(center){
    X <- apply(X,2,function(y) y - median(y));
  }

  if (is.matrix(initialize)) {
    initV <- initialize
  } else if (initialize == "l2pca") {
    mypca <- prcomp(X, center=center)
    initV <- mypca$rotation[nrow(mypca$rotation):1, ]
  }  else if (initialize == "random") {
    initV <- matrix(runif(ncol(X)*ncol(X)), ncol=ncol(X))
  }   else {
    # return an error
  }
  X <- t(X)
  initV <- t(initV)
  pcLength    <- nrow(X) * nrow(X)

  sol <- .C ("l1pcahp", as.double(X), as.integer(dim(X)), as.double(threshold), as.double(initV),loadings=double(pcLength),  PACKAGE="pcaL1")
  
  solution <- new.env()
  solution$loadings <- matrix(sol[["loadings"]], ncol=nrow(X), byrow=TRUE) 
  solution$loadings <-  solution$loadings[nrow(X):1, ]
  solution$loadings <-  t(solution$loadings[1:projDim, ])
  
  if (projDim == 1)
    solution$loadings <- t(solution$loadings)
  
  
  if (scores) {
    solution$scores <- as.matrix(t(X) %*% solution$loadings)
    totalDisp         <- sum(abs(X))
    scoreDisp         <- (apply(abs(solution$scores),2,sum))
    solution$dispExp <- scoreDisp/totalDisp
  }

  if (projPoints) {
    solution$projPoints <- as.matrix(t((solution$loadings) %*% t(solution$loadings) %*% X))
  }

  solution <- as.list(solution)
  class(solution) <- "l1pcahp"
  solution
}
