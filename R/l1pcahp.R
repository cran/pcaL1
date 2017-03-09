l1pcahp <- function (X, projDim=1, center=TRUE, projections="none", initialize="l2pca", threshold=0.0001)
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
    if(nrow(X) < ncol(X))
    {
      mypca <- prcomp(X, center=center)
      initV <- data.frame(mypca$rotation[,ncol(mypca$rotation):1])
      rand  <- matrix(runif((ncol(X) - nrow(X))*ncol(X), min = -1, max = 1), nrow =ncol(X))
      rand  <- data.frame(sweep(rand, 2, sqrt(apply(rand,2, function(x){ sum(x^2)})), FUN="/"))
      initV <- as.matrix(cbind(rand,initV))
    } else {
      mypca <- prcomp(X, center=center)
      initV <- mypca$rotation[,ncol(mypca$rotation):1]
    }
    
  }  else if (initialize == "random") {
    initV <- matrix(runif(ncol(X)*ncol(X), min = -1, max = 1), ncol=ncol(X))
    initV <- sweep(initV, 2, sqrt(apply(initV,2, function(x){ sum(x^2)})), FUN="/")
    initV <- matrix(runif(ncol(X)*ncol(X)), ncol=ncol(X))
  }   else {
    # return an error
  }
  X <- t(X)
  pcLength    <- nrow(X) * nrow(X)
  sol <- .C (C_l1pcahp, as.double(X), as.integer(dim(X)), as.double(threshold), as.double(unlist(initV)),loadings=double(pcLength),  PACKAGE="pcaL1")
  
  solution <- new.env()
  solution$loadings <- matrix(sol[["loadings"]], ncol=nrow(X), byrow=TRUE) 
  solution$loadings <-  solution$loadings[nrow(X):1, ]
  solution$loadings <-  t(solution$loadings[1:projDim, ])
  
  if (projDim == 1)
    solution$loadings <- t(solution$loadings)
  
  
  if (projections=="l1") {
    myl1projection <- l1projection(t(X), as.matrix(solution$loadings))
    #solution$scores <- as.matrix(t(X) %*% solution$loadings)
    solution$scores <- myl1projection$scores

    totalDisp         <- sum(abs(X))
    scoreDisp         <- apply(abs(solution$scores), 2, sum)
    solution$dispExp <- scoreDisp/totalDisp

    #solution$projPoints <- as.matrix(t((solution$loadings) %*% t(solution$loadings) %*% X))
    solution$projPoints <- myl1projection$projPoints
    row.names(solution$projPoints) <- colnames(X)
  } else if (projections=="l2") {
    myl2projection <- l2projection(t(X), as.matrix(solution$loadings))
    #solution$scores <- as.matrix(t(X) %*% solution$loadings)
    solution$scores <- myl2projection$scores

    totalDisp         <- sum(apply(t(X), 2, var))
    scoreDisp         <- apply(solution$scores,2,var)
    solution$dispExp <- scoreDisp/totalDisp

    #solution$projPoints <- as.matrix(t((solution$loadings) %*% t(solution$loadings) %*% X))
    solution$projPoints <- myl1projection$projPoints
    row.names(solution$projPoints) <- colnames(X)
  }

  solution <- as.list(solution)
  class(solution) <- "l1pcahp"
  solution
}
