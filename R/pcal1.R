pcal1 <- function (X, projDim=1, center=TRUE, projections="none", initialize="l2pca")
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

  A <- X
  X <- t(X)
  pcLength <- projDim * (nrow(X))

  initV <- numeric(ncol(X))
  initMethod <- 0
  if (is.numeric(initialize)) {
    initV <- initialize
  }
  else if (initialize == "maxx") {
    initMethod <- 1
  } 
  else if (initialize == "random") {
    initMethod <- 2
  }
  sol <- .C(C_pcal1, as.double(X), as.integer(dim(X)), as.integer(projDim), loadings=double(pcLength), as.integer(initMethod), initV, PACKAGE="pcaL1")
  
  solution <- new.env()
  solution$loadings <- matrix(sol[["loadings"]], ncol=projDim, byrow=FALSE) 

  if (projections == "l1") {
    myl1projection <- l1projection(t(X), as.matrix(solution$loadings[,1:projDim]))
    #solution$projPoints <- as.matrix(t((solution$loadings) %*% t(solution$loadings) %*% X))
    solution$projPoints <- myl1projection$projPoints
    row.names(solution$projPoints) <- colnames(X)
    solution$scores <- myl1projection$scores

    #solution$scores <- as.matrix(A %*% solution$loadings)
    totalDisp         <- sum(abs(A))
    scoreDisp         <- apply(abs(solution$scores), 2, sum)
    solution$dispExp <- scoreDisp/totalDisp
  } else if (projections == "l2") {
    myl2projection <- l2projection(t(X), as.matrix(solution$loadings[,1:projDim]))
    #solution$projPoints <- as.matrix(t((solution$loadings) %*% t(solution$loadings) %*% X))
    solution$projPoints <- myl2projection$projPoints
    row.names(solution$projPoints) <- colnames(X)
    solution$scores <- myl2projection$scores

    #solution$scores <- as.matrix(A %*% solution$loadings)
    totalDisp         <- sum(apply(A,2,var))
    scoreDisp         <- apply(solution$scores,2,var)
    solution$dispExp <- scoreDisp/totalDisp
  }
  

  solution <- as.list(solution)
  class(solution) <- "pcal1"
  solution

}  
