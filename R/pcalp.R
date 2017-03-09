pcalp <- function (X, projDim=1, p = 1.0, center=TRUE, projections="none", initialize="l2pca",solution = "L", epsilon = 0.0000000001, lratio = 0.02)
{
  if (class (X) != "matrix") {
    if (class (X) == "data.frame")
      X <- as.matrix(X)
    else
      X <- matrix(X, ncol = 1)
  }

  if (center) {
    X <- apply(X,2,function(y) y - median(y));
  }

  A <- X
  X <- X[apply(abs(X),1,sum) > 0,] # get rid of origin points for algorithm
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

  solMethod <- 0
  if (solution == "G") {
    solMethod <- 1  
  }
  else if (solution == "L") {
    solMethod <- 2
  }  
  sol <- .C(C_pcalp, as.double(X), as.integer(dim(X)), as.integer(projDim), as.double(p), loadings=double(pcLength), as.integer(initMethod), as.integer(solMethod), initV, as.double(epsilon), as.double(lratio), PACKAGE="pcaL1")
  
  solution <- new.env()
  solution$loadings <- matrix(sol[["loadings"]], ncol=projDim, byrow=FALSE) 

  if (projections == "l1") {
    myl1projection <- l1projection(t(X), as.matrix(solution$loadings[,1:projDim]))
    solution$projPoints <- myl1projection$projPoints
    row.names(solution$projPoints) <- colnames(X)
    solution$scores <- myl1projection$scores

    #solution$scores <- as.matrix(A %*% solution$loadings)
    totalDisp         <- sum(abs(A))
    scoreDisp         <- apply(abs(solution$scores),2,sum)
    solution$dispExp <- scoreDisp/totalDisp
  } else if (projections == "l2") {
    myl2projection <- l2projection(t(X), as.matrix(solution$loadings[,1:projDim]))
    solution$projPoints <- myl2projection$projPoints
    row.names(solution$projPoints) <- colnames(X)
    solution$scores <- myl2projection$scores

    #solution$scores <- as.matrix(A %*% solution$loadings)
    totalDisp         <- sum(apply(A, 2, var))
    scoreDisp         <- apply(solution$scores,2,var)
    solution$dispExp <- scoreDisp/totalDisp
  }
  
  solution <- as.list(solution)
  class(solution) <- "pcalp"
  solution
}  
