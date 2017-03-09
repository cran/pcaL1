l1pcastar <- function (X, projDim=1, center=TRUE, projections="none")
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

  X <- t(X)

  pcLength    <- nrow(X) * nrow(X)
  scoreLength <- 0
  projLength <- 0
  if (projections != "none") {
    scoreLength <- projDim * ncol(X)
  }
  if (projections != "none") {
    projLength  <- nrow(X) * ncol(X)
  }

  #sol <- .C ("l1pcastar", as.double(X), as.integer(dim(X)), as.integer(projDim), as.integer(scores), as.integer(projPoints), loadings=double(pcLength), scores=double(scoreLength), projPoints=double(projLength), PACKAGE="pcaL1")
  # using L1 projection to find projected points
  sol <- .C (C_l1pcastar, as.double(X), as.integer(dim(X)), as.integer(projDim), loadings=double(pcLength), PACKAGE="pcaL1")

  solution <- new.env()
  solution$loadings <- matrix(sol[["loadings"]], ncol=nrow(X), byrow=FALSE)
  
  if (projections == "l1") {
    myl1projection <- l1projection(t(X), as.matrix(solution$loadings[,1:projDim]))
    solution$scores <- myl1projection$scores

    row.names(solution$scores) <- colnames(X)
    totalDisp <- sum(abs(X))
    scoreDisp <- apply(abs(solution$scores), 2, sum)
    solution$dispExp <- scoreDisp/totalDisp
  
    #solution$projPoints <- matrix(sol[["projPoints"]], ncol=nrow(X), byrow=TRUE)  
    solution$projPoints <- myl1projection$projPoints
    row.names(solution$projPoints) <- colnames(X)
  } else if (projections == "l2") {
    myl2projection <- l2projection(t(X), as.matrix(solution$loadings[,1:projDim]))
    solution$scores <- myl2projection$scores

    row.names(solution$scores) <- colnames(X)
    totalDisp <- sum(apply(t(X),2,var))
    scoreDisp <- apply(solution$scores, 2, var)
    solution$dispExp <- scoreDisp/totalDisp
  
    #solution$projPoints <- matrix(sol[["projPoints"]], ncol=nrow(X), byrow=TRUE)  
    solution$projPoints <- myl2projection$projPoints
    row.names(solution$projPoints) <- colnames(X)
  }
  
  solution <- as.list(solution)
  class(solution) <- "l1pcastar"
  solution
}
