sharpel1rs <- function (X, projDim=1, center=TRUE, projections="none") 
{
  if (!inherits(X, "matrix")) {
    if (inherits(X, "data.frame"))
      X <- as.matrix(X)
    else
      X <- matrix(X, ncol = 1)
  }
  if(center){
    X <- apply(X,2,function(y) y - median(y, na.rm = TRUE));
  }

  X <- t(X)

  pcLength    <- nrow(X) * projDim
  print(pcLength)
  scoreLength <- 0
  projLength <- 0
  objLength <- 0
  if (projections != "none") {
    scoreLength <- projDim * ncol(X)
  }
  if (projections != "none") {
    projLength  <- nrow(X) * ncol(X)
  }

  sol <- .C(C_sharpel1rs, as.double(X), as.integer(dim(X)), as.integer(projDim), loadings=double(pcLength), objectives=double(projDim), PACKAGE="pcaL1", NAOK = TRUE)

  solution <- new.env()
  solution$loadings <- matrix(sol[["loadings"]], ncol=projDim, byrow=FALSE)

  solution$minobjectives <- sol[["objectives"]]
  
  if (projections == "l1") {
    myl1projection <- l1projection(t(X), as.matrix(solution$loadings[,1:projDim]))
    solution$projPoints <- myl1projection$projPoints
    #solution$projPoints <- as.matrix(t((solution$loadings) %*% t(solution$loadings) %*% X))
    row.names(solution$projPoints) <- colnames(X)
    #solution$scores <- matrix(sol[["scores"]], ncol=projDim, byrow=FALSE)
    solution$scores <- myl1projection$scores
    rownames(solution$scores) <- colnames(X)
    totalDisp <- sum(abs(X))
    scoreDisp <- apply(abs(solution$scores), 2, sum)
    solution$dispExp <- scoreDisp/totalDisp
  } else if (projections == "l2") {
    myl2projection <- l2projection(t(X), as.matrix(solution$loadings[,1:projDim]))
    solution$projPoints <- myl2projection$projPoints
    #solution$projPoints <- as.matrix(t((solution$loadings) %*% t(solution$loadings) %*% X))
    row.names(solution$projPoints) <- colnames(X)
    #solution$scores <- matrix(sol[["scores"]], ncol=projDim, byrow=FALSE)
    solution$scores <- myl2projection$scores
    rownames(solution$scores) <- colnames(X)
    totalDisp <- sum(apply(t(X),2, var))
    scoreDisp <- apply(solution$scores, 2, var)
    solution$dispExp <- scoreDisp/totalDisp
  }
  
  solution <- as.list(solution)
  class(solution) <- "sharpel1rs"
  solution
}

