sparsel1pca <- function (X, projDim=1, center=TRUE, projections="none", lambda=0) 
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

  X <- t(X)

  pcLength    <- nrow(X) * projDim
  scoreLength <- 0
  projLength <- 0
  objLength <- 0
  if (projections != "none") {
    scoreLength <- projDim * ncol(X)
  }
  if (projections != "none") {
    projLength  <- nrow(X) * ncol(X)
  }

#  sol <- .C ("sparsel1pca", as.double(X), as.integer(dim(X)), as.integer(projDim), as.integer(scores), loadings=double(pcLength), scores=double(scoreLength), objectives=double(projDim), PACKAGE="pcaL1")

  if ((lambda < 0) & (ncol(X) > 100)){
    stop("Too many rows to calculate all lambdas. Specify a particular lambda.")
  } else if (lambda < 0) {
   pcLength    <- 1000*ncol(X)
   lambdas_out <- rep(-1,1000)
  } else {
    pcLength    <- nrow(X) * projDim
    lambdas_out <- c(lambda,0)
  }

  sol <- .C (C_sparsel1pca, as.double(X), as.integer(dim(X)), as.integer(projDim), loadings=double(pcLength), objectives=double(projDim), lambdas_out=as.double(lambdas_out), PACKAGE="pcaL1")
  

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

  if (lambda < 0) {
    num_lambdas <- length(lambdas_out)
    solution$loadings <- matrix(sol[["loadings"]][1:(num_lambdas*nrow(X))], ncol=num_lambdas, byrow=FALSE)
    solution$lambdas <- sol[["lambdas_out"]][1:num_lambdas]
    solution$loadings <- solution$loadings[,order(solution$lambdas)]
    solution$lambdas <- sort(solution$lambdas)
    solution$lambdas <- solution$lambdas[!duplicated(solution$loadings, MARGIN=2)]
    solution$loadings <- unique(solution$loadings, MARGIN=2)
  } else{
    solution$loadings <- matrix(sol[["loadings"]], ncol=projDim, byrow=FALSE)
    solution$minobjectives <- sol[["objectives"]]
  }
  
  solution <- as.list(solution)
  class(solution) <- "sparsel1pca"
  solution
}
