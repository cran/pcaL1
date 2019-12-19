l2projection <- function (X, loadings)
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

  solution <- new.env()
  solution$projPoints <- as.matrix(t(loadings %*% t(loadings) %*% t(X)))
  solution$scores <- as.matrix(X %*% loadings)

  row.names(solution$projPoints) <- rownames(X)
  
  solution <- as.list(solution)
  class(solution) <- "l2projection"
  solution
}
