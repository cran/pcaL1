\name{plot.sparsel1pca}
\alias{plot.sparsel1pca}
\title{Plot a Sparsel1pca Object}
\description{Plots the scores on the first two principal components.}
\usage{
\method{plot}{sparsel1pca}(x, \dots)
}
\arguments{
\item{x}{
an object of class \code{sparsel1pca} with scores for at least the first two dimensions}
\item{\dots}{
arguments to be passed to or from other methods.
}
}
\details{
This function is a method for the generic function \code{plot}, for objects of class \code{sparsel1pca}.
}
\seealso{
\code{\link{sparsel1pca}}
}
\examples{
##for a 100x10 data matrix X, 
## lying (mostly) in the subspace defined by the first 2 unit vectors, 
## projects data into 1 dimension.
X <- matrix(c(runif(100*2, -10, 10), rep(0,100*8)),nrow=100) 
               + matrix(c(rep(0,100*2),rnorm(100*8,0,0.1)),ncol=10)
mysparsel1pca <- sparsel1pca(X)

##projects data into 2 dimensions.
mysparsel1pca <- sparsel1pca(X, projDim=2, center=FALSE, projections="l1")

## plot first two scores
plot(mysparsel1pca$scores)
}
\keyword{multivariate}

