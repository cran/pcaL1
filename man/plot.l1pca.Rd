\name{plot.l1pca}
\alias{plot.l1pca}
\title{Plot an L1pca Object}
\description{Plots the scores on the first two principal components.}
\usage{
\method{plot}{l1pca}(x, \dots)
}
\arguments{
\item{x}{
an object of class \code{l1pca} with scores for at least the first two dimensions}
\item{\dots}{
arguments to be passed to or from other methods.
}
}
\details{
This function is a method for the generic function \code{plot}, for objects of class \code{l1pca}.
}
\seealso{
\code{\link{l1pca}}
}
\examples{
##for a 100x10 data matrix X, 
## lying (mostly) in the subspace defined by the first 2 unit vectors, 
## projects data into 1 dimension.
X <- matrix(c(runif(100*2, -10, 10), rep(0,100*8)),nrow=100) 
              + matrix(c(rep(0,100*2),rnorm(100*8,0,0.1)),ncol=10)
myl1pca <- l1pca(X)

##projects data into 2 dimensions.
myl1pca <- l1pca(X, projDim=2, center=FALSE)

## plot first two scores
plot(myl1pca$scores)
}
\keyword{multivariate}
