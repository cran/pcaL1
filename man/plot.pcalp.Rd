\name{plot.pcalp}
\alias{plot.pcalp}
\title{Plot a Pcalp Object}
\description{Plots the scores on the first two principal components.}
\usage{
\method{plot}{pcalp}(x, \dots)
}
\arguments{
\item{x}{
an object of class \code{pcalp} with scores for at least the first two dimensions}
\item{\dots}{
arguments to be passed to or from other methods.
}
}
\details{
This function is a method for the generic function \code{plot}, for objects of class \code{pcalp}.
}
\seealso{
\code{\link{pcalp}}
}
\examples{
##for a 100x10 data matrix X, 
## lying (mostly) in the subspace defined by the first 2 unit vectors, 
## projects data into 1 dimension.
X <- matrix(c(runif(100*2, -10, 10), rep(0,100*8)),nrow=100) 
               + matrix(c(rep(0,100*2),rnorm(100*8,0,0.1)),ncol=10)
mypcalp <- pcalp(X)

##projects data into 2 dimensions.
mypcalp <- pcalp(X, projDim=2, center=FALSE, projections="l1")

## plot first two scores
plot(mypcalp$scores)
}
\keyword{multivariate}
