\name{plot.pcal1}
\alias{plot.pcal1}
\title{Plot a Pcal1 Object}
\description{Plots the scores on the first two principal components.}
\usage{
\method{plot}{pcal1}(x, \dots)
}
\arguments{
\item{x}{
an object of class \code{pcal1} with scores for at least the first two dimensions}
\item{\dots}{
arguments to be passed to or from other methods.
}
}
\details{
This function is a method for the generic function \code{plot}, for objects of class \code{pcal1}.
}
\seealso{
\code{\link{pcal1}}
}
\examples{
##for a 100x10 data matrix X, 
## lying (mostly) in the subspace defined by the first 2 unit vectors, 
## projects data into 1 dimension.
X <- matrix(c(runif(100*2, -10, 10), rep(0,100*8)),nrow=100) 
               + matrix(c(rep(0,100*2),rnorm(100*8,0,0.1)),ncol=10)
mypcal1 <- pcal1(X)

##projects data into 2 dimensions.
mypcal1 <- pcal1(X, projDim=2, center=FALSE, scores=TRUE)

## plot first two scores
plot(mypcal1$scores)
}
\keyword{multivariate}
