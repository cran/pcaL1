\name{l1pcastar}
\alias{l1pcastar}
\title{L1-PCA*}
\description{Performs a principal component analysis using the algorithm L1-PCA* described by Brooks, Dula, and Boone (2012)}
\usage{
   l1pcastar(X, projDim=1, center=TRUE, scores=FALSE, 
             projPoints=FALSE, dispExp=FALSE)
}
\arguments{
  \item{X}{data, must be in \code{matrix} or table form}
  \item{projDim}{number of dimensions to project data into, must be an integer, default is 1.}
  \item{center}{whether to center the data using the median, default is FALSE}
  \item{scores}{whether to return the matrix of scores, default is FALSE}
  \item{projPoints}{whether to calculate the projected points (reconstruction), default is FALSE}
  \item{dispExp}{whether to return L1 dispersion explained by each component, default is FALSE}
}
\details{The calculation is performed according to the algorithm described by Brooks, Dula, and Boone (2012).  The algorithm finds successive directions of minimum dispersion in the data by finding the L1-norm best-fit hyperplane at each iteration.  Linear programming instances are solved using Clp (http://www.coin-or.org)}
\value{
'l1pcastar' returns a list with class "l1pcastar" containing the following components:
\item{loadings}{the matrix of variable loadings.  The matrix has dimension ncol(X) x ncol(X).  The columns define the projected subspace.} 
\item{scores}{the matrix of projected points.  The matrix has dimension nrow(X) x projDim.}
\item{dispExp}{the proportion of L1 dispersion explained by the loadings vectors.  Calculated as the L1 dispersion of the score on each component divided by the L1 dispersion in the original data.}
\item{projPoints}{the matrix of projected points in terms of the original coordinates.  The matrix has dimension nrow(X) x ncol(X).}
}
\examples{
##for a 100x10 data matrix X, 
## lying (mostly) in the subspace defined by the first 2 unit vectors, 
## projects data into 1 dimension.
X <- matrix(c(runif(100*2, -10, 10), rep(0,100*8)),nrow=100) 
                + matrix(c(rep(0,100*2),rnorm(100*8,0,0.1)),ncol=10)
myl1pcastar <- l1pcastar(X)

##projects data into 2 dimensions.
myl1pcastar <- l1pcastar(X, projDim=2, center=FALSE, scores=TRUE)

## plot first two scores
plot(myl1pcastar$scores)
}
\references{
Brooks J.P., Dula J.H., and Boone E.L. (2012) A pure L1-norm princpal component analysis, in review, available at http://www.optimization-online.org.
}
\keyword{multivariate}
