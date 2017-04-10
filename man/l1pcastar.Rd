\name{l1pcastar}
\alias{l1pcastar}
\title{L1-PCA*}
\description{Performs a principal component analysis using the algorithm L1-PCA* described by Brooks, Dula, and Boone (2013)}
\usage{
   l1pcastar(X, projDim=1, center=TRUE, projections="none")
}
\arguments{
  \item{X}{data, must be in \code{matrix} or table form}
  \item{projDim}{number of dimensions to project data into, must be an integer, default is 1}
  \item{center}{whether to center the data using the median, default is TRUE}
  \item{projections}{whether to calculate reconstructions and scores using the L1 norm ("l1") the L2 norm ("l2") or not at all ("none", default)}
}
\details{The calculation is performed according to the algorithm described by Brooks, Dula, and Boone (2013).  The algorithm finds successive directions of minimum dispersion in the data by finding the L1-norm best-fit hyperplane at each iteration.  Linear programming instances are solved using Clp (http://www.coin-or.org)}
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
X <- matrix(c(runif(100*2, -10, 10), rep(0,100*8)),nrow=100) +
                 matrix(c(rep(0,100*2),rnorm(100*8,0,0.1)),ncol=10)
myl1pcastar <- l1pcastar(X)

##projects data into 2 dimensions.
myl1pcastar <- l1pcastar(X, projDim=2, center=FALSE, projections="l1")

## plot first two scores
plot(myl1pcastar$scores)
}
\references{
  \enumerate{
  \item{Brooks J.P., Dula J.H., and Boone E.L. (2013) A Pure L1-Norm Princpal Component Analysis, \emph{Computational Statistics \& Data Analysis}, 61:83-98. DOI:10.1016/j.csda.2012.11.007}
  \item{Zhou, Y.-H. and Marron, J.S. (2016) Visualization of Robust L1PCA, \emph{Stat}, 5:173-184.  DOI:10.1002/sta4.113}
  }
}
\keyword{multivariate}
