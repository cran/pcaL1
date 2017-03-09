\name{pcal1}
\alias{pcal1}
\title{PCA-L1}
\description{
Performs a principal component analysis using the algorithm PCA-L1 given by Kwak (2008).
}
\usage{
pcal1(X, projDim=1, center=TRUE, projections="none", initialize="l2pca")
}
\arguments{
  \item{X}{data, must be in \code{matrix} or table form}
  \item{projDim}{number of dimensions to project data into, must be an integer, default is 1}
  \item{center}{whether to center the data using the median, default is FALSE}
  \item{projections}{whether to calculate reconstructions and scores using the L1 norm ("l1") the L2 norm ("l2") or not at all ("none", default)}
  \item{initialize}{initial guess for first component.  Options are: "l2pca" - use traditional PCA/SVD, "maxx" - use the point with the largest norm, "random" - use a random vector.  The user can also provide a vector as the initial guess.}
}
\details{The calculation is performed according to the algorithm described by Kwak (2008).  The method is a locally-convergent algorithm for finding successive directions of maximum L1 dispersion.}
\value{
'pcal1' returns a list with class "pcal1" containing the following components:
\item{loadings}{the matrix of variable loadings.  The matrix has dimension ncol(X) x projDim.  The columns define the projected subspace.} 
\item{scores}{the matrix of projected points.  The matrix has dimension nrow(X) x projDim.}
\item{dispExp}{the proportion of L1 dispersion explained by the loadings vectors.  Calculated as the L1 dispersion of the score on each component divided by the L1 dispersion in the original data.}
\item{projPoints}{the matrix of projected points in terms of the original coordinates (reconstructions).  The matrix has dimension nrow(X) x ncol(X).}
}

\examples{
##for 100x10 data matrix X, 
## lying (mostly) in the subspace defined by the first 2 unit vectors, 
## projects data into 1 dimension.
X <- matrix(c(runif(100*2, -10, 10), rep(0,100*8)),nrow=100) 
               + matrix(c(rep(0,100*2),rnorm(100*8,0,0.1)),ncol=10)
mypcal1 <- pcal1(X)

##projects data into 2 dimensions.
mypcal1 <- pcal1(X, projDim=2, center=FALSE, projections="l1")

## plot first two scores
plot(mypcal1$scores)
}
\references{
Kwak N. (2008) Principal component analysis based on L1-norm maximization, \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, 30: 1672-1680. DOI:10.1109/TPAMI.2008.114
}
