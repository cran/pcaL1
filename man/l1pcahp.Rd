\name{l1pcahp}
\alias{l1pcahp}
\title{L1-PCAhp}
\description{Performs a principal component analysis using the algorithm L1-PCAhp described by Visentin, Prestwich and Armagan (2016)}
\usage{
   l1pcahp(X, projDim=1, center=TRUE, projPoints=FALSE, scores=TRUE, 
           initialize="l2pca", threshold=0.0001)
}
\arguments{
  \item{X}{data, must be in \code{matrix} or table form}
  \item{projDim}{number of dimensions to project data into, must be an integer, default is 1}
  \item{center}{whether to center the data using the median, default is FALSE}
  \item{scores}{whether to return the matrix of scores, default is FALSE}
  \item{projPoints}{whether to calculate the projected points (reconstruction), default is FALSE}
  \item{initialize}{initial guess for loadings matrix.  Options are: "l2pca" - use traditional PCA/SVD, "random" - use a randomly-generated matrix.  The user can also specify a matrix as an initial guess.}
  \item{threshold}{sets the convergence threshold for the algorithm, default is 0.001}}


\details{The calculation is performed according to the algorithm described by Visentin, Prestwich and Armagan (2016).  The algorithm computes components iteratively in reverse, using a new heuristic based on Linear Programming.  Linear programming instances are solved using Clp (http://www.coin-or.org)}
\value{
'l1pcahp' returns a list with class "l1pcahp" containing the following components:
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
myl1pcahp <- l1pcahp(X)

##projects data into 2 dimensions.
myl1pcahp <- l1pcahp(X, projDim=2, center=FALSE, scores=TRUE)

## plot first two scores
plot(myl1pcahp$scores)
}
\references{
Visentin A., Prestwich S., and Armagan S. T. (2016) Robust Principal Component Analysis by Reverse Iterative Linear Programming, \emph{Joint European Conference on Machine Learning and Knowledge Discovery in Databases},  593-605. DOI:10.1007/978-3-319-46227-1_37
}
