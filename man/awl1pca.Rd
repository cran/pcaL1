\name{awl1pca}
\alias{awl1pca}
\title{awPCA}
\description{
Performs a principal component analysis using the algorithm awPCA described by Park and Klabjan (2016).
}
\usage{
awl1pca(X, projDim=1, center=TRUE, 
         tolerance=0.001, iterations=200, beta=0.99, gamma=0.1)
}
\arguments{
  \item{X}{data, must be in \code{matrix} or table form}
  \item{projDim}{number of dimensions to project data into, must be an integer, default is 1}
  \item{center}{whether to center the data using the mean, default is FALSE}
  \item{tolerance}{for testing convergence; if the sum of absolute values of loadings vectors is smaller, then the algorithm terminates} 
  \item{iterations}{maximum number of iterations in optimization routine}
  \item{beta}{algorithm parameter to set up bound for weights}
  \item{gamma}{algorithm parameter to determine whether to use approximation formula or prcomp function}
}
\details{The calculation is performed according to the algorithm described by Park and Klabjan (2016).  The method is an  iteratively reweighted least squares algorithm for L1-norm principal component analysis.}
\value{
'awl1pca' returns a list with class "awl1pca" containing the following components:
\item{loadings}{the matrix of variable loadings.  The matrix has dimension ncol(X) x projDim.  The columns define the projected subspace.} 
\item{scores}{the matrix of projected points.  The matrix has dimension nrow(X) x projDim.}
\item{projPoints}{the matrix of projected points in terms of the original coordinates.  The matrix has dimension nrow(X) x ncol(X).}
\item{L1error}{sum of the L1 norm of reconstruction errors.}
\item{nIter}{number of iterations.}
\item{ElapsedTime}{elapsed time.}
}

\examples{
##for 100x10 data matrix X, 
## lying (mostly) in the subspace defined by the first 2 unit vectors, 
## projects data into 1 dimension.
X <- matrix(c(runif(100*2, -10, 10), rep(0,100*8)),nrow=100) 
               + matrix(c(rep(0,100*2),rnorm(100*8,0,0.1)),ncol=10)
myawl1pca <- awl1pca(X)

##projects data into 2 dimensions.
myawl1pca <- awl1pca(X, projDim=2, center=FALSE)

## plot first two scores
plot(myawl1pca$scores)
}
\references{
Park, Y.W. and Klabjan, D. (2016) Iteratively Reweighted Least Squares Algorithms for L1-Norm Principal Component Analysis, \emph{IEEE International Conference on Data Mining (ICDM)}, 2016.
}
