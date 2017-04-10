\name{wl1pca}
\alias{wl1pca}
\title{wPCA}
\description{
Performs a principal component analysis using the algorithm wPCA described by Park and Klabjan (2016).
}
\usage{
wl1pca(X, projDim=1, center=TRUE, projections="l2",
         tolerance=0.001, iterations=200, beta=0.99)
}
\arguments{
  \item{X}{data, must be in \code{matrix} or table form.}
  \item{projDim}{number of dimensions to project data into, must be an integer, default is 1.}
  \item{center}{whether to center the data using the mean, default is TRUE}
  \item{projections}{whether to calculate projections (reconstructions and scores) using the L2 norm ("l2", default) or the L1 norm ("l1").}
  \item{tolerance}{for testing convergence; if the sum of absolute values of loadings vectors is smaller, then the algorithm terminates.} 
  \item{iterations}{maximum number of iterations in optimization routine.}
  \item{beta}{algorithm parameter to set up bound for weights.}
}
\details{The calculation is performed according to the algorithm described by Park and Klabjan (2016).  The method is an  iteratively reweighted least squares algorithm for L1-norm principal component analysis.}
\value{
'wl1pca' returns a list with class "wl1pca" containing the following components:
\item{loadings}{the matrix of variable loadings.  The matrix has dimension ncol(X) x projDim.  The columns define the projected subspace.} 
\item{scores}{the matrix of projected points.  The matrix has dimension nrow(X) x projDim.}
\item{projPoints}{the matrix of L2 projections points on the fitted subspace in terms of the original coordinates.  The matrix has dimension nrow(X) x ncol(X).}
\item{L1error}{sum of the L1 norm of reconstruction errors.}
\item{nIter}{number of iterations.}
\item{ElapsedTime}{elapsed time.}
}

\examples{
##for 100x10 data matrix X, 
## lying (mostly) in the subspace defined by the first 2 unit vectors, 
## projects data into 1 dimension.
X <- matrix(c(runif(100*2, -10, 10), rep(0,100*8)),nrow=100) +
               matrix(c(rep(0,100*2),rnorm(100*8,0,0.1)),ncol=10)
mywl1pca <- wl1pca(X)

##projects data into 2 dimensions.
mywl1pca <- wl1pca(X, projDim=2, center=FALSE)

## plot first two scores
plot(mywl1pca$scores)
}
\references{
Park, Y.W. and Klabjan, D. (2016) Iteratively Reweighted Least Squares Algorithms for L1-Norm Principal Component Analysis, \emph{IEEE International Conference on Data Mining (ICDM)}, 2016.
}
