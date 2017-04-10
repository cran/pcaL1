\name{l1pca}
\alias{l1pca}
\title{L1-PCA}
\description{Performs a principal component analysis using the algorithm L1-PCA given by Ke and Kanade (2005).}
\usage{
   l1pca(X, projDim=1, center=TRUE, projections="l1", 
	 initialize="l2pca", tolerance=0.0001, iterations=10)
}
\arguments{
  \item{X}{data, must be in \code{matrix} or table form.}
  \item{projDim}{number of dimensions to project data into, must be an integer, default is 1.}
  \item{center}{whether to center the data using the median, default is TRUE.}
  \item{projections}{Whether to calculate reconstructions and scores using the L1 ("l1", default) or L2 ("l2") norm.}
  \item{initialize}{initial guess for loadings matrix.  Options are: "l2pca" - use traditional PCA/SVD, "random" - use a randomly-generated matrix.  The user can also provide a matrix as an initial guess.}
  \item{tolerance}{sets the convergence tolerance for the algorithm, default is 0.0001.}
  \item{iterations}{sets the number of iterations to run before returning the result, default is 10.}
}
\details{The calculation is performed according to the linear programming-based algorithm described by Ke and Kanade (2005).  The method is a locally-convergent algorithm for finding the L1-norm best-fit subspace by alternatively optimizing the scores and the loadings matrix at each iteration.  Linear programming instances are solved using Clp (http://www.coin-or.org)
}
\value{
'l1pca' returns a list with class "l1pca" containing the following components:
\item{loadings}{the matrix of variable loadings.  The matrix has dimension ncol(X) x projDim.  The columns defined the projected subspace.} 
\item{scores}{the matrix of projected points.  The matrix has dimension nrow(X) x projDim.}
\item{dispExp}{the proportion of L1 dispersion explained by the loadings vectors.  Calculated as the L1 dispersion of the score on each component divided by the L1 dispersion in the original data.}
\item{projPoints}{the matrix of projected points in terms of the original coordinates (reconstructions).  The matrix has dimension nrow(X) x ncol(X).}
}

\examples{
##for 100x10 data matrix X, 
## lying (mostly) in the subspace defined by the first 2 unit vectors, 
## projects data into 1 dimension.
X <- matrix(c(runif(100*2, -10, 10), rep(0,100*8)),nrow=100) + 
                 matrix(c(rep(0,100*2),rnorm(100*8,0,0.1)),ncol=10)
myl1pca <- l1pca(X)

##projects data into 2 dimensions.
myl1pca <- l1pca(X, projDim=2, center=FALSE, 
                 tolerance=0.00001, iterations=20)

## plot first two scores
plot(myl1pca$scores)
}

\references{
Ke Q. and Kanade T. (2005) Robust L1 norm factorization in the presence of outliers and missing data by alternative convex programming, \emph{IEEE Conference on Computer Vision and Pattern Recognition}. DOI:10.1109/CVPR.2005.309
}
\keyword{multivariate}
