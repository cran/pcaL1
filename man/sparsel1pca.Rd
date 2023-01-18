\name{sparsel1pca}
\alias{sparsel1pca}
\title{SparsEl1-PCA}
\description{L1-norm line fitting with L1-regularization.}
\usage{
   sparsel1pca(X, projDim=1, center=TRUE, projections="none", lambda=0)
}
\arguments{
  \item{X}{data, must be in \code{matrix} or table form.}
  \item{projDim}{number of dimensions to project data into, must be an integer, default is 1.}
  \item{center}{whether to center the data using the median, default is TRUE.}
  \item{projections}{whether to calculate reconstructions and scores using the L1 norm ("l1") the L2 norm ("l2") or not at all ("none", default).}
  \item{lambda}{If negative and number of rows is at most 100, calculates all possible breakpoints for the regularization parameter. Otherwise, fits a regularlized line with lambda set to that value.}
}
\details{The calculation is performed according to the algorithm described by Ling and Brooks (2023, working paper).  The algorithm finds successive, orthogonal fitted lines in the data.}
\value{
'sparsel1pca' returns a list with class "sparsel1pca" containing the following components:
\item{loadings}{the matrix of variable loadings.  The matrix has dimension ncol(X) x projDim.  The columns define the projected subspace.} 
\item{scores}{the matrix of projected points.  The matrix has dimension nrow(X) x projDim.}
\item{dispExp}{the proportion of L1 dispersion explained by the loadings vectors.  Calculated as the L1 dispersion of the score on each component divided by the L1 dispersion in the original data.}
\item{projPoints}{the matrix of projected points in terms of the original coordinates.  The matrix has dimension nrow(X) x ncol(X).}
\item{minobjectives}{the L1 distance of points to their projections in the fitted subspace.}
}
\examples{
##for a 100x10 data matrix X, 
## lying (mostly) in the subspace defined by the first 2 unit vectors, 
## projects data into 1 dimension.
X <- matrix(c(runif(100*2, -10, 10), rep(0,100*8)),nrow=100) +
                matrix(c(rep(0,100*2),rnorm(100*8,0,0.1)),ncol=10)
mysparsel1pca <- sparsel1pca(X, lambda=0.5)

##projects data into 2 dimensions.
mysparsel1pca <- sparsel1pca(X, projDim=2, center=FALSE, projections="l1", lambda=0.5)

## plot first two scores
plot(mysparsel1pca$scores)
}
\references{
Ling, X. and Brooks J.P. (2023) L1-Norm Regularized L1-Norm Best-Fit Lines, working paper.
}
\keyword{multivariate}
