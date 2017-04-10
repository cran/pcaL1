\name{sharpel1pca}
\alias{sharpel1pca}
\title{SharpEl1-PCA}
\description{Performs a principal component analysis using the algorithm SharpEl1-PCA described by Brooks and Dula (2017, submitted)}
\usage{
   sharpel1pca(X, projDim=1, center=TRUE, projections="none")
}
\arguments{
  \item{X}{data, must be in \code{matrix} or table form.}
  \item{projDim}{number of dimensions to project data into, must be an integer, default is 1.}
  \item{center}{whether to center the data using the median, default is TRUE.}
  \item{projections}{whether to calculate reconstructions and scores using the L1 norm ("l1") the L2 norm ("l2") or not at all ("none", default).}
}
\details{The calculation is performed according to the algorithm described by Brooks and Dula (2017, submitted).  The algorithm finds successive, orthogonal fitted lines in the data.}
\value{
'sharpel1pca' returns a list with class "sharpel1pca" containing the following components:
\item{loadings}{the matrix of variable loadings.  The matrix has dimension ncol(X) x projDim.  The columns define the projected subspace.} 
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
mysharpel1pca <- sharpel1pca(X)

##projects data into 2 dimensions.
mysharpel1pca <- sharpel1pca(X, projDim=2, center=FALSE, projections="l1")

## plot first two scores
plot(mysharpel1pca$scores)
}
\references{
Brooks J.P. and Dula J.H. (2017) Estimating L1-Norm Best-Fit Lines, submitted.
}
\keyword{multivariate}
