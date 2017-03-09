\name{pcalp}
\alias{pcalp}
\title{PCA-Lp}
\description{
Performs a principal component analysis using the greedy algorithms PCA-Lp(G) and PCA-Lp(L) given by Kwak (2014). 
}
\usage{
pcalp(X, projDim=1, p = 1.0, center=TRUE, projections="none", 
        initialize="l2pca",solution = "L", 
	epsilon = 0.0000000001, lratio = 0.02)
}
\arguments{
  \item{X}{data, must be in \code{matrix} or table form.}
  \item{projDim}{number of dimensions to project data into, must be an integer, default is 1.}  
  \item{p}{p-norm use to measure the distance between points.}
  \item{center}{whether to center the data using the median, default is FALSE.}
  \item{projections}{whether to calculate reconstructions and scores using the L1 norm ("l1") the L2 norm ("l2") or not at all ("none", default)}
  \item{initialize}{method for initial guess for component.  Options are: "l2pca" - use traditional PCA/SVD, "maxx" - use the point with the largest norm, "random" - use a random vector.}  
  \item{solution}{method projection vector update.  Options are: "G" - PCA-Lp(G) implementation: Gradient search, "L" - PCA-Lp(L) implementation: Lagrangian (default).}
  \item{epsilon}{for checking convergence}
  \item{lratio}{learning ratio, default is 0.02. Suggested value 1/(nr. istances).}
}

\details{The calculation is performed according to the algorithm described by Kwak (2014), an extension of the original Kwak(2008).  The method is a greedy locally-convergent algorithm for finding successive directions of maximum Lp dispersion.}
\value{
'pcalp' returns a list with class "pcalp" containing the following components:
\item{loadings}{the matrix of variable loadings.  The matrix has dimension ncol(X) x projDim.  The columns define the projected subspace.} 
\item{scores}{the matrix of projected points.  The matrix has dimension nrow(X) x projDim.}
\item{dispExp}{the proportion of L1 dispersion explained by the loadings vectors.  Calculated as the L1 dispersion of the score on each component divided by the L1 dispersion in the original data.}
\item{projPoints}{the matrix of projected points in terms of the original coordinates.  The matrix has dimension nrow(X) x ncol(X).}
}

\examples{  
##for 100x10 data matrix X, 
## lying (mostly) in the subspace defined by the first 2 unit vectors, 
## projects data into 1 dimension.
X <- matrix(c(runif(100*2, -10, 10), rep(0,100*8)),nrow=100) 
               + matrix(c(rep(0,100*2),rnorm(100*8,0,0.1)),ncol=10)
mypcalp <- pcalp(X, p = 1.5)

##projects data into 2 dimensions.
mypcalp <- pcalp(X, projDim=2, p = 1.5, center=FALSE, projections="l1")

## plot first two scores
plot(mypcalp$scores)
}
\references{
Kwak N. (2008) Principal component analysis based on L1-norm maximization, \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, 30: 1672-1680. DOI:10.1109/TPAMI.2008.114

Kwak N. (2014). Principal component analysis by Lp-norm maximization. \emph{IEEE transactions on cybernetics, 44(5)}, 594-609. DOI: 10.1109/TCYB.2013.2262936

}
