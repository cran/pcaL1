\name{L2PCA_approx}
\alias{L2PCA_approx}
\alias{l2pcaapprox}
\title{L2PCA_approx}
\description{
Provides an approximation of traditional PCA described by Park and Klabjan (2016) as a subroutine for awl1pca.
}
\usage{
L2PCA_approx(ev.prev, pc.prev, projDim, X.diff)
}
\arguments{
  \item{ev.prev}{matrix of principal component loadings from a previous iteration of awl1pca}
  \item{pc.prev}{vector of eigenvalues from previous iteration of awl1pca}
  \item{projDim}{number of dimensions to project data into, must be an integer}
  \item{X.diff}{The difference between the current weighted matrix estimate and the estimate from the previous iteration}
}
\details{The calculation is performed according to equations (11) and (12) in Park and Klabjan (2016).  The method is an approximation for traditional principal component analysis.}
\value{
'L2PCA_approx' returns a list containing the following components:
\item{eigenvalues}{Estimate of eigenvalues of the covariance matrix.} 
\item{eigenvectors}{Estimate of eigenvectors of the covariance matrix.}
}

\seealso{
  \code{\link{awl1pca}}
}

\references{
Park, Y.W. and Klabjan, D. (2016) Iteratively Reweighted Least Squares Algorithms for L1-Norm Principal Component Analysis, \emph{IEEE International Conference on Data Mining (ICDM)}, 2016.
}
