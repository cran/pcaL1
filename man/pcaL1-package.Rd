\name{pcaL1-package}
\alias{pcaL1-package}
\alias{pcaL1}
\docType{package}
\title{
pcaL1: Three L1-Norm PCA Methods
}
\description{
This package contains three implementations of principal component analysis methods using the L1 norm.  The package depends on COIN-OR Clp version >= 1.12.0.  The methods implemented are PCA-L1 (Kwak 2008), L1-PCA (Ke and Kanade 2003, 2005), and L1-PCA* (Brooks, Dula, and Boone 2013).
}
\details{
\tabular{ll}{
Package: \tab pcaL1\cr
Version: \tab 1.2\cr
Date: \tab 2012-07-05\cr
License: \tab GPL (>=2)\cr
URL: \tab http://www.optimization-online.org/DB_HTML/2012/04/3436.html, http://www.coin-or.org\cr SystemRequirements: \tab COIN-OR Clp (>= 1.12.0)\cr
}

Index:
\preformatted{
l1pca                   L1-PCA
l1pcastar               L1-PCA*
pcaL1-package           pcaL1: Three L1-Norm PCA Methods
pcal1                   PCA-L1
plot.l1pca              Plot an L1pca Object
plot.l1pcastar          Plot an L1pcastar Object
plot.pcal1              Plot a Pcal1 Object
}
}
\author{
Sapan Jot <sapan.madaan@gmail.com> and Paul Brooks <jpbrooks@vcu.edu>.

Maintainer: Paul Brooks <jpbrooks@vcu.edu>
}
\references{
\enumerate{
\item{Kwak N. (2008) Principal Component Analysis Based on L1-Norm Maximization, \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, 30: 1672-1680.}
\item{Ke Q. and Kanade T. (2005) Robust L1 Norm Factorization in the Presence of Outliers and Missing Data by Alternative Convex Programming, \emph{IEEE Conference on Computer Vision and Pattern Recognition}.}
\item{Brooks J.P., Dula J.H., and Boone E.L. (2013) A Pure L1-Norm Princpal Component Analysis, \emph{Computational Statistics \& Dta Analysis}, 61:83-98.}
}
}
\keyword{ package }
