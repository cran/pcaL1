\name{pcaL1-package}
\alias{pcaL1-package}
\alias{pcaL1}
\docType{package}
\title{
pcaL1: Three L1-Norm PCA Methods
}
\description{
This package contains three implementations of principal component analysis methods using the L1 norm.  The package depends on COIN-OR Clp version >= 1.12.0.  The methods implemented are PCA-L1 (Kwak 2008), L1-PCA (Ke and Kanade 2003, 2005), and L1-PCA* (Brooks, Dula, and Boone 2012).
}
\details{
\tabular{ll}{
Package: \tab pcaL1\cr
Version: \tab 1.0\cr
Date: \tab 2012-04-20\cr
License: \tab GPL (>=2)\cr
URL: \tab http://www.optimization-online.org/DB_HTML/2012/04/3436.html, http://www.coin-or.org\cr SystemRequirements: \tab COIN-OR Clp (>= 1.12.0)\cr
Packaged: \tab 2012-04-20 05:20:06 UTC; root\cr
Built: \tab R 2.12.1; x86_64-pc-linux-gnu; 2012-04-20 05:20:18 UTC; unix\cr
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
Sapan Jot <sapan.madaan@gmail.com> and Paul Brooks <jpbrooks@vcu.edu>. The scripts configure and configure.ac and the file INSTALL are based on those in R package clpAPI <geliudie@uni-duesseldorf.de>.

Maintainer: Paul Brooks <jpbrooks@vcu.edu>
}
\references{
\enumerate{
\item{Kwak N. (2008) Principal component analysis based on L1-norm maximization, \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, 30: 1672-1680.}
\item{Ke Q. and Kanade T. (2005) Robust L1 norm factorization in the presence of outliers and missing data by alternative convex programming, \emph{IEEE Conference on Computer Vision and Pattern Recognition}.}
\item{Brooks J.P., Dula J.H., and Boone E.L. (2012) A pure L1-norm princpal component analysis, in review, available at http://www.optimization-online.org.}
}
}
\keyword{ package }
