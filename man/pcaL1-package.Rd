\name{pcaL1-package}
\alias{pcaL1-package}
\alias{pcaL1}
\docType{package}
\title{
pcaL1: L1-Norm PCA Methods
}
\description{
This package contains implementations of six principal component analysis methods using the L1 norm.  The package depends on COIN-OR Clp version >= 1.12.0.  The methods implemented are PCA-L1 (Kwak 2008), L1-PCA (Ke and Kanade 2003, 2005), L1-PCA* (Brooks, Dula, and Boone 2013), L1-PCAhp (Visentin, Prestwich and Armagan 2016), wPCA (Park and Klabjan 2016), and awPCA (Park and Klabjan 2016).
}
\details{
\tabular{ll}{
Package: \tab pcaL1\cr
Version: \tab 1.5.4\cr
Date: \tab 2020-02-10\cr
License: \tab GPL (>=3)\cr
URL: \tab http://www.optimization-online.org/DB_HTML/2012/04/3436.html, http://www.coin-or.org\cr SystemRequirements: \tab COIN-OR Clp (>= 1.12.0)\cr
}

Index:
\preformatted{
awl1pca                 awPCA 
l1pca                   L1-PCA
l1pcahp                 L1-PCAhp
l1pcastar               L1-PCA*
l1projection            L1-Norm Projection on a Subspace
L2PCA_approx            Subroutine for awl1pca
l2projection            L2-Norm Projection on a Subspace
pcal1                   PCA-L1
pcalp                   PCA-Lp
pcaL1-package           pcaL1: L1-Norm PCA Methods
plot.awl1pca            Plot an awl1pca Object
plot.l1pca              Plot an l1pca Object
plot.l1pcahp            Plot an l1pcahp Object
plot.l1pcastar          Plot an l1pcastar Object
plot.pcal1              Plot a pcal1 Object
plot.pcalp              Plot a pcalp Object
plot.wl1pca             Plot an wl1pca Object
plot.sharpel1pca        Plot a sharpel1pca Object
sharpel1pca             SharpeEL1-PCA
wl1pca                  wPCA 
}
}
\author{
Sapan Jot <sapan.madaan@gmail.com>, Paul Brooks <jpbrooks@vcu.edu>, Andrea Visentin <andrea.visentin@insight-centre.org>,Young Woong Park <ywpark@mail.smu.edu>, and Yi-Hui Zhou <yihui_zhou@ncsu.edu>

Maintainer: Paul Brooks <jpbrooks@vcu.edu>
}
\references{
\enumerate{
\item{Brooks and Dula (2017) Estimating L1-Norm Best-Fit Lines, submitted}
\item{Brooks J.P., Dula J.H., and Boone E.L. (2013) A Pure L1-Norm Princpal Component Analysis, \emph{Computational Statistics \& Data Analysis}, 61:83-98. DOI:10.1016/j.csda.2012.11.007}
\item{Ke Q. and Kanade T. (2005) Robust L1 Norm Factorization in the Presence of Outliers and Missing Data by Alternative Convex Programming, \emph{IEEE Conference on Computer Vision and Pattern Recognition}. DOI:10.1109/CVPR.2005.309}
\item{Kwak N. (2008) Principal Component Analysis Based on L1-Norm Maximization, \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, 30: 1672-1680. DOI:10.1109/TPAMI.2008.114}
\item{Kwak N. (2014) Principal Component Analysis by  Lp-Norm Maximization, \emph{IEEE Transactions on Cybernetics}, 44:594-609. DOI:10.1109/TCYB.2013.2262936}
\item{Park, Y.W. and Klabjan, D. (2016) Iteratively Reweighted Least Squares Algorithms for L1-Norm Principal Component Analysis, \emph{IEEE International Conference on Data Mining (ICDM)}. DOI: 10.1109/ICDM.2016.0054}
\item{Visentin A., Prestwich S., and Armagan S. T. (2016) Robust Principal Component Analysis by Reverse Iterative Linear Programming, \emph{Joint European Conference on Machine Learning and Knowledge Discovery in Databases},  593-605. DOI:10.1007/978-3-319-46227-1_37}
\item{Zhou, Y.-H. and Marron, J.S. (2016) Visualization of Robust L1PCA, \emph{Stat}, 5:173-184.  DOI:10.1002/sta4.113}
}
}
\keyword{ package }
