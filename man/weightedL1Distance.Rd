\name{weightedL1Distance}
\alias{weightedL1Distance}
\title{Weighted L1 Distance}
\description{
  Provides the (weighted) L1-norm distances and total distance of points to a subspace.
}
\usage{
weightedL1Distance(X, loadings, weights)
}
\arguments{
  \item{X}{data, in \code{matrix} or table form}
  \item{loadings}{an orthonormal matrix of loadings vectors}
  \item{weights}{a list of weights for loadings vectors}
}
\details{The reconstructions are calculated by solving a linear program.  Then the weights are applied to the distances.}
\value{
  'weightedL1Distance' returns a list containing the following components:
  \item{wDistances}{list of weighted distances}
  \item{totalDistance}{total distance}
}

