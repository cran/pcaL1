\name{l1projection}
\alias{l1projection}
\title{L1 Projection}
\description{
  Provides the L1-norm projection of points on a subspace, including both scores and reconstructions.
}
\usage{
l1projection(X, loadings)
}
\arguments{
  \item{X}{data, in \code{matrix} or table form}
  \item{loadings}{an orthonormal matrix of loadings vectors}
}
\details{The scores and reconstructions are calculated by solving a linear program.}
\value{
  'l1projection' returns a list containing the following components:
  \item{scores}{the matrix of projected points}
  \item{projPoints}{the matrix of projected points in terms of the original coordinates (reconstructions)}
}

