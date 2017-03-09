\name{l2projection}
\alias{l2projection}
\title{L2 Projection}
\description{
  Provides the L2-norm projection of points on a subspace, including both scores and reconstructions.
}
\usage{
l2projection(X, loadings)
}
\arguments{
  \item{X}{data, in \code{matrix} or table form}
  \item{loadings}{an orthonormal matrix of loadings vectors}
}
\details{The scores and reconstructions are calculated by solving a linear program.}
\value{
  'l2projection' returns a list containing the following components:
  \item{scores}{the matrix of projected points}
  \item{projPoints}{the matrix of projected points in terms of the original coordinates (reconstructions)}
}

