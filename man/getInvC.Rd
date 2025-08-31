\name{getInvC}
\alias{getInvC}
\encoding{UTF-8}
\title{Utility function for kriging computation}
\description{Utility function that computes
  (\eqn{\Sigma^{-1} \mathbf{CC}})   where \eqn{\Sigma} is a covariance matrix and 
  \eqn{\mathbf{CC}} is the matrix of of covariances between the observed location and the locations to predict.
  Optionally, the quadratic form
  \eqn{\mathbf{cc}^\top \Sigma^{-1} \mathbf{cc}} can be computed. Both dense and sparse matrix representations are supported.
}
\usage{
getInvC(covmatrix, CC, mse = TRUE)
}
\arguments{
  \item{covmatrix}{A covariance matrix:
    \describe{
      \item{\code{covmatrix}}{A covariance matrix (dense or sparse).}
      \item{\code{sparse}}{Logical scalar indicating whether the matrix
        is stored in sparse format.}
    }
  }
  \item{CC}{Numeric matrix to be multiplied by the inverse of the
    covariance matrix (length must match the dimension of
    \code{covmatrix}).}
  \item{mse}{Logical.  If \code{TRUE} (default) the quadratic form
    \eqn{\mathbf{cc}^\top \Sigma^{-1} \mathbf{cc}} is also returned.}
}
\value{
  A list with components
  \describe{
    \item{\code{a}}{Numeric vector: the product
      \eqn{\Sigma^{-1} \mathbf{cc}}.}
    \item{\code{b}}{Numeric scalar: the quadratic form
      \eqn{\mathbf{cc}^\top \Sigma^{-1} \mathbf{cc}} (only if
      \code{mse = TRUE}).}
  }
}
\details{
  For dense matrices the function uses the Cholesky decomposition provided
  by \code{FastGP::rcppeigen_get_chol}.  
  For sparse matrices (class \code{spam}) the factorisation is performed
  with \code{spam::chol.spam}.  
  If the covariance matrix is not positive definite an error is thrown.
}
\author{
  Moreno Bevilacqua, \email{moreno.bevilacqua@uai.cl}
}

\keyword{Kriging}
