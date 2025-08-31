\name{GeoKriglocWeights}
\alias{GeoKriglocWeights}
\encoding{UTF-8}
\title{
Compute kriging weights (and related quantities) for Gaussian / non-Gaussian
spatial, spatio-temporal or multivariate random fields
}
\description{
Given a set of spatial locations (and possibly temporal instants)
the function returns the vector (or matrix) of kriging weights,
the covariance matrix used in the kriging system and,
optionally, the inverse of the left-hand-side matrix
dor a specified neighborhhod
}
\usage{
GeoKriglocWeights(coordx, coordy=NULL, coordz=NULL, coordt=NULL, coordx_dyn=NULL,
 corrmodel, distance="Eucl", grid=FALSE, loc, neighb=NULL,
 maxdist=NULL, maxtime=NULL, method="cholesky", model="Gaussian",
 n=1, nloc=NULL, param, anisopars=NULL,
 radius=1, sparse=FALSE, time=NULL, which=1,
 copula=NULL, X=NULL, Xloc=NULL, Mloc=NULL,parallel=TRUE)
}
\arguments{
  \item{coordx}{Numeric (\eqn{d\times 2}{d x 2}) or (\eqn{d\times 3}{d x 3}) matrix
    of spatial coordinates.  Coordinates on a sphere are accepted
    (lon/lat in decimal degrees) when \code{distance = "Sphere"}.}
  \item{coordy}{Optional numeric vector giving an additional spatial coordinate
    dimension.  Ignored if \code{coordx} is already a matrix.}
  \item{coordz}{Optional numeric vector giving a third spatial coordinate
    dimension.}
  \item{coordt}{Optional numeric vector of temporal coordinates
    (length \eqn{t}{t}).  If missing, a purely spatial random field is assumed.}
  \item{coordx_dyn}{List of \eqn{m} matrices (\eqn{d_t\times 2}{d_t x 2})
    containing time-varying spatial coordinates.  See \code{\link{GeoKrig}}.}
  \item{corrmodel}{Character string naming a valid correlation model.
    See \code{\link{GeoCovmatrix}} for admissible choices.}
  \item{distance}{Character string specifying the spatial distance.
    Default is \code{"Eucl"} (Euclidean).  See \code{\link{GeoFit}}.}
  \item{grid}{Logical.  If \code{TRUE}, coordinates are interpreted as
    defining a regular grid (see \code{\link{GeoKrig}}).}
  \item{loc}{Numeric (\eqn{n\times 2}{n x 2}) matrix of locations for which
    the kriging weights are required.}
 \item{neighb}{Numeric; an optional positive integer indicating the 
    order of the  neighborhood.}
  \item{maxdist}{Numeric; an optional positive value indicating the distance in the spatial neighborhood.}
  \item{maxtime}{Numeric; an optional positive integer value indicating the order of the temporal neighborhood.}
  \item{method}{Character string indicating the matrix factorisation
    used to solve the kriging system: \code{"cholesky"} (default) or \code{"svd"}.}
  \item{model}{Character string specifying the random field type
    (e.g.\ \code{"Gaussian"}, \code{"SkewGaussian"}, \code{"Gamma"}, …).
    See \code{\link{GeoFit}}.}
  \item{n}{Integer.  Number of trials for Binomial random fields
    (default \code{1}).}
  \item{nloc}{Integer.  Number of trials for the prediction locations
    in Binomial random fields (default \code{1}).}
  \item{param}{Named list of covariance and mean parameters.
    See \code{\link{CorrParam}} and \code{\link{GeoCovmatrix}}.}
  \item{anisopars}{List with components \code{angle} and \code{ratio}
    defining geometric anisotropy (optional).}
  \item{radius}{Positive numeric value: sphere radius when
    coordinates are lon/lat (default \code{1}).}
  \item{sparse}{Logical.  If \code{TRUE}, sparse‐matrix algorithms
    (package \pkg{spam}) are employed.  Only effective with compactly
    supported covariance functions.}
  \item{time}{Numeric vector of length \eqn{m}{m} giving the temporal
    instants for which weights are required.  Ignored if \code{coordt} is
    missing.}
  \item{which}{Integer (\code{1} or \code{2}) selecting the variable whose
    weights are returned in the bivariate case.}
  \item{copula}{Character string naming a copula when a non-Gaussian
    dependence structure is used (\code{"Clayton"} or \code{"Gaussian"}).}
  \item{X}{Numeric matrix of spatio-temporal covariates at data locations.}
  \item{Xloc}{Numeric matrix of spatio-temporal covariates at
    prediction locations.}
   \item{Mloc}{Numeric; Vector of spatio(temporal) estimated means 
     associated to predicted locations.}
   \item{parallel}{Logical; if \code{TRUE} then   parallelization is performed 
      }
}
\details{
The function builds the kriging system
\deqn{\Sigma \mathbf{w} = \boldsymbol{\sigma}_0}{Sigma w = sigma_0}
where \eqn{\Sigma}{Sigma} is the covariance matrix between observed
locations and \eqn{\boldsymbol{\sigma}_0}{sigma_0} the vector of
covariances between observed and prediction locations.  
The solution \eqn{\mathbf{w}}{w} is returned together with \eqn{\Sigma}{Sigma}
and, optionally, \eqn{\Sigma^{-1}}{Sigma^-1}.
Universal kriging with covariates is supported; in that case the
generalised least squares estimate of the mean parameters is appended.
No actual prediction is carried out; for full kriging use
\code{\link{GeoKrig}}.
}
\value{
An object of class \code{KgWeights}, a list containing:
  \item{weights}{Numeric vector (or matrix) of kriging weights.}
  \item{covmat}{Covariance matrix \eqn{\Sigma}{Sigma} used in the kriging system.}
  \item{invcov}{Inverse of \eqn{\Sigma}{Sigma} (only if \code{sparse = FALSE}).}
  \item{beta}{Estimated regression coefficients when covariates are supplied.}
  \item{type}{Character string: \code{"simple"} or \code{"universal"}.}
  \item{model}{Input \code{model}.}
  \item{corrmodel}{Input \code{corrmodel}.}
  \item{spacetime}{Logical: \code{TRUE} for spatio-temporal case,
    \code{FALSE} otherwise.}
  \item{sparse}{Logical: was a sparse algorithm used?}
  \item{loc}{Copy of the \code{loc} argument.}
  \item{time}{Copy of the \code{time} argument (if any).}
}
\references{
  Gaetan, C. and Guyon, X. (2010) \emph{Spatial Statistics and Modeling}.
  Springer-Verlag, New York.
}
\seealso{
  \code{\link{GeoKrig}} for full prediction,
   \code{\link{GeoKrigloc}} for local prediction,
  \code{\link{GeoCovmatrix}} for covariance model details,
}
\author{
  Moreno Bevilacqua, \email{moreno.bevilacqua@uai.cl},  
  Víctor Morales-Oñate, \email{victor.morales@uv.cl},  
  Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl}
}
\keyword{Kriging}
