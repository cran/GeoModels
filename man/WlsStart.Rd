\name{WlsStart}
\alias{WlsStart}
\encoding{UTF-8}
\title{Computes Starting Values based on Weighted Least Squares}
\description{
 Subroutine called by \code{\link{GeoFit}}.
 The function returns appropriate starting values for the
 composite-likelihood fitting procedure based on weighted least squares.
}
\usage{
WlsStart(coordx, coordy, coordz, coordt, coordx_dyn, corrmodel, data,
 distance, fcall, fixed, grid, likelihood, maxdist, neighb,
 maxtime, model, n, param, parscale, paramrange, radius, start,
 taper, tapsep, type, varest, weighted, copula, X, memdist,
 nosym, p_neighb, thin_method)
}
\arguments{
 \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix or (\eqn{d \times 3}{d x 3})-matrix.
 Coordinates on a sphere for a fixed radius \code{radius} are passed in lon/lat format expressed in decimal degrees.}

 \item{coordy}{A numeric vector giving 1-dimension of spatial coordinates; optional argument, default is \code{NULL}.}

 \item{coordz}{A numeric vector giving 1-dimension of spatial coordinates; optional argument, default is \code{NULL}.}

 \item{coordt}{A numeric vector assigning 1-dimension of temporal coordinates.}

 \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d_t x 2})-matrices containing dynamical (in time) spatial coordinates.
 Optional argument, default is \code{NULL}.}

 \item{corrmodel}{String; the name of a correlation model.}

 \item{data}{A numeric vector or a (\eqn{n \times d}{n x d})-matrix or (\eqn{d \times d \times n}{d x d x n})-array of observations.}

 \item{distance}{String; the name of the spatial distance. Default is \code{"Eucl"} (Euclidean distance).}

 \item{fcall}{String; \code{"fitting"} to call the fitting procedure and \code{"simulation"} to call the simulation procedure.}

 \item{fixed}{A named list giving the values of the parameters that will be considered as known values.}

 \item{grid}{Logical; if \code{FALSE} (default) the data are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix;
 if \code{TRUE} then a (\eqn{d \times d \times n}{d x d x n})-array is considered.}

 \item{likelihood}{String; the configuration of the composite likelihood.}

 \item{maxdist}{Numeric; an optional positive value indicating the maximum spatial distance considered in the composite-likelihood computation.}

 \item{neighb}{Numeric; an optional positive integer indicating the order of neighborhood in the composite likelihood computation.}

 \item{maxtime}{Numeric; an optional positive value indicating the maximum temporal separation considered in the composite-likelihood computation.}

 \item{model}{String; the name of the model.}

 \item{n}{Numeric; number of trials in a binomial random field.}

 \item{param}{A numeric vector of parameter values required in the simulation procedure of random fields.}

 \item{parscale}{A numeric vector with scaling values for improving the maximisation routine.}

 \item{paramrange}{A numeric vector with the range of the parameter space.}

 \item{radius}{Numeric; a value indicating the radius of the sphere when using great circle distance.
 Default value is the radius of the earth in km (i.e., \code{6371}).}

 \item{start}{A numeric vector (or list) with starting values.}

 \item{taper}{String; the name of the type of covariance matrix. It can be \code{"Standard"} (default) or
 \code{"Tapering"} for tapered covariance matrices.}

 \item{tapsep}{Numeric; an optional value indicating the separability parameter in the space-time quasi taper (see Details).}

 \item{type}{String; the type of estimation method.}

 \item{varest}{Logical; if \code{TRUE} the estimates variances and standard errors are returned. Default is \code{FALSE}.}

 \item{weighted}{Logical; if \code{TRUE} the likelihood objects are weighted; see \code{\link{GeoFit}}.}

 \item{copula}{String; the type of copula. It can be \code{"Clayton"} or \code{"Gaussian"}.}

 \item{X}{Numeric; matrix of spatio(temporal) covariates in the linear mean specification.}

 \item{memdist}{Logical; if \code{TRUE} then the distances in the composite likelihood are computed before the optimization.}

 \item{nosym}{Logical; if \code{TRUE} symmetric weights are not considered.}

 \item{p_neighb}{Numeric; a value in \eqn{(0,1]} specifying the expected fraction of nearest-neighbor pairs retained through stochastic thinning.
 If \code{1} (default), no thinning is applied and all nearest-neighbor pairs are used. If \code{<1}, pairs are randomly retained using
 independent Bernoulli sampling.}

 \item{thin_method}{String; thinning scheme used when \code{p_neighb<1}. Default is \code{"bernoulli"} (independent Bernoulli thinning).}
}
\value{
 A list with a set of useful information in the estimation procedure.
}
\details{
 Internal function called by \code{\link{GeoFit}}.
}
\seealso{
 \code{\link{GeoFit}}.
}
\author{
Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com}, \url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl}, \url{https://www.researchgate.net/profile/Christian-Caamano}
}
\keyword{LeastSquare}
