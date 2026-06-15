\name{StartParam}
\alias{StartParam}
\encoding{UTF-8}
\title{Initializes the Parameters for Estimation Procedures}
\description{
  Subroutine called by the fitting procedures.
  The procedure initializes the parameters for the fitting procedure.
}
\usage{
StartParam(coordx, coordy, coordz, coordt, coordx_dyn, corrmodel,
           data, distance, fcall, fixed, grid, likelihood,
           maxdist, neighb, maxtime, model, n, param,
           parscale, paramrange, radius, start, taper, tapsep,
           type, typereal, weighted, copula, X, memdist, nosym,
           p_neighb, thin_method)
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

  \item{distance}{String; the name of the spatial distance. Default is \code{"Eucl"} (euclidean distance). See the Section \bold{Details}.}

  \item{fcall}{String; \code{"fitting"} to call the fitting procedure and \code{"simulation"} to call the simulation procedure.}

  \item{fixed}{A named list giving the values of the parameters that will be considered as known values.}

  \item{grid}{Logical; if \code{FALSE} (default) the data are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix;
    if \code{TRUE} then a (\eqn{d \times d \times n}{d x d x n})-array is considered.}

  \item{likelihood}{String; the configuration of the composite likelihood.}

  \item{maxdist}{Numeric; an optional positive value indicating the maximum spatial distance considered in the composite-likelihood computation.}

  \item{neighb}{Numeric; an optional positive integer indicating the order of neighborhood in the composite likelihood computation.
    See the Section \bold{Details} for more information.}

  \item{maxtime}{Numeric; an optional positive value indicating the maximum temporal lag considered in the composite-likelihood computation.}

  \item{model}{String; the density associated to the likelihood objects. \code{"Gaussian"} is the default.}

  \item{n}{Numeric; number of trials for binomial random fields.}

  \item{param}{A numeric vector of parameter values required in the simulation procedure of random fields.}

  \item{parscale}{A numeric vector of scaling factors to improve the maximizing procedure; see \code{\link{optim}}.}

  \item{paramrange}{A numeric vector of parameter ranges; see \code{\link{optim}}.}

  \item{radius}{Numeric; the radius of the sphere in the case of lon-lat coordinates.
    The default is \code{6371}, the radius of the earth.}

  \item{start}{A named list with the initial values of the parameters that are used by the numerical routines in the maximization procedure.}

  \item{taper}{String; the name of the type of covariance matrix. It can be \code{"Standard"} (default) or \code{"Tapering"} for tapered covariance matrices.}

  \item{tapsep}{Numeric; an optional value indicating the separability parameter in the space-time adaptive taper (see \bold{Details}).}

  \item{type}{String; the type of likelihood objects. Temporary value set to \code{"WLeastSquare"} (weighted least-square) in order to compute starting values.}

  \item{typereal}{String; the real type of likelihood objects. See \code{\link{GeoFit}}.}

  \item{weighted}{Logical; if \code{TRUE} the likelihood objects are weighted; see \code{\link{GeoFit}}.}

  \item{copula}{String; the type of copula.}

  \item{X}{Numeric; matrix of space-time covariates.}

  \item{memdist}{Logical; if \code{TRUE} then the distances in the composite likelihood are computed before the optimization.}

  \item{nosym}{Logical; if \code{TRUE} symmetric weights are not considered.}

  \item{p_neighb}{Numeric; a value in \eqn{(0,1]} specifying the expected fraction of nearest-neighbor pairs retained through stochastic thinning.
    If \code{1} (default), no thinning is applied and all nearest-neighbor pairs are used. If \code{<1}, pairs are randomly retained using independent Bernoulli sampling.}

  \item{thin_method}{String; thinning scheme used when \code{p_neighb<1}. Default is \code{"bernoulli"} (independent Bernoulli thinning).}
}
\value{
  A list with a set of useful information in the estimation procedure.
}
\details{
  Internal function called by \code{\link{WlsStart}}.
}
\author{
Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com}, \url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl}, \url{https://www.researchgate.net/profile/Christian-Caamano}
}
\keyword{Composite}
