\name{CompLik}
\alias{CompLik}
\encoding{UTF-8}
\title{Optimizes the Composite Log-likelihood}
\description{
  Subroutine called by \code{GeoFit}.
  The procedure estimates the model parameters by maximization of the
  composite log-likelihood.
}
\usage{
CompLik(copula, bivariate, coordx, coordy, coordz, coordt, 
        coordx_dyn, corrmodel, data, distance, flagcorr, 
        flagnuis, fixed, GPU, grid, likelihood, local, lower, 
        model, n, namescorr, namesnuis, namesparam,
        numparam, numparamcorr, optimizer, 
        onlyvar, parallel, param,
        spacetime, type, upper, varest,
        weigthed, ns, X, sensitivity, MM, aniso)
}
\arguments{
  \item{copula}{String; the type of copula. It can be "Clayton" or "Gaussian".}
  \item{bivariate}{Logical; if \code{TRUE} then the data come from a bivariate random field,
    otherwise from a univariate random field.}
  \item{coordx}{A numeric \eqn{d \times 2} or \eqn{d \times 3} matrix.
    Coordinates on a sphere for a fixed radius \code{radius} 
    are passed in lon/lat format expressed in decimal degrees.}
  \item{coordy}{A numeric vector giving one dimension of
    spatial coordinates; optional, default is \code{NULL}.}
  \item{coordz}{A numeric vector giving one dimension of
    spatial coordinates; optional, default is \code{NULL}.}
  \item{coordt}{A numeric vector giving one dimension of temporal coordinates;
    optional, default is \code{NULL} (in which case a spatial random field is assumed).}
  \item{coordx_dyn}{A list of \eqn{m} numeric \eqn{d_t \times 2} matrices
    containing dynamic (in time) spatial coordinates; optional, default is \code{NULL}.}
  \item{corrmodel}{Numeric; the ID of the correlation model.}
  \item{data}{A numeric vector, or a \eqn{n \times d} matrix, or a \eqn{d \times d \times n} array of observations.}
  \item{distance}{String; the name of the spatial distance.
    Default is \code{"Eucl"} (Euclidean distance). See Details.}
  \item{flagcorr}{Numeric vector of binary values indicating which
    parameters of the correlation function will be estimated.}
  \item{flagnuis}{Numeric vector of binary values indicating which
    nuisance parameters will be estimated.}
  \item{fixed}{Numeric vector of parameters considered as known values.}
  \item{GPU}{Numeric; if \code{NULL} (default), no GPU computation is performed.}
  \item{grid}{Logical; if \code{FALSE} (default), data are interpreted as vector or \eqn{n \times d} matrix;
    if \code{TRUE}, then a \eqn{d \times d \times n} array is considered.}
  \item{likelihood}{String; configuration of the composite likelihood (see \code{\link{GeoFit}}).}
  \item{local}{Numeric; number of local work-items of the GPU.}
  \item{lower}{Named list; optional lower bounds for parameters when using
    optimizers \code{L-BFGS-B}, \code{nlminb}, or \code{optimize}.
    Names must match those in the \code{start} list.}
  \item{model}{Numeric; ID of the density associated with the likelihood objects.}
  \item{n}{Numeric; number of trials in binomial random fields.}
  \item{namescorr}{Character vector; names of the correlation parameters.}
  \item{namesnuis}{Character vector; names of the nuisance parameters.}
  \item{namesparam}{Character vector; names of the parameters to be maximized.}
  \item{numparam}{Numeric; number of parameters to be maximized.}
  \item{numparamcorr}{Numeric; number of correlation parameters.}
  \item{optimizer}{String; optimization algorithm (see \code{\link{optim}}).
    Default is \code{"Nelder-Mead"}. Other options: \code{"nlm"}, \code{"BFGS"}, \code{"L-BFGS-B"}, \code{"nlminb"}.
    For \code{"L-BFGS-B"} and \code{"nlminb"} bounds can be provided.
    For 1D optimization, \code{optimize} is used.}
  \item{onlyvar}{Logical; if \code{TRUE} (and \code{varest} is \code{TRUE}),
    only the variance-covariance matrix is computed without optimizing.
    Default is \code{FALSE}.}
  \item{parallel}{Logical; if \code{TRUE}, optimization uses \code{optimParallel}
    with maximum cores, when \code{optimizer} is \code{"L-BFGS-B"}.
    Default is \code{FALSE}.}
  \item{param}{Numeric vector of parameter values.}
  \item{spacetime}{Logical; if \code{TRUE}, the random field is spatio-temporal,
    otherwise spatial.}
  \item{type}{String; type of likelihood object.
    Default is \code{"Pairwise"} (marginal composite likelihood formed by pairwise marginal likelihoods).}
  \item{upper}{Named list; optional upper bounds for parameters when using
    optimizers \code{L-BFGS-B}, \code{nlminb}, or \code{optimize}.
    Names must match those in the \code{start} list.}
  \item{varest}{Logical; if \code{TRUE}, variance estimates and standard errors are returned.
    Default is \code{FALSE}.}
  \item{weigthed}{Logical; if \code{TRUE}, decreasing weights from a compactly supported
    correlation function with compact support \code{maxdist} (\code{maxtime}) are used.}
  \item{ns}{Numeric; number of (dynamic) temporal instants.}
  \item{X}{Numeric; matrix of space-time covariates in the linear mean specification.}
  \item{sensitivity}{Logical; if \code{TRUE}, the sensitivity matrix is computed.}
  \item{MM}{Numeric; a non-constant fixed mean.}
  \item{aniso}{Logical; whether anisotropy should be considered.}
}
\value{
  Returns a list from an \code{optim} call.
}
\details{
  Subroutine called by \code{GeoFit}.
  The procedure estimates model parameters by maximization of the composite log-likelihood.
}
\seealso{
  \code{\link{GeoFit}}
}
\author{
  Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com}, \url{https://sites.google.com/view/moreno-bevilacqua/home}, \cr
  Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}, \cr
  Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl}, \url{https://www.researchgate.net/profile/Christian-Caamano}
}
\keyword{Composite}
