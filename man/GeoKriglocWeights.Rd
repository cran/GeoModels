\name{GeoKriglocWeights}
\alias{GeoKriglocWeights}
\encoding{UTF-8}
\title{
Compute validated local Gaussian kriging weights
}
\description{
Given a set of spatial locations (and possibly temporal instants), the function
computes local Gaussian kriging weights for each requested prediction task using
a neighborhood selected from the observations.  Non-Gaussian local prediction is
implemented by \code{\link{GeoKrigloc}}, not by this low-level weights helper.
}
\usage{
GeoKriglocWeights(coordx=NULL, coordy=NULL, coordz=NULL, coordt=NULL, coordx_dyn=NULL,
 corrmodel, distance="Eucl", grid=FALSE, loc, neighb=NULL,
 maxdist=NULL, maxtime=NULL, method="cholesky", model="Gaussian",
 n=1, nloc=NULL, param, anisopars=NULL,
 radius=1, sparse=FALSE, time=NULL, which=1,
 copula=NULL, X=NULL, Xloc=NULL, Mloc=NULL,parallel=FALSE,
 ncores=6, compact=FALSE)
}
\arguments{
 \item{coordx}{Numeric (\eqn{d\times 2}{d x 2}) or (\eqn{d\times 3}{d x 3}) matrix
 of spatial coordinates. Coordinates on a sphere are accepted
 (lon/lat in decimal degrees) when \code{distance = "Sphere"}.}
 \item{coordy}{Optional numeric vector giving an additional spatial coordinate
 dimension. Ignored if \code{coordx} is already a matrix.}
 \item{coordz}{Optional numeric vector giving a third spatial coordinate
 dimension.}
 \item{coordt}{Optional numeric vector containing the temporal coordinates of the observations. If missing, a purely spatial random field is assumed. Temporal coordinates may be irregularly spaced; temporal lags are computed from the supplied coordinate values.}
 \item{coordx_dyn}{List of one two- or three-column observation-coordinate matrix per temporal instant.  Blocks are concatenated by time; rows within each block retain their matrix order.  See \code{\link{GeoModels-spacetime-ordering}}.}
 \item{corrmodel}{Character string naming a valid correlation model.
 See \code{\link{GeoCovmatrix}} for admissible choices.}
 \item{distance}{Character string specifying the spatial distance.
 Default is \code{"Eucl"} (Euclidean). See \code{\link{GeoFit}}.}
 \item{grid}{Logical. If \code{TRUE}, coordinates are interpreted as
 defining a regular grid (see \code{\link{GeoKrig}}).}
 \item{loc}{Numeric (\eqn{n\times 2}{n x 2}) matrix of locations for which
 the kriging weights are required.}
 \item{neighb}{Numeric; an optional positive integer indicating the 
 order of the neighborhood.}
 \item{maxdist}{Numeric; an optional positive value indicating the distance in the spatial neighborhood.}
 \item{maxtime}{Numeric; an optional non-negative maximum temporal-distance threshold, expressed in the same units as \code{coordt}, used for the local temporal neighborhood.}
 \item{method}{Character string indicating the matrix factorisation
 used to solve the kriging system: \code{"cholesky"} (default) or \code{"svd"}.}
 \item{model}{Character string.  The validated public implementation currently
 requires \code{"Gaussian"} (the alias \code{"Gauss"} is accepted).}
 \item{n}{Integer. Number of trials for Binomial random fields
 (default \code{1}).}
 \item{nloc}{Integer. Number of trials for the prediction locations
 in Binomial random fields (default \code{1}).}
 \item{param}{Named list of covariance and mean parameters.
 See \code{\link{CorrParam}} and \code{\link{GeoCovmatrix}}.}
 \item{anisopars}{List with components \code{angle} and \code{ratio}
 defining geometric anisotropy (optional).}
 \item{radius}{Positive numeric value: sphere radius when
 coordinates are lon/lat (default \code{1}).}
 \item{sparse}{Logical. If \code{TRUE}, sparse‐matrix algorithms
 (package \pkg{spam}) are employed. Only effective with compactly
 supported covariance functions.}
 \item{time}{Numeric vector giving the temporal instants for which weights are required. Values need not be equally spaced and are interpreted on the same numeric time scale as \code{coordt}. Ignored if \code{coordt} is missing.}
 \item{which}{Integer (\code{1} or \code{2}) selecting the variable whose
 weights are returned in the bivariate case.}
 \item{copula}{Must be \code{NULL}.}
 \item{X}{Numeric design matrix at observation locations.  Fixed-location observations are ordered time then site; dynamic observations are ordered by the temporal blocks of \code{coordx_dyn}.}
 \item{Xloc}{Numeric design matrix at prediction tasks, ordered location then time: all requested times for the first row of \code{loc}, then all requested times for the second row, and so on.}
 \item{Mloc}{Numeric vector of known prediction means in the same location-major order as \code{Xloc}.  Use either \code{Mloc} or \code{Xloc}, not both.}
 \item{parallel}{Logical; default \code{FALSE}. Set \code{TRUE} to compute independent local kriging systems in parallel where supported.}
 \item{ncores}{Positive integer or \code{NULL}; default \code{6}. With \code{parallel=TRUE}, an explicit integer requests that many workers, capped by detected cores and the number of available jobs. Set \code{ncores=NULL} for automatic selection, capped at six workers and normally leaving one detected core free.}
 \item{compact}{Logical. If \code{FALSE} (default), each local element contains
 the full \code{GeoKrigWeights} result for backward compatibility. If
 \code{TRUE}, each non-empty element keeps only the numeric weights and
 \code{neighbor_indices}; this memory-saving representation is used internally
 by \code{GeoSimcond(local=TRUE)}.}
}
\details{
The optional pair cache is bounded by
\code{getOption("GeoModels.local_pair_cache_max_bytes", 256 * 1024^2)}.
This is a per-process budget in bytes for native cache workspace and CSR output,
not a limit on total R memory. The peak during hash growth is included.
A zero budget disables cache construction. If the budget or integer-index
limit would be exceeded, the uncached calculation is used automatically.
Each parallel worker may hold its own R objects.


Mean inputs follow \code{\link{GeoKrigloc}}.  Coefficients \code{mean},
\code{mean1}, and so on correspond in order to the columns of \code{X}; an
intercept-only model uses only \code{mean}.  A vector \code{param$mean} is a
known observation mean and is mutually exclusive with \code{X}.  At local
prediction tasks, use either \code{Xloc} or \code{Mloc}, not both.

For every prediction task the function builds and solves the local kriging
system
\deqn{\Sigma \mathbf{w} = \boldsymbol{\sigma}_0}{Sigma w = sigma_0},
where \eqn{\Sigma}{Sigma} is the Gaussian covariance matrix among the selected
observations and \eqn{\boldsymbol{\sigma}_0}{sigma_0} contains their
covariances with the prediction location.  No actual prediction is carried out;
for model-specific local prediction use \code{\link{GeoKrigloc}}.  Local results are returned directly rather
than through temporary-file references, so their components remain available
after the parallel workers have terminated.
}
\section{Spatio-temporal ordering}{
Observation-level rows follow time-major order.  For fixed sites this is the
order of \code{c(t(data))}; for dynamic sites it is the row-binding of
\code{coordx_dyn} by list element.  Prediction tasks follow location-major
order, so rows of \code{Xloc} and elements of \code{Mloc} enumerate all
requested times for the first prediction location before moving to the next
location.  See \code{\link{GeoModels-spacetime-ordering}}.
}

\value{
An object of class \code{GeoKriglocWeights}.  Its component \code{weights} is a
list ordered by prediction task.  With \code{compact=FALSE}, each
non-\code{NULL} element is the local result returned by
\code{\link{GeoKrigWeights}}, augmented with \code{neighbor_indices}. With
\code{compact=TRUE}, each non-\code{NULL} element contains only
\code{weights} and \code{neighbor_indices}.  An element is \code{NULL} when the requested neighborhood contains no
usable observations.  For purely spatial and bivariate models there is one
list element per row of \code{loc}; for space-time models the order is location
then requested time.  The remaining components record the coordinates,
model, parameters, neighborhood settings, and function call used to construct
the local systems.
}
\references{
 Gaetan, C. and Guyon, X. (2010) \emph{Spatial Statistics and Modeling}.
 Springer-Verlag, New York.
}
\seealso{
\code{\link{GeoKrigloc}} for local kriging prediction,
  \code{\link{GeoKrigWeights}} for global kriging weights,
  \code{\link{GeoCovmatrix}} for covariance matrix construction.
}
\author{
 Moreno Bevilacqua, \email{moreno.bevilacqua@uai.cl}, 
 Víctor Morales-Oñate, \email{victor.morales@uv.cl}, 
 Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl}
}
\keyword{Kriging}
