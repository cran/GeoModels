\name{GeoKrigWeights}
\alias{GeoKrigWeights}
\encoding{UTF-8}
\title{
Compute validated Gaussian kriging weights for spatial and spatio-temporal random fields
}
\description{
Given a set of spatial locations (and possibly temporal instants), the function
returns Gaussian kriging weights and the covariance quantities used to build the
kriging system.  Model-specific non-Gaussian prediction is provided by
\code{\link{GeoKrig}} and \code{\link{GeoKrigloc}} rather than by this low-level
weights helper.
}
\usage{
GeoKrigWeights(
 coordx, coordy = NULL, coordz = NULL, coordt = NULL,
 coordx_dyn = NULL, corrmodel, distance = "Eucl",
 grid = FALSE, loc, method = "cholesky", model = "Gaussian",
 n = 1, nloc = NULL, param, anisopars = NULL, radius = 1,
 sparse = FALSE, time = NULL, which = 1, copula = NULL,
 X = NULL, Xloc = NULL, Mloc=NULL)
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
 \item{copula}{Must be \code{NULL}.  Copula conditional-simulation code requests
 these weights on the latent Gaussian scale.}
 \item{X}{Numeric design matrix at observation locations.  Fixed-location observations are ordered time then site; dynamic observations are ordered by the temporal blocks of \code{coordx_dyn}.}
 \item{Xloc}{Numeric design matrix at prediction tasks, ordered location then time: all requested times for the first row of \code{loc}, then all requested times for the second row, and so on.}
 \item{Mloc}{Numeric vector of known prediction means in the same location-major order as \code{Xloc}.  Use either \code{Mloc} or \code{Xloc}, not both.}
}
\details{

Mean inputs follow \code{\link{GeoKrig}}.  Coefficients \code{mean},
\code{mean1}, and so on correspond in order to the columns of \code{X}; an
intercept-only model uses only \code{mean}.  A vector \code{param$mean} is a
known observation mean and is mutually exclusive with \code{X}.  At
prediction locations, use either \code{Xloc} or \code{Mloc}, not both.

The function builds the kriging system
\deqn{\Sigma \mathbf{w} = \boldsymbol{\sigma}_0}{Sigma w = sigma_0}
where \eqn{\Sigma}{Sigma} is the covariance matrix between observed
locations and \eqn{\boldsymbol{\sigma}_0}{sigma_0} the vector of
covariances between observed and prediction locations. 
The solution \eqn{\mathbf{w}}{w} is returned together with \eqn{\Sigma}{Sigma}
and, optionally, \eqn{\Sigma^{-1}}{Sigma^-1}.
The returned weights solve the Gaussian covariance system only.  Mean inputs are
validated for ordering compatibility, but no non-Gaussian covariance transform is
performed here.  No actual prediction is carried out; for full Gaussian or
non-Gaussian prediction use \code{\link{GeoKrig}}.
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
A list containing:
 \item{weights}{Numeric matrix of Gaussian kriging weights.  The validated
 orientation is observations by prediction tasks, as recorded in
 \code{weights_orientation}.}
 \item{weights_orientation}{Character string describing the weight-matrix
 orientation.}
 \item{covmatrix}{Observed Gaussian covariance matrix.}
 \item{CC}{Observation--prediction Gaussian cross-covariance matrix.}
 \item{model}{The canonical model name (currently \code{"Gaussian"}).}
 \item{corrmodel}{Input correlation model.}
 \item{param}{Parameters used to construct the Gaussian covariance matrix.}
 \item{bivariate}{Logical indicating a bivariate Gaussian correlation model.}
 \item{spacetime}{Logical indicating a space-time correlation model.}
 \item{loc}{Prediction locations.}
 \item{tloc}{Number of requested prediction times.}
 \item{Xloc,Mloc,mean_obs}{Normalized mean inputs retained for alignment and
 conditional-simulation workflows.}
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


\examples{
library(GeoModels)
################################################################

################################################################
###
### Example 1. Spatial kriging weights for
### Gaussian random fields with Gen wendland correlation.
###
################################################################

model="Gaussian"
set.seed(79)
x = runif(300, 0, 1)
y = runif(300, 0, 1)
coords=cbind(x,y)

corrmodel = "GenWend"
mean=0; sill=5; nugget=0
scale=0.2;smooth=0;power2=4

param=list(mean=mean,sill=sill,nugget=nugget,scale=scale,smooth=smooth,power2=power2)

# Simulation of the spatial Gaussian random field:
data = GeoSim(coordx=coords, corrmodel=corrmodel,model=model,
 param=param)$data

xx=seq(0,1,0.25)
loc_to_pred=as.matrix(expand.grid(xx,xx))

W=GeoKrigWeights(,coordx=coords,loc=loc_to_pred,corrmodel=corrmodel,
 model=model,param=param)

dim(W$weights) ### kriging weights
}

\keyword{Kriging}
