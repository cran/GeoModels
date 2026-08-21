\name{GeoVariogram}
\alias{GeoVariogram}
\encoding{UTF-8}
\title{Empirical semivariogram estimation}
\description{
 Computes an empirical estimate of the semivariogram for spatial, spatio-temporal,
 and bivariate random fields.
}
\usage{
GeoVariogram(data, coordx, coordy=NULL, coordz=NULL, coordt=NULL,
 coordx_dyn=NULL, cloud=FALSE, distance="Eucl",
 grid=FALSE, maxdist=NULL, neighb=NULL,
 maxtime=NULL, numbins=NULL,
 radius=1, type='variogram', bivariate=FALSE,
 subsample=1, subsample_t=1, numbins_t=NULL, directed=FALSE)
}
\arguments{
 \item{data}{
 A numeric vector of length \eqn{d}{d} (a single spatial realisation), or a
 \eqn{d \times d}{d x d} matrix (a single realisation on a regular grid), or a
 \eqn{t \times d}{t x d} matrix (a single spatio-temporal realisation), or a
 \eqn{d \times d \times t}{d x d x t} array (a single spatio-temporal realisation on a regular grid).
 See \code{\link{GeoFit}} for details.
 }
 \item{coordx}{
 Spatial coordinates. Either a numeric vector giving the first coordinate, or a
 \eqn{d \times 2}{d x 2} (or \eqn{d \times 3}{d x 3}) matrix of coordinates.
 If \code{distance} refers to great-circle distances, coordinates must be provided in
 lon/lat format (decimal degrees) and the sphere radius is set by \code{radius}.
 }
 \item{coordy}{
 A numeric vector giving the second spatial coordinate.
 Optional, default is \code{NULL}.
 }
 \item{coordz}{
 A numeric vector giving the third spatial coordinate (if needed).
 Optional, default is \code{NULL}.
 }
 \item{coordt}{A numeric vector of temporal coordinates. If \code{NULL} (default), a purely spatial random field is assumed. Temporal coordinates may be irregularly spaced; temporal lags are computed from the supplied coordinate values.}
 \item{coordx_dyn}{
 For dynamic locations, a list with one two- or three-column coordinate matrix per temporal instant.  The list length must match \code{length(coordt)}; \code{data[[t]]} must have one value per row of \code{coordx_dyn[[t]]}, in the same order.
 }
 \item{cloud}{
 Logical; if \code{TRUE} the semivariogram cloud is computed for a univariate
 purely spatial field. If \code{FALSE} (default), a binned empirical semivariogram
 is returned. For large data sets, combine \code{cloud=TRUE} with \code{maxdist},
 \code{neighb}, or \code{subsample} to limit the number of returned pairs.
 }
 \item{distance}{
 String specifying the spatial distance. Default is \code{"Eucl"} (Euclidean distance).
 See \code{\link{GeoFit}} for details.
 }
 \item{grid}{
 Logical; if \code{FALSE} (default) data are interpreted as observations on irregularly
 spaced locations. If \code{TRUE}, data are interpreted as observations on a regular grid.
 }
 \item{maxdist}{
 A positive finite numeric maximum spatial distance. In the bivariate case a
 scalar or a length-3 vector (first marginal, cross, second marginal) is
 accepted. Use \code{NULL} to let the function use the available spatial
 range. See Details.
 }
 \item{neighb}{
 Numeric; an optional positive integer indicating the order of neighborhood
 (useful for large datasets). Neighborhood pair selection is available for
 spatial and spatio-temporal semivariograms. In the bivariate case a length-3
 vector can be used for the first marginal, cross, and second marginal pair sets.
 See Details.
 }
 \item{maxtime}{A positive finite maximum temporal lag, expressed in the same units as \code{coordt}, to be considered for spatio-temporal semivariograms. See Details.}
 \item{numbins}{
 Historical GeoModels argument controlling the spatial bin grid. It is the
 number of spatial bin boundaries, so the number of empirical spatial classes
 is \code{numbins - 1}. The default is 13 boundaries (12 classes).
 }
 \item{numbins_t}{
 Optional positive integer giving the number of temporal classes for irregularly
 spaced \code{coordt}. If \code{NULL}, \code{numbins - 1} classes are used.
 For regularly spaced times the attainable temporal lags are generated directly
 from the common spacing and this argument is not used.
 }
 \item{directed}{
 Logical used only when \code{neighb} is supplied. The default \code{FALSE}
 removes reciprocal nearest-neighbour duplicates so that the empirical
 semivariogram is based on unordered pairs. Set \code{TRUE} to retain the
 directed nearest-neighbour graph used by pairwise composite-likelihood code.
 }
 \item{radius}{
 Numeric; radius of the sphere when using great-circle distances. Default is 1.
 }
 \item{type}{
 String; type of semivariogram. Currently available: \code{"variogram"}.
 }
 \item{bivariate}{
 Logical; if \code{FALSE} (default) data are interpreted as univariate spatial/spatio-temporal
 realisations. If \code{TRUE}, \code{data} is interpreted as a realisation from a bivariate field
 and (cross-)semivariograms are computed.
 }
 \item{subsample}{
 Numeric in \eqn{(0,1]}{(0,1]}. Proportion of spatial locations to be used to compute the semivariogram
 (useful for large datasets). Default is 1 (use all locations).
 }
 \item{subsample_t}{
 Numeric in \eqn{(0,1]}{(0,1]}. Proportion of time points to be used in spatio-temporal settings
 (when \code{coordt} is provided). Default is 1 (use all time points).
 }
}

\details{
 We report the definition of the semivariogram in the spatial case; extensions to spatio-temporal
 and bivariate settings are based on the same principles.

 For a spatial random field \eqn{Z(\cdot)}{Z(.)}, the (classical) binned semivariogram estimator is
 defined as
 \deqn{\hat{\gamma}(h) = \frac{1}{2 |N(h)|}\sum_{(x_i,x_j)\in N(h)} \{Z(x_i)-Z(x_j)\}^2,}
 where \eqn{N(h)}{N(h)} is the set of all sample pairs whose spatial distance falls within a tolerance
 region around lag \eqn{h}{h} (equally spaced intervals are used when \code{cloud=FALSE}).

 The historical \code{numbins} argument sets the number of spatial bin
 boundaries; hence \code{numbins - 1} empirical spatial classes are formed when
 \code{cloud=FALSE}.

 The \code{maxdist} argument sets a strictly positive finite maximum spatial
 distance. If no pair falls below the requested cutoff, a valid empty
 semivariogram object is returned instead of constructing decreasing bins.

 The \code{maxdist} option can be combined with \code{neighb} to reduce the number of pairs when handling
 large datasets, by restricting computations to local neighborhoods. By default
 reciprocal directed nearest-neighbour edges are deduplicated before binning;
 \code{directed=TRUE} retains the directed graph. For chordal and geodesic
 distances, neighbour ordering is obtained from three-dimensional unit-sphere
 coordinates and the reported lags are then evaluated in the requested metric.

 Spatial and temporal bins are left-closed and right-open, except for the last
 bin which is closed on the right. Thus a pair exactly at the maximum retained
 lag is not discarded.

 The \code{maxtime} argument sets the maximum temporal lag considered for
 spatio-temporal semivariograms. For regularly spaced times, attainable temporal
 lags are generated in linear memory/time from the common spacing; no
 \eqn{T \times T}{T x T} matrix of all time differences is formed. For
 irregular times a controlled grid of \code{numbins_t} temporal classes is used.
 The returned spatio-temporal surface is rectangular on \code{centers} by
 \code{centert}; cells with no valid pairs are returned as \code{NA} with count
 zero rather than being removed.

 For dynamic sites, the temporal marginal is the empirical \eqn{\gamma(0,u)}
 and therefore uses only spatially collocated locations (up to numerical
 tolerance). Nearby but non-collocated locations contribute to the positive-
 distance space-time surface, not to the temporal marginal. Consequently, in a
 fully dynamic design with no spatial locations repeated across temporal
 instants, \code{variogramt} is \code{NA} in the affected temporal classes.
 This is an expected property of the sampling design, not a failure of the
 space-time variogram: the positive-distance surface \eqn{\gamma(h,u)} remains
 empirically estimable. In this case \code{plot.GeoVariogram} labels the panel
 ``Temporal marginal unavailable (no repeated spatial locations)''. Repeated
 spatial locations across times are needed only when the empirical temporal
 marginal itself is required.

 In the bivariate case the two marginal semivariograms and the cross-semivariogram
 are evaluated on the same spatial bins. The cross-semivariogram uses the classical
 increment-product estimator on unordered positive-distance pairs. The trivial
 collocated contribution at lag zero is not mixed into the first positive spatial
 class. If two dynamic/support coordinate sets are supplied, they must be aligned
 to define this estimator unambiguously.

 The \code{subsample} and \code{subsample_t} arguments provide additional control
 for large datasets by using only a proportion of spatial locations and/or time
 points. With \code{grid=TRUE}, the grid is first normalized to explicit spatial
 locations, so full-grid and subsampled calculations follow the same path.

 Missing values \code{NA}/\code{NaN} are allowed and pairs involving them are
 skipped consistently. Infinite observations are rejected.
}

\section{Spatio-temporal ordering}{
For fixed sites, \code{data} is a \eqn{T \times N}{T x N} matrix with times in
rows and sites in columns, and is internally read in the order
\code{c(t(data))}.  For dynamic sites, \code{data} and \code{coordx_dyn} are
aligned lists; the function concatenates complete temporal blocks in list
order.  See \code{\link{GeoModels-spacetime-ordering}}.
}

\value{
 Returns an object of class \code{"GeoVariogram"}. The list contains, as applicable:

 \item{bins}{Spatial bin boundaries when \code{cloud=FALSE}; spatial pair distances when \code{cloud=TRUE}.}
 \item{bint}{Temporal lag representatives for a spatio-temporal variogram.}
 \item{bivariate}{Logical indicating a bivariate empirical variogram.}
 \item{cloud}{Logical indicating a variogram cloud.}
 \item{centers}{Spatial bin centers.}
 \item{centert}{Temporal lag representatives used by the rectangular spatio-temporal surface.}
 \item{distance}{Spatial distance type used to construct the empirical variogram.}
 \item{radius}{Sphere radius used for chordal or geodesic distances.}
 \item{grid}{Logical recording whether the original input was supplied as a grid.}
 \item{neighb}{Neighborhood order used for pair selection, or \code{NULL}.}
 \item{directed}{Whether reciprocal nearest-neighbour pairs were retained.}
 \item{lenbins}{Numbers of pairs in the spatial bins. In the bivariate case this is a two-row matrix.}
 \item{lenbinst}{Numbers of pairs in the cross/spatio-temporal bins. For space-time objects this follows the same row-major ordering as \code{variogramst}.}
 \item{lenbint}{Numbers of pairs in the temporal bins.}
 \item{maxdist}{Maximum spatial distance requested by the user.}
 \item{maxtime}{Maximum temporal lag requested by the user.}
 \item{numbins_t}{Requested number of temporal classes for irregular times, or \code{NULL}.}
 \item{regular}{Logical indicating regularly spaced temporal coordinates for a space-time object.}
 \item{time.breaks}{Internal temporal bin boundaries used for space-time binning.}
 \item{spacetime_dyn}{Logical indicating dynamic spatial coordinates.}
 \item{temporal.margin}{For space-time objects, a label indicating that the temporal margin is based on same-site/collocated pairs.}
 \item{subsample}{Spatial subsampling proportion.}
 \item{subsample_t}{Temporal subsampling proportion.}
 \item{variograms}{Empirical spatial semivariogram; a two-row matrix in the bivariate case.}
 \item{variogramst}{Empirical cross-semivariogram in the bivariate case, or the rectangular spatio-temporal surface stored in row-major order. Empty cells are \code{NA}.}
 \item{variogramt}{Empirical temporal marginal semivariogram.}
 \item{type}{Type of empirical variogram.}
}

\references{
 Cressie, N. A. C. (1993) \emph{Statistics for Spatial Data}. New York: Wiley.

 Gaetan, C. and Guyon, X. (2010) \emph{Spatial Statistics and Modeling}. Springer-Verlag, New York.
}

\seealso{
\code{\link{GeoFit}} for model fitting,
  \code{\link{GeoWLS}} for weighted least-squares estimation from empirical variograms,
  \code{\link{GeoCovariogram}} for fitted covariance and variogram values,
  \code{\link{plot.GeoVariogram}} for plotting empirical variograms.
}

\author{
 Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com}, \url{https://sites.google.com/view/moreno-bevilacqua/home},
 Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
 Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl}, \url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{
library(GeoModels)

################################################################
### Example 1. Empirical semivariogram from a spatial Gaussian
### random field with Matérn correlation.
################################################################
set.seed(514)
x = runif(200, 0, 1)
y = runif(200, 0, 1)
coords = cbind(x,y)

corrmodel = "Matern"
mean = 0
sill = 1
nugget = 0
scale = 0.3/3
smooth = 0.5

data = GeoSim(coordx=coords, corrmodel=corrmodel,
 param=list(mean=mean, smooth=smooth, sill=sill,
 nugget=nugget, scale=scale))$data

vario = GeoVariogram(coordx=coords, data=data, maxdist=0.6)
plot(vario, pch=20, ylim=c(0,1), ylab="Semivariogram", xlab="Distance")

################################################################
### Example 2. Empirical semivariogram for a spatio-temporal
### Gaussian random field with Gneiting correlation.
################################################################
set.seed(331)
x = runif(200, 0, 1)
y = runif(200, 0, 1)
coords = cbind(x,y)
times = seq(1,10,1)

data = GeoSim(coordx=coords, coordt=times, corrmodel="gneiting",
 param=list(mean=0, scale_s=0.08, scale_t=0.4, sill=1,
 nugget=0, power_s=1, power_t=1, sep=0.5))$data

vario_st = GeoVariogram(data=data, coordx=coords, coordt=times,
 maxtime=7, maxdist=0.5)
plot(vario_st, pch=20)

################################################################
### Example 3. Empirical (cross-)semivariograms for a bivariate
### Gaussian random field with Bi-Matérn covariance.
################################################################
set.seed(293)
x = runif(400, 0, 1)
y = runif(400, 0, 1)
coords = cbind(x,y)

param = list(mean_1=0, mean_2=0,
 scale_1=0.1/3, scale_2=0.15/3, scale_12=0.15/3,
 sill_1=1, sill_2=1,
 nugget_1=0, nugget_2=0,
 smooth_1=0.5, smooth_12=0.5, smooth_2=0.5,
 pcol=0.3)

data = GeoSim(coordx=coords, corrmodel="Bi_matern", param=param)$data
biv_vario = GeoVariogram(data, coordx=coords, bivariate=TRUE, maxdist=0.5)
plot(biv_vario, pch=20)
}

\keyword{Variogram}
