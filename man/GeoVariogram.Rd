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
             subsample=1, subsample_t=1)
}
\arguments{
  \item{data}{
    A numeric vector of length \eqn{d}{d} (a single spatial realisation),  or a
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
  \item{coordt}{
    A numeric vector of temporal coordinates. If \code{NULL} (default), a purely spatial
    random field is assumed.
  }
  \item{coordx_dyn}{
    A list of \eqn{m}{m} numeric matrices \eqn{d_t \times 2}{d_t x 2} providing time-varying
    spatial coordinates (dynamic locations). Optional, default is \code{NULL}.
  }
  \item{cloud}{
    Logical; if \code{TRUE} the semivariogram cloud is computed. If \code{FALSE} (default),
    a binned empirical semivariogram is returned.
  }
  \item{distance}{
    String specifying the spatial distance. Default is \code{"Eucl"} (Euclidean distance).
    See the \bold{Details} section of \code{\link{GeoFit}}.
  }
  \item{grid}{
    Logical; if \code{FALSE} (default) data are interpreted as observations on irregularly
    spaced locations. If \code{TRUE}, data are interpreted as observations on a regular grid.
  }
  \item{maxdist}{
    Numeric; maximum spatial distance to be considered in semivariogram estimation.
    See \bold{Details}.
  }
  \item{neighb}{
    Numeric; an optional positive integer indicating the order of neighborhood
    (useful for large datasets). See \bold{Details}.
  }
  \item{maxtime}{
    Numeric; maximum temporal lag to be considered for spatio-temporal semivariograms.
    See \bold{Details}.
  }
  \item{numbins}{
    Numeric; number of distance bins used to compute the binned semivariogram.
    See \bold{Details}.
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

  The \code{numbins} argument sets the number of spatial lag bins used when \code{cloud=FALSE}.

  The \code{maxdist} argument sets the maximum spatial distance considered in the estimation.

  The \code{maxdist} option can be combined with \code{neighb} to reduce the number of pairs when handling
  large datasets, by restricting computations to local neighborhoods.

  The \code{maxtime} argument sets the maximum temporal lag considered for spatio-temporal semivariograms.

  The \code{subsample} and \code{subsample_t} arguments provide additional control for large datasets by
  using only a proportion of spatial locations and/or time points.
}

\value{
  Returns an object of class \code{Variogram}.
  An object of class \code{Variogram} is a list containing (at most) the following components:

  \item{bins}{Spatial distance bins if \code{cloud=FALSE}. If \code{cloud=TRUE}, all spatial pairwise distances.}
  \item{bint}{Temporal distance bins if \code{cloud=FALSE}. If \code{cloud=TRUE}, all temporal pairwise distances.}
  \item{cloud}{Logical; \code{TRUE} if the variogram cloud is returned, \code{FALSE} otherwise.}
  \item{centers}{Centers of the spatial bins.}
  \item{distance}{Type of spatial distance.}
  \item{lenbins}{Number of pairs in each spatial bin.}
  \item{lenbinst}{Number of pairs in each spatio-temporal bin.}
  \item{lenbint}{Number of pairs in each temporal bin.}
  \item{maxdist}{Maximum spatial distance used in the estimation; \code{NULL} if not specified.}
  \item{maxtime}{Maximum temporal lag used in the estimation; \code{NULL} if not specified.}
  \item{spacetime_dyn}{Logical; \code{TRUE} if dynamic coordinates (\code{coordx_dyn}) are used.}
  \item{variograms}{Empirical spatial semivariogram.}
  \item{variogramst}{Empirical spatio-temporal semivariogram.}
  \item{variogramt}{Empirical temporal semivariogram.}
  \item{type}{Type of estimated semivariogram.}
}

\references{
  Cressie, N. A. C. (1993) \emph{Statistics for Spatial Data}. New York: Wiley.

  Gaetan, C. and Guyon, X. (2010) \emph{Spatial Statistics and Modeling}. Springer-Verlag, New York.
}

\seealso{\code{\link{GeoFit}}}

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
