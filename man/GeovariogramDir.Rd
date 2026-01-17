\name{GeoVariogramDir}
\alias{GeoVariogramDir}
\encoding{UTF-8}
\title{Empirical directional semivariogram}

\description{
  Computes empirical semivariograms in multiple directions (e.g., 0, 45, 90, 135 degrees) to assess
  spatial anisotropy.
}

\usage{
GeoVariogramDir(data, coordx, coordy = NULL, coordz = NULL,
                directions = c(0, 45, 90, 135), tolerance = 22.5, numbins = 13,
                maxdist = NULL, neighb = NULL, distance = "Eucl",
                subsample = 1)
}

\arguments{
  \item{data}{A numeric vector containing the observed values at each location.}
  \item{coordx}{
    Spatial coordinates. Either a numeric vector giving the first coordinate, or a matrix with 2 (or 3) columns.
    If a matrix is provided, \code{coordy} and \code{coordz} are ignored.
  }
  \item{coordy}{
    A numeric vector of the second coordinate. Optional; defaults to \code{NULL}.
  }
  \item{coordz}{
    A numeric vector of the third coordinate (if needed). Optional; defaults to \code{NULL}.
  }
  \item{directions}{
    Numeric vector giving the principal directions (in degrees) for which the semivariogram is computed
    (default: \code{c(0, 45, 90, 135)}).
  }
  \item{tolerance}{
    Angular tolerance (in degrees) around each direction (default: 22.5).
  }
  \item{numbins}{
    Number of distance bins for the empirical semivariogram (default: 13).
  }
  \item{maxdist}{
    Maximum spatial distance to consider between pairs. If \code{NULL}, an internal default is used.
  }
  \item{neighb}{
    Number of nearest neighbors to use for each location. If \code{NULL}, an internal default is used.
  }
  \item{distance}{
    Type of distance metric to use (default: \code{"Eucl"}). See \code{\link{GeoFit}} for options.
  }
  \item{subsample}{
    Numeric in \eqn{(0,1]}{(0,1]}. Proportion of spatial locations used to compute the directional semivariograms
    (useful for large datasets). Default is 1 (use all locations).
  }
}

\details{
  The function computes empirical semivariograms for several directions by:
  \itemize{
    \item Selecting pairs of points within \code{maxdist} and among the \code{neighb} nearest neighbors using
          \code{\link{GeoNeighIndex}}.
    \item Computing squared differences for each selected pair.
    \item Assigning each pair to a directional class if the vector connecting the pair falls within the
          specified angular tolerance around a given direction.
    \item Binning pairs by distance and computing the average squared difference (semivariogram) within each bin.
  }

  The direction is defined in the \emph{xy-plane} even in 3D. For 2D data, set \code{coordz = NULL}.

  The \code{subsample} argument can be used to reduce computational cost by randomly selecting a subset of
  spatial locations prior to pair construction. This is particularly useful for large datasets.
}

\value{
  A list of class \code{"GeoVariogramDir"} with one element for each direction. Each element is a list with components:
  \item{centers}{Centers of the distance bins.}
  \item{gamma}{Empirical semivariogram values for each bin.}
  \item{npairs}{Number of point pairs in each bin.}
}

\seealso{
  \code{\link{GeoVariogram}}, \code{\link{GeoNeighIndex}}
}

\examples{
require(GeoModels)
set.seed(960)
NN <- 2500
coords <- cbind(runif(NN), runif(NN))
scale <- 0.5/3
param <- list(mean = 0, sill = 1, nugget = 0, scale = scale, smooth = 0.5)
corrmodel <- "Matern"

set.seed(951)
data <- GeoSim(coordx = coords, corrmodel = corrmodel,
               model = "Gaussian", param = param)$data

vario_dir <- GeoVariogramDir(data = data, coordx = coords, maxdist = 0.4,
                             subsample = 0.5)

plot(vario_dir, ylim = c(0,1))
}

\keyword{Variogram}
