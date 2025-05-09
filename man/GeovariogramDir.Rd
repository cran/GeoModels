\name{GeoVariogramDir}
\alias{GeoVariogramDir}
\encoding{UTF-8}
\title{Empirical Directional Semivariogram}

\description{
Computes the empirical semivariogram in multiple directions (e.g., 0, 45, 90, 135 degrees) to assess spatial anisotropy, using only relevant pairs of points  selected using \code{maxdist} and 
\code{neighb} through \code{\link{GeoNeighIndex}}. 
}

\usage{
GeoVariogramDir(data, coordx, coordy = NULL, coordz = NULL,
  directions = c(0, 45, 90, 135), tolerance = 22.5, numbins = 13,
  maxdist = NULL, neighb = NULL, distance = "Eucl")
}

\arguments{
  \item{data}{A numeric vector containing the observed values at each location.}
  \item{coordx}{A numeric vector or matrix of the x-coordinates of the locations. If a matrix with 2 or 3 columns is provided, \code{coordy} and \code{coordz} are ignored.}
  \item{coordy}{A numeric vector of the y-coordinates of the locations. Optional; defaults to \code{NULL} for 1D data or if \code{coordx} is a matrix/data.frame.}
  \item{coordz}{A numeric vector of the z-coordinates of the locations. Optional; defaults to \code{NULL} for 2D data or if \code{coordx} is a matrix/data.frame with 2 columns.}
  \item{directions}{A numeric vector giving the principal directions (in degrees) for which to compute the semivariogram (default: \code{c(0, 45, 90, 135)}).}
  \item{tolerance}{Angular tolerance (in degrees) for each direction (default: 22.5).}
  \item{numbins}{Number of distance bins for the empirical semivariogram (default: 13).}
  \item{maxdist}{Maximum distance to consider between pairs (default: \code{Inf}).}
  \item{neighb}{Number of nearest neighbors to use for each location (default: \code{length(data) - 1}).}
  \item{distance}{Type of distance metric to use (default: \code{"Eucl"}).}
}

\details{
The function computes the empirical semivariogram for several directions by:
\itemize{
  \item Selecting pairs of points within \code{maxdist} and among the \code{neighb} nearest neighbors using \code{\link{GeoNeighIndex}}.
  \item Calculating the squared differences for each pair.
  \item Assigning each pair to a directional bin if the vector connecting the pair falls within the specified angular tolerance of a given direction.
  \item Binning the pairs by distance and computing the average squared difference (semivariogram) for each bin.
}
The direction is defined in the \emph{xy-plane} even in 3D. For 2D data, set \code{coordz = NULL}.

This implementation is optimized: distance bins and directional masks are precomputed for all pairs, minimizing repeated computations for each direction.
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
NN <- 1500
coords <- cbind(runif(NN), runif(NN))
scale <- 0.5/3
param <- list(mean = 0, sill = 1, nugget = 0, scale = scale, smooth = 0.5)
corrmodel <- "Matern"

set.seed(951)
data <- GeoSim(coordx = coords, corrmodel = corrmodel,
               model = "Gaussian", param = param)$data

vario_dir <- GeoVariogramDir(data = data, coordx = coords, maxdist = 0.4)

plot(vario_dir,ylim=c(0,1))
}

\keyword{Variogram}
