\name{GeoModels-3D}
\alias{GeoModels-3D}
\title{Three-dimensional spatial coordinates in GeoModels}
\description{
Conventions and current limitations for purely spatial models whose locations
have three Cartesian coordinates.
}
\details{
For an irregular three-dimensional spatial data set, supply either an
\eqn{N \times 3} numeric matrix in \code{coordx}, or three equal-length
vectors through \code{coordx}, \code{coordy}, and \code{coordz}. Use
\code{grid = FALSE} and \code{distance = "Eucl"}.

Euclidean lags use all three coordinates. The principal fitting, covariance,
variogram, Cholesky simulation, neighborhood, and global/local kriging paths
support this representation.

Regular three-dimensional grids specified by three coordinate axes and
\code{grid = TRUE} are not implemented. Supply their Cartesian product as an
explicit \eqn{N \times 3} matrix and use \code{grid = FALSE}.

For purely spatial univariate fields, spectral turning-bands simulation accepts
explicit irregular three-dimensional Euclidean coordinates through
\code{GeoSimapprox(..., method = "TB", grid = FALSE)}. The same TB simulator
can be used by \code{GeoVarest()} and \code{GeoVarestbootstrap()} for such
fits. Three-dimensional TB does not support bivariate or spatio-temporal
models, and does not accept \code{anisopars}; transform coordinates explicitly
when needed. CE remains restricted to regular two-dimensional grids. Use
\code{method = "cholesky"} for unsupported three-dimensional cases.

Distances \code{"Geod"} and \code{"Chor"} require exactly two columns
(longitude and latitude) and cannot be combined with a third coordinate.

For \code{GeoVariogramDir()}, three-dimensional lag lengths use all three
coordinates, while direction classes are horizontal azimuths computed from
the \eqn{x}- and \eqn{y}-coordinate differences. This is not a full
three-dimensional directional analysis with azimuth and elevation.
}
\seealso{
\code{\link{GeoFit}}, \code{\link{GeoCovmatrix}},
\code{\link{GeoVariogram}}, \code{\link{GeoVariogramDir}},
\code{\link{GeoSim}}, \code{\link{GeoSimapprox}},
\code{\link{GeoKrig}}, \code{\link{GeoKrigloc}}
}
