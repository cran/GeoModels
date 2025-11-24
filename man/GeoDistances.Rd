\name{GeoDistances}
\alias{GeoDistances}
\title{Compute Distance Matrices (Euclidean, Chordal, or Geodesic)}
\description{
Compute a full distance matrix between coordinates in 2D or 3D, using Euclidean, 
Chordal, or Geodesic distance. The computation is done in C for efficiency, 
with parameters passed from R via \code{dotCall64}.
}
\usage{
GeoDistances(coordx = NULL, coordy = NULL, coordz = NULL,
             distance = c("Eucl", "Chor", "Geod"), radius = 1)
}
\arguments{
  \item{coordx}{A numeric matrix of coordinates. Must have 2 or 3 columns 
                (for 2D or 3D) (if coordy and coordz are not provided). Otherwise
                 coordinates can be given via 
                \code{coordx}, \code{coordy}, and optionally \code{coordz}.}
  \item{coordy}{Optional numeric vector of  coordinates.}
  \item{coordz}{Optional numeric vector of  coordinates (for 3D Euclidean).}
  \item{distance}{Type of distance to compute: 
                  \code{"Eucl"} for Euclidean, 
                  \code{"Chor"} for chordal, 
                  \code{"Geod"} for geodesic.}
  \item{radius}{Radius of the sphere for geodesic or chordal distances. 
                Defaults to 1. Use Earth radius (e.g., 6371 Km) for geographic distances.}
}
\details{
\itemize{
  \item Euclidean distance: straight-line distance in 2D or 3D space.
  \item Chordal distance: Euclidean distance between points projected on the 
        unit sphere, scaled by the given radius.
  \item Geodesic distance: shortest path along the surface of a sphere of 
        radius \code{radius}.
}
}
\value{
A symmetric numeric matrix of size \code{n x n}, where \code{n} is the number 
of points. Each entry contains the distance between the corresponding pair 
of coordinates.
}
\examples{
# Example with Euclidean distance
coords <- cbind(c(0, 10), c(0, 10))
GeoDistances(coords, distance = "Eucl")

# Example with geodesic distance (approx distance between Rome and New York)
rome <- c(12.4964, 41.9028)     #
ny   <- c(-74.0060, 40.7128)
coords <- rbind(rome, ny)
GeoDistances(coords, distance = "Geod", radius = 6371)
}

\keyword{Distances}
