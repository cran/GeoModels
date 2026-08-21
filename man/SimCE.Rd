\name{SimCE}
\alias{SimCE}
\encoding{UTF-8}
\title{Circulant embedding simulation on a regular spatial grid}
\description{
Internal subroutine used by \code{GeoSimapprox} to simulate a Gaussian
random field on an exact, regular two-dimensional grid by circulant
embedding.
}
\usage{
SimCE(M, N, x, y, z = NULL, corrmodel, param,
      mean.val = 0, max.ext = 1L)
}
\arguments{
 \item{M}{Integer; number of coordinates in \code{x}.}
 \item{N}{Integer; number of coordinates in \code{y}.}
 \item{x}{A finite, strictly increasing, equally spaced numeric vector of length \code{M}.}
 \item{y}{A finite, strictly increasing, equally spaced numeric vector of length \code{N}.}
 \item{z}{Must be \code{NULL}. Three-dimensional CE grids are not implemented.}
 \item{corrmodel}{String; the name of a spatial correlation model.}
 \item{param}{A list of correlation parameters. \code{sill} is the marginal variance and \code{nugget} is the nugget proportion.}
 \item{mean.val}{A finite scalar or a vector with one value per grid point.}
 \item{max.ext}{Positive integer; maximum number of successive doubling attempts for the embedding.}
}
\details{
The simulated coordinates are exactly \code{expand.grid(x, y)}. The method
uses Euclidean distance. If the embedding still has substantial negative
eigenvalues at the maximum extension, they are truncated to zero and a
warning is issued.
}
\value{
A list containing \code{X}, the simulated values, \code{grid.points}, the
exact grid coordinates, \code{k0}, the selected extension attempt,
\code{approx.flag}, and the minimum embedding eigenvalue.
}
\seealso{\code{\link{GeoSimapprox}}}
\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},
Víctor Morales Oñate, \email{victor.morales@uv.cl},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl}}
\keyword{Simulation}
