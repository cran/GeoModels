\name{GeoNeighIndex}
\alias{GeoNeighIndex}
\encoding{UTF-8}
\title{Spatial or spatiotemporal near neighbour indices.}
\description{
  The function returns the indices associated with a given spatial (temporal) neighbour
  and/or distance. Optionally, a stochastic thinning mechanism can be applied to randomly
  retain only a fraction of the indices.
}
\usage{
GeoNeighIndex(coordx, coordy=NULL, coordz=NULL, coordt=NULL,
              coordx_dyn=NULL, distance="Eucl", neighb=4,
              maxdist=NULL, maxtime=1, radius=1,
              bivariate=FALSE, p_neighb=1, thin_method="iid")
}
\arguments{
  \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix or (\eqn{d \times 3}{d x 3})-matrix.
    Coordinates on a sphere for a fixed radius \code{radius} are passed in lon/lat format
    expressed in decimal degrees.}

  \item{coordy}{A numeric vector giving 1-dimension of spatial coordinates; optional argument,
    default is \code{NULL}.}

  \item{coordz}{A numeric vector giving 1-dimension of spatial coordinates; optional argument,
    default is \code{NULL}.}

  \item{coordt}{A numeric vector assigning 1-dimension of temporal coordinates. Optional argument,
    default is \code{NULL}; if \code{NULL} then a spatial random field is expected.}

  \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d_t x 2})-matrices containing dynamical
    (in time) spatial coordinates. Optional argument, default is \code{NULL}.}

  \item{distance}{String; the name of the spatial distance. Default is \code{"Eucl"} (euclidean distance).
    See the Section \bold{Details} of \code{\link{GeoFit}}.}

  \item{neighb}{Numeric; an optional (vector of) positive integer indicating the order of neighborhood.
    See the Section \bold{Details} for more information.}

  \item{maxdist}{A numeric value denoting the maximum spatial distance; see \bold{Details}.}

  \item{maxtime}{A numeric value denoting the maximum temporal distance; see \bold{Details}.}

  \item{radius}{Numeric; a value indicating the radius of the sphere when using great-circle distances.
    Default value is \code{1}.}

  \item{bivariate}{Logical; if \code{FALSE} (default) the data are interpreted as univariate spatial or
    spatio-temporal realisations. Otherwise, they are interpreted as a realization from a bivariate field.}

  \item{p_neighb}{Numeric; a value in \eqn{(0,1]} specifying the expected fraction of nearest-neighbor pairs
    retained through stochastic thinning. If \code{1} (default), no thinning is applied and all nearest-neighbor
    pairs are used. If \code{<1}, pairs are randomly retained using independent Bernoulli sampling.}

  \item{thin_method}{String; thinning scheme used when \code{p_neighb<1}. Default is \code{"iid"}(independent Bernoulli thinning).}
}
\details{
  The function returns spatial or spatio-temporal indices of the pairs that are neighbours of a given order
  and/or within a given distance threshold.
}
\value{
  Returns a list containing the following components:
  \item{colidx}{First vector of indices.}
  \item{rowidx}{Second vector of indices.}
  \item{lags}{Vector of spatial distances.}
  \item{lagt}{Vector of temporal distances.}
}
\author{
Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com}, \url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl}, \url{https://www.researchgate.net/profile/Christian-Caamano}
}
\examples{
require(GeoModels)
NN <- 400
coords <- cbind(runif(NN), runif(NN))

corrmodel <- "Matern"
scale <- 0.5/3
param <- list(mean=0, sill=1, nugget=0, scale=scale, smooth=0.5)

set.seed(951)
data <- GeoSim(coordx=coords, corrmodel=corrmodel,
              model="Gaussian", param=param)$data

sel <- GeoNeighIndex(coordx=coords, neighb=5)

data1 <- data[sel$colidx]
data2 <- data[sel$rowidx]

## plotting pairs that are neighbours of order 5
plot(data1, data2, xlab="", ylab="",
     main="h-scatterplot, neighb=5")
}
\keyword{Scatterplot}
