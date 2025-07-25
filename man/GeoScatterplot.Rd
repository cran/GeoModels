\name{GeoScatterplot}
\alias{GeoScatterplot}
\encoding{UTF-8}
\title{h-scatterplot for space and space-time data.}
\description{
  The function produces h-scatterplots for the spatial, spatio-temporal and  bivariate setting. 
}
\usage{
GeoScatterplot(data, coordx, coordy=NULL, coordz=NULL, coordt=NULL, coordx_dyn=NULL,
           distance='Eucl', grid=FALSE, maxdist=NULL,neighb=NULL,
           times=NULL, numbins=4, radius=1, bivariate=FALSE,\dots)
}
\arguments{
  \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation) or a (\eqn{n \times d}{n x d})-matrix
   (\eqn{n} iid spatial realisations) or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid)
   or an (\eqn{d \times d \times n}{d x d x n})-array (\eqn{n} iid spatial realisations on regular grid) or a
   (\eqn{t \times d}{t x d})-matrix (a single spatial-temporal realisation) or an (\eqn{t \times d \times n }{t x d x n})-array
   (\eqn{n} iid spatial-temporal realisations) or or an (\eqn{d \times d \times t \times n }{d x d x t})-array
   (a single spatial-temporal realisation on regular grid) or an (\eqn{d \times d \times t \times n }{d x d x t x n})-array
   (\eqn{n} iid spatial-temporal realisations on regular grid). See \code{\link{GeoFit}} for details.}
   \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix  or  (\eqn{d \times 3}{d x 3})-matrix 
     Coordinates on a sphere for a  fixed radius \code{radius} 
    are passed in lon/lat format expressed in decimal degrees.}
  \item{coordy}{A numeric vector giving 1-dimension of
    spatial coordinates; Optional argument, the default is \code{NULL}.}
      \item{coordz}{A numeric vector giving 1-dimension of
    spatial coordinates;  Optional argument, the default is \code{NULL}.}

    \item{coordt}{A numeric vector assigning 1-dimension of
    temporal coordinates. Optional argument, the default is \code{NULL} then a spatial random field is expected.}
   \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d x 2})-matrices
       containing dynamical (in time) spatial coordinates. Optional argument, the default is \code{NULL}
    }   
  \item{distance}{String; the name of the spatial distance. The default
    is \code{Eucl}, the euclidean distance. See the Section
    \bold{Details} of \code{\link{GeoFit}}.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as spatial or spatial-temporal realisations on a set of 
    non-equispaced spatial sites.}
  \item{maxdist}{A numeric value denoting the spatial maximum distance,
   see the Section \bold{Details}.}
    \item{neighb}{Numeric; an optional positive integer indicating the 
    order of neighborhood. See the Section \bold{Details} for more information.}
   \item{times}{A numeric vector denoting the temporal instants  involved \bold{Details}.}
  \item{numbins}{A numeric value denoting the numbers of bins, see the
    Section \bold{Details}.}
  \item{radius}{Numeric; a value indicating  the radius of the sphere when 
  using the great circle distance. Default value is 1.} 
    \item{bivariate}{Logical; if \code{FALSE} (the default) the data
    are interpreted as univariate spatial or spatial-temporal realisations. 
    Otherwise they are intrepreted as a a realization from a bivariate field.}
    \item{\dots}{Optional parameters passed to the plot function.}
}

\value{Produces a plot. No values are returned.}

\details{
  h-scatterplot is the plot of the pair values  that are neighborhood of a certain order or  with distances belonging to a certain interval.
  In the first case a (vector of)  neighborhood must  be specified. In the second case a maximum distance (maxdist) and a number of lag-bins
  (numbins) must be specified.
  The method based on neighborhoods is recommended in particular for large datasets.

}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}


\examples{

library(GeoModels)
set.seed(514)

NN = 600
coords = cbind(runif(NN),runif(NN))

param = list(mean=0,sill=1,nugget=0,power2=4,scale=0.4,smooth=0)

corrmodel = "GenWend"; model = "Gaussian"

data = GeoSim(coordx = coords,corrmodel = corrmodel,
                  model = model,param = param)$data

# h-scatterplots for given a vector of neighborhoods
GeoScatterplot(data,coords,neighb=c(2,4))

}
\keyword{Scatterplot}
