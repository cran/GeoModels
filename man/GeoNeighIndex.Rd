\name{GeoNeighIndex}
\alias{GeoNeighIndex}
\encoding{UTF-8}
\title{Spatial or spatiotemporal  near neighbour indices.}
\description{
  The function returns the indices associated with  a given spatial (temporal) neighbour and/or distance 
}
\usage{
GeoNeighIndex(coordx,coordy=NULL,coordz=NULL,
coordt=NULL,coordx_dyn=NULL,distance="Eucl",neighb=4,maxdist=NULL,
maxtime=1,radius=1,bivariate=FALSE)
}
\arguments{
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
\item{neighb}{Numeric; an optional (vector of) positive integer indicating the 
    order of neighborhood. See the Section \bold{Details} for more information.}
 \item{maxdist}{A numeric value denoting the spatial  distance \bold{Details}.}
   \item{maxtime}{A numeric value denoting the temporal distance \bold{Details}.}
  \item{radius}{Numeric; a value indicating  the radius of the sphere when 
  using the great circle distance. Default value is 1.} 
    \item{bivariate}{Logical; if \code{FALSE} (the default) the data
    are interpreted as univariate spatial or spatial-temporal realisations. 
    Otherwise they are intrepreted as a a realization from a bivariate field.}
}

\details{
 The function returns the  spatial or spatiotemporal indices of the pairs tha are neighboords of a certain order and/or with a certain fixed distance
}


\value{
  Returns a list containing the following components:

  \item{colidx}{First vector of indices}
  \item{rowidx}{Second vector of indices}
  \item{lags}{Vector of spatial distances}
  \item{lagt}{Vector of temporal distances}
}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}


\examples{


require(GeoModels)
NN = 400
coords = cbind(runif(NN),runif(NN))
scale=0.5/3
param = list(mean=0,sill=1,nugget=0,scale=0.5/3,smooth=0.5)
corrmodel = "Matern"; 

param = list(mean=0,sill=1,nugget=0,scale=scale,smooth=0.5)
set.seed(951)
data = GeoSim(coordx = coords,corrmodel = corrmodel,
                  model = "Gaussian",param = param)$data

sel=GeoNeighIndex(coordx=coords,neighb=5)        

data1=data[sel$colidx]; data2=data[sel$rowidx]
## plotting pairs  that are neighboord of order 5
plot(data1,data2,xlab="",ylab="",main="h-scatterplot, neighb=5")

}
\keyword{Scatterplot}
