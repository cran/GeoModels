\name{GeoNeighborhood}
\alias{GeoNeighborhood}
\encoding{UTF-8}
\title{Spatio (temporal) neighborhood selection for local kriging.}
\description{
Given a set of spatio (temporal) locations and data, the procedure selects a spatio (temporal) neighborhood associated to some  given spatio (temporal) locations.
The neighborhood is computed  using a fixed  spatio (temporal) threshold
or considering a fixed number of spatio (temporal) neighbors.
}

\usage{
GeoNeighborhood(data=NULL, coordx, coordy=NULL,coordz=NULL,
coordt=NULL, coordx_dyn=NULL, bivariate=FALSE,
               distance="Eucl", grid=FALSE, 
               loc, neighb=NULL,maxdist=NULL,
               maxtime=NULL, radius=1, time=NULL, 
               X=NULL,M=NULL,spobj=NULL,spdata=NULL,
               parallel=FALSE,ncores=NULL)
}
\arguments{
  \item{data}{An optional \eqn{d}{d}-dimensional vector (a single spatial realisation)  or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid)
   or a
   (\eqn{t \times d}{t x d})-matrix (a single spatial-temporal realisation)   or an (\eqn{d \times d \times t \times n }{d x d x t})-array
   (a single spatial-temporal realisation on regular grid).}
        \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix  or  (\eqn{d \times 3}{d x 3})-matrix 
     Coordinates on a sphere for a  fixed radius \code{radius} 
    are passed in lon/lat format expressed in decimal degrees.}
  \item{coordy}{A numeric vector giving 1-dimension of
    spatial coordinates; Optional argument, the default is \code{NULL}.}
      \item{coordz}{A numeric vector giving 1-dimension of
    spatial coordinates;  Optional argument, the default is \code{NULL}.}
     \item{coordt}{A numeric vector giving 1-dimension of
    temporal coordinates.  Optional argument, the default is \code{NULL}
    then a spatial RF is expected.}
  \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d x 2})-matrices
       containing dynamical (in time) spatial coordinates. Optional argument, the default is \code{NULL}
    }
   \item{bivariate}{If TRUE then data  is considered as spatial  bivariate data.}
  \item{distance}{String; the name of the spatial distance. The default
    is \code{Eucl}, the euclidean distance. See the Section
    \bold{Details}  of \code{\link{GeoFit}}.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as spatial or spatial-temporal realisations on a set
    of non-equispaced spatial sites (irregular grid).}
  \item{loc}{A (\eqn{1 \times 2}{1 x 2})-matrix  giving the spatial coordinate
     of the location for which a neighborhood is computed .}
 \item{neighb}{Numeric; an optional positive integer indicating the 
    order of spatial neighborhood.}
 \item{maxdist}{Numeric; a positive value indicating the maximum
    spatial distance considered in the spatial neighborhood
    selection.}
  \item{maxtime}{Numeric; an optional positive integer indicating the 
    order of temporal  neighborhood.}
    \item{radius}{Numeric; a value indicating  the radius of the sphere when using the great 
    circle distance. Default value is 1.}  
   \item{time}{Numeric; a value  giving the temporal instant for which a neighborhood is computed.}
  \item{X}{Numeric; an optional Matrix of spatio (temporal) covariates.}
   \item{M}{Numeric; an estimated spatio (temporal) mean vector.}
    \item{spobj}{An object of class sp or spacetime}
   \item{spdata}{Character:The name of data in the sp or spacetime object}
    \item{parallel}{Logical; if \code{TRUE} then  parallelizzation is used}
    \item{ncores}{Numeric; number of cores involved in parallelization.} 
}


\value{
  Returns a list  containing the following informations:
  \item{coordx}{A list  of  the  matrix coordinates of the computed spatial neighborhood ;}
  \item{coordt}{A vector  of the computed temporal neighborhood;}
  \item{data}{A list  of the vector of data associated with the spatio  (temporal) neighborhood;}
  \item{distance}{The type of spatial distance;}
  \item{numcoord}{The vector of numbers of location sites involved the spatial neighborhood;}
  \item{numtime}{The vector of numbers of temporal insttants involved the temporal neighborhood;}
  \item{radius}{The radius of the sphere if coordinates are passed in lon/lat format;}
  \item{spacetime}{\code{TRUE} if spatio-temporal and \code{FALSE} if spatial RF;}
  \item{X}{The matrix of spatio  (temporal) covariates associated with the computed spatio  (temporal) neighborhood;}
}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}


\examples{
library(GeoModels)
##########################################
#### Example: spatial neighborhood  ######
##########################################
set.seed(75)
coords=cbind(runif(500),runif(500))

param=list(nugget=0,mean=0,scale=0.2,sill=1,
            power2=4,smooth=1)

data_all = GeoSim(coordx=coords, corrmodel="GenWend", 
                         param=param)$data

plot(coords)
##two locations 
loc_to_pred=matrix(c(0.3,0.5,0.7,0.2),2,2)

points(loc_to_pred,pch=20)
neigh=GeoNeighborhood(data_all, coordx=coords,  
                  loc=loc_to_pred,neighb=8)

# two Neighborhoods 
neigh$coordx
points(neigh$coordx[[1]],pch=20,col="red")
points(neigh$coordx[[2]],pch=20,col="blue")
# associated data
neigh$data


###################################################
#### Example: spatio temporal spatial neighborhood#  
###################################################

set.seed(78)
coords=matrix(runif(80),40,2)
coordt=seq(0,6,0.25)

param=list(nugget=0,mean=0,scale_s=0.2/3,scale_t=0.25/3,sill=2)

data_all = GeoSim(coordx=coords, coordt=coordt,corrmodel="Exp_Exp", 
                         param=param)$data
##  two location to predict
loc_to_pred=matrix(runif(4),2,2)
## three temporal instants to predict
time=c(1,2)

plot(coords,xlim=c(0,1),ylim=c(0,1))
points(loc_to_pred,pch=20)

neigh=GeoNeighborhood(data_all, coordx=coords,  coordt=coordt,
                  loc=loc_to_pred,time=time,neighb=3,maxtime=0.5)

# first spatio-temporal neighborhoods 
# with  associated data
neigh$coordx[[1]]
neigh$coordt[[1]]
neigh$data[[1]]

plot(coords)
points(loc_to_pred,pch=20)
points(neigh$coordx[[1]],col="red",pch=20)

###################################################
#### Example: bivariate  spatial neighborhood #####  
###################################################

set.seed(79)
coords=matrix(runif(100),50,2)

param=list(mean_1=0,mean_2=0,scale=0.12,smooth=0.5,
           sill_1=1,sill_2=1,nugget_1=0,nugget_2=0,pcol=0.5)

data_all = GeoSim(coordx=coords,corrmodel="Bi_matern_sep",
                 param=param)$data
##  two location to predict
loc_to_pred=matrix(runif(4),2,2)

neigh=GeoNeighborhood(data_all, coordx=coords,bivariate=TRUE,
                  loc=loc_to_pred,neighb=5)

plot(coords)
points(loc_to_pred,pch=20)
points(neigh$coordx[[1]],col="red",pch=20)
points(neigh$coordx[[2]],col="red",pch=20)

}
