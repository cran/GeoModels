\name{GeoOutlier}
\alias{GeoOutlier}
\encoding{UTF-8}
\title{Spatio (temporal) outliers detection}
\description{
Given a set of spatio (temporal) locations and data, the procedure select the spatial or spatiotemporal ouliers
using a specific algorithm.
}
\usage{
GeoOutlier(data, coordx, coordy=NULL,coordz=NULL, coordt=NULL, coordx_dyn=NULL, 
             distance="Eucl", grid=FALSE,  neighb=10,alpha=0.001,
             method="Z-Median", radius=1, bivariate=FALSE,X=NULL)
}
\arguments{
  \item{data}{An optional \eqn{d}{d}-dimensional vector (a single spatial realisation)  or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid)
   or a (\eqn{t \times d}{t x d})-matrix (a single spatial-temporal realisation)   or an (\eqn{d \times d \times t \times n }{d x d x t})-array
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

  \item{distance}{String; the name of the spatial distance. The default
    is \code{Eucl}, the euclidean distance. See the Section
    \bold{Details}  of \code{\link{GeoFit}}.}

   \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as spatial or spatial-temporal realisations on a set
    of non-equispaced spatial sites (irregular grid).}

   \item{neighb}{Numeric; an optional positive integer indicating the 
    order of neighborhoodused for Z-Median algorithm.}

    \item{alpha}{Numeric; a numeric value between 0 and 1 used for Z-Median algorithm.}

  \item{method}{String; The name of the algorithm for detecting spatial ouliers. Default is Z-median proposed in Chen et al. (2008)}
\item{radius}{Numeric; a value indicating  the radius of the sphere when using the great 
    circle distance. Default value is 1.}  

  \item{bivariate}{If TRUE then data  is considered as spatial  bivariate data.}

  \item{X}{Numeric; an optional Matrix of spatio (temporal) covariates.}
}


\value{
  Return a  matrix or a list  containing the dected spatial or spatio-temporal outliers
}

\references{
 Chen D, Lu C, Kou Y, Chen F (2008) On detecting spatial outliers. Geoinformatica 12:455–475
 
 Bevilacqua M., Caamaño C., Arellano-Valle R. B., Camilo Gomez C. (2022) A class of random fields with two-piece marginal distributions for modeling point-referenced data with spatial outliers. \emph{Test}  10.1007/s11749-021-00797-5
}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}


\examples{

library(GeoModels)
set.seed(1428)
NN = 1500
coords = cbind(runif(NN),runif(NN))
###
scale=0.5/3
corrmodel = "Matern"; 

param = list(mean=0,sill=1,nugget=0,scale=scale,smooth=0.5,skew=0)
data = GeoSim(coordx = coords,corrmodel = corrmodel,
                  model = "TwoPieceGaussian",param = param)$data

K=15         #parameter for outliers detection alghoritm
alpha=0.005  #parameter for outliers detection alghoritm
outlier=GeoOutlier(data=data, coordx = coords,neighb=K,alpha=alpha)
quilt.plot(coords,data)
for (i in 1:nrow(outlier))  plotrix::draw.circle(outlier[i,1], outlier[i,2],radius=0.02,lwd=2) 
nrow(outlier) # number of outliers

param = list(mean=0,sill=1,nugget=0.4,scale=scale,smooth=0.5)
data = GeoSim(coordx = coords,corrmodel = corrmodel,
                  model = "Gaussian",param = param)$data

K=15         #parameter for outliers detection alghoritm
alpha=0.005  #parameter for outliers detection alghoritm
outlier=GeoOutlier(data=data, coordx = coords,neighb=K,alpha=alpha)
quilt.plot(coords,data)
for (i in 1:nrow(outlier))  plotrix::draw.circle(outlier[i,1], outlier[i,2],radius=0.02,lwd=2)
nrow(outlier) # number of outliers
}
