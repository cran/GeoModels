\name{GeoSim}
\alias{GeoSim}
\alias{print.GeoSim}
\encoding{UTF-8}
\title{Simulation of Gaussian and non Gaussian Random Fields.}
\description{
  Simulation of Gaussian and some non Gaussian spatial, spatio-temporal and spatial bivariate 
 random fields. The function return a realization of a random Field for a given covariance model 
 and covariance parameters.Simulation is based on Cholesky decomposition.
}
\usage{
GeoSim(coordx, coordy=NULL,coordz=NULL, coordt=NULL, coordx_dyn=NULL, corrmodel, 
      distance="Eucl", grid=FALSE, method="cholesky", 
      model='Gaussian', n=1, param,anisopars=NULL,radius=1, 
      sparse=FALSE,X=NULL,spobj=NULL,nrep=1,progress=TRUE)
}
\arguments{
   \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix  or 
   (\eqn{d \times 3}{d x 3})-matrix. Coordinates on a sphere for a  fixed radius \code{radius} 
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
  \item{corrmodel}{String; the name of a correlation model, for the
    description see the Section \bold{Details}.}
  \item{distance}{String; the name of the spatial distance. The default
    is \code{Eucl}, the euclidean distance. See the Section
    \bold{Details}  of \code{\link{GeoFit}}.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as spatial or spatial-temporal realisations on a set
    of non-equispaced spatial sites (irregular grid).}
     \item{method}{String; the type of matrix decomposition used in the simulation. 
       Default is cholesky. The other possible choices is \code{svd}.} 
  \item{model}{String; the type of RF and therefore the densities associated to the likelihood
    objects. \code{Gaussian} is the default, see the Section
    \bold{Details}.}
  \item{n}{Numeric; the number of trials for binomial RFs. 
           The number of successes in the negative Binomial RFs. Default is \eqn{1}.} 
  \item{param}{A list of parameter values required in the
    simulation procedure of RFs, see \bold{Examples}.}
         \item{anisopars}{ A list of two elements "angle" and "ratio" i.e. the anisotropy angle and the anisotropy ratio, respectively.}
    \item{radius}{Numeric; a value indicating  the radius of the sphere when using the great 
    circle distance. Default value is 1.}  
  \item{sparse}{Logical; if \code{TRUE} then  cholesky decomposition is performed
  using sparse matrices algorithms (spam packake). It should be used
  with compactly supported covariance models.\code{FALSE} is the default.}
  \item{X}{Numeric; Matrix of space-time covariates.}
  \item{spobj}{An object of class sp or spacetime}
   \item{nrep}{Numeric; Numbers of indipendent replicates.}
      \item{progress}{Logic; If TRUE then a progress bar is shown.}
  
}


\value{
  Returns an object of class \code{GeoSim}.
  An object of class \code{GeoSim} is a list containing
  at most the following components:
  
  \item{bivariate}{Logical:\code{TRUE} if the Gaussian RF is bivariate, otherwise  \code{FALSE};}
  \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordz}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates;}
  \item{coordx_dyn}{A  list of dynamical (in time) spatial coordinates;}
  \item{corrmodel}{The correlation model; see \code{\link{GeoCovmatrix}}.}
  \item{data}{The vector or matrix or array of data, see
    \code{\link{GeoFit}};}
  \item{distance}{The type of spatial distance;}
    \item{method}{The method of simulation}
  \item{model}{The type of RF, see \code{\link{GeoFit}}.}
  \item{n}{The number of trial for Binomial RFs;the number of successes in a negative Binomial RFs;} 
  \item{numcoord}{The number of spatial coordinates;}
  \item{numtime}{The number the temporal realisations of the RF;}
  \item{param}{The vector of  parameters' estimates;}
  \item{radius}{The radius of the sphere if coordinates are passed in lon/lat format;}
  \item{spacetime}{\code{TRUE} if spatio-temporal and \code{FALSE} if
  spatial RF;}
  \item{nrep}{The number of indipendent replicates;}
}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}


\examples{
library(GeoModels)
library(mapproj)


################################################################
###
### Example 1. Simulation of a spatial Gaussian RF 
### with Matern and Generalized Wendland correlations
###############################################################

# Define the spatial-coordinates of the points:
x <- runif(500);y <- runif(500)
coords=cbind(x,y)
set.seed(261)
# Simulation of a spatial Gaussian RF with Matern correlation function
data1 <- GeoSim(coordx=coords, corrmodel="Matern", param=list(smooth=0.5,
             mean=0,sill=1,scale=0.4/3,nugget=0))$data

set.seed(261)
data2 <- GeoSim(coordx=coords,  corrmodel="GenWend", param=list(smooth=0,
              power2=4,mean=0,sill=1,scale=0.4,nugget=0))$data
opar=par(no.readonly = TRUE)
par(mfrow=c(1,2))
quilt.plot(coords,data1,main="Matern",xlab="",ylab="")
quilt.plot(coords,data2,main="Wendland",xlab="",ylab="")   
par(opar)
             

################################################################
###
### Example 2. Simulation of a spatial geometric RF 
### with  underlying Wend0 correlation
###
################################################################

# Define the spatial-coordinates of the points:
x <- runif(800);y <- runif(800)
coords <- cbind(x,y)
set.seed(251)
# Simulation of a spatial Binomial RF:
sim <- GeoSim(coordx=coords, corrmodel="Wend0",
             model="BinomialNeg",n=1,sparse=TRUE,
             param=list(nugget=0,mean=0,scale=.2,power2=4))

quilt.plot(coords,sim$data,nlevel=max(sim$data),col=terrain.colors(max(sim$data+1))) 

################################################################
###
### Example 3. Simulation of a spatial Weibull  RF
### with  underlying Matern correlation on a regular grid
###
###############################################################
# Define the spatial-coordinates of the points:
x <- seq(0,1,0.032)
y <- seq(0,1,0.032)
set.seed(261)
# Simulation of a spatial Gaussian RF with Matern correlation function
data1 <- GeoSim(x,y,grid=TRUE, corrmodel="Matern",model="Weibull", 
         param=list(shape=1.2,mean=0,scale=0.3/3,nugget=0,smooth=0.5))$data
image.plot(x,y,data1,main="Weibull RF",xlab="",ylab="")

################################################################
###
### Example 4. Simulation of a spatial t  RF
### with  with  underlying Generalized Wendland correlation 
###
###############################################################
# Define the spatial-coordinates of the points:
x <- seq(0,1,0.03)
y <- seq(0,1,0.03)
set.seed(268)
# Simulation of a spatial Gaussian RF with Matern correlation function
data1 <- GeoSim(x,y,grid=TRUE, corrmodel="GenWend",model="StudentT", sparse=TRUE,
         param=list(df=1/4,mean=0,sill=1,scale=0.3,nugget=0,smooth=1,power2=5))$data
image.plot(x,y,data1,col=terrain.colors(100),main="Student-t RF",xlab="",ylab="")


################################################################
###
### Example 5. Simulation of a sinhasinh RF
###   with  underlying Wend0 correlation.
###
###############################################################

# Define the spatial-coordinates of the points:
x <- runif(500, 0, 2)
y <- runif(500, 0, 2)
coords <- cbind(x,y)
set.seed(261)
corrmodel="Wend0"
# Simulation of a spatial Gaussian RF:
param=list(power2=4,skew=0,tail=1,
             mean=0,sill=1,scale=0.2,nugget=0)  ## gaussian case
data0 <- GeoSim(coordx=coords, corrmodel=corrmodel,
               model="SinhAsinh", param=param,sparse=TRUE)$data
plot(density(data0),xlim=c(-7,7))

param=list(power2=4,skew=0,tail=0.7,
             mean=0,sill=1,scale=0.2,nugget=0) ## heavy tails
data1 <- GeoSim(coordx=coords, corrmodel=corrmodel,
               model="SinhAsinh", param=param,sparse=TRUE)$data
lines(density(data1),lty=2)

param=list(power2=4,skew=0.5,tail=1,
             mean=0,sill=1,scale=0.2,nugget=0)  ## asymmetry
data2 <- GeoSim(coordx=coords, corrmodel=corrmodel,
               model="SinhAsinh", param=param,sparse=TRUE)$data
lines(density(data2),lty=3)

################################################################
###
### Example 6. Simulation of a bivariate Gaussian RF
### with  bivariate Matern correlation model
###
###############################################################

# Define the spatial-coordinates of the points:
x <- runif(500, 0, 2)
y <- runif(500, 0, 2)
coords <- cbind(x,y)

# Simulation of a bivariate spatial Gaussian RF:
# with a separable Bivariate Matern
param=list(mean_1=4,mean_2=2,smooth_1=0.5,smooth_2=0.5,smooth_12=0.5,
           scale_1=0.12,scale_2=0.1,scale_12=0.15,
           sill_1=1,sill_2=1,nugget_1=0,nugget_2=0,pcol=0.5)
data <- GeoSim(coordx=coords,corrmodel="Bi_matern",
              param=param)$data
opar=par(no.readonly = TRUE)
par(mfrow=c(1,2))
quilt.plot(coords,data[1,],col=terrain.colors(100),main="1",xlab="",ylab="")
quilt.plot(coords,data[2,],col=terrain.colors(100),main="2",xlab="",ylab="")
par(opar)


################################################################
###
### Example 7. Simulation of a  spatio temporal Gaussian RF.
### observed on  fixed  location sites with double Matern correlation 
###
###############################################################



coordt=1:5

# Define the spatial-coordinates of the points:
x <- runif(50, 0, 2)
y <- runif(50, 0, 2)
coords <- cbind(x,y)

param<-list(nugget=0,mean=0,scale_s=0.2/3,scale_t=2/3,sill=1,smooth_s=0.5,smooth_t=0.5)
data <- GeoSim(coordx=coords, coordt=coordt, corrmodel="Matern_Matern",
                     param=param)$data
dim(data)

################################################################
###
### Example 8. Simulation of a  spatio temporal Gaussian RF.
### observed on  dynamical location sites with double Matern correlation 
###
###############################################################

# Define the dynamical spatial-coordinates of the points:

coordt=1:5
coordx_dyn=list()
maxN=30
set.seed(8)
for(k in 1:length(coordt))
{
NN=sample(1:maxN,size=1)
x <- runif(NN, 0, 1)
y <- runif(NN, 0, 1)
coordx_dyn[[k]]=cbind(x,y)
}
coordx_dyn

param<-list(nugget=0,mean=0,scale_s=0.2/3,scale_t=2/3,sill=1,smooth_s=0.5,smooth_t=0.5)
data <- GeoSim(coordx_dyn=coordx_dyn, coordt=coordt, corrmodel="Matern_Matern",
                     param=param)$data
## spatial realization at first temporal instants
data[[1]]
## spatial realization at third temporal instants
data[[3]]




################################################################
###
### Example 9. Simulation of a Gaussian RF 
###  with a Wend0 correlation in the north emisphere of the planet earth
### using geodesic distance
###############################################################
distance="Geod";radius=6371

NN=3000 ## total point on the sphere on lon/lat format
set.seed(80)
coords=cbind(runif(NN,-180,180),runif(NN,0,90))
## Set the wendland parameters
corrmodel <- "Wend0"
param<-list(mean=0,sill=1,nugget=0,scale=1000,power2=3)
# Simulation of a spatial Gaussian RF on the sphere
#set.seed(2)
data <- GeoSim(coordx=coords,corrmodel=corrmodel,sparse=TRUE,
               distance=distance,radius=radius,param=param)$data
#require(globe)
#globe::globeearth(eye=place("newyorkcity"))
#globe::globepoints(loc=coords,pch=20,col =  cm.colors(length(data),alpha=0.4)[rank(data)])




}
\keyword{Simulation}
