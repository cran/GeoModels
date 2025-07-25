\name{GeoCorrFct}
\alias{GeoCorrFct}
\encoding{UTF-8}
\title{Spatial and Spatio-temporal correlation or  covariance of (non) Gaussian random fields}
\description{
  The function computes the  correlations of 
  a spatial (or spatio-temporal or  bivariate spatial) Gaussian or non Gaussian randomm field
  for a  given  correlation  model and a  set of spatial (temporal) distances.
}
\usage{
GeoCorrFct(x,t=NULL,corrmodel, model="Gaussian",
distance="Eucl", param, radius=6371,n=1,
covariance=FALSE,variogram=FALSE)
}

\arguments{
  \item{x}{A set of spatial distances.}
  \item{t}{A set of (optional) temporal distances.}
  \item{corrmodel}{String; the name of a correlation model, for the
    description see \code{\link{GeoCovmatrix}}.}
  \item{model}{String; the type of RF. See \code{\link{GeoFit}}.}
  \item{distance}{String; the name of the spatial distance. The default
  is \code{Eucl}, the euclidean distance. See \code{\link{GeoFit}}.}
  \item{param}{A list of parameter values required  for the covariance
  model.}
  \item{radius}{Numeric; a value indicating  the radius of the sphere 
     when using covariance models valid using the great circle distance. 
     Default value is the radius of the earth in Km (i.e. 6371)} 
  \item{n}{Numeric; the number of trials in a (negative) binomial random fields.
         Default is \eqn{1}.}
  \item{covariance}{Logic; if TRUE then the covariance is returned. Default is FALSE}
   \item{variogram}{Logic; if FALSE then the covariance/correlation is returned. Otherwise the associated
   semivariogram is returned}
   }


\value{Returns  correlations or covariances values associated to a given parametric spatial and temporal correlation models.}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian, Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{
library(GeoModels)

################################################################
###
### Example 1. Covariance of a Gaussian random field with underlying 
### Matern correlation model with nugget
###
###############################################################
# Define the spatial distances
x = seq(0,1,0.002)
# Correlation Parameters for Matern model 
CorrParam("Matern")
NuisParam("Gaussian")
# Matern Parameters 
param=list(sill=2,smooth=0.5,scale=0.2/3,nugget=0.2,mean=0)
cc= GeoCorrFct(x=x, corrmodel="Matern", covariance=TRUE,
  param=param,model="Gaussian")
plot(cc,ylab="Corr",lwd=2,main="Matern correlation",type="l")

################################################################
###
### Example 2. Covariance of a Gaussian random field with underlying 
### Generalized Wendland-Matern correlation model
###
###############################################################
CorrParam("GenWend_Matern")
NuisParam("Gaussian")
# GenWend Matern Parameters 
param=list(sill=2,smooth=1,scale=0.1,nugget=0,power2=1/4,mean=0)
cc= GeoCorrFct(x=x, corrmodel="GenWend_Matern", param=param,model="Gaussian",covariance=FALSE)
plot(cc,ylab="Cov",lwd=2,,main="GenWend covariance",type="l")

################################################################
###
### Example 3. Semivariogram of a Tukeyh random field with underlying 
### Generalized Wendland correlation model
###
###############################################################
CorrParam("GenWend")
NuisParam("Tukeyh")
x = seq(0,1,0.005)
param=list(sill=1,smooth=1,scale=0.5,nugget=0,power2=5,tail=0.1,mean=0)
cc= GeoCorrFct(x=x, corrmodel="GenWend", param=param,model="Tukeyh",variogram=TRUE)
plot(cc,ylab="Corr",lwd=2,main="Tukey semivariogram",type="l")

################################################################
###
### Example 4. Semi-Variogram of a LoggGaussian random field with underlying 
### Kummer correlation model
###
###############################################################
CorrParam("Kummer")
NuisParam("LogGaussian")
# GenWend Matern Parameters 
param=list(smooth=1,sill=0.5,scale=0.1,nugget=0,power2=1,mean=0)
cc= GeoCorrFct(x=x, corrmodel="Kummer", param=param,model="LogGaussian",
       ,covariance=TRUE,variogram=TRUE)
plot(cc,ylab="Semivario",lwd=2,
  main="LogGaussian semivariogram",type="l")


################################################################
###
### Example 5. Covariance of Poisson random field with underlying 
### Matern correlation model
###
###############################################################
CorrParam("Matern")
NuisParam("Poisson")
x = seq(0,1,0.005)
param=list(scale=0.6/3,nugget=0,smooth=0.5,mean=2)
cc= GeoCorrFct(x=x, corrmodel="Matern", param=param,model="Poisson",covariance=TRUE)
plot(cc,ylab="Cov",lwd=2,
  main="Poisson covariance",type="l")

################################################################
###
### Example 6.  Space time  semivariogram of a Gaussian random field 
### with  separable Matern correlation model
###
###############################################################

## spatial and temporal distances 
h<-seq(0,3,by=0.04)
times<-seq(0,3,by=0.04)

# Correlation Parameters for the space time separable Matern model 
CorrParam("Matern")
NuisParam("Gaussian")
# Matern Parameters 
param=list(sill=1,scale_s=0.6/3,scale_t=0.5,nugget=0,mean=0,smooth_s=1.5,smooth_t=0.5)
cc= GeoCorrFct(x=h,t=times,corrmodel="Matern_Matern", param=param,
        model="Gaussian",variogram=TRUE)
plot(cc,lwd=2,type="l")


################################################################
###
### Example 7. Correlation of a bivariate Gaussian random field 
### with underlying  separable bivariate  Matern correlation model
###
###############################################################
# Define the spatial distances
x = seq(0,1,0.005)
 #Correlation Parameters for the bivariate sep Matern model 
CorrParam("Bi_Matern")
 #Matern Parameters 
param=list(sill_1=1,sill_2=1,smooth_1=0.5,smooth_2=1,smooth_12=0.75,
           scale_1=0.2/3, scale_2=0.2/3, scale_12=0.2/3,
           mean_1=0,mean_2=0,nugget_1=0,nugget_2=0,pcol=-0.2)
cc= GeoCorrFct(x=x, corrmodel="Bi_Matern", param=param,model="Gaussian")
plot(cc,ylab="corr",lwd=2,type="l")
}

