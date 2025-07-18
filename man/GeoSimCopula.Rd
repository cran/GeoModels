\name{GeoSimCopula}
\alias{GeoSimCopula}
\alias{print.GeoSimCopula}
\encoding{UTF-8}
\title{Simulation of Gaussian and non Gaussian Random Fields using copula.}
\description{
  Simulation of Gaussian and some non Gaussian spatial, spatio-temporal and spatial bivariate 
 random fields using Gaussian or Clayton copula.
  The function return a realization of a Random Field for a given covariance model and covariance parameters.
  Simulation is based on Cholesky decomposition.
}
\usage{
GeoSimCopula(coordx, coordy=NULL,coordz=NULL, coordt=NULL, 
coordx_dyn=NULL, corrmodel, distance="Eucl", grid=FALSE, 
      method="cholesky", model='Gaussian', n=1, param,
      anisopars=NULL,radius=1, sparse=FALSE,
      copula="Gaussian",seed=NULL, X=NULL,spobj=NULL,nrep=1)
}
\arguments{
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
  \item{corrmodel}{String; the name of a correlation model, for the
    description see the Section \bold{Details}.}
  \item{distance}{String; the name of the spatial distance. The default
    is \code{Eucl}, the euclidean distance. See the Section
    \bold{Details}  of \code{\link{GeoFit}}.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as spatial or spatial-temporal realisations on a set
    of non-equispaced spatial sites (irregular grid).}
     \item{method}{String; the type of matrix decomposition used in the simulation. Default is cholesky. The other possible choices is \code{svd}.} 
  \item{model}{String; the type of RF and therefore the densities associated to the likelihood
    objects. \code{Gaussian} is the default, see the Section
    \bold{Details}.}
  \item{n}{Numeric; the number of trials for binomial RFs. The number of successes in the negative Binomial RFs. Default is \eqn{1}.} 
  \item{param}{A list of parameter values required in the
    simulation procedure of RFs, see \bold{Examples}.}
      \item{anisopars}{ A list of two elements "angle" and "ratio" i.e. the anisotropy angle and the anisotropy ratio, respectively.}
    \item{radius}{Numeric; a value indicating  the radius of the sphere when using the great 
    circle distance. Default value is 1.}  
  \item{sparse}{Logical; if \code{TRUE} then  cholesky decomposition is performed
  using sparse matrices algorithms (spam packake). It should be used
  with compactly supported covariance models.\code{FALSE} is the default.}
  \item{copula}{String; the type of copula. It can be "Clayton" or "Gaussian"}
  \item{seed}{Numeric; an integer used in set.seed function to reproduce the simulation.}
  \item{X}{Numeric; Matrix of space-time covariates.}
    \item{spobj}{An object of class sp or spacetime}
    \item{nrep}{Numeric; Numbers of indipendent replicates.}
}


\value{
  Returns an object of class \code{GeoSimCopula}.
  An object of class \code{GeoSimCopula} is a list containing
  at most the following components:
  
  \item{bivariate}{Logical:\code{TRUE} if the Gaussian RF is bivariate, otherwise  \code{FALSE};}
  \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
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
  \item{param}{A list of the parameters}
  \item{radius}{The radius of the sphere if coordinates are passed in lon/lat format;}
  \item{randseed}{The seed used for the random simulation;}
  \item{spacetime}{\code{TRUE} if spatio-temporal and \code{FALSE} if
  spatial RF;}
  \item{copula}{The type of copula}
}



\references{

  Bevilacqua M., Alvarado E.,  Caamano C. (2024) A flexible Clayton-like spatial copula with application to bounded support data.  \emph{Journal of Multivariate Analysis}   \bold{201}


}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}


\examples{
library(GeoModels)

################################################################
###
### Example:  Simulation of a reparametrized Beta RF
### for beta regression
### with Gaussian and Clayton Copula 
### with  underlying Wendland correlation.
###
###############################################################
set.seed(261)
NN=1400
x <- runif(NN);y <- runif(NN)
coords=cbind(x,y)

corrmodel="GenWend"
X=cbind(rep(1,NN),runif(NN))

NuisParam("Beta2",num_betas=2,copula="Gaussian")
CorrParam("GenWend")
#### Gaussian copula
param=list(smooth=0,power2=4, min=0,max=1,
             mean=0.1,mean1=0.1,scale=0.3,nugget=0,shape=5)

data <- GeoSimCopula(coordx=coords, corrmodel=corrmodel, model="Beta2",param=param,
  copula="Gaussian",sparse=TRUE,X=X)$data

quilt.plot(coords,data)


#### Clayton copula
NuisParam("Beta2",num_betas=2,copula="Clayton")
CorrParam("GenWend")
param=list(smooth=0,power2=4, min=0,max=1,
             mean=0.2,mean1=0.1,scale=0.3,nugget=0,shape=6,nu=4)
data1 <- GeoSimCopula(coordx=coords, corrmodel=corrmodel, model="Beta2",param=param,
  copula="Clayton",sparse=TRUE,X=X)$data

hist(data1,freq=FALSE)
quilt.plot(coords,data1)

}
\keyword{Simulation}
