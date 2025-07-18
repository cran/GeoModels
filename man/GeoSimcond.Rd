\name{GeoSimcond}
\alias{GeoSimcond}
\title{Conditional (fast) simulation of Gaussian and non Gaussian Random Fields.}
\description{
  Simulates spatial (spatio-temporal) Gaussian and non-Gaussian random fields 
  conditioned on observed data using specified correlation models.
}
\usage{
GeoSimcond(estobj = NULL, data, coordx, coordy = NULL, coordz = NULL, coordt = NULL,
           coordx_dyn = NULL, corrmodel, distance = "Eucl", grid = FALSE, loc,
           maxdist = NULL, maxtime = NULL, method = "Cholesky", model = "Gaussian",
           n = 1, nrep = 1, local = FALSE, L = 1000, neighb = NULL,
           param, anisopars = NULL, radius = 1, sparse = FALSE, time = NULL,
           copula = NULL, X = NULL, Xloc = NULL, Mloc = NULL,
           parallel=TRUE, ncores = NULL)
}
\arguments{
  \item{estobj}{Object of class \code{Geofit} containing model information}
  \item{data}{Numeric vector/matrix/array of observed data}
  \item{coordx}{Numeric matrix of spatial coordinates (d x 2 or d x 3)}
  \item{coordy}{Optional numeric vector of y-coordinates}
  \item{coordz}{Optional numeric vector of z-coordinates}
  \item{coordt}{Optional numeric vector of temporal coordinates}
  \item{coordx_dyn}{Optional list of dynamic spatial coordinates}
  \item{corrmodel}{String specifying correlation model name}
  \item{distance}{String specifying distance metric (default: "Eucl")}
  \item{grid}{Logical for regular grid (default: FALSE)}
  \item{loc}{Numeric matrix of prediction locations (n x 2)}
  \item{maxdist}{Optional maximum distance for local kriging}
  \item{maxtime}{Optional maximum temporal distance}
  \item{method}{String for decomposition method ("Cholesky", "TB", or "CE")}
  \item{model}{String specifying random field type (default: "Gaussian")}
  \item{n}{Number of trials for binomial RFs (default: 1)}
  \item{nrep}{Number of independent replicates (default: 1)}
  \item{local}{Logical for local kriging (default: FALSE)}
  \item{L}{Number of lines for turning bands method (default: 1000)}
  \item{neighb}{Optional neighborhood order for local kriging}
  \item{param}{List of parameter values}
  \item{anisopars}{List with anisotropy angle and ratio}
  \item{radius}{Radius for spherical coordinates (default: Earth's radius)}
  \item{sparse}{Logical for sparse matrix algorithms (default: FALSE)}
  \item{time}{Optional vector of temporal instants to predict}
  \item{copula}{Optional string specifying copula type}
  \item{X}{Optional matrix of spatio-temporal covariates}
  \item{Xloc}{Optional matrix of covariates for prediction locations}
  \item{Mloc}{Optional vector of estimated means for prediction locations}
  \item{parallel}{If TRUE the the computation is parallelized}
  \item{ncores}{Numbers of cores involved in the parallelization}

}
\details{
  For Gaussian RF, performs conditional simulation using three steps:
  \enumerate{
    \item Unconditional simulation at observed and prediction locations
    \item Simple kriging estimates at observed locations
    \item Combination to produce conditional simulations
  }
  For large datasets, approximate methods ("TB" or "CE") are recommended
  coupled with local kriging (local=TRUE and  neighb=k) and using parallelization
  (parallel=T).
}
\value{
  Returns an object of class \code{GeoSimcond} containing:
  \itemize{
    \item Simulated field realizations
    \item Model parameters and settings
    \item Spatial/temporal coordinates
    \item Covariance information
  }
}
\references{
Gaetan, C. and Guyon, X. (2010) \emph{Spatial Statistics and Modelling}.
Springer Verlag, New York.
}
\seealso{
\code{\link[GeoModels]{GeoSim}}, 
\code{\link[GeoModels]{GeoKrig}}
}
\author{
Moreno Bevilacqua \email{moreno.bevilacqua89@gmail.com},\\
Víctor Morales Oñate \email{victor.morales@uv.cl},\\
Christian Caamaño-Carrillo \email{chcaaman@ubiobio.cl}
}
\keyword{spatial}
\keyword{simulation}
\examples{
library(GeoModels)

##############################################
## conditional simulation of a Gaussian rf ###
##############################################
model="Gaussian"
set.seed(79)
### conditioning locations
x = runif(250, 0, 1)
y = runif(250, 0, 1)
coords=cbind(x,y)

# Set the exponential cov parameters:
corrmodel = "GenWend"
mean=0; sill=1; nugget=0
scale=0.2;smooth=0;power2=4

param=list(mean=mean,sill=sill,nugget=nugget,scale=scale,smooth=smooth,power2=power2)

# Simulation 
data = GeoSim(coordx=coords, corrmodel=corrmodel,model=model,
              param=param)$data

## estimation with pairwise likelihood
fixed=list(nugget=nugget,smooth=smooth,power2=power2)
start=list(mean=0,scale=scale,sill=1)
I=Inf
lower=list(mean=-I,scale=0,sill=0)
upper=list(mean= I,scale=I,sill=I)
# Maximum pairwise likelihood fitting :
fit = GeoFit(data, coordx=coords, corrmodel=corrmodel,model=model,
             likelihood='Marginal', type='Pairwise',neighb=3,
             optimizer="nlminb", lower=lower,upper=upper,
             start=start,fixed=fixed)

# locations to simulate 
xx=seq(0,1,0.025)
loc_to_sim=as.matrix(expand.grid(xx,xx))

# Conditional simulation
sim_result <- GeoSimcond(fit,loc = loc_to_sim,nrep=5,parallel=FALSE)$condsim
sim_result <- do.call(rbind, sim_result)

par(mfrow=c(1,2))
quilt.plot(coords, data)
quilt.plot(loc_to_sim, colMeans(sim_result))
par(mfrow=c(1,1))

##############################################
## conditional simulation of a LogGaussian rf 
##############################################
model="LogGaussian"
set.seed(79)
### conditioning locations
x = runif(250, 0, 1)
y = runif(250, 0, 1)
coords=cbind(x,y)

# Set the exponential cov parameters:
corrmodel = "Matern"
mean=0; sill=.1; nugget=0
scale=0.2;smooth=0.5

param=list(mean=mean,sill=sill,nugget=nugget,scale=scale,smooth=smooth)

# Simulation 
data = GeoSim(coordx=coords, corrmodel=corrmodel,model=model,
              param=param)$data

## estimation with pairwise likelihood
fixed=list(nugget=nugget,smooth=smooth)
start=list(mean=0,scale=scale,sill=1)
I=Inf
lower=list(mean=-I,scale=0,sill=0)
upper=list(mean= I,scale=I,sill=I)
# Maximum pairwise likelihood fitting :
fit = GeoFit(data, coordx=coords, corrmodel=corrmodel,model=model,
             likelihood='Marginal', type='Pairwise',neighb=3,
             optimizer="nlminb", lower=lower,upper=upper,
             start=start,fixed=fixed)

# locations to simulate 
xx=seq(0,1,0.025)
loc_to_sim=as.matrix(expand.grid(xx,xx))

# Conditional simulation
sim_result <- GeoSimcond(fit,loc = loc_to_sim,nrep=5,parallel=FALSE)$condsim
sim_result <- do.call(rbind, sim_result)

par(mfrow=c(1,2))
quilt.plot(coords, data)
quilt.plot(loc_to_sim, colMeans(sim_result))
par(mfrow=c(1,1))

}
