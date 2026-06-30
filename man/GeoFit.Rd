\name{GeoFit}
\alias{GeoFit}
\alias{print.GeoFit}
\encoding{UTF-8}
\title{Maximum-Likelihood-Based Fitting of Gaussian and non-Gaussian random fields.}
\description{
 Maximum weighted composite-likelihood fitting for Gaussian and some non-Gaussian
 univariate spatial, spatio-temporal and bivariate spatial random fields.
 The function allows fixing any of the parameters and setting upper/lower bounds in the optimization.
 Different optimization methods can be used.
}
\usage{
GeoFit(data, coordx, coordy=NULL, coordz=NULL, coordt=NULL,
 coordx_dyn=NULL, copula=NULL, corrmodel=NULL, distance="Eucl",
 fixed=NULL, anisopars=NULL, est.aniso=c(FALSE,FALSE),
 grid=FALSE, likelihood="Marginal", lower=NULL, maxdist=Inf,
 neighb=NULL, p_neighb=1, maxtime=Inf, memdist=TRUE,
 method="cholesky", model="Gaussian", n=1, onlyvar=FALSE,
 optimizer="Nelder-Mead", radius=1, score=FALSE,
 sensitivity=FALSE, sparse=FALSE, start=NULL,
 thin_method="bernoulli",type="Pairwise", upper=NULL, 
 varest=FALSE, weighted=FALSE,
 X=NULL, spobj=NULL, spdata=NULL)
}

\arguments{
 \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation) or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid)
 or a (\eqn{t \times d}{t x d})-matrix (a single spatio-temporal realisation) or a (\eqn{d \times d \times t \times n}{d x d x t x n})-array
 (a single spatio-temporal realisation on regular grid).
 See Details for the accepted data layouts.}

 \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix or (\eqn{d \times 3}{d x 3})-matrix.
 Coordinates on a sphere for a fixed radius \code{radius} are passed in lon/lat format expressed in decimal degrees.}

 \item{coordy}{A numeric vector giving 1-dimension of spatial coordinates; optional argument, default is \code{NULL}.}

 \item{coordz}{A numeric vector giving 1-dimension of spatial coordinates; optional argument, default is \code{NULL}.}

 \item{coordt}{A numeric vector assigning 1-dimension of temporal coordinates. Optional argument, default is \code{NULL}:
 if \code{NULL} then a spatial random field is expected.}

 \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d_t x 2})-matrices containing dynamical (in time) spatial coordinates.
 Optional argument, default is \code{NULL}.}

 \item{copula}{String; the type of copula. It can be \code{"Clayton"} or \code{"Gaussian"}.}

 \item{corrmodel}{String; the name of a correlation model; see \code{\link{GeoCovmatrix}} for the list of implemented correlation models.}

 \item{distance}{String; the name of the spatial distance. Default is \code{"Eucl"} (Euclidean distance). See Details for the accepted options.}

 \item{fixed}{An optional named list giving the values of the parameters that will be considered as known values.
 The listed parameters for a given model/correlation function will not be estimated.}

 \item{anisopars}{A list of two elements: \code{"angle"} and \code{"ratio"}, i.e. the anisotropy angle and the anisotropy ratio, respectively.}

 \item{est.aniso}{A bivariate logical vector providing which anisotropy parameters must be estimated.}

 \item{grid}{Logical; if \code{FALSE} (default) the data are interpreted as spatial or spatio-temporal realisations on a set of non-equispaced spatial sites (irregular grid).}

 \item{likelihood}{String; the configuration of the composite likelihood. \code{"Marginal"} is the default; see Details for the accepted options.}

 \item{lower}{An optional named list giving lower bounds for parameters when the optimizer is \code{L-BFGS-B}, \code{nlminb}, \code{bobyqa} or \code{optimize}.
 Names must match those in \code{start}.}

 \item{maxdist}{Numeric; an optional positive value indicating the maximum spatial distance considered in the composite computation.
 See Details for more information.}

 \item{neighb}{Numeric; an optional positive integer indicating the order of neighborhood in the composite likelihood computation.
 See Details for more information.}

 \item{p_neighb}{Numeric; a value in \eqn{(0,1]} specifying the expected fraction of nearest-neighbor pairs retained through stochastic thinning.
 If \code{1} (default), no thinning is applied and all nearest-neighbor pairs are used. If \code{<1}, pairs are randomly retained (Bernoulli sampling).}

 \item{maxtime}{Numeric; an optional positive integer indicating the order of temporal neighborhood in the composite likelihood computation.}

 \item{memdist}{Logical; if \code{TRUE} then all distances useful in the composite likelihood estimation are computed before the optimization.
 \code{FALSE} is deprecated.}

 \item{method}{String; the type of matrix decomposition/linear algebra backend used in likelihood computations.
 Default is \code{"cholesky"}. Another possible choice is \code{"svd"} (when available).}

 \item{model}{String; the type of random field (and associated density) used in the likelihood objects. Default is \code{"Gaussian"}; see Details for the accepted options.}

 \item{n}{Numeric; number of trials in a binomial random field; number of successes in a negative binomial random field.}

 \item{onlyvar}{Logical; if \code{TRUE} (and \code{varest=TRUE}) only the variance-covariance matrix is computed without optimizing.
 Default is \code{FALSE}.}

 \item{optimizer}{String; the optimization algorithm (see \code{\link{optim}} for details). Default is \code{"Nelder-Mead"}.
 Other possible choices are \code{"nlm"}, \code{"BFGS"}, \code{"SANN"}, \code{"L-BFGS-B"}, \code{"nlminb"}, \code{"bobyqa"}.
 For \code{"L-BFGS-B"}, \code{"nlminb"} and \code{"bobyqa"} bounds can be passed via \code{lower} and \code{upper}.
 In the one-dimensional case, \code{optimize} is used.}

 \item{radius}{Numeric; the radius of the sphere in the case of lon-lat coordinates. Default is \code{1}.}

 \item{score}{Logical; if \code{TRUE} the score function is computed. Default is \code{FALSE}.}

 \item{sensitivity}{Logical; if \code{TRUE} the sensitivity matrix is computed. Default is \code{FALSE}.}

 \item{sparse}{Logical; if \code{TRUE} then maximum likelihood / composite likelihood may exploit sparse-matrix algorithms (e.g., \pkg{spam}).
 Typically used with compactly supported covariance models. Default is \code{FALSE}.}

 \item{start}{An optional named list with initial values for parameters to be estimated. Default is \code{NULL} (see Details).}

 \item{thin_method}{String; thinning scheme used when \code{p_neighb < 1}. Default is \code{"bernoulli"} (independent Bernoulli thinning). An alternative thinning scheme is 
 \code{"FixedBudget"}}. 

 \item{type}{String; the type of likelihood objects. If \code{"Pairwise"} (default) the composite likelihood is formed by pairwise components (see Details).}

 \item{upper}{An optional named list giving upper bounds for parameters when the optimizer is \code{L-BFGS-B}, \code{nlminb}, \code{bobyqa} or \code{optimize}.
 Names must match those in \code{start}.}

 \item{varest}{Logical; if \code{TRUE} the estimates variances and standard errors are returned. For composite likelihood estimation it is deprecated.
 Use \code{sensitivity=TRUE} and update the object using \code{GeoVarest}. Default is \code{FALSE}.}

 \item{weighted}{Logical; if \code{TRUE} the likelihood objects are weighted; see Details for the accepted options. Default is \code{FALSE}.}

 \item{X}{Numeric; matrix of spatio(temporal) covariates in the linear mean specification.}

 \item{spobj}{An object of class \code{sp} or \code{spacetime}.}

 \item{spdata}{Character; the name of the data component in the \code{sp} or \code{spacetime} object.}
}
\details{
 \code{GeoFit} provides weighted composite likelihood estimation based on pairs and independence composite likelihood estimation
 for Gaussian and non-Gaussian random fields. Specifically, marginal and conditional pairwise likelihoods are available
 for each type of random field.

 The optimization method is specified using \code{optimizer}. The default method is \code{Nelder-Mead}; other available methods are
 \code{nlm}, \code{BFGS}, \code{SANN}, \code{L-BFGS-B}, \code{bobyqa}, and \code{nlminb}. In the last three cases,
 bounds can be specified using \code{lower} and \code{upper}.

 Depending on the dimension of \code{data} and on the name of the correlation model,
 the observations are assumed to be a realization of a spatial, spatio-temporal or bivariate random field.
 Specifically, with \code{data}, \code{coordx}, \code{coordy}, \code{coordt}:
 \itemize{
 \item If \code{data} is a numeric \eqn{d}{d}-dimensional vector and \code{coordx}, \code{coordy} are two numeric \eqn{d}{d}-dimensional vectors
 (or \code{coordx} is a (\eqn{d \times 2}{d x 2})-matrix and \code{coordy=NULL}), then the data are interpreted as a single spatial realisation observed on \eqn{d}{d} spatial sites;
 \item If \code{data} is a numeric (\eqn{t \times d}{t x d})-matrix and \code{coordt} is a numeric \eqn{t}{t}-dimensional vector, then the data are interpreted as a single spatio-temporal realisation observed on \code{d} sites and \code{t} times;
 \item If \code{data} is a numeric (\eqn{2 \times d}{2 x d})-matrix, then the data are interpreted as a single bivariate spatial realisation observed on \code{d} spatial sites;
 \item If \code{data} is a list, \code{coordx_dyn} is a list and \code{coordt} is a numeric \eqn{t}{t}-dimensional vector, then the data are interpreted as a spatio-temporal realisation observed on dynamical spatial sites (different locations for each time) and for \code{t} times.
 }

 It is also possible to specify a matrix of covariates using \code{X}. Specifically:
 \itemize{
 \item In the spatial case, \code{X} must be a (\eqn{d \times k}{d x k}) matrix associated to \code{data} (a \eqn{d}{d}-vector);
 \item In the spatio-temporal case, \code{X} must be a (\eqn{N \times k}{N x k}) matrix associated to \code{data} (a \eqn{t \times d}{t x d}-matrix), where \eqn{N=t\times d};
 \item In the bivariate case, \code{X} must be a (\eqn{N \times k}{N x k}) matrix associated to \code{data} (a \eqn{2 \times d}{2 x d}-matrix), where \eqn{N=2\times d}.
 }

 The \code{distance} parameter allows different kinds of spatial distances:
 \enumerate{
 \item \code{Eucl}, Euclidean distance (default);
 \item \code{Chor}, chordal distance;
 \item \code{Geod}, geodesic distance.
 }

 The \code{likelihood} parameter represents the composite-likelihood configuration:
 \enumerate{
 \item \code{Conditional}, composite likelihood formed by conditionals;
 \item \code{Marginal}, composite likelihood formed by marginals (default);
 \item \code{Full}, standard likelihood.
 }

 It must be coupled with \code{type}:
 \enumerate{
 \item \code{Pairwise}, composite likelihood based on pairs;
 \item \code{Independence}, composite likelihood based on independence;
 \item \code{Standard}, standard likelihood.
 }

 Stochastic thinning of nearest-neighbor pairs can be enabled via \code{p_neighb<1}. The argument \code{thin_method}
 controls the thinning scheme (default \code{"bernoulli"}).
}
\value{
 Returns an object of class \code{GeoFit}.
 An object of class \code{GeoFit} is a list containing at most the following components:
 \item{bivariate}{Logical: \code{TRUE} if the random field is bivariate, otherwise \code{FALSE}.}
 \item{clic}{The composite information criterion after a \code{GeoVarest} call; if the full likelihood is considered then it coincides with AIC.}
 \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates.}
 \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates.}
 \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates.}
 \item{coordx_dyn}{A list of dynamical (in time) spatial coordinates.}
 \item{conf.int}{Confidence intervals for standard maximum likelihood estimation.}
 \item{convergence}{A string that denotes if convergence is reached.}
 \item{copula}{The type of copula.}
 \item{corrmodel}{The correlation model.}
 \item{data}{The vector/matrix/array (or list) of data.}
 \item{distance}{The type of spatial distance.}
 \item{fixed}{A list of fixed parameters.}
 \item{iterations}{The number of iterations used by the numerical routine.}
 \item{likelihood}{The configuration of the composite likelihood.}
 \item{logCompLik}{The value of the log composite-likelihood at the maximum.}
 \item{maxdist}{The maximum spatial distance used in the weighted composite likelihood (or \code{NULL}).}
 \item{maxtime}{The order of temporal neighborhood in the composite likelihood computation.}
 \item{message}{Extra message passed from the numerical routines.}
 \item{model}{The density associated to the likelihood objects.}
 \item{missp}{\code{TRUE} if a misspecified Gaussian model is used in the composite likelihood.}
 \item{n}{The number of trials/successes for (negative) binomial models.}
 \item{neighb}{The order of spatial neighborhood in the composite likelihood computation.}
 \item{ns}{The number of (different) location sites in the bivariate case.}
 \item{numcoord}{The number of spatial coordinates.}
 \item{numtime}{The number of temporal realisations.}
 \item{param}{A list of parameter estimates.}
 \item{radius}{The radius of the sphere in the case of great-circle distance.}
 \item{stderr}{Standard errors for standard maximum likelihood estimation.}
 \item{sensmat}{The sensitivity matrix.}
 \item{varcov}{The variance-covariance matrix of the estimates.}
 \item{type}{The type of likelihood objects.}
 \item{X}{The matrix of covariates.}
}

\references{


 \bold{General Composite-likelihood}:

 Varin, C., Reid, N. and Firth, D. (2011). An Overview of Composite
 Likelihood Methods. \emph{Statistica Sinica}, \bold{21}, 5--42.

 Varin, C. and Vidoni, P. (2005)
 A Note on Composite Likelihood Inference and Model Selection.
 \emph{Biometrika}, \bold{92}, 519--528.


 \bold{Non Gaussian random fields}: 

 Alegrıa A., Caro S., Bevilacqua M., Porcu E., Clarke J. (2017) \emph{Estimating covariance functions of multivariate skew-Gaussian random fields on the sphere}. Spatial Statistics \bold{22} 388-402

 Alegria A., Bevilacqua, M., Porcu, E. (2016) 
 Likelihood-based inference for multivariate space-time wrapped-Gaussian fields. 
 \emph{Journal of Statistical Computation and Simulation}. \bold{86(13)}, 2583--2597.

 Bevilacqua M., Caamano C., Gaetan C. (2020) On modeling positive continuous data with spatio-temporal dependence. \emph{Environmetrics} \bold{31(7)}

 Bevilacqua M., Caamaño C., Arellano Valle R.B., Morales-Oñate V. (2020) non-Gaussian Geostatistical Modeling using (skew) t Processes. \emph{Scandinavian Journal of Statistics}.

 Blasi F., Caamaño C., Bevilacqua M., Furrer R. (2022) A selective view of climatological data and likelihood estimation \emph{Spatial Statistics} 10.1016/j.spasta.2022.100596

 Bevilacqua M., Caamaño C., Arellano-Valle R. B., Camilo Gomez C. (2022) A class of random fields with two-piece marginal distributions for modeling point-referenced data with spatial outliers. \emph{Test} 10.1007/s11749-021-00797-5
 
 Morales-Navarrete D., Bevilacqua M., Caamaño C., Castro L.M. (2022) Modelling Point Referenced Spatial Count Data: A Poisson Process Approach \emph{TJournal of the American Statistical Association} To appear

 Caamaño C., Bevilacqua M., López C., Morales-Oñate V. (2023) Nearest neighbours weighted composite likelihood based on pairs for (non-)Gaussian massive spatial data with an application to Tukey-hh random fields estimation \emph{Computational Statistics and Data Analysis} To appear

 Bevilacqua M., Alvarado E., Caamaño C., (2023) A flexible Clayton-like spatial copula with application to bounded support data \emph{Journal of Multivariate Analysis} To appear 

 \bold{Weighted Composite-likelihood for (non-)Gaussian random fields}:

 Bevilacqua, M. Gaetan, C., Mateu, J. and Porcu, E. (2012)
 Estimating space and space-time covariance functions for large data
 sets: a weighted composite likelihood approach.
 \emph{Journal of the American Statistical Association, Theory &
 Methods}, \bold{107}, 268--280.

 Bevilacqua, M., Gaetan, C. (2015) Comparing composite likelihood methods based on pairs for spatial Gaussian random fields. \emph{Statistics and Computing}, \bold{25(5)}, 877-892.

 Caamaño C., Bevilacqua M., López C., Morales-Oñate V. (2023) Nearest neighbours weighted composite likelihood based on pairs for (non-)Gaussian massive spatial data with an application to Tukey-hh random fields estimation \emph{Computational Statistics and Data Analysis} To appear

 
}



\seealso{
\code{\link{GeoCovmatrix}} for covariance matrix construction,
  \code{\link{GeoSim}} for simulation,
  \code{\link{GeoKrig}} for prediction,
  \code{\link{GeoVariogram}} and \code{\link{GeoWLS}} for variogram-based tools,
  \code{\link{GeoVarest}} for bootstrap variance estimation.
}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}
 
\examples{
library(GeoModels)

###############################################################
############ Examples of spatial Gaussian random fieldss ################
###############################################################


# Define the spatial-coordinates of the points:
set.seed(3)
N=300 # number of location sites
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)

# Define spatial matrix covariates and regression parameters
X=cbind(rep(1,N),runif(N))
mean <- 0.2
mean1 <- -0.5

# Set the covariance model's parameters:
corrmodel <- "Matern"
sill <- 1
nugget <- 0
scale <- 0.2/3
smooth=0.5


param<-list(mean=mean,mean1=mean1,sill=sill,nugget=nugget,scale=scale,smooth=smooth)

# Simulation of the spatial Gaussian random fields:
data <- GeoSim(coordx=coords,corrmodel=corrmodel, param=param,X=X)$data



################################################################
###
### Example 0. Maximum independence composite likelihood fitting of
### a Gaussian random fields (no dependence parameters)
### 
###############################################################
# setting starting parameters to be estimated
start<-list(mean=mean,mean1=mean1,sill=sill)

fit1 <- GeoFit(data=data,coordx=coords,likelihood="Marginal",
 type="Independence", start=start,X=X)
print(fit1)


################################################################
###
### Example 1. Maximum conditional pairwise likelihood fitting of
### a Gaussian random fields using Nelder-Mead
### 
###############################################################
# setting fixed and starting parameters to be estimated
fixed<-list(nugget=nugget,smooth=smooth)
start<-list(mean=mean,mean1=mean1,scale=scale,sill=sill)

fit1 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, 
 neighb=3,likelihood="Conditional",optimizer="Nelder-Mead",
 type="Pairwise", start=start,fixed=fixed,X=X)
print(fit1)

################################################################
###
### Example 2. Standard Maximum likelihood fitting of
### a Gaussian random fields using nlminb
###
###############################################################
# Define the spatial-coordinates of the points:
set.seed(3)
N=250 # number of location sites
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)

param<-list(mean=mean,sill=sill,nugget=nugget,scale=scale,smooth=smooth)

data <- GeoSim(coordx=coords,corrmodel=corrmodel, param=param)$data

# setting fixed and parameters to be estimated
fixed<-list(nugget=nugget,smooth=smooth)
start<-list(mean=mean,scale=scale,sill=sill)

I=Inf
lower<-list(mean=-I,scale=0,sill=0)
upper<-list(mean=I,scale=I,sill=I)
fit2 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel,
 optimizer="nlminb",upper=upper,lower=lower,
 likelihood="Full",type="Standard", 
 start=start,fixed=fixed)
print(fit2)


###############################################################
############ Examples of spatial non-Gaussian random fieldss #############
###############################################################


################################################################
###
### Example 3. Maximum pairwise likelihood fitting of a Weibull random fields 
### with Generalized Wendland correlation with Nelder-Mead
### 
###############################################################
set.seed(524)
# Define the spatial-coordinates of the points:
N=300
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)
X=cbind(rep(1,N),runif(N))
mean=1; mean1=2 # regression parameters
nugget=0
shape=2
scale=0.2
smooth=0

model="Weibull"
corrmodel="GenWend"
param=list(mean=mean,mean1=mean1,scale=scale,
 shape=shape,nugget=nugget,power2=4,smooth=smooth)
# Simulation of a non stationary weibull random fields:
data <- GeoSim(coordx=coords, corrmodel=corrmodel,model=model,X=X,
 param=param)$data


fixed<-list(nugget=nugget,power2=4,smooth=smooth)
start<-list(mean=mean,mean1=mean1,scale=scale,shape=shape)

# Maximum independence likelihood:
fit <- GeoFit(data=data, coordx=coords, X=X, 
 likelihood="Marginal",type="Independence", corrmodel=corrmodel,
 ,model=model, start=start, fixed=fixed)
print(unlist(fit$param))

## estimating dependence parameter fixing vector mean parameter
Xb=as.numeric(X \%*\% c(mean,mean1))
fixed<-list(nugget=nugget,power2=4,smooth=smooth,mean=Xb)
start<-list(scale=scale,shape=shape)

# Maximum conditional composite-likelihood fitting of the random fields:
fit1 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
 neighb=3,likelihood="Conditional",type="Pairwise",
 optimizer="Nelder-Mead",
 start=start,fixed=fixed)
print(unlist(fit1$param))



### joint estimation of the dependence parameter and mean parameters
fixed<-list(nugget=nugget,power2=4,smooth=smooth)
start<-list(mean=mean,mean1=mean1,scale=scale,shape=shape)
fit2 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
 neighb=3,likelihood="Conditional",type="Pairwise",X=X,
 optimizer="Nelder-Mead",
 start=start,fixed=fixed)
print(unlist(fit2$param))



################################################################
###
### Example 4. Maximum pairwise likelihood fitting of
### a Skew-Gaussian spatial random fields with Wendland correlation
###
###############################################################
set.seed(261)
model="SkewGaussian"
# Define the spatial-coordinates of the points:
x <- runif(500, 0, 1)
y <- runif(500, 0, 1)
coords <- cbind(x,y)

corrmodel="Wend0"
mean=0;nugget=0
sill=1
skew=-4.5
power2=4
c_supp=0.2

# model parameters
param=list(power2=power2,skew=skew,
 mean=mean,sill=sill,scale=c_supp,nugget=nugget)
data <- GeoSim(coordx=coords, corrmodel=corrmodel,model=model, param=param)$data

plot(density(data))
fixed=list(power2=power2,nugget=nugget)
start=list(scale=c_supp,skew=skew,mean=mean,sill=sill)
lower=list(scale=0,skew=-I,mean=-I,sill=0)
upper=list(scale=I,skew=I,mean=I,sill=I)
# Maximum marginal pairwise likelihood:
fit1 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
 neighb=3,likelihood="Marginal",type="Pairwise",
 optimizer="bobyqa",lower=lower,upper=upper,
 start=start,fixed=fixed)
print(unlist(fit1$param))


################################################################
###
### Example 5. Maximum pairwise likelihood fitting of 
### a Bernoulli random fields with exponential correlation 
###
###############################################################

set.seed(422)
N=250
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)
mean=0.1; mean1=0.8; mean2=-0.5 # regression parameters
X=cbind(rep(1,N),runif(N),runif(N)) # marix covariates
corrmodel <- "Wend0"
param=list(mean=mean,mean1=mean1,mean2=mean2,nugget=0,scale=0.2,power2=4)
# Simulation of the spatial Binomial-Gaussian random fields:
data <- GeoSim(coordx=coords, corrmodel=corrmodel, model="Binomial", n=1,X=X,
 param=param)$data


## estimating the marginal parameters using independence cl
fixed <- list(power2=4,scale=0.2,nugget=0)
start <- list(mean=mean,mean1=mean1,mean2=mean2)

# Maximum independence likelihood:
fit <- GeoFit(data=data, coordx=coords,n=1, X=X, 
 likelihood="Marginal",type="Independence",corrmodel=corrmodel,
 ,model="Binomial", start=start, fixed=fixed)
 
print(fit)


## estimating dependence parameter fixing vector mean parameter
Xb=as.numeric(X \%*\% c(mean,mean1,mean2))
fixed <- list(nugget=0,power2=4,mean=Xb)
start <- list(scale=0.2)

# Maximum conditional pairwise likelihood:
fit1 <- GeoFit(data=data, coordx=coords, corrmodel=corrmodel,n=1, 
 likelihood="Conditional",type="Pairwise", neighb=3
 ,model="Binomial", start=start, fixed=fixed)
 
print(fit1)


## estimating jointly marginal and dependnce parameters
fixed <- list(nugget=0,power2=4)
start <- list(mean=mean,mean1=mean1,mean2=mean2,scale=0.2)

# Maximum conditional pairwise likelihood:
fit2 <- GeoFit(data=data, coordx=coords, corrmodel=corrmodel,n=1, X=X, 
 likelihood="Conditional",type="Pairwise", neighb=3
 ,model="Binomial", start=start, fixed=fixed)
 
print(fit2)


###############################################################
######### Examples of Gaussian spatio-temporal random fieldss ###########
###############################################################
set.seed(52)
# Define the temporal sequence:
time <- seq(1, 9, 1)

# Define the spatial-coordinates of the points:
x <- runif(20, 0, 1)
y <- runif(20, 0, 1)
coords=cbind(x,y)

# Set the covariance model's parameters:
scale_s=0.2/3;scale_t=1
smooth_s=0.5;smooth_t=0.5
sill=1
nugget=0
mean=0

param<-list(mean=0,scale_s=scale_s,scale_t=scale_t,
 smooth_t=smooth_t, smooth_s=smooth_s ,sill=sill,nugget=nugget)

# Simulation of the spatio-temporal Gaussian random fields:
data <- GeoSim(coordx=coords,coordt=time,corrmodel="Matern_Matern",
 param=param)$data

################################################################
###
### Example 6. Maximum pairwise likelihood fitting of a
### space time Gaussian random fields with double-exponential correlation
###
###############################################################
# Fixed parameters
fixed<-list(nugget=nugget,smooth_s=smooth_s,smooth_t=smooth_t)
# Starting value for the estimated parameters
start<-list(mean=mean,scale_s=scale_s,scale_t=scale_t,sill=sill)

# Maximum composite-likelihood fitting of the random fields:
fit <- GeoFit(data=data,coordx=coords,coordt=time,
 corrmodel="Matern_Matern",maxtime=1,neighb=3,
 likelihood="Marginal",type="Pairwise",
 start=start,fixed=fixed)
print(fit)



###############################################################
######### Examples of a bivariate Gaussian random fields ###########
###############################################################

################################################################
### Example 7. Maximum pairwise likelihood fitting of a
### bivariate Gaussian random fields with separable Bivariate matern 
### (cross) correlation model 
###############################################################

# Define the spatial-coordinates of the points:
set.seed(89)
x <- runif(300, 0, 1)
y <- runif(300, 0, 1)
coords=cbind(x,y)
# parameters
param=list(mean_1=0,mean_2=0,scale=0.1,smooth=0.5,sill_1=1,sill_2=1,
 nugget_1=0,nugget_2=0,pcol=0.2)

# Simulation of a spatial bivariate Gaussian random fields:
data <- GeoSim(coordx=coords, corrmodel="Bi_Matern_sep", 
 param=param)$data

# selecting fixed and estimated parameters
fixed=list(mean_1=0,mean_2=0,nugget_1=0,nugget_2=0,smooth=0.5)
start=list(sill_1=var(data[1,]),sill_2=var(data[2,]),
 scale=0.1,pcol=cor(data[1,],data[2,]))


# Maximum marginal pairwise likelihood
fitcl<- GeoFit(data=data, coordx=coords, corrmodel="Bi_Matern_sep",
 likelihood="Marginal",type="Pairwise",
 start=start,fixed=fixed,
 neighb=3)
print(fitcl)

}
\keyword{Composite}
