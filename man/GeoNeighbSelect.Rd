\name{GeoNeighbSelect}
\alias{GeoNeighbSelect}
\encoding{UTF-8}
\title{Neighborhood Selection for Pairwise Composite Likelihood Estimation}
\description{
The procedure fits pairwise composite likelihood models over user-specified spatial
or space-time neighborhoods and selects the candidate minimizing the sum of squared
differences between the fitted and empirical semivariograms. A
\code{GeoVariogram} object must be supplied.
}
\usage{
GeoNeighbSelect(data, coordx, coordy=NULL,coordz=NULL, coordt=NULL, coordx_dyn=NULL,
 copula=NULL,corrmodel=NULL, distance="Eucl",fixed=NULL,anisopars=NULL,
 est.aniso=c(FALSE,FALSE), grid=FALSE, likelihood='Marginal',lower=NULL,
 neighb=c(1,2,3,4,5),p_neighb=1,maxtime=Inf, memdist=TRUE,model='Gaussian',
 n=1, ncores=6,optimizer='Nelder-Mead', parallel=FALSE,
 bivariate=FALSE,radius=1, start=NULL,type='Pairwise', upper=NULL,
 weighted=FALSE,X=NULL,spobj=NULL,spdata=NULL,vario=NULL,progress=TRUE,
 check.duplicates=FALSE)
}
\arguments{
 \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation) or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on
 regular grid)
 or a
 (\eqn{t \times d}{t x d})-matrix (a single spatio-temporal realisation) or an (\eqn{d \times d \times t \times n }{d x d x t})-array
 (a single spatio-temporal realisation on regular grid).
 See Details for the accepted data layouts.}
 \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix or (\eqn{d \times 3}{d x 3})-matrix
 Coordinates on a sphere for a fixed radius \code{radius}
 are passed in lon/lat format expressed in decimal degrees.}
 \item{coordy}{A numeric vector giving 1-dimension of
 spatial coordinates; Optional argument, the default is \code{NULL}.}
 \item{coordz}{A numeric vector giving 1-dimension of
 spatial coordinates; Optional argument, the default is \code{NULL}.}

 \item{coordt}{A numeric vector assigning one dimension of temporal coordinates. Optional argument; the default is \code{NULL}, in which case a spatial random field is expected. Temporal coordinates may be irregularly spaced; temporal lags are computed from the supplied coordinate values.}
 \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d x 2})-matrices
 containing dynamical (in time) spatial coordinates. Optional argument, the default is \code{NULL}}
 \item{copula}{String; the type of copula. It can be "Clayton" or "Gaussian"}
 \item{corrmodel}{String; the name of a correlation model, for the
 see \code{\link{GeoCovmatrix}} for the list of implemented correlation models.}
 \item{distance}{String; the name of the spatial distance. The default is \code{Eucl}, the Euclidean distance.
 See Details for the accepted options.}
 \item{fixed}{An optional named list giving the values of the parameters that
 will be considered as known values. The listed parameters for a
 given correlation function will be not estimated.}
 \item{anisopars}{ A list of two elements: "angle" and "ratio" i.e. the anisotropy angle and the anisotropy ratio, respectively.}
 \item{est.aniso}{ A bivariate logical vector providing which anisotropic parameters must be estimated.}


 \item{grid}{Logical; if \code{FALSE} (the default) the data
 are interpreted as spatial or spatio-temporal realisations on a set
 of non-equispaced spatial sites (irregular grid).}

 \item{likelihood}{String; the configuration of the composite
 likelihood. \code{Marginal} is the default (see \code{\link{GeoFit}} for details).}

 \item{lower}{An optional named list giving the values for the lower bound of the space parameter
 when the optimizer is \code{L-BFGS-B} or \code{nlminb} or \code{bobyqa} or \code{optimize}. The
 names of the list must be the same of the names in the \code{start} list.}


 \item{neighb}{Numeric; a vector of positive integers indicating the
 order of neighborhood in the weight function of composite likelihood
 (see \code{\link{GeoFit}} for details).}

 \item{p_neighb}{Numeric; a value in \eqn{(0,1]} specifying the expected fraction
of nearest-neighbor pairs retained through stochastic thinning. If equal to 1
(the default), no thinning is applied. Values below 1 use independent Bernoulli
thinning, so the selected neighborhood depends on the thinning realization. Use a
fixed random seed, or \code{p_neighb=1}, when a stable comparison is required.}

 \item{maxtime}{Numeric vector of non-negative maximum temporal-distance thresholds, expressed in the same units as \code{coordt}, considered when selecting the neighborhood configuration.}


 \item{memdist}{Deprecated logical argument retained for backward compatibility. The selected pair structure is always precomputed and reused. Supplying \code{FALSE} produces a warning and is treated as \code{TRUE}.}

 \item{model}{String; the type of random fields and therefore the densities associated to the likelihood
 objects. \code{Gaussian} is the default, see \code{\link{GeoFit}} for details.}

 \item{n}{Numeric; number of trials in a binomial random fields; number of successes in a negative binomial random fields }
 \item{ncores}{Positive integer or \code{NULL}; default \code{6}. With \code{parallel=TRUE}, an explicit integer requests that many workers, capped by detected cores and the number of candidate fits. Set \code{ncores=NULL} for automatic selection, capped at six workers and normally leaving one detected core free.}

 \item{optimizer}{String; the optimization algorithm
 (see \code{\link{optim}} for details). \code{Nelder-Mead} is the
 default. Other possible choices are \code{nlm}, \code{BFGS}, \code{SANN}, \code{L-BFGS-B} and \code{nlminb} and
 \code{bobyqa}.
 In these last three cases upper and lower bounds can be passed by the user. In the case
 of one-dimensional optimization, the function \code{optimize} is used.}

 \item{parallel}{Logical; default \code{FALSE}. Set \code{TRUE} to evaluate candidate fits in parallel using a \code{future} plan.}
 \item{bivariate}{Logical; if \code{TRUE} the bivariate case is considered.}

 \item{radius}{Numeric; the radius of the sphere in the case of lon-lat coordinates.
 Default value is 1.}

 \item{start}{An optional named list with the initial values of the
 parameters that are used by the numerical routines in maximization
 procedure. \code{NULL} is the default (see \code{\link{GeoFit}} for details).}

 \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
 default) then the marginal composite likelihood is formed by
 pairwise marginal likelihoods (see \code{\link{GeoFit}} for details).}

 \item{upper}{An optional named list giving the values for the upper bound
 of the space parameter when the optimizer is or \code{L-BFGS-B} or \code{bobyqa} or \code{nlminb} or \code{optimize}.
 The names of the list must be the same of the names in the \code{start} list.}

 \item{weighted}{Logical; if \code{TRUE} the likelihood objects are
 weighted (see \code{\link{GeoFit}} for details). If \code{FALSE} (the
 default) the composite likelihood is not weighted.}

 \item{X}{Numeric; Matrix of spatio(temporal)covariates in the linear mean specification.}

 \item{spobj}{An object of class sp or spacetime}
 \item{spdata}{Character:The name of data in the sp or spacetime object}
 \item{vario}{An object of class \code{GeoVariogram} obtained using
\code{\link{GeoVariogram}}. Its spatial, bivariate, or space-time structure must
match the selected correlation model.}
 \item{progress}{Logic; If TRUE then a progress bar is shown.}
\item{check.duplicates}{Logical. If \code{TRUE}, perform a fast exact scan for duplicated observation locations at this user entry point. The default is \code{FALSE}, so no duplicate-location scan is imposed. Internal bootstrap, cross-validation, and refitting calls do not repeat the scan.}
}
\details{
For a spatial univariate model, the criterion compares the fitted and empirical
spatial semivariograms at the empirical bin centers. For a bivariate model, the two
marginal semivariograms and the cross-semivariogram are compared in the order
\eqn{(11,12,22)}. For a space-time model, the spatial and temporal margins and
the complete rectangular space-time semivariogram surface are compared jointly.
Only finite empirical cells contribute to the criterion; a candidate is rejected
when its corresponding fitted values are non-finite or dimensionally incompatible.

The absolute minimum defines the best candidate. The function also reports a more
parsimonious neighborhood using the existing plateau rule: the smallest neighborhood
whose criterion is within 15 percent of the minimum. If every candidate fit fails or
produces a non-finite criterion, the function stops with an error instead of selecting
an arbitrary candidate.
}
\value{
A list with components:
\itemize{
  \item \code{best_neighb}: neighborhood associated with the smallest criterion;
  \item \code{best_maxtime}: selected temporal neighborhood for a space-time model,
  otherwise \code{NULL};
  \item \code{best_estimates}: named parameter estimates for the best candidate;
  \item \code{res}: variogram sum-of-squares criterion for every candidate;
  \item \code{estimates}: matrix of parameter estimates, with columns aligned to
  the names in \code{start};
  \item \code{sugg_neighb}: parsimonious neighborhood selected by the 15 percent
  plateau rule;
  \item \code{sugg_time}: corresponding temporal neighborhood for a space-time
  model, otherwise \code{NULL}.
}
}
\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{

library(GeoModels)


######### spatial case
set.seed(32)
N=150 # number of location sites
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)
mean <- 0.2
# Set the covariance model's parameters:
corrmodel <- "Matern"
sill <- 1;nugget <- 0
scale <- 0.2/3;smooth=0.5

model="Gaussian"
param<-list(mean=mean,sill=sill,nugget=nugget,scale=scale,smooth=smooth)
# Simulation
data <- GeoSim(coordx=coords,corrmodel=corrmodel, param=param,model=model)$data
I=Inf
fixed<-list(nugget=nugget)
start<-list(mean=mean,scale=scale,smooth=smooth,sill=sill)
lower<-list(mean=-I,scale=0,sill=0,smooth=0)
upper<-list(mean=I,scale=I,sill=I,smooth=I)

vario = GeoVariogram(coordx=coords,data=data,maxdist=0.3,numbins=15)

neighb=c(1,2,3,4) ## trying different neighbs
selK <- GeoNeighbSelect(vario=vario,data=data,coordx=coords,corrmodel=corrmodel,
 model=model,neighb=neighb,
 likelihood="Conditional",type="Pairwise",parallel=FALSE,
 optimizer="nlminb",lower=lower,upper=upper,
 start=start,fixed=fixed)
print(selK$best_neighb) ## selected neighbor




}
\keyword{Composite}
