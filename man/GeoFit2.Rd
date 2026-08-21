\name{GeoFit2}
\alias{GeoFit2}
\encoding{UTF-8}
\title{Fitting Gaussian and Non-Gaussian Random Fields with Automatic Marginal Initialization}
\description{
 A univariate starting-value wrapper around \code{\link{GeoFit}}. The function optionally
 performs a preliminary fit under spatial independence to estimate available
 marginal parameters, uses those estimates to update the starting values, and
 then calls \code{GeoFit} for the requested final fit. If the preliminary
 independence fit is unavailable or fails, the original starting values are
 retained. Bivariate correlation models are not supported.
}
\usage{
GeoFit2(data, coordx = NULL, coordy = NULL, coordz = NULL, coordt = NULL,
 coordx_dyn = NULL, copula = NULL, corrmodel = NULL,
 distance = "Eucl", fixed = NULL, anisopars = NULL,
 est.aniso = c(FALSE, FALSE), grid = FALSE,
 likelihood = "Marginal", lower = NULL, maxdist = Inf,
 neighb = NULL, p_neighb = 1, maxtime = Inf, memdist = TRUE,
 method = "cholesky", model = "Gaussian", n = 1,
 onlyvar = FALSE, optimizer = "Nelder-Mead", radius = 1,
 score = FALSE, sensitivity = FALSE, sparse = FALSE,
 start = NULL, thin_method = "bernoulli", type = "Pairwise",
 upper = NULL, varest = FALSE, weighted = FALSE, X = NULL,
 spobj = NULL, spdata = NULL, independence_start = TRUE,
 independence_optimizer = "Nelder-Mead",
 warn_independence_failure = TRUE, check.duplicates = FALSE)
}
\arguments{
 \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation) or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid)
 or a (\eqn{t \times d}{t x d})-matrix (a single spatio-temporal realisation) or an (\eqn{d \times d \times t \times n}{d x d x t x n})-array
 (a single spatio-temporal realisation on regular grid).
 See Details for the accepted data layouts.}

 \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix or (\eqn{d \times 3}{d x 3})-matrix.
 Coordinates on a sphere for a fixed radius \code{radius} are passed in lon/lat format expressed in decimal degrees.
 It may be omitted when coordinates are supplied through \code{coordx_dyn} or \code{spobj}.}

 \item{coordy}{A numeric vector giving 1-dimension of spatial coordinates; optional argument, default is \code{NULL}.}

 \item{coordz}{A numeric vector giving 1-dimension of spatial coordinates; optional argument, default is \code{NULL}.}

 \item{coordt}{A numeric vector assigning one dimension of the observation-time coordinates. Optional argument, default is \code{NULL}; if \code{NULL}, a spatial random field is expected. Temporal coordinates may be irregularly spaced; temporal lags are computed from the supplied coordinate values.}

 \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d_t x 2})-matrices containing dynamical (in time) spatial coordinates.
 Optional argument, default is \code{NULL}.}

 \item{copula}{String; the type of copula. It can be \code{"Clayton"} or \code{"Gaussian"}.}

 \item{corrmodel}{String; the name of a correlation model; see \code{\link{GeoCovmatrix}} for the list of implemented correlation models.}

 \item{distance}{String; the name of the spatial distance. Default is \code{"Eucl"} (Euclidean distance). See Details for the accepted options.}

 \item{fixed}{An optional named list giving the values of the parameters that will be considered as known values.
 The listed parameters for a given correlation function will not be estimated.}

 \item{anisopars}{A list of two elements: \code{"angle"} and \code{"ratio"}, i.e. the anisotropy angle and the anisotropy ratio, respectively. Geometric anisotropy requires \code{distance="Eucl"}.}

 \item{est.aniso}{A bivariate logical vector providing which anisotropic parameters must be estimated.}

 \item{grid}{Logical; if \code{FALSE} (default) the data are interpreted as spatial or spatio-temporal realisations on a set of non-equispaced spatial sites. If \code{TRUE}, the grid axes and data dimensions are handled by \code{GeoFit} and converted internally to the canonical explicit-coordinate layout before fitting.}

 \item{likelihood}{String; the configuration of the composite likelihood. \code{"Marginal"} is the default; see Details for the accepted options.}

 \item{lower}{An optional named list giving lower bounds for parameters when the optimizer is \code{L-BFGS-B}, \code{nlminb}, \code{bobyqa} or \code{optimize}.
 Names must match parameters that are being estimated.}

 \item{maxdist}{Numeric; an optional positive value indicating the maximum spatial distance considered in the composite likelihood computation.
 See Details for more information.}

 \item{neighb}{Numeric; an optional positive integer indicating the order of neighborhood in the composite likelihood computation.
 See Details for more information.}

 \item{p_neighb}{Numeric scalar in \eqn{(0,1]}. If \code{1} (default), no thinning is applied. For \code{thin_method="bernoulli"}, it is the candidate-pair inclusion probability in the constant-weight design and controls the retained size in expectation. For \code{thin_method="FixedBudget"}, it defines the exact global retained-pair budget \eqn{K=\mathrm{round}(p_{neighb}d)}, where \eqn{d} is the candidate-pair count.}

 \item{maxtime}{Numeric; an optional non-negative maximum temporal-distance threshold, expressed in the same units as \code{coordt}, used to select pairwise composite-likelihood contributions.}

 \item{memdist}{Deprecated logical argument retained for backward compatibility. The selected pair structure is always precomputed and reused during composite-likelihood optimization. Supplying \code{FALSE} produces a warning and is treated as \code{TRUE}.}

 \item{method}{String; the type of matrix decomposition used in the likelihood computation. Default is \code{"cholesky"}.
 Another possible choice is \code{"svd"} (when available).}

 \item{model}{String; the type of RF and therefore the densities associated to the likelihood objects.
 \code{"Gaussian"} is the default; see Details for the accepted options.}

 \item{n}{Positive integer scalar, or one positive integer per observation, for models using a Binomial/Negative-Binomial trial or success count.}

 \item{onlyvar}{Logical; if \code{TRUE} (and \code{varest=TRUE}) only the variance-covariance matrix is computed without optimizing.
 Default is \code{FALSE}.}

 \item{optimizer}{String; the optimization algorithm (see \code{\link{optim}} for details). \code{"Nelder-Mead"} is the default.
 Other possible choices are \code{"nlm"}, \code{"BFGS"}, \code{"SANN"}, \code{"L-BFGS-B"}, \code{"nlminb"}, \code{"bobyqa"}.
 In these last three cases upper and lower bounds can be passed by the user. In the one-dimensional case, \code{optimize} is used.}

 \item{radius}{Numeric; the radius of the sphere in the case of lon-lat coordinates. Default value is \code{1}.}

 \item{score}{Logical; should score function be computed? Default is \code{FALSE}.}

 \item{sensitivity}{Logical; if \code{TRUE} then the sensitivity matrix is computed.}

 \item{sparse}{Logical; if \code{TRUE} then maximum likelihood is computed using sparse matrix algorithms (e.g., \pkg{spam}).
 It should be used with compactly supported covariance models. Default is \code{FALSE}.}

 \item{start}{An optional named list with initial values for parameters used by the numerical routines in the maximization procedure.
 Default is \code{NULL}; omitted parameters receive internal starting values (see Details).}

 \item{thin_method}{String; thinning scheme used when \code{p_neighb < 1}. Default is \code{"bernoulli"} (independent Bernoulli thinning). An alternative thinning scheme is \code{"FixedBudget"}. The legacy name \code{"TargetBalanced"} is accepted as an alias for \code{"FixedBudget"}.}

 \item{type}{String; the type of the likelihood objects. If \code{"Pairwise"} (default) then the marginal composite likelihood is formed by pairwise marginal likelihoods (see Details).}

 \item{upper}{An optional named list giving upper bounds for parameters when the optimizer is \code{L-BFGS-B}, \code{nlminb}, \code{bobyqa} or \code{optimize}.
 Names must match parameters that are being estimated.}

 \item{varest}{Logical; if \code{TRUE} the estimates' variances and standard errors are returned. For composite likelihood estimation it is deprecated.
 Use \code{sensitivity=TRUE} and update the object using \code{GeoVarest}. Default is \code{FALSE}.}

 \item{weighted}{Logical; if \code{TRUE} the likelihood objects are weighted; see Details for the accepted options. Default is \code{FALSE}.}

 \item{X}{Numeric; matrix of spatio(temporal) covariates in the linear mean specification.}

 \item{spobj}{An object of class \code{sp} or \code{spacetime}.}

 \item{spdata}{Character; the name of data in the \code{sp} or \code{spacetime} object.}


 \item{independence_start}{Logical; if \code{TRUE} (default), \code{GeoFit2}
 attempts a preliminary marginal fit under spatial independence before the final
 fit. The preliminary fit is skipped when the requested fit already has
 \code{type="Independence"}. Bivariate correlation models are not supported by
 \code{GeoFit2}; use \code{GeoFit} directly for those models.}

 \item{independence_optimizer}{String; optimizer used only for the preliminary
 independence fit. Default is \code{"Nelder-Mead"}. The optimizer used for the
 final fit is still specified by \code{optimizer}.}

 \item{warn_independence_failure}{Logical; if \code{TRUE} (default), a warning
 is issued when the preliminary independence fit is unavailable or fails. In
 that case the supplied \code{start} values are used unchanged.}
\item{check.duplicates}{Logical. If \code{TRUE}, perform a fast exact scan for duplicated observation locations at this user entry point. The default is \code{FALSE}, so no duplicate-location scan is imposed. Internal bootstrap, cross-validation, and refitting calls do not repeat the scan.}
}
\details{
 \code{GeoFit2} does not implement a second fitting engine. It is a wrapper
 around the canonical \code{\link{GeoFit}} function.

 When \code{independence_start=TRUE}, the function first constructs a
 preliminary call to \code{GeoFit} with \code{likelihood="Marginal"} and
 \code{type="Independence"}. The preliminary fit estimates the marginal
 parameters available for the selected model. When \code{start} is supplied,
matching marginal entries and omitted mean coefficients are updated as before.
When \code{start=NULL}, all eligible marginal estimates from the independence
fit are used, while dependence parameters use the conservative automatic
initial values prepared by \code{GeoFit}. Parameters supplied in \code{fixed}
always take precedence.

 The final model is then fitted by calling \code{GeoFit} with the original
 requested likelihood, likelihood-object type, optimizer, bounds, pair-selection
 settings, and the updated starting values. Consequently, the objective function,
 parameter constraints, and returned standard \code{GeoFit} components are the
 same as in a direct call to \code{GeoFit}; only the starting values may differ.

 For Gaussian models, \code{sill} is the total marginal variance and
 \code{nugget} attenuates off-diagonal correlation. Therefore an independence
 estimate may initialize \code{sill}, while a user-supplied \code{nugget} is
 preserved unchanged and is never added to or subtracted from \code{sill}.

 Bivariate correlation models are not supported by \code{GeoFit2}; use
 \code{GeoFit} directly. The preliminary initialization is skipped when
 \code{type="Independence"}. If an independence likelihood is not
 implemented for the selected model, or if the preliminary optimization fails,
 the final fit still proceeds using the supplied starting values or, when
\code{start=NULL}, the internal starting values prepared by \code{GeoFit}. This
behavior can be controlled with \code{warn_independence_failure}.

 Stochastic thinning of nearest-neighbor pairs in the final fit can be enabled
 through \code{p_neighb < 1}; \code{thin_method} specifies the thinning
 scheme.
}
\value{
 Returns an object of class \code{GeoFit}.
 An object of class \code{GeoFit} is a list containing at most the following components:
 \item{bivariate}{Logical: \code{TRUE} if the Gaussian RF is bivariate, otherwise \code{FALSE}.}
 \item{clic}{The composite information criterion after a \code{GeoVarest} call; if the full likelihood is considered then it coincides with AIC.}
 \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates.}
 \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates.}
 \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates.}
 \item{coordx_dyn}{A list of dynamical (in time) spatial coordinates.}
 \item{conf.int}{Confidence intervals for standard maximum likelihood estimation.}
 \item{convergence}{A string that denotes if convergence is reached.}
 \item{copula}{The type of copula.}
 \item{corrmodel}{The correlation model.}
 \item{data}{The vector or matrix or array (or list) of data.}
 \item{distance}{The type of spatial distance.}
 \item{fixed}{A list of fixed parameters.}
 \item{iterations}{The number of iterations used by the numerical routine.}
 \item{likelihood}{The configuration of the composite likelihood.}
 \item{logCompLik}{The value of the log composite-likelihood at the maximum.}
 \item{maxdist}{The maximum spatial distance used in the weighted composite likelihood (or \code{NULL}).}
 \item{maxtime}{The maximum temporal-distance threshold used in the composite likelihood, expressed in the same units as \code{coordt}.}
 \item{message}{Extra message passed from the numerical routines.}
 \item{model}{The density associated to the likelihood objects.}
 \item{missp}{\code{TRUE} if a misspecified Gaussian model is used in the composite likelihood.}
 \item{n}{The number of trials in a binomial RF; the number of successes in a negative binomial random field.}
 \item{neighb}{The order of spatial neighborhood in the composite likelihood computation.}
 \item{ns}{The number of (different) location sites in the bivariate case.}
 \item{numcoord}{The number of spatial coordinates.}
 \item{numtime}{The number of temporal realisations of the random field.}
 \item{param}{A list of parameter estimates.}
 \item{radius}{The radius of the sphere in the case of great-circle distance.}
 \item{stderr}{Standard errors for standard maximum likelihood estimation.}
 \item{sensmat}{The sensitivity matrix.}
 \item{varcov}{The variance-covariance matrix of the estimates.}
 \item{type}{The type of the likelihood objects.}
 \item{X}{The matrix of covariates.}
 \item{start_original}{The named list of starting values supplied by the user.}
 \item{start_used}{The named list of starting values passed to the final
 \code{GeoFit} call after any successful marginal initialization.}
 \item{independence_fit}{The preliminary independence \code{GeoFit} object, or
 \code{NULL} when the preliminary fit was skipped or failed.}
 \item{independence_start_message}{A character string describing why the
 preliminary initialization was unavailable or failed, or \code{NULL} when no
 such message was generated.}

}

\references{
 The methodological references for maximum weighted composite-likelihood
 fitting of Gaussian and non-Gaussian random fields are reported in
 \code{\link{GeoFit}}.
}

\seealso{
\code{\link{GeoCovmatrix}} for covariance matrix construction,
  \code{\link{GeoSim}} for simulation,
  \code{\link{GeoKrig}} for prediction,
  \code{\link{GeoVariogram}} and \code{\link{GeoWLS}} for variogram-based tools,
  \code{\link{GeoVarest}} for bootstrap variance estimation.
}

\author{
Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com}, \url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl}, \url{https://www.researchgate.net/profile/Christian-Caamano}
}
\examples{


library(GeoModels)


###############################################################
############ Examples of spatial Gaussian random fields ################
###############################################################

################################################################
###
### Example 1 : Maximum pairwise conditional likelihood fitting 
### of a Gaussian RF with Matern correlation
###
###############################################################
model="Gaussian"
# Define the spatial-coordinates of the points:
set.seed(3)
N=400 # number of location sites
x <- runif(N, 0, 1)
set.seed(6)
y <- runif(N, 0, 1)
coords <- cbind(x,y)

# Define spatial matrix covariates
X=cbind(rep(1,N),runif(N))

# Set the covariance model's parameters:
corrmodel <- "Matern"
mean <- 0.2
mean1 <- -0.5
sill <- 1
nugget <- 0
scale <- 0.2/3
smooth=0.5
param<-list(mean=mean,mean1=mean1,sill=sill,nugget=nugget,scale=scale,smooth=smooth)

# Simulation of the spatial Gaussian RF:
data <- GeoSim(coordx=coords,model=model,corrmodel=corrmodel, param=param,X=X)$data

fixed<-list(nugget=nugget,smooth=smooth)
start<-list(mean=mean,mean1=mean1,scale=scale,sill=sill)

################################################################
###
### Maximum pairwise likelihood fitting of
### Gaussian random fields with exponential correlation.
### 
###############################################################
fit1 <- GeoFit2(data=data,coordx=coords,corrmodel=corrmodel, 
 neighb=3,likelihood="Conditional",
 type="Pairwise", start=start,fixed=fixed,X=X)
print(fit1)




###############################################################
############ Examples of spatial non-Gaussian random fields #############
###############################################################


################################################################
###
### Example 2. Maximum pairwise likelihood fitting of 
### a LogGaussian RF with Generalized Wendland correlation
### 
###############################################################
set.seed(524)
# Define the spatial-coordinates of the points:
N=500
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)
X=cbind(rep(1,N),runif(N))
mean=1; mean1=2 # regression parameters
nugget=0
sill=0.5
scale=0.2
smooth=0

model="LogGaussian"
corrmodel="GenWend"
param=list(mean=mean,mean1=mean1,sill=sill,scale=scale,
 nugget=nugget,power2=4,smooth=smooth)
# Simulation of a non stationary LogGaussian RF:
data <- GeoSim(coordx=coords, corrmodel=corrmodel,model=model,X=X,
 param=param)$data

fixed<-list(nugget=nugget,power2=4,smooth=smooth)
start<-list(mean=mean,mean1=mean1,scale=scale,sill=sill)
I=Inf
lower<-list(mean=-I,mean1=-I,scale=0,sill=0)
upper<-list(mean= I,mean1= I,scale=I,sill=I)

# Maximum pairwise composite-likelihood fitting of the RF:
fit <- GeoFit2(data=data,coordx=coords,corrmodel=corrmodel, model=model,
 neighb=3,likelihood="Conditional",type="Pairwise",X=X,
 optimizer="nlminb",lower=lower,upper=upper,
 start=start,fixed=fixed)
print(unlist(fit$param))


################################################################
###
### Example 3. Maximum pairwise likelihood fitting of
### SinhAsinh random fields with Wendland0 correlation
###
###############################################################
set.seed(261)
model="SinhAsinh"
# Define the spatial-coordinates of the points:
x <- runif(500, 0, 1)
y <- runif(500, 0, 1)
coords <- cbind(x,y)

corrmodel="Wend0"
mean=0;nugget=0
sill=1
skew=-0.5
tail=1.5
power2=4
c_supp=0.2

# model parameters
param=list(power2=power2,skew=skew,tail=tail,
 mean=mean,sill=sill,scale=c_supp,nugget=nugget)
data <- GeoSim(coordx=coords, corrmodel=corrmodel,model=model, param=param)$data

plot(density(data))
fixed=list(power2=power2,nugget=nugget)
start=list(scale=c_supp,skew=skew,tail=tail,mean=mean,sill=sill)
# Maximum pairwise likelihood:
fit1 <- GeoFit2(data=data,coordx=coords,corrmodel=corrmodel, model=model,
 neighb=3,likelihood="Marginal",type="Pairwise",
 start=start,fixed=fixed)
print(unlist(fit1$param))



}
\keyword{Composite}
