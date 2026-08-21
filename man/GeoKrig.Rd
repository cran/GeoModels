\name{GeoKrig}
\alias{GeoKrig}
\encoding{UTF-8}
\title{Spatial (bivariate) and spatio temporal optimal linear prediction 
 for Gaussian and non-Gaussian random fields.}
\description{
  For a given set of spatial location sites (and temporal instants),
the function computes optimal linear prediction and associated mean square error 
for the Gaussian and non-Gaussian case.
}
\usage{
GeoKrig(estobj=NULL,data, coordx=NULL, coordy=NULL, coordz=NULL, coordt=NULL, 
coordx_dyn=NULL, corrmodel,distance="Eucl",
 grid=FALSE, loc, 
 method="cholesky", model="Gaussian", n=1,nloc=NULL,mse=FALSE, 
 param, anisopars=NULL,radius=1, sparse=FALSE,
 time=NULL, type_krig="Simple",weigthed=TRUE,which=1,
 copula=NULL, X=NULL,Xloc=NULL,Mloc=NULL,spobj=NULL,spdata=NULL,varcov=NULL,
 progress=FALSE,check.duplicates=FALSE)
}
\arguments{
 \item{estobj}{An object of class Geofit that includes information about data, model and estimates.}
 \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation) or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid) or a
 (\eqn{t \times d}{t x d})-matrix (a single spatio-temporal realisation) or an (\eqn{d \times d \times t \times n }{d x d x t})-array
 (a single spatio-temporal realisation on regular grid) giving the data used for prediction.}
\item{coordx}{A numeric \eqn{d \times 2}{d x 2} or \eqn{d \times 3}{d x 3}
coordinate matrix.  It may be omitted when \code{coordx_dyn} or \code{spobj}
is supplied.  Coordinates on a sphere for a fixed radius \code{radius} are
passed in lon/lat format expressed in decimal degrees.}
 \item{coordy}{A numeric vector giving 1-dimension of
 spatial coordinates; Optional argument, the default is \code{NULL}.}
 \item{coordz}{A numeric vector giving 1-dimension of
 spatial coordinates; Optional argument, the default is \code{NULL}.}
 \item{coordt}{A numeric vector giving the temporal coordinates of the observations. The default is \code{NULL}, in which case a spatial random field is expected. Temporal coordinates may be irregularly spaced; temporal lags are computed from the supplied coordinate values.}
 \item{coordx_dyn}{For dynamic observation locations, a list with one two- or three-column coordinate matrix per element of \code{coordt}.  The rows of \code{coordx_dyn[[t]]} correspond to \code{data[[t]]} and to the same temporal block of \code{X}.  Dynamic observation locations do not alter the ordering of prediction tasks.  See \code{\link{GeoModels-spacetime-ordering}}.}
 \item{corrmodel}{String; the name of a correlation model, for the
 see \code{\link{GeoCovmatrix}} for the list of implemented correlation models.}
 \item{distance}{String; the name of the spatial distance. The default
 is \code{Eucl}, the Euclidean distance. See \code{\link{GeoFit}} for details.}
 \item{grid}{Logical; if \code{FALSE} (the default) the data used for prediction
 are interpreted as spatial or spatio-temporal realisations on a set
 of non-equispaced spatial sites (irregular grid).}
 \item{loc}{A numeric (\eqn{n \times 2}{n x 2})-matrix (where
 \code{n} is the number of spatial sites) giving 2-dimensions of spatial coordinates to be predicted.}
 \item{method}{String; matrix decomposition used to solve the kriging
system.  The choices are \code{cholesky} (default) and \code{svd}.
Only \code{cholesky} is available when \code{sparse=TRUE}.}
 \item{n}{Positive integer size parameter. For Binomial it may be scalar or contain one value per observation. For direct Negative Binomial it is the single common number \eqn{r} of successes; Geometric corresponds to \eqn{r=1}.}
 \item{nloc}{Positive integer size at prediction tasks. For Binomial, if observation-side \code{n} is location-specific then \code{nloc} is required and is never inferred by averaging \code{n}; otherwise a scalar \code{n} is reused when \code{nloc=NULL}. For direct Negative Binomial, \code{nloc}, if supplied, must equal the common \eqn{r}; for Geometric it must equal 1.}
 \item{mse}{Logical; if \code{TRUE}, compute the MSE of the kriging predictor.
 The default is \code{FALSE}; in that case the returned \code{mse} component is
 \code{NULL}, including for fixed- and dynamic-support space-time prediction.} 
 \item{model}{String; the type of RF and therefore the densities associated to the likelihood
 objects. \code{Gaussian} is the default, see the Section
 Details.} 
 \item{param}{A list of parameter values required for the correlation model. See Details for the accepted options.}
 \item{anisopars}{ A list of two elements: "angle" and "ratio" i.e. the anisotropy angle and the anisotropy ratio, respectively.}
 \item{radius}{Numeric: the radius of the sphere if coordinates are passed in lon/lat format;Default value is 1.}
 \item{sparse}{Logical; if \code{TRUE} kriging is computed with sparse matrices algorithms 
 using spam package. Default is FALSE. It should be used with compactly supported covariances.} 
 \item{time}{A numeric vector giving the temporal instants to be predicted. Values need not be equally spaced and are interpreted on the same numeric time scale as \code{coordt}. The default is \code{NULL}, in which case only spatial prediction is performed.}
 \item{type_krig}{String; \code{Simple} (default) treats the supplied or
fitted mean coefficients as known.  \code{Universal} uses the same plug-in
coefficients and, when \code{mse=TRUE}, adds their estimation uncertainty
from \code{varcov} to the prediction MSE.}
 \item{weigthed}{Logical legacy argument retained for backward
compatibility.  It has no effect in standard kriging.}
 \item{which}{Numeric; In the case of bivariate cokriging it indicates which variable to predict.
 	 It can be 1 or 2}
 \item{copula}{String; optional copula specification. Linear prediction is implemented for \code{"Gaussian"}, the constructive Clayton-like copula \code{"Clayton"}, and \code{"SkewGaussian"}.}
 \item{X}{Numeric design matrix at the observations.  For fixed-location space-time data, rows follow \code{c(t(data))}, i.e. time then site.  For dynamic observations, use a stacked matrix ordered by temporal blocks, or a list aligned with \code{coordx_dyn}.}
 \item{Xloc}{Numeric design matrix at prediction tasks.  For space-time prediction its rows are location-major: all requested times for \code{loc[1, ]}, then all requested times for \code{loc[2, ]}, and so on.}
 \item{Mloc}{Numeric vector giving the known marginal location predictor at prediction tasks, in the same location-major order as \code{Xloc}. For Gaussian additive margins this is the mean itself; for models with link-scale parametrizations (for example Gamma/Weibull/LogGaussian or Beta2) it is on the same scale as \code{param$mean}. Use either \code{Mloc} or \code{Xloc}, not both.}
 \item{spobj}{An object of class \code{sp} or \code{spacetime}. For space-time objects, the current \code{sp2Geo()} conversion uses sequential temporal indices and does not preserve irregular original time spacing; to retain irregular temporal distances, use the explicit-coordinate interface with numeric \code{coordt}.}
 \item{spdata}{Character:The name of data in the sp or spacetime object}
 \item{varcov}{Covariance matrix of the estimated parameters.  It is
required for \code{type_krig="Universal"} when \code{mse=TRUE}.  For
composite-likelihood fits this is the sandwich covariance, i.e. the inverse
Godambe information, usually taken automatically from a \code{GeoFit}
object.  Mean coefficients absent from \code{varcov} are treated as fixed.}
 \item{progress}{Logical; if \code{TRUE}, show progress while the optimized
blocked prediction path processes multiple prediction blocks.  The default is
\code{FALSE}.  No progress bar is shown when the calculation uses a single
block or a legacy one-shot prediction path.}
\item{check.duplicates}{Logical. If \code{TRUE}, perform a fast exact scan for duplicated observation locations at this user entry point. The default is \code{FALSE}, so no duplicate-location scan is imposed. Internal bootstrap, cross-validation, and refitting calls do not repeat the scan.}
}

\details{

For univariate kriging, the mean at the observations is either
\eqn{X\beta}, where the entries of \eqn{\beta} are represented by
\code{mean}, \code{mean1}, and so on, or a site-specific vector supplied as
\code{param$mean}.  With \code{X=NULL}, the
model is intercept-only and uses the scalar coefficient \code{mean}.

The prediction mean is specified independently: use \code{Xloc} to obtain
\eqn{X_{loc}\beta}, or use \code{Mloc} to supply the prediction mean directly.
\code{Mloc} and \code{Xloc} are mutually exclusive.  If the observation mean
is site-specific, \code{Mloc} is required.  When \code{X} is used at the
observations, either \code{Xloc} or \code{Mloc} must be supplied; the latter
allows a known prediction mean that is not generated by the fitted design.

For \code{type_krig="Universal"}, prediction uses the same fitted mean
coefficients as \code{type_krig="Simple"}.  Thus coefficients estimated by
maximum likelihood, REML, or composite likelihood are not recomputed inside
\code{GeoKrig}.  When \code{mse=TRUE}, the MSE is augmented by
\eqn{(X_{loc}-\lambda^T X)\,Var(\widehat\beta)\,
(X_{loc}-\lambda^T X)^T}, using the mean-parameter block of \code{varcov}.
For composite likelihood, \code{varcov} should be the sandwich covariance
(the inverse Godambe information).  Coefficients absent from that matrix are
regarded as fixed and contribute zero uncertainty.  External mean vectors
supplied through \code{param$mean} or \code{Mloc} are already known, so a
requested universal prediction is treated as simple kriging with a warning.


For prediction, names beginning with \code{Gaussian_misp_} are interpreted as
estimation specifications rather than new data-generating models.  They are
automatically mapped to the corresponding response model before covariance and
prediction calculations (for example \code{Gaussian_misp_Poisson} to
\code{Poisson}).  Without a copula, prediction is available only for models with
a validated model-specific covariance/predictor implementation.  In particular,
plain \code{Logistic}, \code{Beta2}, \code{Kumaraswamy}, and
\code{Kumaraswamy2} are rejected rather than falling through an incomplete branch;
these margins are available in the validated continuous-copula prediction path.

For discrete models with nonlinear marginal means, the Universal-kriging MSE uses
the derivative of the marginal mean with respect to the regression predictor when
propagating the mean-parameter block of \code{varcov}.  The MSE is evaluated from
its diagonal directly and does not construct a prediction-by-prediction matrix.

For \code{copula="Gaussian"}, \code{copula="Clayton"}, and \code{copula="SkewGaussian"}, prediction is performed on the observed marginal scale using marginal means, marginal variances, and the copula-induced cross-covariances. The supported continuous margins are \code{"Gaussian"}, \code{"StudentT"}, \code{"LogGaussian"}, \code{"Gamma"}, \code{"Weibull"}, \code{"Beta"}, \code{"Beta2"}, \code{"Kumaraswamy"}, \code{"Kumaraswamy2"}, \code{"Logistic"}, \code{"SkewLaplace"}, \code{"Tukeyh"}, \code{"Tukeyh2"}, and \code{"SinhAsinh"}. The count margins \code{"Poisson"}, \code{"Binomial"}, and \code{"BinomialNeg"} are also supported by global \code{GeoKrig} for Gaussian, Clayton-like, and skew-Gaussian copulas; the Binomial copula path currently requires a common trial count \code{n}. For a Gaussian copula, the last three margins use their exact transformed-Gaussian covariance rather than a truncated Hermite approximation. For location-dependent margins such as \code{"Beta2"}, \code{"Kumaraswamy2"}, and the positive-scale models, the covariance between two sites uses both supplied location predictors. Thus \code{X}/\code{Xloc} or external means are propagated into the covariance consistently. The resulting predictor is the optimal linear predictor on the observed scale; except when the joint field is Gaussian it is generally different from the full conditional mean. For \code{type_krig="Universal"}, the MSE correction uses the derivative of the marginal mean with respect to the regression predictor (delta-method correction for nonlinear marginal means). For the Clayton-like and skew-Gaussian copulas, the expensive covariance ingredients are precomputed and cached. Clayton-like prediction uses deterministic copula quadrature tables, while skew-Gaussian prediction uses a cached bivariate Hermite expansion of the latent Gaussian representation; neither path performs adaptive two-dimensional integration separately for each covariance entry. For \code{copula="SkewGaussian"}, \code{param$nu} is the paper's \eqn{\eta\in(-1,1)}; when \code{nu=0} the covariance calculation is routed exactly to the Gaussian-copula implementation.

For the Clayton-like copula, direct adaptive two-dimensional integration for every covariance entry would be prohibitively expensive. GeoModels therefore uses a deterministic cached covariance engine: the Clayton-like copula density is evaluated once on a Gauss--Legendre quadrature grid over a compact grid of the underlying squared correlations, the discrete joint distributions are numerically balanced to preserve uniform margins, and covariance values are obtained by cached matrix products and interpolation. The expensive copula kernel is reused across prediction locations and across calls in the same R session. Positive-scale margins exploit their exact multiplicative scaling, while \code{"Beta2"} and \code{"Kumaraswamy2"} use an additional compact cached grid on the location-predictor scale, built only over the range required by the current kriging problem. This approximation is designed for kriging covariance construction and avoids one numerical double integral per pair.

For global prediction with a supported continuous copula, \code{GeoKrig()} constructs observation--prediction cross-covariances in bounded blocks of prediction locations. The observation covariance matrix is factorized once and the same Cholesky (or SVD, when requested) factorization is reused for every block. This changes only the evaluation order: predictions, MSE values, and kriging weights use the same covariance equations as the one-shot calculation. The historical full weights matrix is retained in the returned object, so blocking reduces the peak memory associated with the cross-covariance matrix without changing the returned results. The ordinary Gaussian model without a copula uses the same blocked path. When \code{progress=TRUE} and more than one block is required, progress is updated once per completed prediction block; the numerical calculation and block sizes are unchanged.

 Best linear unbiased predictor and associated mean square error is computed 
 for Gaussian and some non-Gaussian cases.
 Specifically, for a spatial or spatio-temporal or spatial bivariate dataset, given a set of spatial locations and
 temporal istants and a correlation model
\code{corrmodel} with some fixed parameters and given the type of RF (\code{model}) the function computes
 simple or universal kriging, for the specified spatial locations
 \code{loc} and temporal instants \code{time}, providing also the respective mean square error.
 For the choice of the spatial or spatio temporal correlation model see details in \code{\link{GeoCovmatrix}} function.
 The list \code{param} specifies mean and covariance parameters, see \code{\link{CorrParam}} and \code{\link{GeoCovmatrix}} for details. The \code{type_krig} parameter indicates the type of kriging. In the
case of simple kriging, the known mean can be specified by the parameter
 \code{mean} in the list \code{param} (See examples).
}

\section{Spatio-temporal ordering}{
Observed fixed-location data use time-major order: a \eqn{T \times N}{T x N}
matrix is vectorized as \code{c(t(data))}.  Dynamic observations are supplied
as aligned lists \code{data[[t]]} and \code{coordx_dyn[[t]]}, concatenated by
time.  The rows of \code{X} and any known observation mean use the same
observation order.

Prediction tasks use the different, location-major order
\preformatted{
loc[1, ] at time[1], ..., loc[1, ] at time[Tloc],
loc[2, ] at time[1], ..., loc[2, ] at time[Tloc], ...
}
Rows of \code{Xloc} and elements of \code{Mloc} must follow this order.  The
returned \code{pred} and \code{mse} objects are \eqn{Tloc \times Nloc}{Tloc x Nloc}
matrices with prediction times in rows and prediction locations in columns.
See \code{\link{GeoModels-spacetime-ordering}}.
}

\value{
 Returns an object of class \code{Kg}.
 An object of class \code{Kg} is a list containing
 at most the following components:
 \item{bivariate}{\code{TRUE} if spatial bivariate cokriging is performed, otherwise \code{FALSE};}
 \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates used for prediction;}
 \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates used for prediction;}
 \item{coordz}{A \eqn{d}{d}-dimensional vector of spatial coordinates used for prediction;}
 \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates used for prediction;}
 \item{corrmodel}{String: the correlation model;}
 \item{covmatrix}{The covariance matrix. }
 \item{data}{The vector or matrix or array of data used for prediction}
 \item{distance}{String: the type of spatial distance;}
 \item{grid}{\code{TRUE} if the spatial data used for prediction are observed in a regular grid, otherwise \code{FALSE};}
 \item{loc}{A (\eqn{n \times 2}{n x 2})-matrix of spatial locations to be predicted.}
 \item{n}{The Binomial number of trials, or the common Negative-Binomial number \eqn{r} of successes.} 
 \item{nozero}{In the case of sparse simple kriging the percentage of non zero values in the covariance matrix.
 Otherwise is NULL.}
 \item{numcoord}{Numeric:he number \eqn{d}{d} of spatial coordinates used for prediction;}
 \item{numloc}{Numeric: the number \eqn{n}{n} of spatial coordinates to be predicted;}
 \item{numtime}{Numeric: the number \eqn{d}{d} of the temporal instants used for prediction;}
 \item{numt}{Numeric: the number \eqn{m}{m} of the temporal instants to be predicted;}
 \item{model}{The response model used for prediction after canonicalizing any
 misspecified-Gaussian fitting name.}
 \item{fit_model}{The model name supplied by the caller or stored in the input
 \code{GeoFit} object before prediction canonicalization.}
 \item{param}{The original parameter list supplied to prediction, before any
 internal standardization used to construct covariance matrices;}
 \item{pred}{For spatio-temporal prediction, a \eqn{m \times n}{m x n} matrix with rows corresponding to \code{time} and columns corresponding to rows of \code{loc}; for spatial prediction, a numeric vector.}
 \item{radius}{Numeric: the radius of the sphere if coordinates are pssed in lon/lat format;}
 \item{spacetime}{\code{TRUE} if spatio-temporal kriging and \code{FALSE} if
 spatial kriging;}
 \item{time}{A \eqn{m}{m}-dimensional vector of temporal coordinates to be predicted;}
 \item{type}{String: the type of kriging (Standard).}
 \item{type_krig}{String: the type of kriging (simple or universal)}
 \item{mse}{When \code{mse=TRUE}, the prediction MSE (a time-by-location matrix
 for space-time prediction and a vector for spatial prediction); otherwise
 \code{NULL}.}

}


\references{
 Gaetan, C. and Guyon, X. (2010) \emph{Spatial Statistics and Modelling}.
 Springer-Verlag, New York.
}

\seealso{
\code{\link{GeoCovmatrix}} for covariance matrix construction,
  \code{\link{GeoKrigloc}} and \code{\link{GeoKriglocWeights}} for local kriging,
  \code{\link{GeoFit}} for parameter estimation.
}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{

library(GeoModels)
################################################################
########### Examples of spatial kriging ############
################################################################

################################################################
###
### Example 1. Spatial kriging of a
### Gaussian random fields with Gen wendland correlation.
###
################################################################

model="Gaussian"
set.seed(79)
x = runif(300, 0, 1)
y = runif(300, 0, 1)
coords=cbind(x,y)
# Set the exponential cov parameters:
corrmodel = "GenWend"
mean=0; sill=5; nugget=0
scale=0.2;smooth=0;power2=4

param=list(mean=mean,sill=sill,nugget=nugget,scale=scale,smooth=smooth,power2=power2)

# Simulation of the spatial Gaussian random field:
data = GeoSim(coordx=coords, corrmodel=corrmodel,model=model,
 param=param)$data

## estimation with pairwise likelihood
fixed=list(nugget=nugget,smooth=0,power2=power2)
start=list(mean=0,scale=scale,sill=1)
I=Inf
lower=list(mean=-I,scale=0,sill=0)
upper=list(mean= I,scale=I,sill=I)
# Maximum pairwise likelihood fitting :
fit = GeoFit(data, coordx=coords, corrmodel=corrmodel,model=model,
 likelihood='Marginal', type='Pairwise',neighb=3,
 optimizer="nlminb", lower=lower,upper=upper,
 start=start,fixed=fixed)

# locations to predict
xx=seq(0,1,0.03)
loc_to_pred=as.matrix(expand.grid(xx,xx))

## first option
#param=append(fit$param,fit$fixed)
#pr=GeoKrig(loc=loc_to_pred,coordx=coords,corrmodel=corrmodel,
# model=model,param=param,data=data,mse=TRUE)

## second option using object GeoFit
pr=GeoKrig(fit,loc=loc_to_pred,mse=TRUE)


colour = rainbow(100)

opar=par(no.readonly = TRUE)
par(mfrow=c(1,3))
if (requireNamespace("fields", quietly = TRUE)) {
 fields::quilt.plot(coords, data, col = colour)
}
# simple kriging map prediction
if (requireNamespace("fields", quietly = TRUE)) {
 fields::image.plot(
  xx, xx, matrix(pr$pred, ncol = length(xx)), col = colour,
  xlab = "", ylab = "", main = " Kriging "
 )
}

# simple kriging MSE map prediction variance
if (requireNamespace("fields", quietly = TRUE)) {
 fields::image.plot(
  xx, xx, matrix(pr$mse, ncol = length(xx)), col = colour,
  xlab = "", ylab = "", main = "Std error"
 )
}
par(opar)

################################################################
###
### Example 2. Spatial kriging of a Skew
### Gaussian random fields with Matern correlation.
###
################################################################
model="SkewGaussian"
set.seed(79)
x = runif(300, 0, 1)
y = runif(300, 0, 1)
coords=cbind(x,y)
# Set the exponential cov parameters:
corrmodel = "Matern"
mean=0
sill=2
nugget=0
scale=0.1
smooth=0.5
skew=3
param=list(mean=mean,sill=sill,nugget=nugget,scale=scale,smooth=smooth,skew=skew)

# Simulation of the spatial skew Gaussian random field:
data = GeoSim(coordx=coords, corrmodel=corrmodel,model=model,
 param=param)$data

fixed=list(nugget=nugget,smooth=smooth)
start=list(mean=0,scale=scale,sill=1,skew=skew)
I=Inf
lower=list(mean=-I,scale=0,sill=0,skew=-I)
upper=list(mean= I,scale=I,sill=I,skew=I)
# Maximum pairwise likelihood fitting :
fit = GeoFit2(data, coordx=coords, corrmodel=corrmodel,model=model,
 likelihood='Marginal', type='Pairwise',neighb=3,
 optimizer="nlminb", lower=lower,upper=upper,
 start=start,fixed=fixed)

# locations to predict
xx=seq(0,1,0.03)
loc_to_pred=as.matrix(expand.grid(xx,xx))
## optimal linear kriging
pr=GeoKrig(fit,loc=loc_to_pred,mse=TRUE)

colour = rainbow(100)

opar=par(no.readonly = TRUE)
par(mfrow=c(1,3))
if (requireNamespace("fields", quietly = TRUE)) {
 fields::quilt.plot(coords, data, col = colour)
}
# simple kriging map prediction
if (requireNamespace("fields", quietly = TRUE)) {
 fields::image.plot(
  xx, xx, matrix(pr$pred, ncol = length(xx)), col = colour,
  xlab = "", ylab = "", main = " Kriging "
 )
}

# simple kriging MSE map prediction variance
if (requireNamespace("fields", quietly = TRUE)) {
 fields::image.plot(
  xx, xx, matrix(pr$mse, ncol = length(xx)), col = colour,
  xlab = "", ylab = "", main = "Std error"
 )
}
par(opar)

################################################################
###
### Example 3. Spatial kriging of a 
### Gamma random field with mean spatial regression
###
###############################################################
set.seed(312)
model="Gamma"
corrmodel = "GenWend" 
# Define the spatial-coordinates of the points:
NN=300
coords=cbind(runif(NN),runif(NN))
## matrix covariates
a0=rep(1,NN)
a1=runif(NN,0,1)
X=cbind(a0,a1)
##Set model parameters
shape=2
## regression parameters
mean = 1;mean1= -0.2
# correlation parameters
nugget = 0;power2=4
scale = 0.3;smooth=0 

## simulation
param=list(shape=shape,nugget=nugget,mean=mean,mean1=mean1, 
 scale=scale,power2=power2,smooth=smooth)
data = GeoSim(coordx=coords,corrmodel=corrmodel, param=param,
 model=model,X=X)$data

#####starting and fixed parameters
fixed=list(nugget=nugget,power2=power2,smooth=smooth)
start=list(mean=mean,mean1=mean1, scale=scale,shape=shape)

## estimation with pairwise likelihood
fit2 = GeoFit(data=data,coordx=coords,corrmodel=corrmodel,X=X,
 neighb=3,likelihood="Conditional",type="Pairwise",
 start=start,fixed=fixed, model = model)

# locations to predict with associated covariates
xx=seq(0,1,0.03)
loc_to_pred=as.matrix(expand.grid(xx,xx))
NP=nrow(loc_to_pred)
a0=rep(1,NP)
a1=runif(NP,0,1)
Xloc=cbind(a0,a1)

#optimal linear kriging 
pr=GeoKrig(fit2,loc=loc_to_pred,Xloc=Xloc,sparse=TRUE,mse=TRUE)

## map 
opar=par(no.readonly = TRUE)
par(mfrow=c(1,3))
if (requireNamespace("fields", quietly = TRUE)) {
 fields::quilt.plot(coords, data, main = "Data")
}
map=matrix(pr$pred,ncol=length(xx))
mapmse=matrix(pr$mse,ncol=length(xx))
if (requireNamespace("fields", quietly = TRUE)) {
 fields::image.plot(xx, xx, map, xlab = "", ylab = "", main = "Kriging ")
}

if (requireNamespace("fields", quietly = TRUE)) {
 fields::image.plot(xx, xx, mapmse, xlab = "", ylab = "", main = "MSE")
}
par(opar)


################################################################
########### Examples of spatio temporal kriging ############
################################################################

################################################################
###
### Example 4. Spatio temporal simple kriging of n locations
### sites and m temporal instants for a Gaussian random fields
### with estimated double Wendland correlation.
###
###############################################################
model="Gaussian"
# Define the spatial-coordinates of the points:
x = runif(300, 0, 1)
y = runif(300, 0, 1)
coords=cbind(x,y)
times=1:4

# Define model correlation modek and associated parameters
corrmodel="Wend0_Wend0"
param=list(nugget=0,mean=0,power2_s=4,power2_t=4,
 scale_s=0.2,scale_t=2,sill=1)

# Simulation of the space time Gaussian random field:
set.seed(31)
data=GeoSim(coordx=coords,coordt=times,corrmodel=corrmodel,sparse=TRUE,
 param=param)$data

# Maximum pairwise likelihood fitting of the space time random field:
start = list(scale_s=0.15,scale_t=2,sill=1,mean=0)
fixed = list(nugget=0,power2_s=4,power2_t=4)

fit = GeoFit(data, coordx=coords, coordt=times, model=model, corrmodel=corrmodel, 
 likelihood='Conditional', type='Pairwise',start=start,fixed=fixed,
 neighb=3,maxtime=1)

# locations to predict
xx=seq(0,1,0.04)
loc_to_pred=as.matrix(expand.grid(xx,xx))
# Define the times to predict
times_to_pred=2

pr=GeoKrig(fit,loc=loc_to_pred,time=times_to_pred,sparse=TRUE,mse=TRUE)

opar=par(no.readonly = TRUE)
par(mfrow=c(1,3))
zlim=c(-2.5,2.5)
colour = rainbow(100)


if (requireNamespace("fields", quietly = TRUE)) {
 fields::quilt.plot(
  coords, data[2, ], col = colour, main = " data at Time 2"
 )
}
if (requireNamespace("fields", quietly = TRUE)) {
 fields::image.plot(
  xx, xx, matrix(pr$pred, ncol = length(xx)), col = colour,
  main = " Kriging at Time 2", ylab = ""
 )
}
if (requireNamespace("fields", quietly = TRUE)) {
 fields::image.plot(
  xx, xx, matrix(pr$mse, ncol = length(xx)), col = colour,
  main = "Std err Time at time 2", ylab = ""
 )
}


par(opar)


################################################################
###
### Example r. Spatial bivariate simple cokriging of n locations
### sites for a bivariate Gaussian random fields
### with estimated Matern correlation.
###
###############################################################
#set.seed(6)
#NN=1500 # number of spatial locations 
#x = runif(NN, 0, 1); 
#y = runif(NN, 0, 1) 
#coords=cbind(x,y) 

## setting parameters
#mean_1 = 2; mean_2= -1
#nugget_1 =0;nugget_2=0
#sill_1 =0.5; sill_2 =1; 

### correlation parameters
#CorrParam("Bi_Matern")
#scale_1=0.2/3; scale_2=0.15/3; scale_12=0.5*(scale_2+scale_1) 
#smooth_1=smooth_2=smooth_12=0.5
#pcol = -0.4
#param= list(nugget_1=nugget_1,nugget_2=nugget_2,
# sill_1=sill_1,sill_2=sill_2,
# mean_1=mean_1,mean_2=mean_2,
# smooth_1=smooth_1, smooth_2=smooth_2,smooth_12=smooth_12,
# scale_1=scale_1, scale_2=scale_2,scale_12=scale_12,
# pcol=pcol)

## simulation
#data = GeoSim(coordx=coords, corrmodel="Bi_Matern",model=model,param=param)$data

#fixed=list(mean_1=mean_1,mean_2=mean_2, nugget_1=nugget_1,nugget_2=nugget_2, 
# smooth_1=smooth_1, smooth_2=smooth_2,smooth_12=smooth_12)

#start=list( sill_1=sill_1,sill_2=sill_2,
# scale_1=scale_1,scale_2=scale_2,scale_12=scale_12, pcol=pcol)

## estimation with maximum likelihood 
#fit = GeoFit(data=data,coordx=coords, corrmodel="Bi_Matern",
 #likelihood="Marginal",type="Pairwise",optimizer="BFGS",neighb=5,
 #start=start,fixed=fixed)

###### co-kriging for the fist component ##############
#xx=seq(0,1,0.022)
#loc_to_pred=as.matrix(expand.grid(xx,xx))
#pr1 = GeoKrig(fit,which=1,mse=TRUE,loc=loc_to_pred)

#opar=par(no.readonly = TRUE)
#par(mfrow=c(1,2))
#zlim=c(-2.5,2.5)
#colour = rainbow(100)
#fields::quilt.plot(coords,data[1,] ,col=colour,main = paste(" Fist component")) 
#fields::quilt.plot(loc_to_pred,pr1$pred,col=colour,
# main = paste(" Kriging first component"),ylab="")
#par(opar)


}
\keyword{Kriging}
