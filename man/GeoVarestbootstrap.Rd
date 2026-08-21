\name{GeoVarestbootstrap} 
\alias{GeoVarestbootstrap}
\encoding{UTF-8}
\title{Update a \code{GeoFit} object using parametric bootstrap for std error estimation.}
\description{
 The procedure updates a \code{GeoFit} object using a
 classical parametric bootstrap. The fitted model is simulated and refitted
 for each bootstrap replication; the empirical covariance matrix of the
 successful parameter estimates is used as the bootstrap covariance estimate.
 A faster score-bootstrap alternative is \code{GeoVarest}. Full/Standard
 likelihood fits are not accepted; obtain their Hessian-based standard errors
 directly with \code{GeoFit(..., varest = TRUE)}.}
\usage{
GeoVarestbootstrap(
  fit, K = 100, sparse = FALSE,
  optimizer = NULL, lower = NULL, upper = NULL,
  method = "cholesky", alpha = 0.95, L = 10000,
  parallel = FALSE, ncores = 6, progress = TRUE,
  seed = NULL, min_success_rate = 0.8
)
}
\arguments{
 \item{fit}{A fitted object obtained from the
 \code{\link{GeoFit}}.}
 \item{K}{The number of simulations in the parametric bootstrap.}
 \item{sparse}{Logical; if \code{TRUE} then cholesky decomposition is performed
 using sparse matrices algorithms (spam packake).}
 \item{optimizer}{The type of optimization algorithm (see \code{\link{GeoFit}} for details). If NULL then the optimization algorithm of 
 the object fit is chosen.}
 \item{lower}{An optional named list giving the values for the lower bound of the space parameter
 when the optimizer is \code{L-BFGS-B} or \code{nlminb} or \code{optimize}.}
 \item{upper}{An optional named list giving the values for the upper bound of the space parameter
 when the optimizer is \code{L-BFGS-B} or \code{nlminb} or \code{optimize}.}
 \item{method}{String; The method of simulation. Default is \code{cholesky}. For large data set three options are \code{TB} or \code{CE} (see the \code{GeoSimapprox}) function. } 
 \item{alpha}{ Numeric; The level of the confidence interval. }
 \item{L}{Numeric; the number of lines in the turning band method.} 
 \item{parallel}{Logical; default \code{FALSE}. If \code{TRUE} then 
 the estimation step is parallelized 
 }
 \item{ncores}{Positive integer or \code{NULL}; default \code{6}. With \code{parallel=TRUE}, an explicit integer requests that many workers, capped by detected cores and the number of bootstrap jobs. Set \code{ncores=NULL} for automatic selection, capped at six workers and normally leaving one detected core free.}
 \item{progress}{Logic; If TRUE then a progress bar is shown.}
  \item{seed}{Optional integer seed for reproducibility of the simulated
 samples.}
  \item{min_success_rate}{Minimum fraction of successful bootstrap refits
 required to return variance estimates. The default is 0.8.}
 }
\details{
For spatio-temporal fits, the original numeric \code{coordt} values stored in \code{fit} are reused. Irregularly spaced times are supported with \code{method = "cholesky"}; approximate simulation methods retain the restrictions documented in \code{\link{GeoSimapprox}}.

 The function updates a \code{GeoFit} object by refitting each simulated
data set. A sensitivity matrix is not required to estimate bootstrap standard
errors and confidence intervals. If \code{fit$sensmat} is available, CLIC and
CLBIC penalties are also computed. For stochastic nearest-neighbor fits, the
realized retained-pair graph from the original fit is reused unchanged across
all bootstrap refits.

Parallel multisession workers are started with the package-library paths of the calling R session, including the library containing the loaded GeoModels installation.

For \code{method = "TB"} without a copula, parallel workers are persistent and replications are streamed one at a time (simulate, evaluate, release). This keeps peak memory bounded when both \code{L} and \code{K} are large while still reusing the selected pair graph and static fitting context within each worker. For fitted copula models, \code{method="TB"} is simulated through the turning-bands backend of \code{GeoSimCopula}; the copula path is currently not streamed replicate-by-replicate. On platforms where \pkg{future} reports forked multicore execution as safe, GeoModels uses it automatically to reduce worker startup and serialization overhead; otherwise it falls back to multisession with the current GeoModels library path propagated explicitly.
}
\value{
 Returns an updated object of class \code{GeoFit}. The main updated components
 include \code{stderr}, \code{varcov}, \code{godambe} when invertible,
\code{conf.int}, \code{pvalues}, and \code{estimates}. When a compatible
\code{sensmat} is available, \code{claic}/\code{clic}, \code{clbic},
and \code{clic_penalty} are also updated.
}



\section{Three-dimensional coordinates}{
For a purely spatial univariate fit with explicit irregular three-dimensional
Euclidean coordinates, \code{method = "TB"} is supported when no
\code{anisopars} were used. Method \code{"CE"} remains restricted to
regular two-dimensional grids. Use \code{method = "cholesky"} for
bivariate, spatio-temporal, anisotropic, or other unsupported
three-dimensional cases.
}
\seealso{
\code{\link{GeoFit}} for the fitted objects used as input
}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{


library(GeoModels)

################################################################
###
### Example 1. Test on the parameter
### of a regression model using conditional composite likelihood
###
###############################################################
set.seed(342)
model="Gaussian" 
# Define the spatial-coordinates of the points:
NN=3500
x = runif(NN, 0, 1)
y = runif(NN, 0, 1)
coords = cbind(x,y)
# Parameters
mean=1; mean1=-1.25; # regression parameters
 sill=1 # variance

# matrix covariates
X=cbind(rep(1,nrow(coords)),runif(nrow(coords)))

# model correlation 
corrmodel="Matern"
smooth=0.5;scale=0.1; nugget=0;

# simulation
param=list(smooth=smooth,mean=mean,mean1=mean1,
 sill=sill,scale=scale,nugget=nugget)
data = GeoSim(coordx=coords, corrmodel=corrmodel,
 model=model, param=param,X=X)$data

I=Inf

fixed=list(nugget=nugget,smooth=smooth)
start=list(mean=mean,mean1=mean1,scale=scale,sill=sill)

lower=list(mean=-I,mean1=-I,scale=0,sill=0)
upper=list(mean=I,mean1=I,scale=I,sill=I)
# Maximum pairwise composite-likelihood fitting of the RF:
fit = GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
 likelihood="Conditional",type="Pairwise",sensitivity=TRUE,
 lower=lower,upper=upper,neighb=3,
 optimizer="nlminb",X=X,
 start=start,fixed=fixed)

unlist(fit$param)


#fit_update=GeoVarestbootstrap(fit,K=100,parallel=TRUE)
#fit_update$stderr
#fit_update$conf.int
#fit_update$pvalues

}


\keyword{Composite}
