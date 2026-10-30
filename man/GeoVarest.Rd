\name{GeoVarest}
\alias{GeoVarest}
\encoding{UTF-8}
\title{Score-based variance estimation for \code{GeoFit} objects}

\description{
 The function updates a fitted \code{GeoFit} object by estimating the
 variability matrix of the composite likelihood score through parametric simulation.
 The fitted model is used to generate \code{K} independent datasets. For each simulated dataset, the composite likelihood score is evaluated at the original estimate \eqn{\hat\theta}, without refitting the model. The empirical
 variance of these simulated scores provides an estimate of the variability matrix \eqn{J}. Together with the sensitivity matrix \eqn{H}, computed by
 \code{GeoFit} when \code{sensitivity = TRUE}, this yields the Godambe
 sandwich covariance matrix
 \deqn{
 G^{-1} = H^{-1} J H^{-1}.
 }
 The updated object contains standard errors, Wald confidence intervals,
 p-values, the estimated matrices \eqn{J}, \eqn{H^{-1}} and \eqn{G^{-1}},
 and composite likelihood information criteria based on the penalty
 \eqn{\mathrm{tr}(H^{-1}J)}.
}
\usage{
GeoVarest(fit, K = 100, sparse = FALSE,
 method = c("cholesky", "TB", "CE"),
 alpha = 0.95, L = 10000,
 parallel = FALSE, ncores = 6, progress = TRUE, seed = NULL,
 min_success_rate = 0.8)
}
\arguments{
 \item{fit}{A fitted object obtained from \code{\link{GeoFit}}. The object
 must contain the sensitivity matrix, hence \code{GeoFit} should be called
 with \code{sensitivity = TRUE}. Full/Standard likelihood fits are not accepted;
 use \code{GeoFit(..., likelihood = "Full", type = "Standard", varest = TRUE)}
 for Hessian-based standard errors.}

 \item{K}{The number of simulations used in the parametric score bootstrap.}

 \item{sparse}{Logical; if \code{TRUE}, then Cholesky decomposition is
 performed using sparse matrix algorithms.}

 \item{method}{String; the method of simulation. The default is
 \code{"cholesky"}. For large data sets the options \code{"TB"} and
 \code{"CE"} call approximate simulation methods; see
 \code{\link{GeoSimapprox}}.}

 \item{alpha}{Numeric; the level of the confidence interval.}

 \item{L}{Numeric; the number of lines in the turning bands method.}

 \item{parallel}{Logical; default \code{FALSE}. If \code{TRUE}, the score evaluation step is
 parallelized.}

 \item{ncores}{Positive integer or \code{NULL}; default \code{6}. With \code{parallel=TRUE}, an explicit integer requests that many workers, capped by detected cores and the number of bootstrap jobs. Set \code{ncores=NULL} for automatic selection, capped at six workers and normally leaving one detected core free.}

 \item{progress}{Logical; if \code{TRUE}, progress information is shown.}

 \item{seed}{Optional integer seed for reproducibility of the simulated
 samples.}
 \item{min_success_rate}{Minimum fraction of successful score-bootstrap
 replications required to return variance estimates. The default is 0.8.}
}
\details{
For spatio-temporal fits, the original numeric \code{coordt} values stored in \code{fit} are reused. Irregularly spaced times are supported with \code{method = "cholesky"}; approximate simulation methods retain the restrictions documented in \code{\link{GeoSimapprox}}.


 Let \eqn{cl(\theta)} denote the composite log-likelihood and let
 \eqn{U(\theta) = \nabla cl(\theta)} be the corresponding composite score.
 The function simulates \code{K} data sets from the fitted model and evaluates
 the composite score at the fitted parameter value \eqn{\hat\theta}. The
 variability matrix is estimated as
 \deqn{\hat J = Var\{U_1(\hat\theta), \ldots, U_K(\hat\theta)\}.}

 If \eqn{H} is the sensitivity matrix stored in \code{fit$sensmat}, the
 inverse Godambe matrix is estimated by
 \deqn{\widehat{G^{-1}} = H^{-1} \hat J H^{-1}.}
 Standard errors are obtained from the square root of the diagonal of
 \eqn{\widehat{G^{-1}}}.

 For composite likelihoods, the penalty used in the information criterion is
 \deqn{tr(H^{-1}\hat J),}
 and the composite likelihood information criterion is computed as
 \deqn{-2 cl(\hat\theta) + 2 tr(H^{-1}\hat J).}

 Differently from \code{\link{GeoVarestbootstrap}}, this function does not
 refit the model for each simulated data set. It estimates the variability
 matrix of the score and then computes the sandwich/Godambe covariance
 matrix.

 For stochastic nearest-neighbor fits, the realized retained-pair graph stored
 in the original fit is reused unchanged for every simulated data set. Thus the
 score bootstrap estimates variability conditional on the selected pair graph,
 in agreement with the current implementation described in the stochastic-NN
 methodology.


Parallel multisession workers are started with the package-library paths of the calling R session, including the library containing the loaded GeoModels installation.

For \code{method = "TB"} without a copula, parallel workers are persistent and replications are streamed one at a time (simulate, evaluate, release). This keeps peak memory bounded when both \code{L} and \code{K} are large while still reusing the selected pair graph and static fitting context within each worker. For fitted copula models, \code{method="TB"} is simulated through the turning-bands backend of \code{GeoSimCopula}; the copula path is currently not streamed replicate-by-replicate. On platforms where \pkg{future} reports forked multicore execution as safe, GeoModels uses it automatically to reduce worker startup and serialization overhead; otherwise it falls back to multisession with the current GeoModels library path propagated explicitly.
}
\value{
 Returns an updated object of class \code{GeoFit}. The following components
 are added or updated:
 \item{stderr}{Estimated standard errors obtained from the inverse Godambe
 matrix.}
 \item{varcov}{Estimated inverse Godambe matrix \eqn{\widehat{G^{-1}}}.}
 \item{godambe}{Estimated Godambe matrix.}
 \item{Jmat}{Estimated variability matrix of the composite score.}
 \item{Hinv}{Inverse, or generalized inverse, of the sensitivity matrix.}
 \item{claic}{Composite likelihood AIC-type criterion.}
 \item{clic}{Same value as \code{claic}.}
 \item{clbic}{Composite likelihood BIC-type criterion.}
 \item{clic_penalty}{Penalty term \eqn{tr(H^{-1}\hat J)}.}
 \item{conf.int}{Wald-type confidence intervals based on the estimated
 standard errors.}
 \item{pvalues}{Wald-type p-values.}
 \item{scores}{Matrix of successful bootstrap score evaluations.}
 \item{score_logCompLik}{Composite log-likelihood values corresponding to the
 successful score evaluations.}
 \item{score_failures}{Data frame with failed score evaluations, if any.}
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
\author{
Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},
\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl},
\url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},
\url{https://www.researchgate.net/profile/Christian-Caamano}
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
NN=150
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


#fit_update=GeoVarest(fit,K=100,parallel=TRUE)
#fit_update$stderr
#fit_update$conf.int
#fit_update$pvalues

}
\keyword{Composite}
