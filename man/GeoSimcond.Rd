\name{GeoSimcond}
\alias{GeoSimcond}
\title{Conditional simulation of spatial Gaussian and non-Gaussian random fields}
\description{
 Performs global or nearest-neighbour local conditional simulation for univariate
 spatial random fields at explicit prediction coordinates. Cholesky simulation is
 available together with the approximate turning-bands method when supported by
 the selected correlation model.  The combination \code{local=TRUE} and
 \code{method="TB"} avoids global observation covariance matrices and is intended
 for large spatial datasets.
}
\usage{
GeoSimcond(estobj = NULL, data, coordx, coordy = NULL, coordz = NULL, coordt = NULL,
 coordx_dyn = NULL, corrmodel, distance = "Eucl", grid = FALSE, loc,
 maxdist = NULL, maxtime = NULL, method = "Cholesky", model = "Gaussian",
 n = 1, nrep = 1, local = FALSE, L = 1000, neighb = NULL,
 param, anisopars = NULL, radius = 1, sparse = FALSE, time = NULL,
 copula = NULL, X = NULL, Xloc = NULL, Mloc = NULL,
 parallel=FALSE, ncores = 6, progress=FALSE, n_iter=25L,
 check.duplicates=FALSE, nloc=NULL, mcmc_thin=1L)
}
\arguments{
 \item{estobj}{Object of class \code{GeoFit} containing model information}
 \item{data}{Numeric vector/matrix/array of observed data}
 \item{coordx}{Numeric matrix with one row per observed location and two or three spatial coordinate columns}
 \item{coordy}{Optional numeric vector of y-coordinates}
 \item{coordz}{Optional numeric vector of z-coordinates}
 \item{coordt}{Retained for API compatibility. Space-time conditional simulation is not currently implemented.}
 \item{coordx_dyn}{Retained for API compatibility. Dynamic coordinates are not currently implemented.}
 \item{corrmodel}{String specifying correlation model name}
 \item{distance}{String specifying distance metric (default: "Eucl")}
 \item{grid}{Must currently be \code{FALSE}; supply explicit coordinates.}
 \item{loc}{Numeric matrix of prediction locations (n x 2)}
 \item{maxdist}{Optional maximum distance for local kriging}
 \item{maxtime}{Optional maximum temporal distance}
 \item{method}{Unconditional simulation method. Currently \code{"Cholesky"} and \code{"TB"} are available; \code{"CE"} is retained for API compatibility but rejected by the explicit-coordinate conditional workflow.}
 \item{model}{String specifying random field type (default: "Gaussian")}
 \item{n}{For direct \code{Binomial}, the number of trials, either scalar or one positive integer per observed location. For direct \code{BinomialNeg}, the common positive integer number \eqn{r} of successes.}
 \item{nrep}{Number of retained conditional simulation replicates (default: 1). For direct count models these are successive thinned draws from the latent Gibbs chain.}
 \item{local}{Logical. If \code{FALSE} (default), global conditioning is used. If
 \code{TRUE}, each prediction location is conditioned only on the observations
 selected by \code{neighb} and/or \code{maxdist}.}
 \item{L}{Number of lines for turning bands method (default: 1000)}
 \item{neighb}{Optional positive integer giving the number of nearest observed
 locations used for local conditional simulation.  When \code{local=TRUE}, at
 least one of \code{neighb} or \code{maxdist} must be supplied.}
 \item{param}{List of parameter values}
 \item{anisopars}{List with anisotropy angle and ratio}
 \item{radius}{Radius used by spherical distance calculations (default: 1)}
 \item{sparse}{Must currently be \code{FALSE}.}
 \item{time}{Retained for API compatibility. Space-time conditional simulation is not currently implemented.}
 \item{copula}{Optional string specifying copula type}
 \item{X}{Optional design matrix for the marginal location at the observed locations.}
 \item{Xloc}{Optional design matrix for the marginal location at prediction locations. If \code{X} is supplied and the margin uses a location parameter, either \code{Xloc} or \code{Mloc} must be supplied.}
 \item{Mloc}{Optional vector of known/fitted marginal location values at prediction locations; an alternative to \code{Xloc}.}
 \item{parallel}{Logical; default \code{FALSE}. If \code{TRUE}, supported expensive stages use parallel workers. For \code{method="TB"} or \code{"CE"}, the setting is passed to \code{GeoSimapprox}; for \code{local=TRUE}, local kriging weights may also be computed in parallel. The final local substitution is evaluated in the main R process using a compact sparse weight operator and bounded batches, avoiding duplication of large simulation and weight objects across workers.}
 \item{ncores}{Positive integer or \code{NULL}; default \code{6}. With \code{parallel=TRUE}, an explicit integer requests that many workers, capped by detected cores and available jobs. Set \code{ncores=NULL} for automatic selection, capped at six workers and normally leaving one detected core free. A value of 1 is treated as serial execution.}
 \item{progress}{If TRUE then a progress bar is shown.}
 \item{n_iter}{Positive integer; number of Gibbs sweeps used by latent conditional samplers. For direct \code{Binomial} and \code{BinomialNeg}, it is the burn-in length of the count-constrained Gibbs chain; retained simulations are then separated by \code{mcmc_thin} sweeps. The default for these direct count models is 25. For the other latent Gibbs engines, omitting \code{n_iter} retains the historical default of 1000 sweeps. An explicitly supplied value is always respected. It is ignored by the Gaussian copula conditional simulator.}

\item{check.duplicates}{Logical. If \code{TRUE}, perform a fast exact scan for duplicated observation locations at this user entry point. The default is \code{FALSE}, so no duplicate-location scan is imposed. Internal bootstrap, cross-validation, and refitting calls do not repeat the scan.}
 \item{nloc}{Prediction-side count size. For direct \code{Binomial}, a scalar or one positive integer per prediction location; it is required when observation-side \code{n} varies by site. For direct \code{BinomialNeg}, if supplied it must equal the common \eqn{r}.}
 \item{mcmc_thin}{Positive integer thinning interval for the direct Binomial and Negative-Binomial count-constrained latent Gibbs sampler. The default is 1.}
}
\details{
The currently validated domain is univariate, purely spatial, non-dynamic and
non-grid conditional simulation. Both global and nearest-neighbour local
conditioning are available for the Gaussian latent-substitution paths described
below. Unsupported domain combinations are rejected before covariance matrices
or unconditional simulations are constructed.

Without a copula, the currently supported models are
\code{Gaussian}, \code{LogGaussian}, \code{Tukeyh}, \code{Tukeyh2},
\code{SinhAsinh}, \code{SkewGaussian}, \code{Gamma}, \code{Weibull}, \code{Binomial}, and \code{BinomialNeg}.  Direct non-copula \code{Gamma}
conditional simulation uses the finite Gaussian-square construction and therefore
requires \code{param$shape} to be a positive integer; non-integer values are
rejected rather than rounded. For direct and Gaussian-copula Tukey margins,
\code{Tukeyh} requires \code{0 <= tail < 0.5}; \code{Tukeyh2} requires
\code{0 <= tail1 < 0.5} and \code{0 <= tail2 < 0.5}, with \code{tail1}
the right-tail and \code{tail2} the left-tail parameter. Zero is valid and
recovers the Gaussian transformation on that side. \code{SinhAsinh} requires
a strictly positive \code{tail}.

For the direct \code{Binomial} and \code{BinomialNeg} random fields, conditional
simulation uses the repeated latent Gaussian threshold construction defining the
models.  At each observed site the latent Gaussian copies are updated as one
count-constrained block.  Their signs are sampled exactly from a conditional
Bernoulli distribution by a log-scale dynamic program, after which the latent
Gaussian magnitudes are sampled from the corresponding univariate truncated
normal full conditionals.  For \code{Binomial}, the constraint is
\deqn{\sum_{l=1}^{n_i} I\{Z_l(s_i)>0\}=y_i.}
For \code{BinomialNeg}, writing \eqn{t_i=y_i+r}, the constraint is
\deqn{\sum_{l=1}^{t_i-1} I\{Z_l(s_i)>0\}=r-1, \qquad Z_{t_i}(s_i)>0.}
The native sampler maintains the dense Gaussian precision update in compiled
code and directly generates the prediction-side latent fields, so the full
\eqn{L\times n_{obs}\times n_{rep}} latent array is not returned to R.
The first \code{n_iter} full site sweeps are discarded and retained conditional
realisations are separated by \code{mcmc_thin} sweeps.  This direct count path
currently requires \code{method="Cholesky"} and \code{local=FALSE}.

With \code{copula="Gaussian"}, conditional simulation is implemented for the
continuous margins \code{Gaussian}, \code{StudentT}, \code{LogGaussian},
\code{Gamma}, \code{Weibull}, \code{Beta}, \code{Beta2},
\code{Kumaraswamy}, \code{Kumaraswamy2}, \code{Logistic},
\code{SkewLaplace}, \code{Tukeyh}, \code{Tukeyh2}, and
\code{SinhAsinh}. With \code{copula="Clayton"} or
\code{copula="SkewGaussian"}, the currently validated continuous margins are
the same list except for \code{Tukeyh}, \code{Tukeyh2}, and
\code{SinhAsinh}.
For copula-based \code{Gamma}, \code{Weibull}, \code{Beta},
\code{Beta2}, \code{Kumaraswamy}, and \code{Kumaraswamy2} margins,
\code{param$sill} is not required: their marginal dispersion is determined by
the corresponding shape parameters.  The \code{sill} parameter is required
only for copula margins whose marginal scale is explicitly parameterized by it.
For the skew--Gaussian copula, \code{param$nu} is the bounded asymmetry
parameter \eqn{\eta\in(-1,1)}; the latent Gibbs sampler uses
\eqn{\gamma_\eta=\eta/\sqrt{1-\eta^2}}.

For the constructive Clayton-like copula, \code{param$nu} must be a positive
integer.  If \eqn{u_i} is the observed copula uniform, the sampler uses
\eqn{r_i=u_i^{2/\nu}} and the latent representation
\deqn{r_i=\frac{\sum_{k=1}^{\nu} Z_{ki}^2}{\sum_{k=1}^{\nu} Z_{ki}^2+W_{1i}^2+W_{2i}^2}.}
At each observed site the \eqn{\nu+2} Gaussian latent variables are updated
with a Gibbs kernel that preserves this ratio exactly.  Directional full
conditionals are von Mises--Fisher and the common latent radius is updated from
a one-dimensional log-concave full conditional.  The retained finite-sweep
state is therefore a Monte Carlo approximation to the Clayton latent
conditional distribution after \code{n_iter} sweeps.  Conditional Gaussian
simulation of the latent fields at \code{loc} is then followed by the Clayton
ratio reconstruction and the location-specific marginal quantile.

Discrete copula margins are not routed through the former mid-PIT
approximation, because that approximation is not exact conditioning on the
latent intervals.

For Gaussian fields, conditional simulations use the substitution method
\deqn{Z_c(s_0)=E\{Z(s_0)\mid Z(s)=z\}+Z^*(s_0)-WZ^*(s),}
where the unconditional residual field \eqn{Z^*} is generated with zero mean
and the covariance specified by \code{param}. The declared orientation of the
kriging weights is used explicitly.

With \code{local=TRUE}, the same substitution identity is applied separately at
each prediction location using only its selected neighborhood \eqn{N_j}:
\deqn{Z_c(s_j)=\hat Z_j+Z^*(s_j)-w_j^T Z^*_{N_j}.}
The local systems are solved once, and GeoSimcond keeps a compact map of
neighbor indices and weights.  When possible this map is assembled once as a
\code{spam} sparse matrix with \eqn{n_{loc}} rows and \eqn{n_{obs}} columns.
Conditional corrections for multiple realizations are then evaluated in bounded
matrix batches rather than by repeating the local graph reduction for every
replicate.  The sparse representation contains only the selected local weights;
GeoSimcond does not retain one local covariance matrix per prediction site.
Prediction locations with no observation inside a \code{maxdist}-only
neighborhood receive zero conditioning correction, so their local conditional
approximation reduces to the marginal mean plus the unconditional residual at
that location.

Local conditioning is currently validated for \code{Gaussian},
\code{LogGaussian}, \code{Tukeyh}, \code{Tukeyh2}, and \code{SinhAsinh}
without a copula, and for all supported continuous margins with
\code{copula="Gaussian"}.  In particular, \code{Tukeyh}, \code{Tukeyh2}
and \code{SinhAsinh} are also supported as continuous margins under the
Gaussian copula and use the same latent-Gaussian fast path.  The Gibbs-based constructive
\code{SkewGaussian}, \code{Gamma}, \code{Weibull}, \code{Binomial}, and \code{BinomialNeg} paths and the
\code{Clayton} and \code{SkewGaussian} copulas still require
\code{local=FALSE}.

For the sinh--arcsinh model, the implementation uses the monotone transformation
\deqn{Y(s)=\mu(s)+\sqrt{\sigma^2}\,\sinh\{[\operatorname{asinh}\{Z(s)\}+\eta]/\nu\},}
and conditions on the Gaussian scale using its exact inverse
\deqn{Z(s)=\sinh\{\nu\,\operatorname{asinh}[(Y(s)-\mu(s))/\sqrt{\sigma^2}]-\eta\}.}
Thus \code{LogGaussian}, \code{Tukeyh}, \code{Tukeyh2} and
\code{SinhAsinh} require only one latent Gaussian conditional simulation per
replicate; no non-Gaussian Gibbs step is introduced.

For margins with a location parameter, location-specific covariates are preserved:
if \code{X} is supplied for the observations, prediction requires either
\code{Xloc} or \code{Mloc}; the function does not silently replace a varying
prediction-side location by the intercept.

For continuous copula margins, observations are transformed by their fitted
marginal CDF to the copula scale, conditional simulation is performed on the
latent copula construction, and each realization is transformed back with the
location-specific marginal quantile.  For the Gaussian copula this is standard
Gaussian conditional simulation.  The skew--Gaussian case uses its latent Gibbs
sampler, while the Clayton-like case uses the Gaussian-square ratio Gibbs
sampler described above.  In both MCMC cases, a retained finite-sweep draw is
a Monte Carlo approximation to the target conditional distribution after
\code{n_iter} sweeps.

\code{method="TB"} makes the unconditional simulation step approximate.  For
large datasets it can be combined with \code{local=TRUE}; the turning-bands
simulation then avoids a dense covariance factorization for the unconditional
field, while local conditioning avoids the global observation covariance system.
For a fixed neighborhood size \eqn{m}, local preprocessing consists of small
\eqn{m \times m} kriging systems and the substitution correction uses only the
retained local weights. For multiple conditional realizations, the same sparse
local weight operator is reused and applied in bounded matrix batches. This
optimization is shared automatically by the one-Gaussian monotone transformations
and by Gaussian-copula margins. When \code{parallel=TRUE}, approximate
unconditional simulation and local-weight construction honor \code{ncores}; the
final sparse substitution remains in the main R process to avoid copying the
unconditional simulations and weight map to additional worker processes. Spatial circulant
embedding is not currently available
because it requires a regular grid, whereas the validated conditional workflow
combines explicit observation and prediction coordinates. The nugget, distance,
radius and anisotropy settings are propagated to the latent Gaussian simulations.
}
\value{
An object of class \code{GeoSimcond} containing \code{condsim}, a list with one
element per successful conditional replicate. Each element is a numeric vector of
length \code{numloc}, in prediction-location order. Thus, if a dense matrix is
needed for post-processing, \code{do.call(rbind, lapply(x$condsim, as.numeric))}
produces a matrix with replications in rows and prediction locations in columns.
The object also contains \code{cond_mean} and \code{cond_var}, computed across
replicates using a streaming update without materializing an additional full
replicate-by-location matrix, together with the model, parameter, conditioning
(\code{local}, \code{neighb}) and simulation-method metadata. With a single
replicate, \code{cond_var} contains \code{NA} values because a sample variance
cannot be estimated.
}
\references{
Gaetan, C. and Guyon, X. (2010) \emph{Spatial Statistics and Modelling}.
Springer Verlag, New York.

Blasi, F., Caamaño-Carrillo, C., Bevilacqua, M. and Furrer, R. (2022).
A selective view of climatological data and likelihood estimation.
\emph{Spatial Statistics}, 50, 100596.

Caamaño-Carrillo, C., Bevilacqua, M., López, C. and Morales-Oñate, V. (2024).
Nearest neighbors weighted composite likelihood based on pairs for (non-)Gaussian
massive spatial data with an application to Tukey-hh random fields estimation.
\emph{Computational Statistics and Data Analysis}, 191, 107887.
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
sim_result <- GeoSimcond(fit,loc = loc_to_sim,nrep=50)

cond_mean=sim_result$cond_mean # conditional mean
cond_var =sim_result$cond_var # conditional var

# Empirical pointwise intervals from the conditional simulations.
# Rows are replications and columns are prediction locations.
sim_mat <- do.call(rbind, lapply(sim_result$condsim, as.numeric))
interval_summary <- t(apply(sim_mat, 2, quantile,
                            probs = c(0.025, 0.5, 0.975)))
colnames(interval_summary) <- c("q025", "q500", "q975")
head(interval_summary)

par(mfrow=c(1,3))
if (requireNamespace("fields", quietly = TRUE)) fields::quilt.plot(coords, data)
if (requireNamespace("fields", quietly = TRUE)) fields::quilt.plot(loc_to_sim, cond_mean)
if (requireNamespace("fields", quietly = TRUE)) fields::quilt.plot(loc_to_sim, cond_var)
par(mfrow=c(1,1))

##############################################
## conditional simulation of a LogGaussian rf 
##############################################
model="LogGaussian"
set.seed(79)
### conditioning locations
x = runif(500, 0, 1)
y = runif(500, 0, 1)
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
sim_result <- GeoSimcond(fit,loc = loc_to_sim,nrep=50)

cond_mean=sim_result$cond_mean # conditional mean
cond_var =sim_result$cond_var # conditional var

par(mfrow=c(1,3))
if (requireNamespace("fields", quietly = TRUE)) fields::quilt.plot(coords, data)
if (requireNamespace("fields", quietly = TRUE)) fields::quilt.plot(loc_to_sim,cond_mean)
if (requireNamespace("fields", quietly = TRUE)) fields::quilt.plot(loc_to_sim,cond_var)
par(mfrow=c(1,1))

}
