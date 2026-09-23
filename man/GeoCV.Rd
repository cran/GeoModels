\name{GeoCV}
\alias{GeoCV}
\encoding{UTF-8}
\title{Repeated holdout kriging cross-validation}
\description{
The procedure uses the \code{\link{GeoKrig}} or \code{\link{GeoKrigloc}}
function to compute repeated holdout kriging cross-validation using information from a
\code{\link{GeoFit}} object. The function returns prediction scores.
}

\usage{
GeoCV(
  fit, K = 100, estimation = TRUE,
  optimizer = NULL, lower = NULL, upper = NULL,
  n.fold = 0.05, local = FALSE,
  neighb = NULL, maxdist = NULL, maxtime = NULL,
  sparse = FALSE, type_krig = "Simple", which = 1,
  parallel = FALSE, ncores = 6, progress = TRUE,
  seed = NULL, predictor = "linear", conditional_nrep = 1000L,
  conditional_n_iter = 25L, conditional_mcmc_thin = 1L
)
}

\arguments{
 \item{fit}{An object of class \code{\link{GeoFit}}.}

 \item{K}{Positive integer greater than or equal to 2 giving the number of cross-validation iterations. Non-integer values are rejected rather than rounded.}

 \item{estimation}{Logical; if \code{TRUE}, the model is re-estimated on the training observations at each iteration and the fold-specific estimates are used for prediction. If \code{FALSE}, the parameter estimates in \code{fit} are reused; this is therefore a conditional predictive assessment given parameters estimated from the complete dataset, not a fully refitted out-of-sample cross-validation.}

 \item{optimizer}{The type of optimization algorithm if \code{estimation} is
 \code{TRUE}. See \code{\link{GeoFit}} for details. If \code{NULL}, then the
 optimization algorithm stored in \code{fit} is used.}

 \item{lower}{An optional named list giving the values for the lower bounds of
 the parameters when bounded optimization is used and \code{estimation} is
 \code{TRUE}.}

 \item{upper}{An optional named list giving the values for the upper bounds of
 the parameters when bounded optimization is used and \code{estimation} is
 \code{TRUE}.}

 \item{n.fold}{Numeric; the fraction of observations randomly deleted and
 predicted in each cross-validation iteration. At least one training and one
 prediction observation are retained. In the space-time case, sampling also
 retains at least one training observation at every observed time.}

 \item{local}{Logical; if \code{TRUE}, local kriging is performed. The default
 is \code{FALSE}.}

 \item{neighb}{Numeric; an optional positive integer indicating the order of
 neighborhood if local kriging is performed.}

 \item{maxdist}{Numeric; an optional positive value indicating the spatial
 neighborhood distance if local kriging is performed.}

 \item{maxtime}{Numeric; an optional non-negative temporal-distance threshold, expressed in the same units as the fitted \code{coordt}, when local kriging is performed.}

 \item{sparse}{Logical; if \code{TRUE}, kriging and simulation are computed with
 sparse matrix algorithms using the \pkg{spam} package. The default is
 \code{FALSE}. It should be used with compactly supported covariance models.}

 \item{type_krig}{String; the type of kriging. If \code{"Simple"}, the default, the fitted mean coefficients are treated as plug-in values. If \code{"Universal"}, the prediction MSE is additionally corrected using the covariance matrix of the estimated mean coefficients. With \code{estimation = TRUE}, Universal cross-validation is supported for full-likelihood refits, for which \code{varest = TRUE} is requested automatically in each fold. It is deliberately rejected for composite-likelihood refits because each fold would require its own Godambe covariance matrix. With \code{estimation = FALSE}, a stored \code{fit$varcov} may be used; for composite likelihood it can be obtained with \code{GeoVarest}.}

 \item{which}{Numeric; in the case of bivariate cokriging, it indicates which
 variable to predict. It can be \code{1} or \code{2}.}


 \item{predictor}{Character vector selecting the prediction method. The default
 \code{"linear"} preserves the historical behaviour based on \code{GeoKrig} or
 \code{GeoKrigloc}. For direct spatial \code{Binomial} and \code{BinomialNeg}
 random fields without a copula, \code{"conditional"} uses the conditional mean
 estimated by \code{GeoSimcond}. Supplying
 \code{c("linear", "conditional")} evaluates both predictors on exactly the same
 held-out observations and, when \code{estimation=TRUE}, from the same fold-specific
 refit. Conditional prediction is currently restricted to purely spatial,
 univariate, non-copula direct count models with \code{local=FALSE} and
 \code{sparse=FALSE}.}

 \item{conditional_nrep}{Positive integer; number of retained conditional
 simulations used by \code{GeoSimcond} to estimate the conditional mean within
 each cross-validation iteration. Used only when
 \code{predictor} includes \code{"conditional"}. Default is 1000.}

 \item{conditional_n_iter}{Positive integer; burn-in length, in full Gibbs
 sweeps, for the direct count conditional sampler used within each
 cross-validation iteration. Used only when \code{predictor} includes
 \code{"conditional"}. Default is 25.}

 \item{conditional_mcmc_thin}{Positive integer; thinning interval for retained
 direct-count conditional simulations within each cross-validation iteration.
 Used only when \code{predictor} includes \code{"conditional"}. Default is 1.}

 \item{parallel}{Logical; default \code{FALSE}. Set \code{TRUE} to evaluate the cross-validation iterations in parallel.}

 \item{ncores}{Positive integer or \code{NULL}; default \code{6}. With \code{parallel=TRUE}, an explicit integer requests that many workers, capped by detected cores and the number of cross-validation jobs. Set \code{ncores=NULL} for automatic selection, capped at six workers and normally leaving one detected core free. Non-integer values are rejected.}

 \item{progress}{Logical; if \code{TRUE}, a progress bar is shown.}

 \item{seed}{Optional finite integer seed used to make the random selection of folds and any stochastic refitting step reproducible. Non-integer values are rejected. If \code{NULL}, the current random number generator state is left unchanged. When supplied, the previous \code{.Random.seed} is restored on exit.}
}

\details{
The option \code{GeoModels.cv_strict = TRUE} rejects any excluded observation,
non-finite prediction, or invalid Gaussian predictive MSE. The default is
\code{FALSE}: exclusions are counted explicitly, and a fold with no usable
predictions returns \code{NA} scores rather than an unreported change of denominator.
Warnings raised by iterations are collected and summarized in the calling R
session; \code{options(warn = 2)} continues to treat warnings as errors.

For a spatio-temporal \code{GeoFit} object, the stored temporal coordinates are reused without imposing equal spacing. Local temporal neighborhoods interpret \code{maxtime} as a distance threshold in the same units as those coordinates.


The function randomly removes a fraction \code{n.fold} of the observations at
each iteration, predicts the removed observations, and computes predictive
scores. With the default \code{predictor="linear"}, prediction uses kriging as
in previous versions. For supported direct binomial and negative-binomial random
fields, \code{predictor="conditional"} instead estimates the conditional mean
from repeated calls to the count-constrained conditional simulator. When both
predictors are requested, the fold is constructed and the model is refitted
only once, so the comparison is paired within each cross-validation iteration.

If \code{estimation = TRUE}, the model is re-estimated at each cross-validation iteration before prediction. Refits preserve the fitted likelihood settings, including pair weighting, pair thinning, distance-memory handling and anisotropy. For misspecified estimators, the response model stored in \code{fit$model} remains distinct from the working likelihood stored in \code{fit$estimation_model}; the latter is used for each refit. Estimated anisotropy parameters are passed once through \code{anisopars}, avoiding duplicate \code{angle}/\code{ratio} entries in the starting or fixed parameter lists.


For conditional prediction of direct \code{Binomial} and \code{BinomialNeg}
random fields, each fold-specific fitted object is passed to
\code{GeoSimcond}. The conditional predictor is the empirical mean of
\code{conditional_nrep} retained conditional simulations after
\code{conditional_n_iter} burn-in sweeps, with thinning
\code{conditional_mcmc_thin}. Parallelisation, when requested, is performed over
cross-validation iterations; the nested \code{GeoSimcond} call is deliberately
run serially within each worker to avoid nested parallel execution.

For Universal kriging with fold-specific re-estimation, full-likelihood fits request \code{varest = TRUE} within each fold. Composite-likelihood refits are not assigned the Godambe covariance matrix from the full dataset: \code{GeoCV} instead rejects \code{estimation = TRUE, type_krig = "Universal"} for composite likelihood, because a statistically coherent analysis would require a fold-specific \code{GeoVarest} calculation.

If \code{estimation = FALSE}, the parameter estimates from the original fit are reused after removing the validation observations from the conditioning data. Consequently, the reported errors assess prediction conditional on parameters that were estimated using the complete dataset and may be more optimistic than a fully refitted out-of-sample cross-validation.

For \code{Poisson}, \code{Binomial}, and \code{BinomialNeg} with \code{copula="SkewGaussian"}, global cross-validation (\code{local=FALSE}) supports both fold refitting and fixed-parameter prediction. Discrete Clayton-like copula models still require \code{estimation=FALSE} until their pairwise likelihood is implemented.

Regular-grid fits are converted to explicit point coordinates after observations
are removed, because each training sample is no longer a complete Cartesian
grid. Space-time covariates may be stored either as one matrix in observation
order or as a list containing one matrix per time.

For Gaussian fitted models all documented scores are computed from the Gaussian
predictive mean and MSE. For non-Gaussian fitted models, only RMSE, MAE and MAD
are returned; \code{brie}, \code{crps}, \code{lscore}, \code{pit},
\code{intscore} and \code{coverage} are set to \code{NA}, because an exact
non-Gaussian predictive distribution is not available from only a mean and MSE.

When \code{seed} is supplied, the cross-validation samples are reproducible. The
function preserves and restores the user's random number generator state.
}

\value{
Returns a list containing predictive score vectors. With the default
\code{predictor="linear"}, the historical top-level structure is unchanged.
When conditional prediction is requested, components \code{linear} and/or
\code{conditional} contain the corresponding score lists. For the conditional
predictor, RMSE, MAE and MAD are populated and the Gaussian-only probability
scores remain \code{NA}. For backward compatibility, the top-level score
aliases refer to the linear predictor when it is requested, and otherwise to
the conditional predictor.
 \item{linear}{When requested, a list of score vectors for the linear predictor.}
 \item{conditional}{When requested, a list of score vectors for the conditional-mean predictor.}
 \item{rmse}{The vector of root mean squared errors for the primary predictor.}
 \item{mae}{The vector of mean absolute errors for the primary predictor.}
 \item{mad}{The vector of median absolute errors for the primary predictor.}
 \item{brie}{The vector of Brier scores, or \code{NA} for non-Gaussian fits.}
 \item{crps}{The vector of continuous ranked probability scores, or
 \code{NA} for non-Gaussian fits.}
 \item{lscore}{The vector of log-scores, or \code{NA} for non-Gaussian fits.}
 \item{pit_mean}{The vector of mean probability integral transform values within each holdout fold, or \code{NA} for non-Gaussian fits. This is a descriptive PIT summary, not a calibration score: a mean near 0.5 does not by itself imply a uniform PIT distribution.}
 \item{pit}{Backward-compatible alias of \code{pit_mean}.}
 \item{intscore}{The vector of interval scores, or \code{NA} for
 non-Gaussian fits.}
 \item{coverage}{The vector of empirical coverage values, or \code{NA} for
 non-Gaussian fits.}
 \item{n_requested}{Number of requested predictions in each iteration.}
 \item{n_valid}{Number of finite observation/prediction pairs used for error scores.}
 \item{n_failed}{Number of excluded pairs; \code{n_requested - n_valid}.}
 \item{n_invalid_mse}{Among valid pairs, the number with non-finite or nonpositive
 predictive MSE. \code{NA} when Gaussian probability scores are not applicable.}
 \item{n_valid_prob}{Number of pairs used for Gaussian probability scores;
 \code{NA} when not applicable.}
 \item{warnings}{A frequency table of warning messages from all iterations,
 or \code{NULL} when none were raised.}
 \item{seed}{The seed used for reproducibility, or \code{NULL}.}
 \item{predictor}{The prediction method or methods requested.}
 \item{conditional_nrep}{When conditional prediction is requested, the number of retained conditional simulations per fold.}
 \item{conditional_n_iter}{When conditional prediction is requested, the Gibbs burn-in length per fold.}
 \item{conditional_mcmc_thin}{When conditional prediction is requested, the Gibbs thinning interval.}
}

\seealso{\code{\link{GeoKrig}}, \code{\link{GeoKrigloc}}, \code{\link{GeoSimcond}}, \code{\link{GeoFit}}.}

\author{
Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},
\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl},
\url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},
\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{
\donttest{
library(GeoModels)

###############################################################
### Example of spatial kriging cross-validation
###############################################################

model <- "Gaussian"
set.seed(79)
x <- runif(400, 0, 1)
y <- runif(400, 0, 1)
coords <- cbind(x, y)

corrmodel <- "GenWend"
mean <- 0
sill <- 5
nugget <- 0
scale <- 0.2
smooth <- 0
power2 <- 4

param <- list(
  mean = mean, sill = sill, nugget = nugget,
  scale = scale, smooth = smooth, power2 = power2
)

data <- GeoSim(coordx = coords, corrmodel = corrmodel,
               param = param)$data

fixed <- list(nugget = nugget, smooth = 0, power2 = power2)
start <- list(mean = 0, scale = scale, sill = 1)
I <- Inf
lower <- list(mean = -I, scale = 0, sill = 0)
upper <- list(mean = I, scale = I, sill = I)

fit <- GeoFit(
  data, coordx = coords, corrmodel = corrmodel,
  model = model, likelihood = "Marginal", type = "Pairwise",
  neighb = 3, optimizer = "nlminb", lower = lower,
  upper = upper, start = start, fixed = fixed
)

#a <- GeoCV(fit, K = 100, estimation = TRUE,
#           parallel = TRUE, seed = 123)
#mean(a$rmse)
}
}

\keyword{Composite}
