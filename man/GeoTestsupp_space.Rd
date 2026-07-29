\name{GeoTestsupp_space}
\alias{GeoTestsupp_space}
\title{Parametric Bootstrap Test for Spatial Support}
\description{
Performs a parametric bootstrap likelihood-ratio test for the spatial
\emph{scale}, interpreted as the compact-support radius of the
\code{"GenWend"} correlation model. The null model either constrains the
scale by \eqn{\mathrm{scale} \le h_0}{scale <= h0} or, when \code{h0 = NULL},
represents spatial independence by fixing the support below the minimum
observed inter-site distance. The alternative estimates the scale without the
null upper constraint.
}
\usage{
GeoTestsupp_space(data, coordx,
 start, fixed,
 model = "Gaussian",
 h0 = NULL,
 optimizer = "bobyqa",
 lower = NULL, upper = NULL,
 neighb = 5,
 B = 10000,
 likelihood = NULL,
 type = NULL,
 method = c("cholesky", "TB"),
 L = 10000,
 parallel = TRUE,
 ncores = NULL,
 progress = TRUE,
 seed = NULL)
}
\arguments{
  \item{data}{Numeric vector containing the spatial observations. Missing or
  non-finite values are not allowed.}

  \item{coordx}{Numeric matrix whose rows contain the spatial coordinates.
  Its number of rows must equal \code{length(data)}. Duplicated locations are
  not allowed because the minimum inter-site distance must be positive.}

  \item{start}{Non-empty named list containing parameters to be estimated.
  Together with \code{fixed}, it must contain exactly the parameters returned
  by \code{CorrParam("GenWend")} and \code{NuisParam(model)}.}

  \item{fixed}{Named list containing parameters kept fixed, or \code{NULL}.
  A parameter cannot appear in both \code{start} and \code{fixed}. If
  \code{scale} is supplied in \code{fixed}, it is internally moved to the
  estimated-parameter list for construction of the unrestricted model.}

  \item{model}{Marginal model passed to \code{GeoFit}. Default is
  \code{"Gaussian"}.}

  \item{h0}{Positive finite threshold defining
  \eqn{H_0:\mathrm{scale}\le h_0}{H0: scale <= h0}. If \code{NULL}, the
  support is fixed at \code{max(1e-6, d_min - 1e-6)}, where
  \eqn{d_{min}}{d_min} is the minimum inter-site distance. This removes
  dependence among the sampled locations; for Gaussian random fields it is
  equivalent to independence of the sampled observations. It should not be
  confused with the process-level null hypothesis tested by
  \code{GeoTestIndependence}.}

  \item{optimizer}{Optimization method passed to \code{GeoFit}. Default is
  \code{"bobyqa"}.}

  \item{lower, upper}{Named lists containing bounds for parameters in
  \code{start}. Missing bounds are replaced internally by finite fallback
  bounds. The null upper bound for \code{scale} is always reimposed after any
  numerical widening of the parameter box.}

  \item{neighb}{Number of nearest neighbors used for pairwise composite
  likelihood. It is ignored for full likelihood. When full likelihood is
  selected automatically, it is internally set to \code{NULL}.}

  \item{B}{Positive integer giving the requested number of parametric
  bootstrap replications. Default is \code{10000}. Values below 99 are
  allowed but generate a warning because they give a coarse bootstrap
  \emph{p}-value.}

  \item{likelihood}{Likelihood passed to \code{GeoFit}. Supported values are
  \code{"Full"}, \code{"Marginal"}, and \code{"Conditional"}. Either both
  \code{likelihood} and \code{type} must be supplied, or both must be
  \code{NULL}.}

  \item{type}{Likelihood type passed to \code{GeoFit}. Full likelihood must
  use \code{"Standard"}; marginal or conditional composite likelihood must
  use \code{"Pairwise"}.}

  \item{method}{Simulation method used under the fitted null model.
  \code{"cholesky"} calls \code{GeoSim}; \code{"TB"} calls
  \code{GeoSimapprox}.}

  \item{L}{Number of spectral components used by \code{GeoSimapprox} when
  \code{method = "TB"}. Default is \code{10000}. Ignored for Cholesky
  simulation.}

  \item{parallel}{Logical value indicating whether bootstrap fits should be
  evaluated in parallel when the required packages are available.}

  \item{ncores}{Number of workers. If \code{NULL}, the value is selected
  automatically from the available cores and the number of bootstrap
  replications.}

  \item{progress}{Logical value indicating whether progress information should
  be displayed through \pkg{progressr}.}

  \item{seed}{Optional finite numeric seed. When supplied, the test is
  reproducible and the user's previous random-number-generator state is
  restored when the function exits.}
}
\details{
The correlation model is fixed internally to \code{"GenWend"}. Therefore,
the parameter called \code{scale} is interpreted as the compact-support
radius.

\strong{Likelihood selection.}
When both \code{likelihood} and \code{type} are \code{NULL}, the function uses
full likelihood with \code{type = "Standard"} for at most 2000 observations,
and marginal pairwise composite likelihood for larger datasets. User-supplied
combinations are checked for consistency.

\strong{Null and alternative models.}
For a numeric \code{h0}, the null model estimates \code{scale} subject to
\eqn{\mathrm{scale}\le h_0}{scale <= h0}, while the unrestricted model uses a
data-dependent upper bound strictly larger than \code{h0}. For
\code{h0 = NULL}, the null model fixes the support radius just below the
minimum observed inter-site distance. Since the generalized Wendland
correlation is compactly supported, all distinct observations are then
uncorrelated under the null model.

The lower bound under the alternative must be smaller than \code{h0}, and its
upper bound must be larger than \code{h0}, so that the null parameter space is
nested in the alternative. The function stops if these conditions are not
satisfied.

\strong{Likelihood-ratio statistic.}
Let \eqn{\ell_0}{l0} and \eqn{\ell_1}{l1} denote the fitted log-likelihood or
log-composite-likelihood values under the null and unrestricted models. The
observed statistic is
\deqn{\Lambda_{\mathrm{obs}} = 2(\ell_1-\ell_0).}{Lambda_obs = 2*(l1-l0).}
Tiny negative differences attributable to numerical tolerance are set to zero.
A materially smaller unrestricted likelihood is treated as an optimization
failure rather than automatically converted to a zero statistic.

For pairwise composite likelihood, the function verifies that the null and
unrestricted fits use identical \code{rowidx} and \code{colidx} pair sets. A
likelihood-ratio comparison based on different pair sets is not accepted.

\strong{Adaptive unrestricted-fit fallback.}
The unrestricted model is first fitted once using the primary starting values.
Only if this fit fails or returns a likelihood smaller than the null likelihood
beyond numerical tolerance is a second fit attempted. The fallback starts from
the fitted null parameters and places \code{scale} slightly above \code{h0}.
The better valid unrestricted fit is retained. The same adaptive rule is used
for each bootstrap dataset.

\strong{Parametric bootstrap.}
Bootstrap datasets are simulated from the fitted null model. Each dataset is
refitted under both hypotheses, using the observed null and unrestricted
estimates as starting values. Replications with failed fits, incompatible pair
sets, non-finite likelihoods, or an unresolved ordering
\eqn{\ell_1 < \ell_0}{l1 < l0} are discarded. If \eqn{B_{\mathrm{valid}}}{B_valid}
replications are successful, the bootstrap \emph{p}-value is
\deqn{
\widehat p =
\frac{1 + \sum_{b=1}^{B_{\mathrm{valid}}}
I(\Lambda_b \ge \Lambda_{\mathrm{obs}})}
{B_{\mathrm{valid}} + 1}.
}{
p_hat = (1 + sum(Lambda_b >= Lambda_obs)) / (B_valid + 1).
}
At least \code{min(20, B)} valid replications are required. A warning is issued
when fewer than 80 percent of the requested replications are successful.

\strong{Early exit.}
If the unrestricted estimate itself satisfies the null constraint, the
likelihood-ratio statistic is zero and the function returns \code{pvalue = 1}
without running the bootstrap.

\strong{Parallel computation.}
For \code{method = "TB"}, simulation may use parallel computation inside
\code{GeoSimapprox}. Bootstrap refits are performed in a separate phase. In
the parallel bootstrap, each simulated dataset is written to a separate RDS
file so that each worker reads only its own replicate. The user's
\pkg{future} plan, \pkg{progressr} handlers, RNG state, and temporary files are
restored or removed on exit.
}
\value{
An invisible list. When the bootstrap is run, it contains:

  \item{d_min}{Minimum inter-site distance.}
  \item{h0}{Threshold actually used under the null hypothesis.}
  \item{lambda_obs}{Observed likelihood-ratio statistic.}
  \item{pvalue}{Parametric-bootstrap \emph{p}-value.}
  \item{seed}{Seed supplied by the user, or \code{NULL}.}
  \item{method}{Simulation method used.}
  \item{L}{Number of turning-bands components when \code{method = "TB"},
  otherwise \code{NULL}.}
  \item{B_requested}{Requested number of bootstrap replications.}
  \item{bootstrap_successful}{Number of valid bootstrap replications.}
  \item{bootstrap_failed}{Number of failed or discarded bootstrap
  replications.}
  \item{bootstrap_success_rate}{Fraction of requested replications that were
  valid.}
  \item{fallback_observed}{Logical value indicating whether the adaptive
  unrestricted-fit fallback was used for the observed data.}
  \item{fallback_bootstrap}{Number of bootstrap replications in which the
  adaptive unrestricted-fit fallback was attempted.}
  \item{B_rep}{Vector of valid bootstrap likelihood-ratio statistics.}
  \item{fit_H0}{Fitted \code{GeoFit} object under the null hypothesis.}
  \item{fit_H1}{Fitted \code{GeoFit} object under the unrestricted
  alternative.}

For an early exit, \code{lambda_obs = 0}, \code{pvalue = 1},
\code{bootstrap_successful = 0}, \code{bootstrap_success_rate = NA}, and
\code{B_rep} is empty. In the current implementation, \code{L} and
\code{bootstrap_failed} are omitted from this early-return object.
}
\note{
The test concerns the compact-support radius of the generalized Wendland
correlation model and is not a generic range-parameter test for arbitrary
correlation functions. For final inference, a substantially larger value of
\code{B} than that used in quick examples is recommended.
}
\seealso{
\code{\link{GeoFit}}, \code{\link{GeoSim}},
\code{\link{GeoSimapprox}}, \code{\link{GeoTestIndependence}},
\code{\link{GeoTestIsotropy}}, \code{\link{GeoTests}}
}
\examples{
\dontrun{
library(GeoModels)

set.seed(123)
n <- 300
coords <- cbind(runif(n), runif(n))

param <- list(
  mean = 0,
  nugget = 0,
  sill = 1,
  scale = 0.20,
  smooth = 0.5,
  power2 = 4
)

z <- GeoSim(
  coordx = coords,
  corrmodel = "GenWend",
  model = "Gaussian",
  param = param
)$data

start <- list(
  mean = 0,
  sill = 1,
  scale = 0.15
)
fixed <- list(
  nugget = 0,
  smooth = 0.5,
  power2 = 4
)

lower <- list(mean = -Inf, sill = 1e-6, scale = 1e-6)
upper <- list(mean = Inf, sill = Inf, scale = 1)

## Test H0: scale <= 0.10 against scale > 0.10.
ans <- GeoTestsupp_space(
  data = z,
  coordx = coords,
  start = start,
  fixed = fixed,
  h0 = 0.10,
  lower = lower,
  upper = upper,
  B = 199,
  likelihood = "Full",
  type = "Standard",
  parallel = FALSE,
  seed = 321
)

ans$pvalue
ans$lambda_obs
ans$bootstrap_success_rate

## Independence test: support fixed below the minimum inter-site distance.
ans_ind <- GeoTestsupp_space(
  data = z,
  coordx = coords,
  start = start,
  fixed = fixed,
  h0 = NULL,
  lower = lower,
  upper = upper,
  B = 199,
  likelihood = "Full",
  type = "Standard",
  parallel = FALSE,
  seed = 654
)
}
}
\author{
Moreno Bevilacqua
\email{moreno.bevilacqua89@gmail.com}
\url{https://sites.google.com/view/moreno-bevilacqua/home}
\cr
Victor Morales Onate
\email{victor.morales@uv.cl}
\url{https://sites.google.com/site/moralesonatevictor/}
\cr
Christian Caamano-Carrillo
\email{chcaaman@ubiobio.cl}
\url{https://www.researchgate.net/profile/Christian-Caamano}
}
\keyword{Bootstrap}
\keyword{Composite}
\keyword{HypothesisTest}
\keyword{Spatial}
