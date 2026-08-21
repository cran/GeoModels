\name{GeoTestIsotropy}
\alias{GeoTestIsotropy}
\title{Parametric Bootstrap Tests for Spatial Isotropy}
\description{
Performs two parametric-bootstrap tests of spatial isotropy by comparing an
isotropic model under the null hypothesis with a geometrically anisotropic
model under the alternative. The first test (MWCLT) is based on the estimated
anisotropy ratio. The second test (WCLRT) is based on a likelihood-ratio or
composite-likelihood-ratio statistic. Their null distributions are estimated
by simulation from the fitted isotropic model.
}

\usage{
GeoTestIsotropy(data, coordx,
 start, fixed,
 optimizer = "bobyqa",
 model = "Gaussian",
 corrmodel = "Matern",
 lower = NULL, upper = NULL,
 B = 1000,
 likelihood = NULL,
 type = NULL,
 copula = NULL,
 neighb = 5,
 method = c("cholesky", "TB"),
 L = 10000,
 n = 1,
 angle_starts = c(pi/4, 3*pi/4),
 batch_size = NULL,
 max_attempts = NULL,
 parallel = FALSE,
 ncores = 6,
 progress = TRUE,
 seed = NULL)
}

\arguments{
 \item{data}{Numeric vector of observations.}
 \item{coordx}{Numeric matrix with exactly two columns containing the spatial
 coordinates. The number of rows must equal \code{length(data)}.}
 \item{start}{Named list containing the covariance and marginal parameters to
 be estimated. Parameters \code{angle} and \code{ratio} must not be included,
 because they are handled internally by the isotropy test.}
 \item{fixed}{Named list containing the covariance and marginal parameters kept
 fixed. Together, \code{start} and \code{fixed} must contain exactly the
 parameters required by \code{CorrParam(corrmodel)} and
 \code{NuisParam(model)}.}
 \item{optimizer}{Optimization algorithm passed to \code{GeoFit}.}
 \item{model}{Marginal random-field model.}
 \item{corrmodel}{Spatial correlation model.}
 \item{lower, upper}{Named lists of lower and upper bounds for the parameters
 in \code{start}. Missing bounds are replaced by \code{-Inf} and \code{Inf},
 respectively.}
 \item{B}{Number of valid parametric-bootstrap replications used to calibrate
 each test. Failed fitting attempts are replaced by newly simulated samples.}
 \item{likelihood}{Likelihood type passed to \code{GeoFit}. If both
 \code{likelihood} and \code{type} are \code{NULL}, full likelihood is used for
 at most 2000 observations and marginal pairwise composite likelihood is used
 otherwise.}
 \item{type}{Likelihood type passed to \code{GeoFit}; typically
 \code{"Standard"} for full likelihood or \code{"Pairwise"} for composite
 likelihood. It must be supplied together with \code{likelihood}.}
 \item{copula}{Optional copula specification. With a supported purely spatial
 correlation model, both \code{method="cholesky"} and \code{method="TB"} are
 available.}
 \item{neighb}{Number of nearest neighbours used when
 \code{type = "Pairwise"}. Ignored otherwise.}
 \item{method}{Simulation method under the null hypothesis. Use
 \code{"cholesky"} for exact simulation or \code{"TB"} for approximate
 turning-bands simulation. When \code{copula} is supplied, the corresponding
 latent engine is selected internally by \code{GeoSimCopula}.}
 \item{L}{Number of spectral components used when \code{method = "TB"}.}
 \item{n}{Common number of trials/successes used by models that require an
 integer \code{n}, such as Binomial-type models. It can be a scalar or, where
 supported, a vector with one value per observation. For \code{"Bernoulli"}
 and \code{"Binary"}, \code{n} is internally fixed to 1.}
 \item{angle_starts}{Finite numeric vector of starting angles used for the
 anisotropic fit. All supplied angles are tried (modulo \eqn{\pi}) and the fit
 with the largest objective value is retained. The default uses two starts,
 \eqn{\pi/4} and \eqn{3\pi/4}, for both the observed fit and every bootstrap
 refit.}
 \item{batch_size}{Number of bootstrap datasets simulated and processed at one
 time. If \code{NULL}, a modest batch size is selected automatically. Batching
 limits memory use for large spatial datasets.}
 \item{max_attempts}{Maximum number of simulated bootstrap datasets allowed in
 order to collect exactly \code{B} valid replications. If \code{NULL}, the
 default is \code{max(B + 20, ceiling(1.5 * B))}. The function stops if this
 limit is reached before \code{B} valid replications are obtained.}
 \item{parallel}{Logical; if \code{TRUE}, bootstrap refits are distributed
 across multiple R sessions. With \code{method="TB"}, bootstrap simulation
 may also use multiple workers when supported by the simulation routine.}
 \item{ncores}{Positive integer or \code{NULL}. With \code{parallel=TRUE}, an
 explicit integer requests that many workers, capped by detected cores and the
 current number of jobs. Set \code{ncores=NULL} for automatic selection.}
 \item{progress}{Logical; if \code{TRUE}, display the same
 \pkg{progressr} text progress bar used by other parallel \pkg{GeoModels}
 functions. A single bar is shown for the whole bootstrap and counts valid
 bootstrap replications in both sequential and parallel execution.}
 \item{seed}{Optional finite numeric seed. The previous random-number-generator
 state is restored when the function exits.}
}

\details{
The null and alternative hypotheses are
\deqn{H_0: \mathrm{ratio}=1,}
and
\deqn{H_1: \mathrm{ratio}>1,}
with the anisotropy angle estimated under the alternative. Under \eqn{H_0},
the angle is not identifiable and the ratio is on the boundary of its
parameter space. Consequently, a standard chi-squared reference distribution
is not assumed; both tests are calibrated by parametric bootstrap.

The MWCLT statistic is the unconstrained anisotropy-ratio estimate
\eqn{\widehat\tau}. If \eqn{\widehat\tau_b^*} denotes the corresponding
estimate in bootstrap replication \eqn{b}, the upper-tail bootstrap p-value is
\deqn{\widehat p_\tau = \frac{1 + \sum_{b=1}^B I(\widehat\tau_b^* \geq \widehat\tau)}{B+1}.}

The WCLRT statistic is
\deqn{T_{\mathrm{obs}}=\max\{0,2(\ell_1-\ell_0)\},}
where \eqn{\ell_0} and \eqn{\ell_1} are the maximized likelihoods or composite
likelihoods under the isotropic and anisotropic models. Its bootstrap p-value
is
\deqn{\widehat p_T = \frac{1 + \sum_{b=1}^B I(T_b \geq T_{\mathrm{obs}})}{B+1}.}
The returned components \code{statistic} and \code{pvalue} remain aliases for
this WCLRT statistic and p-value for backward compatibility.

For every observed or bootstrap anisotropic fit, all angles in
\code{angle_starts} are used as starting values and the valid fit with the
largest objective value is retained. During bootstrap refitting, the fitted
observed anisotropy angle is also included as a starting angle. If the best
numerical anisotropic fit has an objective value below the isotropic fit by
more than numerical tolerance, the boundary solution is selected, giving
ratio 1 and WCLRT statistic 0.

For pairwise composite likelihood, the isotropic and anisotropic fits must use
the same retained pairs. The function checks this condition before computing
the WCLRT statistic.

Bootstrap datasets are generated from the fitted isotropic model in batches.
Failed bootstrap fits do not enter the empirical null distribution: new
samples are generated until exactly \code{B} valid replications have been
collected, or until \code{max_attempts} is reached. This makes the denominator
of both bootstrap p-values exactly \code{B+1}.

The argument \code{n} is propagated consistently to \code{GeoFit},
\code{GeoSim}, \code{GeoSimapprox}, and \code{GeoSimCopula}. Direct latent
Gaussian simulation of some non-Gaussian models imposes additional structural
constraints. In particular, direct \code{"Gamma"} simulation requires the
fitted \code{shape} to be a positive integer because it is the number of
squared latent Gaussian fields. For the direct Gamma model used in the latent
construction, it is therefore recommended to keep \code{shape} fixed at a
positive integer during the isotropy test. This integer restriction does not
apply merely because the marginal distribution is Gamma when a supported
copula construction is used.
}

\value{
An invisible list containing, among others:
\item{statistic, pvalue}{Backward-compatible aliases for the observed WCLRT
 statistic and its bootstrap p-value.}
\item{statistic_lrt, pvalue_lrt}{Observed WCLRT statistic and p-value.}
\item{statistic_ratio, pvalue_ratio}{Observed MWCLT anisotropy-ratio statistic
 and p-value.}
\item{pvalues}{Named vector containing the MWCLT and WCLRT p-values.}
\item{ratio_hat, angle_hat}{Estimated anisotropy ratio and angle. The angle is
 \code{NA} when the boundary solution is selected.}
\item{parallel, ncores}{Whether parallel bootstrap refitting was used and the
number of workers selected by the function.}
\item{fit_H0, fit_H1}{Fitted isotropic and best anisotropic \code{GeoFit}
 objects for the observed data.}
\item{B_rep, B_rep_lrt}{The \code{B} valid bootstrap WCLRT statistics.}
\item{B_rep_ratio}{The \code{B} valid bootstrap anisotropy-ratio estimates.}
\item{B_rep_angle}{Bootstrap angle estimates; boundary solutions are
 \code{NA}.}
\item{bootstrap_B}{Requested and obtained number of valid bootstrap
 replications.}
\item{bootstrap_total_attempts}{Total number of simulated datasets needed to
 obtain \code{B} valid replications.}
\item{bootstrap_failed_attempts}{Number of failed attempts that were replaced.}
\item{bootstrap_failures}{Data frame describing failed attempts.}
\item{bootstrap_boundary_solutions}{Number of valid bootstrap replications for
 which the boundary solution was selected.}
\item{observed_multistart_angles}{Starting angles used for the observed H1 fit.}
\item{batch_size, max_attempts}{Bootstrap batching and retry settings used by
 the function.}
}

\examples{
\dontrun{
## Isotropy tests for a Gaussian random field simulated under H0
set.seed(123)
nsite <- 500
coords <- cbind(runif(nsite), runif(nsite))

param <- list(
  mean = 0,
  nugget = 0,
  sill = 1,
  scale = 0.20,
  smooth = 0.5
)

z <- GeoSim(
  coordx = coords,
  corrmodel = "Matern",
  model = "Gaussian",
  param = param
)$data

start <- list(
  mean = 0,
  sill = 1,
  scale = 0.20,
  smooth = 0.5
)
fixed <- list(nugget = 0)

lower <- list(mean = -5, sill = 0.01, scale = 0.01, smooth = 0.05)
upper <- list(mean =  5, sill = 5,    scale = 1,    smooth = 2.5)

iso_test <- GeoTestIsotropy(
  data = z,
  coordx = coords,
  start = start,
  fixed = fixed,
  model = "Gaussian",
  corrmodel = "Matern",
  lower = lower,
  upper = upper,
  B = 99,
  likelihood = "Full",
  type = "Standard",
  method = "cholesky",
  seed = 321
)

iso_test$pvalue_ratio  # MWCLT
iso_test$pvalue_lrt    # WCLRT
}
}

\seealso{
\code{\link{GeoFit}}, \code{\link{GeoSim}},
\code{\link{GeoSimapprox}}, \code{\link{GeoVariogramDir}}
}

\author{
Moreno Bevilacqua \email{moreno.bevilacqua89@gmail.com}
\url{https://sites.google.com/view/moreno-bevilacqua/home}\cr
Victor Morales Onate \email{victor.morales@uv.cl}
\url{https://sites.google.com/site/moralesonatevictor/}\cr
Christian Caamano-Carrillo \email{chcaaman@ubiobio.cl}
\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\keyword{Bootstrap}
\keyword{Composite}
\keyword{HypothesisTest}
\keyword{Spatial}
