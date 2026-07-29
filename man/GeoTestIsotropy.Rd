\name{GeoTestIsotropy}
\alias{GeoTestIsotropy}
\title{Parametric Bootstrap Test for Spatial Isotropy}
\description{
Performs a parametric-bootstrap test of spatial isotropy by comparing an
isotropic model under the null hypothesis with a geometrically anisotropic
model under the alternative. The test statistic is a likelihood-ratio or
composite-likelihood-ratio statistic, and its null distribution is estimated
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
 parallel = TRUE,
 ncores = NULL,
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
 \item{B}{Number of parametric-bootstrap replications.}
 \item{likelihood}{Likelihood type passed to \code{GeoFit}. If both
 \code{likelihood} and \code{type} are \code{NULL}, full likelihood is used for
 at most 2000 observations and marginal pairwise composite likelihood is used
 otherwise.}
 \item{type}{Likelihood type passed to \code{GeoFit}; typically
 \code{"Standard"} for full likelihood or \code{"Pairwise"} for composite
 likelihood. It must be supplied together with \code{likelihood}.}
 \item{copula}{Optional copula specification. When a copula is supplied, only
 \code{method = "cholesky"} is currently supported.}
 \item{neighb}{Number of nearest neighbours used when
 \code{type = "Pairwise"}. Ignored otherwise.}
 \item{method}{Simulation method under the null hypothesis. Use
 \code{"cholesky"} for exact simulation through \code{GeoSim}, or \code{"TB"}
 for approximate turning-bands simulation through \code{GeoSimapprox}.}
 \item{L}{Number of spectral components used when \code{method = "TB"}.}
 \item{parallel}{Logical; if \code{TRUE}, bootstrap refits may be distributed
 across multiple R sessions.}
 \item{ncores}{Number of parallel workers. If \code{NULL}, it is selected
 automatically.}
 \item{progress}{Logical; if \code{TRUE}, display progress information through
 \pkg{progressr}.}
 \item{seed}{Optional finite numeric seed. The previous random-number-generator
 state is restored when the function exits.}
}

\details{
The null and alternative hypotheses are
\deqn{H_0: \mathrm{ratio}=1,}
and
\deqn{H_1: \mathrm{ratio}>1,}
with the anisotropy angle estimated under the alternative. Under \eqn{H_0},
the angle is not identifiable. Consequently, the usual chi-squared reference
distribution for the likelihood-ratio statistic is not used. Instead, the
null distribution is estimated by parametric bootstrap.

The observed statistic is
\deqn{T_{\mathrm{obs}}=\max\{0,2(\ell_1-\ell_0)\},}
where \eqn{\ell_0} and \eqn{\ell_1} are the maximized likelihoods or composite
likelihoods under the isotropic and anisotropic models. Bootstrap datasets are
simulated from the fitted isotropic model, and both models are refitted to each
simulated dataset. The bootstrap p-value is
\deqn{\frac{1+\sum_b I(T_b \geq T_{\mathrm{obs}})}{1+B_{\mathrm{valid}}}.}

The anisotropic model is normally fitted once, starting from
\code{angle = 0} and \code{ratio = 1.2}. An adaptive multi-start fallback is
used only if that fit fails or gives a smaller likelihood than the isotropic
fit. In that case, additional starting angles are tried and the fit with the
largest likelihood is retained. The same estimation rule is used for the
observed data and for every bootstrap replication.

For pairwise composite likelihood, the isotropic and anisotropic fits must use
the same retained pairs. The function checks this condition before computing
the likelihood-ratio statistic.

The estimates obtained from the observed fits are used as starting values in
the bootstrap refits. Failed bootstrap replications are removed from the
empirical null distribution and are reported in the returned object. A warning
is issued when fewer than 80 percent of the replications are successful or
when fewer than 20 valid bootstrap statistics are available.
}

\value{
An invisible list with components:
\item{statistic}{Observed likelihood-ratio statistic.}
\item{pvalue}{Parametric-bootstrap p-value.}
\item{ratio_hat}{Estimated anisotropy ratio under the alternative.}
\item{angle_hat}{Estimated anisotropy angle under the alternative.}
\item{corrmodel}{Correlation model used in the test.}
\item{model}{Marginal model used in the test.}
\item{likelihood}{Likelihood type used in the fits.}
\item{type}{Likelihood construction used in the fits.}
\item{method}{Simulation method used under the null hypothesis.}
\item{L}{Number of turning-bands components, or \code{NULL} when not used.}
\item{seed}{Seed supplied by the user, or \code{NULL}.}
\item{fit_H0}{Fitted isotropic \code{GeoFit} object.}
\item{fit_H1}{Fitted anisotropic \code{GeoFit} object.}
\item{observed_fallback_used}{Logical indicating whether the adaptive
 multi-start fallback was used for the observed anisotropic fit.}
\item{observed_H1_attempts}{Number of anisotropic fitting attempts for the
 observed data.}
\item{B_rep}{Vector of valid bootstrap likelihood-ratio statistics.}
\item{bootstrap_B}{Requested number of bootstrap replications.}
\item{bootstrap_successful}{Number of valid bootstrap replications.}
\item{bootstrap_success_rate}{Proportion of valid bootstrap replications.}
\item{bootstrap_fallbacks}{Number of bootstrap replications in which the
 adaptive multi-start fallback was activated.}
\item{bootstrap_attempts}{Number of anisotropic fitting attempts used in each
 bootstrap replication.}
\item{bootstrap_failures}{Data frame containing the indices and error messages
 of failed bootstrap replications.}
}

\examples{
\dontrun{
## Isotropy test for a Gaussian random field simulated under H0
set.seed(123)
n <- 500
coords <- cbind(runif(n), runif(n))

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

lower <- list(
  mean = -5,
  sill = 0.01,
  scale = 0.01,
  smooth = 0.05
)
upper <- list(
  mean = 5,
  sill = 5,
  scale = 1,
  smooth = 2.5
)

#iso_test <- GeoTestIsotropy(
#  data = z,
#  coordx = coords,
#  start = start,
#  fixed = fixed,
#  optimizer = "bobyqa",
#  model = "Gaussian",
#  corrmodel = "Matern",
#  lower = lower,
#  upper = upper,
#  B = 99,
#  likelihood = "Full",
#  type = "Standard",
#  method = "cholesky",
#  parallel =TRUE,
#  seed = 321)

#iso_test$statistic
#iso_test$pvalue
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
