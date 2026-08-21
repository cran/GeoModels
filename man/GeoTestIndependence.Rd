\name{GeoTestIndependence}
\alias{GeoTestIndependence}
\title{Parametric bootstrap test for absence of spatial dependence}
\description{
Performs a parametric bootstrap likelihood-ratio test of spatial independence
against a spatially dependent random-field model. Under the null hypothesis,
the observations at distinct sites are fitted using the marginal independence
likelihood. Under the alternative hypothesis, a full spatial likelihood is
fitted with the nugget proportion estimated below its independence boundary.
The bootstrap is required because the null hypothesis corresponds to a boundary
point and the spatial correlation parameters are not identified under the null.
}
\usage{
GeoTestIndependence(data, coordx,
 start, fixed = NULL,
 corrmodel = "Matern",
 model = "Gaussian",
 optimizer = "bobyqa",
 optimizer_ind = "Nelder-Mead",
 lower = NULL, upper = NULL,
 X = NULL, n = 1,
 distance = "Eucl", radius = 1,
 anisopars = NULL,
 est.aniso = c(FALSE, FALSE),
 sparse = FALSE,
 B = 1000,
 parallel = FALSE,
 ncores = 6,
 progress = TRUE,
 seed = NULL)
}
\arguments{
 \item{data}{Numeric vector of observations.}

 \item{coordx}{Matrix or object coercible to a matrix containing the spatial
 coordinates. The number of rows must equal \code{length(data)}. Duplicate
 locations are not allowed.}

 \item{start}{Non-empty named list containing the parameters to be estimated
 under the spatial alternative. It must include, either directly or through
 \code{fixed}, all parameters required by \code{model} and \code{corrmodel}.
 The parameter \code{nugget} is always estimated under the alternative; if it
 is supplied in \code{fixed}, it is automatically moved to \code{start}.}

 \item{fixed}{Optional named list containing parameters kept fixed under the
 spatial alternative. Parameters cannot occur in both \code{start} and
 \code{fixed}.}

 \item{corrmodel}{Correlation model used under the spatial alternative. The
 default is \code{"Matern"}.}

 \item{model}{Marginal model for the observations. The default is
 \code{"Gaussian"}.}

 \item{optimizer}{Optimization method passed to \code{GeoFit} for the full
 spatial fit under \eqn{H_1}. The default is \code{"bobyqa"}.}

 \item{optimizer_ind}{Optimization method used for the marginal independence
 fit under \eqn{H_0}. The default is \code{"Nelder-Mead"}. If the requested
 method fails, the function also tries \code{"Nelder-Mead"} as a fallback and
 retains the successful fit with the largest finite log-likelihood.}

 \item{lower, upper}{Optional named lists of lower and upper bounds for the
 parameters estimated under the spatial alternative. Names must refer only to
 parameters in \code{start}, after the possible promotion of \code{nugget}
 from \code{fixed}. The upper bound for \code{nugget} must be one, or omitted,
 so that the independence model belongs to the closure of the alternative.}

 \item{X}{Optional design matrix for the mean or marginal regression model.
 Its number of rows must equal \code{length(data)}.}

 \item{n}{Model-specific parameter passed to \code{GeoFit} and \code{GeoSim}.
 The default is \code{1}.}

 \item{distance}{Distance metric passed to \code{GeoFit}, \code{GeoDistances},
 and \code{GeoSim}. The default is \code{"Eucl"}.}

 \item{radius}{Positive radius used for spherical or geodesic distances. The
 default is \code{1}.}

 \item{anisopars}{Optional list of anisotropy parameters passed to the spatial
 fit under the alternative.}

 \item{est.aniso}{Logical vector of length two indicating whether the
 anisotropy angle and ratio are estimated under the alternative. The default is
 \code{c(FALSE, FALSE)}.}

 \item{sparse}{Logical value passed to the full-likelihood spatial fit. The
 default is \code{FALSE}.}

 \item{B}{Number of parametric bootstrap replications. The default is
 \code{1000}. Values below \code{99} produce a coarse bootstrap \emph{p}-value
 and generate a warning.}

 \item{parallel}{Logical; default \code{FALSE}. Set \code{TRUE} to evaluate bootstrap fits in parallel when the required packages are available.}

 \item{ncores}{Positive integer or \code{NULL}; default \code{6}. With \code{parallel=TRUE}, an explicit integer requests that many workers, capped by detected cores and the number of bootstrap replications. Set \code{ncores=NULL} for automatic selection, capped at six workers and normally leaving one detected core free.}

 \item{progress}{Logical value indicating whether progress bars should be
 displayed through \pkg{progressr}.}

 \item{seed}{Optional finite numeric seed. When supplied, the previous global
 random-number-generator state is restored when the function exits.}
}
\details{
The test compares
\deqn{H_0: \text{absence of spatial dependence}}
with
\deqn{H_1: \text{a spatially dependent random-field model}.}

Under \eqn{H_0}, the function fits
\code{GeoFit(..., likelihood = "Marginal", type = "Independence")}. Only
marginal and regression parameters are estimated. Spatial correlation
parameters and the nugget proportion are absent from the null fit because they
are not identified under independence.
The null fit uses \code{optimizer_ind}; if necessary, a
\code{"Nelder-Mead"} fallback is attempted. The optimizer associated with the
selected null fit is returned in \code{optimizer_independence}.

Under \eqn{H_1}, the function fits
\code{GeoFit(..., likelihood = "Full", type = "Standard")} using
\code{corrmodel}. The nugget parameter is estimated on
\eqn{[0,1)}. Independence is represented by the boundary
\eqn{\texttt{nugget}=1}, where the spatial correlation contribution vanishes.
The upper optimization bound is internally set to a value immediately below
one because \code{GeoFit} requires \code{nugget < 1}.

The observed likelihood-ratio statistic is
\deqn{\Lambda_{obs}=2\{\ell_1-\ell_0\},}
where \eqn{\ell_0} is the maximized independence log-likelihood and
\eqn{\ell_1} is the maximized full spatial log-likelihood. If the unrestricted
spatial fit reaches the independence boundary, the statistic is set to zero
and the function returns \code{pvalue = 1} without running the bootstrap.

The usual chi-squared likelihood-ratio approximation is not used because the
null value is on the boundary and the spatial correlation parameters are not
identified under \eqn{H_0}. Instead, the null distribution is approximated by
parametric bootstrap. The bootstrap \emph{p}-value is
\deqn{
\widehat p =
\frac{1+\sum_{b=1}^{B_{valid}}
I(\Lambda_b^*\geq\Lambda_{obs})}
{B_{valid}+1}.
}

Bootstrap datasets are generated exactly independently at the observed sites.
Because \code{GeoSim} does not accept \code{nugget = 1}, the implementation
uses an auxiliary \code{GenWend} correlation model with compact support equal
to one half of the minimum inter-site distance. Consequently, all off-diagonal
correlations at the observed locations are exactly zero. This auxiliary model
is used only as a simulation device; the fitted null model remains the marginal
independence model.

The observed parameter estimates are used as starting values in the bootstrap
refits. If the first spatial fit fails to attain at least the independence
log-likelihood, adaptive fallback fits are attempted using alternative initial
values for \code{nugget}. When the best spatial fit is numerically equivalent
to the null fit, the bootstrap statistic is treated as a boundary value equal
to zero.

When \code{parallel = TRUE}, bootstrap fits are evaluated with
\pkg{future.apply}. Each simulated dataset is written temporarily to a separate
RDS file so that workers do not receive the complete simulation object. The
user's previous \pkg{future} plan, \pkg{progressr} handlers, and random-number-
generator state are restored on exit.

The current implementation is restricted to univariate models with a single
nugget proportion named \code{nugget}. It compares the independence likelihood
only with the full spatial likelihood. It does not directly compare the
independence likelihood with a pairwise composite likelihood, because the two
criteria are not on the same scale.
}
\value{
An invisible list with components:

\item{lambda_obs}{Observed likelihood-ratio statistic.}

\item{pvalue}{Parametric bootstrap \emph{p}-value.}

\item{nugget_hat}{Estimated nugget proportion under the spatial alternative.}

\item{boundary_observed}{Logical value indicating whether the unrestricted
spatial fit was classified as an independence-boundary solution.}

\item{optimizer_independence}{Optimizer used by the selected independence fit
under \eqn{H_0}. This can differ from \code{optimizer_ind} when the
\code{"Nelder-Mead"} fallback is selected.}

\item{optimizer_spatial}{Optimizer used for the full spatial fit under
\eqn{H_1}.}

\item{seed}{Seed supplied by the user, or \code{NULL}.}

\item{B_requested}{Requested number of bootstrap replications.}

\item{bootstrap_successful}{Number of valid bootstrap replications.}

\item{bootstrap_failed}{Number of failed bootstrap replications.}

\item{bootstrap_success_rate}{Fraction of successful bootstrap replications.
This is \code{NA} when the observed statistic is zero and the bootstrap is not
run.}

\item{fallback_observed}{Logical value indicating whether adaptive starting
values were required for the observed spatial fit.}

\item{fallback_bootstrap}{Number of bootstrap replications in which adaptive
starting values were used for the spatial fit.}

\item{boundary_bootstrap}{Number of valid bootstrap replications classified as
independence-boundary solutions.}

\item{simulation_support}{Compact-support value used by the auxiliary
\code{GenWend} model to simulate independent observations. This component is
absent when the observed statistic is zero and the bootstrap is skipped.}

\item{B_rep}{Vector of valid bootstrap likelihood-ratio statistics.}

\item{fit_H0}{\code{GeoFit} object fitted under spatial independence.}

\item{fit_H1}{\code{GeoFit} object fitted under the spatial alternative.}
}
\note{
For Gaussian random fields, the null hypothesis is equivalent to mutual
independence of the observations at distinct spatial locations. For
non-Gaussian models, the interpretation as independence requires that zero
latent spatial dependence imply independence in the corresponding observed
random-field construction.

A small bootstrap \emph{p}-value indicates evidence of spatial dependence. For
final analyses, a value such as \code{B = 999} or larger is recommended.
}
\seealso{
\code{\link{GeoFit}},
\code{\link{GeoSim}},
\code{\link{GeoDistances}},
\code{\link{GeoTestsupp_space}},
\code{\link{GeoTestIsotropy}}
}
\examples{
\dontrun{
library(GeoModels)

set.seed(3)
N <- 300
coords <- cbind(runif(N), runif(N))

X <- cbind(1, runif(N))

param <- list(
  mean = 0.2,
  mean1 = -0.5,
  sill = 1,
  nugget = 0,
  scale = 0.2 / 3,
  smooth = 0.5
)

data <- GeoSim(
  coordx = coords,
  corrmodel = "Matern",
  model = "Gaussian",
  param = param,
  X = X
)$data

start <- list(
  mean = 0,
  mean1 = 0,
  sill = 1,
  nugget = 0.1,
  scale = 0.1
)

fixed <- list(smooth = 0.5)

lower <- list(
  mean = -Inf,
  mean1 = -Inf,
  sill = 1e-6,
  nugget = 0,
  scale = 1e-6
)

upper <- list(
  mean = Inf,
  mean1 = Inf,
  sill = Inf,
  nugget = 1,
  scale = 1
)

ans <- GeoTestIndependence(
  data = data,
  coordx = coords,
  X = X,
  start = start,
  fixed = fixed,
  corrmodel = "Matern",
  model = "Gaussian",
  optimizer = "bobyqa",
  optimizer_ind = "Nelder-Mead",
  lower = lower,
  upper = upper,
  B = 199,
  parallel = TRUE,
  seed = 123
)

ans$pvalue
ans$lambda_obs
ans$nugget_hat
ans$bootstrap_success_rate
}
}
\author{
Moreno Bevilacqua
\email{moreno.bevilacqua89@gmail.com}
\url{https://sites.google.com/view/moreno-bevilacqua/home}
}
\keyword{Bootstrap}
\keyword{Likelihood}
\keyword{HypothesisTest}
\keyword{Spatial}
