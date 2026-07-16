\name{GeoTestsupp_space}
\alias{GeoTestsupp_space}
\title{Parametric bootstrap test for spatial support (scale)}
\description{
Performs a parametric bootstrap likelihood-ratio test on the spatial \emph{scale} parameter using the compactly supported \code{GenWend} correlation model.
Two hypotheses are compared: \strong{H0} imposes a maximum admissible scale (or independence), while \strong{H1} estimates scale freely.
The test statistic is the likelihood ratio and the \emph{p}-value is obtained by parametric bootstrap under \strong{H0}.
}
\usage{
GeoTestsupp_space(data, coordx,
 start, fixed,
 model = "Gaussian",
 h0 = NULL,
 optimizer = "bobyqa",
 lower = NULL, upper = NULL,
 neighb = 5,
 B = 1000,
 likelihood = NULL,
 type = NULL,
 method = c("cholesky", "TB"),
 L = 1000,
 parallel = TRUE,
 ncores = NULL,
 progress = TRUE,
 seed = NULL)
}
\arguments{
 \item{data}{Numeric vector of observations of length \code{n}.}
 \item{coordx}{\code{n x d} matrix with spatial coordinates.}
 \item{start}{\emph{Named} list of parameters \strong{to be estimated}. Names must match the required correlation and nuisance parameters not supplied in \code{fixed}.}
 \item{fixed}{\emph{Named} list of parameters \strong{kept fixed}. Together with \code{start}, it must exactly cover the required set of correlation and nuisance parameters.}
 \item{model}{Marginal model for the data. Default is \code{"Gaussian"}.}
 \item{h0}{Positive threshold for the \emph{scale} parameter under \strong{H0}. If \code{NULL}, \strong{H0} represents spatial independence with \code{scale = h0_val} set to the adjusted minimum inter-point distance.}
 \item{optimizer}{Optimization method passed to \code{GeoFit}, for example \code{"bobyqa"}.}
 \item{lower, upper}{\emph{Named} lists with bounds for parameters in \code{start} only. \code{NULL} means no extra bounds beyond the internal defaults. Bounds are internally widened if too narrow.}
 \item{neighb}{Number of neighbors for composite likelihood when applicable. Ignored when the full likelihood is used.}
 \item{B}{Number of parametric bootstrap replications. Default is \code{1000}.}
 \item{likelihood}{Type of likelihood to pass to \code{GeoFit}, for example \code{"Marginal"}. Either specify both \code{likelihood} and \code{type}, or leave both \code{NULL} for automatic selection.}
 \item{type}{Likelihood type to pass to \code{GeoFit}, for example \code{"Pairwise"}. Either specify both \code{likelihood} and \code{type}, or leave both \code{NULL} for automatic selection.}
 \item{method}{Simulation method for bootstrap samples under \strong{H0}. Use \code{"cholesky"} for \code{GeoSim} or \code{"TB"} for turning-bands simulation through \code{GeoSimapprox}.}
 \item{L}{Number of spectral components used by \code{GeoSimapprox} when \code{method = "TB"}. Default is \code{1000}. Ignored when \code{method = "cholesky"}.}
 \item{parallel}{If \code{TRUE}, parallel computation is used when available. For \code{method = "TB"}, \code{GeoSimapprox} may parallelize the simulation of the \code{B} bootstrap datasets; the subsequent bootstrap fits may also be parallelized in a separate phase.}
 \item{ncores}{Number of cores for parallelization. If \code{NULL}, it is selected automatically up to \code{B} or \code{parallel::detectCores()-1}.}
 \item{progress}{If \code{TRUE}, progress bars are shown via \pkg{progressr}.}
 \item{seed}{Optional integer seed used to initialize the random number generator before the parametric bootstrap simulation. If \code{NULL}, the current random state is left unchanged. When supplied, the previous \code{.Random.seed} is restored on exit.}
}
\details{
The function fixes the correlation model to \code{corrmodel = "GenWend"}, a compactly supported generalized Wendland model.
Input validation checks consistency of parameter names between \code{start} and \code{fixed} with respect to \code{CorrParam("GenWend")} and \code{NuisParam(model)}.

\strong{Reproducibility.}
If \code{seed} is supplied, \code{set.seed(seed)} is called before simulation/bootstrap operations. The user's previous random-number-generator state is restored when the function exits, so setting \code{seed} makes the test reproducible without permanently changing the global RNG state.

\strong{Likelihood selection.}
If both \code{likelihood} and \code{type} are \code{NULL}, the function selects \code{likelihood = "Marginal"}, \code{type = "Pairwise"} for large datasets and \code{likelihood = "Full"}, \code{type = "Standard"} otherwise. If only one of \code{likelihood} or \code{type} is supplied, an error is raised.

\strong{Threshold \code{h0} and minimum distance.}
Let \eqn{d_{min}}{d_min} be the minimum inter-point distance, computed with \pkg{nabor}. If \code{h0 = NULL}, the test defaults to an independence test with \code{h0_val = max(1e-6, d_min - 1e-6)}.
Under \strong{H0}, the support parameter is constrained by \code{scale <= h0}, or fixed to \code{h0_val} in the independence case. Under \strong{H1}, the scale is estimated freely, with a large data-driven upper bound.

\strong{Likelihood-ratio statistic.}
Let \eqn{\Lambda_{obs} = \max\{0, 2(\ell_1 - \ell_0)\}}{Lambda_obs = max(0, 2*(l1 - l0))}, where \eqn{\ell_0} and \eqn{\ell_1} are the log-likelihood or log-composite-likelihood values under \strong{H0} and \strong{H1}. The null distribution of \eqn{\Lambda} is approximated by parametric bootstrap with \code{B} datasets simulated under \strong{H0}. The \emph{p}-value is computed as
\deqn{(1 + \#\{\Lambda_b \ge \Lambda_{obs}\}) / (1 + B_{valid}).}

\strong{Simulation.}
For \code{method = "cholesky"}, bootstrap datasets are simulated with \code{GeoSim}. For \code{method = "TB"}, they are simulated with \code{GeoSimapprox}, using \code{L} spectral components. When \code{parallel = TRUE}, \code{GeoSimapprox} may parallelize over the \code{B} independent bootstrap replications. The function saves and restores the user's \pkg{future} plan and \pkg{progressr} handlers.

\strong{Early exit.}
If the unconstrained fit under \strong{H1} already satisfies \code{scale <= h0} or \code{scale <= h0_val}, the function returns \code{lambda_obs = 0} and \code{pvalue = 1} without running the bootstrap.

\strong{Robust bounds.}
The helper routines widen degenerate or ultra-narrow parameter boxes and enforce positivity for \code{sill} and \code{scale} where needed.
}
\value{
An invisible list with components:
\item{d_min}{Minimum inter-point distance \eqn{d_{min}}.}
\item{h0}{Threshold actually used under \strong{H0}.}
\item{lambda_obs}{Observed likelihood-ratio statistic \eqn{\Lambda_{obs}}.}
\item{pvalue}{Bootstrap \emph{p}-value.}
\item{B_rep}{Vector of bootstrap statistics \eqn{\Lambda_b} for valid replications.}
\item{fit_H0}{\code{GeoFit} object fitted under \strong{H0}.}
\item{fit_H1}{\code{GeoFit} object fitted under \strong{H1}.}
\item{seed}{Seed supplied by the user, or \code{NULL}.}
}
\note{
The function uses \code{corrmodel = "GenWend"} internally. For large spatial datasets, \code{method = "TB"} may be preferable to \code{method = "cholesky"}, especially when used with a suitable value of \code{L}.
}

\seealso{
\code{\link{GeoFit}}, \code{\link{GeoSim}}, \code{\link{GeoSimapprox}}, \code{\link{GeoVarest}}
}

\author{
Moreno Bevilacqua \email{moreno.bevilacqua89@gmail.com} \url{https://sites.google.com/view/moreno-bevilacqua/home}\cr
Víctor Morales Oñate \email{victor.morales@uv.cl} \url{https://sites.google.com/site/moralesonatevictor/}\cr
Christian Caamaño-Carrillo \email{chcaaman@ubiobio.cl} \url{https://www.researchgate.net/profile/Christian-Caamano}
}

\keyword{Bootstrap}
\keyword{Composite}
\keyword{HypothesisTest}
