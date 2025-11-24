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
                  method = "Cholesky",
                  parallel = TRUE,
                  ncores = NULL,
                  progress = TRUE)
}
\arguments{
  \item{data}{Numeric vector of observations of length \code{n}.}
  \item{coordx}{\code{n x d} matrix with spatial coordinates (at least 2 columns).}
  \item{start}{\emph{Named} list of parameters \strong{to be estimated}. Names must match \code{CorrParam(corrmodel)} and \code{NuisParam(model)}.}
  \item{fixed}{\emph{Named} list of parameters \strong{kept fixed}. Together with \code{start} it must exactly cover the required set.}
  \item{model}{Marginal model for the data (default \code{"Gaussian"}).}
  \item{h0}{Positive threshold for the \emph{scale} parameter under \strong{H0}. If \code{NULL} (default), \strong{H0} represents \emph{spatial independence} with \code{scale = h0_val} set to the (adjusted) minimum inter-point distance.}
  \item{optimizer}{Optimization method passed to \code{GeoFit} (e.g., \code{"bobyqa"}).}
  \item{lower, upper}{\emph{Named} lists with bounds for parameters in \code{start} only. \code{NULL} means no extra bounds beyond defaults. Bounds are internally widened if too narrow.}
  \item{neighb}{Number of neighbors for composite likelihood when applicable. Ignored when the full likelihood is used.}
  \item{B}{Number of parametric bootstrap replications (default \code{1000}).}
  \item{likelihood}{Type of (composite) likelihood to pass to \code{GeoFit} (e.g., \code{"Marginal"}).}
  \item{type}{Composite-likelihood type (e.g., \code{"Pairwise"}).}
  \item{method}{Simulation method for bootstrap under \strong{H0}: \code{"Cholesky"} or \code{"TB"} (turning bands via \code{GeoSimapprox}).}
  \item{parallel}{If \code{TRUE}, bootstrap estimations are parallelized using \pkg{future}.}
  \item{ncores}{Number of cores for parallelization. If \code{NULL}, it is selected automatically up to \code{B} or \code{parallel::detectCores()-1}.}
  \item{progress}{If \code{TRUE}, progress bars are shown via \pkg{progressr}.}
}
\details{
The function fixes the correlation model to \code{corrmodel = "GenWend"} (generalized Wendland, compact support).
Input validation ensures consistency of parameter names between \code{start} and \code{fixed} w.r.t. \code{CorrParam(corrmodel)} and \code{NuisParam(model)}.

\strong{Automatic likelihood selection.}
If \eqn{n > 10{,}000}, a composite-likelihood is used (\code{likelihood = "Marginal"}, \code{type = "Pairwise"}) with \code{neighb} neighbors; otherwise the full likelihood is used (\code{likelihood = "Full"}, \code{type = "Standard"}).

\strong{Threshold \code{h0} and minimum distance.}
Let \eqn{d_{min}}{d_min} be the minimum inter-point distance, computed with \pkg{nabor}. If \code{h0 = NULL}, the test defaults to independence with \code{h0\_val = max(1e-6, d_min - 1e-6)}.
Under \strong{H0} one enforces \code{scale <= h0} (or \code{scale = h0\_val} for independence), while under \strong{H1} the scale is estimated freely but capped above by a large, data-driven \code{U\_big}.

\strong{Likelihood-ratio statistic.}
Let \eqn{\Lambda_{obs} = \max\{0, 2(\ell_1 - \ell_0)\}}{Lambda_obs = max(0, 2*(l1 - l0))} where \eqn{\ell_0} and \eqn{\ell_1} are the log-(composite) likelihoods under \strong{H0} and \strong{H1}. The null distribution of \eqn{\Lambda} is approximated via parametric bootstrap with \code{B} datasets simulated under \strong{H0}; the \emph{p}-value is computed as
\deqn{(1 + \#\{\Lambda_b \ge \Lambda_{obs}\}) / (1 + B_{valid}).}

\strong{Early exit.} If the unconstrained fit under \strong{H1} already satisfies \code{scale <= h0} (or \code{<= h0\_val}), the function returns \code{lambda\_obs = 0} and \code{pvalue = 1} without running the bootstrap.

\strong{Robust bounds.} The helper automatically widens degenerate or ultra-narrow boxes and enforces positivity for \code{sill} and \code{scale}.
}
\value{
An (invisible) list with components:
\item{d_min}{Minimum inter-point distance \eqn{d_{min}}.}
\item{h0}{Threshold actually used under \strong{H0}.}
\item{lambda_obs}{Observed likelihood-ratio statistic \eqn{\Lambda_{obs}}.}
\item{pvalue}{Bootstrap \emph{p}-value.}
\item{B_rep}{Vector of bootstrap statistics \eqn{\Lambda_b} for valid replications.}
\item{fit_H0}{\code{GeoFit} object fitted under \strong{H0}.}
\item{fit_H1}{\code{GeoFit} object fitted under \strong{H1}.}
}
\note{
The function uses \code{corrmodel = "GenWend"} internally.
For very large \code{n}, \code{method = "TB"} may be preferable for speed.
}

\seealso{
\code{\link{GeoFit}}, \code{\link{GeoSim}}, \code{\link{GeoSimapprox}}, \code{\link{GeoVarestbootstrap}}
}

\author{
Moreno Bevilacqua \email{moreno.bevilacqua89@gmail.com} \url{https://sites.google.com/view/moreno-bevilacqua/home}\cr
Víctor Morales Oñate \email{victor.morales@uv.cl} \url{https://sites.google.com/site/moralesonatevictor/}\cr
Christian Caamaño-Carrillo \email{chcaaman@ubiobio.cl} \url{https://www.researchgate.net/profile/Christian-Caamano}
}

\keyword{Bootstrap}
\keyword{Composite}
\keyword{HypothesisTest}
