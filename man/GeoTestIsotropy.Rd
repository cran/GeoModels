\name{GeoTestIsotropy}
\alias{GeoTestIsotropy}
\title{Parametric bootstrap test for isotropy}
\description{
Performs an isotropy test based on \emph{parametric bootstrap} by comparing an isotropic model (H0) with an anisotropic model (H1).
A user-chosen spatial correlation model is fitted to spatial data, and the test statistic is a likelihood ratio
between H1 (with anisotropy estimated) and H0 (anisotropy fixed at ratio = 1 and angle = 0). The \emph{p}-value is computed via parametric bootstrap.
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
 parallel = TRUE,
 ncores = NULL,
 progress = TRUE)
}
\arguments{
 \item{data}{Numeric vector of observations (length \code{n}).}
 \item{coordx}{\code{n x d} matrix with the spatial coordinates of the sites (at least 2 columns for anisotropy analysis).}
 \item{start}{\emph{Named} list with the parameters \strong{to be estimated}; names must match those expected by \code{CorrParam(corrmodel)} and \code{NuisParam(model)}.}
 \item{fixed}{\emph{Named} list with the parameters \strong{kept fixed} during estimation; together with \code{start} it must exactly cover the required set of parameters.}
 \item{optimizer}{Optimization algorithm to be used in \code{GeoFit} (e.g., \code{"bobyqa"}).}
 \item{model}{Marginal model for the data (e.g., \code{"Gaussian"}).}
 \item{corrmodel}{Spatial correlation model (e.g., \code{"Matern"}).}
 \item{lower, upper}{\emph{Named} lists with lower/upper bounds for parameters in \code{start} only. If \code{NULL}, no extra bounds beyond defaults.}
 \item{B}{Number of parametric bootstrap replications (default \code{1000}).}
 \item{likelihood}{Type of (composite) likelihood to pass to \code{GeoFit} (e.g., \code{"Marginal"}).}
 \item{type}{Composite-likelihood type (e.g., \code{"Pairwise"}).}
 \item{copula}{Copula object for simulation via \code{GeoSimCopula}; if \code{NULL} Gaussian simulation is used.}
 \item{neighb}{Number of neighbors for composite-likelihood estimation, where applicable.}
 \item{parallel}{If \code{TRUE}, bootstrap fits are distributed in parallel through \pkg{future}.}
 \item{ncores}{Number of cores to use in parallel; if \code{NULL} it is chosen automatically up to \code{B} or \code{parallel::detectCores()-1}.}
 \item{progress}{If \code{TRUE}, show progress bars using \pkg{progressr}.}
}
\details{
\code{GeoTestIsotropy} implements a \emph{likelihood ratio test} (LRT) between:
\itemize{
 \item \strong{H0}: isotropic model with anisotropy parameters fixed at \code{angle = 0}, \code{ratio = 1};
 \item \strong{H1}: anisotropic model with \code{angle} and \code{ratio} estimated.
}
The observed statistic is \eqn{\mathrm{LRT}_{\text{obs}} = 2(\ell_1 - \ell_0)}, truncated at zero if negative.
The null distribution is approximated by simulating \code{B} datasets under H0 with parameters equal to the H0 fit, and recomputing the LRT.
The \emph{p}-value is computed as \eqn{(1 + \#\{\mathrm{LRT}_b \ge \mathrm{LRT}_{\text{obs}}\}) / (1 + B_{\text{valid}})}.
The function includes robust checks on parameter names, bounds, simulation method choice (Cholesky vs TB), and parallelization.
}
\value{
An (invisible) list with components:
\item{statistic}{Observed LRT statistic.}
\item{pvalue}{Bootstrap \emph{p}-value.}
\item{ratio_hat}{Estimated anisotropy ratio (\code{> 1} indicates stretching along the principal axis).}
\item{angle_hat}{Estimated anisotropy angle (radians/degrees according to \code{GeoFit} convention).}
\item{corrmodel}{Correlation model used.}
\item{model}{Marginal model used.}
\item{fit_H0}{\code{GeoFit} object fitted under H0.}
\item{fit_H1}{\code{GeoFit} object fitted under H1.}
\item{B_rep}{Vector of LRT values obtained on the valid bootstrap replications.}
}


\seealso{
\code{\link{GeoVariogram}} for empirical variograms,
  \code{\link{GeoVariogramDir}} for directional variograms,
  \code{\link{GeoFit}} for covariance-model fitting.
}

\author{
Moreno Bevilacqua \email{moreno.bevilacqua89@gmail.com} \url{https://sites.google.com/view/moreno-bevilacqua/home}\cr
Víctor Morales Oñate \email{victor.morales@uv.cl} \url{https://sites.google.com/site/moralesonatevictor/}\cr
Christian Caamaño-Carrillo \email{chcaaman@ubiobio.cl} \url{https://www.researchgate.net/profile/Christian-Caamano}
}

\keyword{Bootstrap}
\keyword{Composite}
\keyword{HypothesisTest}
