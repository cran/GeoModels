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
 alpha = 0.95, L = 1000,
 parallel = TRUE, ncores = NULL, progress = TRUE,
 score_method = c("geofit_score", "finite"),
 eps = 1e-5, seed = NULL)
}
\arguments{
 \item{fit}{A fitted object obtained from \code{\link{GeoFit}}. The object
 must contain the sensitivity matrix, hence \code{GeoFit} should be called
 with \code{sensitivity = TRUE}.}

 \item{K}{The number of simulations used in the parametric score bootstrap.}

 \item{sparse}{Logical; if \code{TRUE}, then Cholesky decomposition is
 performed using sparse matrix algorithms.}

 \item{method}{String; the method of simulation. The default is
 \code{"cholesky"}. For large data sets the options \code{"TB"} and
 \code{"CE"} call approximate simulation methods; see
 \code{\link{GeoSimapprox}}.}

 \item{alpha}{Numeric; the level of the confidence interval.}

 \item{L}{Numeric; the number of lines in the turning bands method.}

 \item{parallel}{Logical; if \code{TRUE}, the score evaluation step is
 parallelized.}

 \item{ncores}{Numeric; number of cores involved in the parallel score
 evaluation. If \code{NULL}, the number of cores is chosen automatically.}

 \item{progress}{Logical; if \code{TRUE}, progress information is shown.}

 \item{score_method}{String; the method used to compute the composite score
 at the fitted parameter value. If \code{"geofit_score"}, the score returned
 by \code{GeoFit(..., onlyvar = TRUE, score = TRUE)} is used. If
 \code{"finite"}, the score is computed by external finite differences of
 the composite log-likelihood evaluated through \code{GeoFit}.}

 \item{eps}{Numeric; finite difference step multiplier used only when
 \code{score_method = "finite"}.}

 \item{seed}{Optional integer seed for reproducibility of the simulated
 samples.}
}
\details{
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
 matrix. As a consequence, the component \code{estimates} returned by this
 function contains simulated scores evaluated at \eqn{\hat\theta}, not
 bootstrap parameter estimates.
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
 \item{estimates}{Matrix containing the successful scores and the associated
 composite log-likelihood values. These are not refitted parameter
 estimates.}
}
\seealso{
\code{\link{GeoFit}} for the fitted objects used as input,
  \code{\link{GeoVarestbootstrap}} for bootstrap variance estimation,
  \code{\link{CompLik}} for composite-likelihood calculations.
}
\author{
Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},
\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl},
\url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},
\url{https://www.researchgate.net/profile/Christian-Caamano}
}
\keyword{Composite}
