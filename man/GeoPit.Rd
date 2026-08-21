\name{GeoPit}
\alias{GeoPit}
\encoding{UTF-8}
\title{Probability Integral Transform for Fitted GeoModels Objects}
\description{
Transforms observations or prediction targets using the fitted marginal
cumulative distribution function. The output can be represented on the
uniform probability-integral-transform scale or on the standard Gaussian-score
scale.
}
\usage{
GeoPit(object, type = c("Uniform", "Gaussian"), data_to_pred = NULL)
}
\arguments{
  \item{object}{An object of class \code{GeoFit}, \code{GeoKrig}, or
  \code{GeoKrigloc}. The spatial bivariate case is not currently supported.}
  \item{type}{Character string. \code{"Uniform"} returns fitted PIT values
  \eqn{F_i(Y_i)}. \code{"Gaussian"} returns
  \eqn{\Phi^{-1}\{F_i(Y_i)\}}. Partial matching is performed by
  \code{match.arg}.}
  \item{data_to_pred}{Optional observed values at the prediction locations when
  \code{object} is a \code{GeoKrig} or \code{GeoKrigloc} object. This is
  required unless the object already contains a \code{data_to_pred} component.}
}
\details{

For a \code{GeoFit} object, marginal transformations use the complete fitted
mean: a site-specific fixed vector when present, otherwise \eqn{X\beta} with
coefficients \code{mean}, \code{mean1}, and so on matched to the columns of
\code{X}.  The intercept-only case uses scalar \code{mean}.

For a \code{GeoFit} object, site-specific fitted means are reconstructed from
the design matrix and fitted regression coefficients. An externally supplied
fixed mean vector takes precedence. The fitted marginal CDF is then evaluated
at each observation. Copula parameters do not enter the PIT because the
transformation is marginal. The Student-t and skew Student-t transformations
use the fitted degrees of freedom \eqn{1/df} without integer rounding. The
continuous PIT implementation includes the bounded Kumaraswamy families,
Logistic, and SkewLaplace marginals used by the copula diagnostics.

For \code{type = "Gaussian"}, PIT values are clipped only for numerical
protection before applying the standard-normal quantile function. For discrete
marginal models, the current transformation is the ordinary non-randomized PIT
and therefore is not exactly uniform. For \code{PoissonGamma}, the fitted
marginal CDF is the negative-binomial CDF with mean \code{exp(mean)} (or the
site-specific fitted mean) and size parameter \code{shape}. Pearson-residual objects returned by
\code{GeoResiduals()} are rejected; PIT diagnostics for discrete models must use
the original fitted object.

For \code{GeoKrig} and \code{GeoKrigloc} objects, predictive PIT is currently
available only for Gaussian models and requires predictive MSE values.  Linear
kriging for the non-Gaussian models does not by itself define the full conditional
predictive distribution needed for a valid PIT.
}
\value{
The input object with its \code{data} component replaced by the transformed
values.
}
\examples{
\dontrun{
fit_u <- GeoPit(fit, type = "Uniform")
fit_z <- GeoPit(fit, type = "Gaussian")
}
}
\keyword{spatial}
\keyword{distribution}
