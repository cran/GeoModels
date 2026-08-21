\name{GeoCorrFct_Cop}
\alias{GeoCorrFct_Cop}
\encoding{UTF-8}
\title{Correlation, Covariance, and Variogram Functions for Spatial Copula Models}
\description{
Computes the observed-scale correlation, covariance, or semivariogram of a
univariate spatial or spatio-temporal random field constructed with a Gaussian,
Clayton-like, or skew-Gaussian copula. The function uses the same numerical
second-order covariance engines as \code{GeoCovmatrix()}, \code{GeoCovariogram()},
and \code{GeoKrig()}.
}
\usage{
GeoCorrFct_Cop(x, t = NULL, corrmodel,
model = "Gaussian", copula = "Gaussian",
distance = "Eucl", param, radius = 6371,
n = 1, covariance = FALSE, variogram = FALSE)
}

\arguments{
 \item{x}{Numeric vector of non-negative spatial distances.}
 \item{t}{Optional numeric vector of non-negative temporal distances for a
 spatio-temporal correlation model.}
 \item{corrmodel}{String giving the latent Gaussian correlation model. See
 \code{\link{GeoCovmatrix}} for implemented correlation models.}
 \item{model}{String giving the marginal distribution. Copula covariance is
 implemented for \code{Gaussian}, \code{StudentT}, \code{LogGaussian},
 \code{Gamma}, \code{Weibull}, \code{Beta}, \code{Beta2},
 \code{Kumaraswamy}, \code{Kumaraswamy2}, \code{Logistic},
 \code{SkewLaplace}, \code{Tukeyh}, \code{Tukeyh2}, and
 \code{SinhAsinh}.}
 \item{copula}{String giving the copula. Supported values are
 \code{Gaussian}, \code{Clayton}, and \code{SkewGaussian}. For the
 Clayton-like copula, \code{param$nu} must be a positive integer. For the
 skew-Gaussian copula, \code{param$nu} is the asymmetry parameter and must
 belong to \eqn{(-1,1)}.}
 \item{distance}{String giving the spatial distance. The default is
 \code{Eucl}. The supplied \code{x} values are already distances; this
 argument is retained for consistency with the other GeoModels functions.}
 \item{param}{List of correlation, marginal, nugget, and copula parameters.
 Since this is a lag-only function and has no design matrix, location-dependent
 margins use an intercept-only marginal predictor \code{mean}.}
 \item{radius}{Numeric radius of the sphere when great-circle distances are
 used. The default is 6371 km.}
 \item{n}{Numeric number of trials, retained for interface compatibility.}
 \item{covariance}{Logical. If \code{TRUE} and \code{variogram=FALSE}, return
 observed-scale covariance. If \code{FALSE}, return observed-scale
 correlation.}
 \item{variogram}{Logical. If \code{TRUE}, return a semivariogram. With
 \code{covariance=TRUE}, this is the observed-scale semivariogram
 \eqn{V_Y-C_Y(h)}. With \code{covariance=FALSE}, it is the normalized
 semivariogram \eqn{1-\rho_Y(h)}.}
}

\details{
The spatial dependence supplied by \code{corrmodel} is interpreted as the
correlation \eqn{\rho(h)} of the latent Gaussian random field underlying the
copula. The marginal covariance is then obtained on the observed scale using
the selected copula and marginal quantile transformation.

For Gaussian copulas, the implementation uses the Gaussian-copula covariance
engine. The \code{Tukeyh}, \code{Tukeyh2}, and \code{SinhAsinh} margins are
monotone transforms of one Gaussian field, so their Gaussian-copula covariance
is evaluated with the corresponding exact transformed-Gaussian covariance rather
than a truncated Hermite expansion. Clayton-like covariances use the cached
deterministic quadrature/interpolation engine, avoiding adaptive two-dimensional
integration at every lag. Skew-Gaussian covariances use the cached bivariate
Hermite representation. Consequently, repeated calls with the same copula and
marginal parameters can reuse cached numerical objects.

The nugget is applied to the latent correlation for non-zero lags. At the exact
zero spatial and temporal lag, the function represents the same random variable
and therefore returns correlation one, the marginal variance for covariance,
and zero for the semivariogram.
}

\value{
An object of class \code{GeoCorrFct}. Its \code{corr} component contains the
requested correlation, covariance, or semivariogram values. The object also
contains the spatial distances, temporal distances, marginal model, copula,
parameters, and flags describing the requested scale.
}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{
library(GeoModels)

################################################################
### Correlation of a mean-reparametrized Beta random field
### with a Matern latent correlation model.
################################################################

x <- seq(0, 0.4, 0.02)

param <- list(smooth = 0.5, scale = 0.2 / 3, nugget = 0,
              mean = 0, min = 0, max = 1, shape = 2)

corr_gauss <- GeoCorrFct_Cop(
    x = x, corrmodel = "Matern", param = param,
    copula = "Gaussian", model = "Beta2"
)
plot(corr_gauss, ylab = "Correlation", lwd = 2)

## Clayton-like and skew-Gaussian copulas are selected through param$nu.
## Their fast covariance engines are also used by GeoCovmatrix and GeoKrig.
}
