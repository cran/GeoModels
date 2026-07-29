\name{GeoTests}
\alias{GeoTests}
\encoding{UTF-8}
\title{Statistical Hypothesis Tests for Nested Models}
\description{
Performs Wald and likelihood-ratio tests for successively nested
\code{GeoFit} models. Standard likelihood fits can be compared using
Wald or Wilks tests, whereas composite likelihood fits can be compared
using Wald or Satterthwaite-adjusted composite likelihood-ratio tests.
}

\usage{
GeoTests(object1, object2, \dots, statistic)
}

\arguments{
  \item{object1}{An unrestricted fitted model of class \code{GeoFit}.}
  \item{object2}{A restricted fitted model of class \code{GeoFit}, nested
  within \code{object1}.}
  \item{\dots}{Further fitted models of class \code{GeoFit}, ordered from
  less restricted to more restricted. Adjacent models are compared.}
  \item{statistic}{A character string specifying the test statistic.
  Available choices are \code{"Wald"}, \code{"Wilks"}, and
  \code{"WilksS"}. See Details.}
}

\details{
The models must be supplied from the unrestricted model to successively
more restricted models. For each adjacent pair, the tested parameters are
the parameters estimated in the unrestricted model but fixed in the
restricted model. Each tested parameter must therefore occur explicitly
in the \code{fixed} component of the restricted fit.

The function checks that adjacent models use the same data, coordinates,
design matrix, marginal model, correlation model, likelihood type, pair
construction, and other relevant fitting options. Common fixed parameters
must have the same values in the two models.

The following statistics are available:

\itemize{
  \item \code{"Wald"}: available for both standard and composite
  likelihood fits. The statistic is
  \deqn{
  W_{\mathrm{Wald}} =
  (\widehat{\boldsymbol\psi}-\boldsymbol\psi_0)^\top
  \widehat V_{\psi\psi}^{-1}
  (\widehat{\boldsymbol\psi}-\boldsymbol\psi_0),
  }
  where \eqn{\boldsymbol\psi} is the vector of tested parameters and
  \eqn{\widehat V_{\psi\psi}} is the corresponding block of the estimated
  covariance matrix of the unrestricted estimator. Under the null
  hypothesis, the statistic is asymptotically chi-squared with degrees of
  freedom equal to the number of restrictions. For composite likelihood
  fits, \code{varcov} must contain the Godambe covariance matrix, usually
  obtained by applying \code{GeoVarest} to the unrestricted fit.

  \item \code{"Wilks"}: available only for standard full-likelihood fits.
  The statistic is
  \deqn{
  W_{\mathrm{Wilks}} =
  2\{\ell(\widehat{\boldsymbol\theta})-
  \ell(\widetilde{\boldsymbol\theta})\},
  }
  where \eqn{\widehat{\boldsymbol\theta}} and
  \eqn{\widetilde{\boldsymbol\theta}} are the unrestricted and restricted
  maximum likelihood estimates, respectively. Under standard regularity
  conditions, the statistic is asymptotically chi-squared with degrees of
  freedom equal to the number of restrictions.

  \item \code{"WilksS"}: available only for composite likelihood fits,
  including marginal, conditional, and difference composite likelihoods.
  Let
  \deqn{
  W_{CL}=2\{\ell_{CL}(\widehat{\boldsymbol\theta})-
  \ell_{CL}(\widetilde{\boldsymbol\theta})\}.
  }
  Its asymptotic null distribution is generally a weighted sum
  \deqn{
  W_{CL}\ \dot{\sim}\ \sum_{j=1}^{r}\lambda_j\chi^2_{1,j},
  }
  where \eqn{r} is the number of restrictions. In the implementation,
  \eqn{\lambda_1,\ldots,\lambda_r} are the eigenvalues of
  \deqn{
  H_{\mathrm{eff}}V_{\psi\psi}, \qquad
  H_{\mathrm{eff}}=\{(H^{-1})_{\psi\psi}\}^{-1},
  }
  with \eqn{H} the sensitivity matrix and
  \eqn{V=H^{-1}JH^{-1}} the Godambe covariance matrix from the unrestricted
  fit. The Satterthwaite moment-matching constants are
  \deqn{
  \nu=\frac{(\sum_j\lambda_j)^2}{\sum_j\lambda_j^2},
  \qquad
  \kappa=\frac{\sum_j\lambda_j^2}{\sum_j\lambda_j},
  }
  and the adjusted statistic
  \deqn{
  W_S=\frac{W_{CL}}{\kappa}
  }
  is compared with a chi-squared distribution having \eqn{\nu} degrees of
  freedom. The unrestricted fit must contain both \code{sensmat} and the
  Godambe covariance matrix \code{varcov}; these can be obtained using
  \code{sensitivity = TRUE} in \code{GeoFit}, followed by \code{GeoVarest}.
}

For a composite likelihood-ratio comparison, the unrestricted and
restricted fits must be based on exactly the same retained pairs. This is
automatic for deterministic pair constructions. With stochastic thinning,
the two fits must use the same thinning realization; \code{GeoTests} checks
that their \code{rowidx} and \code{colidx} components are identical and
stops otherwise.

The matrices used by \code{"WilksS"} are evaluated at the unrestricted
fit. Under the null hypothesis, this is a consistent plug-in evaluation
because the unrestricted and restricted estimators converge to the same
parameter value.
}

\value{
A data frame containing one row for each supplied model. The first row
corresponds to the unrestricted model and contains missing test results.
Each subsequent row reports the comparison between that model and the
model in the preceding row. The columns are:

\describe{
  \item{Num.Par}{Number of estimated parameters in the model.}
  \item{Diff.Par}{Number of restrictions relative to the preceding model.}
  \item{Df}{Degrees of freedom used for the reference chi-squared
  distribution. For \code{"WilksS"}, this is the Satterthwaite effective
  degrees of freedom and need not be an integer.}
  \item{Chisq}{Observed value of the test statistic.}
  \item{Pr(>chisq)}{Upper-tail p-value from the chi-squared reference
  distribution with the reported degrees of freedom.}
}
}

\references{

Kent, J. T. (1982). Robust Properties of Likelihood Ratio Tests.
\emph{Biometrika}, \bold{69}, 19--27.

Satterthwaite, F. E. (1946). An Approximate Distribution of Estimates of
Variance Components. \emph{Biometrics Bulletin}, \bold{2}, 110--114.

Varin, C., Reid, N. and Firth, D. (2011). An Overview of Composite
Likelihood Methods. \emph{Statistica Sinica}, \bold{21}, 5--42.
}

\author{
Moreno Bevilacqua,
\email{moreno.bevilacqua89@gmail.com},
\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate,
\email{victor.morales@uv.cl},
\url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo,
\email{chcaaman@ubiobio.cl},
\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\seealso{
\code{\link{GeoFit}}, \code{\link{GeoVarest}}.
}

\examples{
library(GeoModels)

################################################################
### Example 1. Test of a regression coefficient using
### conditional pairwise composite likelihood
################################################################
set.seed(342)
model <- "Gaussian"

NN <- 1500
x <- runif(NN, 0, 1)
y <- runif(NN, 0, 1)
coords <- cbind(x, y)

mean <- 1
mean1 <- -1.25
nugget <- 0
sill <- 1

X <- cbind(1, runif(nrow(coords)))

corrmodel <- "Wend0"
power2 <- 4
c_supp <- 0.15

param <- list(
  power2 = power2, mean = mean, mean1 = mean1,
  sill = sill, scale = c_supp, nugget = nugget
)

data <- GeoSim(
  coordx = coords, corrmodel = corrmodel,
  model = model, param = param, X = X
)$data

I <- Inf

## H1: regression mean with a nonzero slope
fixed <- list(nugget = nugget, power2 = power2)
start <- list(mean = mean, mean1 = mean1, scale = c_supp, sill = sill)
lower <- list(mean = -I, mean1 = -I, scale = 0, sill = 0)
upper <- list(mean = I, mean1 = I, scale = I, sill = I)

fitH1 <- GeoFit(
  data = data, coordx = coords, corrmodel = corrmodel,
  model = model, likelihood = "Conditional", type = "Pairwise",
  sensitivity = TRUE, lower = lower, upper = upper, neighb = 3,
  optimizer = "nlminb", X = X, start = start, fixed = fixed
)

## H0: constant mean, that is, mean1 = 0
fixed0 <- list(power2 = power2, nugget = nugget, mean1 = 0)
start0 <- list(mean = mean, scale = c_supp, sill = sill)
lower0 <- list(mean = -I, scale = 0, sill = 0)
upper0 <- list(mean = I, scale = I, sill = I)

fitH0 <- GeoFit(
  data = data, coordx = coords, corrmodel = corrmodel,
  model = model, likelihood = "Conditional", type = "Pairwise",
  sensitivity = TRUE, lower = lower0, upper = upper0, neighb = 3,
  optimizer = "nlminb", X = X, start = start0, fixed = fixed0
)

\dontrun{
## GeoVarest is needed on the unrestricted fit to estimate its
## Godambe covariance matrix. Increase K for final analyses.
fitH1 <- GeoVarest(fitH1, K = 100, seed = 1234)

GeoTests(fitH1, fitH0, statistic = "Wald")
GeoTests(fitH1, fitH0, statistic = "WilksS")
}

################################################################
### Example 2. Test of Gaussianity for a Sinh--arcsinh random field
### using the standard full likelihood
################################################################
set.seed(99)
model <- "SinhAsinh"

NN <- 200
x <- runif(NN, 0, 1)
y <- runif(NN, 0, 1)
coords <- cbind(x, y)

mean <- 0
nugget <- 0
sill <- 1
skew <- 0
tail <- 1

corrmodel <- "Wend0"
power2 <- 4
c_supp <- 0.2

param <- list(
  power2 = power2, skew = skew, tail = tail,
  mean = mean, sill = sill, scale = c_supp, nugget = nugget
)

data <- GeoSim(
  coordx = coords, corrmodel = corrmodel,
  model = model, param = param
)$data

## H1: unrestricted Sinh--arcsinh model
fixed <- list(power2 = power2, nugget = nugget, mean = mean)
start <- list(scale = c_supp, skew = skew, tail = tail, sill = sill)
lower <- list(scale = 0, skew = -I, tail = 0, sill = 0)
upper <- list(scale = I, skew = I, tail = I, sill = I)

fitH1 <- GeoFit2(
  data = data, coordx = coords, corrmodel = corrmodel,
  model = model, likelihood = "Full", type = "Standard",
  varest = TRUE, lower = lower, upper = upper,
  optimizer = "nlminb", start = start, fixed = fixed
)

## H0: Gaussianity, that is, skew = 0 and tail = 1
fixed0 <- list(
  power2 = power2, nugget = nugget, mean = mean,
  tail = 1, skew = 0
)
start0 <- list(scale = c_supp, sill = sill)
lower0 <- list(scale = 0, sill = 0)
upper0 <- list(scale = 2, sill = 5)

fitH0 <- GeoFit(
  data = data, coordx = coords, corrmodel = corrmodel,
  model = model, likelihood = "Full", type = "Standard",
  varest = TRUE, lower = lower0, upper = upper0,
  optimizer = "nlminb", start = start0, fixed = fixed0
)

GeoTests(fitH1, fitH0, statistic = "Wald")
GeoTests(fitH1, fitH0, statistic = "Wilks")
}

\keyword{spatial}
\keyword{htest}
