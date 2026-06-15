\name{GeoScores}
\alias{GeoScores}
\encoding{UTF-8}

\title{Computation of predictive scores}

\description{
  The function computes predictive scores for observed validation values and
  predicted values, optionally using predictive variances or 95 percent
  prediction intervals.
}

\usage{
GeoScores(data_to_pred,
          probject = NULL,
          pred = NULL,
          mse = NULL,
          score = c("pe", "crps", "intscore", "coverage"),
          lower95 = NULL,
          upper95 = NULL,
          threshold = 0.5,
          na.rm = TRUE)
}

\arguments{
  \item{data_to_pred}{
    Numeric vector, matrix or array containing the observed validation values.
    Values are internally coerced to a numeric vector.
  }

  \item{probject}{
    Optional object of class \code{GeoKrig} or \code{GeoKrigloc}. If supplied,
    point predictions and prediction variances are extracted from
    \code{probject$pred} and \code{probject$mse}, respectively.
  }

  \item{pred}{
    Numeric vector, matrix or array of point predictions. This argument is
    required when \code{probject} is not supplied.
  }

  \item{mse}{
    Optional numeric vector, matrix or array of prediction variances. Predictive
    standard errors are computed as \code{sqrt(mse)}. This argument is required
    for Gaussian-distribution scores unless both \code{lower95} and
    \code{upper95} are supplied.
  }

  \item{score}{
    Character vector specifying which predictive scores should be computed.
    Possible values are \code{"brie"}, \code{"crps"}, \code{"lscore"},
    \code{"pit"}, \code{"pe"}, \code{"intscore"}, and \code{"coverage"}.
    Several values can be supplied.
  }

  \item{lower95}{
    Optional numeric vector, matrix or array containing the lower bounds of the
    95 percent prediction intervals. If supplied together with
    \code{upper95}, these bounds are used directly for \code{"intscore"} and
    \code{"coverage"}. If \code{mse} is not supplied, they are also used to
    approximate predictive standard errors for Gaussian-distribution scores.
  }

  \item{upper95}{
    Optional numeric vector, matrix or array containing the upper bounds of the
    95 percent prediction intervals. See \code{lower95}.
  }

  \item{threshold}{
    Numeric threshold used to compute the Brier score \code{"brie"}. The binary
    event is \eqn{Y > threshold}.
  }

  \item{na.rm}{
    Logical. If \code{TRUE}, non-finite observations, predictions, invalid
    intervals and non-positive predictive variances or standard errors are
    removed before computing the scores. If \code{FALSE}, such values produce
    an error.
  }
}

\details{
  \code{GeoScores} computes point prediction scores and probabilistic
  predictive scores. Predictions and prediction variances can be supplied
  directly through \code{pred} and \code{mse}, or indirectly through a
  \code{GeoKrig} or \code{GeoKrigloc} object supplied to \code{probject}.

  The point prediction option \code{"pe"} returns the mean absolute error and
  the root mean squared prediction error.

  The scores \code{"brie"}, \code{"lscore"}, \code{"pit"} and \code{"crps"}
  are computed under the assumption that the predictive distribution is
  Gaussian, with mean equal to the point prediction and standard deviation equal
  to the predictive standard error. If \code{mse} is supplied, the predictive
  standard error is \code{sqrt(mse)}. If \code{mse} is absent and both
  \code{lower95} and \code{upper95} are supplied, the predictive standard error
  is approximated as
  \deqn{
    \frac{U - L}{2 \Phi^{-1}(0.975)},
  }
  where \eqn{L} and \eqn{U} are the lower and upper bounds of the 95 percent
  prediction interval.

  The Brier score is computed for the binary event \eqn{Y > threshold}, using
  the Gaussian predictive probability
  \eqn{1 - \Phi((threshold - \hat{Y}) / s)}, where \eqn{\hat{Y}} is the
  prediction and \eqn{s} is the predictive standard error.

  The logarithmic score is returned as the negative log predictive density,
  averaged over validation observations:
  \deqn{
    \frac{1}{2}\left\{\log(2\pi\sigma^2) +
    \left(\frac{y - \hat{y}}{\sigma}\right)^2\right\}.
  }

  The CRPS is computed using the closed-form expression for a Gaussian
  predictive distribution.

  For \code{"pit"}, if \code{probject} is a \code{GeoKrig} or
  \code{GeoKrigloc} object and \code{GeoPit} is available, the PIT values are
  obtained through \code{GeoPit(probject, type = "Uniform")$data}; otherwise
  they are computed as \code{pnorm(data_to_pred, mean = pred, sd = sqrt(mse))}.

  The 95 percent interval score is computed as
  \deqn{
    (u - l) +
    \frac{2}{\alpha}(l - y) I(y < l) +
    \frac{2}{\alpha}(y - u) I(y > u),
  }
  with \eqn{\alpha = 0.05}. If \code{lower95} and \code{upper95} are supplied,
  they are used directly as \eqn{l} and \eqn{u}; otherwise the interval is
  computed from the Gaussian predictive distribution using \code{pred} and
  \code{mse}. The empirical coverage is the proportion of validation values
  falling inside the corresponding 95 percent intervals.
}

\value{
  A list containing only the components requested through \code{score}. Possible
  components are:
  \item{MAE}{
    Mean absolute error. Returned when \code{"pe"} is requested.
  }
  \item{RMSPE}{
    Root mean squared prediction error. Returned when \code{"pe"} is requested.
  }
  \item{Brier}{
    Brier score for the binary event \eqn{Y > threshold}. Returned when
    \code{"brie"} is requested.
  }
  \item{LogScore}{
    Mean negative log predictive density under a Gaussian predictive
    distribution. Returned when \code{"lscore"} is requested.
  }
  \item{PIT}{
    Probability integral transform values. Returned when \code{"pit"} is
    requested.
  }
  \item{CRPS}{
    Mean continuous ranked probability score under a Gaussian predictive
    distribution. Returned when \code{"crps"} is requested.
  }
  \item{IS95}{
    Mean interval score for the 95 percent prediction intervals. Returned when
    \code{"intscore"} is requested.
  }
  \item{Cvg95}{
    Empirical coverage of the 95 percent prediction intervals. Returned when
    \code{"coverage"} is requested.
  }
}

\references{
  Gneiting, T. and Raftery, A. E. (2007).
  \emph{Strictly Proper Scoring Rules, Prediction, and Estimation}.
  Journal of the American Statistical Association, \bold{102}, 359--378.

  Heaton, M. J., Datta, A., Finley, A. O., Furrer, R., Guinness, J.,
  Guhaniyogi, R., Gerber, F., Gramacy, R. B., Hammerling, D.,
  Katzfuss, M., Lindgren, F., Nychka, D. W., Sun, F.,
  and Zammit-Mangion, A. (2019).
  \emph{A Case Study Competition Among Methods for Analyzing Large Spatial Data}.
  Journal of Agricultural, Biological, and Environmental Statistics,
  \bold{24}, 398--425.
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

\examples{
library(GeoModels)

################################################################
######### Example of predictive score computation  #############
################################################################

model <- "Gaussian"
set.seed(79)
N <- 1000

x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x, y)

# Set covariance parameters
corrmodel <- "GenWend"
mean <- 0
sill <- 5
nugget <- 0
scale <- 0.2
smooth <- 0
power2 <- 4

param <- list(mean = mean, sill = sill, nugget = nugget,
              scale = scale, smooth = smooth, power2 = power2)

# Simulation of the spatial Gaussian random field
data <- GeoSim(coordx = coords, corrmodel = corrmodel,
               param = param)$data

# Training and validation split
sel <- sample(1:N, N * 0.8)

coords_est <- coords[sel, ]
coords_to_pred <- coords[-sel, ]

data_est <- data[sel]
data_to_pred <- data[-sel]

# Pairwise likelihood fitting
fixed <- list(nugget = nugget, smooth = smooth, power2 = power2)

start <- list(mean = 0, scale = scale, sill = 1)

I <- Inf
lower <- list(mean = -I, scale = 0, sill = 0)
upper <- list(mean =  I, scale = I, sill = I)

fit <- GeoFit(data_est, coordx = coords_est,
              corrmodel = corrmodel, model = model,
              likelihood = "Marginal", type = "Pairwise",
              neighb = 3, optimizer = "nlminb",
              lower = lower, upper = upper,
              start = start, fixed = fixed)

# Prediction at validation locations
pr <- GeoKrig(fit,loc = coords_to_pred,  data = data_est, mse = TRUE)

# Predictive scores from pred and mse
Pr_scores <- GeoScores(data_to_pred, pred = pr$pred, mse = pr$mse,
                       score = c("pe", "brie", "crps", "lscore",
                                 "pit", "intscore", "coverage"),
                       threshold = 0)

Pr_scores$MAE
Pr_scores$RMSPE
Pr_scores$Brier
Pr_scores$CRPS
Pr_scores$LogScore
Pr_scores$IS95
Pr_scores$Cvg95

}

\keyword{Predictive scores}
