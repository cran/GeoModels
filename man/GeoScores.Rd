\name{GeoScores}
\alias{GeoScores}
\encoding{UTF-8}

\title{Computation of predictive scores}

\description{
 The function computes predictive scores for observed validation values from
 point predictions, kriging prediction objects, or conditional-simulation
 objects. Gaussian predictive scores are available from \code{GeoKrig}/\code{GeoKrigloc}
 through the \code{pred}/\code{mse} convention, while \code{GeoSimcond}
 objects are scored from their empirical conditional predictive distribution.
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
 Optional object of class \code{GeoKrig}, \code{GeoKrigloc}, or \code{GeoSimcond}.
 For \code{GeoKrig}/\code{GeoKrigloc}, point predictions and prediction variances
 are extracted from \code{probject$pred} and \code{probject$mse}. Requested
 probabilistic scores are computed from the Gaussian predictive convention
 \eqn{N(\widehat Y, MSE)} for every marginal model. For a non-Gaussian model
 this is a Gaussian approximation based on the optimal linear predictor and its
 MSE; \code{GeoScores} reports it without deciding whether that approximation is
 appropriate for a particular analysis.

 For \code{GeoSimcond}, the conditional simulations in \code{probject$condsim}
 define an empirical predictive distribution. Point scores use
 \code{probject$cond_mean}; CRPS, PIT, Brier score, prediction intervals,
 interval score, and coverage are computed directly from the conditional draws,
 without a Gaussian approximation.
 }

 \item{pred}{
 Numeric vector, matrix or array of point predictions. This argument is
 required when \code{probject} is not supplied.
 }

 \item{mse}{
 Optional numeric vector, matrix or array of prediction variances. Predictive
 standard errors are computed as \code{sqrt(mse)}. This argument is not needed
 for \code{"pe"}, for interval score or coverage when \code{lower95} and
 \code{upper95} are supplied, or when \code{probject} is a
 \code{GeoSimcond} object. It is required for probabilistic scores based on the
 Gaussian \code{pred}/\code{mse} convention unless both interval bounds are
 supplied and used to infer a Gaussian standard error.
 }

 \item{score}{
 Character vector specifying which predictive scores should be computed.
 Possible values are \code{"brie"}, \code{"crps"}, \code{"lscore"},
 \code{"pit"}, \code{"pe"}, \code{"intscore"}, and \code{"coverage"}.
 Several values can be supplied.
 }

 \item{lower95}{
 Optional numeric vector, matrix or array containing the lower bounds of the
 95 percent prediction intervals. If supplied together with \code{upper95},
 these bounds are used directly for \code{"intscore"} and \code{"coverage"}.
 For a \code{GeoSimcond} object they override the default empirical 2.5 and
 97.5 percent conditional-simulation quantiles. If \code{mse} is not supplied
 in the Gaussian \code{pred}/\code{mse} path, they can also be used to infer a
 Gaussian predictive standard error.
 }

 \item{upper95}{
 Optional numeric vector, matrix or array containing the upper bounds of the
 95 percent prediction intervals. See \code{lower95}.
 }

 \item{threshold}{
 Finite numeric scalar used to compute the Brier score \code{"brie"}. The binary event is \eqn{Y > threshold}.
 }

 \item{na.rm}{
 Logical. If \code{TRUE}, invalid observations are removed separately for each requested score class. Thus an invalid predictive variance does not remove an otherwise valid observation from MAE/RMSPE, and an invalid interval does not affect point scores. If \code{FALSE}, invalid values relevant to any requested score produce an error.
 }
}

\details{
 \code{GeoScores} dispatches predictive scoring according to the information
 contained in \code{probject}.

 For a \code{GeoKrig} or \code{GeoKrigloc} object, \code{pred} is the point
 predictor and \code{mse} its mean squared prediction error. The option
 \code{"pe"} returns mean absolute error and root mean squared prediction
 error. Whenever probabilistic scores are requested, \code{GeoScores} uses the
 Gaussian predictive convention
 \deqn{
 Y_0\mid Y \ \dot\sim\ N\{\widehat Y_0, MSE_0\}.
 }
 This convention is used for every marginal model. For a Gaussian random field
 it is the natural conditional predictive distribution (conditional on fitted
 parameters). For a non-Gaussian model, \code{GeoKrig}/\code{GeoKrigloc}
 still returns an optimal linear predictor, so CRPS, LogScore, PIT, Brier score,
 and Gaussian intervals obtained this way must be interpreted as Gaussian
 predictive approximations based on the OLP and its MSE. \code{GeoScores} does
 not suppress those quantities; their scientific interpretation is left to the
 user.

 Predictions and prediction variances can equivalently be supplied directly
 through \code{pred} and \code{mse}; this uses the same Gaussian predictive
 convention.

 For a \code{GeoSimcond} object, the replications stored in
 \code{probject$condsim} are treated as an empirical conditional predictive
 distribution at each prediction location. No Gaussian approximation is used.
 The point forecast for \code{"pe"} is \code{probject$cond_mean} (or the
 empirical conditional mean reconstructed from the simulations if needed).
 The Brier probability for the event \eqn{Y>threshold} is the empirical
 conditional exceedance probability. PIT is the empirical conditional CDF at
 the validation observation. The default 95 percent predictive interval is
 given by the empirical 0.025 and 0.975 quantiles, unless explicit
 \code{lower95}/\code{upper95} values are supplied.

 The empirical CRPS from conditional draws
 \eqn{x_1,\ldots,x_B} is
 \deqn{
 \widehat{CRPS}(F,y)=
 \frac{1}{B}\sum_{b=1}^B |x_b-y|-
 \frac{1}{2B^2}\sum_{b=1}^B\sum_{c=1}^B |x_b-x_c|.
 }
 It is evaluated internally with an equivalent sorted-sample formula requiring
 \eqn{O(B\log B)} rather than \eqn{O(B^2)} work per prediction location.

 A logarithmic score is not computed automatically from a
 \code{GeoSimcond} object. Conditional draws identify the predictive
 distribution empirically but do not by themselves provide a predictive
 density evaluated exactly at the observed value. If \code{"lscore"} is
 requested from a \code{GeoSimcond} object, \code{LogScore = NA} is returned
 with a warning; no kernel-density approximation is introduced silently.

 In the Gaussian \code{pred}/\code{mse} path, the Brier score uses the
 Gaussian predictive probability of \eqn{Y>threshold}; the logarithmic score
 is the mean negative Gaussian log predictive density; PIT uses the Gaussian
 CDF; and CRPS uses the closed-form Gaussian expression.

 The 95 percent interval score is
 \deqn{
 (u-l)+\frac{2}{\alpha}(l-y)I(y<l)+
 \frac{2}{\alpha}(y-u)I(y>u), \qquad \alpha=0.05.
 }
 For \code{GeoKrig}/\code{GeoKrigloc} without explicit interval bounds,
 \eqn{l} and \eqn{u} are Gaussian intervals based on \code{pred} and
 \code{mse}. For \code{GeoSimcond}, they are empirical conditional-simulation
 quantiles. Empirical coverage is the proportion of validation observations
 lying in the corresponding intervals.
}
\value{
 A list containing only the components requested through \code{score}:
 \item{MAE}{Mean absolute error. For \code{GeoSimcond}, the empirical
 conditional mean is used as the point forecast.}
 \item{RMSPE}{Root mean squared prediction error.}
 \item{Brier}{Brier score for the event \eqn{Y>threshold}. It is Gaussian
 probability based for \code{GeoKrig}/\code{GeoKrigloc} and empirical
 conditional-simulation based for \code{GeoSimcond}.}
 \item{LogScore}{Mean negative Gaussian log predictive density for the
 \code{pred}/\code{mse} or \code{GeoKrig}/\code{GeoKrigloc} path. For a
 \code{GeoSimcond} object this component is \code{NA} if requested.}
 \item{PIT}{Gaussian PIT values for the \code{pred}/\code{mse} or kriging
 path, or empirical conditional CDF values for \code{GeoSimcond}.}
 \item{CRPS}{Mean Gaussian CRPS for the kriging/direct-prediction path, or
 mean empirical CRPS from conditional simulations for \code{GeoSimcond}.}
 \item{IS95}{Mean interval score for the 95 percent prediction intervals.}
 \item{Cvg95}{Empirical coverage of the 95 percent prediction intervals.}
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
######### Example of predictive score computation #############
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
upper <- list(mean = I, scale = I, sill = I)

fit <- GeoFit(data_est, coordx = coords_est,
 corrmodel = corrmodel, model = model,
 likelihood = "Marginal", type = "Pairwise",
 neighb = 3, optimizer = "nlminb",
 lower = lower, upper = upper,
 start = start, fixed = fixed)

# Prediction at validation locations
pr <- GeoKrig(fit,loc = coords_to_pred, data = data_est, mse = TRUE)

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

# Conditional-simulation objects can be scored directly. For example,
# after obtaining cs <- GeoSimcond(..., nrep = 500), use
# GeoScores(data_to_pred, cs)
# to compute point scores, empirical CRPS, empirical 95 percent interval
# score, and empirical coverage from cs$condsim.

}

\keyword{Predictive scores}
