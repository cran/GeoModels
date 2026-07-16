\name{GeoCV}
\alias{GeoCV}
\encoding{UTF-8}
\title{n-fold kriging Cross-validation}
\description{
The procedure uses the \code{\link{GeoKrig}} or \code{\link{GeoKrigloc}}
function to compute n-fold kriging cross-validation using information from a
\code{\link{GeoFit}} object. The function returns prediction scores.
}

\usage{
GeoCV(
  fit, K = 100, estimation = TRUE,
  optimizer = NULL, lower = NULL, upper = NULL,
  n.fold = 0.05, local = FALSE,
  neighb = NULL, maxdist = NULL, maxtime = NULL,
  sparse = FALSE, type_krig = "Simple", which = 1,
  parallel = TRUE, ncores = NULL, progress = TRUE,
  seed = NULL
)
}

\arguments{
 \item{fit}{An object of class \code{\link{GeoFit}}.}

 \item{K}{The number of iterations in cross-validation.}

 \item{estimation}{Logical; if \code{TRUE}, an estimation is performed at each
 iteration and the estimates are used in the prediction. Otherwise, the estimates
 in the object \code{fit} are used.}

 \item{optimizer}{The type of optimization algorithm if \code{estimation} is
 \code{TRUE}. See \code{\link{GeoFit}} for details. If \code{NULL}, then the
 optimization algorithm stored in \code{fit} is used.}

 \item{lower}{An optional named list giving the values for the lower bounds of
 the parameters when bounded optimization is used and \code{estimation} is
 \code{TRUE}.}

 \item{upper}{An optional named list giving the values for the upper bounds of
 the parameters when bounded optimization is used and \code{estimation} is
 \code{TRUE}.}

 \item{n.fold}{Numeric; the percentage of data to be deleted and predicted in
 each cross-validation iteration.}

 \item{local}{Logical; if \code{TRUE}, local kriging is performed. The default
 is \code{FALSE}.}

 \item{neighb}{Numeric; an optional positive integer indicating the order of
 neighborhood if local kriging is performed.}

 \item{maxdist}{Numeric; an optional positive value indicating the spatial
 neighborhood distance if local kriging is performed.}

 \item{maxtime}{Numeric; an optional positive value indicating the temporal
 neighborhood distance if local kriging is performed.}

 \item{sparse}{Logical; if \code{TRUE}, kriging and simulation are computed with
 sparse matrix algorithms using the \pkg{spam} package. The default is
 \code{FALSE}. It should be used with compactly supported covariance models.}

 \item{type_krig}{String; the type of kriging. If \code{"Simple"}, the default,
 simple kriging is performed. If \code{"Optim"}, optimal kriging is performed for
 some non-Gaussian random fields.}

 \item{which}{Numeric; in the case of bivariate cokriging, it indicates which
 variable to predict. It can be \code{1} or \code{2}.}

 \item{parallel}{Logical; if \code{TRUE}, the cross-validation iterations are
 parallelized.}

 \item{ncores}{Numeric; number of cores involved in parallelization.}

 \item{progress}{Logical; if \code{TRUE}, a progress bar is shown.}

 \item{seed}{Optional integer seed used to make the random selection of folds and
 any stochastic refitting step reproducible. If \code{NULL}, the current random
 number generator state is left unchanged. When supplied, the previous
 \code{.Random.seed} is restored on exit.}
}

\details{
The function randomly removes a fraction \code{n.fold} of the observations at
each iteration, predicts the removed observations using kriging, and computes a
set of predictive scores.

If \code{estimation = TRUE}, the model is re-estimated at each cross-validation
iteration before prediction. Otherwise, the parameter estimates contained in
\code{fit} are used.

When \code{seed} is supplied, the cross-validation folds are reproducible. The
function preserves and restores the user's random number generator state.
}

\value{
Returns a list containing the following components:
 \item{rmse}{The vector of root mean squared errors.}
 \item{mae}{The vector of mean absolute errors.}
 \item{mad}{The vector of median absolute errors.}
 \item{brie}{The vector of Brier scores.}
 \item{crps}{The vector of continuous ranked probability scores.}
 \item{lscore}{The vector of log-scores.}
 \item{pit}{The vector of average probability integral transform values.}
 \item{intscore}{The vector of interval scores.}
 \item{coverage}{The vector of empirical coverage values.}
 \item{seed}{The seed used for reproducibility, or \code{NULL}.}
}

\seealso{\code{\link{GeoKrig}}, \code{\link{GeoKrigloc}}, \code{\link{GeoFit}}.}

\author{
Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},
\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl},
\url{https://sites.google.com/site/moralesonatevictor/},
Christian Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},
\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{
\donttest{
library(GeoModels)

###############################################################
### Example of spatial kriging cross-validation
###############################################################

model <- "Gaussian"
set.seed(79)
x <- runif(400, 0, 1)
y <- runif(400, 0, 1)
coords <- cbind(x, y)

corrmodel <- "GenWend"
mean <- 0
sill <- 5
nugget <- 0
scale <- 0.2
smooth <- 0
power2 <- 4

param <- list(
  mean = mean, sill = sill, nugget = nugget,
  scale = scale, smooth = smooth, power2 = power2
)

data <- GeoSim(coordx = coords, corrmodel = corrmodel,
               param = param)$data

fixed <- list(nugget = nugget, smooth = 0, power2 = power2)
start <- list(mean = 0, scale = scale, sill = 1)
I <- Inf
lower <- list(mean = -I, scale = 0, sill = 0)
upper <- list(mean = I, scale = I, sill = I)

fit <- GeoFit(
  data, coordx = coords, corrmodel = corrmodel,
  model = model, likelihood = "Marginal", type = "Pairwise",
  neighb = 3, optimizer = "nlminb", lower = lower,
  upper = upper, start = start, fixed = fixed
)

#a <- GeoCV(fit, K = 100, estimation = TRUE,
#           parallel = TRUE, seed = 123)
#mean(a$rmse)
}
}

\keyword{Composite}
