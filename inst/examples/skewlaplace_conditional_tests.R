## Targeted post-installation tests of both SkewLaplace conditional engines.
## Test 1: direct non-copula interpolation at observed sites.
## Test 2: exact conditional mean for one observed location.
## Test 3: Gaussian-copula mode still executes and returns finite draws.
library(GeoModels)
set.seed(303)
p <- 0.75
sigma <- 1.2
eta1 <- 0.3
eta2 <- -0.2
r <- 0.6

## A two-site pure Gaussian correlation model; the exact cross covariance is r.
## Use the Gaussian covariance matrix to compute the latent correlation,
## avoiding any assumption about the numeric interpretation of 'scale'.
obs <- matrix(c(0.0, 0.0), nrow = 1)
pred <- matrix(c(0.25, 0.0), nrow = 1)
param <- list(mean = eta1, skew = p, sill = sigma^2,
              nugget = 0, scale = 0.5)
gpars <- list(mean = 0, sill = 1, nugget = 0, scale = 0.5)
R2 <- GeoCovmatrix(coordx = rbind(obs, pred), corrmodel = "Exp",
                   model = "Gaussian", param = gpars)$covmatrix
r <- as.numeric(R2[1, 2])

## Choose a response value deliberately away from its unconditional mean.
yobs <- 1.8
nrep <- 2500
set.seed(11)
direct <- GeoSimcond(data = yobs, coordx = obs, loc = pred,
                     corrmodel = "Exp", model = "SkewLaplace",
                     param = param, Mloc = eta2, copula = NULL,
                     nrep = nrep, n_iter = 30, mcmc_thin = 1)

## With exactly ONE conditioning site, the direct AL-RF has the identity:
## E(Y2 | Y1=y) = eta2 + r^2*(y-eta1) +
##   (1-r^2)*sigma*(1/p-1/(1-p)).
mean_exact <- eta2 + r^2 * (yobs - eta1) +
              (1-r^2) * sigma * (1/p - 1/(1-p))
mean_mc <- mean(vapply(direct$condsim, function(x) x[1], numeric(1)))
cat(sprintf("Direct conditional E[Y2|Y1]: MC = %.4f ; theory = %.4f\n",
            mean_mc, mean_exact))
stopifnot(abs(mean_mc - mean_exact) < 0.5)

## This model also has an explicit one-site conditional variance formula.
y_standard <- (yobs - eta1) / sigma
r2 <- r^2
v <- 1 - r2
var_exact <- sigma^2 * (
  v^2 * (1/p^2 + 1/(1-p)^2) +
  2 * r2 * v * ((1 + max(y_standard, 0))/p +
                (1 + max(-y_standard, 0))/(1-p))
)
var_mc <- var(vapply(direct$condsim, function(x) x[1], numeric(1)))
cat(sprintf("Direct conditional Var[Y2|Y1]: MC = %.4f ; theory = %.4f\n",
            var_mc, var_exact))
stopifnot(abs(var_mc - var_exact) < max(2.5, 0.2 * var_exact))

## Exact interpolation: prediction at training site is deterministic.
same_site <- GeoSimcond(data = yobs, coordx = obs, loc = obs,
                        corrmodel = "Exp", model = "SkewLaplace",
                        param = param, Mloc = eta1,
                        copula = NULL, nrep = 4)
stopifnot(all(vapply(same_site$condsim,
                     function(x) abs(x[1] - yobs) < 1e-11, logical(1))))

## Copula Gaussian case must remain available and produce finite values.
set.seed(12)
cop <- GeoSimcond(data = yobs, coordx = obs, loc = pred,
                  corrmodel = "Exp", model = "SkewLaplace",
                  param = param, Mloc = eta2,
                  copula = "Gaussian", nrep = 500)
stopifnot(length(cop$condsim) == 500, all(is.finite(unlist(cop$condsim))))
cat("BOTH SkewLaplace conditional modes: checks passed.\n")
