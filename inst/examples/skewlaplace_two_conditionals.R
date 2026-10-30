## GeoModels 2.3.0 - SkewLaplace direct vs Gaussian-copula conditional simulation
## Reproducible illustrative example. Run AFTER installing the patched package.
## This deliberately uses the SAME observed data for the two conditional models:
## the Gaussian copula is then a competing predictive assumption, not the true
## data-generating model.
library(GeoModels)

set.seed(20261027)
NN <- 100L
coords <- cbind(runif(NN), runif(NN))
X <- cbind(1, runif(NN, -1, 1))
loc <- cbind(runif(12), runif(12))
Xloc <- cbind(1, runif(nrow(loc), -1, 1))

param <- list(mean = 0.1, mean1 = -0.25, skew = 0.75,
              sill = 1, nugget = 0, scale = 0.15,
              smooth = 0, power2 = 4)

simulated <- GeoSim(coordx = coords, X = X, corrmodel = "GenWend",
                    model = "SkewLaplace", param = param)
data <- as.numeric(simulated$data)

## (1) DIRECT AL random field of the paper (four Gaussian latent fields):
set.seed(123)
direct <- GeoSimcond(data = data, coordx = coords, loc = loc, X = X,
                     Xloc = Xloc, corrmodel = "GenWend", param = param,
                     model = "SkewLaplace", copula = NULL,
                     nrep = 250, n_iter = 1000, mcmc_thin = 5)

## (2) Gaussian copula with the very same AL marginal distribution:
set.seed(123)
cop <- GeoSimcond(data = data, coordx = coords, loc = loc, X = X,
                  Xloc = Xloc, corrmodel = "GenWend", param = param,
                  model = "SkewLaplace", copula = "Gaussian", nrep = 250)

comparison <- data.frame(
  s1 = loc[, 1], s2 = loc[, 2],
  direct_mean = direct$cond_mean,
  gaussian_copula_mean = cop$cond_mean,
  direct_q25 = apply(do.call(rbind, direct$condsim), 2, quantile, 0.25),
  copula_q25 = apply(do.call(rbind, cop$condsim), 2, quantile, 0.25)
)
print(comparison)

## Interpolation test for the DIRECT model at observed locations:
## Even with different covariates, the conditional residual is fixed exactly.
matched <- GeoSimcond(data = data, coordx = coords, loc = coords[1:4, , drop = FALSE],
                       X = X, Xloc = X[1:4, , drop = FALSE],
                       corrmodel = "GenWend", param = param,
                       model = "SkewLaplace", copula = NULL, nrep = 3)
stopifnot(all(vapply(matched$condsim,
                     function(z) isTRUE(all.equal(z, data[1:4], tolerance = 1e-10)),
                     logical(1))))
