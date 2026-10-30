# Internal diagnostic interface for the direct RF density (not a copula density).
.GeoSkewLaplaceLogDensity <- function(x, y, rho, skew, sill = 1, log = TRUE) {
  x <- as.double(x); y <- as.double(y); rho <- as.double(rho)
  if (!(length(y) %in% c(1L, length(x))) ||
      !(length(rho) %in% c(1L, length(x))))
    stop("y and rho must have length one or length(x).", call. = FALSE)
  if (length(skew) != 1L || !is.finite(skew) || skew <= 0 || skew >= 1 ||
      length(sill) != 1L || !is.finite(sill) || sill <= 0 ||
      any(!is.finite(c(x, y, rho))) || any(abs(rho) >= 1))
    stop("Invalid direct SkewLaplace density arguments.", call. = FALSE)
  ans <- .Call("geo_skewlaplace_logdensity", x, y, rho,
               as.double(skew), as.double(sill), PACKAGE = "GeoModels")
  if (isTRUE(log)) ans else exp(ans)
}

# Coordinate identity map: matching observations are conditioned exactly;
# repeated prediction sites share the same random value, rather than a new draw.
.GeoSimcondPredictionMap <- function(coord_obs, loc) {
  key <- function(z) apply(z, 1L, function(r)
    paste(sprintf("%.17g", r), collapse = "|"))
  ko <- key(coord_obs); kl <- key(loc)
  obs_match <- match(kl, ko)
  unknown <- which(is.na(obs_match))
  unique_positions <- unknown[!duplicated(kl[unknown])]
  list(obs_match = obs_match, unknown = unknown,
       unique_positions = unique_positions,
       loc_map = match(kl, kl[unique_positions]))
}

# Prepare a dense GLOBAL Gaussian conditional law once per GeoSimcond call.
# Cholesky of observations is shared by the solve and (if needed) the precision.
# Only the prediction conditional covariance is factored for new simulations.
.GeoSimcondGaussianPlan <- function(coord_obs, loc, corrmodel, param,
                                    distance, radius, anisopars,
                                    precision = FALSE) {
  nobs <- nrow(coord_obs); nloc <- nrow(loc)
  cov <- GeoCovmatrix(
    coordx = rbind(coord_obs, loc), corrmodel = corrmodel,
    model = "Gaussian", param = param, distance = distance, radius = radius,
    anisopars = anisopars, sparse = FALSE, copula = NULL,
    check.duplicates = FALSE
  )$covmatrix
  cov <- as.matrix(cov)
  if (!identical(dim(cov), c(nobs + nloc, nobs + nloc)) ||
      any(!is.finite(cov)))
    stop("Invalid joint latent Gaussian covariance.", call. = FALSE)
  cov <- 0.5 * (cov + t(cov))
  io <- seq_len(nobs); il <- nobs + seq_len(nloc)
  Roo <- cov[io, io, drop = FALSE]
  Rol <- cov[io, il, drop = FALSE]
  Rll <- cov[il, il, drop = FALSE]
  rm(cov)
  U <- tryCatch(chol(Roo), error = function(e) NULL)
  if (is.null(U))
    stop("Latent observation covariance is not positive definite; check duplicate sites and covariance parameters.",
         call. = FALSE)
  Wobs <- backsolve(U, forwardsolve(t(U), Rol))
  if (any(!is.finite(Wobs)))
    stop("Non-finite latent Gaussian conditional weights.", call. = FALSE)
  C <- Rll - crossprod(Rol, Wobs)
  C <- 0.5 * (C + t(C))
  V <- tryCatch(chol(C), error = function(e) NULL)
  factor_method <- "Cholesky"
  if (!is.null(V)) {
    L <- t(V)
  } else {
    # A semidefinite conditional law can occur at a numerical rank boundary.
    # Do not add an artificial nugget. Clip ONLY roundoff-scale negative roots.
    ee <- eigen(C, symmetric = TRUE)
    tol <- 100 * .Machine$double.eps * max(1L, nloc) * max(1, max(abs(Rll)))
    if (min(ee$values) < -tol)
      stop("Latent prediction conditional covariance is not positive semidefinite.",
           call. = FALSE)
    L <- sweep(ee$vectors, 2L, sqrt(pmax(ee$values, 0)), "*")
    factor_method <- "PSD eigenfactor (roundoff clipping only)"
  }
  Q <- if (isTRUE(precision)) chol2inv(U) else NULL
  if (!is.null(Q)) {
    Q <- 0.5 * (Q + t(Q))
    if (any(!is.finite(Q)) || any(diag(Q) <= 0))
      stop("Invalid latent Gaussian precision matrix.", call. = FALSE)
  }
  list(W = t(Wobs), L = L, Q = Q, conditional_factor = factor_method)
}

# Exact global Gaussian simulation for the SkewLaplace Gaussian-copula route.
# Does not use a Gibbs sampler or the direct RF's four latent fields.
.GeoSimcondGaussianExact <- function(data, coord_obs, loc, corrmodel, param,
                                     nrep, distance, radius, anisopars,
                                     progress = FALSE) {
  coord_obs <- as.matrix(coord_obs); loc <- as.matrix(loc)
  mp <- .GeoSimcondPredictionMap(coord_obs, loc)
  nloc <- nrow(loc); nu <- length(mp$unique_positions)
  observed <- which(!is.na(mp$obs_match))
  base <- numeric(nloc)
  base[observed] <- data[mp$obs_match[observed]]
  if (nu == 0L) return(rep(list(base), nrep))
  if (isTRUE(progress)) message("Gaussian copula: preparing the conditional Gaussian factor once ...")
  plan <- .GeoSimcondGaussianPlan(
    coord_obs, loc[mp$unique_positions, , drop = FALSE], corrmodel, param,
    distance, radius, anisopars, precision = FALSE
  )
  mu <- as.numeric(plan$W %*% data)
  batch <- max(1L, min(nrep, floor(64 * 1024^2 / (8 * nu))))
  ans <- vector("list", nrep)
  for (first in seq.int(1L, nrep, by = batch)) {
    last <- min(nrep, first + batch - 1L)
    nb <- last - first + 1L
    z <- plan$L %*% matrix(stats::rnorm(nu * nb), nrow = nu, ncol = nb)
    z <- sweep(z, 1L, mu, "+")
    for (j in seq_len(nb)) {
      value <- base
      value[mp$unknown] <- z[mp$loc_map[mp$unknown], j]
      ans[[first + j - 1L]] <- value
    }
  }
  attr(ans, "simulation_details") <- list(
    engine = "Gaussian copula: cached exact Gaussian conditional factor",
    observation_factorizations = 1L, prediction_factorizations = 1L,
    conditional_factor = plan$conditional_factor, mcmc = FALSE
  )
  ans
}
