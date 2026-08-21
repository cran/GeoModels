####################################################
### File name: GeoCorrFct_Cop.R
####################################################

GeoCorrFct_Cop <- function(x, t = NULL, corrmodel, model = "Gaussian", copula = "Gaussian",
                           distance = "Eucl", param, radius = 6371, n = 1,
                           covariance = FALSE, variogram = FALSE)
{
  ############################################################################
  ## Latent Gaussian correlation from the native correlation engine
  ############################################################################
  CorrelationFct <- function(corrmodel, lags, lagt, numlags, numlagt,
                             mu, paramcorr, N) {
    p <- dotCall64::.C64(
      "VectCorrelation",
      SIGNATURE = c("double", "integer", "double", "integer", "integer",
                    "double", "integer", "double", "double", "double", "integer"),
      corr = dotCall64::vector_dc("double", numlags * numlagt),
      corrmodel, lags, numlags, numlagt, mu, CkModel("Gaussian"),
      c(nugget = 0, sill = 1), paramcorr, lagt, N,
      INTENT = c("rw", rep("r", 10)),
      PACKAGE = "GeoModels", VERBOSE = 0, NAOK = TRUE
    )
    p$corr
  }

  ############################################################################
  ## Checks
  ############################################################################
  call <- match.call()

  corr_id <- CkCorrModel(corrmodel)
  model_id <- CkModel(model)
  if (is.null(corr_id))
    stop("The name of the correlation model is not correct\n")
  if (is.null(model_id))
    stop("The name of the model is not correct\n")
  if (!is.numeric(x) || any(!is.finite(x)))
    stop("Distances must be finite numeric values\n")
  if (any(x < 0))
    stop("Distances must be positive\n")
  if (!is.list(param))
    stop("param must be a list\n")

  copula <- as.character(copula)[1L]
  if (is.na(copula) || !(copula %in% c("Gaussian", "Clayton", "SkewGaussian"))) {
    stop("GeoCorrFct_Cop supports copula='Gaussian', 'Clayton', or 'SkewGaussian'.",
         call. = FALSE)
  }
  if (!(as.integer(model_id) %in% .GeoCopulaPredictionModels)) {
    stop("GeoCorrFct_Cop supports the validated copula prediction margins, including Poisson, Binomial, and BinomialNeg.",
         call. = FALSE)
  }

  spacetime <- CheckST(corr_id)
  bivariate <- CheckBiv(corr_id)
  if (bivariate)
    stop("GeoCorrFct_Cop is implemented only for univariate copula margins.", call. = FALSE)

  ## GeoCorrFct_Cop is a lag-only function and has no design matrix.  Therefore
  ## location-dependent margins can only use an intercept-only marginal predictor.
  supplied_mean <- .GeoMean_supplied_names(param)
  if (length(setdiff(supplied_mean, "mean"))) {
    stop("GeoCorrFct_Cop() has no design matrix X; only the intercept parameter 'mean' is allowed.",
         call. = FALSE)
  }

  param_work <- param
  if (is.null(param_work$mean)) param_work$mean <- 0
  if (length(param_work$mean) != 1L || !is.numeric(param_work$mean) ||
      !is.finite(param_work$mean)) {
    stop("GeoCorrFct_Cop() requires a finite scalar intercept 'mean'.", call. = FALSE)
  }
  ## Keep the same convenience convention as GeoCorrFct(): when a scale-based
  ## continuous margin needs sill and it is omitted, use unit sill.
  if (is.null(param_work$sill)) param_work$sill <- 1
  if (is.null(param_work$nugget)) param_work$nugget <- 0

  nugget <- suppressWarnings(as.numeric(param_work$nugget)[1L])
  if (!is.finite(nugget) || nugget < 0 || nugget >= 1)
    stop("GeoCorrFct_Cop requires nugget in [0,1).", call. = FALSE)

  if (identical(copula, "Clayton")) {
    nu <- suppressWarnings(as.numeric(param_work$nu)[1L])
    if (!is.finite(nu) || nu <= 0 ||
        abs(nu - round(nu)) > sqrt(.Machine$double.eps)) {
      stop("The GeoModels Clayton-like copula requires a finite positive integer 'nu'.",
           call. = FALSE)
    }
  }
  if (identical(copula, "SkewGaussian")) {
    eta_cop <- suppressWarnings(as.numeric(param_work$nu)[1L])
    if (!is.finite(eta_cop) || abs(eta_cop) >= 1)
      stop("The skew-Gaussian copula requires nu (eta) in (-1,1).", call. = FALSE)
  }

  ############################################################################
  ## Setup and latent correlation
  ############################################################################
  x <- as.numeric(x)
  nx <- length(x)
  if (is.null(t)) t <- 0
  if (!is.numeric(t) || any(!is.finite(t)) || any(t < 0))
    stop("Temporal distances must be finite non-negative numeric values.", call. = FALSE)
  t <- as.numeric(t)

  if (spacetime) {
    nt <- length(t)
    lag_grid <- as.matrix(expand.grid(x, t))
  } else {
    t <- 0
    nt <- 1L
    lag_grid <- cbind(x, 0)
  }

  corr_names <- CorrelationPar(corr_id)
  if (is.null(corr_names))
    stop("Correlation parameters are not defined for this correlation model.", call. = FALSE)
  missing_corr <- setdiff(corr_names, names(param_work))
  if (length(missing_corr)) {
    stop("Missing correlation parameter(s): ", paste(missing_corr, collapse = ", "), ".",
         call. = FALSE)
  }
  parcorr <- suppressWarnings(as.numeric(unlist(param_work[corr_names], use.names = FALSE)))
  if (length(parcorr) != length(corr_names) || any(!is.finite(parcorr)))
    stop("Correlation parameters must be finite numeric scalars.", call. = FALSE)

  eta <- as.numeric(param_work$mean)

  ## Pass only marginal/copula nuisance parameters to the covariance engine.
  ## Correlation parameters do not affect the copula transformation itself and
  ## including them would unnecessarily fragment the numerical caches.
  cop_names <- NuisParam2(model, FALSE, 1L, copula = copula)
  params_cop <- unlist(param_work[cop_names], use.names = TRUE)
  if (as.integer(model_id) %in% c(11L, 16L)) {
    nn <- suppressWarnings(as.numeric(n))
    if (length(nn) != 1L || !is.finite(nn) || nn < 1 ||
        abs(nn - round(nn)) > sqrt(.Machine$double.eps))
      stop("GeoCorrFct_Cop requires a positive integer n for Binomial/BinomialNeg margins.", call. = FALSE)
    params_cop["n"] <- nn
  }

  latent_rho <- CorrelationFct(corr_id, x, t, nx, nt, eta, parcorr, n)
  latent_rho <- as.numeric(latent_rho)
  if (length(latent_rho) != nrow(lag_grid) || any(!is.finite(latent_rho)))
    stop("Non-finite latent correlations returned by the correlation model.", call. = FALSE)
  if (any(abs(latent_rho) > 1 + 1e-8))
    stop("The latent correlation model returned values outside [-1,1].", call. = FALSE)
  latent_rho <- pmax(-1, pmin(1, latent_rho))

  ## The nugget acts on the latent correlation for distinct observations, exactly
  ## as in Cmatrix_copula().  At the exact zero lag GeoCorrFct_Cop represents the
  ## same random variable, hence its covariance is the marginal variance.
  rho_cop <- (1 - nugget) * latent_rho
  zero_lag <- lag_grid[, 1L] == 0 & lag_grid[, 2L] == 0
  rho_cop[zero_lag] <- 1

  ############################################################################
  ## Observed-scale covariance using the common three-copula dispatcher
  ############################################################################
  nout <- length(rho_cop)
  cova <- .GeoCopulaCrossCovFast(
    copula = copula,
    rho = rho_cop,
    model_type = model_id,
    params = params_cop,
    eta1 = rep(eta, nout),
    eta2 = rep(eta, nout)
  )

  moments <- .GeoCopulaMarginalMoments(model_id, params_cop, eta)
  marginal_variance <- as.numeric(moments$variance)[1L]
  if (!is.finite(marginal_variance) || marginal_variance <= 0)
    stop("Invalid marginal variance in GeoCorrFct_Cop().", call. = FALSE)

  ## Enforce the exact same-variable endpoint.  The fast numerical engines are
  ## calibrated to this limit, but setting it explicitly avoids tiny quadrature
  ## or truncation errors in plotted correlation/variogram curves.
  cova[zero_lag] <- marginal_variance
  corr_obs <- cova / marginal_variance
  corr_obs[zero_lag] <- 1

  if (isTRUE(variogram)) {
    if (isTRUE(covariance)) {
      res <- marginal_variance - cova
    } else {
      res <- 1 - corr_obs
    }
    res[zero_lag] <- 0
  } else if (isTRUE(covariance)) {
    res <- cova
  } else {
    res <- corr_obs
  }

  out <- list(
    corr = res,
    distances = x,
    times = t,
    model = model,
    copula = copula,
    distance = distance,
    param = param_work,
    radius = radius,
    n = n,
    covariance = covariance,
    variogram = variogram,
    spacetime = spacetime,
    bivariate = FALSE
  )

  structure(c(out, call = call), class = c("GeoCorrFct"))
}
