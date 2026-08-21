GeoKrigWeights <- function(
  coordx, coordy = NULL, coordz = NULL, coordt = NULL, coordx_dyn = NULL,
  corrmodel, distance = "Eucl", grid = FALSE, loc,
  method = "cholesky", model = "Gaussian", n = 1, nloc = NULL,
  param, anisopars = NULL, radius = 1, sparse = FALSE, time = NULL,
  which = 1, copula = NULL, X = NULL, Xloc = NULL, Mloc = NULL
) {
  ############################################################################
  ## GeoKrigWeights: global kriging weights
  ##
  ## Conservative updated version.
  ## Main reliable target: Gaussian latent/global weights.
  ## The default orientation is kept as nobs x nloc for compatibility with
  ## existing GeoSimcond code using crossprod(W, x).
  ############################################################################

  call <- match.call()

  ## Internal defaults kept for compatibility with previous internal code.
  ## They are intentionally not formal arguments, so the public API matches
  ## the existing Rd documentation.
  return_covmatrix <- TRUE
  return_CC <- TRUE
  weights_orientation <- "nobs_by_nloc"

  method <- gsub("[[:blank:]]", "", method)
  model <- .GeoPredictionCanonicalModel(gsub("[[:blank:]]", "", model))
  corrmodel <- gsub("[[:blank:]]", "", corrmodel)
  distance <- gsub("[[:blank:]]", "", distance)

  if (is.null(CkModel(model))) {
    stop("The name of the model is not correct\n")
  }
  corrmodel_code <- CkCorrModel(corrmodel)
  if (!is.character(corrmodel) || is.null(corrmodel_code)) {
    stop("the name of the correlation model is wrong")
  }

  bivariate <- CheckBiv(corrmodel_code)
  spacetime <- CheckST(corrmodel_code)

  ## This version is intentionally conservative.  GeoSimcond uses this function
  ## for Gaussian latent weights, so do not silently return approximate/wrong
  ## weights for models whose cross-covariance requires model-specific code.
  if (!identical(model, "Gaussian")) {
    stop("GeoKrigWeights currently returns validated weights only for model = 'Gaussian'. ",
         "For non-Gaussian prediction use GeoKrig() or GeoKrigloc().", call. = FALSE)
  }
  if (!is.null(copula)) {
    stop("GeoKrigWeights currently requires copula = NULL; conditional-simulation ",
         "code must request weights on the latent Gaussian scale.", call. = FALSE)
  }
  if (bivariate && !(which == 1 || which == 2)) {
    stop("which parameter must be 1 or 2")
  }

  if (spacetime && !is.null(coordx_dyn)) {
    coordx_dyn <- .GeoNormalizeDynamicCoordinates(
      coordx_dyn, coordt = coordt, context = "GeoKrigWeights"
    )
    .GeoValidateDynamicRows(X, coordx_dyn, name = "X",
                            context = "GeoKrigWeights")
  }

  biv_info <- NULL
  if (bivariate) {
    biv_info <- .GeoBivariateCoordinateInfo(
      coordx = coordx, coordy = coordy, coordz = coordz,
      coordx_dyn = coordx_dyn, grid = grid,
      context = "GeoKrigWeights"
    )
    if (biv_info$dynamic) {
      coordx_dyn <- biv_info$blocks
      coordx <- coordy <- coordz <- NULL
    }
    if (!is.null(X)) {
      X <- .GeoBivariateDesign(X, biv_info$ns,
                               context = "GeoKrigWeights")
    }
  }

  ## Location matrix
  if (is.vector(loc)) {
    loc <- matrix(loc, nrow = 1)
  } else if (!is.matrix(loc)) {
    loc <- as.matrix(loc)
  }
  if (!(ncol(loc) %in% c(2L, 3L))) {
    stop("loc parameter must be a matrix N x 2 or N x 3")
  }


  .GeoValidate3D(
    coordx = coordx, coordy = coordy, coordz = coordz,
    coordx_dyn = coordx_dyn, distance = distance, grid = grid,
    context = "GeoKrigWeights"
  )
  .GeoValidate3D(
    coordx = loc, distance = distance, grid = FALSE,
    context = "GeoKrigWeights prediction locations"
  )
  .GeoValidateCoordinateCompatibility(
    coordx = coordx, coordy = coordy, coordz = coordz,
    coordx_dyn = coordx_dyn, loc = loc, context = "GeoKrigWeights"
  )

  if (spacetime && is.null(time)) {
    stop("At least one temporal instant is needed for space-time kriging\n")
  }
  if (is.null(time)) time <- 0

  ## Mean inputs are normalized after dimat/dimat2 are known.
  if (!is.null(X) && !is.matrix(X) && !is.list(X)) X <- as.matrix(X)
  Mobs <- if (!bivariate) .GeoMean_external(param) else NULL
  if (!is.null(Mobs) && !is.null(X))
    stop("A fixed observation mean and X cannot be used together.")
  if (!is.null(Mloc) && !is.null(Xloc))
    stop("Use either Mloc or Xloc, not both.")
  if (!is.null(Mobs) && is.null(Mloc))
    stop("Mloc is required when the observation mean is supplied as a vector.")
  if (!is.null(X) && is.null(Xloc) && is.null(Mloc))
    stop("Xloc or Mloc is required when X is supplied.")
  if (is.null(X) && !is.null(Xloc))
    stop("X must be supplied when Xloc is used.")
  if (is.null(X) && is.null(Mobs) &&
      length(setdiff(.GeoMean_supplied_names(param), "mean")))
    stop("Parameters mean1, mean2, ... require an explicit design matrix X; without X the mean is intercept-only.")

  ## Prediction locations in the metric used by Corr_c.
  ## GeoCovmatrix handles the observed covariance.  For the cross-covariance we
  ## apply the same anisotropic map to ccc and loc below.  If the caller already
  ## passes anisotropized coordinates, anisopars must be NULL.
  loc_metric <- loc
  if (!is.null(anisopars)) {
    loc_metric <- GeoAniso(loc_metric, c(anisopars$angle, anisopars$ratio))
  }

  numloc <- nrow(loc_metric)
  tloc <- length(time)
  if (!tloc) tloc <- 1L
  if (bivariate) tloc <- 1L

  if (ncol(loc_metric) == 2L) {
    locx <- loc_metric[, 1L]
    locy <- loc_metric[, 2L]
    locz <- rep(0, numloc)
  } else {
    locx <- loc_metric[, 1L]
    locy <- loc_metric[, 2L]
    locz <- loc_metric[, 3L]
  }

  ############################################################################
  ## Build covariance matrix at observed sites
  ############################################################################
  param_cov <- param
  if (!is.null(Mobs)) param_cov$mean <- 0
  covmatrix <- GeoCovmatrix(
    coordx = coordx, coordy = coordy, coordz = coordz, coordt = coordt,
    coordx_dyn = coordx_dyn, corrmodel = corrmodel, distance = distance,
    grid = grid, model = model, n = n, param = param_cov,
    anisopars = anisopars, radius = radius, sparse = sparse,
    copula = copula, X = X
  )
  covmatrix$param <- unlist(covmatrix$param)

  spacetime_dyn <- !is.null(covmatrix$coordx_dyn)
  if (!spacetime_dyn) {
    dimat <- covmatrix$numcoord * covmatrix$numtime
  } else {
    dimat <- sum(covmatrix$ns)
  }
  dimat2 <- numloc * tloc

  if (!bivariate) {
    if (is.null(X)) X <- covmatrix$X
    X <- .GeoMean_design(X, n = dimat, name = "X")
    p_mean <- ncol(X)
    invisible(.GeoMean_beta(param_cov, p_mean))
    if (!is.null(Mobs) && length(Mobs) != dimat)
      stop("The fixed observation mean must have one value per observation.")

    if (!is.null(Mloc)) {
      if (!is.numeric(Mloc) || any(!is.finite(Mloc)))
        stop("Mloc must be finite and numeric.")
      Mloc <- as.numeric(Mloc)
      if (length(Mloc) == 1L) Mloc <- rep(Mloc, dimat2)
      if (length(Mloc) != dimat2)
        stop("Mloc must have one value per prediction point/time combination.")
      Xloc <- NULL
    } else {
      Xloc <- .GeoMean_design(Xloc, n = dimat2, p = p_mean, name = "Xloc")
    }
  } else {
    X <- .GeoBivariateDesign(
      if (is.null(X)) covmatrix$X else X,
      covmatrix$ns, context = "GeoKrigWeights"
    )
    p_mean <- ncol(X)
    invisible(.GeoBivariateValidateMeanParameters(
      covmatrix$param, p_mean, context = "GeoKrigWeights"
    ))
    if (!is.null(Mloc)) {
      if (!is.numeric(Mloc) || any(!is.finite(Mloc)))
        stop("Mloc must be finite and numeric.")
      Mloc <- as.numeric(Mloc)
      if (length(Mloc) == 1L) Mloc <- rep(Mloc, dimat2)
      if (length(Mloc) != dimat2)
        stop("Mloc must have one value per prediction location.")
      Xloc <- NULL
    } else {
      Xloc <- .GeoMean_design(Xloc, n = dimat2, p = p_mean, name = "Xloc")
    }
  }

  ############################################################################
  ## Observed coordinates stacked as required by Corr_c
  ############################################################################
  NS <- 0L
  if (spacetime || bivariate) {
    NS <- cumsum(covmatrix$ns)
    NS <- c(0L, NS)[-(length(covmatrix$ns) + 1L)]
  }

  if (!grid) {
    if (!spacetime && !bivariate) {
      if (is.null(covmatrix$coordz)) {
        ccc <- cbind(covmatrix$coordx, covmatrix$coordy, 0)
      } else {
        ccc <- cbind(covmatrix$coordx, covmatrix$coordy, covmatrix$coordz)
      }
    }
  }

  if (grid) {
    if (is.null(coordz)) {
      ccc <- as.matrix(expand.grid(covmatrix$coordx, covmatrix$coordy))
      ccc <- cbind(ccc, 0)
    } else {
      ccc <- as.matrix(expand.grid(covmatrix$coordx, covmatrix$coordy, covmatrix$coordz))
    }
    grid <- FALSE
  } else {
    if ((spacetime || bivariate) && !spacetime_dyn) {
      if (!is.null(covmatrix$coordz)) {
        ccc <- cbind(rep(covmatrix$coordx, covmatrix$numtime),
                     rep(covmatrix$coordy, covmatrix$numtime),
                     rep(covmatrix$coordz, covmatrix$numtime))
      } else {
        ccc <- cbind(rep(covmatrix$coordx, covmatrix$numtime),
                     rep(covmatrix$coordy, covmatrix$numtime), 0)
      }
    }
    if ((spacetime || bivariate) && spacetime_dyn) {
      ccc <- do.call(rbind, args = c(covmatrix$coordx_dyn))
      if (ncol(ccc) == 2L) ccc <- cbind(ccc, 0)
    }
  }

  if (!exists("ccc")) stop("Could not build observed coordinate matrix")
  ccc <- as.matrix(ccc)
  if (ncol(ccc) == 2L) ccc <- cbind(ccc, 0)
  if (nrow(ccc) != dimat) stop("Observed coordinate matrix has wrong number of rows")

  if (!is.null(anisopars)) {
    ccc <- GeoAniso(ccc, c(anisopars$angle, anisopars$ratio))
  }

  ############################################################################
  ## Cross-correlation and cross-covariance CC
  ############################################################################
  corrmodel_num <- CkCorrModel(covmatrix$corrmodel)
  distance_num <- CheckDistance(covmatrix$distance)
  corrparam <- covmatrix$param[covmatrix$namescorr]

  if (!(covmatrix$model %in% c(1))) {
    stop("Internal error: GeoKrigWeights is restricted to Gaussian model code.", call. = FALSE)
  }

  cc_call <- dotCall64::.C64(
    "Corr_c",
    SIGNATURE = c("double", "double", "double", "double", "double", "integer",
                  "integer", "double", "double", "double", "integer", "integer",
                  "integer", "integer", "integer", "integer", "double", "integer",
                  "integer", "double", "integer", "integer", "double"),
    corri = dotCall64::vector_dc("double", dimat * dimat2),
    ccc[, 1L], ccc[, 2L], ccc[, 3L], covmatrix$coordt,
    corrmodel_num, as.integer(0),
    locx, locy, locz,
    covmatrix$numcoord, numloc, tloc,
    covmatrix$ns, NS, covmatrix$numtime,
    corrparam,
    covmatrix$spacetime, covmatrix$bivariate,
    time, distance_num, as.integer(which - 1L), covmatrix$radius,
    INTENT = c("w", rep("r", 22)),
    NAOK = TRUE, PACKAGE = "GeoModels"
  )

  rho <- cc_call$corri

  if (bivariate) {
    ## Corr_c evaluates the bivariate covariance model directly.  The returned
    ## vector already includes marginal sills and cross-covariance scaling.
    corri <- rho
  } else {
    nug <- as.numeric(covmatrix$param["nugget"])
    sill <- as.numeric(covmatrix$param["sill"])
    if (is.na(nug)) nug <- 0
    if (is.na(sill)) stop("sill parameter is missing")
    corri <- (1 - nug) * rho
    vvar <- sill
  }

  CC <- if (bivariate) {
    matrix(corri, nrow = dimat, ncol = dimat2)
  } else {
    matrix(corri * as.numeric(vvar), nrow = dimat, ncol = dimat2)
  }
  if (nrow(CC) != dimat || ncol(CC) != dimat2) {
    stop("Wrong dimension of cross-covariance matrix CC")
  }

  ############################################################################
  ## Final weights
  ############################################################################
  weights_result <- getInvC(covmatrix, CC, mse = FALSE)
  krig_weights <- weights_result$a

  ## Default is the historical orientation used by GeoSimcond:
  ## kriging vector = crossprod(weights, data), i.e. weights is nobs x nloc.
  if (weights_orientation == "nloc_by_nobs") {
    krig_weights <- t(krig_weights)
  }

  result <- list(
    weights = krig_weights,
    weights_orientation = weights_orientation,
    model = model,
    corrmodel = corrmodel,
    coordx = covmatrix$coordx,
    coordy = covmatrix$coordy,
    coordz = covmatrix$coordz,
    coordx_dyn = covmatrix$coordx_dyn,
    bivariate = bivariate,
    spacetime = covmatrix$spacetime,
    param = covmatrix$param,
    xloc = Xloc,
    Mloc = Mloc,
    mean_obs = Mobs,
    tloc = tloc,
    loc = loc,
    loc_metric = loc_metric,
    call = call
  )

  if (return_covmatrix) result$covmatrix <- covmatrix$covmatrix
  if (return_CC) result$CC <- CC

  return(result)
}