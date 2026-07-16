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
  model <- gsub("[[:blank:]]", "", model)
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
    stop("GeoKrigWeights_updated currently supports model = 'Gaussian' reliably. ",
         "Use the original full model-specific version for non-Gaussian weights, ",
         "or extend this function model by model.")
  }
  if (!is.null(copula)) {
    stop("GeoKrigWeights_updated currently supports copula = NULL only. ",
         "For conditional simulation with copulas, call it on the latent Gaussian scale.")
  }
  if (bivariate && !(which == 1 || which == 2)) {
    stop("which parameter must be 1 or 2")
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

  if (spacetime && is.null(time)) {
    stop("At least one temporal instant is needed for space-time kriging\n")
  }
  if (is.null(time)) time <- 0

  ## X/Xloc checks.  They are mostly irrelevant for simple Gaussian covariance
  ## weights, but keeping them coherent avoids hidden dimension bugs when the
  ## object is passed downstream.
  if (!is.null(X) && !is.matrix(X)) X <- as.matrix(X)
  if (!is.null(Xloc)) {
    if (is.vector(Xloc)) Xloc <- matrix(Xloc, nrow = 1) else Xloc <- as.matrix(Xloc)
  }

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
  covmatrix <- GeoCovmatrix(
    coordx = coordx, coordy = coordy, coordz = coordz, coordt = coordt,
    coordx_dyn = coordx_dyn, corrmodel = corrmodel, distance = distance,
    grid = grid, model = model, n = n, param = param,
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

  if (is.null(X)) {
    X <- covmatrix$X
  } else {
    X <- as.matrix(X)
    if (nrow(X) != dimat) stop("X has wrong number of rows")
  }

  if (is.null(Xloc)) {
    if (!is.null(X) && ncol(as.matrix(X)) > 1L) {
      stop("Xloc is missing but X has more than one column")
    }
    Xloc <- matrix(1, nrow = dimat2, ncol = 1L)
  } else {
    Xloc <- as.matrix(Xloc)
    if (nrow(Xloc) != dimat2) stop("Xloc has wrong number of rows")
    if (!is.null(X) && ncol(Xloc) != ncol(as.matrix(X))) {
      stop("Xloc and X must have the same number of columns")
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
    stop("Internal error: this updated function is restricted to Gaussian model code")
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
    vvar <- if (which == 1L) {
      covmatrix$param["sill_1"] + covmatrix$param["nugget_1"]
    } else {
      covmatrix$param["sill_2"] + covmatrix$param["nugget_2"]
    }
    corri <- rho
  } else {
    nug <- as.numeric(covmatrix$param["nugget"])
    sill <- as.numeric(covmatrix$param["sill"])
    if (is.na(nug)) nug <- 0
    if (is.na(sill)) stop("sill parameter is missing")
    corri <- (1 - nug) * rho
    vvar <- sill
  }

  CC <- matrix(corri * as.numeric(vvar), nrow = dimat, ncol = dimat2)
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
    bivariate = bivariate,
    spacetime = covmatrix$spacetime,
    param = covmatrix$param,
    xloc = Xloc,
    tloc = tloc,
    loc = loc,
    loc_metric = loc_metric,
    call = call
  )

  if (return_covmatrix) result$covmatrix <- covmatrix$covmatrix
  if (return_CC) result$CC <- CC

  return(result)
}