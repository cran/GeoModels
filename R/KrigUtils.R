####################################################
### Internal utilities for global and local kriging
####################################################

.GeoKrigMatchMethod <- function(method = "cholesky", sparse = FALSE,
                                context = "Kriging") {
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop(context, ": method must be 'cholesky' or 'svd'.", call. = FALSE)
  }
  method <- tolower(gsub("[[:blank:]]", "", method))
  if (!(method %in% c("cholesky", "svd"))) {
    stop(context, ": method must be 'cholesky' or 'svd'.", call. = FALSE)
  }
  if (isTRUE(sparse) && identical(method, "svd")) {
    stop(context, ": method = 'svd' is not available with sparse = TRUE; ",
         "use method = 'cholesky'.", call. = FALSE)
  }
  method
}

.GeoKrigCoordinateMatrix <- function(coordx = NULL, coordy = NULL,
                                     coordz = NULL, grid = FALSE,
                                     context = "Kriging") {
  if (is.null(coordx)) {
    stop(context, ": observed coordinates are missing.", call. = FALSE)
  }

  if ((is.matrix(coordx) || is.data.frame(coordx)) &&
      is.null(coordy) && is.null(coordz)) {
    coords <- as.matrix(coordx)
  } else {
    if (is.null(coordy)) {
      stop(context, ": coordy is required when coordx is not a coordinate matrix.",
           call. = FALSE)
    }
    if (!is.numeric(coordx) || !is.numeric(coordy) ||
        any(!is.finite(coordx)) || any(!is.finite(coordy))) {
      stop(context, ": coordx and coordy must be finite numeric vectors.",
           call. = FALSE)
    }
    if (!is.null(coordz) &&
        (!is.numeric(coordz) || any(!is.finite(coordz)))) {
      stop(context, ": coordz must be a finite numeric vector.", call. = FALSE)
    }

    if (isTRUE(grid)) {
      coords <- if (is.null(coordz)) {
        as.matrix(expand.grid(coordx, coordy))
      } else {
        as.matrix(expand.grid(coordx, coordy, coordz))
      }
    } else {
      if (length(coordx) != length(coordy)) {
        stop(context, ": coordx and coordy must have the same length.",
             call. = FALSE)
      }
      if (!is.null(coordz) && length(coordz) != length(coordx)) {
        stop(context, ": coordz must have the same length as coordx and coordy.",
             call. = FALSE)
      }
      coords <- if (is.null(coordz)) {
        cbind(coordx, coordy)
      } else {
        cbind(coordx, coordy, coordz)
      }
    }
  }

  if (!is.matrix(coords) || !is.numeric(coords) ||
      nrow(coords) < 1L || !(ncol(coords) %in% c(2L, 3L)) ||
      any(!is.finite(coords))) {
    stop(context, ": coordinates must be a finite numeric N x 2 or N x 3 matrix.",
         call. = FALSE)
  }
  storage.mode(coords) <- "double"
  unname(coords)
}

.GeoKrigNormalizeData <- function(data, nspace, coordt = NULL,
                                  grid = FALSE,
                                  context = "Kriging") {
  if (is.null(data)) stop(context, ": data are missing.", call. = FALSE)
  if (!is.numeric(data) || any(!is.finite(data))) {
    stop(context, ": data must be finite and numeric.", call. = FALSE)
  }
  if (!is.numeric(nspace) || length(nspace) != 1L ||
      !is.finite(nspace) || nspace < 1L || nspace != as.integer(nspace)) {
    stop(context, ": invalid number of spatial observations.", call. = FALSE)
  }
  nspace <- as.integer(nspace)

  if (is.null(coordt)) {
    if (length(data) != nspace) {
      stop(context, ": data must contain one value per spatial coordinate (",
           nspace, " values expected).", call. = FALSE)
    }
    ## For regular grids this is the order of expand.grid(coordx, coordy):
    ## the first coordinate varies fastest.  It is also the order returned by
    ## GeoSim for a spatial grid.
    return(unname(as.numeric(data)))
  }

  if (!is.numeric(coordt) || !length(coordt) || any(!is.finite(coordt))) {
    stop(context, ": coordt must be a non-empty finite numeric vector.",
         call. = FALSE)
  }
  nt <- length(coordt)
  total <- nspace * nt

  if (is.matrix(data)) {
    dd <- dim(data)
    if (identical(dd, c(nt, nspace))) {
      ## Rows are times, columns are spatial sites.
      return(unname(as.numeric(t(data))))
    }
    if (identical(dd, c(nspace, nt))) {
      ## Columns are times; column-major vectorisation is already time-major.
      return(unname(as.numeric(data)))
    }
  }

  if (length(dim(data)) >= 3L) {
    dd <- dim(data)
    if (!isTRUE(grid) || utils::tail(dd, 1L) != nt || prod(utils::head(dd, -1L)) != nspace) {
      stop(context, ": a regular-grid space-time array must have the last ",
           "dimension equal to length(coordt) and the remaining dimensions ",
           "equal to the spatial grid size.", call. = FALSE)
    }
    return(unname(as.numeric(data)))
  }

  if (is.vector(data) && length(data) == total) {
    return(unname(as.numeric(data)))
  }

  stop(context, ": data have an incompatible format for fixed-location ",
       "space-time observations; expected a length-", total,
       " vector, a ", nt, " x ", nspace, " matrix, a ", nspace,
       " x ", nt, " matrix, or a compatible regular-grid array.",
       call. = FALSE)
}

.GeoKrigObservedCoordinates <- function(covmatrix, dimat,
                                        context = "Kriging") {
  if (!is.null(covmatrix$coordx_dyn)) {
    coords <- do.call(rbind, lapply(covmatrix$coordx_dyn, as.matrix))
  } else {
    base <- if (is.null(covmatrix$coordz)) {
      cbind(covmatrix$coordx, covmatrix$coordy)
    } else {
      cbind(covmatrix$coordx, covmatrix$coordy, covmatrix$coordz)
    }
    repeats <- if (isTRUE(covmatrix$spacetime) || isTRUE(covmatrix$bivariate)) {
      as.integer(covmatrix$numtime)
    } else {
      1L
    }
    coords <- base[rep(seq_len(nrow(base)), times = repeats), , drop = FALSE]
  }

  coords <- as.matrix(coords)
  if (ncol(coords) == 2L) coords <- cbind(coords, 0)
  if (ncol(coords) != 3L || nrow(coords) != dimat ||
      !is.numeric(coords) || any(!is.finite(coords))) {
    stop(context, ": internal observed-coordinate ordering is inconsistent ",
         "with the covariance matrix.", call. = FALSE)
  }
  storage.mode(coords) <- "double"
  unname(coords)
}

.GeoKrigFactor <- function(covmatrix, method = "cholesky",
                           context = "Kriging") {
  method <- .GeoKrigMatchMethod(method, sparse = isTRUE(covmatrix$sparse),
                                context = context)
  C <- covmatrix$covmatrix
  nobs <- nrow(C)

  if (isTRUE(covmatrix$sparse)) {
    Cspam <- if (spam::is.spam(C)) C else spam::as.spam(C)
    U <- tryCatch(spam::chol.spam(Cspam), error = function(e) NULL)
    if (is.null(U)) {
      stop(context, ": covariance matrix is not positive definite.",
           call. = FALSE)
    }
    return(structure(
      list(method = "cholesky", sparse = TRUE, factor = U, nobs = nobs),
      class = "GeoKrigFactor"
    ))
  }

  C <- as.matrix(C)
  if (identical(method, "cholesky")) {
    U <- tryCatch(FastGP::rcppeigen_get_chol(C), error = function(e) NULL)
    if (is.null(U)) {
      stop(context, ": covariance matrix is not positive definite.",
           call. = FALSE)
    }
    return(structure(
      list(method = "cholesky", sparse = FALSE, factor = U, nobs = nobs),
      class = "GeoKrigFactor"
    ))
  }

  dec <- tryCatch(svd(C), error = function(e) NULL)
  if (is.null(dec) || !length(dec$d) || any(!is.finite(dec$d))) {
    stop(context, ": SVD decomposition of the covariance matrix failed.",
         call. = FALSE)
  }
  tol <- max(dim(C)) * max(dec$d) * .Machine$double.eps
  if (any(dec$d <= tol)) {
    stop(context, ": covariance matrix is numerically singular under SVD ",
         "(smallest singular value = ", format(min(dec$d)), ").",
         call. = FALSE)
  }
  structure(
    list(method = "svd", sparse = FALSE, factor = dec, nobs = nobs),
    class = "GeoKrigFactor"
  )
}

.GeoKrigSolveFactor <- function(factor, rhs, context = "Kriging") {
  if (!inherits(factor, "GeoKrigFactor")) {
    stop(context, ": invalid cached kriging factorization.", call. = FALSE)
  }
  rhs <- as.matrix(rhs)
  if (!is.numeric(rhs) || any(!is.finite(rhs))) {
    stop(context, ": the right-hand side of the kriging system must be finite.",
         call. = FALSE)
  }
  if (nrow(rhs) != factor$nobs) {
    stop(context, ": incompatible kriging-system dimensions.", call. = FALSE)
  }

  if (isTRUE(factor$sparse)) {
    tmp <- spam::forwardsolve(factor$factor, rhs)
    return(as.matrix(spam::backsolve(factor$factor, tmp)))
  }

  if (identical(factor$method, "cholesky")) {
    tmp <- forwardsolve(factor$factor, rhs)
    return(as.matrix(forwardsolve(factor$factor, tmp, transpose = TRUE)))
  }

  dec <- factor$factor
  as.matrix(dec$v %*% ((crossprod(dec$u, rhs)) / dec$d))
}

.GeoKrigSolve <- function(covmatrix, rhs, method = "cholesky",
                          context = "Kriging") {
  factor <- .GeoKrigFactor(covmatrix, method = method, context = context)
  .GeoKrigSolveFactor(factor, rhs, context = context)
}

.GeoKrigDenseSPDSolve <- function(C, rhs, context = "Kriging") {
  C <- as.matrix(C)
  rhs <- as.matrix(rhs)
  n <- nrow(C)
  nrhs <- ncol(rhs)

  if (!is.numeric(C) || n < 1L || ncol(C) != n || any(!is.finite(C))) {
    stop(context, ": covariance matrix must be a finite square numeric matrix.",
         call. = FALSE)
  }
  if (!is.numeric(rhs) || nrhs < 1L || nrow(rhs) != n ||
      any(!is.finite(rhs))) {
    stop(context, ": incompatible or non-finite kriging right-hand side.",
         call. = FALSE)
  }

  ## dposv performs Cholesky factorization and all triangular solves in one
  ## native LAPACK call.  Private copies are passed because LAPACK overwrites
  ## both the covariance matrix and the right-hand sides.  This is especially
  ## useful for GeoKrigloc, where thousands of small dense systems would
  ## otherwise cross the R/native boundary several times per neighborhood.
  z <- dotCall64::.C64(
    "GeoDenseSPDSolve",
    SIGNATURE = c("double", "double", "integer", "integer", "integer"),
    A = as.double(C),
    B = as.double(rhs),
    as.integer(n), as.integer(nrhs), info = as.integer(0L),
    INTENT = c("rw", "rw", "r", "r", "w"),
    PACKAGE = "GeoModels", VERBOSE = 0, NAOK = TRUE
  )
  info <- as.integer(z$info)[1L]
  if (!is.finite(info) || info != 0L) {
    if (is.finite(info) && info > 0L) {
      stop(context, ": covariance matrix is not positive definite.",
           call. = FALSE)
    }
    stop(context, ": native dense Cholesky solve failed (LAPACK info = ",
         info, ").", call. = FALSE)
  }
  out <- matrix(z$B, nrow = n, ncol = nrhs)
  if (any(!is.finite(out))) {
    stop(context, ": native dense Cholesky solve returned non-finite values.",
         call. = FALSE)
  }
  out
}

.GeoKrigLocationBlocks <- function(nobs, nloc, tloc = 1L,
                                   context = "Kriging") {
  nobs <- as.integer(nobs)
  nloc <- as.integer(nloc)
  tloc <- as.integer(tloc)
  if (any(!is.finite(c(nobs, nloc, tloc))) || nobs < 1L || nloc < 1L ||
      tloc < 1L) {
    stop(context, ": invalid dimensions for blocked prediction.", call. = FALSE)
  }

  ## The option is intentionally internal.  It controls the approximate size
  ## of one dense cross-covariance matrix, not the total memory footprint.
  max_bytes <- getOption("GeoModels.krig.block_bytes", 32 * 1024^2)
  max_bytes <- suppressWarnings(as.numeric(max_bytes)[1L])
  if (!is.finite(max_bytes) || max_bytes <= 0) max_bytes <- 32 * 1024^2

  bytes_per_location <- 8 * as.double(nobs) * as.double(tloc)
  block_nloc <- max(1L, min(nloc, floor(max_bytes / bytes_per_location)))
  starts <- seq.int(1L, nloc, by = block_nloc)
  lapply(starts, function(i) seq.int(i, min(nloc, i + block_nloc - 1L)))
}

.GeoKrigCorrelationBlock <- function(covmatrix, observed_coords,
                                     corrmodel, corrparam, distance,
                                     loc_block, time, tloc, NS,
                                     which = 1L, context = "Kriging") {
  loc_block <- as.matrix(loc_block)
  if (!is.numeric(loc_block) || nrow(loc_block) < 1L ||
      !(ncol(loc_block) %in% c(2L, 3L)) || any(!is.finite(loc_block))) {
    stop(context, ": invalid prediction-coordinate block.", call. = FALSE)
  }
  if (ncol(loc_block) == 2L) loc_block <- cbind(loc_block, 0)
  observed_coords <- as.matrix(observed_coords)
  if (ncol(observed_coords) == 2L) observed_coords <- cbind(observed_coords, 0)

  nobs <- nrow(observed_coords)
  nloc_block <- nrow(loc_block)
  ntarget <- nloc_block * as.integer(tloc)
  out <- dotCall64::.C64(
    "Corr_c",
    SIGNATURE = c(
      "double", "double", "double", "double", "double", "integer",
      "integer", "double", "double", "double", "integer", "integer",
      "integer", "integer", "integer", "integer", "double", "integer",
      "integer", "double", "integer", "integer", "double"
    ),
    corri = dotCall64::vector_dc("double", nobs * ntarget),
    observed_coords[, 1L], observed_coords[, 2L], observed_coords[, 3L],
    covmatrix$coordt, corrmodel, 0,
    loc_block[, 1L], loc_block[, 2L], loc_block[, 3L],
    covmatrix$numcoord, nloc_block, as.integer(tloc),
    covmatrix$ns, NS, covmatrix$numtime, corrparam,
    covmatrix$spacetime, covmatrix$bivariate, time, distance,
    as.integer(which) - 1L, covmatrix$radius,
    INTENT = c("w", rep("r", 22)), NAOK = TRUE, PACKAGE = "GeoModels"
  )$corri
  if (length(out) != nobs * ntarget || any(!is.finite(out))) {
    stop(context, ": non-finite cross-correlation block.", call. = FALSE)
  }
  out
}

.GeoKrigMeanVarcov <- function(varcov, mean_names,
                                context = "Kriging") {
  if (is.null(varcov)) {
    stop(context, ": type_krig = 'Universal' with mse = TRUE requires ",
         "varcov.  For a GeoFit object, compute/update its parameter ",
         "covariance matrix before kriging.", call. = FALSE)
  }

  V <- as.matrix(varcov)
  if (!is.numeric(V) || nrow(V) != ncol(V) || any(!is.finite(V))) {
    stop(context, ": varcov must be a finite square numeric matrix.",
         call. = FALSE)
  }
  rn <- rownames(V)
  cn <- colnames(V)
  if (is.null(rn) || is.null(cn) || anyDuplicated(rn) || anyDuplicated(cn)) {
    stop(context, ": varcov must have unique row and column names.",
         call. = FALSE)
  }
  if (!setequal(rn, cn)) {
    stop(context, ": row and column names of varcov must identify the same ",
         "parameters.", call. = FALSE)
  }

  out <- matrix(0, nrow = length(mean_names), ncol = length(mean_names),
                dimnames = list(mean_names, mean_names))
  available <- mean_names[mean_names %in% rn & mean_names %in% cn]
  if (length(available)) {
    block <- V[available, available, drop = FALSE]
    asym <- max(abs(block - t(block)))
    tol <- 1e-8 * max(1, max(abs(block)))
    if (!is.finite(asym) || asym > tol) {
      stop(context, ": the mean-parameter block of varcov is not symmetric.",
           call. = FALSE)
    }
    out[available, available] <- (block + t(block)) / 2
  }
  out
}


.GeoKrigSystem <- function(covmatrix, CC, method = "cholesky",
                           type_krig = "Simple", X = NULL, Xloc = NULL,
                           target_variance = NULL, mse = FALSE,
                           varcov = NULL, mean_names = NULL, factor = NULL,
                           weights = NULL, context = "Kriging") {
  type_key <- tolower(gsub("[[:blank:]]", "", type_krig))
  if (!(type_key %in% c("simple", "universal"))) {
    stop(context, ": type_krig must be 'Simple' or 'Universal'.",
         call. = FALSE)
  }

  CC <- as.matrix(CC)
  nobs <- nrow(CC)
  npred <- ncol(CC)
  if (nobs != nrow(covmatrix$covmatrix)) {
    stop(context, ": cross-covariance matrix has the wrong number of rows.",
         call. = FALSE)
  }

  ## Both Simple and Universal use the mean coefficients supplied by the
  ## fitted model (ML, REML, composite likelihood, or fixed by the user).
  ## Universal differs only in the prediction-MSE correction below.
  if (is.null(weights)) {
    weights <- if (is.null(factor)) {
      .GeoKrigSolve(covmatrix, CC, method = method, context = context)
    } else {
      .GeoKrigSolveFactor(factor, CC, context = context)
    }
  } else {
    weights <- as.matrix(weights)
    if (!is.numeric(weights) || any(!is.finite(weights)) ||
        nrow(weights) != nobs || ncol(weights) != npred) {
      stop(context, ": invalid precomputed kriging weights.", call. = FALSE)
    }
  }
  quad <- colSums(CC * weights)

  mse_value <- NULL
  if (isTRUE(mse)) {
    if (is.null(target_variance)) {
      stop(context, ": target variance is required to compute the MSE.",
           call. = FALSE)
    }
    mse_value <- rep(as.numeric(target_variance), length.out = npred) - quad

    if (identical(type_key, "universal")) {
      X <- as.matrix(X)
      Xloc <- as.matrix(Xloc)
      if (!is.numeric(X) || any(!is.finite(X)) || nrow(X) != nobs ||
          ncol(X) < 1L) {
        stop(context, ": X must be a finite design matrix with one row per ",
             "observation for universal-MSE correction.", call. = FALSE)
      }
      if (!is.numeric(Xloc) || any(!is.finite(Xloc)) ||
          nrow(Xloc) != npred || ncol(Xloc) != ncol(X)) {
        stop(context, ": Xloc must have one row per prediction task and the ",
             "same number of columns as X for universal-MSE correction.",
             call. = FALSE)
      }
      if (is.null(mean_names)) {
        mean_names <- colnames(X)
        if (is.null(mean_names) || any(!nzchar(mean_names))) {
          mean_names <- .GeoMean_names(ncol(X))
        }
      }
      if (length(mean_names) != ncol(X) || anyDuplicated(mean_names)) {
        stop(context, ": invalid mean-parameter names for universal kriging.",
             call. = FALSE)
      }

      var_mean <- .GeoKrigMeanVarcov(
        varcov, mean_names = mean_names, context = context
      )
      trend_gap <- Xloc - t(weights) %*% X
      correction <- rowSums((trend_gap %*% var_mean) * trend_gap)
      mse_value <- mse_value + correction
    }
  }

  list(weights = weights, mse = mse_value,
       type = if (identical(type_key, "universal")) "Universal" else "Simple")
}

################################################################################
## Wrapped-Gaussian circular component helpers                                 ##
################################################################################

.GeoWrappedMeanDirection <- function(eta) {
  (2 * atan(as.numeric(eta)) + pi) %% (2 * pi)
}

.GeoWrappedSinVariance <- function(sill) {
  sill <- as.numeric(sill)[1L]
  if(!is.finite(sill) || sill <= 0)
    stop("Wrapped Gaussian prediction requires sill > 0.", call. = FALSE)
  0.5 * (-expm1(-2 * sill))
}

.GeoWrappedCosVariance <- function(sill) {
  sill <- as.numeric(sill)[1L]
  if(!is.finite(sill) || sill <= 0)
    stop("Wrapped Gaussian prediction requires sill > 0.", call. = FALSE)
  0.5 * (1 - exp(-sill))^2
}

.GeoWrappedSinCovariance <- function(rho, sill) {
  sill <- as.numeric(sill)[1L]
  if(!is.finite(sill) || sill <= 0)
    stop("Wrapped Gaussian prediction requires sill > 0.", call. = FALSE)
  rho <- pmax(-1, pmin(1, as.numeric(rho)))
  0.5 * (exp(-sill * (1 - rho)) - exp(-sill * (1 + rho)))
}

.GeoWrappedCosCovariance <- function(rho, sill) {
  sill <- as.numeric(sill)[1L]
  if(!is.finite(sill) || sill <= 0)
    stop("Wrapped Gaussian prediction requires sill > 0.", call. = FALSE)
  rho <- pmax(-1, pmin(1, as.numeric(rho)))
  0.5 * (exp(-sill * (1 - rho)) + exp(-sill * (1 + rho))) - exp(-sill)
}

.GeoWrappedCircularCorrelation <- function(rho, sill) {
  ## Stable evaluation of sinh(sill * rho) / sinh(sill).
  .GeoWrappedSinCovariance(rho, sill) / .GeoWrappedSinVariance(sill)
}

################################################################################
## Prediction model/counter helpers                                            ##
################################################################################
.GeoPredictionCanonicalModel <- function(model) {
  if (!is.character(model) || length(model) != 1L || is.na(model)) return(model)
  model <- gsub("[[:blank:]]", "", model)
  aliases <- c(
    Gauss = "Gaussian", LogGauss = "LogGaussian", SkewGauss = "SkewGaussian",
    TwoPieceGauss = "TwoPieceGaussian", tukeyh = "Tukeyh", tukeyh2 = "Tukeyh2",
    poisson = "Poisson", poissongamma = "PoissonGamma",
    Binomiallogistic = "BinomialLogistic",
    Gaussian_misp_StudentT = "StudentT",
    Gaussian_misp_Poisson = "Poisson",
    Gaussian_misp_SkewStudentT = "SkewStudentT",
    Gaussian_misp_Tukeygh = "Tukeygh",
    Gaussian_misp_PoissonZIP = "PoissonZIP",
    Gaussian_misp_PoissonGamma = "PoissonGamma",
    Gaussian_misp_Binomial = "Binomial",
    Gaussian_misp_BinomialNeg = "BinomialNeg"
  )
  if (model %in% names(aliases)) unname(aliases[[model]]) else model
}

.GeoKrigSupportedPlainModels <- c(
  "Gaussian", "SkewGaussian", "SkewLaplace", "StudentT", "SkewStudentT",
  "Gamma", "Weibull", "LogLogistic", "LogGaussian", "TwoPieceStudentT",
  "Beta", "TwoPieceGaussian", "TwoPieceTukeyh", "TwoPieceBimodal",
  "Tukeygh", "Tukeyh", "Tukeyh2", "SinhAsinh", "Wrapped",
  "Binary", "Bernoulli", "Binomial", "Geom", "Geometric", "BinomialNeg",
  "Poisson", "PoissonZIP", "BinomialNegZINB",
  "PoissonGamma", "PoissonGammaZIP", "PoissonGammaZIP1"
)

.GeoKrigValidatePredictionModel <- function(model, copula = NULL, bivariate = FALSE,
                                            context = "GeoKrig") {
  if (isTRUE(bivariate) && !identical(model, "Gaussian")) {
    stop(context, ": bivariate prediction is currently implemented only for model = 'Gaussian'.",
         call. = FALSE)
  }
  if (is.null(copula)) {
    if (!(model %in% .GeoKrigSupportedPlainModels)) {
      stop(context, ": prediction is not implemented for model = '", model,
           "' without a copula.", call. = FALSE)
    }
  } else {
    model_code <- as.integer(CkModel(model))
    supported_codes <- .GeoGaussianCopulaContinuousModels
    if (identical(copula, "Gaussian") || identical(context, "GeoKrig"))
      supported_codes <- c(supported_codes, .GeoGaussianCopulaDiscreteModels)
    if (!(model_code %in% supported_codes)) {
      if (model_code %in% .GeoGaussianCopulaDiscreteModels &&
          copula %in% c("Clayton", "SkewGaussian") && !identical(context, "GeoKrig")) {
        stop(context, ": ", copula,
             " copula with Poisson/Binomial/BinomialNeg is currently available in global GeoKrig/GeoCV only; local GeoKrigloc is not yet implemented.",
             call. = FALSE)
      }
      stop(context, ": copula-based prediction is not implemented for copula = '",
           copula, "' and model = '", model, "'.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

.GeoKrigPositiveIntegerVector <- function(x, name, n = NULL) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
      any(x < 1) || any(abs(x - round(x)) > sqrt(.Machine$double.eps))) {
    stop(name, " must contain positive integers.", call. = FALSE)
  }
  x <- as.integer(round(x))
  if (!is.null(n)) {
    if (length(x) == 1L) x <- rep.int(x, n)
    if (length(x) != n) {
      stop(name, " must be scalar or contain one value per corresponding location.",
           call. = FALSE)
    }
  }
  x
}

.GeoKrigCountSizes <- function(model, n, nloc, nobs, npred,
                               context = "GeoKrig") {
  model <- .GeoPredictionCanonicalModel(model)
  ## Models whose response mean/variance depends on a size/trial parameter.
  binomial_models <- c("Binomial", "Binomial2")
  negbin_models <- c("BinomialNeg", "BinomialNegZINB")
  if (model %in% c("Binary", "Bernoulli")) {
    if(!is.numeric(n) || length(n) != 1L || !is.finite(n) || as.integer(round(n)) != 1L)
      stop(context, ": Binary/Bernoulli is the one-trial field; n must equal 1.", call.=FALSE)
    if(!is.null(nloc) && any(as.numeric(nloc) != 1))
      stop(context, ": Binary/Bernoulli prediction requires nloc=1 when nloc is supplied.", call.=FALSE)
    return(list(obs = rep.int(1L, nobs), pred = rep.int(1L, npred)))
  }
  if (model %in% c("Geom", "Geometric")) {
    if(!is.numeric(n) || length(n) != 1L || !is.finite(n) || as.integer(round(n)) != 1L)
      stop(context, ": Geometric is the r = 1 Negative-Binomial field; n must equal 1.", call.=FALSE)
    if(!is.null(nloc) && any(as.numeric(nloc) != 1))
      stop(context, ": Geometric is the r = 1 Negative-Binomial field; nloc, when supplied, must equal 1.", call.=FALSE)
    return(list(obs = rep.int(1L, nobs), pred = rep.int(1L, npred)))
  }
  if (!(model %in% c(binomial_models, negbin_models))) {
    ## Corr_c_bin still expects integer vectors; ones are neutral for Poisson-type models.
    return(list(obs = rep.int(1L, nobs), pred = rep.int(1L, npred)))
  }

  if(model %in% negbin_models) {
    r <- .GeoKrigPositiveIntegerVector(n, paste0(context, ": n"), n = NULL)
    if(length(r) != 1L)
      stop(context, ": for BinomialNeg, n is the common number r of successes and must be a positive integer scalar.", call.=FALSE)
    if(!is.null(nloc)) {
      rloc <- .GeoKrigPositiveIntegerVector(nloc, paste0(context, ": nloc"), n = NULL)
      if(any(rloc != r[1L]))
        stop(context, ": for BinomialNeg, r is a common field parameter; nloc cannot differ from n.", call.=FALSE)
    }
    return(list(obs = rep.int(r[1L], nobs), pred = rep.int(r[1L], npred)))
  }

  obs <- .GeoKrigPositiveIntegerVector(n, paste0(context, ": n"), n = nobs)
  if (is.null(nloc)) {
    if (length(n) != 1L) {
      stop(context, ": nloc is required when Binomial n is location-specific.", call. = FALSE)
    }
    pred <- rep.int(as.integer(round(n[1L])), npred)
  } else {
    pred <- .GeoKrigPositiveIntegerVector(nloc, paste0(context, ": nloc"), n = npred)
  }
  list(obs = obs, pred = pred)
}

.GeoKrigCountMeanDerivative <- function(model_code, param, eta, size) {
  eta <- as.numeric(eta)
  size <- rep(as.numeric(size), length.out = length(eta))
  code <- as.integer(model_code)
  if (code == 30L) return(exp(eta))                         # Poisson
  if (code %in% c(46L, 57L, 58L)) {                        # Poisson-Gamma (+ ZIP)
    d <- exp(eta)
    if (code %in% c(57L, 58L)) d <- (1 - pnorm(as.numeric(param["pmu"]))) * d
    return(d)
  }
  if (code == 43L) {                                       # Poisson ZIP
    return((1 - pnorm(as.numeric(param["pmu"]))) * exp(eta))
  }
  if (code %in% c(2L, 11L, 19L)) {                         # Binary/Binomial
    return(size * dnorm(eta))
  }
  if (code %in% c(14L, 16L)) {                             # Geometric/NegBin
    p <- pnorm(eta)
    return(-size * dnorm(eta) / p^2)
  }
  if (code == 45L) {                                       # zero-inflated NegBin
    p <- pnorm(eta)
    z <- pnorm(as.numeric(param["pmu"]))
    return(-(1 - z) * size * dnorm(eta) / p^2)
  }
  stop("GeoKrig: no marginal-mean derivative is implemented for this count model.",
       call. = FALSE)
}

###########################################################################
## Sparse observation-pair correlation cache for local kriging           ##
###########################################################################

.GeoLocalPairOccurrenceCount <- function(indices) {
  if (!is.list(indices) || !length(indices)) return(0)
  sum(vapply(indices, function(idx) {
    k <- length(idx)
    if (k < 2L) 0 else as.double(k) * as.double(k - 1L) / 2
  }, numeric(1L)))
}

.GeoPairCacheOption <- function(auto = FALSE) {
  opt <- getOption("GeoModels.local_pair_cache", NULL)
  if (is.null(opt) || identical(opt, "auto")) return(isTRUE(auto))
  if (!is.logical(opt) || length(opt) != 1L || is.na(opt)) {
    stop("Option 'GeoModels.local_pair_cache' must be TRUE, FALSE, 'auto', or NULL.",
         call. = FALSE)
  }
  isTRUE(opt)
}

.GeoPairCacheNear <- function(x, values) {
  x <- suppressWarnings(as.numeric(x)[1L])
  is.finite(x) && any(abs(x - values) <= 100 * .Machine$double.eps * max(1, abs(x)))
}

.GeoPairCacheSpatialAuto <- function(corrmodel, param, weights = FALSE) {
  if (isTRUE(weights)) return(TRUE)
  code <- as.integer(CkCorrModel(corrmodel))
  if (!length(code) || is.na(code)) return(FALSE)

  ## Always expensive: Kummer / hypergeometric / hole-effect families.
  if (code %in% c(21L, 22L, 23L, 24L, 25L, 26L, 27L, 29L, 30L)) return(TRUE)

  ## General Matern is expensive except the closed half-integer cases.
  if (code == 14L) {
    return(!.GeoPairCacheNear(param[["smooth"]], c(0.5, 1.5, 2.5, 3.5)))
  }

  ## GenWend and its Matern reparameterizations share CorFunW_gen().
  ## smooth 0:3 use the cheap closed forms measured in the local benchmark.
  if (code %in% c(6L, 7L, 19L)) {
    return(!.GeoPairCacheNear(param[["smooth"]], 0:3))
  }
  FALSE
}

.GeoPairCacheSpacetimeAuto <- function(corrmodel, param, weights = FALSE) {
  if (isTRUE(weights)) return(TRUE)
  code <- as.integer(CkCorrModel(corrmodel))
  if (!length(code) || is.na(code)) return(FALSE)

  ## Space-time families containing a general Matern or generalized Wendland
  ## evaluation.  Their pair correlation depends only on fixed (h_ij, u_ij)
  ## during one local-kriging call, so it is safe to reuse across targets.
  code %in% c(61L, 62L, 78L, 85L, 86L, 87L, 88L, 89L, 95L)
}

.GeoPairCacheEnabled <- function(corrmodel, param, weights = FALSE,
                                 spacetime = FALSE) {
  auto <- if (isTRUE(spacetime)) {
    .GeoPairCacheSpacetimeAuto(corrmodel, param, weights = weights)
  } else {
    .GeoPairCacheSpatialAuto(corrmodel, param, weights = weights)
  }
  .GeoPairCacheOption(auto)
}

.GeoPairCacheThresholdOK <- function(indices) {
  total_pairs <- .GeoLocalPairOccurrenceCount(indices)
  min_pairs <- suppressWarnings(as.numeric(
    getOption("GeoModels.local_pair_cache_min_pairs", 1e6)
  )[1L])
  if (!is.finite(min_pairs) || min_pairs < 0) min_pairs <- 1e6
  is.finite(total_pairs) && total_pairs >= min_pairs
}

.GeoValidatePairCacheGeometry <- function(coords, times = NULL,
                                          context = "local kriging") {
  coords <- as.matrix(coords)
  if (!is.numeric(coords) || nrow(coords) < 1L ||
      !(ncol(coords) %in% c(2L, 3L)) || any(!is.finite(coords))) {
    stop(context, ": invalid coordinates for pair-cache construction.",
         call. = FALSE)
  }
  storage.mode(coords) <- "double"
  if (!is.null(times)) {
    times <- as.numeric(times)
    if (length(times) != nrow(coords) || any(!is.finite(times))) {
      stop(context, ": invalid times for pair-cache construction.", call. = FALSE)
    }
  }
  list(coords = coords, times = times)
}

.GeoPairCacheBudgetOK <- function(indices, nobs) {
  budget <- getOption("GeoModels.local_pair_cache_max_bytes", 256 * 1024^2)
  if (!is.numeric(budget) || length(budget) != 1L ||
      !is.finite(budget) || budget < 0)
    stop("GeoModels.local_pair_cache_max_bytes must be a finite nonnegative number.",
         call. = FALSE)
  if (!is.list(indices)) stop("Neighborhood indices must be a list.", call. = FALSE)
  ## These arrays alone must fit, even before accounting for the hash table.
  k <- if (length(indices)) max(lengths(indices)) else 0
  lower_bound <- 65536 + 4 * (as.double(nobs) + 1) + 12 * k * (k - 1) / 2
  budget >= lower_bound
}

.GeoBuildSpatialPairCache <- function(indices, coords, corrmodel_code,
                                      corrparam, distance_code, radius,
                                      context = "local spatial kriging") {
  if (!.GeoPairCacheThresholdOK(indices)) return(NULL)
  geom <- .GeoValidatePairCacheGeometry(coords, context = context)
  if (!.GeoPairCacheBudgetOK(indices, nrow(geom$coords))) return(NULL)
  corrparam <- as.double(corrparam)
  if (!length(corrparam) || any(!is.finite(corrparam))) {
    stop(context, ": invalid correlation parameters for pair cache.",
         call. = FALSE)
  }

  cache <- .Call(
    "GeoBuildSpatialPairCache",
    indices, geom$coords, as.integer(corrmodel_code), corrparam,
    as.integer(distance_code), as.double(radius),
    PACKAGE = "GeoModels"
  )
  .GeoValidatePairCache(cache, context)
}

.GeoBuildSpacetimePairCache <- function(indices, coords, times,
                                        corrmodel_code, corrparam,
                                        distance_code, radius,
                                        context = "local space-time kriging") {
  if (!.GeoPairCacheThresholdOK(indices)) return(NULL)
  geom <- .GeoValidatePairCacheGeometry(coords, times, context)
  if (!.GeoPairCacheBudgetOK(indices, nrow(geom$coords))) return(NULL)
  corrparam <- as.double(corrparam)
  if (!length(corrparam) || any(!is.finite(corrparam))) {
    stop(context, ": invalid correlation parameters for pair cache.",
         call. = FALSE)
  }

  cache <- .Call(
    "GeoBuildSpacetimePairCache",
    indices, geom$coords, geom$times, as.integer(corrmodel_code), corrparam,
    as.integer(distance_code), as.double(radius),
    PACKAGE = "GeoModels"
  )
  .GeoValidatePairCache(cache, context)
}

.GeoValidatePairCache <- function(cache, context) {
  ## NULL is an intentional memory-budget/allocation fallback, not corrupt data.
  if (is.null(cache)) return(NULL)
  if (!is.list(cache) || is.null(cache$row_ptr) || is.null(cache$col_idx) ||
      is.null(cache$rho) || length(cache$col_idx) != length(cache$rho)) {
    stop(context, ": native pair-cache construction returned invalid data.",
         call. = FALSE)
  }
  cache
}

.GeoPairCachePacked <- function(indices, cache,
                                context = "local kriging") {
  if (is.null(cache)) return(NULL)
  idx <- as.integer(indices)
  if (length(idx) < 2L) return(numeric(0L))
  out <- .Call(
    "GeoPairCachePacked",
    idx, cache$row_ptr, cache$col_idx, cache$rho,
    PACKAGE = "GeoModels"
  )
  if (length(out) != as.double(length(idx)) * as.double(length(idx) - 1L) / 2 ||
      any(!is.finite(out))) {
    stop(context, ": invalid pair-cache lookup result.", call. = FALSE)
  }
  as.numeric(out)
}

.GeoSpacetimeObservationGeometry <- function(coords, coordx_dyn, coordt,
                                             context = "local space-time kriging") {
  coordt <- as.numeric(coordt)
  if (!length(coordt) || any(!is.finite(coordt))) {
    stop(context, ": invalid observation times.", call. = FALSE)
  }
  if (!is.null(coordx_dyn)) {
    if (!is.list(coordx_dyn) || length(coordx_dyn) != length(coordt)) {
      stop(context, ": dynamic coordinates and times are inconsistent.", call. = FALSE)
    }
    sizes <- vapply(coordx_dyn, nrow, integer(1L))
    flat <- do.call(rbind, lapply(coordx_dyn, as.matrix))
    times <- rep(coordt, sizes)
  } else {
    coords <- as.matrix(coords)
    ns <- nrow(coords)
    flat <- coords[rep(seq_len(ns), times = length(coordt)), , drop = FALSE]
    times <- rep(coordt, each = ns)
  }
  geom <- .GeoValidatePairCacheGeometry(flat, times, context)
  list(coords = geom$coords, times = geom$times)
}

.GeoPairCrossCorrelation <- function(indices, coords, times = NULL, loc,
                                     pred_time = 0, corrmodel_code, corrparam,
                                     distance_code, radius,
                                     context = "local kriging") {
  geom <- .GeoValidatePairCacheGeometry(coords, times, context)
  .GeoPairCrossCorrelationValidated(
    indices, geom, loc, pred_time, corrmodel_code, corrparam,
    distance_code, radius, context
  )
}

## geom is constructed once by .GeoValidatePairCacheGeometry().  Keep all
## per-target validation O(k), never rescan the N global observations here.
.GeoPairCrossCorrelationValidated <- function(indices, geom, loc,
                                               pred_time = 0, corrmodel_code,
                                               corrparam, distance_code, radius,
                                               context = "local kriging") {
  if (!is.numeric(indices) || any(!is.finite(indices)) ||
      any(indices != floor(indices)) || any(indices < 1) ||
      any(indices > nrow(geom$coords)))
    stop(context, ": invalid neighborhood indices.", call. = FALSE)
  if (!is.numeric(pred_time) || length(pred_time) != 1L || !is.finite(pred_time))
    stop(context, ": invalid prediction time.", call. = FALSE)
  loc <- as.numeric(loc)
  if (length(loc) != ncol(geom$coords) || any(!is.finite(loc))) {
    stop(context, ": invalid prediction coordinates for pair correlation.",
         call. = FALSE)
  }
  out <- .Call(
    "GeoPairCrossCorrelation",
    as.integer(indices), geom$coords, geom$times, loc, as.double(pred_time),
    as.integer(corrmodel_code), as.double(corrparam),
    as.integer(distance_code), as.double(radius), PACKAGE = "GeoModels"
  )
  if (length(out) != length(indices) || any(!is.finite(out))) {
    stop(context, ": invalid indexed cross-correlation result.", call. = FALSE)
  }
  as.numeric(out)
}

.GeoPairCorrelationParam <- function(corrmodel, param,
                                     context = "local kriging") {
  .GeoValidateCorrelationParameters(corrmodel, param, context = context)
}

.GeoGaussianPairCovariance <- function(indices, cache, param,
                                       context = "local kriging") {
  idx <- as.integer(indices)
  k <- length(idx)
  if (k < 1L) return(matrix(numeric(0L), 0L, 0L))
  packed <- .GeoPairCachePacked(idx, cache, context = context)
  R <- diag(1, k)
  if (k > 1L) {
    R[lower.tri(R)] <- packed
    R <- R + t(R) - diag(k)
  }
  nugget <- suppressWarnings(as.numeric(param[["nugget"]])[1L])
  sill <- suppressWarnings(as.numeric(param[["sill"]])[1L])
  if (!is.finite(nugget)) nugget <- 0
  if (!is.finite(sill) || sill <= 0) {
    stop(context, ": invalid Gaussian sill for pair-cache covariance.",
         call. = FALSE)
  }
  sill * ((1 - nugget) * R + nugget * diag(k))
}

.GeoSpatialCrossCorrelation <- function(coords_obs, loc,
                                        corrmodel_code, corrparam,
                                        distance_code, radius,
                                        context = "local spatial kriging") {
  coords_obs <- as.matrix(coords_obs)
  loc <- as.numeric(loc)
  if (ncol(coords_obs) == 2L) coords_obs <- cbind(coords_obs, 0)
  if (length(loc) == 2L) loc <- c(loc, 0)
  if (ncol(coords_obs) != 3L || length(loc) != 3L) {
    stop(context, ": incompatible coordinates for cross-correlation.",
         call. = FALSE)
  }
  storage.mode(coords_obs) <- "double"
  nobs <- nrow(coords_obs)
  ans <- dotCall64::.C64(
    "Corr_c",
    SIGNATURE = c(
      "double", "double", "double", "double", "double", "integer",
      "integer", "double", "double", "double", "integer", "integer",
      "integer", "integer", "integer", "integer", "double", "integer",
      "integer", "double", "integer", "integer", "double"
    ),
    corri = dotCall64::vector_dc("double", nobs),
    coords_obs[, 1L], coords_obs[, 2L], coords_obs[, 3L], 0,
    as.integer(corrmodel_code), 0L,
    as.double(loc[1L]), as.double(loc[2L]), as.double(loc[3L]),
    as.integer(nobs), 1L, 1L, as.integer(nobs), 0L, 1L,
    as.double(corrparam), 0L, 0L, 0,
    as.integer(distance_code), 0L, as.double(radius),
    INTENT = c("w", rep("r", 22)), NAOK = TRUE, PACKAGE = "GeoModels"
  )$corri
  if (length(ans) != nobs || any(!is.finite(ans))) {
    stop(context, ": non-finite spatial cross-correlation.", call. = FALSE)
  }
  as.numeric(ans)
}
