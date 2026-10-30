####################################################
### File name: GeoKrigloc.r
####################################################


###########################################################################
##  Future backend helpers                                                ##
###########################################################################

.GeoKriglocSetFuturePlan <- function(workers) {
  backend <- tolower(as.character(
    getOption("GeoModels.GeoKrigloc_future_backend", "auto")
  )[1L])
  if (!(backend %in% c("auto", "multicore", "multisession"))) {
    stop(
      "GeoKrigloc: option 'GeoModels.GeoKrigloc_future_backend' must be one of 'auto', 'multicore', or 'multisession'.",
      call. = FALSE
    )
  }

  can_multicore <- tryCatch(
    isTRUE(future::supportsMulticore()),
    error = function(e) FALSE
  )
  if ((identical(backend, "multicore") ||
       (identical(backend, "auto") && can_multicore)) && can_multicore) {
    future::plan(future::multicore, workers = workers)
    return("multicore")
  }
  if (identical(backend, "multicore") && !can_multicore) {
    warning(
      "GeoKrigloc: multicore futures are not supported in this R session; using multisession.",
      call. = FALSE
    )
  }

  libs <- .libPaths()
  pkg_dir <- tryCatch(find.package("GeoModels", quiet = TRUE),
                      error = function(e) character(0))
  if (length(pkg_dir) == 1L && nzchar(pkg_dir)) libs <- c(dirname(pkg_dir), libs)
  libs <- unique(vapply(
    libs[nzchar(libs)], normalizePath, character(1L),
    winslash = "/", mustWork = FALSE
  ))
  future::plan(future::multisession, workers = workers, rscript_libs = libs)
  "multisession"
}

###########################################################################
##  Top-level worker — defined OUTSIDE GeoKrigloc 
###########################################################################

.perform_kriging_worker <- function(k, neigh, param, corrmodel, distance,
                                    n, nloc, model, type_krig, sparse, anisopars,
                                    radius, mse, copula, varcov, method,
                                    loc, time, Xloc, Mloc, M,
                                    task_indices, is_space, spacetime,
                                    which, factor_cache = NULL, factor_keys = NULL) {
  .dup_old <- .GeoSuppressDuplicateChecksForNestedCalls()
  on.exit(.GeoRestoreDuplicateCheckOption(.dup_old), add = TRUE)

  if (is.environment(factor_cache) && length(factor_keys) >= k &&
      !is.na(factor_keys[k])) {
    old_factor_options <- options(
      GeoModels.local_krig_factor_cache = factor_cache,
      GeoModels.local_krig_factor_key = factor_keys[k])
    on.exit(options(old_factor_options), add = TRUE)
  }


  if (!is.null(neigh$numpoints) && length(neigh$numpoints) >= k &&
      neigh$numpoints[k] == 0L) {
    out <- list(indices = as.integer(k), pred = NA_real_,
                mse = if (mse) NA_real_ else NULL)
    if (identical(model, "Wrapped") && is.null(copula)) {
      out$wrapped_mse_sin <- if (mse) NA_real_ else NULL
      out$wrapped_mse_cos <- if (mse) NA_real_ else NULL
      out$wrapped_resultant <- NA_real_
    }
    return(out)
  }

  local_param <- param
  if (!is.null(M)) local_param$mean <- neigh$M[[k]]

  Xloc_k <- if (!is.null(Xloc)) Xloc[k, , drop = FALSE] else NULL

  Mloc_k <- if (!is.null(Mloc)) Mloc[k] else NULL
  nloc_k <- if (is.null(nloc) || length(nloc) == 1L) nloc else nloc[k]
  n_k <- n
  if (length(n) > 1L && !is.null(neigh$indices) && length(neigh$indices) >= k) {
    idx_n <- neigh$indices[[k]]
    if (length(n) >= max(idx_n, 0L)) n_k <- n[idx_n]
  }

  if (is_space) {
    pr <- GeoKrig(
      loc       = loc[k, ],
      data      = neigh$data[[k]],
      coordx    = neigh$coordx[[k]],
      corrmodel = corrmodel,
      distance  = distance,
      n         = n_k,
      nloc      = nloc_k,
      method    = method,
      X         = neigh$X[[k]],
      Xloc      = Xloc_k,
      Mloc      = Mloc_k,
      type_krig = type_krig,
      sparse    = sparse,
      model     = model,
      param     = local_param,
      anisopars = anisopars,
      radius    = radius,
      mse       = mse,
      copula    = copula,
      varcov    = varcov
    )

  } else {
    i <- task_indices[k, 1L]
    j <- task_indices[k, 2L]

    if (spacetime) {
      coordx_dyn_k <- if (!is.null(neigh$coordx_dyn)) neigh$coordx_dyn[[k]] else NULL
      pr <- GeoKrig(
        data      = neigh$data[[k]],
        coordx    = if (is.null(coordx_dyn_k)) neigh$coordx[[k]] else NULL,
        coordt    = neigh$coordt[[k]],
        coordx_dyn = coordx_dyn_k,
        loc       = loc[i, ],
        time      = time[j],
        X         = neigh$X[[k]],
        Mloc      = Mloc_k,
        Xloc      = Xloc_k,
        type_krig = type_krig,
        sparse    = sparse,
        corrmodel = corrmodel,
        distance  = distance,
        model     = model,
        param     = local_param,
        anisopars = anisopars,
        radius    = radius,
        mse       = mse,
        varcov    = varcov,
        copula    = copula,
        n         = n_k,
        nloc      = nloc_k,
        method    = method
      )
    } else {
      # bivariate: preserve the two variable-specific coordinate blocks.
      coordx_dyn_k <- if (!is.null(neigh$coordx_dyn)) neigh$coordx_dyn[[k]] else NULL
      if (!is.null(neigh$block_sizes) && any(neigh$block_sizes[[k]] == 0L)) {
        out <- list(indices = as.integer(k), pred = NA_real_,
                    mse = if (mse) NA_real_ else NULL)
        if (identical(model, "Wrapped") && is.null(copula)) {
          out$wrapped_mse_sin <- if (mse) NA_real_ else NULL
          out$wrapped_mse_cos <- if (mse) NA_real_ else NULL
          out$wrapped_resultant <- NA_real_
        }
        return(out)
      }
      pr <- GeoKrig(
        data      = neigh$data[[k]],
        coordx    = if (is.null(coordx_dyn_k)) neigh$coordx[[k]] else NULL,
        coordx_dyn = coordx_dyn_k,
        loc       = loc[i, ],
        X         = neigh$X[[k]],
        Mloc      = Mloc_k,
        Xloc      = Xloc_k,
        type_krig = type_krig,
        sparse    = sparse,
        corrmodel = corrmodel,
        distance  = distance,
        model     = model,
        param     = local_param,
        varcov    = varcov,
        anisopars = anisopars,
        radius    = radius,
        mse       = mse,
        copula    = copula,
        n         = n_k,
        nloc      = nloc_k,
        method    = method,
        which     = which
      )
    }
  }

  out <- list(indices = as.integer(k), pred = pr$pred,
              mse = if (mse) pr$mse else NULL)
  if (identical(model, "Wrapped") && is.null(copula)) {
    out$wrapped_mse_sin <- if (mse) pr$wrapped_mse_sin else NULL
    out$wrapped_mse_cos <- if (mse) pr$wrapped_mse_cos else NULL
    out$wrapped_resultant <- pr$wrapped_resultant
  }
  out
}


###########################################################################
##  Main function                                                        ##
###########################################################################

GeoKrigloc = function(estobj=NULL, data, coordx=NULL, coordy=NULL, coordz=NULL, coordt=NULL, coordx_dyn=NULL,
                      corrmodel, distance="Eucl", grid=FALSE, loc, neighb=NULL,
                      maxdist=NULL, maxtime=NULL, method="cholesky", model="Gaussian",
                      n=1, nloc=NULL, mse=FALSE, param, anisopars=NULL,
                      radius=1, sparse=FALSE, time=NULL, type="Standard",
                      type_krig="Simple", weigthed=TRUE, which=1,
                      copula=NULL, X=NULL, Xloc=NULL, Mloc=NULL, varcov=NULL,
                      spobj=NULL, spdata=NULL, parallel=FALSE, ncores=6, progress=TRUE,
                      check.duplicates=FALSE)
{


  ###########################################################################
  ##  0. Parallel helpers                                                  ##
  ###########################################################################

  setup_parallel <- function(workers) {
    old_opts <- options(
      warn = -1,
      verbose = FALSE,
      future.globals.maxSize = 4 * 1024^3
    )
    old_plan <- future::plan()
    .GeoKriglocSetFuturePlan(workers)
    list(old_opts = old_opts, old_plan = old_plan)
  }

  cleanup_parallel <- function(cfg) {
    options(cfg$old_opts)
    try(future::plan(cfg$old_plan), silent = TRUE)
    invisible(gc(verbose = FALSE))
  }

  choose_workers_safe <- function(neigh, max_cores, n_jobs,
                                      payload = neigh) {
    max_cores <- min(as.integer(max_cores), max(1L, as.integer(n_jobs)))
    if (max_cores <= 1L) return(1L)

    nn <- if (!is.null(neigh$numpoints)) {
      suppressWarnings(max(as.integer(neigh$numpoints), na.rm = TRUE))
    } else if (!is.null(neigh$indices)) {
      suppressWarnings(max(vapply(neigh$indices, length, integer(1L)),
                           na.rm = TRUE))
    } else if (!is.null(neigh$coordx)) {
      suppressWarnings(max(vapply(neigh$coordx, function(z) {
        if (is.null(z)) 0L else nrow(as.matrix(z))
      }, integer(1L)), na.rm = TRUE))
    } else 0L
    if (!is.finite(nn) || nn < 1L) nn <- 1L

    p <- if (!is.null(neigh$X)) {
      suppressWarnings(max(vapply(neigh$X, function(z) {
        if (is.null(z)) 0L else ncol(as.matrix(z))
      }, integer(1L)), na.rm = TRUE))
    } else 0L
    if (!is.finite(p) || p < 0L) p <- 0L

    ## Conservative dense-local estimate: covariance, factor/solve workspaces,
    ## cross-covariance/design temporaries and R object overhead.
    nn <- as.double(nn)
    p <- as.double(p)
    bytes_per_job <- 8 * (nn^2 + nn * max(1, p) + 8 * nn) * 4
    if (!is.finite(bytes_per_job) || bytes_per_job <= 0) return(1L)

    ## Multisession workers each receive their own neighborhood object;
    ## multicore workers can initially share it through copy-on-write.
    backend <- tolower(as.character(
      getOption("GeoModels.GeoKrigloc_future_backend", "auto")
    )[1L])
    can_multicore <- tryCatch(
      isTRUE(future::supportsMulticore()), error = function(e) FALSE
    )
    multisession <- identical(backend, "multisession") ||
      (identical(backend, "auto") && !can_multicore)
    if (multisession) {
      bytes_per_job <- bytes_per_job + as.numeric(object.size(payload))
    }

    free_ram <- tryCatch({
      if (requireNamespace("memuse", quietly = TRUE)) {
        as.numeric(memuse::Sys.meminfo()["freeram"], units = "b")
      } else 2 * 1024^3
    }, error = function(e) 2 * 1024^3)

    safe <- max(1L, as.integer(floor(0.70 * free_ram / bytes_per_job)))
    min(max_cores, safe)
  }


  ###########################################################################
  ##  1. Input processing                                                  ##
  ###########################################################################

  call = match.call()

  if (!is.null(estobj)) {
    if (!inherits(estobj, "GeoFit")) stop("need a 'GeoFit' object as input\n")
    data = estobj$data
    if (!estobj$grid) {
      if (!estobj$bivariate) {
        if (is.null(estobj$coordx_dyn))
          coordx = cbind(estobj$coordx, estobj$coordy, estobj$coordz)
        else
          coordx_dyn = estobj$coordx_dyn
      } else {
        if (is.null(estobj$coordx_dyn)) {
          coordx = estobj$coordx[1:estobj$ns[1]]
          coordy = estobj$coordy[1:estobj$ns[2]]
        } else {
          coordx_dyn = estobj$coordx_dyn
        }
      }
    } else {
      coordx = estobj$coordx
      coordy = estobj$coordy
      coordz = estobj$coordz
    }
    if (isTRUE(estobj$spacetime)) coordt = estobj$coordt else coordt = NULL
    coordx_dyn = estobj$coordx_dyn
    corrmodel  = estobj$corrmodel
    model      = estobj$model
    distance   = estobj$distance
    grid       = estobj$grid
    n          = estobj$n
    param      = append(estobj$param, estobj$fixed)
    radius     = estobj$radius
    copula     = estobj$copula
    anisopars  = estobj$anisopars
    X          = estobj$X
    varcov     = estobj$varcov
    if (!is.null(estobj$sparse)) sparse = estobj$sparse
  }

  model_requested <- gsub("[[:blank:]]", "", as.character(model)[1L])
  model <- .GeoPredictionCanonicalModel(model_requested)
  if (is.null(CkModel(model))) stop("The name of the model is not correct")
  if (!is.character(corrmodel) || is.null(CkCorrModel(corrmodel)))
    stop("The name of the correlation model is wrong")
  corrmodel <- gsub("[[:blank:]]", "", corrmodel)
  if (!is.null(copula)) {
    copula <- as.character(copula)[1L]
    if (is.na(copula) || !nzchar(gsub("[[:blank:]]", "", copula))) {
      copula <- NULL
    } else {
      copula <- gsub("[[:blank:]]", "", copula)
    }
  }
  if (!is.null(copula) && !(copula %in% c("Gaussian", "Clayton", "SkewGaussian")))
    stop("GeoKrigloc currently implements copula-based linear prediction for 'Gaussian', 'Clayton', and 'SkewGaussian'.", call. = FALSE)
  distance <- gsub("[[:blank:]]", "", distance)

  bivariate = CheckBiv(CkCorrModel(corrmodel))
  spacetime = CheckST(CkCorrModel(corrmodel))
  space     = !spacetime && !bivariate
  .GeoKrigValidatePredictionModel(model, copula = copula, bivariate = bivariate, context = "GeoKrigloc")

  method <- .GeoKrigMatchMethod(method, sparse = sparse, context = "GeoKrigloc")
  type_key <- tolower(gsub("[[:blank:]]", "", type))
  if (!identical(type_key, "standard")) {
    warning("GeoKrigloc currently computes the Standard local kriging system; ",
            "the supplied type is ignored.", call. = FALSE)
  }
  type <- "Standard"
  krig_key <- tolower(gsub("[[:blank:]]", "", type_krig))
  if (!(krig_key %in% c("simple", "universal")))
    stop("type_krig must be 'Simple' or 'Universal'")
  type_krig <- if (krig_key == "simple") "Simple" else "Universal"
  if (identical(model, "Wrapped") && identical(type_krig, "Universal")) {
    warning("For model = 'Wrapped', prediction uses componentwise Simple kriging of sine/cosine residuals; type_krig = 'Universal' is replaced by 'Simple'.",
            call. = FALSE)
    type_krig <- "Simple"
  }
  if (!is.logical(weigthed) || length(weigthed) != 1L || is.na(weigthed))
    stop("weigthed must be TRUE or FALSE")
  if (!isTRUE(weigthed))
    warning("weigthed is retained for backward compatibility and has no effect for type = 'Standard'.", call. = FALSE)

  ## GeoFit stores a single spatial realization as a 1 x n matrix.
  ## GeoNeighborhood expects spatial data as a vector (or an n-row object).
  ## Normalize only this representation; space-time and bivariate formats
  ## must retain their matrix/list structure.
  if (space && is.matrix(data) && nrow(data) == 1L) {
    data <- as.numeric(data[1L, , drop = TRUE])
  }

  if (!is.null(spobj)) {
    a = sp2Geo(spobj, spdata); coordx = a$coords
    if (!a$pj && distance != "Chor") distance = "Geod"
    if (spacetime) coordt = a$coordt
    if (!is.null(a$Y)) data <- a$Y
   if (!is.null(a$X)) X <- a$X
  }

  if (spacetime && !is.null(coordx_dyn)) {
    coordx_dyn <- .GeoNormalizeDynamicCoordinates(
      coordx_dyn, coordt = coordt, context = "GeoKrigloc"
    )
    .GeoValidateDynamicData(data, coordx_dyn, context = "GeoKrigloc")
    .GeoValidateDynamicRows(X, coordx_dyn, name = "X", context = "GeoKrigloc")
  }

  biv_info <- NULL
  if (bivariate) {
    biv_info <- .GeoBivariateCoordinateInfo(
      coordx = coordx, coordy = coordy, coordz = coordz,
      coordx_dyn = coordx_dyn, grid = grid,
      context = "GeoKrigloc"
    )
    data <- .GeoBivariateDataVector(data, biv_info$ns,
                                    context = "GeoKrigloc")
    if (!is.null(X)) {
      X <- .GeoBivariateDesign(X, biv_info$ns,
                               context = "GeoKrigloc")
    }
    if (biv_info$dynamic) {
      coordx_dyn <- biv_info$blocks
      coordx <- coordy <- coordz <- NULL
    }
  }

  if (is.null(coordx_dyn)) {
    coords <- .GeoKrigCoordinateMatrix(
      coordx = coordx, coordy = coordy, coordz = coordz, grid = grid,
      context = "GeoKrigloc"
    )
    if (!bivariate) {
      data <- .GeoKrigNormalizeData(
        data, nspace = nrow(coords),
        coordt = if (spacetime) coordt else NULL,
        grid = grid, context = "GeoKrigloc"
      )
    }
  } else {
    coordx = coordy = coordz = coords = NULL
  }

  .GeoCheckDuplicateLocations(
    coordx = if (is.null(coordx_dyn)) coords else NULL,
    coordt = coordt, coordx_dyn = coordx_dyn,
    spacetime = spacetime, bivariate = bivariate, grid = FALSE,
    context = "GeoKrigloc", check.duplicates = check.duplicates
  )
  .dup_old <- .GeoSuppressDuplicateChecksForNestedCalls()
  on.exit(.GeoRestoreDuplicateCheckOption(.dup_old), add = TRUE)

  if (is.vector(loc) && !is.list(loc)) loc <- matrix(loc, nrow = 1L)
  loc <- as.matrix(loc)
  loc_orig <- loc
  if (!is.numeric(loc) || !nrow(loc) || any(!is.finite(loc)))
    stop("loc must be a finite numeric N x 2 or N x 3 matrix")
  if (spacetime) {
    if (is.null(time) || !is.numeric(time) || !length(time) || any(!is.finite(time)))
      stop("time must be a non-empty finite numeric vector for space-time kriging")
    time <- as.numeric(time)
  }

  .st_old <- .GeoSTEnter(corrmodel,
    coordx = if (is.null(coordx_dyn)) coords else NULL,
    coordx_dyn = coordx_dyn, distance = distance, context = "GeoKrigloc")
  on.exit(.GeoSTRestore(.st_old), add = TRUE)
  .GeoSTValidate(corrmodel, param, context = "GeoKrigloc")

  .GeoValidate3D(
    coordx = if (is.null(coordx_dyn)) coords else NULL,
    coordx_dyn = coordx_dyn, distance = distance, grid = grid,
    context = "GeoKrigloc"
  )
  .GeoValidate3D(
    coordx = loc, distance = distance, grid = FALSE,
    context = "GeoKrigloc prediction locations"
  )
  .GeoValidateCoordinateCompatibility(
    coordx = if (is.null(coordx_dyn)) coords else NULL,
    coordx_dyn = coordx_dyn, loc = loc, context = "GeoKrigloc"
  )

  ## Select neighbors and compute covariances in exactly the same metric.
  coords_metric <- coords
  coordx_dyn_metric <- coordx_dyn
  loc_metric <- loc
  worker_anisopars <- anisopars
  if (!is.null(anisopars)) {
    if (!identical(distance, "Eucl"))
      stop("GeoKrigloc: anisopars requires distance = 'Eucl'.")
    if (!is.null(coords_metric))
      coords_metric <- GeoAniso(coords_metric, c(anisopars$angle, anisopars$ratio))
    if (!is.null(coordx_dyn_metric))
      coordx_dyn_metric <- lapply(coordx_dyn_metric, function(z)
        GeoAniso(as.matrix(z), c(anisopars$angle, anisopars$ratio)))
    loc_metric <- GeoAniso(loc_metric, c(anisopars$angle, anisopars$ratio))
    ## Coordinates sent to GeoKrig are already transformed.
    worker_anisopars <- NULL
    param$angle <- NULL
    param$ratio <- NULL
  }

  Nloc = nrow(loc)
  if (is.null(Nloc)) Nloc = 1
  if (spacetime && is.null(time))
    stop("At least one prediction time is required for space-time kriging.")
  Tloc = length(time)
  if (bivariate) Tloc = 1

  ## Normalize the same three mean specifications used by GeoKrig():
  ## intercept only, X beta, or an external observation mean plus Mloc.
  if (!bivariate) {
    nobs_mean <- length(unlist(data, use.names = FALSE))
    M <- .GeoMean_external(param, n = nobs_mean)
    if (!is.null(M) && !is.null(X))
      stop("A fixed observation mean and X cannot be used together.")
    if (!is.null(Mloc) && !is.null(Xloc))
      stop("Use either Mloc or Xloc, not both.")
    if (!is.null(M) && is.null(Mloc))
      stop("Mloc is required when the observation mean is supplied as a vector.")
    if (!is.null(X)) {
      X <- .GeoMean_design(X, n = nobs_mean, name = "X")
      p_mean <- ncol(X)
      invisible(.GeoMean_beta(param, p_mean))
      if (is.null(Xloc) && is.null(Mloc))
        stop("Xloc or Mloc is required when X is supplied.")
    } else {
      p_mean <- 1L
      if (!is.null(Xloc)) stop("X must be supplied when Xloc is used.")
      if (is.null(M)) {
        if (length(setdiff(.GeoMean_supplied_names(param), "mean")))
          stop("Parameters mean1, mean2, ... require an explicit design matrix X; without X the mean is intercept-only.")
        invisible(.GeoMean_beta(param, 1L))
      }
    }

    total_mean_tasks <- if (space) Nloc else Nloc * Tloc
    if (!is.null(Mloc)) {
      if (!is.numeric(Mloc) || any(!is.finite(Mloc)))
        stop("Mloc must be finite and numeric.")
      Mloc <- as.numeric(Mloc)
      if (length(Mloc) == 1L) Mloc <- rep(Mloc, total_mean_tasks)
      if (length(Mloc) != total_mean_tasks)
        stop("Mloc must have one value per local prediction task.")
      Xloc <- NULL
    } else if (!is.null(Xloc)) {
      Xloc <- .GeoMean_design(Xloc, n = total_mean_tasks, p = p_mean, name = "Xloc")
    }
  } else {
    M <- NULL
    p_mean <- if (is.null(X)) 1L else ncol(X)
    invisible(.GeoBivariateValidateMeanParameters(
      param, p_mean, context = "GeoKrigloc"
    ))
    if (!is.null(Mloc) && !is.null(Xloc))
      stop("Use either Mloc or Xloc, not both.")
    if (!is.null(X) && is.null(Xloc) && is.null(Mloc))
      stop("Xloc or Mloc is required when X is supplied.")
    if (is.null(X) && !is.null(Xloc))
      stop("X must be supplied when Xloc is used.")
    if (!is.null(Mloc)) {
      if (!is.numeric(Mloc) || any(!is.finite(Mloc)))
        stop("Mloc must be finite and numeric.")
      Mloc <- as.numeric(Mloc)
      if (length(Mloc) == 1L) Mloc <- rep(Mloc, Nloc)
      if (length(Mloc) != Nloc)
        stop("Mloc must have one value per prediction location.")
      Xloc <- NULL
    } else if (!is.null(Xloc)) {
      Xloc <- .GeoMean_design(Xloc, n = Nloc, p = p_mean, name = "Xloc")
    }
  }

  if (identical(type_krig, "Universal") && (!is.null(M) || !is.null(Mloc))) {
    warning("External observation/prediction means are treated as known; ",
            "type_krig is therefore set to 'Simple'.", call. = FALSE)
    type_krig <- "Simple"
  }
  if (identical(type_krig, "Universal") && isTRUE(mse) && is.null(varcov)) {
    stop("Universal kriging with mse = TRUE requires varcov (the covariance matrix of the estimated parameters).")
  }


  ###########################################################################
  ##  2. Pre-compute task indices for spacetime / bivariate               ##
  ###########################################################################

  task_indices <- if (spacetime || bivariate) {
    idx <- matrix(0L, nrow = Nloc * Tloc, ncol = 2L)
    k   <- 1L
    for (i in seq_len(Nloc)) {
      for (j in seq_len(Tloc)) {
        idx[k, 1L] <- i
        idx[k, 2L] <- j
        k <- k + 1L
      }
    }
    idx
  } else {
    NULL
  }


  ###########################################################################
  ##  3. Neighbourhood setup                                               ##
  ###########################################################################

  if (space) {
    neigh       = GeoNeighborhood(data, coordx = coords_metric, distance = distance,
                                  loc = loc_metric, neighb = neighb, maxdist = maxdist,
                                  X = X, M = M, parallel = FALSE, ncores = ncores)
    total_tasks <- Nloc

  } else if (spacetime) {
    neigh       = GeoNeighborhood(data, coordx = coords_metric, coordt = coordt,
                                  distance = distance, neighb = neighb,
                                  coordx_dyn = coordx_dyn_metric, loc = loc_metric, time = time,
                                  maxdist = maxdist, maxtime = maxtime,
                                  X = X, M = M, parallel = FALSE, ncores = ncores)
    total_tasks <- Nloc * Tloc

  } else {   # bivariate
    neigh       = GeoNeighborhood(data, coordx = coords_metric, distance = distance,
                                  neighb = neighb, coordx_dyn = coordx_dyn_metric,
                                  bivariate = TRUE,
                                  loc = loc_metric, maxdist = maxdist,
                                  X = X, M = M, parallel = FALSE, ncores = ncores)
    total_tasks <- Nloc * Tloc
  }


  ###########################################################################
  ##  4. Progress bar setup                                                ##
  ###########################################################################

  if (!is.logical(progress)) stop("progress must be logical (TRUE/FALSE)")

  if (progress) {
    restore_progress_handlers <- .GeoProgressHandlersPush(
      handler = "txtprogressbar", global = TRUE
    )
    on.exit(restore_progress_handlers(), add = TRUE)
    pb <- progressr::progressor(
      along = seq_len(total_tasks), on_exit = FALSE, auto_finish = FALSE
    )
    ## Finish before restoring the temporary global handler, including on
    ## errors or user interrupts.
    on.exit(pb(type = "finish"), add = TRUE, after = FALSE)
  } else {
    pb <- function(...) {}
  }

  # Pre-allocate result vectors
  res1 <- numeric(total_tasks)
  res2 <- if (mse) numeric(total_tasks) else NULL
  wrapped_local <- identical(model, "Wrapped") && is.null(copula)
  res_wrap_sin <- if (wrapped_local && mse) numeric(total_tasks) else NULL
  res_wrap_cos <- if (wrapped_local && mse) numeric(total_tasks) else NULL
  res_wrap_resultant <- if (wrapped_local) numeric(total_tasks) else NULL


  ###########################################################################
  ##  5. Kriging jobs                                                      ##
  ###########################################################################

  ## GeoKrigloc deliberately uses GeoKrig as its single local prediction
  ## engine.  GeoNeighborhood defines each local data set; GeoKrig then
  ## performs exactly the same validated kriging calculation used elsewhere
  ## in the package.  Keeping one numerical path avoids duplicating covariance,
  ## copula, mean, MSE, and solver logic inside GeoKrigloc.
  jobs <- seq_len(total_tasks)
  factor_cache <- .GeoKrigNewFactorCache()
  factor_keys <- NULL
  if (!is.null(factor_cache) && length(neigh$indices) == total_tasks) {
    factor_keys <- vapply(neigh$indices, function(ii) paste(ii, collapse = ","), character(1L))
    repeated <- duplicated(factor_keys) | duplicated(factor_keys, fromLast = TRUE)
    factor_keys[!repeated | !nzchar(factor_keys)] <- NA_character_
    if (all(is.na(factor_keys))) factor_cache <- NULL
  } else {
    factor_cache <- NULL
  }
  worker_args <- list(
    neigh = neigh, param = param, corrmodel = corrmodel,
    distance = distance, n = n, nloc = nloc, model = model,
    type_krig = type_krig, sparse = sparse,
    anisopars = worker_anisopars, radius = radius, mse = mse,
    copula = copula, varcov = varcov, method = method,
    loc = loc_metric, time = time, Xloc = Xloc, Mloc = Mloc, M = M,
    task_indices = task_indices, is_space = space,
    spacetime = spacetime, which = which,
    factor_cache = factor_cache, factor_keys = factor_keys
  )


  ###########################################################################
  ##  6. Execute jobs                                                       ##
  ###########################################################################

  parallel_exec <- isTRUE(parallel) && total_tasks > 1L
  workers <- 1L
  if (parallel_exec) {
    requested <- .GeoResolveWorkers(TRUE, ncores, n_jobs = total_tasks)

    ## A numeric ncores value is an explicit worker request.  Respect it up to
    ## the number of detected cores and available prediction tasks.  The RAM
    ## safety heuristic is reserved for ncores = NULL, which remains available
    ## as an opt-in automatic/safe mode.
    if (is.null(ncores)) {
      workers <- choose_workers_safe(
        neigh, requested, total_tasks, payload = worker_args
      )
      if (workers < requested) {
        warning("Reducing cores from ", requested, " to ", workers,
                " to avoid memory pressure.")
      }
    } else {
      workers <- requested
    }
    parallel_exec <- workers > 1L
  }

  if (parallel_exec) {
    config <- setup_parallel(workers)
    on.exit(cleanup_parallel(config), add = TRUE)
    if (progress) cat("Performing local kriging using", workers, "cores...\n")

    res <- future.apply::future_lapply(
      jobs,
      function(k) {
        worker <- get(
          ".perform_kriging_worker",
          envir = asNamespace("GeoModels"),
          inherits = FALSE
        )
        ans <- do.call(worker, c(list(k = as.integer(k)), worker_args))
        if (progress) pb(sprintf("k=%d", k))
        ans
      },
      future.globals = c("worker_args", "progress", "pb"),
      future.packages = "GeoModels",
      future.seed = TRUE,
      future.stdout = FALSE,
      future.conditions = if (progress) "condition" else "none"
    )
  } else {
    res <- vector("list", total_tasks)
    for (k in jobs) {
      ans <- do.call(
        .perform_kriging_worker,
        c(list(k = as.integer(k)), worker_args)
      )
      res[[k]] <- ans
      if (progress) pb(sprintf("k=%d", k))
    }
  }

  for (ans in res) {
    res1[ans$indices] <- as.numeric(ans$pred)
    if (mse) res2[ans$indices] <- as.numeric(ans$mse)
    if (wrapped_local) {
      res_wrap_resultant[ans$indices] <- as.numeric(ans$wrapped_resultant)
      if (mse) {
        res_wrap_sin[ans$indices] <- as.numeric(ans$wrapped_mse_sin)
        res_wrap_cos[ans$indices] <- as.numeric(ans$wrapped_mse_cos)
      }
    }
  }


  ###########################################################################
  ##  7. Format output                                                     ##
  ###########################################################################

  varpred <- NULL
  wrapped_mse_sin <- wrapped_mse_cos <- wrapped_resultant <- NULL

  if (spacetime || bivariate) {
    pred <- matrix(res1, nrow = Tloc, ncol = Nloc)
    if (mse) varpred <- matrix(res2, nrow = Tloc, ncol = Nloc)
    if (wrapped_local) {
      wrapped_resultant <- matrix(res_wrap_resultant, nrow = Tloc, ncol = Nloc)
      if (mse) {
        wrapped_mse_sin <- matrix(res_wrap_sin, nrow = Tloc, ncol = Nloc)
        wrapped_mse_cos <- matrix(res_wrap_cos, nrow = Tloc, ncol = Nloc)
      }
    }
  } else {
    pred <- res1
    if (mse) varpred <- res2
    if (wrapped_local) {
      wrapped_resultant <- res_wrap_resultant
      if (mse) {
        wrapped_mse_sin <- res_wrap_sin
        wrapped_mse_cos <- res_wrap_cos
      }
    }
  }

  if (Tloc == 1) {
    pred    <- c(pred)
    varpred <- c(varpred)
    wrapped_resultant <- c(wrapped_resultant)
    wrapped_mse_sin <- c(wrapped_mse_sin)
    wrapped_mse_cos <- c(wrapped_mse_cos)
  }

  GeoKrigloc_out <- list(
    bivariate = bivariate,
    coordx    = coordx,
    coordy    = coordy,
    coordz    = coordz,
    coordt    = coordt,
    coordx_dyn = coordx_dyn,
    corrmodel = corrmodel,
    data      = data,
    distance  = distance,
    grid      = grid,
    loc       = loc_orig,
    copula    = copula,
    numcoord  = if (exists("coords", inherits = FALSE)) nrow(coords) else NA,
    numloc    = Nloc,
    numtime   = length(coordt),
    numt      = Tloc,
    maxdist   = maxdist,
    maxtime   = maxtime,
    model     = model,
    fit_model = model_requested,
    n         = n,
    param     = param,
    pred      = pred,
    radius    = radius,
    spacetime = spacetime,
    time      = time,
    type       = type,
    type_krig = type_krig,
    method     = method,
    weigthed   = weigthed,
    which     = which,
    mse       = varpred,
    varcov    = varcov
  )

  if (wrapped_local) {
    GeoKrigloc_out$wrapped_mse_sin <- wrapped_mse_sin
    GeoKrigloc_out$wrapped_mse_cos <- wrapped_mse_cos
    GeoKrigloc_out$wrapped_resultant <- wrapped_resultant
    GeoKrigloc_out$wrapped_prediction_method <-
      "componentwise linear prediction of sine/cosine residuals"
  }

  structure(c(GeoKrigloc_out, call = call), class = c("GeoKrigloc"))
}
