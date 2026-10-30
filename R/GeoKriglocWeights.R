###########################################################################
##  Internal covariance helpers used by GeoKriglocWeights               ##
###########################################################################

###########################################################################
##  Shared dense spatial covariance helpers                              ##
###########################################################################

.GeoKriglocPrepareSpatialCovariance <- function(neigh, coordx_global,
                                                 X_global, param, corrmodel,
                                                 distance, n, model,
                                                 radius, copula, M = NULL) {
  if (is.null(neigh$indices)) return(NULL)
  nonempty <- which(vapply(neigh$indices, length, integer(1L)) >= 2L)
  if (!length(nonempty)) return(NULL)

  k0 <- nonempty[1L]
  idx0 <- as.integer(neigh$indices[[k0]])
  take <- idx0[seq_len(min(2L, length(idx0)))]
  coords0 <- as.matrix(coordx_global[take, , drop = FALSE])
  X0 <- if (!is.null(X_global)) {
    as.matrix(X_global[take, , drop = FALSE])
  } else NULL

  param0 <- param
  ## External observation means are location-specific and are supplied to the
  ## direct covariance assembler at run time.  The one-time template only
  ## needs a valid scalar mean so GeoCovmatrix can normalize metadata exactly
  ## as in the historical path.
  if (!is.null(M)) param0$mean <- 0

  n0 <- n
  if (length(n) > 1L && length(take) && length(n) >= max(take)) {
    n0 <- n[take]
  }

  template <- GeoCovmatrix(
    coordx = coords0, corrmodel = corrmodel, distance = distance,
    model = model, n = n0, param = param0, anisopars = NULL,
    radius = radius, sparse = FALSE, copula = copula, X = X0
  )
  if (!is.null(template$copula)) {
    cp <- as.character(template$copula)[1L]
    if (is.na(cp) || !nzchar(gsub("[[:blank:]]", "", cp))) {
      template$copula <- NULL
    }
  }
  template$param <- unlist(template$param)
  template$local_cov_radius <- if (CheckSph(CkCorrModel(corrmodel))) 1 else radius
  template
}

.GeoKriglocLocalParamVector <- function(template, param) {
  out <- template$param
  z <- unlist(param)
  if (length(z) && !is.null(names(z))) {
    keep <- intersect(names(z), names(out))
    if (length(keep)) out[keep] <- z[keep]
  }
  out
}

## Assemble the Gaussian covariance from the already cached packed pairs.
## This helper is deliberately Gaussian-only, matching the public weights API.
.GeoKriglocSpatialDenseCovariance <- function(nobs, nuisance, raw_corr,
                                               context = "GeoKriglocWeights") {
  npairs <- as.double(nobs) * as.double(nobs - 1L) / 2
  if (length(raw_corr) != npairs || any(!is.finite(raw_corr)))
    stop(context, ": invalid packed local correlation vector.", call. = FALSE)
  nugget <- as.numeric(nuisance["nugget"])[1L]
  sill <- as.numeric(nuisance["sill"])[1L]
  if (!is.finite(nugget)) nugget <- 0
  if (!is.finite(sill) || sill <= 0)
    stop(context, ": invalid Gaussian sill parameter.", call. = FALSE)
  corr <- as.numeric(raw_corr) * (1 - nugget)
  C <- diag(nobs)
  if (length(corr)) {
    C[lower.tri(C)] <- corr
    C <- t(C)
    C[lower.tri(C)] <- corr
  }
  C * sill
}

## Fast compact spatial Gaussian weights using the shared observation-pair
## correlation cache.  The cross-correlation to the prediction location is
## still evaluated directly because those pairs are location-specific.
.GeoKriglocWeightsCachedSpatial <- function(idx, i, vars) {
  coords_obs <- as.matrix(vars$coords_global[idx, , drop = FALSE])
  raw_corr <- .GeoPairCachePacked(
    idx, vars$pair_cache, context = "GeoKriglocWeights"
  )
  nuisance <- vars$cache_param_vec[vars$cov_template$namesnuis]
  nuisance <- nuisance[!is.na(names(nuisance))]
  C <- .GeoKriglocSpatialDenseCovariance(
    nobs = length(idx), nuisance = nuisance, raw_corr = raw_corr,
    context = "GeoKriglocWeights"
  )

  cross_corr <- .GeoSpatialCrossCorrelation(
    coords_obs, vars$loc_metric[i, ],
    corrmodel_code = vars$corrmodel_code,
    corrparam = vars$cache_param_vec[vars$cov_template$namescorr],
    distance_code = vars$distance_code,
    radius = vars$cov_template$local_cov_radius,
    context = "GeoKriglocWeights"
  )
  nugget <- as.numeric(nuisance["nugget"])[1L]
  sill <- as.numeric(nuisance["sill"])[1L]
  if (!is.finite(nugget)) nugget <- 0
  if (!is.finite(sill))
    stop("GeoKriglocWeights: invalid sill in cached spatial core.", call. = FALSE)
  CC <- matrix((1 - nugget) * sill * cross_corr, ncol = 1L)
  as.numeric(.GeoKrigDenseSPDSolve(C, CC, context = "GeoKriglocWeights"))
}

## Compact Gaussian space-time weights from the shared pair cache.  Global
## observation coordinates/times are flattened once in neighborhood-index
## order, so static and dynamic space-time data use the same kernel.
.GeoKriglocWeightsCachedSpacetime <- function(idx, loc, pred_time, vars) {
  C <- .GeoGaussianPairCovariance(
    idx, vars$pair_cache, vars$param, context = "GeoKriglocWeights"
  )
  cross_corr <- .GeoPairCrossCorrelationValidated(
    idx, vars$st_geometry,
    loc = loc, pred_time = pred_time,
    corrmodel_code = vars$corrmodel_code,
    corrparam = vars$corrparam,
    distance_code = vars$distance_code,
    radius = vars$radius,
    context = "GeoKriglocWeights"
  )
  nugget <- suppressWarnings(as.numeric(vars$param[["nugget"]])[1L])
  sill <- suppressWarnings(as.numeric(vars$param[["sill"]])[1L])
  if (!is.finite(nugget)) nugget <- 0
  if (!is.finite(sill) || sill <= 0)
    stop("GeoKriglocWeights: invalid sill in cached space-time core.", call. = FALSE)
  CC <- matrix((1 - nugget) * sill * cross_corr, ncol = 1L)
  as.numeric(.GeoKrigDenseSPDSolve(C, CC, context = "GeoKriglocWeights"))
}

## -------------------------------------------------------------------------
## Internal worker shared by serial and parallel spatial/bivariate weights.
## The neighborhood payload is intentionally compact: only indices and small
## per-task metadata are retained; coordinates/design rows are reconstructed
## from the global objects once inside the worker.
## -------------------------------------------------------------------------
.GeoKriglocWeightsTask <- function(i, vars) {
  idx <- vars$neigh$indices[[i]]
  if (is.null(idx) || !length(idx)) return(NULL)
  idx <- as.integer(idx)

  if (identical(vars$mode, "space") && isTRUE(vars$compact) &&
      !is.null(vars$pair_cache) && !is.null(vars$cov_template)) {
    return(list(
      weights = .GeoKriglocWeightsCachedSpatial(idx, i, vars),
      neighbor_indices = idx
    ))
  }

  param_local <- vars$param
  X_neigh <- NULL
  coordx_arg <- NULL
  coordx_dyn_arg <- NULL

  if (identical(vars$mode, "space")) {
    coordx_arg <- vars$coords_global[idx, , drop = FALSE]
    if (!is.null(vars$X_global))
      X_neigh <- vars$X_global[idx, , drop = FALSE]
    if (!is.null(vars$M_global))
      param_local$mean <- vars$M_global[idx]
  } else {
    idx_by_var <- vars$neigh$indices_by_variable[[i]]
    if (is.null(idx_by_var) || length(idx_by_var) != 2L ||
        any(vapply(idx_by_var, length, integer(1L)) == 0L)) return(NULL)

    if (!is.null(vars$coordx_dyn_global)) {
      coordx_dyn_arg <- Map(
        function(z, ii) z[ii, , drop = FALSE],
        vars$coordx_dyn_global, idx_by_var
      )
    } else {
      coordx_dyn_arg <- lapply(
        idx_by_var,
        function(ii) vars$coords_global[ii, , drop = FALSE]
      )
    }
    if (!is.null(vars$X_global))
      X_neigh <- vars$X_global[idx, , drop = FALSE]
  }

  xloc_row <- if (!is.null(vars$Xloc)) vars$Xloc[i, , drop = FALSE] else NULL

  ans <- GeoKrigWeights(
    coordx = coordx_arg,
    coordx_dyn = coordx_dyn_arg,
    loc = matrix(vars$loc_metric[i, ], nrow = 1L),
    corrmodel = vars$corrmodel,
    distance = vars$distance,
    n = vars$n,
    X = X_neigh,
    Xloc = xloc_row,
    Mloc = if (!is.null(vars$Mloc)) vars$Mloc[i] else NULL,
    sparse = vars$sparse,
    model = vars$model,
    param = param_local,
    anisopars = NULL,
    radius = vars$radius,
    copula = vars$copula,
    method = vars$method,
    which = vars$which
  )

  if (isTRUE(vars$compact)) {
    return(list(
      weights = as.numeric(ans$weights),
      neighbor_indices = idx
    ))
  }
  ans$neighbor_indices <- idx
  ans
}

.GeoKriglocWeightsRunTasks <- function(n_tasks, vars, parallel, workers) {
  ids <- seq_len(n_tasks)
  if (!isTRUE(parallel) || workers <= 1L) {
    return(lapply(ids, .GeoKriglocWeightsTask, vars = vars))
  }
  future.apply::future_lapply(
    ids, .GeoKriglocWeightsTask, vars = vars,
    future.seed = TRUE, future.globals = TRUE
  )
}

GeoKriglocWeights <- function(coordx=NULL, coordy=NULL, coordz=NULL, coordt=NULL, coordx_dyn=NULL,
                             corrmodel, distance="Eucl", grid=FALSE, loc, neighb=NULL,
                             maxdist=NULL, maxtime=NULL, method="cholesky", model="Gaussian",
                             n=1, nloc=NULL, param, anisopars=NULL,
                             radius=1, sparse=FALSE, time=NULL, which=1,
                             copula=NULL, X=NULL, Xloc=NULL, Mloc=NULL,
                             parallel=FALSE, ncores=6, compact=FALSE)  
{
  ###########################################################################
  ##  Validazione input
  ###########################################################################
  call = match.call()
  
  if(is.null(CkModel(model))) stop("The name of the model is not correct\n")
  if(!is.character(corrmodel) || is.null(CkCorrModel(corrmodel))) stop("the name of the correlation model is wrong")
  
  corrmodel <- gsub("[[:blank:]]", "", corrmodel)
  model <- .GeoPredictionCanonicalModel(gsub("[[:blank:]]", "", model))
  distance <- gsub("[[:blank:]]", "", distance)
  method <- .GeoKrigMatchMethod(method, sparse = isTRUE(sparse),
                                context = "GeoKriglocWeights")

  bivariate = CheckBiv(CkCorrModel(corrmodel))
  spacetime = CheckST(CkCorrModel(corrmodel))
  space = !spacetime && !bivariate

  if (!identical(model, "Gaussian"))
    stop("GeoKriglocWeights currently returns validated local weights only for model = 'Gaussian'.", call. = FALSE)
  if (!is.null(copula))
    stop("GeoKriglocWeights currently requires copula = NULL.", call. = FALSE)

  if (spacetime && !is.null(coordx_dyn)) {
    coordx_dyn <- .GeoNormalizeDynamicCoordinates(
      coordx_dyn, coordt = coordt, context = "GeoKriglocWeights"
    )
    .GeoValidateDynamicRows(X, coordx_dyn, name = "X",
                            context = "GeoKriglocWeights")
  }

  biv_info <- NULL
  if (bivariate) {
    biv_info <- .GeoBivariateCoordinateInfo(
      coordx = coordx, coordy = coordy, coordz = coordz,
      coordx_dyn = coordx_dyn, grid = grid,
      context = "GeoKriglocWeights"
    )
    if (!is.null(X)) {
      X <- .GeoBivariateDesign(X, biv_info$ns,
                               context = "GeoKriglocWeights")
    }
    if (biv_info$dynamic) {
      coordx_dyn <- biv_info$blocks
      coordx <- coordy <- coordz <- NULL
    }
  }

  # coordinate osservate
  if(is.null(coordx_dyn)){
    if(is.null(coordy)) {
      coords = as.matrix(coordx)
    } else {
      if(!grid) {
        coords = cbind(coordx, coordy)
        if(!is.null(coordz)) coords = cbind(coords, coordz)
      } else {
        if(!is.null(coordz)) {
          coords = as.matrix(expand.grid(coordx, coordy, coordz))
        } else {
          coords = as.matrix(expand.grid(coordx, coordy))
        }
      }
    }
  } else {
    coords = NULL
  }

  if(is.vector(loc)) loc = matrix(loc, nrow=1)
  loc = as.matrix(loc)

  .st_old <- .GeoSTEnter(corrmodel,
    coordx = if (is.null(coordx_dyn)) coords else NULL,
    coordx_dyn = coordx_dyn, distance = distance, context = "GeoKriglocWeights")
  on.exit(.GeoSTRestore(.st_old), add = TRUE)
  .GeoSTValidate(corrmodel, param, context = "GeoKriglocWeights")

  .GeoValidate3D(
    coordx = if (is.null(coordx_dyn)) coords else NULL,
    coordx_dyn = coordx_dyn, distance = distance, grid = grid,
    context = "GeoKriglocWeights"
  )
  .GeoValidate3D(
    coordx = loc, distance = distance, grid = FALSE,
    context = "GeoKriglocWeights prediction locations"
  )
  .GeoValidateCoordinateCompatibility(
    coordx = if (is.null(coordx_dyn)) coords else NULL,
    coordx_dyn = coordx_dyn, loc = loc, context = "GeoKriglocWeights"
  )
  Nloc = nrow(loc)
  
  if(is.null(time)) time = 0
  Tloc = length(time)
  if(bivariate) Tloc = 1
  
  ## Mean specification for local weights.  The weights themselves do not
  ## depend on the mean in the Gaussian case, but X/Xloc and M/Mloc must be
  ## aligned with the same local tasks used by GeoKrigloc().
  if (!bivariate) {
    nobs_mean <- if (!is.null(coordx_dyn)) {
      sum(vapply(coordx_dyn, nrow, integer(1L)))
    } else if (spacetime) {
      nrow(coords) * length(coordt)
    } else {
      nrow(coords)
    }
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
      param, p_mean, context = "GeoKriglocWeights"
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

  ###########################################################################
  ##  SUGG. 1: anisotropia anche per la selezione dei vicini
  ###########################################################################
  coords_aniso <- coords
  coordx_dyn_aniso <- coordx_dyn
  loc_aniso    <- loc
  if (!is.null(anisopars)) {
    if (!is.null(coords)) coords_aniso <- GeoAniso(coords, c(anisopars$angle, anisopars$ratio))
    if (!is.null(coordx_dyn)) {
      if (ncol(coordx_dyn[[1L]]) != 2L)
        stop("GeoKriglocWeights: anisotropy for dynamic three-dimensional coordinates is not implemented.", call. = FALSE)
      coordx_dyn_aniso <- lapply(coordx_dyn, function(z)
        GeoAniso(z, c(anisopars$angle, anisopars$ratio)))
    }
    loc_aniso <- GeoAniso(loc, c(anisopars$angle, anisopars$ratio))
  }

  ###########################################################################
  ##  Memory + centralized parallel plan
  ###########################################################################
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel))
    stop("parallel must be TRUE or FALSE.", call. = FALSE)
  if (!is.logical(compact) || length(compact) != 1L || is.na(compact))
    stop("compact must be TRUE or FALSE.", call. = FALSE)

  total_parallel_tasks <- if (spacetime) Nloc * Tloc else Nloc
  workers <- .GeoResolveWorkers(parallel, ncores, n_jobs = total_parallel_tasks)
  parallel <- isTRUE(parallel) && workers > 1L

  if (parallel && !spacetime) {
    old_limit <- getOption("future.globals.maxSize")
    old_plan <- future::plan()
    options(future.globals.maxSize = 4 * 1024^3)
    .GeoKriglocSetFuturePlan(workers)
    on.exit(options(future.globals.maxSize = old_limit), add = TRUE)
    on.exit(try(future::plan(old_plan), silent = TRUE), add = TRUE)
  }

  ###########################################################################
  ##  Loop Kriging - SPACE
  ###########################################################################
  if(space){
    neigh = GeoNeighborhood(coordx = coords_aniso, distance = distance, loc = loc_aniso,
                            neighb = neighb, maxdist = maxdist, X = X, M = M,
                            parallel = FALSE, ncores = NULL)

    compact_neigh <- list(indices = neigh$indices)
    cov_template <- pair_cache <- cache_param_vec <- NULL
    if (isTRUE(compact) && !isTRUE(sparse) && identical(method, "cholesky")) {
      cov_template <- .GeoKriglocPrepareSpatialCovariance(
        neigh = compact_neigh, coordx_global = coords_aniso,
        X_global = X, param = param, corrmodel = corrmodel,
        distance = distance, n = n, model = model, radius = radius,
        copula = NULL, M = M
      )
      if (!is.null(cov_template) &&
          .GeoPairCacheEnabled(corrmodel, param, weights = TRUE, spacetime = FALSE)) {
        cache_param_vec <- .GeoKriglocLocalParamVector(cov_template, param)
        pair_cache <- .GeoBuildSpatialPairCache(
          compact_neigh$indices, coords_aniso,
          corrmodel_code = CkCorrModel(corrmodel),
          corrparam = cache_param_vec[cov_template$namescorr],
          distance_code = CheckDistance(distance),
          radius = cov_template$local_cov_radius,
          context = "GeoKriglocWeights"
        )
      }
    }
    task_vars <- list(
      mode = "space",
      neigh = compact_neigh,
      coords_global = coords_aniso,
      coordx_dyn_global = NULL,
      X_global = X,
      M_global = M,
      loc_metric = loc_aniso,
      Xloc = Xloc,
      Mloc = Mloc,
      param = param,
      corrmodel = corrmodel,
      distance = distance,
      n = n,
      sparse = sparse,
      model = model,
      radius = radius,
      copula = copula,
      method = method,
      which = which,
      compact = compact,
      pair_cache = pair_cache,
      cov_template = cov_template,
      cache_param_vec = cache_param_vec,
      corrmodel_code = CkCorrModel(corrmodel),
      distance_code = CheckDistance(distance)
    )
    rm(neigh, compact_neigh)
    invisible(gc(verbose = FALSE))
    weights_list <- .GeoKriglocWeightsRunTasks(
      Nloc, task_vars, parallel = parallel, workers = workers
    )
  }

  ###########################################################################
  ##  Loop Kriging - SPACETIME
  ###########################################################################
  if(spacetime){
    neigh = GeoNeighborhood(coordx = coords_aniso, coordt = coordt,
                            distance = distance, neighb = neighb,
                            coordx_dyn = coordx_dyn_aniso, loc = loc_aniso, time = time,
                            maxdist = maxdist, maxtime = maxtime, X = X, M = M,
                            parallel = FALSE, ncores = NULL)
    total_tasks = Nloc * Tloc
    weights_list = vector("list", total_tasks)

    st_pair_cache <- st_geometry <- corrparam_st <- NULL
    cache_st <- isTRUE(compact) && !isTRUE(sparse) &&
      identical(tolower(method), "cholesky") &&
      .GeoPairCacheEnabled(corrmodel, param, weights = TRUE, spacetime = TRUE)
    if (cache_st) {
      st_geometry <- .GeoSpacetimeObservationGeometry(
        coords_aniso, coordx_dyn_aniso, coordt, context = "GeoKriglocWeights"
      )
      corrparam_st <- .GeoPairCorrelationParam(
        corrmodel, param, context = "GeoKriglocWeights"
      )
      st_pair_cache <- .GeoBuildSpacetimePairCache(
        neigh$indices, st_geometry$coords, st_geometry$times,
        corrmodel_code = CkCorrModel(corrmodel), corrparam = corrparam_st,
        distance_code = CheckDistance(distance), radius = radius,
        context = "GeoKriglocWeights"
      )
    }
    st_cache_vars <- if (!is.null(st_pair_cache)) list(
      pair_cache = st_pair_cache, st_geometry = st_geometry, param = param,
      corrmodel_code = CkCorrModel(corrmodel), corrparam = corrparam_st,
      distance_code = CheckDistance(distance), radius = radius
    ) else NULL

    k = 1
    for(i in 1:Nloc){
      for(j in 1:Tloc){
        param_local = param
        if(!is.null(M) && !is.null(neigh$M)) param_local$mean = neigh$M[[k]]
        coords_neigh = neigh$coordx[[k]]
        # SUGG. 4: check robusto
        if (!is.null(coords_neigh) && !(is.matrix(coords_neigh) && nrow(coords_neigh) == 0) && length(coords_neigh) > 0){
          # SUGG. 3: Xloc riga robusta con indice spacetime
          idx_st <- k
          xloc_row <- if (!is.null(Xloc)) Xloc[idx_st, , drop = FALSE] else NULL

          if (!is.null(st_cache_vars)) {
            weights_list[[k]] <- list(
              weights = .GeoKriglocWeightsCachedSpacetime(
                as.integer(neigh$indices[[k]]), loc_aniso[i, ], time[j],
                st_cache_vars
              ),
              neighbor_indices = as.integer(neigh$indices[[k]])
            )
            k = k + 1L
            next
          }

          coordx_dyn_neigh = if (!is.null(neigh$coordx_dyn)) neigh$coordx_dyn[[k]] else NULL
          weights_result = GeoKrigWeights(
            coordx = if (is.null(coordx_dyn_neigh)) coords_neigh else NULL,
            coordt = neigh$coordt[[k]],
            coordx_dyn = coordx_dyn_neigh,
            loc = matrix(loc_aniso[i, ], nrow = 1),
            time = time[j], 
            X = if(!is.null(neigh$X)) neigh$X[[k]] else NULL,
            Mloc = if(!is.null(Mloc)) Mloc[idx_st] else NULL,
            Xloc = xloc_row,
            sparse = sparse, 
            corrmodel = corrmodel,
            distance = distance, 
            model = model, 
            param = param_local,
            anisopars = NULL,
            radius = radius,
            copula = copula, 
            n = n, 
            method = method,
            which = which
          )
          if (isTRUE(compact)) {
            weights_list[[k]] <- list(
              weights = as.numeric(weights_result$weights),
              neighbor_indices = as.integer(neigh$indices[[k]])
            )
          } else {
            weights_result$neighbor_indices = neigh$indices[[k]]
            weights_list[[k]] = weights_result
          }
        } else {
          weights_list[k] <- list(NULL)
        }
        k = k + 1
      }
    }
  }

  ###########################################################################
  ##  Loop Kriging - BIVARIATE
  ###########################################################################
  if(bivariate){
    neigh = GeoNeighborhood(coordx = coords_aniso, distance = distance,
                            neighb = neighb, coordx_dyn = coordx_dyn_aniso,
                            bivariate = TRUE,
                            loc = loc_aniso, maxdist = maxdist, X = X, M = M,
                            parallel = FALSE, ncores = NULL)

    compact_neigh <- list(
      indices = neigh$indices,
      indices_by_variable = neigh$indices_by_time
    )
    task_vars <- list(
      mode = "bivariate",
      neigh = compact_neigh,
      coords_global = coords_aniso,
      coordx_dyn_global = coordx_dyn_aniso,
      X_global = X,
      M_global = NULL,
      loc_metric = loc_aniso,
      Xloc = Xloc,
      Mloc = Mloc,
      param = param,
      corrmodel = corrmodel,
      distance = distance,
      n = n,
      sparse = sparse,
      model = model,
      radius = radius,
      copula = copula,
      method = method,
      which = which,
      compact = compact
    )
    rm(neigh, compact_neigh)
    invisible(gc(verbose = FALSE))
    weights_list <- .GeoKriglocWeightsRunTasks(
      Nloc, task_vars, parallel = parallel, workers = workers
    )
  }

  ###########################################################################
  ##  Output
  ###########################################################################
  GeoKriglocWeights_out = list(
    bivariate = bivariate,
    coordx = coordx,
    coordy = coordy,
    coordz = coordz,
    coordt = coordt,
    coordx_dyn = coordx_dyn,
    corrmodel = corrmodel,
    data = NULL,
    distance = distance,
    grid = grid,
    loc = loc,
    copula = copula,
    numcoord = if(!is.null(coords)) nrow(coords) else length(coordx),
    numloc = Nloc,
    numtime = if(!is.null(coordt)) length(coordt) else 0,
    numt = Tloc,
    maxdist = maxdist,
    maxtime = maxtime,
    model = model,
    n = n,
    param = param,
    weights = weights_list,
    radius = radius,
    spacetime = spacetime,
    time = time,
    which = which,
    neighb = neighb,
    compact = compact,
    ncores = if (parallel) workers else 1L
  )
  
  structure(c(GeoKriglocWeights_out, call = call), class = c("GeoKriglocWeights"))
}
