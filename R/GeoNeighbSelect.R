GeoNeighbSelect <- function(
  data, coordx, coordy=NULL, coordz=NULL, coordt=NULL, coordx_dyn=NULL,
  copula=NULL, corrmodel=NULL, distance="Eucl", fixed=NULL, anisopars=NULL,
  est.aniso=c(FALSE,FALSE), grid=FALSE, likelihood='Marginal',
  lower=NULL, neighb=c(1,2,3,4,5), p_neighb=1,
  maxtime=Inf, memdist=TRUE, model='Gaussian', n=1, ncores=6,
  optimizer='Nelder-Mead', parallel=FALSE, bivariate=FALSE, radius=1, start=NULL,
  type='Pairwise', upper=NULL, weighted=FALSE,
  X=NULL, spobj=NULL, spdata=NULL, vario=NULL, progress=TRUE, check.duplicates=FALSE
)
{
  if (!isTRUE(memdist)) {
    warning(
      "'memdist = FALSE' is deprecated and ignored; GeoNeighbSelect now always precomputes the selected pair structure.",
      call. = FALSE
    )
  }
  memdist <- TRUE

  if (!is.numeric(neighb) || length(neighb) == 0L)
    stop("neighb must be a non-empty numeric vector")

  if (anyNA(neighb) || any(!is.finite(neighb)))
    stop("neighb must contain only finite, non-missing values")

  if (any(neighb <= 0))
    stop("neighb must contain only positive values")

  tol <- sqrt(.Machine$double.eps)
  if (any(abs(neighb - round(neighb)) > tol))
    stop("neighb must be a positive integer numeric vector")

  neighb <- as.integer(round(neighb))

  if (!is.numeric(p_neighb) || length(p_neighb) != 1L || !is.finite(p_neighb))
    stop("p_neighb must be a finite numeric scalar.")
  if (p_neighb <= 0 || p_neighb > 1)
    stop("p_neighb must be in (0,1].")
  if (p_neighb < 1) {
    warning(
      "With 'p_neighb < 1', neighborhood selection depends on the stochastic thinning realization; use a fixed seed or p_neighb = 1 for a more stable comparison.",
      call. = FALSE
    )
  }

  tau <- 0.15

  corrmodel_valid <- try(CkCorrModel(corrmodel), silent = TRUE)
  if (inherits(corrmodel_valid, "try-error") || is.null(corrmodel_valid))
    stop("correlation model is not valid\n")

  bivariate <- CheckBiv(corrmodel_valid)
  spacetime <- CheckST(corrmodel_valid)
  space <- !spacetime && !bivariate

  ## Validate the observation design once before fitting all candidate
  ## neighborhoods. Subsequent GeoFit calls reuse subsets/identical locations
  ## and must not rescan a potentially massive coordinate set.
  dup_coordx <- coordx; dup_coordy <- coordy; dup_coordz <- coordz
  dup_coordt <- coordt; dup_coordx_dyn <- coordx_dyn
  if (!is.null(spobj)) {
    aa_dup <- sp2Geo(spobj, spdata)
    dup_coordx <- aa_dup$coords
    if (spacetime) dup_coordt <- aa_dup$coordt
  }
  .GeoCheckDuplicateLocations(
    coordx = dup_coordx, coordy = dup_coordy, coordz = dup_coordz,
    coordt = dup_coordt, coordx_dyn = dup_coordx_dyn,
    spacetime = spacetime, bivariate = bivariate, grid = grid,
    context = "GeoNeighbSelect", check.duplicates = check.duplicates
  )
  .dup_old <- .GeoSuppressDuplicateChecksForNestedCalls()
  on.exit(.GeoRestoreDuplicateCheckOption(.dup_old), add = TRUE)

  if (!is.null(vario)) {
    if (!inherits(vario, "GeoVariogram"))
      stop("A GeoVariogram object is needed as input for vario\n")
    vario_bivariate <- isTRUE(vario$bivariate)
    vario_spacetime <- !is.null(vario$bint)
    ok_type <- (bivariate && vario_bivariate && !vario_spacetime) ||
      (spacetime && vario_spacetime && !vario_bivariate) ||
      (space && !vario_spacetime && !vario_bivariate)
    if (!ok_type)
      stop("The GeoVariogram object is not of the same type of the correlation model\n")
    semiv <- vario
  } else {
    stop("A GeoVariogram object is needed as input for vario\n")
  }

  if (is.null(start) || length(start) == 0L)
    stop("start parameter must be provided\n")
  start_names <- names(start)
  if (is.null(start_names) || any(!nzchar(start_names)) || anyDuplicated(start_names))
    stop("start must be a uniquely named list.\n")

  # Build one empirical layout and use the same ordering for every candidate fit.
  build_variogram_layout <- function() {
    spatial_lags <- as.numeric(semiv$centers)
    if (!length(spatial_lags) || anyNA(spatial_lags) || any(!is.finite(spatial_lags)))
      stop("The GeoVariogram object does not contain valid spatial bin centers.\n")

    if (bivariate) {
      marginal <- semiv$variograms
      cross <- as.numeric(semiv$variogramst)
      if (is.null(dim(marginal)) || nrow(marginal) != 2L ||
          ncol(marginal) != length(spatial_lags) ||
          length(cross) != length(spatial_lags)) {
        stop(
          "The bivariate GeoVariogram object must contain two marginal variograms and one cross-variogram evaluated at the same spatial bin centers.\n"
        )
      }
      empirical <- rbind(
        as.numeric(marginal[1L, ]),
        cross,
        as.numeric(marginal[2L, ])
      )
      return(list(
        kind = "bivariate",
        x = spatial_lags,
        t = NULL,
        empirical = empirical
      ))
    }

    if (spacetime) {
      temporal_lags <- semiv$centert
      if (is.null(temporal_lags) || length(temporal_lags) == 0L)
        temporal_lags <- semiv$bint
      temporal_lags <- as.numeric(temporal_lags)

      nx <- length(spatial_lags)
      nt <- length(temporal_lags)
      spatial_vario <- as.numeric(semiv$variograms)
      temporal_vario <- as.numeric(semiv$variogramt)
      cross_vario <- as.numeric(semiv$variogramst)

      if (!nt || anyNA(temporal_lags) || any(!is.finite(temporal_lags)) ||
          length(spatial_vario) != nx || length(temporal_vario) != nt ||
          length(cross_vario) != nx * nt) {
        stop(
          "The space-time GeoVariogram object must contain compatible spatial and temporal margins and a complete rectangular space-time variogram surface.\n"
        )
      }

      cross_matrix <- matrix(cross_vario, nrow = nx, ncol = nt, byrow = TRUE)
      empirical <- rbind(
        c(0, temporal_vario),
        cbind(spatial_vario, cross_matrix)
      )
      return(list(
        kind = "spacetime",
        x = c(0, spatial_lags),
        t = c(0, temporal_lags),
        empirical = empirical
      ))
    }

    empirical <- as.numeric(semiv$variograms)
    if (length(empirical) != length(spatial_lags))
      stop("The GeoVariogram values and spatial bin centers have incompatible lengths.\n")

    list(
      kind = "spatial",
      x = spatial_lags,
      t = NULL,
      empirical = empirical
    )
  }

  variogram_layout <- build_variogram_layout()

  format_theoretical_variogram <- function(cc) {
    if (variogram_layout$kind == "bivariate") {
      cc <- as.matrix(cc)
      if (!identical(dim(cc), dim(variogram_layout$empirical)))
        return(NULL)
      return(cc)
    }

    if (variogram_layout$kind == "spacetime") {
      nr <- length(variogram_layout$x)
      nc <- length(variogram_layout$t)
      if (length(cc) != nr * nc)
        return(NULL)
      return(matrix(as.numeric(cc), nrow = nr, ncol = nc))
    }

    cc <- as.numeric(cc)
    if (length(cc) != length(variogram_layout$empirical))
      return(NULL)
    cc
  }

  variogram_sse <- function(theoretical) {
    empirical <- variogram_layout$empirical
    if (is.null(theoretical) || !identical(dim(empirical), dim(theoretical))) {
      return(Inf)
    }

    empirical <- as.numeric(empirical)
    theoretical <- as.numeric(theoretical)
    empirical_ok <- is.finite(empirical)
    if (!any(empirical_ok) || any(!is.finite(theoretical[empirical_ok])))
      return(Inf)

    sum((empirical[empirical_ok] - theoretical[empirical_ok])^2)
  }

  align_estimates <- function(param) {
    ans <- rep(NA_real_, length(start_names))
    names(ans) <- start_names
    values <- unlist(param, use.names = TRUE)

    if (!length(values))
      return(ans)

    if (!is.null(names(values)) && any(nzchar(names(values)))) {
      common <- intersect(start_names, names(values))
      ans[common] <- suppressWarnings(as.numeric(values[common]))
    } else if (length(values) == length(ans)) {
      ans[] <- suppressWarnings(as.numeric(values))
    }
    ans
  }

  if (spacetime) {
    if (!is.numeric(maxtime) || length(maxtime) == 0L)
      stop("maxtime must be a non-empty numeric vector")
    candidate_grid <- expand.grid(
      neighb = neighb, maxtime = maxtime,
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
  } else {
    candidate_grid <- data.frame(
      neighb = neighb, maxtime = rep(Inf, length(neighb))
    )
  }

  total_iter <- nrow(candidate_grid)
  estimates <- matrix(
    NA_real_, nrow = total_iter, ncol = length(start_names),
    dimnames = list(NULL, start_names)
  )
  res <- rep(NA_real_, total_iter)

  common_params <- list(
    data = data, coordx = coordx, coordy = coordy, coordz = coordz, coordt = coordt,
    coordx_dyn = coordx_dyn, copula = copula, corrmodel = corrmodel, distance = distance,
    fixed = fixed, anisopars = anisopars, est.aniso = est.aniso, grid = grid,
    likelihood = likelihood, lower = lower, memdist = memdist, model = model, n = n,
    optimizer = optimizer, radius = radius, start = start, type = type, upper = upper,
    weighted = weighted, X = X, spobj = spobj, spdata = spdata,
    p_neighb = p_neighb
  )

  compute_fit <- function(neighb_val, maxtime_val = Inf) {
    params <- common_params
    params$neighb <- neighb_val
    params$maxtime <- maxtime_val

    aa <- try(suppressWarnings(do.call(GeoFit, params)), silent = TRUE)
    if (inherits(aa, "try-error") || is.null(aa)) {
      return(list(estimates = align_estimates(NULL), res = Inf))
    }

    estimates_val <- align_estimates(aa$param)
    valid_loglik <- is.numeric(aa$logCompLik) && length(aa$logCompLik) == 1L &&
      is.finite(aa$logCompLik)
    if (!valid_loglik)
      return(list(estimates = estimates_val, res = Inf))

    corr_params <- c(aa$param, aa$fixed)
    corr_call <- if (is.null(copula)) {
      try(
        suppressWarnings(
          GeoCorrFct(
            variogram_layout$x, t = variogram_layout$t,
            corrmodel = corrmodel, model = model, distance = distance,
            param = corr_params, radius = radius, n = n,
            covariance = TRUE, variogram = TRUE
          )
        ),
        silent = TRUE
      )
    } else {
      try(
        suppressWarnings(
          GeoCorrFct_Cop(
            variogram_layout$x, t = variogram_layout$t,
            corrmodel = corrmodel, model = model, distance = distance,
            copula = copula, param = corr_params, radius = radius, n = n,
            covariance = TRUE, variogram = TRUE
          )
        ),
        silent = TRUE
      )
    }

    if (inherits(corr_call, "try-error") || is.null(corr_call$corr))
      return(list(estimates = estimates_val, res = Inf))

    theoretical <- format_theoretical_variogram(corr_call$corr)
    list(estimates = estimates_val, res = variogram_sse(theoretical))
  }

  task_label <- function(i) {
    if (spacetime) {
      sprintf(
        "neighb=%g, maxtime=%g",
        candidate_grid$neighb[i], candidate_grid$maxtime[i]
      )
    } else {
      sprintf("neighb=%g", candidate_grid$neighb[i])
    }
  }

  run_candidate <- function(i) {
    out <- compute_fit(candidate_grid$neighb[i], candidate_grid$maxtime[i])
    c(out$estimates, res = out$res)
  }

  if (progress) {
    restore_progress_handlers <- .GeoProgressHandlersPush(
      handler = "txtprogressbar", global = TRUE
    )
    on.exit(restore_progress_handlers(), add = TRUE)
  }

  resolved_workers <- .GeoResolveWorkers(parallel, ncores, n_jobs = total_iter)
  if (resolved_workers <= 1L) {
    if (progress) {
      pb <- progressr::progressor(
        along = seq_len(total_iter), on_exit = FALSE, auto_finish = FALSE
      )
      ## The serial branch owns this progressor directly.  Finish it before
      ## restoring the temporary global handler on every exit path.
      on.exit(pb(type = "finish"), add = TRUE, after = FALSE)
    }

    n_params <- length(start_names)
    for (i in seq_len(total_iter)) {
      if (progress)
        pb(task_label(i))
      out <- run_candidate(i)
      estimates[i, ] <- out[seq_len(n_params)]
      res[i] <- unname(out[n_params + 1L])
    }
  } else {
    n.cores <- resolved_workers

    cat("Performing", total_iter, "estimations using", n.cores, "cores...\n")

    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = n.cores)

    if (progress) {
      results <- progressr::with_progress({
        pb <- progressr::progressor(along = seq_len(total_iter))
        future.apply::future_lapply(
          seq_len(total_iter),
          function(i) {
            out <- run_candidate(i)
            pb(task_label(i))
            out
          },
          future.seed = TRUE, future.stdout = FALSE, future.conditions = "none"
        )
      })
    } else {
      results <- future.apply::future_lapply(
        seq_len(total_iter), run_candidate,
        future.seed = TRUE, future.stdout = FALSE, future.conditions = "none"
      )
    }

    n_params <- length(start_names)
    estimates <- do.call(rbind, lapply(results, function(x) x[seq_len(n_params)]))
    colnames(estimates) <- start_names
    res <- vapply(results, function(x) unname(x[n_params + 1L]), numeric(1))
  }

  finite_ok <- is.finite(res)
  if (!any(finite_ok)) {
    stop(
      "All candidate neighborhood fits failed or produced a non-finite variogram criterion.",
      call. = FALSE
    )
  }

  indexmin <- which.min(res)
  bestK <- candidate_grid$neighb[indexmin]
  best_T <- if (spacetime) candidate_grid$maxtime[indexmin] else NULL
  best_estimates <- estimates[indexmin, ]

  minres <- min(res[finite_ok])
  idx_ok <- which(finite_ok & res <= (1 + tau) * minres)
  sugg_neighb <- min(candidate_grid$neighb[idx_ok])

  if (spacetime) {
    idx_m <- idx_ok[candidate_grid$neighb[idx_ok] == sugg_neighb]
    idx_s <- idx_m[which.min(res[idx_m])]
    sugg_time <- candidate_grid$maxtime[idx_s]
  } else {
    sugg_time <- NULL
  }

  list(
    best_neighb = as.numeric(bestK),
    best_maxtime = best_T,
    best_estimates = best_estimates,
    res = unname(res),
    estimates = estimates,
    sugg_neighb = as.numeric(sugg_neighb),
    sugg_time = sugg_time
  )
}
