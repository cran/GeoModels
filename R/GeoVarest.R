# Estimates J = Var{score(theta_hat)} using GeoFit(score=TRUE),
# then G^{-1} = H^{-1} J H^{-1}.

# Library paths to use for multisession workers.  In addition to the
# current .libPaths(), explicitly prepend the library that contains the
# GeoModels namespace currently in use.  This makes bootstrap parallelism
# robust when GeoModels is installed in a personal/non-standard library and
# also avoids relying on version-specific future inheritance of .libPaths().
.GeoVarest_future_libpaths <- function() {
  libs <- .libPaths()

  pkg_dir <- tryCatch(
    find.package("GeoModels", quiet = TRUE),
    error = function(e) character(0)
  )

  if (length(pkg_dir) == 1L && nzchar(pkg_dir)) {
    libs <- c(dirname(pkg_dir), libs)
  }

  libs <- libs[nzchar(libs)]
  libs <- vapply(
    libs,
    function(z) normalizePath(z, winslash = "/", mustWork = FALSE),
    character(1L),
    USE.NAMES = FALSE
  )
  unique(libs)
}

# Configure the future backend used by GeoVarest bootstrap work.  On Unix-like
# systems a forked multicore backend avoids serializing the fitted model and
# reloading GeoModels in every worker, which is materially faster for large TB
# jobs.  future::supportsMulticore() already disables this route in environments
# where forking is unsafe (for example some IDE configurations).  Windows and
# unsupported environments transparently fall back to multisession with the
# explicit library-path propagation required by GeoModels.
.GeoVarest_set_future_plan <- function(workers) {
  backend <- tolower(as.character(
    getOption("GeoModels.GeoVarest_future_backend", "auto")
  )[1L])

  if (!(backend %in% c("auto", "multicore", "multisession"))) {
    stop(
      "GeoVarest: option 'GeoModels.GeoVarest_future_backend' must be one of 'auto', 'multicore', or 'multisession'.",
      call. = FALSE
    )
  }

  can_multicore <- tryCatch(
    isTRUE(future::supportsMulticore()),
    error = function(e) FALSE
  )

  use_multicore <- identical(backend, "multicore") ||
    (identical(backend, "auto") && can_multicore)

  if (use_multicore && can_multicore) {
    future::plan(future::multicore, workers = workers)
    return("multicore")
  }

  if (identical(backend, "multicore") && !can_multicore) {
    warning(
      "GeoVarest: multicore futures are not supported in this R session; using multisession.",
      call. = FALSE
    )
  }

  future::plan(
    future::multisession,
    workers = workers,
    rscript_libs = .GeoVarest_future_libpaths()
  )
  "multisession"
}


# Return the exact pair graph stored by the original pairwise fit.  Reusing
# it avoids rebuilding nearest-neighbour indices for every bootstrap sample.
# For old deterministic fits that do not contain valid indices, the caller can
# safely fall back to rebuilding the deterministic graph.  Missing indices are
# an error only for stochastic thinning, because rebuilding would change the
# estimator conditionally on the selected graph.
.GeoVarest_selected_pairs <- function(fit) {
  thin_method <- if (is.null(fit$thin_method)) "bernoulli" else fit$thin_method
  thin_method <- tolower(as.character(thin_method))[1L]
  p_neighb <- if (is.null(fit$p_neighb)) 1 else fit$p_neighb
  p_neighb <- as.numeric(p_neighb)[1L]
  stochastic <- identical(thin_method, "match") ||
    (thin_method %in% c("bernoulli", "fixedbudget", "fixed_budget", "targetbalanced") &&
       is.finite(p_neighb) && p_neighb < 1)

  if (!identical(as.character(fit$type)[1L], "Pairwise")) return(NULL)

  if (!is.null(fit$memdist) && !isTRUE(fit$memdist) && stochastic) {
    stop(
      "This GeoFit object was created with the deprecated memdist = FALSE path and does not contain a reusable stochastic-thinning graph. Refit the model with the current GeoFit implementation before calling GeoVarest or GeoVarestbootstrap.",
      call. = FALSE
    )
  }

  colidx <- fit$colidx
  rowidx <- fit$rowidx
  valid <- is.numeric(colidx) && is.numeric(rowidx) &&
    length(colidx) == length(rowidx) && length(colidx) >= 1L &&
    all(is.finite(colidx)) && all(is.finite(rowidx))

  if (!valid) {
    if (stochastic) {
      stop(
        "The original stochastic-thinning fit does not contain valid pair indices.",
        call. = FALSE
      )
    }
    return(NULL)
  }

  if (!is.integer(colidx)) colidx <- as.integer(colidx)
  if (!is.integer(rowidx)) rowidx <- as.integer(rowidx)

  list(colidx = colidx, rowidx = rowidx)
}

.GeoVarest_with_selected_pairs <- function(pair_index, expr, pair_info = NULL,
                                                fit_context = NULL) {
  run <- function() .GeoWithDuplicateCheckSuppressed(force(expr))
  if (is.null(pair_index) && is.null(fit_context)) return(run())

  old_options <- options(
    GeoModels.selected_pair_index = pair_index,
    GeoModels.selected_pair_info = pair_info,
    GeoModels.fit_context_cache = fit_context
  )
  on.exit(options(old_options), add = TRUE)
  run()
}

# Precompute transient pair distances once per worker for the common purely
# spatial univariate pairwise case.  The cache lives only for the worker call;
# it is never stored in the returned GeoFit object.  More complicated layouts
# retain exact pair-index reuse and let StartParam rebuild their distances.
.GeoVarest_prepare_pair_info <- function(pair_index, coords, fit_spec) {
  if (is.null(pair_index) ||
      !identical(as.character(fit_spec$type)[1L], "Pairwise") ||
      isTRUE(fit_spec$bivariate) || isTRUE(fit_spec$spacetime) ||
      (!is.null(fit_spec$coordt) && length(fit_spec$coordt) > 1L) ||
      !is.null(fit_spec$coordx_dyn)) {
    return(NULL)
  }

  coords <- as.matrix(coords)
  if (nrow(coords) < 1L ||
      max(c(pair_index$colidx, pair_index$rowidx)) > nrow(coords)) {
    return(NULL)
  }

  distance_code <- switch(
    if (is.null(fit_spec$distance)) "Eucl" else as.character(fit_spec$distance)[1L],
    eucl = 0L, Eucl = 0L,
    chor = 1L, Chor = 1L,
    geod = 2L, Geod = 2L,
    proj = 3L, Proj = 3L,
    NULL
  )
  if (is.null(distance_code) || identical(distance_code, 3L)) return(NULL)

  pair_fun <- getFromNamespace(".GeoSelectedPairInfo", "GeoModels")
  tryCatch(
    pair_fun(
      pair_index = pair_index,
      coords = coords,
      distance = distance_code,
      radius = if (is.null(fit_spec$radius)) 1 else fit_spec$radius,
      spacetime = FALSE,
      bivariate = FALSE,
      ns = NULL,
      coordt = NULL
    ),
    error = function(e) NULL
  )
}

.GeoVarest_release_fit_context <- function(fit_spec) {
  if (identical(as.character(fit_spec$type)[1L], "Pairwise")) {
    try(
      .C("DeleteGlobalVar2", PACKAGE = "GeoModels", NAOK = TRUE),
      silent = TRUE
    )
  }
  invisible(NULL)
}


.GeoVarest_estimation_model <- function(fit, context = "GeoVarest") {
  stored <- fit$estimation_model
  if (!is.null(stored) && length(stored) == 1L && !is.na(stored) && nzchar(stored)) {
    return(as.character(stored))
  }

  model <- as.character(fit$model)[1L]
  if (!isTRUE(fit$missp)) return(model)

  ## Backward compatibility for GeoFit objects created before
  ## estimation_model was stored explicitly.  BinomialNeg is ambiguous
  ## because both Gaussian- and binary-working misspecifications exist; such
  ## old objects must be refitted so that the inferential criterion is known.
  fallback <- c(
    StudentT = "Gaussian_misp_StudentT",
    Poisson = "Gaussian_misp_Poisson",
    PoissonZIP = "Gaussian_misp_PoissonZIP",
    PoissonGamma = "Gaussian_misp_PoissonGamma",
    Binomial = "Gaussian_misp_Binomial",
    SkewStudentT = "Gaussian_misp_SkewStudentT",
    Tukeygh = "Gaussian_misp_Tukeygh"
  )

  if (identical(model, "BinomialNeg")) {
    stop(
      context, ": this older misspecified BinomialNeg fit does not record the working likelihood. Refit it with the current GeoFit before bootstrap variance estimation.",
      call. = FALSE
    )
  }
  if (model %in% names(fallback)) return(unname(fallback[[model]]))

  stop(
    context, ": cannot recover the misspecified estimation model from this older GeoFit object. Refit it with the current GeoFit implementation.",
    call. = FALSE
  )
}

.GeoVarest_bootstrap_parameters <- function(fit) {
  start <- fit$param
  fixed <- fit$fixed
  anisopars <- fit$anisopars

  if (!is.null(anisopars)) {
    anames <- intersect(c("angle", "ratio"), names(unlist(anisopars, use.names = TRUE)))
    if (length(anames)) {
      if (!is.null(start)) start[intersect(names(start), anames)] <- NULL
      if (!is.null(fixed)) fixed[intersect(names(fixed), anames)] <- NULL
    }
  }

  list(
    start = start,
    fixed = fixed,
    sim_param = append(start, fixed),
    anisopars = anisopars
  )
}

.GeoVarest_validate_controls <- function(K, parallel, ncores, min_success_rate, context) {
  tol <- sqrt(.Machine$double.eps)
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) || K < 2 ||
      abs(K - round(K)) > tol) {
    stop(context, ": K must be an integer >= 2", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop(context, ": parallel must be TRUE or FALSE", call. = FALSE)
  }
  ncores <- .GeoValidateNcores(ncores, paste0(context, ": ncores"))
  if (!is.numeric(min_success_rate) || length(min_success_rate) != 1L ||
      !is.finite(min_success_rate) || min_success_rate <= 0 || min_success_rate > 1) {
    stop(context, ": min_success_rate must be in (0,1]", call. = FALSE)
  }

  list(
    K = as.integer(round(K)),
    ncores = ncores,
    min_success_rate = as.numeric(min_success_rate)
  )
}

.GeoVarest_check_success <- function(n_successful, K, min_success_rate, context) {
  min_success <- max(2L, as.integer(ceiling(min_success_rate * K)))
  if (n_successful < min_success) {
    stop(
      sprintf(
        "%s: only %d of %d bootstrap replications succeeded (%.1f%%); at least %d (%.1f%%) are required.",
        context, n_successful, K, 100 * n_successful / K,
        min_success, 100 * min_success_rate
      ),
      call. = FALSE
    )
  }
  if (n_successful < K) {
    warning(
      sprintf("%s: %d of %d bootstrap replications failed; variance estimates use the %d successful replications.",
              context, K - n_successful, K, n_successful),
      call. = FALSE
    )
  }
  invisible(min_success)
}

.GeoVarest_preflight_simulation <- function(fit, method, model_sim, sim_param,
                                             context = "GeoVarest") {
  if (!is.null(fit$copula)) {
    if (!(method %in% c("cholesky", "TB"))) {
      stop(context, ": copula bootstrap simulation supports method='cholesky' or method='TB'.",
           call. = FALSE)
    }
    ## Copula simulations use their own latent construction. Validate the
    ## marginal at the copula layer and, for TB, validate the latent backend.
    model_check <- .GeoSimulationCanonicalModel(model_sim)
    tryCatch({
      .GeoValidateCopulaSimulationModel(model_check)
      .GeoValidateCopulaMarginal(model_check, sim_param)
      if (identical(method, "TB")) {
        if (isTRUE(fit$bivariate))
          stop("turning-bands copula simulation is univariate only.", call. = FALSE)
        if (isTRUE(fit$spacetime))
          stop("turning-bands copula simulation is not implemented for spatio-temporal models.", call. = FALSE)
        .GeoValidateTBDistance(
          if (is.null(fit$distance)) "Eucl" else fit$distance,
          context = context
        )
        .GeoValidateTBCorrelationModel(.GeoSimulationCanonicalCorrModel(fit$corrmodel), bivariate = FALSE)
      }
    }, error = function(e) {
      stop(context, ": fitted copula model cannot be simulated: ",
           conditionMessage(e), call. = FALSE)
    })
    return(invisible(TRUE))
  }

  simulator <- if (identical(method, "cholesky")) "GeoSim" else "GeoSimapprox"
  model_check <- .GeoSimulationCanonicalModel(model_sim)
  tryCatch({
    .GeoValidateDirectSimulationModel(model_check, simulator = simulator)
    .GeoValidateLatentIntegerParameters(model_check, sim_param)
    .GeoValidateSimulationMarginal(model_check, sim_param)
  }, error = function(e) {
    stop(context, ": fitted model cannot be simulated by ", simulator, ": ", conditionMessage(e), call. = FALSE)
  })
  invisible(TRUE)
}

.GeoVarest_reject_full <- function(fit, context) {
  if (identical(as.character(fit$likelihood)[1L], "Full") &&
      identical(as.character(fit$type)[1L], "Standard")) {
    stop(
      context,
      ": Full/Standard likelihood fits do not use bootstrap Godambe variance estimation. Refit with GeoFit(..., likelihood='Full', type='Standard', varest=TRUE) to obtain Hessian-based standard errors.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

## Evaluate an expression under a deterministic scalar seed without changing
## the RNG state seen by the caller.  GeoVarest uses this for per-replication
## bootstrap streams so results do not depend on worker scheduling.
.GeoVarest_with_rng_kind <- function(rng_kind, expr) {
  if (is.null(rng_kind) || !length(rng_kind)) return(force(expr))

  rng_kind <- as.character(rng_kind)
  if (length(rng_kind) < 2L || length(rng_kind) > 3L || any(!nzchar(rng_kind))) {
    stop("GeoVarest: invalid internal RNG kind specification.", call. = FALSE)
  }

  old_kind <- RNGkind()
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL

  restore <- function() {
    do.call(RNGkind, as.list(old_kind))
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }
  on.exit(restore(), add = TRUE)

  do.call(RNGkind, as.list(rng_kind))
  force(expr)
}

.GeoVarest_seeded_eval <- function(seed_value, expr, rng_kind = NULL) {
  if (is.null(seed_value) || length(seed_value) != 1L || !is.finite(seed_value)) {
    return(.GeoVarest_with_rng_kind(rng_kind, force(expr)))
  }

  seed_value <- as.integer(seed_value)

  .GeoVarest_with_rng_kind(rng_kind, {
    set.seed(seed_value)
    force(expr)
  })
}

## Simulate all TB replications assigned to one worker in one GeoSimapprox
## call.  GeoSimapprox receives the original per-replication seeds through an
## internal option and reseeds immediately around each do_one_sim() call.
## Therefore batching removes repeated setup/validation overhead without
## changing the simulated fields.  The hidden option below is useful for
## regression tests and emergency fallback; it is not part of the public API.
.GeoVarest_tb_simulate_batch <- function(ids, sim_spec, sim_seeds, L, rng_kind = NULL) {
  ids <- as.integer(ids)
  if (!length(ids)) return(list())

  GeoSimapprox_fun <- getFromNamespace("GeoSimapprox", "GeoModels")
  use_batch <- !identical(
    getOption("GeoModels.GeoVarest_tb_batch", FALSE),
    FALSE
  )

  make_args <- function(nrep) {
    args <- sim_spec
    args$nrep <- as.integer(nrep)
    args$progress <- FALSE
    args$method <- "TB"
    args$L <- L
    ## GeoVarest already parallelizes across persistent worker chunks.
    ## Avoid nested futures inside each batch.
    args$parallel <- FALSE
    args$ncores <- NULL
    args
  }

  normalize_data <- function(sim_obj, nrep) {
    if (nrep == 1L) return(list(sim_obj$data))
    if (!is.list(sim_obj$data) || length(sim_obj$data) != nrep) {
      stop("GeoSimapprox returned an unexpected batched data layout.", call. = FALSE)
    }
    sim_obj$data
  }

  if (use_batch) {
    old_opt <- options(
      GeoModels.GeoSimapprox_rep_seeds = as.integer(sim_seeds[ids])
    )
    on.exit(options(old_opt), add = TRUE)

    batch_try <- tryCatch({
      sim_obj <- .GeoVarest_with_rng_kind(
        rng_kind,
        do.call(
          GeoSimapprox_fun,
          make_args(length(ids))
        )
      )
      normalize_data(sim_obj, length(ids))
    }, error = function(e) e)

    ## Restore immediately before a possible legacy fallback; the on.exit
    ## above remains as protection against interrupts.
    options(old_opt)

    if (!inherits(batch_try, "error")) return(batch_try)
    ## Preserve the old per-replication failure behavior if an unexpected
    ## model/layout cannot use the batch path.
  }

  lapply(ids, function(k) {
    .GeoVarest_seeded_eval(
      sim_seeds[[k]],
      do.call(GeoSimapprox_fun, make_args(1L)),
      rng_kind = rng_kind
    )$data
  })
}

# Simulate exactly one TB replication under its deterministic replication seed.
# This is the default large-job path used by GeoVarest: a worker simulates one
# field, scores it immediately, releases it, and only then advances to the next
# replication.  Peak memory is therefore O(n) per worker rather than O(K*n),
# which is important when L is large (e.g. L >= 10000).
.GeoVarest_tb_simulate_one <- function(k, sim_spec, sim_seeds, L, rng_kind = NULL) {
  k <- as.integer(k)[1L]
  GeoSimapprox_fun <- getFromNamespace("GeoSimapprox", "GeoModels")

  args <- sim_spec
  args$nrep <- 1L
  args$progress <- FALSE
  args$method <- "TB"
  args$L <- L
  ## Outer GeoVarest parallelism already occupies the available workers.
  ## Nested TB futures would oversubscribe the machine.
  args$parallel <- FALSE
  args$ncores <- NULL

  .GeoVarest_seeded_eval(
    sim_seeds[[k]],
    do.call(GeoSimapprox_fun, args),
    rng_kind = rng_kind
  )$data
}

# Internal worker for TB-based GeoVarest calculations.
# Workers are persistent, but bootstrap fields are streamed one at a time:
# simulate -> score -> release.  This keeps memory bounded for large L/K while
# still reusing the selected pair graph, pair distances and GeoFit static setup
# within each worker.
.GeoVarest_score_worker <- function(ids, mode = c("tb", "file"),
                                      sim_spec, fit_spec, theta_hat,
                                      data_files = NULL, sim_seeds = NULL,
                                      L = NULL, rng_kind = NULL,
                                      progressor = NULL) {
  .dup_old <- .GeoSuppressDuplicateChecksForNestedCalls()
  on.exit(.GeoRestoreDuplicateCheckOption(.dup_old), add = TRUE)

  mode <- match.arg(mode)
  ids <- as.integer(ids)
  target_names <- names(theta_hat)
  ctx <- .GeoVarest_worker_context(fit_spec, sim_spec$coordx)
  on.exit(.GeoVarest_release_fit_context(fit_spec), add = TRUE)

  use_stream <- identical(mode, "tb") && !identical(
    getOption("GeoModels.GeoVarest_tb_stream", TRUE), FALSE
  )
  sim_data <- NULL
  if (identical(mode, "tb") && !use_stream) {
    sim_data <- .GeoVarest_tb_simulate_batch(
      ids, sim_spec, sim_seeds, L, rng_kind = rng_kind
    )
  }

  start <- .GeoVarest_start_from_theta(theta_hat, fit_spec$aniso_names)
  out <- vector("list", length(ids))
  names(out) <- as.character(ids)

  for (ii in seq_along(ids)) {
    k <- ids[[ii]]
    current_file <- if (identical(mode, "file")) data_files[[k]] else NULL

    out[[ii]] <- tryCatch({
      current_data <- if (identical(mode, "file")) {
        readRDS(current_file)
      } else if (use_stream) {
        .GeoVarest_tb_simulate_one(k, sim_spec, sim_seeds, L, rng_kind = rng_kind)
      } else {
        sim_data[[ii]]
      }

      ff <- .GeoVarest_run_fit(
        current_data, sim_spec$coordx, fit_spec, start,
        ctx$pair_info, ctx$fit_context, score = TRUE
      )
      .GeoVarest_extract_score(ff, target_names)
    }, error = function(e) {
      .GeoVarest_failed_result(target_names, e)
    }, finally = {
      if (!is.null(current_file)) try(unlink(current_file), silent = TRUE)
    })

    gc_every <- if (identical(mode, "tb")) 25L else 5L
    if (ii %% gc_every == 0L) gc(verbose = FALSE)
    if (!is.null(progressor)) {
      progressor(sprintf(if (identical(mode, "tb")) "replication %d" else "score %d", k))
    }
  }

  if (!is.null(sim_data)) rm(sim_data)
  gc(verbose = FALSE, full = identical(mode, "file"))
  out
}



GeoVarest <- function(fit, K = 100, sparse = FALSE,
                      method = c("cholesky", "TB", "CE"),
                      alpha = 0.95, L = 10000,
                      parallel = FALSE, ncores = 6, progress = TRUE,
                      seed = NULL, min_success_rate = 0.8) {
  .dup_old <- .GeoSuppressDuplicateChecksForNestedCalls()
  on.exit(.GeoRestoreDuplicateCheckOption(.dup_old), add = TRUE)


  `%||%` <- function(a, b) if (!is.null(a)) a else b

  method <- match.arg(method)

  .GeoVarest_reject_full(fit, "GeoVarest")
  ctl <- .GeoVarest_validate_controls(
    K, parallel, ncores, min_success_rate, "GeoVarest"
  )
  K <- ctl$K
  ncores <- ctl$ncores
  min_success_rate <- ctl$min_success_rate

  .GeoValidate3D(
    coordx = fit$coordx, coordy = fit$coordy, coordz = fit$coordz,
    coordx_dyn = fit$coordx_dyn,
    distance = if (is.null(fit$distance)) "Eucl" else fit$distance,
    grid = isTRUE(fit$grid),
    approximate = !identical(method, "cholesky"), method = method,
    context = "GeoVarest"
  )

  fit_dimension <- .GeoCoordinateDimension(
    fit$coordx, fit$coordy, fit$coordz, fit$coordx_dyn
  )
  if (identical(fit_dimension, 3L) && identical(method, "TB")) {
    if (isTRUE(fit$bivariate) || isTRUE(fit$spacetime)) {
      stop("GeoVarest: three-dimensional TB simulation is supported only for purely spatial univariate fits.",
           call. = FALSE)
    }
    if (!is.null(fit$anisopars)) {
      stop("GeoVarest: three-dimensional TB simulation does not support anisopars; transform the coordinates explicitly.",
           call. = FALSE)
    }
  }

  ## ---------------------------------------------------------------
  ## 0. Checks
  ## ---------------------------------------------------------------
  if (!is.logical(progress) || length(progress) != 1L) {
    stop("progress must be logical", call. = FALSE)
  }

  if (!(is.numeric(alpha) && length(alpha) == 1L && alpha > 0 && alpha < 1)) {
    stop("alpha must be a single numeric in (0,1)", call. = FALSE)
  }


  if (is.null(fit$sensmat)) {
    stop("Sensitivity matrix is missing: use sensitivity = TRUE in GeoFit", call. = FALSE)
  }

  restore_rng <- .GeoVarest_preserve_seed(seed)
  on.exit(restore_rng(), add = TRUE)

  ## ---------------------------------------------------------------
  ## Progress handlers
  ## ---------------------------------------------------------------
  use_progressr <- isTRUE(progress) && requireNamespace("progressr", quietly = TRUE)

  if (isTRUE(progress) && !use_progressr) {
    warning(
      "progress=TRUE but progressr is not available; using text messages only.",
      call. = FALSE
    )
  }

  if (use_progressr) {
    restore_progress_handlers <- .GeoProgressHandlersPush(
      handler = "txtprogressbar", global = TRUE
    )
    on.exit(restore_progress_handlers(), add = TRUE)
  }

  gm_fun <- function(name) {
    if (requireNamespace("GeoModels", quietly = TRUE) &&
        exists(name, envir = asNamespace("GeoModels"), inherits = FALSE)) {
      return(getFromNamespace(name, "GeoModels"))
    }
    if (exists(name, mode = "function", inherits = TRUE)) {
      return(get(name, mode = "function", inherits = TRUE))
    }
    stop("Function ", name, " is not available. Load GeoModels first.", call. = FALSE)
  }


  set_proj_env_for_workers <- function() {
    get_proj_path <- function() {
      path <- Sys.getenv("PROJ_DATA", unset = "")

      if (!nzchar(path)) {
        path <- Sys.getenv("PROJ_LIB", unset = "")
      }

     # if (!nzchar(path) && requireNamespace("sf", quietly = TRUE)) {
     #   path <- tryCatch(sf::sf_proj_info("path"), error = function(e) "")
     # }

      if (length(path) != 1L || !nzchar(path)) {
        return("")
      }

      normalizePath(path, winslash = "/", mustWork = FALSE)
    }

    proj_path <- get_proj_path()

    if (!nzchar(proj_path)) {
      return(character(0))
    }

    Sys.setenv(PROJ_DATA = proj_path)
    Sys.setenv(PROJ_LIB  = proj_path)

    c(PROJ_DATA = proj_path, PROJ_LIB = proj_path)
  }

  GeoFit_fun <- gm_fun("GeoFit")
  GeoSim_fun <- gm_fun("GeoSim")

  ## Use an internal representation for purely spatial models without
  ## modifying the GeoFit object returned to the user. GeoFit may store
  ## absence of time as the scalar 0, whereas GeoSim/GeoFit internal calls
  ## should receive NULL.
  coordt_use <- fit$coordt
  if (is.null(coordt_use) || length(coordt_use) <= 1L) {
    coordt_use <- NULL
  }

  fit$thin_method <- fit$thin_method %||% "bernoulli"
  fit$p_neighb    <- fit$p_neighb %||% 1
  selected_pair_index <- .GeoVarest_selected_pairs(fit)

  optimizer <- fit$optimizer %||% "Nelder-Mead"
  lower     <- fit$lower
  upper     <- fit$upper

  theta_hat <- as.numeric(unlist(fit$param, use.names = TRUE))
  names(theta_hat) <- names(unlist(fit$param, use.names = TRUE))

  num_params <- length(theta_hat)
  if (num_params < 1L) {
    stop("fit$param is empty", call. = FALSE)
  }

  H <- as.matrix(fit$sensmat)
  storage.mode(H) <- "double"

  if (!all(dim(H) == c(num_params, num_params))) {
    stop("fit$sensmat has incompatible dimension", call. = FALSE)
  }

  H <- (H + t(H)) / 2
  rownames(H) <- colnames(H) <- names(theta_hat)

  safe_inverse <- function(A, label) {
    A <- as.matrix(A)
    storage.mode(A) <- "double"
    A <- (A + t(A)) / 2

    out <- try(solve(A), silent = TRUE)

    if (!inherits(out, "try-error") && all(is.finite(out))) {
      return((out + t(out)) / 2)
    }

    ee <- eigen(A, symmetric = TRUE)
    tol <- max(dim(A)) * max(abs(ee$values)) * .Machine$double.eps
    keep <- abs(ee$values) > tol

    if (!any(keep)) {
      stop("Cannot invert ", label, call. = FALSE)
    }

    out <- ee$vectors[, keep, drop = FALSE] %*%
      diag(1 / ee$values[keep], nrow = sum(keep)) %*%
      t(ee$vectors[, keep, drop = FALSE])

    warning("Used spectral generalized inverse for ", label, call. = FALSE)
    (out + t(out)) / 2
  }

  geometry <- .GeoVarest_prepare_geometry(fit, coordt_use)
  dimat <- geometry$dimat
  X_use <- geometry$X
  coords <- geometry$coords

  selected_pair_info <- .GeoVarest_prepare_pair_info(
    selected_pair_index,
    coords,
    list(
      type = fit$type,
      bivariate = fit$bivariate,
      spacetime = fit$spacetime,
      coordx_dyn = fit$coordx_dyn,
      distance = fit$distance,
      radius = fit$radius
    )
  )
  fit_context_main <- if (as.character(fit$type)[1L] %in% c("Pairwise", "Independence")) new.env(parent = emptyenv()) else NULL
  on.exit(.GeoVarest_release_fit_context(list(type = fit$type)), add = TRUE)

  ## Data-generating model and inferential criterion are distinct for
  ## misspecified fits.  GeoFit now records both explicitly.
  model_sim <- as.character(fit$model)[1L]
  model_est <- .GeoVarest_estimation_model(fit, "GeoVarest")
  boot_par <- .GeoVarest_bootstrap_parameters(fit)

  .GeoVarest_preflight_simulation(
    fit, method, model_sim, boot_par$sim_param, context = "GeoVarest"
  )

  estimated_size_mb <- (K * dimat * 8) / (1024^2)

  ## For TB simulations, never send GeoSimapprox's local do_one_sim(),
  ## RFfct1(), and simu_approx() closures to multisession workers.  Those
  ## closures retain the complete simulation frame and can become hundreds
  ## of MiB each.  Instead, stream complete simulation/score replications
  ## through a small namespace-level worker.  Streaming is also activated
  ## sequentially when retaining all K simulated fields would exceed 500 MiB.
  tb_stream_mode <- identical(method, "TB") && is.null(fit$copula) &&
    ((isTRUE(parallel) && K > 1L) || estimated_size_mb > 500)

  if (estimated_size_mb > 500) {
    warning(
      sprintf("Estimated simulated dataset size: %.1f MB", estimated_size_mb),
      call. = FALSE
    )
  }

  ## ---------------------------------------------------------------
  ## 1--3. Streaming TB simulation and score evaluation
  ## ---------------------------------------------------------------
  if (tb_stream_mode) {
    sim_spec <- list(
      coordx     = coords,
      coordt     = coordt_use,
      coordx_dyn = fit$coordx_dyn,
      anisopars  = boot_par$anisopars,
      corrmodel  = fit$corrmodel,
      model      = model_sim,
      param      = boot_par$sim_param,
      grid       = fit$grid,
      X          = fit$X,
      n          = fit$n,
      distance   = fit$distance,
      radius     = fit$radius
    )

    fit_spec <- .GeoVarest_make_fit_spec(
      fit, boot_par, model_est, X_use, optimizer, lower, upper,
      selected_pair_index, coordt_use, memdist = TRUE,
      weighted = fit$weighted %||% FALSE
    )

    ## One deterministic simulation seed per replication.  The result is
    ## reproducible for a fixed GeoVarest(seed=...) independently of worker
    ## scheduling and of the number of score workers.
    sim_seeds <- sample.int(.Machine$integer.max, K, replace = FALSE)
    ## future.apply uses L'Ecuyer-CMRG streams in multisession workers when
    ## future.seed=TRUE.  The bootstrap's explicit scalar replicate seeds must
    ## nevertheless be interpreted under the caller's RNG kind; otherwise the
    ## same scalar seed produces different TB fields in sequential and parallel
    ## execution.  Carry the caller RNG kind explicitly to every worker.
    sim_rng_kind <- RNGkind()

    workers <- .GeoResolveWorkers(parallel, ncores, n_jobs = K)

    use_outer_stream <- isTRUE(parallel) && workers > 1L &&
      requireNamespace("future", quietly = TRUE) &&
      requireNamespace("future.apply", quietly = TRUE)

    if (isTRUE(parallel) && !use_outer_stream) {
      warning(
        "Parallel TB evaluation requested but future/future.apply is unavailable; using sequential streaming.",
        call. = FALSE
      )
    }

    worker_fun <- if (requireNamespace("GeoModels", quietly = TRUE) &&
                      exists(".GeoVarest_score_worker",
                             envir = asNamespace("GeoModels"), inherits = FALSE)) {
      getFromNamespace(".GeoVarest_score_worker", "GeoModels")
    } else {
      .GeoVarest_score_worker
    }

    ids <- seq_len(K)

    if (use_outer_stream) {
      old_plan <- future::plan()
      on.exit(try(future::plan(old_plan), silent = TRUE), add = TRUE)

      old_proj_data <- Sys.getenv("PROJ_DATA", unset = NA)
      old_proj_lib  <- Sys.getenv("PROJ_LIB", unset = NA)
      on.exit({
        if (is.na(old_proj_data)) Sys.unsetenv("PROJ_DATA") else
          Sys.setenv(PROJ_DATA = old_proj_data)
        if (is.na(old_proj_lib)) Sys.unsetenv("PROJ_LIB") else
          Sys.setenv(PROJ_LIB = old_proj_lib)
      }, add = TRUE)

      set_proj_env_for_workers()
      future_backend <- .GeoVarest_set_future_plan(workers)

      if (progress) {
        cat("Performing", K,
            "TB simulations and scores on", workers, "cores...\n")
      }

      chunks <- split(ids, rep(seq_len(workers), length.out = K))
      explicit_payload <- sum(vapply(
        list(sim_spec, fit_spec, theta_hat, sim_seeds, sim_rng_kind),
        function(x) as.numeric(utils::object.size(x)),
        numeric(1)
      ))
      old_max_size <- getOption("future.globals.maxSize")
      on.exit(options(future.globals.maxSize = old_max_size), add = TRUE)
      options(
        future.globals.maxSize = max(500 * 1024^2, 1.5 * explicit_payload)
      )

      run_chunks <- function(progressor = NULL) {
        future.apply::future_lapply(
          chunks,
          FUN = worker_fun,
          mode = "tb",
          sim_spec = sim_spec,
          fit_spec = fit_spec,
          theta_hat = theta_hat,
          sim_seeds = sim_seeds,
          L = L,
          rng_kind = sim_rng_kind,
          progressor = progressor,
          future.seed = TRUE,
          future.scheduling = 1,
          future.globals = FALSE,
          future.packages = "GeoModels"
        )
      }

      result_chunks <- if (use_progressr) {
        progressr::with_progress({
          pb <- progressr::progressor(steps = K)
          pb_worker <- if (as.numeric(utils::object.size(pb)) <= 1024^2) pb else NULL
          run_chunks(pb_worker)
        })
      } else {
        run_chunks(NULL)
      }

      score_list <- vector("list", K)
      for (chunk_result in result_chunks) {
        ids_chunk <- as.integer(names(chunk_result))
        score_list[ids_chunk] <- unname(chunk_result)
      }
    } else {
      if (progress) cat("Performing", K, "TB simulations and scores sequentially...\n")
      if (use_progressr) {
        score_list <- progressr::with_progress({
          pb <- progressr::progressor(steps = K)
          unname(worker_fun(
            ids = ids,
            mode = "tb",
            sim_spec = sim_spec,
            fit_spec = fit_spec,
            theta_hat = theta_hat,
            sim_seeds = sim_seeds,
            L = L,
            rng_kind = sim_rng_kind,
            progressor = pb
          ))
        })
      } else {
        score_list <- unname(worker_fun(
          ids = ids,
          mode = "tb",
          sim_spec = sim_spec,
          fit_spec = fit_spec,
          theta_hat = theta_hat,
          sim_seeds = sim_seeds,
          L = L,
          rng_kind = sim_rng_kind,
          progressor = NULL
        ))
      }
    }

  } else {

  ## ---------------------------------------------------------------
  ## 1. Simulate K datasets
  ## ---------------------------------------------------------------
  sim_args <- list(
    coordx     = coords,
    coordt     = coordt_use,
    coordx_dyn = fit$coordx_dyn,
    anisopars  = boot_par$anisopars,
    corrmodel  = fit$corrmodel,
    model      = model_sim,
    param      = boot_par$sim_param,
    grid       = fit$grid,
    X          = fit$X,
    n          = fit$n,
    distance   = fit$distance,
    radius     = fit$radius,
    nrep       = K,
    progress   = progress
  )

  if (is.null(fit$copula)) {
    if (method == "cholesky") {
      data_sim_full <- do.call(
        GeoSim_fun,
        c(sim_args, list(sparse = sparse, method = method))
      )
    } else {
      GeoSimapprox_fun <- gm_fun("GeoSimapprox")

      ## GeoSimapprox may parallelize over nrep = K. This simulation phase
      ## is completed before score evaluation starts, so the score
      ## parallelism below is not nested. We defensively restore the current
      ## future plan before entering the score phase.
      old_plan_before_sim <- NULL
      if (isTRUE(parallel) && requireNamespace("future", quietly = TRUE)) {
        old_plan_before_sim <- future::plan()
      }

      data_sim_full <- do.call(
        GeoSimapprox_fun,
        c(
          sim_args,
          list(
            method   = method,
            L        = L,
            parallel = parallel,
            ncores   = ncores
          )
        )
      )

      if (!is.null(old_plan_before_sim)) {
        try(future::plan(old_plan_before_sim), silent = TRUE)
      }
    }
  } else {
    if (!(method %in% c("cholesky", "TB"))) {
      stop("Copula simulation supports method = 'cholesky' or method = 'TB'.",
           call. = FALSE)
    }

    GeoSimCopula_fun <- gm_fun("GeoSimCopula")
    copula_args <- list(copula = fit$copula, sparse = sparse, method = method)
    if (identical(method, "TB")) {
      copula_args <- c(copula_args, list(L = L, parallel = parallel, ncores = ncores))
    }

    data_sim_full <- do.call(
      GeoSimCopula_fun,
      c(sim_args, copula_args)
    )
  }

  data_sim <- data_sim_full$data
  rm(data_sim_full)
  gc(verbose = FALSE, full = TRUE)

  ## ---------------------------------------------------------------
  ## 2. Evaluate composite likelihood and score through GeoFit
  ## ---------------------------------------------------------------
  make_start <- function(par) {
    par <- as.numeric(par)
    names(par) <- names(theta_hat)
    .GeoVarest_start_from_theta(par, if (is.null(boot_par$anisopars)) character(0) else c("angle", "ratio"))
  }

  geofit_onlyvar <- function(current_data, par, want_score) {
    call_it <- function() {
      GeoFit_fun(
        data        = current_data,
        start       = make_start(par),
        fixed       = boot_par$fixed,
        coordx      = coords,
        coordt      = coordt_use,
        coordx_dyn  = fit$coordx_dyn,
        copula      = fit$copula,
        anisopars   = boot_par$anisopars,
        est.aniso   = fit$est.aniso,
        thin_method = fit$thin_method,
        lower       = lower,
        upper       = upper,
        neighb      = fit$neighb,
        p_neighb    = fit$p_neighb,
        corrmodel   = fit$corrmodel,
        model       = model_est,
        sparse      = FALSE,
        n           = fit$n,
        maxdist     = fit$maxdist,
        maxtime     = fit$maxtime,
        memdist     = TRUE,
        optimizer   = optimizer,
        grid        = fit$grid,
        likelihood  = fit$likelihood,
        type        = fit$type,
        X           = X_use,
        distance    = fit$distance,
        radius      = fit$radius,
        onlyvar     = TRUE,
        score       = want_score,
        sensitivity = FALSE,
        varest      = FALSE,
        weighted    = fit$weighted %||% FALSE
      )
    }

    .GeoVarest_with_selected_pairs(
      selected_pair_index,
      call_it(),
      pair_info = selected_pair_info,
      fit_context = fit_context_main
    )
  }

  geofit_score <- function(current_data, par) {
    capture.output({
      ff <- geofit_onlyvar(current_data, par, TRUE)
    }, file = nullfile())
    ans <- .GeoVarest_extract_score(ff, names(theta_hat))
    list(score = ans[seq_along(theta_hat)], logCompLik = unname(ans[[length(ans)]]))
  }

  make_score <- function(k, current_data) {
    tryCatch({
      out <- geofit_score(current_data, theta_hat)

      c(as.numeric(out$score), logCompLik = as.numeric(out$logCompLik))

    }, error = function(e) {
      z <- rep(NA_real_, num_params + 1L)
      attr(z, "error") <- conditionMessage(e)
      z
    })
  }

  ## ---------------------------------------------------------------
  ## 3. Score loop. Parallelization is only here.
  ## ---------------------------------------------------------------
  if (progress) {
    cat("Computing", K, "scores ...\n")
  }

  use_parallel <- isTRUE(parallel) && K > 1L &&
    requireNamespace("future", quietly = TRUE) &&
    requireNamespace("future.apply", quietly = TRUE)

  if (isTRUE(parallel) && !use_parallel) {
    warning(
      "Parallel evaluation requested but future/future.apply is unavailable; using sequential evaluation.",
      call. = FALSE
    )
  }

  if (use_parallel) {
    ncores <- .GeoResolveWorkers(TRUE, ncores, n_jobs = K)
    use_parallel <- ncores > 1L
  }

  if (use_parallel) {
    old_plan <- future::plan()
    on.exit(try(future::plan(old_plan), silent = TRUE), add = TRUE)

    ## Workers receive an explicit compact payload below.  Avoid measuring
    ## nested closures here: their captured environments can be very large and
    ## are not exported by the current implementation.
    old_max_size <- getOption("future.globals.maxSize")
    on.exit(options(future.globals.maxSize = old_max_size), add = TRUE)
    safety_margin <- 1.5
    default_cap <- 500 * 1024^2

    old_proj_data <- Sys.getenv("PROJ_DATA", unset = NA)
    old_proj_lib  <- Sys.getenv("PROJ_LIB",  unset = NA)

    on.exit({
      if (is.na(old_proj_data)) {
        Sys.unsetenv("PROJ_DATA")
      } else {
        Sys.setenv(PROJ_DATA = old_proj_data)
      }
      if (is.na(old_proj_lib)) {
        Sys.unsetenv("PROJ_LIB")
      } else {
        Sys.setenv(PROJ_LIB = old_proj_lib)
      }
    }, add = TRUE)

    rscript_envs <- set_proj_env_for_workers()

    future::plan(
      future::multisession,
      workers = ncores,
      rscript_libs = .GeoVarest_future_libpaths()
    )

    temp_dir <- tempdir()
    tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

    data_files <- vapply(seq_len(K), function(k) {
      fp <- file.path(temp_dir, sprintf("scoreJ_%s_%04d.rds", tag, k))
      saveRDS(data_sim[[k]], fp, compress = FALSE)
      fp
    }, character(1))

    on.exit(unlink(data_files), add = TRUE)

    rm(data_sim)
    gc(verbose = FALSE, full = TRUE)

    ## Use the package-level worker with an explicit payload.  This avoids
    ## exporting the nested GeoFit closures (and their captured fit object) to
    ## every multisession worker.  It is especially important for large data
    ## and also keeps deterministic/Independence fits on their original path.
    score_worker <- if (requireNamespace("GeoModels", quietly = TRUE) &&
                        exists(".GeoVarest_score_worker",
                               envir = asNamespace("GeoModels"), inherits = FALSE)) {
      getFromNamespace(".GeoVarest_score_worker", "GeoModels")
    } else {
      .GeoVarest_score_worker
    }

    score_sim_spec <- list(
      coordx = coords,
      coordt = coordt_use,
      coordx_dyn = fit$coordx_dyn
    )

    score_fit_spec <- .GeoVarest_make_fit_spec(
      fit, boot_par, model_est, X_use, optimizer, lower, upper,
      selected_pair_index, coordt_use, memdist = TRUE,
      weighted = fit$weighted %||% FALSE
    )

    chunks <- split(seq_len(K), rep(seq_len(ncores), length.out = K))

    explicit_payload <- sum(vapply(
      list(data_files, score_sim_spec, score_fit_spec, theta_hat),
      function(x) as.numeric(utils::object.size(x)),
      numeric(1)
    ))
    options(
      future.globals.maxSize = max(default_cap, explicit_payload * safety_margin)
    )

    run_score_chunks <- function(progressor = NULL) {
      future.apply::future_lapply(
        chunks,
        FUN = score_worker,
        mode = "file",
        data_files = data_files,
        sim_spec = score_sim_spec,
        fit_spec = score_fit_spec,
        theta_hat = theta_hat,
        progressor = progressor,
        future.seed = TRUE,
        future.globals = FALSE,
        future.packages = "GeoModels"
      )
    }

    result_chunks <- if (use_progressr) {
      progressr::with_progress({
        pb <- progressr::progressor(steps = K)
        pb_worker <- if (as.numeric(utils::object.size(pb)) <= 1024^2) pb else NULL
        run_score_chunks(pb_worker)
      })
    } else {
      run_score_chunks(NULL)
    }

    score_list <- vector("list", K)
    for (chunk_result in result_chunks) {
      ids_chunk <- as.integer(names(chunk_result))
      score_list[ids_chunk] <- unname(chunk_result)
    }

  } else {
    if (use_progressr) {
      score_list <- progressr::with_progress({
        pb <- progressr::progressor(along = seq_len(K))
        out <- vector("list", K)

        for (k in seq_len(K)) {
          out[[k]] <- make_score(k, data_sim[[k]])
          if (k %% 5L == 0L) gc(verbose = FALSE)
          pb(sprintf("score %d", k))
        }

        out
      })
    } else {
      score_list <- vector("list", K)

      for (k in seq_len(K)) {
        if (progress && (k == 1L || k %% 10L == 0L || k == K)) {
          cat("score", k, "of", K, "\n")
        }

        score_list[[k]] <- make_score(k, data_sim[[k]])
        if (k %% 5L == 0L) gc(verbose = FALSE)
      }
    }

    rm(data_sim)
    gc(verbose = FALSE, full = TRUE)
  }

  } # end non-streaming branch

  errors <- vapply(
    score_list,
    function(x) attr(x, "error") %||% "",
    character(1)
  )

  score_raw <- do.call(rbind, score_list)
  colnames(score_raw) <- c(names(theta_hat), "logCompLik")

  valid <- apply(
    score_raw[, seq_len(num_params), drop = FALSE],
    1L,
    function(z) all(is.finite(z))
  ) & is.finite(score_raw[, num_params + 1L])

  score_out <- score_raw[valid, , drop = FALSE]
  n_successful <- nrow(score_out)

  if (n_successful < max(2L, as.integer(ceiling(min_success_rate * K)))) {
    shown <- utils::head(errors[nzchar(errors)], 8L)
    msg <- sprintf(
      "GeoVarest: only %d of %d score-bootstrap replications succeeded (%.1f%%); at least %.1f%% are required.",
      n_successful, K, 100 * n_successful / K, 100 * min_success_rate
    )
    if (length(shown) > 0L) {
      msg <- paste0(msg, "\nFirst errors:\n- ", paste(shown, collapse = "\n- "))
    }
    stop(msg, call. = FALSE)
  }
  .GeoVarest_check_success(n_successful, K, min_success_rate, "GeoVarest")

  if (progress) {
    cat("Successful score evaluations:", n_successful, "out of", K, "\n")
  }

  ## ---------------------------------------------------------------
  ## 4. Sandwich/Godambe and criteria
  ## ---------------------------------------------------------------
  score_mat <- score_out[, seq_len(num_params), drop = FALSE]

  J <- stats::var(score_mat)
  J <- (J + t(J)) / 2
  rownames(J) <- colnames(J) <- names(theta_hat)

  ## ---------------------------------------------------------------
  ## Numerically stable sandwich through diagonal rescaling
  ##
  ## Let D = diag{1 / sqrt(diag(H))}. Then
  ##   Hs = D H D,  Js = D J D,
  ## and, in exact arithmetic,
  ##   H^{-1} J H^{-1}
  ##     = D Hs^{-1} Js Hs^{-1} D.
  ##
  ## This is algebraically identical to the unscaled sandwich, but it
  ## avoids unstable inversion when parameters have very different scales.
  ## ---------------------------------------------------------------
  hdiag <- diag(H)

  if (any(!is.finite(hdiag)) || any(hdiag <= 0)) {
    stop(
      "The sensitivity matrix has non-positive or non-finite diagonal elements.",
      call. = FALSE
    )
  }

  scale_vec <- 1 / sqrt(hdiag)
  names(scale_vec) <- names(theta_hat)

  D <- diag(scale_vec, nrow = num_params, ncol = num_params)
  rownames(D) <- colnames(D) <- names(theta_hat)

  H_scaled <- D %*% H %*% D
  H_scaled <- (H_scaled + t(H_scaled)) / 2
  rownames(H_scaled) <- colnames(H_scaled) <- names(theta_hat)

  J_scaled <- D %*% J %*% D
  J_scaled <- (J_scaled + t(J_scaled)) / 2
  rownames(J_scaled) <- colnames(J_scaled) <- names(theta_hat)

  H_scaled_inv <- safe_inverse(H_scaled, "scaled H/sensmat")
  rownames(H_scaled_inv) <- colnames(H_scaled_inv) <- names(theta_hat)

  Ginv_scaled <- H_scaled_inv %*% J_scaled %*% H_scaled_inv
  Ginv_scaled <- (Ginv_scaled + t(Ginv_scaled)) / 2
  rownames(Ginv_scaled) <- colnames(Ginv_scaled) <- names(theta_hat)

  ## Transform the covariance matrix back to the original parameter scale.
  Ginv <- D %*% Ginv_scaled %*% D
  Ginv <- (Ginv + t(Ginv)) / 2
  rownames(Ginv) <- colnames(Ginv) <- names(theta_hat)

  ## H^{-1} on the original parameter scale, retained for compatibility.
  Hinv <- D %*% H_scaled_inv %*% D
  Hinv <- (Hinv + t(Hinv)) / 2
  rownames(Hinv) <- colnames(Hinv) <- names(theta_hat)

  ## Invert the covariance on the scaled parameter system.  Directly
  ## inverting Ginv on the original scale can recreate the same purely
  ## unit-driven ill-conditioning that motivated the scaling of H.
  G_scaled <- safe_inverse(Ginv_scaled, "scaled Ginv/varcov")
  rownames(G_scaled) <- colnames(G_scaled) <- names(theta_hat)

  Dinv <- diag(1 / scale_vec, nrow = num_params, ncol = num_params)
  rownames(Dinv) <- colnames(Dinv) <- names(theta_hat)

  G <- Dinv %*% G_scaled %*% Dinv
  G <- (G + t(G)) / 2
  rownames(G) <- colnames(G) <- names(theta_hat)

  stderr <- sqrt(pmax(0, diag(Ginv)))
  names(stderr) <- names(theta_hat)

  ## Equivalent to tr(H^{-1} J), evaluated on the stable scaled system.
  penalty <- sum(diag(H_scaled_inv %*% J_scaled))
  penalty_alt <- sum(diag(H_scaled %*% Ginv_scaled))

  fit_loglik <- if (!is.null(fit$logCompLik)) fit$logCompLik else fit$logLik
  fit_loglik <- as.numeric(fit_loglik)

  if (length(fit_loglik) != 1L || !is.finite(fit_loglik)) {
    stop("fit does not contain a finite logCompLik/logLik value", call. = FALSE)
  }

  lik <- as.character(fit$likelihood)
  typ <- as.character(fit$type)

  if ((lik == "Marginal" && typ %in% c("Independence", "Pairwise")) ||
      (lik == "Conditional" && typ == "Pairwise")) {

    claic <- -2 * fit_loglik + 2 * penalty
    clbic <- -2 * fit_loglik + log(dimat) * penalty
    fit$varimat <- H %*% Ginv %*% H

  } else {
    claic <- clbic <- NA_real_
  }

  z_alpha <- stats::qnorm(1 - (1 - alpha) / 2)

  fit$stderr <- stderr
  fit$varcov <- Ginv
  fit$godambe <- G
  fit$Jmat <- J
  fit$Hinv <- Hinv

  ## Numerical diagnostics for the rescaled sensitivity matrix.
  fit$H_scaled <- H_scaled
  fit$J_scaled <- J_scaled
  fit$Ginv_scaled <- Ginv_scaled
  fit$godambe_scaled <- G_scaled
  fit$parameter_scaling <- scale_vec
  fit$H_scaled_condition <- kappa(H_scaled)
  fit$H_scaled_eigenvalues <- eigen(
    H_scaled,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  fit$claic <- claic
  fit$clic <- claic
  fit$clbic <- clbic
  fit$clic_penalty <- penalty
  fit$clic_penalty_alt <- penalty_alt

  fit$conf.int <- rbind(
    theta_hat - z_alpha * stderr,
    theta_hat + z_alpha * stderr
  )
  colnames(fit$conf.int) <- names(theta_hat)
  rownames(fit$conf.int) <- c("Lower", "Upper")

  fit$pvalues_godambe <- 2 * stats::pnorm(-abs(theta_hat / stderr))
  names(fit$pvalues_godambe) <- names(theta_hat)

  ## Keep the standard component name for compatibility with GeoModels.
  fit$pvalues <- fit$pvalues_godambe
  fit$pvalues_type <- "Wald p-values based on score-bootstrap Godambe variance"

  fit$scores <- score_mat
  fit$score_logCompLik <- score_out[, "logCompLik"]

  fit$score_failures <- data.frame(
    iteration = which(!valid),
    error = errors[!valid],
    stringsAsFactors = FALSE
  )

  fit$bootstrap_type <- "parametric_score_J_geofit_score"
  fit$bootstrap_K <- K
  fit$bootstrap_successful <- n_successful
  fit$bootstrap_success_rate <- n_successful / K
  fit$bootstrap_min_success_rate <- min_success_rate
  fit$bootstrap_pair_graph <- if (is.null(selected_pair_index)) "deterministic" else "reused_original"
  fit$bootstrap_thin_method <- fit$thin_method
  fit$bootstrap_p_neighb <- fit$p_neighb
  fit$bootstrap_thinning_seed_rule <- if (is.null(selected_pair_index)) {
    "not applicable"
  } else {
    "not used; the original selected-pair graph is reused"
  }

  fit
}