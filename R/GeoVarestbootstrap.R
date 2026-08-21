.GeoVarestbootstrap_coordt_use <- function(coordt, spacetime = NULL,
                                                coordx_dyn = NULL) {
  ## Older GeoFit objects may encode the absence of a temporal dimension as
  ## a scalar (typically zero).  Passing that scalar back to GeoFit makes a
  ## purely spatial refit look space-time and fails CkInput().  Respect an
  ## explicit spacetime flag; for older objects without the flag, infer a
  ## temporal fit only from dynamic coordinates or more than one time value.
  if (!is.null(coordx_dyn)) return(coordt)
  if (!is.null(spacetime)) {
    if (!isTRUE(spacetime)) return(NULL)
    return(coordt)
  }
  if (is.numeric(coordt) && length(coordt) > 1L) coordt else NULL
}

.GeoVarestbootstrap_refit_worker <- function(ids, mode = c("tb", "file"),
                                               fit_spec, coords, target_names,
                                               data_files = NULL, sim_spec = NULL,
                                               sim_seeds = NULL, L = NULL,
                                               rng_kind = NULL, progressor = NULL) {
  .dup_old <- .GeoSuppressDuplicateChecksForNestedCalls()
  on.exit(.GeoRestoreDuplicateCheckOption(.dup_old), add = TRUE)

  mode <- match.arg(mode)
  ids <- as.integer(ids)
  ctx <- .GeoVarest_worker_context(fit_spec, coords)
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

      fit_result <- .GeoVarest_run_fit(
        current_data, coords, fit_spec, fit_spec$param,
        ctx$pair_info, ctx$fit_context, score = FALSE
      )
      .GeoVarest_extract_refit(fit_result, target_names)
    }, error = function(e) {
      .GeoVarest_failed_result(target_names, e)
    }, finally = {
      if (!is.null(current_file)) unlink(current_file, force = TRUE)
    })

    gc_every <- if (identical(mode, "tb")) 5L else 25L
    if (ii %% gc_every == 0L) gc(verbose = FALSE)
    if (!is.null(progressor)) progressor(sprintf("k=%d", k))
  }

  if (!is.null(sim_data)) rm(sim_data)
  gc(verbose = FALSE, full = identical(mode, "file"))
  out
}


GeoVarestbootstrap <- function(fit, K = 100, sparse = FALSE,
                               optimizer = NULL, lower = NULL, upper = NULL,
                               method = "cholesky", alpha = 0.95,
                               L = 10000, parallel = FALSE, ncores = 6,
                               progress = TRUE, seed = NULL,
                               min_success_rate = 0.8) {
  .dup_old <- .GeoSuppressDuplicateChecksForNestedCalls()
  on.exit(.GeoRestoreDuplicateCheckOption(.dup_old), add = TRUE)


  .GeoVarest_reject_full(fit, "GeoVarestbootstrap")
  ctl <- .GeoVarest_validate_controls(
    K, parallel, ncores, min_success_rate, "GeoVarestbootstrap"
  )
  K <- ctl$K
  ncores <- ctl$ncores
  min_success_rate <- ctl$min_success_rate

  ## ------------------------------------------------------------------
  ## Optional RNG seed, restoring user RNG state on exit
  ## ------------------------------------------------------------------
  restore_rng <- .GeoVarest_preserve_seed(seed)
  on.exit(restore_rng(), add = TRUE)

  ## ------------------------------------------------------------------
  ## Setup cleanup handlers
  ## ------------------------------------------------------------------
  future_plan_original <- future::plan()

  cleanup_resources <- function() {
    tryCatch({
      future::plan(future_plan_original)
    }, error = function(e) {
      tryCatch(
        future::plan(future::sequential),
        error = function(e2) invisible(NULL)
      )
    })
    gc(verbose = FALSE, full = TRUE)
  }
  on.exit(cleanup_resources(), add = TRUE)

  ## ------------------------------------------------------------------
  ## Progress handlers
  ## ------------------------------------------------------------------
  if (!is.logical(progress)) stop("progress must be logical (TRUE or FALSE)")
  restore_progress_handlers <- .GeoProgressHandlersPush(
    handler = if (progress) "txtprogressbar" else "void",
    global = if (progress) TRUE else NULL
  )
  on.exit(restore_progress_handlers(), add = TRUE)

  ## ------------------------------------------------------------------
  ## 1. Check input
  ## ------------------------------------------------------------------
  if (!(method %in% c("cholesky", "TB", "CE"))) {
    stop("The method of simulation is not correct")
  }

  coordt_use <- .GeoVarestbootstrap_coordt_use(
    fit$coordt, fit$spacetime, fit$coordx_dyn
  )

  .GeoValidate3D(
    coordx = fit$coordx, coordy = fit$coordy, coordz = fit$coordz,
    coordx_dyn = fit$coordx_dyn,
    distance = if (is.null(fit$distance)) "Eucl" else fit$distance,
    grid = isTRUE(fit$grid),
    approximate = !identical(method, "cholesky"), method = method,
    context = "GeoVarestbootstrap"
  )

  fit_dimension <- .GeoCoordinateDimension(
    fit$coordx, fit$coordy, fit$coordz, fit$coordx_dyn
  )
  if (identical(fit_dimension, 3L) && identical(method, "TB")) {
    if (isTRUE(fit$bivariate) || isTRUE(fit$spacetime)) {
      stop("GeoVarestbootstrap: three-dimensional TB simulation is supported only for purely spatial univariate fits.",
           call. = FALSE)
    }
    if (!is.null(fit$anisopars)) {
      stop("GeoVarestbootstrap: three-dimensional TB simulation does not support anisopars; transform the coordinates explicitly.",
           call. = FALSE)
    }
  }

  if (!(is.numeric(alpha) && length(alpha) == 1 && alpha > 0 && alpha < 1)) {
    stop("alpha must be a single numeric in (0,1)")
  }

  if (is.null(optimizer)) {
    optimizer <- fit$optimizer
    lower     <- fit$lower
    upper     <- fit$upper
  }

if (is.null(fit$thin_method)) {
  fit$thin_method <- "bernoulli"
}

thin_method_canonical <- function(x) {
  x <- tolower(as.character(x))
  if (x %in% c("bernoulli")) {
    return("bernoulli")
  }
  if (x %in% c("fixedbudget", "fixed_budget", "targetbalanced")) {
    return("FixedBudget")
  }
  if (x %in% c("match")) {
    return("match")
  }
  stop("Unsupported thin_method: ", x, call. = FALSE)
}

fit$thin_method <- thin_method_canonical(fit$thin_method)
  if (is.null(fit$p_neighb)) {
    fit$p_neighb <- 1
  }
  fit_memdist <- TRUE
  fit_weighted <- if (is.null(fit$weighted)) FALSE else isTRUE(fit$weighted)
  selected_pair_index <- .GeoVarest_selected_pairs(fit)
  ## Local reference only: it points to the two index vectors already stored
  ## in fit and is removed before returning the updated object.
  fit$pair_index <- selected_pair_index

  model <- fit$model
  cat("Parametric bootstrap can be time consuming ...\n")

  ## ------------------------------------------------------------------
  ## 2. Data-generating versus estimation model
  ## ------------------------------------------------------------------
  model_sim <- as.character(model)[1L]
  model_est <- .GeoVarest_estimation_model(fit, "GeoVarestbootstrap")
  boot_par <- .GeoVarest_bootstrap_parameters(fit)

  .GeoVarest_preflight_simulation(
    fit, method, model_sim, boot_par$sim_param,
    context = "GeoVarestbootstrap"
  )

  ## ------------------------------------------------------------------
  ## 3. Geometry and design matrix
  ## ------------------------------------------------------------------
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

  ## For TB without a copula, stream each simulation directly into its
  ## bootstrap refit.  No temporary field files are created.
  tb_stream_mode <- identical(method, "TB") && is.null(fit$copula)

  ## ------------------------------------------------------------------
  ## 3.1. Memory pre-check
  ## ------------------------------------------------------------------
  estimated_size_mb <- (K * fit$numtime * fit$numcoord * 8) / (1024^2)
  if (!tb_stream_mode && estimated_size_mb > 500) {
    warning(sprintf(
      "Estimated simulated dataset size: %.1f MB (>500MB). Consider reducing K or using sparse methods.",
      estimated_size_mb
    ))
  }

  ## ------------------------------------------------------------------
  ## 4. Simulation
  ## ------------------------------------------------------------------
  if (!tb_stream_mode) {
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
        data_sim_full <- do.call(GeoSim, c(sim_args, list(sparse = sparse, method = method)))
      } else if (method %in% c("TB", "CE")) {
        ## GeoSimapprox may parallelize over nrep = K. This phase is
        ## completed before the bootstrap refits start, so the refit
        ## parallelism below is not nested.
        old_plan_after_sim <- NULL

        if (isTRUE(parallel) && requireNamespace("future", quietly = TRUE)) {
          old_plan_after_sim <- future::plan()
        }

        data_sim_full <- do.call(
          GeoSimapprox,
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

        if (!is.null(old_plan_after_sim)) {
          try(future::plan(old_plan_after_sim), silent = TRUE)
        }
      } else {
        stop("Unsupported method for simulation")
      }
    } else {
      if (method %in% c("cholesky", "TB")) {
        copula_args <- list(copula = fit$copula, sparse = sparse, method = method)
        if (identical(method, "TB")) {
          copula_args <- c(copula_args, list(L = L, parallel = parallel, ncores = ncores))
        }
        data_sim_full <- do.call(
          GeoSimCopula,
          c(sim_args, copula_args)
        )
      } else {
        stop("Copula simulation supports method = 'cholesky' or method = 'TB'.")
      }
    }

    data_sim <- data_sim_full$data
    rm(data_sim_full)
    gc(verbose = FALSE, full = TRUE)

  }

  ## ------------------------------------------------------------------
  ## 5. Setting parallelization
  ## ------------------------------------------------------------------
  ncores <- .GeoResolveWorkers(parallel, ncores, n_jobs = K)
  parallel <- isTRUE(parallel) && ncores > 1L

  ## ------------------------------------------------------------------
  ## 6. Parameter names and sequential refit helper
  ## ------------------------------------------------------------------
  target_names <- names(fit$param)
  if (is.null(target_names) || length(target_names) != length(fit$param)) {
    target_names <- names(unlist(fit$param, use.names = TRUE))
  }
  if (is.null(target_names) || length(target_names) != length(fit$param)) {
    target_names <- paste0("par", seq_along(fit$param))
  }

  main_fit_spec <- .GeoVarest_make_fit_spec(
    fit, boot_par, model_est, X_use, optimizer, lower, upper,
    selected_pair_index, coordt_use, memdist = fit_memdist,
    weighted = fit_weighted, include_start = TRUE
  )

  num_params <- length(fit$param)

  estimate_fun <- function(k, current_data) {
    result <- tryCatch({
      capture.output({
        fit_result <- .GeoVarest_run_fit(
          current_data, coords, main_fit_spec, main_fit_spec$param,
          selected_pair_info, fit_context_main, score = FALSE
        )
      }, file = nullfile())

      out <- .GeoVarest_extract_refit(fit_result, target_names)
      rm(fit_result)
      if (k %% 5L == 0L) gc(verbose = FALSE)
      out
    }, error = function(e) {
      if (k <= 3) {
        cat("Warning: iteration", k, "failed:", conditionMessage(e), "\n")
      }
      rep(NA_real_, num_params + 1L)
    })

    result
  }

  ## ------------------------------------------------------------------
  ## 9. Bootstrap
  ## ------------------------------------------------------------------
  param_cols <- seq_len(num_params)

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
      selected_pair_index, coordt_use, memdist = fit_memdist,
      weighted = fit_weighted, include_start = TRUE
    )

    sim_seeds <- sample.int(.Machine$integer.max, K, replace = FALSE)
    ## Preserve the caller's RNG kind across future multisession workers.
    ## future.seed=TRUE switches workers to L'Ecuyer-CMRG; without this
    ## explicit metadata the same scalar replicate seed generates different
    ## TB bootstrap fields in sequential and parallel execution.
    sim_rng_kind <- RNGkind()
    worker_fun <- if (requireNamespace("GeoModels", quietly = TRUE) &&
                      exists(".GeoVarestbootstrap_refit_worker",
                             envir = asNamespace("GeoModels"), inherits = FALSE)) {
      getFromNamespace(".GeoVarestbootstrap_refit_worker", "GeoModels")
    } else {
      .GeoVarestbootstrap_refit_worker
    }

    workers <- if (isTRUE(parallel)) ncores else 1L
    use_outer_stream <- isTRUE(parallel) && workers > 1L &&
      requireNamespace("future", quietly = TRUE) &&
      requireNamespace("future.apply", quietly = TRUE)

    if (isTRUE(parallel) && !use_outer_stream) {
      warning(
        "Parallel TB bootstrap requested but future/future.apply is unavailable; using sequential streaming.",
        call. = FALSE
      )
    }

    if (use_outer_stream) {
      future_backend <- .GeoVarest_set_future_plan(workers)
      chunks <- split(seq_len(K), rep(seq_len(workers), length.out = K))

      payload_size <- as.numeric(utils::object.size(list(
        sim_spec = sim_spec,
        fit_spec = fit_spec,
        coords = coords,
        target_names = target_names,
        sim_seeds = sim_seeds,
        sim_rng_kind = sim_rng_kind
      )))
      old_limit <- getOption("future.globals.maxSize")
      on.exit(options(future.globals.maxSize = old_limit), add = TRUE)
      options(future.globals.maxSize = max(500 * 1024^2, 1.5 * payload_size))

      cat("Performing", K, "TB simulations and bootstrap refits on",
          workers, "cores...\n")

      run_chunks <- function(progressor = NULL) {
        future.apply::future_lapply(
          X = chunks,
          FUN = worker_fun,
          mode = "tb",
          sim_spec = sim_spec,
          fit_spec = fit_spec,
          coords = coords,
          target_names = target_names,
          sim_seeds = sim_seeds,
          L = L,
          rng_kind = sim_rng_kind,
          progressor = progressor,
          future.seed = TRUE,
          future.scheduling = 1,
          future.packages = "GeoModels",
          future.globals = FALSE
        )
      }

      chunk_results <- if (progress) {
        progressr::with_progress({
          pb <- progressr::progressor(steps = K)
          pb_worker <- if (as.numeric(utils::object.size(pb)) <= 1024^2) pb else NULL
          run_chunks(pb_worker)
        })
      } else {
        run_chunks(NULL)
      }

      res_list <- do.call(c, chunk_results)
      if (length(res_list) > 0L && !is.null(names(res_list))) {
        ids <- suppressWarnings(as.integer(names(res_list)))
        if (all(is.finite(ids))) res_list <- res_list[order(ids)]
      }
    } else {
      cat("Performing", K, "TB simulations and bootstrap refits sequentially...\n")
      if (progress) {
        res_list <- progressr::with_progress({
          pb <- progressr::progressor(steps = K)
          unname(worker_fun(
            ids = seq_len(K),
            mode = "tb",
            sim_spec = sim_spec,
            fit_spec = fit_spec,
            coords = coords,
            target_names = target_names,
            sim_seeds = sim_seeds,
            L = L,
            rng_kind = sim_rng_kind,
            progressor = pb
          ))
        })
      } else {
        res_list <- unname(worker_fun(
          ids = seq_len(K),
          mode = "tb",
          sim_spec = sim_spec,
          fit_spec = fit_spec,
          coords = coords,
          target_names = target_names,
          sim_seeds = sim_seeds,
          L = L,
          rng_kind = sim_rng_kind,
          progressor = NULL
        ))
      }
    }

    for (k in seq_len(min(3L, length(res_list)))) {
      err <- attr(res_list[[k]], "error")
      if (!is.null(err)) cat("Warning: iteration", k, "failed:", err, "\n")
    }

    res_raw <- do.call(rbind, res_list)
    if (is.null(res_raw)) {
      res_raw <- matrix(NA_real_, nrow = 0L, ncol = num_params + 1L)
    }

    if (exists("chunk_results", inherits = FALSE)) rm(chunk_results)
    rm(res_list, sim_spec, fit_spec, sim_seeds)
    gc(verbose = FALSE, full = TRUE)

  } else if (!parallel) {
    cat("Performing", K, "estimations sequentially...\n")
    res_list <- vector("list", K)

    if (progress) {
      progressr::with_progress({
        pb <- progressr::progressor(along = seq_len(K))
        for (k in seq_len(K)) {
          pb(sprintf("k=%d", k))
          res_list[[k]] <- suppressWarnings(estimate_fun(k, data_sim[[k]]))
        }
      })
    } else {
      for (k in seq_len(K)) {
        res_list[[k]] <- suppressWarnings(estimate_fun(k, data_sim[[k]]))
      }
    }

    res_raw <- do.call(rbind, res_list)
    if (is.null(res_raw)) {
      res_raw <- matrix(NA_real_, nrow = 0, ncol = num_params + 1L)
    }

  } else {
    cat("Performing", K, "estimations using", ncores, "cores...\n")

    ## Save one dataset per file. Workers receive only paths and a compact
    ## fit specification; no local GeoFit/extraction closures are exported.
    temp_root <- tempfile("GeoVarestbootstrap_")
    dir.create(temp_root, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(temp_root, recursive = TRUE, force = TRUE), add = TRUE)
    data_files <- file.path(temp_root, sprintf("boot_%05d.rds", seq_len(K)))
    for (k in seq_len(K)) {
      saveRDS(data_sim[[k]], data_files[[k]], compress = FALSE)
    }

    rm(data_sim)
    gc(verbose = FALSE, full = TRUE)

    fit_spec <- .GeoVarest_make_fit_spec(
      fit, boot_par, model_est, X_use, optimizer, lower, upper,
      selected_pair_index, coordt_use, memdist = fit_memdist,
      weighted = fit_weighted, include_start = TRUE
    )

    worker_fun <- if (requireNamespace("GeoModels", quietly = TRUE) &&
                      exists(".GeoVarestbootstrap_refit_worker",
                             envir = asNamespace("GeoModels"),
                             inherits = FALSE)) {
      getFromNamespace(".GeoVarestbootstrap_refit_worker", "GeoModels")
    } else {
      .GeoVarestbootstrap_refit_worker
    }

    future::plan(
      future::multisession,
      workers = ncores,
      rscript_libs = .GeoVarest_future_libpaths()
    )
    chunks <- split(seq_len(K), rep(seq_len(ncores), length.out = K))

    payload_size <- as.numeric(utils::object.size(list(
      data_files = data_files,
      fit_spec = fit_spec,
      coords = coords,
      target_names = target_names
    )))
    old_limit <- getOption("future.globals.maxSize")
    on.exit(options(future.globals.maxSize = old_limit), add = TRUE)
    options(future.globals.maxSize = max(500 * 1024^2, 1.5 * payload_size))

    run_refits <- function(progressor = NULL) {
      future.apply::future_lapply(
        X = chunks,
        FUN = worker_fun,
        mode = "file",
        data_files = data_files,
        fit_spec = fit_spec,
        coords = coords,
        target_names = target_names,
        progressor = progressor,
        future.seed = TRUE,
        future.packages = "GeoModels",
        future.globals = FALSE
      )
    }

    chunk_results <- if (progress) {
      progressr::with_progress({
        pb_fit <- progressr::progressor(steps = K)
        pb_worker <- if (as.numeric(utils::object.size(pb_fit)) <= 1024^2) {
          pb_fit
        } else {
          NULL
        }
        run_refits(pb_worker)
      })
    } else {
      run_refits(NULL)
    }

    res_list <- do.call(c, chunk_results)
    if (length(res_list) > 0L && !is.null(names(res_list))) {
      ids <- suppressWarnings(as.integer(names(res_list)))
      if (all(is.finite(ids))) res_list <- res_list[order(ids)]
    }

    for (k in seq_len(min(3L, length(res_list)))) {
      err <- attr(res_list[[k]], "error")
      if (!is.null(err)) {
        cat("Warning: iteration", k, "failed:", err, "\n")
      }
    }

    res_raw <- do.call(rbind, res_list)
    if (is.null(res_raw)) {
      res_raw <- matrix(NA_real_, nrow = 0L, ncol = num_params + 1L)
    }

    rm(res_list, chunk_results, fit_spec)
    gc(verbose = FALSE, full = TRUE)
  }

  ## Filter rows with all NA parameter estimates.
  if (nrow(res_raw) > 0) {
    valid <- !apply(
      res_raw[, param_cols, drop = FALSE],
      1L,
      function(r) all(is.na(r))
    )
    res <- res_raw[valid, , drop = FALSE]
  } else {
    res <- res_raw
  }

  if (ncol(res) == num_params + 1L) {
    colnames(res) <- c(target_names, "logCompLik")
  }

  ## ------------------------------------------------------------------
  ## Clean memory
  ## ------------------------------------------------------------------
  if (exists("data_sim", inherits = FALSE)) rm(data_sim)
  if (exists("res_list", inherits = FALSE)) rm(res_list)
  if (exists("res_raw", inherits = FALSE)) rm(res_raw)
  gc(verbose = FALSE, full = TRUE)

    ## ------------------------------------------------------------------
  ## 10. Post-processing
  ## ------------------------------------------------------------------
  n_successful <- nrow(res)
  .GeoVarest_check_success(
    n_successful, K, min_success_rate, "GeoVarestbootstrap"
  )

  cat("Successful bootstrap iterations:", n_successful, "out of", K, "\n")

  param_mat <- res[, param_cols, drop = FALSE]
  invG      <- stats::var(param_mat)
  G         <- try(solve(invG), silent = TRUE)
  if (!is.matrix(G)) {
    warning("Bootstrap estimated Godambe matrix is singular")
  }

  stderr <- sqrt(pmax(0, diag(invG)))
  names(stderr) <- target_names

  numparam <- length(fit$param)
  penalty <- NA_real_
  claic <- clbic <- NA_real_
  fit$varimat <- NULL
  fit$Jmat <- NULL

  if (!is.null(fit$sensmat)) {
    H <- as.matrix(fit$sensmat)
    storage.mode(H) <- "double"
    if (all(dim(H) == c(numparam, numparam)) && all(is.finite(H))) {
      H <- (H + t(H)) / 2
      penalty <- sum(diag(H %*% invG))
      claic   <- -2 * fit$logCompLik + 2 * penalty
      clbic   <- -2 * fit$logCompLik + log(dimat) * penalty
      ## Under Var(theta_hat)=H^{-1} J H^{-1}, recover J for
      ## compatibility with score-bootstrap output.
      fit$varimat <- H %*% invG %*% H
      fit$Jmat <- fit$varimat
    } else {
      warning(
        "GeoVarestbootstrap: fit$sensmat is incompatible with the bootstrap covariance; CLIC/CLBIC were not computed.",
        call. = FALSE
      )
    }
  }

  ## ------------------------------------------------------------------
  ## 10.1 Bootstrap percentile confidence intervals
  ## ------------------------------------------------------------------
  theta_hat <- as.numeric(fit$param)
  names(theta_hat) <- target_names
  probs <- c((1 - alpha) / 2, 1 - (1 - alpha) / 2)

  q_boot <- t(apply(param_mat, 2L, stats::quantile, probs = probs, na.rm = TRUE))
  colnames(q_boot) <- c("Lower_perc", "Upper_perc")
  fit$conf.int.bootstrap_percentile <- t(q_boot)
  rownames(fit$conf.int.bootstrap_percentile) <- c("Lower", "Upper")

  ## ------------------------------------------------------------------
  ## 11. Final output
  ## ------------------------------------------------------------------
  fit$claic     <- claic
  fit$clic      <- claic
  fit$clbic     <- clbic
  fit$clic_penalty <- penalty
  fit$stderr    <- stderr
  fit$varcov    <- invG
  fit$godambe   <- if (is.matrix(G)) (G + t(G)) / 2 else NULL
  fit$estimates <- res

  z_alpha <- stats::qnorm(1 - (1 - alpha) / 2)
  aa      <- z_alpha * stderr
  pp      <- as.numeric(fit$param)

  fit$conf.int <- rbind(pp - aa, pp + aa)
  colnames(fit$conf.int) <- target_names
  rownames(fit$conf.int) <- c("Lower", "Upper")

  fit$pvalues <- 2 * stats::pnorm(-abs(pp / stderr))
  names(fit$pvalues) <- target_names
  fit$pvalues_type <- "Wald p-values based on parametric-bootstrap variance"

  fit$bootstrap_type <- "parametric"
  fit$bootstrap_K <- K
  fit$bootstrap_successful <- n_successful
  fit$bootstrap_success_rate <- n_successful / K
  fit$bootstrap_min_success_rate <- min_success_rate
  fit$bootstrap_failures <- data.frame(
    iteration = which(!valid),
    error = rep(NA_character_, sum(!valid)),
    stringsAsFactors = FALSE
  )
  fit$pair_index <- NULL
  fit$bootstrap_pair_graph <- if (is.null(selected_pair_index)) "deterministic" else "reused_original"
  fit$bootstrap_thin_method <- fit$thin_method
  fit$bootstrap_p_neighb <- fit$p_neighb
  fit$bootstrap_seed <- seed
  fit$bootstrap_thinning_seed_rule <- if (is.null(selected_pair_index)) {
    "not applicable"
  } else {
    "not used; the original selected-pair graph is reused"
  }
  fit$bootstrap_simulation_method <- method
  fit$bootstrap_memdist <- fit_memdist
  fit$bootstrap_weighted <- fit_weighted
  fit$bootstrap_L <- if (identical(method, "TB")) L else NULL
  fit$bootstrap_workflow <- if (tb_stream_mode) "one-pass TB simulation/refit" else "standard"
  if (tb_stream_mode) {
    fit$bootstrap_execution_mode <- if (isTRUE(use_outer_stream)) {
      "parallel one-pass TB simulation and full refitting"
    } else {
      "sequential one-pass TB simulation and full refitting"
    }
    fit$bootstrap_parallel_workers <- if (isTRUE(use_outer_stream)) workers else 1L
    fit$bootstrap_simulation_seed_rule <-
      "explicit per-replication seeds generated from the optional seed argument"
  }

  return(fit)
}