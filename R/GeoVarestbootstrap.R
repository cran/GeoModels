GeoVarestbootstrap <- function(fit, K = 100, sparse = FALSE,
                               optimizer = NULL, lower = NULL, upper = NULL,
                               method = "cholesky", alpha = 0.95,
                               L = 1000, parallel = TRUE, ncores = NULL,
                               progress = TRUE) {

  ## ------------------------------------------------------------------
  ## Setup cleanup handlers
  ## ------------------------------------------------------------------
  future_plan_original <- future::plan()

  cleanup_resources <- function() {
    tryCatch({
      future::plan(future_plan_original)
    }, error = function(e) {
      tryCatch(future::plan(future::sequential),
               error = function(e2) invisible(NULL))
    })
    gc(verbose = FALSE, full = TRUE)
  }
  on.exit(cleanup_resources(), add = TRUE)

  ## ------------------------------------------------------------------
  ## Progress handlers
  ## FIX 7: restore correttamente anche quando old_handlers è lista vuota
  ## ------------------------------------------------------------------
  old_handlers <- progressr::handlers()
  if (!is.logical(progress)) stop("progress must be logical (TRUE or FALSE)")
  if (!progress) progressr::handlers("void")
  on.exit({
    if (length(old_handlers) > 0) progressr::handlers(old_handlers)
    else progressr::handlers("default")
  }, add = TRUE)

  ## ------------------------------------------------------------------
  ## 1. Check input
  ## ------------------------------------------------------------------
  if (length(fit$coordt) == 1) fit$coordt <- NULL
  if (is.null(fit$sensmat)) stop("Sensitivity matrix is missing: use sensitivity = TRUE in GeoFit")
  if (!(method %in% c("cholesky", "TB", "CE"))) stop("The method of simulation is not correct")
  if (!is.numeric(K) || K < 1) stop("K must be a positive integer")
  K <- as.integer(K)
  if (!(is.numeric(alpha) && length(alpha) == 1 && alpha > 0 && alpha < 1))
    stop("alpha must be a single numeric in (0,1)")

  if (is.null(optimizer)) {
    optimizer <- fit$optimizer
    lower     <- fit$lower
    upper     <- fit$upper
  }

  model <- fit$model
  cat("Parametric bootstrap can be time consuming ...\n")

  ## ------------------------------------------------------------------
  ## 2. Misspecification mapping
  ## ------------------------------------------------------------------
  model_sim <- model
  model_est <- model

  if (isTRUE(fit$missp)) {
    model_map <- c(
      StudentT     = "Gaussian_misp_StudentT",
      Poisson      = "Gaussian_misp_Poisson",
      PoissonZIP   = "Gaussian_misp_PoissonZIP",
      SkewStudentT = "Gaussian_misp_SkewStudentT",
      Tukeygh      = "Gaussian_misp_Tukeygh"
    )
    if (model %in% names(model_map)) {
      model_est <- model_map[[model]]
    }
  }

  ## ------------------------------------------------------------------
  ## 3. X matrix computing
  ## ------------------------------------------------------------------
  ## FIX 9: dimat corretto per caso dinamico (coordx_dyn)
  if (!is.null(fit$coordx_dyn)) {
    dimat <- sum(fit$ns)
  } else {
    dimat <- fit$numtime * fit$numcoord
  }

  if (!is.null(fit$X) && ncol(fit$X) == 1) {
    ncheck <- min(NROW(fit$X), dimat)
    X_use <- if (ncheck > 0 && all(fit$X[seq_len(ncheck), 1] == 1)) NULL else fit$X
  } else {
    X_use <- fit$X
  }

  coords <- cbind(fit$coordx, fit$coordy)
  if (isTRUE(fit$bivariate) && is.null(fit$coordx_dyn)) {
    if (nrow(coords) %% 2L != 0L) stop("bivariate=TRUE but odd number of coordinates")
    coords <- coords[seq_len(nrow(coords) / 2L), , drop = FALSE]
  }

  ## ------------------------------------------------------------------
  ## 3.1. Memory pre-check
  ## ------------------------------------------------------------------
  estimated_size_mb <- (K * fit$numtime * fit$numcoord * 8) / (1024^2)
  if (estimated_size_mb > 500) {
    warning(sprintf(
      "Estimated simulated dataset size: %.1f MB (>500MB). Consider reducing K or using sparse methods.",
      estimated_size_mb
    ))
  }

  ## ------------------------------------------------------------------
  ## 4. Simulation
  ## ------------------------------------------------------------------
  sim_args <- list(
    coordx     = coords,
    coordt     = fit$coordt,
    coordx_dyn = fit$coordx_dyn,
    anisopars  = fit$anisopars,
    corrmodel  = fit$corrmodel,
    model      = model_sim,
    param      = append(fit$param, fit$fixed),
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
      data_sim_full <- do.call(GeoSimapprox,
                          c(sim_args, list(method = method, L = L, parallel = FALSE)))
    } else {
      stop("Unsupported method for simulation")
    }
  } else {
    if (method == "cholesky") {
      data_sim_full <- do.call(GeoSimCopula,
                          c(sim_args, list(copula = fit$copula, sparse = sparse, method = method)))
    } else {
      stop("Unsupported method for copula simulation")
    }
  }

  data_sim <- data_sim_full$data
  rm(data_sim_full)
  gc(verbose = FALSE, full = TRUE)

  ## ------------------------------------------------------------------
  ## 5. Setting parallelization
  ## ------------------------------------------------------------------
  coremax <- parallel::detectCores()
  if (is.na(coremax) || coremax <= 1) {
    parallel <- FALSE
    ncores   <- 1
  } else {
    if (is.null(parallel)) parallel <- TRUE
    ncores <- max(1L, min(if (is.null(ncores)) coremax - 1L else ncores, K))
  }

  ## ------------------------------------------------------------------
  ## 6. Helper: estrai solo param — UNICA FUNZIONE, usata ovunque
  ## FIX 4: logica di validazione unificata, elimina divergenza
  ##        tra percorso sequenziale e worker parallelo
  ## ------------------------------------------------------------------
  extract_param <- function(res_est, num_params) {
    if (is.null(res_est)) return(rep(NA_real_, num_params + 1L))

    conv <- res_est$convergence
    if (is.null(conv) || !is.character(conv) || conv != "Successful")
      return(rep(NA_real_, num_params + 1L))

    ll <- if (!is.null(res_est$logCompLik)) res_est$logCompLik
          else if (!is.null(res_est$logLik)) res_est$logLik
          else return(rep(NA_real_, num_params + 1L))

    if (!is.finite(ll) || ll == 0 || abs(ll) >= 1e10)
      return(rep(NA_real_, num_params + 1L))

    th <- unlist(res_est$param)
    if (is.null(th) || length(th) != num_params || any(!is.finite(th)))
      return(rep(NA_real_, num_params + 1L))

    c(th, logCompLik = ll)
  }

  ## ------------------------------------------------------------------
  ## 7. Helper: chiama GeoFit evitando duplicazione
  ## FIX 3: unica funzione, il seed viene applicato esternamente se serve
  ## ------------------------------------------------------------------
  call_geofit <- function(current_data, fit_ess, coords, X_use,
                          lower, upper, optimizer, model_est) {
    GeoFit(
      data       = current_data,
      start      = fit_ess$param,
      fixed      = fit_ess$fixed,
      coordx     = coords,
      coordt     = fit_ess$coordt,
      coordx_dyn = fit_ess$coordx_dyn,
      copula     = fit_ess$copula,
      anisopars  = fit_ess$anisopars,
      est.aniso  = fit_ess$est.aniso,
      lower      = lower,
      upper      = upper,
      neighb     = fit_ess$neighb,
      p_neighb   = fit_ess$p_neighb,
      corrmodel  = fit_ess$corrmodel,
      model      = model_est,
      sparse     = FALSE,
      n          = fit_ess$n,
      maxdist    = fit_ess$maxdist,
      maxtime    = fit_ess$maxtime,
      optimizer  = optimizer,
      grid       = fit_ess$grid,
      likelihood = fit_ess$likelihood,
      type       = fit_ess$type,
      X          = X_use,
      distance   = fit_ess$distance,
      radius     = fit_ess$radius
    )
  }

  ## ------------------------------------------------------------------
  ## 8. Estimation function (sequenziale)
  ## ------------------------------------------------------------------
  num_params <- length(fit$param)

  estimate_fun <- function(k, current_data) {
    result <- tryCatch({
      capture.output({
        seed_thin <- as.integer(1234567L + k * 9999L)
        fit_result <- if (isTRUE(fit$p_neighb < 1)) {
          withr::with_seed(seed_thin,
            call_geofit(current_data, fit, coords, X_use,
                        lower, upper, optimizer, model_est))
        } else {
          call_geofit(current_data, fit, coords, X_use,
                      lower, upper, optimizer, model_est)
        }
      }, file = nullfile())

      out <- extract_param(fit_result, num_params)
      rm(fit_result)
      ## FIX 8: gc senza full=TRUE nei percorsi interni
      if (k %% 10 == 0) gc(verbose = FALSE)
      out
    },
    error = function(e) {
      if (k <= 3) cat("Warning: iteration", k, "failed:", conditionMessage(e), "\n")
      rep(NA_real_, num_params + 1L)
    })
    result
  }

  ## ------------------------------------------------------------------
  ## 9. Bootstrap
  ## ------------------------------------------------------------------
  param_cols <- seq_len(num_params)

  if (!parallel) {
    cat("Performing", K, "estimations sequentially...\n")
    res_list <- vector("list", K)

    if (progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
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
    if (is.null(res_raw)) res_raw <- matrix(NA_real_, nrow = 0, ncol = num_params + 1L)

  } else {
    cat("Performing", K, "estimations using", ncores, "cores...\n")

    ## ------------------------------------------------------------------
    ## FIX 1: salva K file separati — ogni worker legge solo il suo
    ## ------------------------------------------------------------------
    temp_dir  <- tempdir()
    ts_tag    <- format(Sys.time(), "%Y%m%d_%H%M%S")
    data_files <- vapply(seq_len(K), function(k) {
      fp <- file.path(temp_dir, sprintf("boot_%s_%04d.rds", ts_tag, k))
      saveRDS(data_sim[[k]], fp, compress = FALSE)
      fp
    }, character(1))
    on.exit(unlink(data_files), add = TRUE)

    rm(data_sim)
    gc(verbose = FALSE, full = TRUE)

    fit_essentials <- list(
      param      = fit$param,
      fixed      = fit$fixed,
      coordt     = fit$coordt,
      coordx_dyn = fit$coordx_dyn,
      copula     = fit$copula,
      anisopars  = fit$anisopars,
      est.aniso  = fit$est.aniso,
      neighb     = fit$neighb,
      p_neighb   = fit$p_neighb,
      corrmodel  = fit$corrmodel,
      n          = fit$n,
      maxdist    = fit$maxdist,
      maxtime    = fit$maxtime,
      grid       = fit$grid,
      likelihood = fit$likelihood,
      type       = fit$type,
      distance   = fit$distance,
      radius     = fit$radius
    )

    ## Worker ultra-leggero: legge UN solo file
    estimate_worker <- function(k, data_files, fit_ess, coords,
                                X_use, lower, upper, optimizer,
                                model_est, num_params,
                                call_geofit, extract_param) {
      current_data <- readRDS(data_files[[k]])

      result <- tryCatch({
        capture.output({
          seed_thin  <- as.integer(1234567L + k * 9999L)
          fit_result <- if (isTRUE(fit_ess$p_neighb < 1)) {
            withr::with_seed(seed_thin,
              call_geofit(current_data, fit_ess, coords, X_use,
                          lower, upper, optimizer, model_est))
          } else {
            call_geofit(current_data, fit_ess, coords, X_use,
                        lower, upper, optimizer, model_est)
          }
        }, file = nullfile())

        out <- extract_param(fit_result, num_params)
        rm(fit_result, current_data)
        ## FIX 8: gc leggero nel worker
        gc(verbose = FALSE)
        out
      },
      error = function(e) {
        if (k <= 3) cat("Warning: iteration", k, "failed:", conditionMessage(e), "\n")
        rep(NA_real_, num_params + 1L)
      })
      result
    }

    old_limit <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = 200 * 1024^2)
    on.exit(options(future.globals.maxSize = old_limit), add = TRUE)

    future::plan(future::multisession, workers = ncores)

    if (progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
    }

    xx <- tryCatch({
      if (progress) {
        progressr::with_progress({
          pb <- progressr::progressor(along = seq_len(K))
          foreach::foreach(
            k = seq_len(K),
            .combine  = rbind,
            .errorhandling = "pass",
            .options.future = list(
              packages = c("GeoModels", "withr"),
              seed     = TRUE,
              globals  = structure(TRUE, add = c(
                "data_files", "fit_essentials", "coords", "X_use",
                "lower", "upper", "optimizer", "model_est",
                "num_params", "estimate_worker", "call_geofit", "extract_param"
              ))
            )
          ) %dofuture% {
            pb(sprintf("k=%d", k))
            estimate_worker(k, data_files, fit_essentials, coords,
                            X_use, lower, upper, optimizer,
                            model_est, num_params, call_geofit, extract_param)
          }
        })
      } else {
        foreach::foreach(
          k = seq_len(K),
          .combine  = rbind,
          .errorhandling = "pass",
          .options.future = list(
            packages = c("GeoModels", "withr"),
            seed     = TRUE,
            globals  = structure(TRUE, add = c(
              "data_files", "fit_essentials", "coords", "X_use",
              "lower", "upper", "optimizer", "model_est",
              "num_params", "estimate_worker", "call_geofit", "extract_param"
            ))
          )
        ) %dofuture% {
          estimate_worker(k, data_files, fit_essentials, coords,
                          X_use, lower, upper, optimizer,
                          model_est, num_params, call_geofit, extract_param)
        }
      }
    }, error = function(e) {
      warning("Parallel bootstrap failed: ", e$message)
      matrix(NA_real_, nrow = 0, ncol = num_params + 1L)
    })

    res_raw <- if (!is.matrix(xx) || nrow(xx) == 0)
      matrix(NA_real_, nrow = 0, ncol = num_params + 1L)
    else xx
  }

  ## Filtra righe con tutti NA sui parametri
  if (nrow(res_raw) > 0) {
    valid <- !apply(res_raw[, param_cols, drop = FALSE], 1L,
                    function(r) all(is.na(r)))
    res <- res_raw[valid, , drop = FALSE]
  } else {
    res <- res_raw
  }

  if (ncol(res) == num_params + 1L)
    colnames(res) <- c(names(fit$param), "logCompLik")

  ## ------------------------------------------------------------------
  ## Clean memory
  ## ------------------------------------------------------------------
  if (exists("data_sim",  inherits = FALSE)) rm(data_sim)
  if (exists("xx",        inherits = FALSE)) rm(xx)
  if (exists("res_list",  inherits = FALSE)) rm(res_list)
  if (exists("res_raw",   inherits = FALSE)) rm(res_raw)
  gc(verbose = FALSE, full = TRUE)

  ## ------------------------------------------------------------------
  ## 10. Post-processing
  ## ------------------------------------------------------------------
  n_successful <- nrow(res)
  if (n_successful < 2)
    stop("Insufficient successful bootstrap iterations for variance estimation")
  cat("Successful bootstrap iterations:", n_successful, "out of", K, "\n")

  param_mat <- res[, param_cols, drop = FALSE]
  invG      <- var(param_mat)          # stima bootstrap di G⁻¹ (paper, sez. 2)
  G         <- try(solve(invG), silent = TRUE)
  if (!is.matrix(G)) warning("Bootstrap estimated Godambe matrix is singular")

  stderr <- sqrt(pmax(0, diag(invG)))
  names(stderr) <- names(fit$param)

  numparam <- length(fit$param)
  if ((fit$likelihood == "Marginal" && fit$type %in% c("Independence", "Pairwise")) ||
      (fit$likelihood == "Conditional" && fit$type == "Pairwise")) {
    H       <- fit$sensmat
    penalty <- sum(diag(H %*% invG))   # tr(H · G⁻¹), eq. (6) del paper
    claic   <- -2 * fit$logCompLik + 2 * penalty
    clbic   <- -2 * fit$logCompLik + log(dimat) * penalty
    fit$varimat <- H %*% invG %*% H    # stima di J_a
  } else if (fit$likelihood == "Full" && fit$type == "Standard") {
    claic <- -2 * fit$logCompLik + 2 * numparam
    clbic <- -2 * fit$logCompLik + log(dimat) * numparam
  } else {
    claic <- clbic <- NA_real_
  }

  ## ------------------------------------------------------------------
  ## 10.1 Confidence intervals e bootstrap p-values
  ## FIX 5: p-value calcolato solo per parametri che ammettono H₀:θ=0
  ##        (regressione); per parametri di scala/dipendenza viene
  ##        restituito NA con un avviso.
  ## ------------------------------------------------------------------
  theta_hat <- as.numeric(fit$param)
  names(theta_hat) <- names(fit$param)
  probs <- c((1 - alpha) / 2, 1 - (1 - alpha) / 2)

  q_boot <- t(apply(param_mat, 2L, quantile, probs = probs, na.rm = TRUE))
  colnames(q_boot) <- c("Lower_perc", "Upper_perc")
  fit$conf.int.bootstrap_percentile <- t(q_boot)
  rownames(fit$conf.int.bootstrap_percentile) <- c("Lower", "Upper")

  ## Parametri positivi per costruzione: test contro 0 non ha senso
  positive_params <- c("sill", "sill_1", "sill_2", "nugget", "nugget_1",
                       "nugget_2", "scale", "scale_1", "scale_2",
                       "smooth", "shape", "df", "sigma2", "tail",
                       "tail1", "tail2")
  pval_boot <- vapply(seq_along(theta_hat), function(j) {
    nm <- names(theta_hat)[j]
    if (any(vapply(positive_params, function(p) grepl(p, nm), logical(1)))) {
      return(NA_real_)
    }
    bj      <- param_mat[, j]
    pj_low  <- mean(bj <= 0, na.rm = TRUE)
    pj_high <- mean(bj >= 0, na.rm = TRUE)
    pmin(1, 2 * min(pj_low, pj_high))
  }, numeric(1))
  names(pval_boot) <- names(theta_hat)

  fit$pvalues_bootstrap <- pval_boot

  ## ------------------------------------------------------------------
  ## 11. Final output
  ## ------------------------------------------------------------------
  fit$claic     <- claic
  fit$clbic     <- clbic
  fit$stderr    <- stderr
  fit$varcov    <- invG
  fit$estimates <- res

  z_alpha <- qnorm(1 - (1 - alpha) / 2)
  aa      <- z_alpha * stderr
  pp      <- as.numeric(fit$param)
  fit$conf.int <- rbind(pp - aa, pp + aa)
  colnames(fit$conf.int) <- names(fit$param)
  rownames(fit$conf.int) <- c("Lower", "Upper")

  fit$pvalues <- 2 * pnorm(-abs(pp / stderr))
  names(fit$pvalues) <- names(fit$param)

  return(fit)
}