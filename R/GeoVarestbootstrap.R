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
  ## ------------------------------------------------------------------
  old_handlers <- progressr::handlers()
  if (!is.logical(progress)) stop("progress must be logical (TRUE or FALSE)")
  if (!progress) progressr::handlers("void")
  on.exit(progressr::handlers(old_handlers), add = TRUE)

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
  dimat <- fit$numtime * fit$numcoord
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
  ## 3.1. Memory pre-check (estimate before simulation)
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
      data_sim <- do.call(GeoSim, c(sim_args, list(sparse = sparse, method = method)))
    } else if (method %in% c("TB", "CE")) {
      data_sim <- do.call(GeoSimapprox,
                          c(sim_args, list(method = method, L = L, parallel = FALSE)))
    } else {
      stop("Unsupported method for simulation")
    }
  } else {
    if (method == "cholesky") {
      data_sim <- do.call(GeoSimCopula,
                          c(sim_args, list(copula = fit$copula, sparse = sparse, method = method)))
    } else {
      stop("Unsupported method for copula simulation")
    }
  }

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
  ## 6. Helper function to validate likelihood 
  ## ------------------------------------------------------------------
  is_valid_fit <- function(res_est) {
    if (is.null(res_est)) return(FALSE)

    # Check convergence
    conv <- res_est$convergence
    if (is.null(conv) || !is.character(conv)) return(FALSE)
    if (!(conv %in% c("Successful"))) return(FALSE)

    # Check loglik
    ll <- if (!is.null(res_est$logCompLik)) res_est$logCompLik else res_est$logLik
    if (is.null(ll) || !is.finite(ll)) return(FALSE)
    if (ll == 0 || abs(ll) >= 1e10) return(FALSE)

    # Check parameters
    th <- unlist(res_est$param)
    if (is.null(th) || any(!is.finite(th))) return(FALSE)
    
    TRUE
  }

  ## ------------------------------------------------------------------
  ## 7. Estimation function 
  ## ------------------------------------------------------------------
  estimate_fun <- function(k) {
    result <- tryCatch(
      {
        capture.output({
          # Simple seed based on iteration k for reproducibility
          seed_thin <- as.integer(1234567L + k * 9999L)
          
          fit_result <- if (isTRUE(fit$p_neighb < 1)) {
            withr::with_seed(seed_thin, GeoFit(
              data       = data_sim$data[[k]],
              start      = fit$param,
              fixed      = fit$fixed,
              coordx     = coords,
              coordt     = fit$coordt,
              coordx_dyn = fit$coordx_dyn,
              copula     = fit$copula,
              anisopars  = fit$anisopars,
              est.aniso  = fit$est.aniso,
              lower      = lower,
              upper      = upper,
              neighb     = fit$neighb,
              p_neighb   = fit$p_neighb,
              corrmodel  = fit$corrmodel,
              model      = model_est,
              sparse     = FALSE,
              n          = fit$n,
              maxdist    = fit$maxdist,
              maxtime    = fit$maxtime,
              optimizer  = optimizer,
              grid       = fit$grid,
              likelihood = fit$likelihood,
              type       = fit$type,
              X          = X_use,
              distance   = fit$distance,
              radius     = fit$radius
            ))
          } else {
            GeoFit(
              data       = data_sim$data[[k]],
              start      = fit$param,
              fixed      = fit$fixed,
              coordx     = coords,
              coordt     = fit$coordt,
              coordx_dyn = fit$coordx_dyn,
              copula     = fit$copula,
              anisopars  = fit$anisopars,
              est.aniso  = fit$est.aniso,
              lower      = lower,
              upper      = upper,
              neighb     = fit$neighb,
              p_neighb   = fit$p_neighb,
              corrmodel  = fit$corrmodel,
              model      = model_est,
              sparse     = FALSE,
              n          = fit$n,
              maxdist    = fit$maxdist,
              maxtime    = fit$maxtime,
              optimizer  = optimizer,
              grid       = fit$grid,
              likelihood = fit$likelihood,
              type       = fit$type,
              X          = X_use,
              distance   = fit$distance,
              radius     = fit$radius
            )
          }
        }, file = nullfile())

        if (is_valid_fit(fit_result)) {
          ll <- if (!is.null(fit_result$logCompLik)) fit_result$logCompLik else NA_real_
          c(unlist(fit_result$param), logCompLik = ll)
        } else {
          NULL
        }
      },
      error = function(e) {
        if (k <= 3) cat("Warning: iteration", k, "failed:", conditionMessage(e), "\n")
        NULL
      }
    )
    if (k %% 10 == 0) gc(verbose = FALSE)
    result
  }

  ## ------------------------------------------------------------------
  ## 8. Bootstrap 
  ## ------------------------------------------------------------------
  num_params <- length(fit$param)
  param_cols <- seq_len(num_params)

  if (!parallel) {
    cat("Performing", K, "estimations sequentially...\n")
    if (progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      res_list <- vector("list", K)
      progressr::with_progress({
        pb <- progressr::progressor(along = seq_len(K))
        for (k in seq_len(K)) {
          pb(sprintf("k=%d", k))
          res_list[[k]] <- suppressWarnings(estimate_fun(k))
        }
      })
    } else {
      res_list <- vector("list", K)
      for (k in seq_len(K)) {
        res_list[[k]] <- suppressWarnings(estimate_fun(k))
      }
    }
    res <- do.call(rbind, Filter(Negate(is.null), res_list))
    if (is.null(res)) res <- matrix(NA_real_, nrow = 0, ncol = num_params + 1L)

  } else {
    cat("Performing", K, "estimations using", ncores, "cores...\n")

    data_size_mb <- as.numeric(object.size(data_sim)) / 1024^2
    future_limit_mb <- max(600, ceiling(data_size_mb * 2))
    old_limit <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = future_limit_mb * 1024^2)
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
            .combine = rbind,
            .errorhandling = "pass",
            .options.future = list(
              packages = c("GeoModels"),
              seed = TRUE
            )
          ) %dofuture% {
            pb(sprintf("k=%d", k))
            res_est <- suppressWarnings(estimate_fun(k))
            if (is.null(res_est)) {
              rep(NA_real_, num_params + 1L)
            } else {
              res_est
            }
          }
        })
      } else {
        foreach::foreach(
          k = seq_len(K),
          .combine = rbind,
          .errorhandling = "pass",
          .options.future = list(
            packages = c("GeoModels"),
            seed = TRUE
          )
        ) %dofuture% {
          res_est <- suppressWarnings(estimate_fun(k))
          if (is.null(res_est)) {
            rep(NA_real_, num_params + 1L)
          } else {
            res_est
          }
        }
      }
    }, error = function(e) {
      warning("Parallel bootstrap failed: ", e$message)
      matrix(NA_real_, nrow = 0, ncol = num_params + 1L)
    })

    # More robust check for parallel results
    if (!is.matrix(xx) || is.null(dim(xx)) || nrow(xx) == 0) {
      res <- matrix(NA_real_, nrow = 0, ncol = num_params + 1L)
    } else {
      valid <- !apply(xx[, param_cols, drop = FALSE], 1L, function(row) all(is.na(row)))
      res   <- xx[valid, , drop = FALSE]
    }
  }

  if (ncol(res) == num_params + 1L) {
    colnames(res) <- c(names(fit$param), "logCompLik")
  }

  ## ------------------------------------------------------------------
  ## Clean memory 
  ## ------------------------------------------------------------------
  rm(data_sim)
  if (exists("xx", inherits = FALSE)) rm(xx)
  if (exists("res_list", inherits = FALSE)) rm(res_list)
  gc(verbose = FALSE, full = TRUE)

  ## ------------------------------------------------------------------
  ## 9. Post-processing 
  ## ------------------------------------------------------------------
  n_successful <- nrow(res)
  if (n_successful < 2) stop("Insufficient successful bootstrap iterations for variance estimation")
  cat("Successful bootstrap iterations:", n_successful, "out of", K, "\n")

  param_mat <- res[, param_cols, drop = FALSE]
  invG      <- var(param_mat)
  G         <- try(solve(invG), silent = TRUE)
  if (!is.matrix(G)) warning("Bootstrap estimated Godambe matrix is singular")

  stderr <- sqrt(pmax(0, diag(invG)))
  names(stderr) <- names(fit$param)

  numparam <- length(fit$param)
  if ((fit$likelihood == "Marginal" && fit$type %in% c("Independence", "Pairwise")) ||
      (fit$likelihood == "Conditional" && fit$type == "Pairwise")) {
    H       <- fit$sensmat
    penalty <- sum(diag(H %*% invG))
    claic   <- -2 * fit$logCompLik + 2 * penalty
    clbic   <- -2 * fit$logCompLik + log(dimat) * penalty
    fit$varimat <- H %*% invG %*% H
  } else if (fit$likelihood == "Full" && fit$type == "Standard") {
    claic <- -2 * fit$logCompLik + 2 * numparam
    clbic <- -2 * fit$logCompLik + log(dimat) * numparam
  } else {
    claic <- clbic <- NA_real_
  }

  ## ------------------------------------------------------------------
  ## 9.1 Confidence intervals and bootstrap p-values
  ## ------------------------------------------------------------------
  theta_hat <- as.numeric(fit$param)
  names(theta_hat) <- names(fit$param)
  probs <- c((1 - alpha)/2, 1 - (1 - alpha)/2)

  q_boot <- t(apply(param_mat, 2L, quantile, probs = probs, na.rm = TRUE))
  colnames(q_boot) <- c("Lower_perc", "Upper_perc")
  confint_boot_perc <- q_boot

  # Two-tailed p-values against zero
  theta0 <- rep(0, length(theta_hat))
  names(theta0) <- names(theta_hat)
  pval_boot <- vapply(seq_along(theta_hat), function(j) {
    bj <- param_mat[, j]
    pj_low  <- mean(bj <= theta0[j], na.rm = TRUE)
    pj_high <- mean(bj >= theta0[j], na.rm = TRUE)
    2 * min(pj_low, pj_high)
  }, numeric(1))
  names(pval_boot) <- names(theta_hat)

  fit$conf.int.bootstrap_percentile <- t(confint_boot_perc)
  rownames(fit$conf.int.bootstrap_percentile) <- c("Lower", "Upper")
  fit$pvalues_bootstrap <- pmin(1, pval_boot)

  ## ------------------------------------------------------------------
  ## 10. Final output
  ## ------------------------------------------------------------------
  fit$claic      <- claic
  fit$clbic      <- clbic
  fit$stderr     <- stderr
  fit$varcov     <- invG
  fit$estimates  <- res

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