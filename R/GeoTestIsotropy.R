GeoTestIsotropy <- function(data, coordx,
                            start, fixed,
                            optimizer = "bobyqa",
                            model = "Gaussian",
                            corrmodel = "Matern",
                            lower = NULL, upper = NULL,
                            B = 1000,
                            likelihood = NULL,
                            type = NULL,
                            copula = NULL,
                            neighb = 5,
                            method = c("cholesky", "TB"),
                            L = 1000,
                            parallel = TRUE,
                            ncores = NULL,
                            progress = TRUE,
                            seed = NULL) {

  ## ====== CLEANUP: FUTURE PLAN AND PROGRESS HANDLERS ======
  future_plan_original <- future::plan()
  old_handlers <- progressr::handlers()

  cleanup_resources <- function() {
    tryCatch({
      future::plan(future_plan_original)
    }, error = function(e) invisible(NULL))

    tryCatch({
      if (length(old_handlers) > 0L) {
        progressr::handlers(old_handlers)
      } else {
        progressr::handlers("default")
      }
    }, error = function(e) invisible(NULL))

    gc(verbose = FALSE, full = TRUE)
  }
  on.exit(cleanup_resources(), add = TRUE)

  ## ====== RNG SEED, OPTIONAL AND NON-INTRUSIVE ======
  ## If seed is supplied, make simulations/bootstrap reproducible,
  ## then restore the user's original .Random.seed on exit.
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("seed must be NULL or a single finite numeric value", call. = FALSE)
    }

    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }

    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)

    set.seed(as.integer(seed))
  }

  set_progress_handlers <- function() {
    if (isTRUE(progress)) {
      progressr::handlers(global = TRUE)
      progressr::handlers(progressr::handler_txtprogressbar(clear = TRUE))
    } else {
      progressr::handlers("void")
    }
  }

  ## ====== VALIDAZIONE INPUT INIZIALE ======
  method <- match.arg(method)

  if (!is.matrix(coordx)) coordx <- as.matrix(coordx)
  if (!is.numeric(data))  data   <- as.numeric(data)

  if (nrow(coordx) != length(data)) {
    stop("coordx rows (", nrow(coordx), ") must match data length (", length(data), ")")
  }
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 1) {
    stop("B must be a positive integer")
  }
  B <- as.integer(B)

  if (!is.numeric(L) || length(L) != 1L || !is.finite(L) || L < 1) {
    stop("L must be a positive integer")
  }
  L <- as.integer(L)

  if (!is.numeric(neighb) || length(neighb) != 1L || !is.finite(neighb) || neighb < 1) {
    stop("neighb must be positive")
  }

  if (!is.logical(progress) || length(progress) != 1L) {
    stop("progress must be logical (TRUE or FALSE)")
  }

  ## Numero di punti
  NN <- nrow(coordx)

  ## ====== GESTIONE CORES ROBUSTA ======
  coremax <- parallel::detectCores()
  if (is.na(coremax) || coremax <= 1L) {
    parallel <- FALSE
    ncores   <- 1L
  } else {
    if (is.null(ncores)) {
      ncores <- min(coremax - 1L, B)
    }
    if (!is.numeric(ncores) || length(ncores) != 1L || !is.finite(ncores)) {
      stop("ncores must be a positive integer or NULL")
    }
    ncores <- as.integer(max(1L, min(ncores, B, coremax)))
    if (ncores == 1L) parallel <- FALSE
  }

  ## ====== AUTOMATIC LIKELIHOOD SELECTION ======
  if (is.null(likelihood) && is.null(type)) {
    if (NN > 2000) {
      likelihood <- "Marginal"
      type <- "Pairwise"
    } else {
      likelihood <- "Full"
      type <- "Standard"
      neighb <- NULL
    }
  } else if (xor(is.null(likelihood), is.null(type))) {
    stop("Either specify both 'likelihood' and 'type', or leave both NULL for automatic selection.")
  }

  set_progress_handlers()

  corrmodel_used <- corrmodel
  model_used <- model

  ## ====== VALIDAZIONE PARAMETRI ======
  corr_names <- CorrParam(corrmodel_used)
  nuis_names <- NuisParam(model_used)
  expected <- sort(unique(c(corr_names, nuis_names)))
  provided <- sort(unique(c(names(start), names(fixed))))
  if (!identical(expected, provided)) {
    missing <- setdiff(expected, provided)
    extra   <- setdiff(provided, expected)
    msg <- c()
    if (length(missing)) msg <- c(msg, paste0("missing: ", paste(missing, collapse = ", ")))
    if (length(extra))   msg <- c(msg, paste0("extra: ",   paste(extra,   collapse = ", ")))
    stop("Parameter names in start+fixed must match CorrParam(corrmodel)+NuisParam(model). ",
         paste(msg, collapse = " | "))
  }

  ## Bounds solo per 'start'
  start_names <- sort(names(start))
  norm_bounds <- function(b, which = c("lower", "upper")) {
    which <- match.arg(which)
    if (is.null(b)) b <- list()
    b_names <- names(b)
    if (length(b) && is.null(b_names)) {
      stop(which, " must be a *named* list with names matching start")
    }
    if (length(setdiff(b_names, start_names))) {
      stop(which, " contains names not in start: ",
           paste(setdiff(b_names, start_names), collapse = ", "))
    }
    out <- setNames(rep(if (which == "lower") -Inf else Inf, length(start_names)), start_names)
    if (length(b_names)) out[b_names] <- unlist(b, use.names = FALSE)
    as.list(out)
  }
  lower <- norm_bounds(lower, "lower")
  upper <- norm_bounds(upper, "upper")

  ## ====== ANISOTROPIA ======
  aniso_H0    <- list(angle = 0, ratio = 1)
  aniso_start <- list(angle = 0, ratio = 1.2)

  ## ====== HELPERS DI FIT ======
  fit_H0 <- function(Z) {
    GeoFit(
      data = Z, coordx = coordx, corrmodel = corrmodel_used,
      likelihood = likelihood, optimizer = optimizer,
      model = model_used, type = type, neighb = neighb,
      start = start, fixed = fixed,
      copula = copula,
      lower = lower, upper = upper,
      anisopars = aniso_H0
    )
  }

  fit_H1 <- function(Z) {
    GeoFit(
      data = Z, coordx = coordx, corrmodel = corrmodel_used,
      likelihood = likelihood, optimizer = optimizer,
      model = model_used, type = type, neighb = neighb,
      start = start, fixed = fixed,
      copula = copula,
      lower = lower, upper = upper,
      anisopars = aniso_start,
      est.aniso = c(TRUE, TRUE)
    )
  }

  ## ====== FIT OSSERVATO ======
  fit0 <- fit_H0(data)
  fit1 <- fit_H1(data)

  LL0 <- if (!is.null(fit0$logCompLik)) fit0$logCompLik else fit0$logLik
  LL1 <- if (!is.null(fit1$logCompLik)) fit1$logCompLik else fit1$logLik
  LL0 <- as.numeric(LL0)
  LL1 <- as.numeric(LL1)

  if (!is.finite(LL0) || !is.finite(LL1)) {
    stop("Observed H0/H1 fit returned non-finite log-likelihood values")
  }

  if (LL1 < LL0) {
    message("Note: H1 log-likelihood (", round(LL1, 3),
            ") < H0 log-likelihood (", round(LL0, 3), "), LRT set to 0")
  }
  LRT_obs <- max(0, 2 * (LL1 - LL0))

  ## ====== SIMULAZIONI H0 ======
  par0 <- c(fit0$param, fit0$fixed)
  param_H0 <- as.list(par0[expected])

  if (is.null(copula)) {
    if (method == "cholesky") {
      sim <- GeoSim(
        coordx = coordx, corrmodel = corrmodel_used,
        param = param_H0, model = model_used,
        nrep = B, anisopars = aniso_H0, progress = FALSE
      )
    } else if (method == "TB") {
      ## GeoSimapprox may parallelize over nrep = B. This phase is completed
      ## before bootstrap refits start, so the refit parallelism below is not nested.
      plan_before_sim <- future::plan()

      sim <- GeoSimapprox(
        coordx = coordx, corrmodel = corrmodel_used, model = model_used,
        method = "TB", L = L,
        param = param_H0, nrep = B, anisopars = aniso_H0,
        progress = FALSE,
        parallel = parallel,
        ncores = ncores
      )

      try(future::plan(plan_before_sim), silent = TRUE)
    } else {
      stop("method must be 'cholesky' or 'TB'.")
    }
  } else {
    if (method != "cholesky") {
      stop("Unsupported method for copula simulation: use method = 'cholesky'")
    }

    sim <- GeoSimCopula(
      coordx = coordx, corrmodel = corrmodel_used,
      method = "Cholesky",
      param = param_H0, model = model_used, copula = copula,
      nrep = B, anisopars = aniso_H0
    )
  }

  get_rep <- function(simobj, b) {
    X <- simobj$data
    if (is.list(X)) X[[b]] else X[, b]
  }

  loglik_value <- function(fit_obj) {
    val <- if (!is.null(fit_obj$logCompLik)) fit_obj$logCompLik else fit_obj$logLik
    as.numeric(val)
  }

  ## ====== BOOTSTRAP ======
  if (!parallel) {
    message(sprintf("Performing %d estimations under H0 and H1 sequentially...", B))

    set_progress_handlers()
    pb <- if (isTRUE(progress)) progressr::progressor(along = seq_len(B)) else NULL

    LRT_boot <- vector("numeric", B)
    for (b in seq_len(B)) {
      if (!is.null(pb)) pb(sprintf("b=%d", b))

      Zb  <- get_rep(sim, b)
      f0b <- suppressWarnings(try(fit_H0(Zb), silent = TRUE))
      f1b <- suppressWarnings(try(fit_H1(Zb), silent = TRUE))

      if (inherits(f0b, "try-error") || inherits(f1b, "try-error")) {
        LRT_boot[b] <- NA_real_
      } else {
        ll0b <- loglik_value(f0b)
        ll1b <- loglik_value(f1b)
        if (!is.finite(ll0b) || !is.finite(ll1b)) {
          LRT_boot[b] <- NA_real_
        } else {
          LRT_boot[b] <- max(0, 2 * (ll1b - ll0b))
        }
      }
    }
    LRT_boot <- LRT_boot[is.finite(LRT_boot)]

  } else {
    message(sprintf("Performing %d estimations under H0 and H1 using %d cores...", B, ncores))

    future::plan(future::multisession, workers = ncores)

    set_progress_handlers()
    pb <- if (isTRUE(progress)) progressr::progressor(along = seq_len(B)) else NULL

    xx <- foreach::foreach(
      b = seq_len(B),
      .combine = c,
      .options.future = list(
        packages   = "GeoModels",
        seed       = TRUE,
        stdout     = FALSE,
        conditions = "none"
      )
    ) %dofuture% {
      if (!is.null(pb)) pb(sprintf("b=%d", b))

      Zb  <- get_rep(sim, b)
      f0b <- suppressWarnings(try(fit_H0(Zb), silent = TRUE))
      f1b <- suppressWarnings(try(fit_H1(Zb), silent = TRUE))

      if (inherits(f0b, "try-error") || inherits(f1b, "try-error")) {
        as.numeric(NA)
      } else {
        ll0b <- if (!is.null(f0b$logCompLik)) f0b$logCompLik else f0b$logLik
        ll1b <- if (!is.null(f1b$logCompLik)) f1b$logCompLik else f1b$logLik
        ll0b <- as.numeric(ll0b)
        ll1b <- as.numeric(ll1b)

        if (!is.finite(ll0b) || !is.finite(ll1b)) {
          as.numeric(NA)
        } else {
          as.numeric(max(0, 2 * (ll1b - ll0b)))
        }
      }
    }

    LRT_boot <- xx[is.finite(xx)]
  }

  ## ====== RISULTATI ======
  failed <- B - length(LRT_boot)
  if (failed > 0) {
    warning(sprintf("%d/%d bootstrap replications failed (%.1f%%)",
                    failed, B, 100 * failed / B))
  }

  if (!length(LRT_boot)) {
    warning("No valid bootstrap replications; p-value set to NA")
    pval <- NA_real_
  } else {
    pval <- (1 + sum(LRT_boot >= LRT_obs)) / (length(LRT_boot) + 1)
  }

  ratio_hat <- fit1$anisopars$ratio
  angle_hat <- fit1$anisopars$angle

  cat("\n=== FINAL RESULTS ===\n")
  cat("Observed LRT statistic   =", round(LRT_obs, 5), "\n")
  cat("Valid bootstrap reps     =", length(LRT_boot), "/", B, "\n")
  cat("Estimated anisotropy: ratio =", round(ratio_hat, 4),
      "| angle =", round(angle_hat, 4), "\n")
  cat("p-value                  =", round(pval, 5), "\n\n")

  concl <- if (!is.na(pval) && pval < 0.05) {
    "Reject H0 => significant anisotropy detected"
  } else {
    "Do not reject H0 => data consistent with isotropy"
  }
  cat("Conclusion:", concl, "\n\n")

  invisible(list(
    statistic = LRT_obs,
    pvalue = pval,
    ratio_hat = ratio_hat,
    angle_hat = angle_hat,
    corrmodel = corrmodel,
    model = model_used,
    method = method,
    L = if (method == "TB") L else NULL,
    seed = seed,
    fit_H0 = fit0,
    fit_H1 = fit1,
    B_rep = LRT_boot
  ))
}
