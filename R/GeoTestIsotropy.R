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
                            L = 10000,
                            parallel = TRUE,
                            ncores = NULL,
                            progress = TRUE,
                            seed = NULL) {

  `%||%` <- function(x, y) if (!is.null(x)) x else y

  ## ====== CLEANUP: FUTURE PLAN AND PROGRESS HANDLERS ======
  future_plan_original <- future::plan()
  old_handlers <- progressr::handlers()

  cleanup_resources <- function() {
    tryCatch(
      future::plan(future_plan_original),
      error = function(e) invisible(NULL)
    )

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

  ## ====== INITIAL INPUT VALIDATION ======
  method <- match.arg(method)

  if (!is.matrix(coordx)) coordx <- as.matrix(coordx)
  if (!is.numeric(data)) data <- as.numeric(data)

  if (ncol(coordx) != 2L) {
    stop("GeoTestIsotropy currently requires two-dimensional coordinates", call. = FALSE)
  }
  if (any(!is.finite(coordx))) {
    stop("coordx contains non-finite values", call. = FALSE)
  }
  if (any(!is.finite(data))) {
    stop("data contains non-finite values", call. = FALSE)
  }
  if (nrow(coordx) != length(data)) {
    stop(
      "coordx rows (", nrow(coordx), ") must match data length (",
      length(data), ")",
      call. = FALSE
    )
  }

  if (!is.list(start) || is.null(names(start)) || any(names(start) == "")) {
    stop("start must be a named list", call. = FALSE)
  }
  if (is.null(fixed)) fixed <- list()
  if (!is.list(fixed) || (length(fixed) > 0L &&
                          (is.null(names(fixed)) || any(names(fixed) == "")))) {
    stop("fixed must be a named list", call. = FALSE)
  }

  duplicated_parameters <- intersect(names(start), names(fixed))
  if (length(duplicated_parameters) > 0L) {
    stop(
      "Parameters cannot be present in both start and fixed: ",
      paste(duplicated_parameters, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 1) {
    stop("B must be a positive integer", call. = FALSE)
  }
  B <- as.integer(B)

  if (!is.numeric(L) || length(L) != 1L || !is.finite(L) || L < 1) {
    stop("L must be a positive integer", call. = FALSE)
  }
  L <- as.integer(L)

  if (!is.logical(parallel) || length(parallel) != 1L) {
    stop("parallel must be logical", call. = FALSE)
  }
  if (!is.logical(progress) || length(progress) != 1L) {
    stop("progress must be logical", call. = FALSE)
  }

  NN <- nrow(coordx)

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
    stop(
      "Either specify both 'likelihood' and 'type', or leave both NULL for automatic selection",
      call. = FALSE
    )
  }

  if (identical(type, "Pairwise")) {
    if (is.null(neighb) || !is.numeric(neighb) || length(neighb) != 1L ||
        !is.finite(neighb) || neighb < 1) {
      stop("neighb must be a positive integer for pairwise likelihood", call. = FALSE)
    }
    neighb <- as.integer(round(neighb))
  } else {
    neighb <- NULL
  }

  if (!is.null(copula) && method != "cholesky") {
    stop("Copula simulation currently requires method = 'cholesky'", call. = FALSE)
  }

  ## ====== ROBUST CORE SELECTION ======
  coremax <- parallel::detectCores()
  if (is.na(coremax) || coremax <= 1L) {
    parallel <- FALSE
    ncores <- 1L
  } else {
    if (is.null(ncores)) ncores <- min(coremax - 1L, B)
    if (!is.numeric(ncores) || length(ncores) != 1L || !is.finite(ncores)) {
      stop("ncores must be a positive integer or NULL", call. = FALSE)
    }
    ncores <- as.integer(max(1L, min(ncores, B, coremax)))
    if (ncores == 1L) parallel <- FALSE
  }

  set_progress_handlers()

  corrmodel_used <- corrmodel
  model_used <- model

  ## ====== PARAMETER VALIDATION ======
  corr_names <- CorrParam(corrmodel_used)
  nuis_names <- NuisParam(model_used)
  expected <- sort(unique(c(corr_names, nuis_names)))
  provided <- sort(c(names(start), names(fixed)))

  if (!identical(expected, provided)) {
    missing_parameters <- setdiff(expected, provided)
    extra_parameters <- setdiff(provided, expected)
    msg <- character(0)

    if (length(missing_parameters) > 0L) {
      msg <- c(msg, paste0("missing: ", paste(missing_parameters, collapse = ", ")))
    }
    if (length(extra_parameters) > 0L) {
      msg <- c(msg, paste0("extra: ", paste(extra_parameters, collapse = ", ")))
    }

    stop(
      "Parameter names in start+fixed must match CorrParam(corrmodel)+NuisParam(model). ",
      paste(msg, collapse = " | "),
      call. = FALSE
    )
  }

  ## Bounds apply only to parameters in start.
  start_names <- sort(names(start))

  norm_bounds <- function(b, which = c("lower", "upper")) {
    which <- match.arg(which)
    if (is.null(b)) b <- list()

    b_names <- names(b)
    if (length(b) > 0L && is.null(b_names)) {
      stop(which, " must be a named list with names matching start", call. = FALSE)
    }
    if (length(setdiff(b_names, start_names)) > 0L) {
      stop(
        which, " contains names not in start: ",
        paste(setdiff(b_names, start_names), collapse = ", "),
        call. = FALSE
      )
    }

    out <- setNames(
      rep(if (which == "lower") -Inf else Inf, length(start_names)),
      start_names
    )
    if (length(b_names) > 0L) out[b_names] <- unlist(b, use.names = FALSE)
    as.list(out)
  }

  lower <- norm_bounds(lower, "lower")
  upper <- norm_bounds(upper, "upper")

  lower_vec <- unlist(lower, use.names = TRUE)
  upper_vec <- unlist(upper, use.names = TRUE)
  start_vec <- unlist(start, use.names = TRUE)[names(lower_vec)]

  if (any(lower_vec >= upper_vec)) {
    stop("Each lower bound must be smaller than its upper bound", call. = FALSE)
  }
  if (any(start_vec < lower_vec | start_vec > upper_vec)) {
    stop("Starting values must lie within the specified bounds", call. = FALSE)
  }

  ## ====== ANISOTROPY AND FIT HELPERS ======
  aniso_H0 <- list(angle = 0, ratio = 1)
  aniso_start <- list(angle = 0, ratio = 1.2)
  fallback_angles <- c(0, pi / 4, pi / 2, 3 * pi / 4)

  loglik_value <- function(fit_obj) {
    if (inherits(fit_obj, "try-error") || is.null(fit_obj)) return(NA_real_)
    val <- if (!is.null(fit_obj$logCompLik)) fit_obj$logCompLik else fit_obj$logLik
    val <- as.numeric(val)
    if (length(val) != 1L || !is.finite(val)) NA_real_ else val
  }

  strip_anisotropy <- function(par) {
    par <- as.list(par)
    par[setdiff(names(par), c("angle", "ratio"))]
  }

  fit_H0 <- function(Z, start_values = start) {
    GeoFit(
      data = Z,
      coordx = coordx,
      corrmodel = corrmodel_used,
      likelihood = likelihood,
      optimizer = optimizer,
      model = model_used,
      type = type,
      neighb = neighb,
      start = start_values,
      fixed = fixed,
      copula = copula,
      lower = lower,
      upper = upper,
      anisopars = aniso_H0
    )
  }

  fit_H1_once <- function(Z, start_values = start,
                          anisopars_values = aniso_start) {
    GeoFit(
      data = Z,
      coordx = coordx,
      corrmodel = corrmodel_used,
      likelihood = likelihood,
      optimizer = optimizer,
      model = model_used,
      type = type,
      neighb = neighb,
      start = start_values,
      fixed = fixed,
      copula = copula,
      lower = lower,
      upper = upper,
      anisopars = anisopars_values,
      est.aniso = c(TRUE, TRUE)
    )
  }

  same_pairs <- function(fit0_obj, fit1_obj) {
    if (!identical(type, "Pairwise")) return(TRUE)

    !is.null(fit0_obj$rowidx) && !is.null(fit0_obj$colidx) &&
      !is.null(fit1_obj$rowidx) && !is.null(fit1_obj$colidx) &&
      identical(fit0_obj$rowidx, fit1_obj$rowidx) &&
      identical(fit0_obj$colidx, fit1_obj$colidx)
  }

  fit_H1_adaptive <- function(Z, fit0_obj, start_values = start,
                              primary_anisopars = aniso_start) {
    ll0 <- loglik_value(fit0_obj)
    if (!is.finite(ll0)) {
      return(list(
        fit = NULL,
        used_fallback = FALSE,
        attempts = 0L,
        error = "H0 fit has a non-finite likelihood"
      ))
    }

    primary_angle <- as.numeric(primary_anisopars$angle %||% 0)
    primary_ratio <- as.numeric(primary_anisopars$ratio %||% 1.2)

    if (!is.finite(primary_angle)) primary_angle <- 0
    if (!is.finite(primary_ratio) || primary_ratio <= 1) primary_ratio <- 1.2

    primary_angle <- primary_angle %% pi
    primary_anisopars <- list(angle = primary_angle, ratio = primary_ratio)

    fits <- list(
      suppressWarnings(try(
        fit_H1_once(
          Z,
          start_values = start_values,
          anisopars_values = primary_anisopars
        ),
        silent = TRUE
      ))
    )
    attempts <- 1L

    ll_primary <- loglik_value(fits[[1L]])
    tol_ll <- 1e-8 * max(1, abs(ll0))
    need_fallback <- !is.finite(ll_primary) || ll_primary < ll0 - tol_ll

    if (need_fallback) {
      additional_angles <- fallback_angles[
        abs(fallback_angles - primary_angle) > sqrt(.Machine$double.eps)
      ]

      for (angle_start in additional_angles) {
        fits[[length(fits) + 1L]] <- suppressWarnings(try(
          fit_H1_once(
            Z,
            start_values = start_values,
            anisopars_values = list(angle = angle_start, ratio = 1.2)
          ),
          silent = TRUE
        ))
        attempts <- attempts + 1L
      }
    }

    ll_values <- vapply(fits, loglik_value, numeric(1))
    valid <- is.finite(ll_values)

    if (!any(valid)) {
      return(list(
        fit = NULL,
        used_fallback = need_fallback,
        attempts = attempts,
        error = "All anisotropic fits failed"
      ))
    }

    valid_indices <- which(valid)
    best_index <- valid_indices[which.max(ll_values[valid])]
    best_fit <- fits[[best_index]]

    if (!same_pairs(fit0_obj, best_fit)) {
      return(list(
        fit = NULL,
        used_fallback = need_fallback,
        attempts = attempts,
        error = "H0 and H1 do not use the same pair set"
      ))
    }

    list(
      fit = best_fit,
      used_fallback = need_fallback,
      attempts = attempts,
      error = NULL
    )
  }

  ## ====== OBSERVED FITS ======
  fit0 <- fit_H0(data)
  fit1_result <- fit_H1_adaptive(
    data,
    fit0_obj = fit0,
    start_values = start,
    primary_anisopars = aniso_start
  )

  if (is.null(fit1_result$fit)) {
    stop("Observed anisotropic fit failed: ", fit1_result$error, call. = FALSE)
  }

  fit1 <- fit1_result$fit
  LL0 <- loglik_value(fit0)
  LL1 <- loglik_value(fit1)

  if (!is.finite(LL0) || !is.finite(LL1)) {
    stop("Observed H0/H1 fit returned non-finite log-likelihood values", call. = FALSE)
  }

  if (!same_pairs(fit0, fit1)) {
    stop(
      "Observed H0 and H1 fits must use the same pairs for the composite-likelihood ratio",
      call. = FALSE
    )
  }

  if (fit1_result$used_fallback) {
    message(
      "Adaptive angle fallback used for the observed H1 fit (",
      fit1_result$attempts, " attempts)"
    )
  }

  tol_obs <- 1e-8 * max(1, abs(LL0))
  if (LL1 < LL0 - tol_obs) {
    warning(
      "The best anisotropic fit has a lower likelihood than the isotropic fit; ",
      "the observed LRT is set to zero",
      call. = FALSE
    )
  }
  LRT_obs <- max(0, 2 * (LL1 - LL0))

  ## Use observed estimates as bootstrap starting values.  Angle and ratio
  ## remain in anisopars and are therefore removed from start_H1_boot.
  start_H0_boot <- strip_anisotropy(fit0$param)
  start_H1_boot <- strip_anisotropy(fit1$param)

  primary_aniso_boot <- list(
    angle = as.numeric(fit1$anisopars$angle %||% 0),
    ratio = as.numeric(fit1$anisopars$ratio %||% 1.2)
  )

  ## ====== SIMULATIONS UNDER H0 ======
  par0 <- c(fit0$param, fit0$fixed)
  param_H0 <- as.list(par0[expected])

  if (is.null(copula)) {
    if (method == "cholesky") {
      sim <- GeoSim(
        coordx = coordx,
        corrmodel = corrmodel_used,
        param = param_H0,
        model = model_used,
        nrep = B,
        anisopars = aniso_H0,
        progress = FALSE
      )
    } else {
      plan_before_sim <- future::plan()

      sim <- GeoSimapprox(
        coordx = coordx,
        corrmodel = corrmodel_used,
        model = model_used,
        method = "TB",
        L = L,
        param = param_H0,
        nrep = B,
        anisopars = aniso_H0,
        progress = FALSE,
        parallel = parallel,
        ncores = ncores
      )

      try(future::plan(plan_before_sim), silent = TRUE)
    }
  } else {
    sim <- GeoSimCopula(
      coordx = coordx,
      corrmodel = corrmodel_used,
      method = "Cholesky",
      param = param_H0,
      model = model_used,
      copula = copula,
      nrep = B,
      anisopars = aniso_H0
    )
  }

  get_rep <- function(simobj, b) {
    X <- simobj$data
    if (is.list(X)) X[[b]] else X[, b]
  }

  bootstrap_one <- function(Zb) {
    f0b <- suppressWarnings(try(
      fit_H0(Zb, start_values = start_H0_boot),
      silent = TRUE
    ))

    if (inherits(f0b, "try-error") || !is.finite(loglik_value(f0b))) {
      return(list(
        lrt = NA_real_,
        used_fallback = FALSE,
        attempts = 0L,
        error = "H0 fit failed"
      ))
    }

    f1_result <- fit_H1_adaptive(
      Zb,
      fit0_obj = f0b,
      start_values = start_H1_boot,
      primary_anisopars = primary_aniso_boot
    )

    if (is.null(f1_result$fit)) {
      return(list(
        lrt = NA_real_,
        used_fallback = f1_result$used_fallback,
        attempts = f1_result$attempts,
        error = f1_result$error
      ))
    }

    ll0b <- loglik_value(f0b)
    ll1b <- loglik_value(f1_result$fit)

    if (!is.finite(ll0b) || !is.finite(ll1b)) {
      return(list(
        lrt = NA_real_,
        used_fallback = f1_result$used_fallback,
        attempts = f1_result$attempts,
        error = "Non-finite bootstrap likelihood"
      ))
    }

    list(
      lrt = max(0, 2 * (ll1b - ll0b)),
      used_fallback = f1_result$used_fallback,
      attempts = f1_result$attempts,
      error = NULL
    )
  }

  ## ====== BOOTSTRAP ======
  if (!parallel) {
    message(sprintf("Performing %d bootstrap refits sequentially...", B))

    run_sequential <- function() {
      out <- vector("list", B)
      pb <- if (isTRUE(progress)) progressr::progressor(along = seq_len(B)) else NULL

      for (b in seq_len(B)) {
        out[[b]] <- bootstrap_one(get_rep(sim, b))
        if (!is.null(pb)) pb(sprintf("b=%d", b))
      }
      out
    }

    bootstrap_results <- if (isTRUE(progress)) {
      progressr::with_progress(run_sequential())
    } else {
      run_sequential()
    }

  } else {
    message(sprintf(
      "Performing %d bootstrap refits using %d cores...",
      B, ncores
    ))

    future::plan(future::multisession, workers = ncores)

    temp_dir <- tempdir()
    tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
    data_files <- vapply(seq_len(B), function(b) {
      path <- file.path(temp_dir, sprintf("isotropy_%s_%05d.rds", tag, b))
      saveRDS(get_rep(sim, b), path, compress = FALSE)
      path
    }, character(1))

    on.exit(unlink(data_files), add = TRUE)
    rm(sim)
    gc(verbose = FALSE, full = TRUE)

    run_parallel <- function() {
      pb <- if (isTRUE(progress)) progressr::progressor(along = seq_len(B)) else NULL

      future.apply::future_lapply(
        seq_len(B),
        function(b) {
          out <- bootstrap_one(readRDS(data_files[[b]]))
          if (!is.null(pb)) pb(sprintf("b=%d", b))
          out
        },
        future.seed = TRUE,
        future.packages = "GeoModels"
      )
    }

    bootstrap_results <- if (isTRUE(progress)) {
      progressr::with_progress(run_parallel())
    } else {
      run_parallel()
    }
  }

  LRT_all <- vapply(bootstrap_results, function(x) x$lrt, numeric(1))
  valid <- is.finite(LRT_all)
  LRT_boot <- LRT_all[valid]
  errors <- vapply(
    bootstrap_results,
    function(x) x$error %||% "",
    character(1)
  )
  fallback_used <- vapply(
    bootstrap_results,
    function(x) isTRUE(x$used_fallback),
    logical(1)
  )
  attempts <- vapply(
    bootstrap_results,
    function(x) as.integer(x$attempts %||% 0L),
    integer(1)
  )

  ## ====== RESULTS ======
  n_valid <- length(LRT_boot)
  failed <- B - n_valid
  success_rate <- n_valid / B

  if (failed > 0L) {
    warning(
      sprintf(
        "%d/%d bootstrap replications failed (%.1f%%)",
        failed, B, 100 * failed / B
      ),
      call. = FALSE
    )
  }

  if (success_rate < 0.8) {
    warning(
      sprintf(
        "Only %.1f%% of bootstrap replications were successful; the p-value may be unreliable",
        100 * success_rate
      ),
      call. = FALSE
    )
  }

  if (n_valid < 20L) {
    warning(
      "Fewer than 20 valid bootstrap replications are available; the p-value is very imprecise",
      call. = FALSE
    )
  }

  if (n_valid == 0L) {
    pval <- NA_real_
  } else {
    pval <- (1 + sum(LRT_boot >= LRT_obs)) / (n_valid + 1)
  }

  ratio_hat <- as.numeric(fit1$anisopars$ratio)
  angle_hat <- as.numeric(fit1$anisopars$angle)

  cat("\n=== FINAL RESULTS ===\n")
  cat("Observed LRT statistic   =", round(LRT_obs, 5), "\n")
  cat("Valid bootstrap reps     =", n_valid, "/", B, "\n")
  cat("Adaptive fallbacks       =", sum(fallback_used), "/", B, "\n")
  cat(
    "Estimated anisotropy: ratio =", round(ratio_hat, 4),
    "| angle =", round(angle_hat, 4), "\n"
  )
  cat("p-value                  =", round(pval, 5), "\n\n")

  conclusion <- if (!is.na(pval) && pval < 0.05) {
    "Reject H0 => significant anisotropy detected"
  } else {
    "Do not reject H0 => data consistent with isotropy"
  }
  cat("Conclusion:", conclusion, "\n\n")

  failure_table <- data.frame(
    iteration = which(!valid),
    error = errors[!valid],
    stringsAsFactors = FALSE
  )

  invisible(list(
    statistic = LRT_obs,
    pvalue = pval,
    ratio_hat = ratio_hat,
    angle_hat = angle_hat,
    corrmodel = corrmodel_used,
    model = model_used,
    likelihood = likelihood,
    type = type,
    method = method,
    L = if (method == "TB") L else NULL,
    seed = seed,
    fit_H0 = fit0,
    fit_H1 = fit1,
    observed_fallback_used = fit1_result$used_fallback,
    observed_H1_attempts = fit1_result$attempts,
    B_rep = LRT_boot,
    bootstrap_B = B,
    bootstrap_successful = n_valid,
    bootstrap_success_rate = success_rate,
    bootstrap_fallbacks = sum(fallback_used),
    bootstrap_attempts = attempts,
    bootstrap_failures = failure_table
  ))
}