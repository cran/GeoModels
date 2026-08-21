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
                            n = 1,
                            angle_starts = c(pi / 4, 3 * pi / 4),
                            batch_size = NULL,
                            max_attempts = NULL,
                            parallel = FALSE,
                            ncores = 6,
                            progress = TRUE,
                            seed = NULL) {

  `%||%` <- function(x, y) if (!is.null(x)) x else y

  ## ====== preserve runtime settings and optional RNG seed ======
  ## Do not use .GeoTestPreserveRuntime() here.  That helper inspects and
  ## restores global progressr handlers, which is unsafe when this function
  ## is called inside tryCatch()/withCallingHandlers() and can trigger
  ## "should not be called with handlers on the stack" even when
  ## progress = FALSE.  GeoTestIsotropy manages its own progress state below.
  old_future_plan <- if (requireNamespace("future", quietly = TRUE)) future::plan() else NULL
  on.exit({
    if (!is.null(old_future_plan)) try(future::plan(old_future_plan), silent = TRUE)
    gc(verbose = FALSE, full = TRUE)
  }, add = TRUE)
  restore_seed <- .GeoTestSetSeed(seed)
  on.exit(restore_seed(), add = TRUE)

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

  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 1 ||
      abs(B - round(B)) > sqrt(.Machine$double.eps)) {
    stop("B must be a positive integer", call. = FALSE)
  }
  B <- as.integer(round(B))

  if (!is.numeric(L) || length(L) != 1L || !is.finite(L) || L < 1 ||
      abs(L - round(L)) > sqrt(.Machine$double.eps)) {
    stop("L must be a positive integer", call. = FALSE)
  }
  L <- as.integer(round(L))

  if (!is.logical(parallel) || length(parallel) != 1L) {
    stop("parallel must be logical", call. = FALSE)
  }
  if (!is.logical(progress) || length(progress) != 1L) {
    stop("progress must be logical", call. = FALSE)
  }

  if (!is.numeric(angle_starts) || length(angle_starts) < 1L ||
      any(!is.finite(angle_starts))) {
    stop("angle_starts must contain at least one finite numeric angle", call. = FALSE)
  }
  angle_starts <- as.numeric(angle_starts) %% pi
  angle_starts <- angle_starts[!duplicated(round(angle_starts, 12L))]

  NN <- nrow(coordx)

  model_canonical <- .GeoSimulationCanonicalModel(model)
  n_used <- .GeoSimulationEffectiveN(model, model_canonical, n)
  if (model_canonical %in% c("Binomial", "BinomialLogistic", "BinomialNeg", "BinomialNegZINB")) {
    n_used <- .GeoValidateIntegerVector(n_used, "n", nobs = NN, minimum = 1L)
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

  if (!is.null(copula) && identical(method, "TB")) {
    .GeoValidateTBDistance("Eucl", context = "GeoTestIsotropy")
    .GeoValidateTBCorrelationModel(
      .GeoSimulationCanonicalCorrModel(corrmodel),
      bivariate = FALSE
    )
  }

  ## ====== ROBUST CORE SELECTION ======
  core_cfg <- .GeoTestCoreConfig(parallel, ncores, B, strict_positive = FALSE)
  parallel <- core_cfg$parallel
  ncores <- core_cfg$ncores

  ## Use the same progress-bar machinery adopted by the other parallel
  ## GeoModels functions (e.g. GeoNeighbSelect and GeoKrigloc).  In
  ## particular, install one temporary txtprogressbar handler and keep a
  ## single progressor owned by the main R process.  Do not use
  ## .GeoTestConfigureProgress() here: changing global progressr handlers
  ## from inside nested calling handlers can trigger
  ## "should not be called with handlers on the stack".
  use_progressr <- isTRUE(progress) && requireNamespace("progressr", quietly = TRUE)
  if (isTRUE(progress) && !use_progressr) {
    warning(
      "progress=TRUE but progressr is not available; progress bar disabled.",
      call. = FALSE
    )
  }
  if (use_progressr) {
    restore_progress_handlers <- .GeoProgressHandlersPush(
      handler = "txtprogressbar", global = TRUE
    )
    on.exit(restore_progress_handlers(), add = TRUE)
  }

  if (is.null(batch_size)) {
    batch_size <- if (parallel) max(2L, 2L * ncores) else min(10L, B)
  }
  if (!is.numeric(batch_size) || length(batch_size) != 1L ||
      !is.finite(batch_size) || batch_size < 1 ||
      abs(batch_size - round(batch_size)) > sqrt(.Machine$double.eps)) {
    stop("batch_size must be NULL or a positive integer", call. = FALSE)
  }
  batch_size <- min(B, as.integer(round(batch_size)))

  if (is.null(max_attempts)) {
    max_attempts <- max(B + 20L, as.integer(ceiling(1.5 * B)))
  }
  if (!is.numeric(max_attempts) || length(max_attempts) != 1L ||
      !is.finite(max_attempts) || max_attempts < B ||
      abs(max_attempts - round(max_attempts)) > sqrt(.Machine$double.eps)) {
    stop("max_attempts must be NULL or an integer greater than or equal to B", call. = FALSE)
  }
  max_attempts <- as.integer(round(max_attempts))

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
    .GeoTestNormalizeBounds(b, which, start_names, require_list = FALSE)
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
  aniso_ratio_start <- 1.2

  loglik_value <- .GeoTestLogLikValue

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
      n = n_used,
      anisopars = aniso_H0
    )
  }

  fit_H1_once <- function(Z, start_values = start,
                          anisopars_values = list(angle = 0, ratio = 1.2)) {
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
      n = n_used,
      anisopars = anisopars_values,
      est.aniso = c(TRUE, TRUE)
    )
  }

  same_pairs <- function(fit0_obj, fit1_obj) {
    .GeoTestSamePairs(fit0_obj, fit1_obj, type, require_indices = TRUE)
  }

  fit_H1_multistart <- function(Z, fit0_obj, start_values = start,
                                angles = angle_starts,
                                ratio_start = aniso_ratio_start) {
    ll0 <- loglik_value(fit0_obj)
    if (!is.finite(ll0)) {
      return(list(
        fit = NULL,
        attempts = 0L,
        valid_attempts = 0L,
        error = "H0 fit has a non-finite likelihood"
      ))
    }

    ratio_start <- as.numeric(ratio_start)[1L]
    if (!is.finite(ratio_start) || ratio_start <= 1) ratio_start <- 1.2

    angles <- as.numeric(angles) %% pi
    angles <- angles[is.finite(angles)]
    angles <- angles[!duplicated(round(angles, 12L))]
    if (length(angles) < 1L) angles <- 0

    fits <- lapply(angles, function(angle_start) {
      suppressWarnings(try(
        fit_H1_once(
          Z,
          start_values = start_values,
          anisopars_values = list(angle = angle_start, ratio = ratio_start)
        ),
        silent = TRUE
      ))
    })

    ll_values <- vapply(fits, loglik_value, numeric(1))
    valid <- is.finite(ll_values)

    if (!any(valid)) {
      return(list(
        fit = NULL,
        attempts = length(fits),
        valid_attempts = 0L,
        error = "All anisotropic multi-start fits failed"
      ))
    }

    valid_indices <- which(valid)
    best_index <- valid_indices[which.max(ll_values[valid])]
    best_fit <- fits[[best_index]]

    if (!same_pairs(fit0_obj, best_fit)) {
      return(list(
        fit = NULL,
        attempts = length(fits),
        valid_attempts = sum(valid),
        error = "H0 and H1 do not use the same pair set"
      ))
    }

    list(
      fit = best_fit,
      attempts = length(fits),
      valid_attempts = sum(valid),
      starting_angle = angles[best_index],
      error = NULL
    )
  }

  ## ====== OBSERVED FITS ======
  fit0 <- fit_H0(data)
  if (!is.finite(loglik_value(fit0))) {
    stop("Observed isotropic fit returned a non-finite log-likelihood", call. = FALSE)
  }

  fit1_result <- fit_H1_multistart(
    data,
    fit0_obj = fit0,
    start_values = start,
    angles = angle_starts,
    ratio_start = aniso_ratio_start
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

  tol_obs <- 1e-8 * max(1, abs(LL0))
  observed_boundary_selected <- LL1 < LL0 - tol_obs
  LRT_obs <- if (observed_boundary_selected) 0 else max(0, 2 * (LL1 - LL0))

  ratio_fit <- as.numeric(fit1$anisopars$ratio %||% NA_real_)
  angle_fit <- as.numeric(fit1$anisopars$angle %||% NA_real_)
  ratio_hat <- if (observed_boundary_selected) 1 else ratio_fit
  angle_hat <- if (observed_boundary_selected) NA_real_ else angle_fit

  if (!is.finite(ratio_hat) || ratio_hat < 1) {
    stop("Observed anisotropy-ratio estimate is invalid", call. = FALSE)
  }

  ## Use observed estimates as bootstrap starting values. Angle and ratio
  ## remain in anisopars and are therefore removed from start_H1_boot.
  start_H0_boot <- strip_anisotropy(fit0$param)
  start_H1_boot <- strip_anisotropy(fit1$param)

  ## Use the same H1 multi-start angles for the observed fit and every
  ## bootstrap refit. Keeping this set fixed makes the numerical search
  ## comparable across the observed and bootstrap samples and controls cost.
  bootstrap_angle_starts <- angle_starts
  bootstrap_ratio_start <- if (is.finite(ratio_fit) && ratio_fit > 1) ratio_fit else 1.2

  ## ====== NULL MODEL FOR PARAMETRIC BOOTSTRAP ======
  par0 <- c(fit0$param, fit0$fixed)
  param_H0 <- as.list(par0[expected])

  if (is.null(copula)) {
    latent_check <- suppressWarnings(try(
      .GeoValidateLatentIntegerParameters(model_canonical, param_H0),
      silent = TRUE
    ))
    if (inherits(latent_check, "try-error")) {
      cond <- attr(latent_check, "condition")
      msg <- if (!is.null(cond)) conditionMessage(cond) else as.character(latent_check)
      extra <- if (identical(model_canonical, "Gamma")) {
        paste0(
          " For the direct Gamma latent-squares model used in the paper, ",
          "fix 'shape' to a positive integer (put it in 'fixed') so the fitted ",
          "null model can be simulated exactly by the latent-Gaussian construction."
        )
      } else {
        ""
      }
      stop(
        "The fitted null model cannot be simulated by the direct latent-Gaussian bootstrap: ",
        msg, extra,
        call. = FALSE
      )
    }
  }

  simulate_null <- function(nrep, workers = ncores) {
    if (is.null(copula)) {
      if (method == "cholesky") {
        ## GeoSim() currently prints an informational line for Cholesky
        ## simulation even when progress=FALSE. Capture that output here so
        ## GeoTestIsotropy() can present a single bootstrap progress bar.
        sim_out <- NULL
        invisible(utils::capture.output(
          sim_out <- GeoSim(
            coordx = coordx,
            corrmodel = corrmodel_used,
            param = param_H0,
            model = model_used,
            n = n_used,
            nrep = nrep,
            anisopars = aniso_H0,
            progress = FALSE
          ),
          type = "output"
        ))
        sim_out
      } else {
        plan_before_sim <- future::plan()
        on.exit(try(future::plan(plan_before_sim), silent = TRUE), add = TRUE)
        GeoSimapprox(
          coordx = coordx,
          corrmodel = corrmodel_used,
          model = model_used,
          method = "TB",
          L = L,
          n = n_used,
          param = param_H0,
          nrep = nrep,
          anisopars = aniso_H0,
          progress = FALSE,
          parallel = parallel && workers > 1L,
          ncores = workers
        )
      }
    } else {
      GeoSimCopula(
        coordx = coordx,
        corrmodel = corrmodel_used,
        method = method,
        param = param_H0,
        model = model_used,
        n = n_used,
        copula = copula,
        nrep = nrep,
        anisopars = aniso_H0,
        L = L,
        parallel = parallel && workers > 1L,
        ncores = workers,
        progress = FALSE
      )
    }
  }

  bootstrap_one <- function(Zb) {
    f0b <- suppressWarnings(try(
      fit_H0(Zb, start_values = start_H0_boot),
      silent = TRUE
    ))

    if (inherits(f0b, "try-error") || !is.finite(loglik_value(f0b))) {
      return(list(
        valid = FALSE,
        lrt = NA_real_,
        ratio = NA_real_,
        angle = NA_real_,
        boundary_selected = FALSE,
        h1_attempts = 0L,
        h1_valid_attempts = 0L,
        error = "H0 fit failed"
      ))
    }

    f1_result <- fit_H1_multistart(
      Zb,
      fit0_obj = f0b,
      start_values = start_H1_boot,
      angles = bootstrap_angle_starts,
      ratio_start = bootstrap_ratio_start
    )

    if (is.null(f1_result$fit)) {
      return(list(
        valid = FALSE,
        lrt = NA_real_,
        ratio = NA_real_,
        angle = NA_real_,
        boundary_selected = FALSE,
        h1_attempts = f1_result$attempts,
        h1_valid_attempts = f1_result$valid_attempts,
        error = f1_result$error
      ))
    }

    ll0b <- loglik_value(f0b)
    ll1b <- loglik_value(f1_result$fit)

    if (!is.finite(ll0b) || !is.finite(ll1b)) {
      return(list(
        valid = FALSE,
        lrt = NA_real_,
        ratio = NA_real_,
        angle = NA_real_,
        boundary_selected = FALSE,
        h1_attempts = f1_result$attempts,
        h1_valid_attempts = f1_result$valid_attempts,
        error = "Non-finite bootstrap likelihood"
      ))
    }

    tol_b <- 1e-8 * max(1, abs(ll0b))
    boundary_selected <- ll1b < ll0b - tol_b
    ratio_b <- if (boundary_selected) {
      1
    } else {
      as.numeric(f1_result$fit$anisopars$ratio %||% NA_real_)
    }
    angle_b <- if (boundary_selected) {
      NA_real_
    } else {
      as.numeric(f1_result$fit$anisopars$angle %||% NA_real_)
    }

    if (!is.finite(ratio_b) || ratio_b < 1) {
      return(list(
        valid = FALSE,
        lrt = NA_real_,
        ratio = NA_real_,
        angle = NA_real_,
        boundary_selected = boundary_selected,
        h1_attempts = f1_result$attempts,
        h1_valid_attempts = f1_result$valid_attempts,
        error = "Invalid bootstrap anisotropy-ratio estimate"
      ))
    }

    list(
      valid = TRUE,
      lrt = if (boundary_selected) 0 else max(0, 2 * (ll1b - ll0b)),
      ratio = ratio_b,
      angle = angle_b,
      boundary_selected = boundary_selected,
      h1_attempts = f1_result$attempts,
      h1_valid_attempts = f1_result$valid_attempts,
      error = NULL
    )
  }

  process_batch <- function(sim, batch_n, batch_id) {
    ## Keep one package-standard progress bar for the whole bootstrap.  In the
    ## sequential branch the master updates it immediately after each valid
    ## refit.  In the parallel branch the master polls individual futures and
    ## updates it as each valid refit actually completes, avoiding buffered or
    ## batchwise jumps in the display.
    if (!parallel) {
      out <- vector("list", batch_n)
      for (b in seq_len(batch_n)) {
        out[[b]] <- bootstrap_one(.GeoTestGetRep(sim, b))
        if (!is.null(pb) && isTRUE(out[[b]]$valid) &&
            is.finite(out[[b]]$lrt) && is.finite(out[[b]]$ratio)) {
          pb(sprintf("b=%d", b))
        }
      }
      out
    } else {
      data_files <- .GeoTestSaveReplicates(
        sim, batch_n, sprintf("isotropy_batch%04d", batch_id)
      )
      on.exit(unlink(data_files), add = TRUE)
      rm(sim)
      gc(verbose = FALSE, full = TRUE)
      .GeoTestParallelFileApply(
        data_files = data_files,
        FUN = bootstrap_one,
        ncores = min(ncores, batch_n),
        progressor = pb,
        future_packages = "GeoModels"
      )
    }
  }

  ## ====== BOOTSTRAP: collect exactly B valid replications ======
  ## One progressor is created for the complete bootstrap, using the same
  ## txtprogressbar handler as the other parallel GeoModels functions.  The
  ## bar counts valid replications and is advanced by the main R process when
  ## each refit has genuinely completed.
  if (use_progressr) {
    pb <- progressr::progressor(
      steps = B, on_exit = FALSE, auto_finish = FALSE
    )
    on.exit(pb(type = "finish"), add = TRUE, after = FALSE)
  } else {
    pb <- NULL
  }

  run_bootstrap <- function() {
    valid_results <- vector("list", 0L)
    failure_records <- vector("list", 0L)
    total_attempts <- 0L
    batch_id <- 0L

    while (length(valid_results) < B && total_attempts < max_attempts) {
      outstanding <- B - length(valid_results)
      remaining_attempts <- max_attempts - total_attempts
      batch_n <- min(batch_size, outstanding, remaining_attempts)
      if (batch_n < 1L) break

      batch_id <- batch_id + 1L
      workers_batch <- min(ncores, batch_n)
      sim <- simulate_null(batch_n, workers = workers_batch)
      batch_results <- process_batch(sim, batch_n, batch_id)
      if (exists("sim", inherits = FALSE)) rm(sim)
      gc(verbose = FALSE, full = TRUE)

      for (j in seq_along(batch_results)) {
        total_attempts <- total_attempts + 1L
        ans <- batch_results[[j]]

        if (isTRUE(ans$valid) && is.finite(ans$lrt) && is.finite(ans$ratio)) {
          valid_results[[length(valid_results) + 1L]] <- ans
        } else {
          failure_records[[length(failure_records) + 1L]] <- data.frame(
            attempt = total_attempts,
            error = ans$error %||% "Unknown bootstrap failure",
            h1_attempts = as.integer(ans$h1_attempts %||% 0L),
            h1_valid_attempts = as.integer(ans$h1_valid_attempts %||% 0L),
            stringsAsFactors = FALSE
          )
        }
      }

    }

    list(
      valid_results = valid_results,
      failure_records = failure_records,
      total_attempts = total_attempts,
      batch_id = batch_id
    )
  }

  bootstrap_state <- run_bootstrap()

  valid_results <- bootstrap_state$valid_results
  failure_records <- bootstrap_state$failure_records
  total_attempts <- bootstrap_state$total_attempts
  batch_id <- bootstrap_state$batch_id

  n_valid <- length(valid_results)
  if (n_valid < B) {
    stop(
      "GeoTestIsotropy obtained only ", n_valid, " valid bootstrap replications after ",
      total_attempts, " attempts (requested B = ", B, "). Increase 'max_attempts', ",
      "improve starting values/bounds, or inspect optimization failures.",
      call. = FALSE
    )
  }

  LRT_boot <- vapply(valid_results, function(x) x$lrt, numeric(1))
  ratio_boot <- vapply(valid_results, function(x) x$ratio, numeric(1))
  angle_boot <- vapply(valid_results, function(x) x$angle, numeric(1))
  boundary_boot <- vapply(
    valid_results, function(x) isTRUE(x$boundary_selected), logical(1)
  )
  h1_attempts <- vapply(
    valid_results, function(x) as.integer(x$h1_attempts %||% 0L), integer(1)
  )
  h1_valid_attempts <- vapply(
    valid_results, function(x) as.integer(x$h1_valid_attempts %||% 0L), integer(1)
  )

  pval_lrt <- (1 + sum(LRT_boot >= LRT_obs)) / (B + 1)
  pval_ratio <- (1 + sum(ratio_boot >= ratio_hat)) / (B + 1)

  if (!is.null(pb)) pb(type = "finish")

  if (length(failure_records) > 0L) {
    failure_table <- do.call(rbind, failure_records)
    rownames(failure_table) <- NULL
  } else {
    failure_table <- data.frame(
      attempt = integer(0),
      error = character(0),
      h1_attempts = integer(0),
      h1_valid_attempts = integer(0),
      stringsAsFactors = FALSE
    )
  }

  ## ====== RESULTS ======
  cat("\n=== FINAL RESULTS ===\n")
  cat("Observed WCLRT/LRT statistic =", round(LRT_obs, 5), "\n")
  cat("Observed anisotropy ratio    =", round(ratio_hat, 5), "\n")
  cat("Valid bootstrap reps         =", B, "\n")
  cat("Total bootstrap attempts     =", total_attempts, "\n")
  cat("Failed attempts replaced     =", nrow(failure_table), "\n")
  cat("H1 multi-start angles        =", length(angle_starts), "\n")
  if (is.finite(angle_hat)) {
    cat("Estimated anisotropy angle   =", round(angle_hat, 5), "\n")
  } else {
    cat("Estimated anisotropy angle   = NA (boundary solution)\n")
  }
  cat("MWCLT p-value (ratio)        =", round(pval_ratio, 5), "\n")
  cat("WCLRT p-value                =", round(pval_lrt, 5), "\n\n")

  conclusion_ratio <- if (pval_ratio < 0.05) {
    "Reject H0 => significant anisotropy detected"
  } else {
    "Do not reject H0 => data consistent with isotropy"
  }
  conclusion_lrt <- if (pval_lrt < 0.05) {
    "Reject H0 => significant anisotropy detected"
  } else {
    "Do not reject H0 => data consistent with isotropy"
  }

  cat("MWCLT conclusion:", conclusion_ratio, "\n")
  cat("WCLRT conclusion:", conclusion_lrt, "\n\n")

  invisible(list(
    ## Backward-compatible primary result: WCLRT/LRT.
    statistic = LRT_obs,
    pvalue = pval_lrt,
    ## Explicit results for both tests described in the paper.
    statistic_lrt = LRT_obs,
    pvalue_lrt = pval_lrt,
    statistic_ratio = ratio_hat,
    pvalue_ratio = pval_ratio,
    pvalues = c(MWCLT = pval_ratio, WCLRT = pval_lrt),
    ratio_hat = ratio_hat,
    angle_hat = angle_hat,
    corrmodel = corrmodel_used,
    model = model_used,
    likelihood = likelihood,
    type = type,
    method = method,
    parallel = parallel,
    ncores = if (parallel) ncores else 1L,
    L = if (method == "TB") L else NULL,
    n = n_used,
    seed = seed,
    fit_H0 = fit0,
    fit_H1 = fit1,
    observed_boundary_selected = observed_boundary_selected,
    observed_multistart_angles = angle_starts,
    observed_H1_attempts = fit1_result$attempts,
    observed_H1_valid_attempts = fit1_result$valid_attempts,
    ## Backward-compatible alias for the LRT bootstrap statistics.
    B_rep = LRT_boot,
    B_rep_lrt = LRT_boot,
    B_rep_ratio = ratio_boot,
    B_rep_angle = angle_boot,
    bootstrap_B = B,
    bootstrap_successful = B,
    bootstrap_success_rate = 1,
    bootstrap_total_attempts = total_attempts,
    bootstrap_failed_attempts = nrow(failure_table),
    bootstrap_boundary_solutions = sum(boundary_boot),
    bootstrap_H1_attempts = h1_attempts,
    bootstrap_H1_valid_attempts = h1_valid_attempts,
    bootstrap_failures = failure_table,
    batch_size = batch_size,
    max_attempts = max_attempts,
    ## Retained for compatibility with objects produced by earlier versions.
    observed_fallback_used = FALSE,
    bootstrap_fallbacks = 0L,
    bootstrap_attempts = h1_attempts
  ))
}
