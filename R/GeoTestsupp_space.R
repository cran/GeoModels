########################################
## Parametric bootstrap test for spatial support
########################################
GeoTestsupp_space <- function(data, coordx,
                              start, fixed,
                              model = "Gaussian",
                              h0 = NULL,
                              optimizer = "bobyqa",
                              lower = NULL, upper = NULL,
                              neighb = 5,
                              B = 10000,
                              likelihood = NULL,
                              type = NULL,
                              method = c("cholesky", "TB"),
                              L = 10000,
                              parallel = FALSE,
                              ncores = 6,
                              progress = TRUE,
                              seed = NULL) {

  ## ====== preserve runtime settings and optional RNG seed ======
  cleanup_resources <- .GeoTestPreserveRuntime(gc_on_exit = FALSE)
  on.exit(cleanup_resources(), add = TRUE)
  restore_seed <- .GeoTestSetSeed(seed)
  on.exit(restore_seed(), add = TRUE)

  ## ====== input validation ======
  method <- match.arg(method)

  if (!is.matrix(coordx)) coordx <- as.matrix(coordx)
  if (!is.numeric(data)) data <- as.numeric(data)

  if (!is.numeric(coordx) || any(!is.finite(coordx))) {
    stop("coordx must contain only finite numeric values", call. = FALSE)
  }
  if (any(!is.finite(data))) {
    stop("data must contain only finite numeric values", call. = FALSE)
  }
  if (nrow(coordx) != length(data)) {
    stop(
      "coordx rows (", nrow(coordx), ") must match data length (",
      length(data), ")",
      call. = FALSE
    )
  }
  if (nrow(coordx) < 2L) {
    stop("At least two spatial observations are required", call. = FALSE)
  }

  if (!is.list(start) || length(start) < 1L ||
      is.null(names(start)) || any(names(start) == "")) {
    stop("start must be a non-empty named list", call. = FALSE)
  }

  if (is.null(fixed)) fixed <- list()
  if (!is.list(fixed)) {
    stop("fixed must be NULL or a named list", call. = FALSE)
  }
  if (length(fixed) > 0L &&
      (is.null(names(fixed)) || any(names(fixed) == ""))) {
    stop("fixed must be NULL or a named list", call. = FALSE)
  }

  if (anyDuplicated(names(start)) || anyDuplicated(names(fixed))) {
    stop("start and fixed cannot contain duplicated parameter names", call. = FALSE)
  }

  overlap <- intersect(names(start), names(fixed))
  if (length(overlap) > 0L) {
    stop(
      "Parameters cannot be present in both start and fixed: ",
      paste(overlap, collapse = ", "),
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

  if (!is.logical(progress) || length(progress) != 1L) {
    stop("progress must be logical (TRUE or FALSE)", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1L) {
    stop("parallel must be logical (TRUE or FALSE)", call. = FALSE)
  }
  if (!is.null(neighb) &&
      (!is.numeric(neighb) || length(neighb) != 1L ||
       !is.finite(neighb) || neighb < 1)) {
    stop("neighb must be a positive integer or NULL", call. = FALSE)
  }
  if (!is.null(neighb)) neighb <- as.integer(round(neighb))

  if (!is.null(h0) &&
      (!is.numeric(h0) || length(h0) != 1L || !is.finite(h0) || h0 <= 0)) {
    stop("h0 must be NULL or a positive finite number", call. = FALSE)
  }
  if (xor(is.null(likelihood), is.null(type))) {
    stop(
      "Either specify both likelihood and type, or leave both NULL for automatic selection.",
      call. = FALSE
    )
  }

  ## ====== constants ======
  max_scale_cap <- 1e6
  min_box_width <- 1e-3
  corrmodel <- "GenWend"
  eps <- 1e-6
  NN <- nrow(coordx)

  ## ====== progress handlers ======
  use_progressr <- .GeoTestConfigureProgress(progress)

  ## ====== automatic likelihood selection ======
  if (is.null(likelihood) && is.null(type)) {
    if (NN > 2000L) {
      likelihood <- "Marginal"
      type <- "Pairwise"
    } else {
      likelihood <- "Full"
      type <- "Standard"
      neighb <- NULL
    }
  }

  if (!(likelihood %in% c("Full", "Marginal", "Conditional"))) {
    stop("Unsupported likelihood", call. = FALSE)
  }
  if (identical(likelihood, "Full") && !identical(type, "Standard")) {
    stop("Full likelihood must be coupled with type = 'Standard'", call. = FALSE)
  }
  if (!identical(likelihood, "Full") && !identical(type, "Pairwise")) {
    stop("Composite likelihood must be coupled with type = 'Pairwise'", call. = FALSE)
  }

  ## ====== robust core handling ======
  core_cfg <- .GeoTestCoreConfig(parallel, ncores, B, strict_positive = FALSE)
  parallel <- core_cfg$parallel
  ncores <- core_cfg$ncores

  ## ====== compute d_min and h0 ======
  nn <- nabor::knn(coordx, coordx, k = 2)
  d_min <- min(nn$nn.dists[, 2])

  if (!is.finite(d_min) || d_min <= 0) {
    stop(
      "The minimum inter-site distance is not positive; check duplicated coordinates.",
      call. = FALSE
    )
  }

  h0_val <- if (is.null(h0)) max(eps, d_min - eps) else as.numeric(h0)
  if (!is.finite(h0_val) || h0_val <= eps) {
    stop("Invalid h0: it must be greater than the numerical tolerance", call. = FALSE)
  }

  ## ====== finite fallback bounds ======
  dat_sd <- try(stats::sd(data), silent = TRUE)
  if (inherits(dat_sd, "try-error") || !is.finite(dat_sd)) dat_sd <- 1
  dat_sd <- max(dat_sd, 1e-6)
  BIG_M <- max(1e3, 100 * dat_sd)
  SILL_MIN <- 1e-8

  mins <- apply(coordx, 2, min)
  maxs <- apply(coordx, 2, max)
  bbox_diag <- sqrt(sum((maxs - mins)^2))
  U_big <- min(max(10 * h0_val, 5 * bbox_diag, 1e-4), max_scale_cap)

  if (U_big <= h0_val) {
    stop(
      "The practical upper bound for scale must be greater than h0. ",
      "Reduce h0 or increase max_scale_cap in the implementation.",
      call. = FALSE
    )
  }

  msg <- sprintf("Testing H0: compact support <= %.5f", h0_val)
  if (is.null(h0)) msg <- paste0(msg, " (independence test)")
  message(msg)

  ## ====== parameter validation ======
  corr_names <- CorrParam(corrmodel)
  nuis_names <- NuisParam(model)
  expected <- sort(unique(c(corr_names, nuis_names)))
  provided <- sort(c(names(start), names(fixed)))

  if (!identical(expected, provided)) {
    missing <- setdiff(expected, provided)
    extra <- setdiff(provided, expected)
    msg_parts <- character(0)
    if (length(missing)) {
      msg_parts <- c(msg_parts, paste0("missing: ", paste(missing, collapse = ", ")))
    }
    if (length(extra)) {
      msg_parts <- c(msg_parts, paste0("extra: ", paste(extra, collapse = ", ")))
    }
    stop(
      "Parameter names in start+fixed must match CorrParam(corrmodel)+NuisParam(model). ",
      paste(msg_parts, collapse = " | "),
      call. = FALSE
    )
  }

  if (!("scale" %in% expected)) {
    stop("The selected correlation model does not contain a scale parameter", call. = FALSE)
  }

  ## The unrestricted model must estimate scale. If the user supplied it
  ## in fixed, use that value as its starting value and remove it from fixed.
  base_start <- start
  base_fixed <- fixed
  if (!("scale" %in% names(base_start))) {
    if (!("scale" %in% names(base_fixed))) {
      stop("scale must be supplied in either start or fixed", call. = FALSE)
    }
    base_start$scale <- base_fixed$scale
    base_fixed$scale <- NULL
  }

  ## ====== bounds utilities ======
  norm_bounds <- function(b, which = c("lower", "upper"), start_names) {
    .GeoTestNormalizeBounds(b, which, start_names, require_list = TRUE)
  }

  widen_boxes <- function(lo, up, starts = NULL, positive = character()) {
    nms <- union(names(lo), names(up))

    for (nm in nms) {
      if (is.null(lo[[nm]]) || !is.finite(lo[[nm]])) {
        lo[[nm]] <- if (nm %in% positive) SILL_MIN else -BIG_M
      }
      if (is.null(up[[nm]]) || !is.finite(up[[nm]])) {
        up[[nm]] <- BIG_M
      }

      w <- up[[nm]] - lo[[nm]]
      if (!is.finite(w) || w <= 0) {
        mid <- if (!is.null(starts) && !is.null(starts[[nm]]) &&
                   is.finite(starts[[nm]])) {
          starts[[nm]]
        } else {
          0
        }
        lo[[nm]] <- mid - min_box_width / 2
        up[[nm]] <- mid + min_box_width / 2
      } else if (w < min_box_width) {
        extra <- (min_box_width - w) / 2
        lo[[nm]] <- lo[[nm]] - extra
        up[[nm]] <- up[[nm]] + extra
      }

      if (nm %in% positive) {
        if (lo[[nm]] < eps) lo[[nm]] <- eps
        if (up[[nm]] <= lo[[nm]]) up[[nm]] <- lo[[nm]] + min_box_width
      }
    }

    list(lower = lo, upper = up)
  }

  clamp_start <- function(starts, lo, up) {
    for (nm in names(starts)) {
      x <- as.numeric(starts[[nm]])
      l <- as.numeric(lo[[nm]])
      u <- as.numeric(up[[nm]])

      if (!is.finite(x)) x <- (l + u) / 2

      gap <- max(1e-10, 1e-8 * max(1, abs(l), abs(u)))
      if (x <= l) x <- l + gap
      if (x >= u) x <- u - gap

      if (!(x > l && x < u)) x <- (l + u) / 2
      starts[[nm]] <- x
    }
    starts
  }

  ## ====== fit helpers ======
  do_fit <- function(data_, starts, fixeds, lows, ups) {
    fixed_arg <- if (length(fixeds) > 0L) fixeds else NULL

    GeoFit(
      data = data_,
      coordx = coordx,
      corrmodel = corrmodel,
      model = model,
      start = starts,
      fixed = fixed_arg,
      optimizer = optimizer,
      lower = lows,
      upper = ups,
      neighb = neighb,
      likelihood = likelihood,
      type = type
    )
  }

  loglik_value <- .GeoTestLogLikValue

  same_pairs <- function(f0, f1) .GeoTestSamePairs(f0, f1, type)

  make_start_from_fit <- function(fit_obj, template, lo, up,
                                  force_scale_above_h0 = FALSE) {
    ans <- template
    est <- unlist(fit_obj$param, use.names = TRUE)
    common <- intersect(names(ans), names(est))
    ans[common] <- as.list(est[common])

    if (force_scale_above_h0 && "scale" %in% names(ans)) {
      scale_try <- max(
        as.numeric(ans$scale),
        h0_val * (1 + 1e-3),
        as.numeric(lo$scale) * (1 + 1e-3)
      )
      ans$scale <- min(scale_try, as.numeric(up$scale) * (1 - 1e-8))
    }

    clamp_start(ans, lo, up)
  }

  fit_unrestricted <- function(data_, start_primary, fixed1, lower1, upper1,
                               fit0_ref = NULL, ll0_ref = NA_real_) {
    candidates <- list()

    first <- suppressWarnings(
      try(do_fit(data_, start_primary, fixed1, lower1, upper1), silent = TRUE)
    )
    if (is.finite(loglik_value(first))) candidates[[length(candidates) + 1L]] <- first

    best_ll <- if (length(candidates)) {
      max(vapply(candidates, loglik_value, numeric(1)))
    } else {
      -Inf
    }

    tol <- 1e-8 * max(1, abs(ll0_ref), abs(best_ll))
    need_fallback <- !is.finite(best_ll) ||
      (is.finite(ll0_ref) && best_ll < ll0_ref - tol)

    fallback_used <- FALSE

    if (need_fallback && !is.null(fit0_ref)) {
      fallback_used <- TRUE
      start_fallback <- make_start_from_fit(
        fit0_ref,
        template = start_primary,
        lo = lower1,
        up = upper1,
        force_scale_above_h0 = TRUE
      )

      retry <- suppressWarnings(
        try(do_fit(data_, start_fallback, fixed1, lower1, upper1), silent = TRUE)
      )
      if (is.finite(loglik_value(retry))) {
        candidates[[length(candidates) + 1L]] <- retry
      }
    }

    if (!length(candidates)) {
      return(list(fit = NULL, fallback = fallback_used, valid = FALSE))
    }

    ll <- vapply(candidates, loglik_value, numeric(1))
    best <- candidates[[which.max(ll)]]
    best_ll <- max(ll)
    tol <- 1e-8 * max(1, abs(ll0_ref), abs(best_ll))

    valid <- !is.finite(ll0_ref) || best_ll >= ll0_ref - tol

    list(fit = best, fallback = fallback_used, valid = valid)
  }

  ## ====== construct H0 and H1 parameterizations ======
  start1 <- base_start
  fixed1 <- base_fixed

  lower1 <- norm_bounds(lower, "lower", start_names = sort(names(start1)))
  upper1 <- norm_bounds(upper, "upper", start_names = sort(names(start1)))

  if (is.null(lower1$scale) || !is.finite(lower1$scale) || lower1$scale <= 0) {
    lower1$scale <- eps
  }
  if (is.null(upper1$scale) || !is.finite(upper1$scale) || upper1$scale <= 0) {
    upper1$scale <- U_big
  }
  upper1$scale <- min(as.numeric(upper1$scale), U_big)

  wb1 <- widen_boxes(
    lower1,
    upper1,
    starts = start1,
    positive = intersect(c("sill", "scale"), names(start1))
  )
  lower1 <- wb1$lower
  upper1 <- wb1$upper

  lower1$scale <- max(as.numeric(lower1$scale), eps)
  upper1$scale <- min(as.numeric(upper1$scale), U_big)

  if (lower1$scale >= h0_val) {
    stop(
      "The lower bound for scale under H1 must be smaller than h0 so that H0 is nested in H1.",
      call. = FALSE
    )
  }
  if (upper1$scale <= h0_val) {
    stop(
      "The upper bound for scale under H1 must be greater than h0.",
      call. = FALSE
    )
  }

  start1 <- clamp_start(start1, lower1, upper1)

  start0 <- start1
  fixed0 <- fixed1

  if (is.null(h0)) {
    start0$scale <- NULL
    fixed0$scale <- h0_val

    lower0 <- lower1[names(start0)]
    upper0 <- upper1[names(start0)]
  } else {
    lower0 <- lower1
    upper0 <- upper1
    upper0$scale <- min(as.numeric(upper0$scale), h0_val)
    lower0$scale <- max(as.numeric(lower0$scale), eps)

    if (lower0$scale >= upper0$scale) {
      stop(
        "The bounds under H0 do not define a non-empty interval for scale.",
        call. = FALSE
      )
    }

    start0 <- clamp_start(start0, lower0, upper0)
  }

  ## ====== observed fits ======
  fit0 <- suppressWarnings(
    try(do_fit(data, start0, fixed0, lower0, upper0), silent = TRUE)
  )
  ll0 <- loglik_value(fit0)
  if (!is.finite(ll0)) {
    stop("The observed fit under H0 failed", call. = FALSE)
  }

  unrestricted_obs <- fit_unrestricted(
    data_ = data,
    start_primary = start1,
    fixed1 = fixed1,
    lower1 = lower1,
    upper1 = upper1,
    fit0_ref = fit0,
    ll0_ref = ll0
  )

  if (!isTRUE(unrestricted_obs$valid) || is.null(unrestricted_obs$fit)) {
    stop(
      "The observed unrestricted fit failed to attain a log-likelihood at least as large as H0.",
      call. = FALSE
    )
  }

  fit1 <- unrestricted_obs$fit
  ll1 <- loglik_value(fit1)

  if (!same_pairs(fit0, fit1)) {
    stop(
      "H0 and H1 must use the same pairs for the composite-likelihood ratio.",
      call. = FALSE
    )
  }

  ll_tol_obs <- 1e-8 * max(1, abs(ll0), abs(ll1))
  Lambda_obs <- if (ll1 < ll0 && ll1 >= ll0 - ll_tol_obs) {
    0
  } else {
    2 * (ll1 - ll0)
  }

  if (!is.finite(Lambda_obs) || Lambda_obs < -ll_tol_obs) {
    stop("Invalid observed likelihood-ratio statistic", call. = FALSE)
  }
  Lambda_obs <- max(0, Lambda_obs)

  ## ====== early exit if unrestricted estimate satisfies H0 ======
  scale_hat <- unlist(c(fit1$param, fit1$fixed), use.names = TRUE)["scale"]
  if (length(scale_hat) == 1L && is.finite(scale_hat) &&
      scale_hat <= h0_val + 1e-10 * max(1, h0_val)) {
    cat("Unrestricted fit satisfies H0 (bootstrap not needed)\n")

    return(invisible(list(
      d_min = round(d_min, 5),
      h0 = round(h0_val, 5),
      lambda_obs = 0,
      pvalue = 1,
      seed = seed,
      method = method,
      B_requested = B,
      bootstrap_successful = 0L,
      bootstrap_success_rate = NA_real_,
      fallback_observed = unrestricted_obs$fallback,
      fallback_bootstrap = 0L,
      B_rep = numeric(0),
      fit_H0 = fit0,
      fit_H1 = fit1
    )))
  }

  ## ====== H0 parameters for simulation ======
  message("Running parametric bootstrap (B = ", B, ") ...")
  parH0 <- c(fit0$param, fit0$fixed)
  param_H0 <- as.list(parH0[expected])

  ## Use observed estimates as bootstrap starting values.
  start0_boot <- as.list(fit0$param)
  start1_boot <- as.list(fit1$param)
  start0_boot <- clamp_start(start0_boot, lower0, upper0)
  start1_boot <- clamp_start(start1_boot, lower1, upper1)

  ## ====== simulations ======
  if (method == "cholesky") {
    data_sim <- GeoModels::GeoSim(
      coordx = coordx,
      corrmodel = corrmodel,
      param = param_H0,
      model = model,
      sparse = TRUE,
      nrep = B,
      progress = FALSE
    )
  } else {
    old_plan_after_sim <- NULL
    if (isTRUE(parallel) && requireNamespace("future", quietly = TRUE)) {
      old_plan_after_sim <- future::plan()
    }

    data_sim <- GeoModels::GeoSimapprox(
      coordx = coordx,
      corrmodel = corrmodel,
      param = param_H0,
      model = model,
      method = "TB",
      L = L,
      nrep = B,
      parallel = parallel,
      ncores = ncores,
      progress = FALSE
    )

    if (!is.null(old_plan_after_sim)) {
      try(future::plan(old_plan_after_sim), silent = TRUE)
    }
  }

  ## ====== bootstrap helper ======
  estimate_fun <- function(Zb) {
    f0 <- suppressWarnings(
      try(do_fit(Zb, start0_boot, fixed0, lower0, upper0), silent = TRUE)
    )
    ll0b <- loglik_value(f0)

    if (!is.finite(ll0b)) {
      return(c(lambda = NA_real_, fallback = 0))
    }

    unrestricted <- fit_unrestricted(
      data_ = Zb,
      start_primary = start1_boot,
      fixed1 = fixed1,
      lower1 = lower1,
      upper1 = upper1,
      fit0_ref = f0,
      ll0_ref = ll0b
    )

    if (!isTRUE(unrestricted$valid) || is.null(unrestricted$fit)) {
      return(c(lambda = NA_real_, fallback = as.numeric(unrestricted$fallback)))
    }

    f1 <- unrestricted$fit
    ll1b <- loglik_value(f1)

    if (!same_pairs(f0, f1)) {
      return(c(lambda = NA_real_, fallback = as.numeric(unrestricted$fallback)))
    }

    tol <- 1e-8 * max(1, abs(ll0b), abs(ll1b))
    if (ll1b < ll0b - tol) {
      return(c(lambda = NA_real_, fallback = as.numeric(unrestricted$fallback)))
    }

    lambda <- if (ll1b < ll0b) 0 else 2 * (ll1b - ll0b)
    c(lambda = max(0, lambda), fallback = as.numeric(unrestricted$fallback))
  }

  ## ====== sequential or parallel bootstrap ======
  if (!parallel) {
    message(sprintf("Performing %d estimations under H0 and H1 sequentially...", B))
    result_list <- .GeoTestSequentialApply(
      data_sim, B, estimate_fun, progress = use_progressr
    )
  } else if (!requireNamespace("future", quietly = TRUE) ||
             !requireNamespace("future.apply", quietly = TRUE)) {
    warning(
      "Parallel bootstrap requested but future/future.apply is unavailable; using sequential evaluation.",
      call. = FALSE
    )
    result_list <- .GeoTestSequentialApply(
      data_sim, B, estimate_fun, progress = FALSE
    )
  } else {
    message(sprintf(
      "Performing %d estimations under H0 and H1 using %d cores...",
      B, ncores
    ))
    data_files <- .GeoTestSaveReplicates(data_sim, B, "supp_space")
    on.exit(unlink(data_files), add = TRUE)
    rm(data_sim)
    gc(verbose = FALSE, full = TRUE)
    pb_boot <- if (use_progressr) {
      progressr::progressor(steps = B)
    } else {
      NULL
    }
    result_list <- .GeoTestParallelFileApply(
      data_files = data_files,
      FUN = estimate_fun,
      ncores = ncores,
      progressor = pb_boot,
      globals = list(coordx = coordx, estimate_fun = estimate_fun),
      progress_every = TRUE,
      error_value = c(lambda = NA_real_, fallback = 0)
    )
  }

  result_mat <- do.call(rbind, result_list)
  Lambda_all <- result_mat[, "lambda"]
  fallback_all <- result_mat[, "fallback"]
  valid <- is.finite(Lambda_all)
  Lambda_boot <- Lambda_all[valid]

  ## ====== results and diagnostics ======
  n_valid <- length(Lambda_boot)
  failed <- B - n_valid
  success_rate <- n_valid / B
  fallback_boot <- sum(fallback_all > 0, na.rm = TRUE)

  if (failed > 0L) {
    warning(sprintf(
      "%d/%d bootstrap replications failed (%.1f%%)",
      failed, B, 100 * failed / B
    ), call. = FALSE)
  }

  min_required <- min(20L, B)
  if (n_valid < min_required) {
    stop(sprintf(
      "Only %d/%d bootstrap replications were successful; at least %d are required.",
      n_valid, B, min_required
    ), call. = FALSE)
  }

  if (success_rate < 0.8) {
    warning(
      sprintf(
        "Only %.1f%% of bootstrap replications were successful; results may be unreliable.",
        100 * success_rate
      ),
      call. = FALSE
    )
  }

  if (B < 99L) {
    warning(
      "B < 99 gives a coarse bootstrap p-value; use a larger B for final analyses.",
      call. = FALSE
    )
  }

  pval <- (1 + sum(Lambda_boot >= Lambda_obs)) / (n_valid + 1)

  res <- list(
    d_min = round(d_min, 5),
    h0 = round(h0_val, 5),
    lambda_obs = round(Lambda_obs, 5),
    pvalue = round(pval, 5),
    seed = seed,
    method = method,
    L = if (method == "TB") L else NULL,
    B_requested = B,
    bootstrap_successful = n_valid,
    bootstrap_failed = failed,
    bootstrap_success_rate = success_rate,
    fallback_observed = unrestricted_obs$fallback,
    fallback_bootstrap = fallback_boot,
    B_rep = Lambda_boot,
    fit_H0 = fit0,
    fit_H1 = fit1
  )

  cat("\n=== FINAL RESULTS ===\n")
  cat("Minimum distance (d_min)           =", res$d_min, "\n")
  cat("H0 threshold (compact support)    <=", res$h0, "\n")
  cat("Observed LRT statistic             =", res$lambda_obs, "\n")
  cat("Valid bootstrap reps               =", n_valid, "/", B, "\n")
  cat("Adaptive H1 fallbacks              =", fallback_boot, "bootstrap reps\n")
  cat("p-value                            =", res$pvalue, "\n\n")

  concl <- if (is.null(h0)) {
    if (res$pvalue < 0.05) {
      "Reject H0 => significant spatial dependence detected"
    } else {
      "Do not reject H0 => data consistent with spatial independence"
    }
  } else {
    if (res$pvalue < 0.05) {
      sprintf(
        "Reject H0 => compact support significantly > %.4f (stronger dependence)",
        h0
      )
    } else {
      sprintf(
        "Do not reject H0 => insufficient evidence that scale > %.4f",
        h0
      )
    }
  }
  cat("Conclusion:", concl, "\n\n")

  invisible(res)
}