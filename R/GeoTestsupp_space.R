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
                               B = 1000,
                               likelihood = NULL,
                               type = NULL,
                               method = "Cholesky",
                               parallel = TRUE,
                               ncores = NULL,
                               progress = TRUE) {

  ## ====== input validations  ======
  if (!is.matrix(coordx)) coordx <- as.matrix(coordx)
  if (!is.numeric(data)) data <- as.numeric(data)
  
  if (nrow(coordx) != length(data)) 
    stop("coordx rows (", nrow(coordx), ") must match data length (", length(data), ")")
  if (B < 1) 
    stop("B must be at least 1")
  if (!is.null(neighb) && neighb < 1) 
    stop("neighb must be positive")
  if (!is.null(h0) && (h0 <= 0 || !is.finite(h0))) 
    stop("h0 must be positive and finite")
  
  ## ====== COSTANTI ======
  max_scale_cap <- 1e6
  min_box_width <- 1e-3
  corrmodel <- "GenWend"
  eps <- 1e-6
  NN <- nrow(coordx)

  ## ====== SETUP PROGRESS HANDLERS ======
  if (progress) {
    progressr::handlers(global = TRUE)
    progressr::handlers(progressr::handler_txtprogressbar(clear = TRUE))
  } else {
    progressr::handlers("void")
  }

  ## ====== AUTOMATIC LIKELIHOOD SELECTION ======
  if (is.null(likelihood) || is.null(type)) {
    if (NN > 2000) {
      likelihood <- "Marginal"
      type <- "Pairwise"
    } else {
      likelihood <- "Full"
      type <- "Standard"
      neighb <- NULL
    }
  }

  ## ====== COMPUTE d_min E h0 ======
  nn <- nabor::knn(coordx, coordx, k = 2)
  d_min <- min(nn$nn.dists[, 2])
  h0_val <- if (is.null(h0)) max(eps, d_min - eps) else h0
  if (h0_val <= eps) stop("Invalid h0: must be positive.")

  ## ====== BOUND FINITI PER mean/sill ======
  dat_sd <- try(sd(as.numeric(data), na.rm = TRUE), silent = TRUE)
  if (inherits(dat_sd, "try-error") || !is.finite(dat_sd)) dat_sd <- 1
  dat_sd <- max(dat_sd, 1e-6)
  BIG_M <- max(1e3, 100 * dat_sd)
  SILL_MIN <- 1e-8

  ## ====== U_big PER scale (upper in H1) ======
  mins <- apply(coordx, 2, min)
  maxs <- apply(coordx, 2, max)
  bbox_diag <- sqrt(sum((maxs - mins)^2))
  U_big <- max(10 * h0_val, 5 * bbox_diag, 1e-4)
  U_big <- min(U_big, max_scale_cap)

  if (!is.null(h0) && h0 > U_big) {
    warning(sprintf("Provided h0 (%.4f) exceeds practical upper bound (%.4f); results may be unstable.", h0, U_big))
  }

  msg <- sprintf("Testing H0: compact support <= %.5f  ", h0_val)
  if (is.null(h0)) msg <- paste0(msg, "  (independence test)")
  message(msg)

  ## ====== parameters validations  ======
  corr_names <- CorrParam(corrmodel)
  nuis_names <- NuisParam(model)
  expected <- sort(unique(c(corr_names, nuis_names)))
  provided <- sort(unique(c(names(start), names(fixed))))
  
  if (!identical(expected, provided)) {
    missing <- setdiff(expected, provided)
    extra <- setdiff(provided, expected)
    msg_parts <- c()
    if (length(missing)) msg_parts <- c(msg_parts, paste0("missing: ", paste(missing, collapse = ", ")))
    if (length(extra)) msg_parts <- c(msg_parts, paste0("extra: ", paste(extra, collapse = ", ")))
    stop("Parameter names in start+fixed must match CorrParam(corrmodel)+NuisParam(model). ",
         paste(msg_parts, collapse = " | "))
  }

  ## ====== UTILITY: NORMALIZZA BOUNDS ======
  start_names_ref <- sort(names(start))
  norm_bounds <- function(b, which = c("lower", "upper"), start_names = start_names_ref) {
    which <- match.arg(which)
    if (is.null(b)) b <- list()
    if (length(b) && is.null(names(b))) {
      stop(which, " must be a *named* list with names matching start")
    }
    if (length(setdiff(names(b), start_names))) {
      stop(which, " contains names not in start: ",
           paste(setdiff(names(b), start_names), collapse = ", "))
    }
    base_val <- if (which == "lower") -Inf else Inf
    out <- setNames(rep(base_val, length(start_names)), start_names)
    if (length(b)) out[names(b)] <- unlist(b, use.names = FALSE)
    as.list(out)
  }

  widen_boxes <- function(lo, up, starts = NULL, positive = character()) {
    nms <- union(names(lo), names(up))
    for (nm in nms) {
      if (is.null(lo[[nm]]) || !is.finite(lo[[nm]])) 
        lo[[nm]] <- if (nm %in% positive) SILL_MIN else -BIG_M
      if (is.null(up[[nm]]) || !is.finite(up[[nm]])) 
        up[[nm]] <- BIG_M
      
      w <- up[[nm]] - lo[[nm]]
      if (!is.finite(w) || w <= 0) {
        mid <- if (!is.null(starts) && !is.null(starts[[nm]]) && is.finite(starts[[nm]])) 
          starts[[nm]] else 0
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

  ## ====== HELPER: DO FIT (passa sempre il dataset) ======
  do_fit <- function(data_, starts, fixeds, lows, ups) {
    GeoFit(
      data = data_, coordx = coordx, corrmodel = corrmodel,
      start = starts, fixed = fixeds,
      optimizer = optimizer,
      lower = lows, upper = ups,
      neighb = neighb, likelihood = likelihood, type = type
    )
  }

  ## ====== FIT H0 ======
  start0 <- start
  fixed0 <- fixed

  if (is.null(h0)) {
    # Test di indipendenza: scale fissata a h0_val
    if ("scale" %in% names(start0)) {
      start0$scale <- NULL
    }
    fixed0$scale <- h0_val
  }

  # Bounds per H0, solo sui nomi start0
  lower_h0 <- lower
  upper_h0 <- upper
  if (is.null(h0)) {
    if (!is.null(lower_h0) && "scale" %in% names(lower_h0)) lower_h0$scale <- NULL
    if (!is.null(upper_h0) && "scale" %in% names(upper_h0)) upper_h0$scale <- NULL
  }
  
  lower0 <- norm_bounds(lower_h0, "lower", start_names = sort(names(start0)))
  upper0 <- norm_bounds(upper_h0, "upper", start_names = sort(names(start0)))

  if (!is.null(h0) && "scale" %in% names(start0)) {
    upper0$scale <- min(
      if (is.null(upper0$scale) || !is.finite(upper0$scale)) h0_val else upper0$scale, 
      h0_val
    )
    if (is.null(lower0$scale) || !is.finite(lower0$scale) || lower0$scale <= 0) 
      lower0$scale <- eps
  }

  wb0 <- widen_boxes(lower0, upper0, starts = start0,
                     positive = intersect(c("sill", "scale"), names(start0)))
  lower0 <- wb0$lower
  upper0 <- wb0$upper

  fit0 <- do_fit(data, start0, fixed0, lower0, upper0)

  ## ====== FIT H1 ======
  start1 <- start
  fixed1 <- fixed
  if (!("scale" %in% names(start1)) && "scale" %in% names(fixed1)) {
    start1$scale <- fixed1$scale
    fixed1$scale <- NULL
  }
  
  if ("scale" %in% names(start1)) {
    if (!is.finite(start1$scale) || start1$scale <= 0) {
      start1$scale <- max(h0_val * 1.5, eps)
    }
    start1$scale <- min(start1$scale, 0.8 * U_big)
  }

  lower1 <- norm_bounds(lower, "lower", start_names = sort(names(start1)))
  upper1 <- norm_bounds(upper, "upper", start_names = sort(names(start1)))

  if ("scale" %in% names(start1)) {
    if (is.null(lower1$scale) || !is.finite(lower1$scale) || lower1$scale <= 0) 
      lower1$scale <- eps
    if (is.null(upper1$scale) || !is.finite(upper1$scale) || upper1$scale <= 0) 
      upper1$scale <- U_big
    upper1$scale <- min(upper1$scale, U_big)
  }

  wb1 <- widen_boxes(lower1, upper1, starts = start1,
                     positive = intersect(c("sill", "scale"), names(start1)))
  lower1 <- wb1$lower
  upper1 <- wb1$upper

  fit1 <- do_fit(data, start1, fixed1, lower1, upper1)

  ## ====== LRT OSSERVATO (mantiene logCompLik) ======
  ll0 <- fit0$logCompLik
  ll1 <- fit1$logCompLik
  if (ll1 < ll0) {
    message("Note: H1 log-likelihood (", round(ll1, 3), 
            ") < H0 log-likelihood (", round(ll0, 3), "), LRT set to 0")
  }
  Lambda_obs <- max(0, 2 * (ll1 - ll0))

  ## ====== EARLY EXIT SE H1 SODDISFA H0 ======
  if (!is.null(fit1$param$scale) && is.finite(fit1$param$scale) && 
      fit1$param$scale <= h0_val + 1e-12) {
    cat("Unconstrained fit satisfies H0 (bootstrap not needed)\n")
    return(list(
      d_min = round(d_min, 5), 
      h0 = round(h0_val, 5),
      lambda_obs = 0, 
      pvalue = 1, 
      fit_H0 = fit0, 
      fit_H1 = fit1
    ))
  }

  ## ====== PARAMETRI H0 PER SIMULAZIONE ======
  message("Running parametric bootstrap (B = ", B, ") ...")
  parH0 <- c(fit0$param, fit0$fixed)
  param_H0 <- as.list(parH0[expected])

  ## ====== SIMULAZIONI ======
  if (method == "Cholesky") {
    data_sim <- GeoModels::GeoSim(
      coordx = coordx, 
      corrmodel = corrmodel,
      param = param_H0, 
      model = model,
      sparse = TRUE, 
      nrep = B,
      progress = FALSE
    )
  } else if (method == "TB") {
    data_sim <- GeoModels::GeoSimapprox(
      coordx = coordx, 
      corrmodel = corrmodel,
      param = param_H0,
      model = model,
      nrep = B,
      progress = FALSE
    )
  } else {
    stop("method must be 'Cholesky' or 'TB'.")
  }

  get_rep <- function(sim, b) {
    X <- sim$data
    if (is.list(X)) X[[b]] else X[, b]
  }

  ## ====== GESTIONE CORES ROBUSTA ======
  coremax <- parallel::detectCores()
  if (is.na(coremax) || coremax <= 1) {
    parallel <- FALSE
    ncores <- 1
  } else {
    if (is.null(ncores)) ncores <- min(coremax - 1, B)
    ncores <- max(1, min(ncores, B, coremax))
    if (ncores == 1) parallel <- FALSE
  }

  ## ====== HELPER BOOTSTRAP (usa do_fit con data passati) ======
  estimate_fun <- function(Zb) {
    f0 <- suppressWarnings(try(do_fit(Zb, start0, fixed0, lower0, upper0), silent = TRUE))
    f1 <- suppressWarnings(try(do_fit(Zb, start1, fixed1, lower1, upper1), silent = TRUE))
    if (inherits(f0, "try-error") || inherits(f1, "try-error") ||
        !is.finite(f0$logCompLik) || !is.finite(f1$logCompLik)) {
      return(NA_real_)
    }
    max(0, 2 * (f1$logCompLik - f0$logCompLik))
  }

  ## ====== BOOTSTRAP SEQUENZIALE O PARALLELO ======
  if (!parallel) {
    message(sprintf("Performing %d estimations under H0 and H1 sequentially...", B))
    if (progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
    }
    pb <- if (progress) progressr::progressor(along = seq_len(B)) else NULL
    
    Lambda_boot <- numeric(B)
    ok <- rep(TRUE, B)
    for (b in seq_len(B)) {
      if (progress) pb(sprintf("b=%d", b))
      Zb <- get_rep(data_sim, b)
      val <- estimate_fun(Zb)
      if (!is.finite(val)) { ok[b] <- FALSE } else { Lambda_boot[b] <- val }
    }
    Lambda_boot <- Lambda_boot[ok]
    
  } else {
    message(sprintf("Performing %d estimations under H0 and H1 using %d cores...", B, ncores))
    future::plan(future::multisession, workers = ncores)
    on.exit(future::plan(future::sequential), add = TRUE)
    if (progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = seq_len(B))
    }
    b <- NULL
    Lambda_boot <- foreach::foreach(
      b = seq_len(B),
      .combine = c,
      .options.future = list(
        packages = "GeoModels",
        seed = TRUE,
        stdout = FALSE,
        conditions = "none"
      )
    ) %dofuture% {
      if (progress) pb(sprintf("b=%d", b))
      Zb <- get_rep(data_sim, b)
      f0 <- suppressWarnings(try(
        GeoFit(
          data = Zb, coordx = coordx, corrmodel = corrmodel,
          start = start0, fixed = fixed0, optimizer = optimizer,
          lower = lower0, upper = upper0,
          neighb = neighb, likelihood = likelihood, type = type
        ), silent = TRUE))
      f1 <- suppressWarnings(try(
        GeoFit(
          data = Zb, coordx = coordx, corrmodel = corrmodel,
          start = start1, fixed = fixed1, optimizer = optimizer,
          lower = lower1, upper = upper1,
          neighb = neighb, likelihood = likelihood, type = type
        ), silent = TRUE))
      if (inherits(f0, "try-error") || inherits(f1, "try-error") ||
          !is.finite(f0$logCompLik) || !is.finite(f1$logCompLik)) {
        as.numeric(NA)
      } else {
        as.numeric(max(0, 2 * (f1$logCompLik - f0$logCompLik)))
      }
    }
    Lambda_boot <- Lambda_boot[is.finite(Lambda_boot)]
  }

  ## ====== RISULTATI CON WARNING INFORMATIVO ======
  failed <- B - length(Lambda_boot)
  if (failed > 0) {
    warning(sprintf("%d/%d bootstrap replications failed (%.1f%%)", 
                    failed, B, 100 * failed / B))
  }
  if (B > 0 && failed / B > 0.7) {
    warning("High failure rate in bootstrap fits; consider adjusting starts/bounds or reducing model complexity.")
  }
  
  if (!length(Lambda_boot)) {
    warning("No valid bootstrap replications; p-value set to NA")
    pval <- NA_real_
  } else {
    pval <- (1 + sum(Lambda_boot >= Lambda_obs)) / (length(Lambda_boot) + 1)
  }

  ## ====== OUTPUT FINALE ======
  res <- list(
    d_min = round(d_min, 5),
    h0 = round(h0_val, 5),
    lambda_obs = round(Lambda_obs, 5),
    pvalue = if (is.na(pval)) NA_real_ else round(pval, 5),
    B_rep = Lambda_boot,
    fit_H0 = fit0,
    fit_H1 = fit1
  )

  cat("\n=== FINAL RESULTS ===\n")
  cat("Minimum distance (d_min)           =", res$d_min, "\n")
  cat("H0 threshold (compact support)    <=", res$h0, "\n")
  cat("Observed LRT statistic             =", res$lambda_obs, "\n")
  cat("Valid bootstrap reps               =", length(res$B_rep), "/", B, "\n")
  cat("p-value                            =", if (is.na(res$pvalue)) "NA" else res$pvalue, "\n\n")

  concl <- if (is.na(res$pvalue)) {
    "Bootstrap failed (no valid replications)"
  } else if (is.null(h0)) {
    if (res$pvalue < 0.05) 
      "Reject H0 => significant spatial dependence detected"
    else 
      "Do not reject H0 => data consistent with spatial independence"
  } else {
    if (res$pvalue < 0.05) 
      sprintf("Reject H0 => compact support significantly > %.4f (stronger dependence)", h0)
    else 
      sprintf("Do not reject H0 => insufficient evidence that scale > %.4f", h0)
  }
  cat("Conclusion:", concl, "\n\n")

  invisible(res)
}
