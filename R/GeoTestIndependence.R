########################################
## Parametric bootstrap test for absence of spatial dependence
########################################
GeoTestIndependence <- function(data, coordx,
                                start, fixed = NULL,
                                corrmodel = "Matern",
                                model = "Gaussian",
                                optimizer = "bobyqa",
                                optimizer_ind = "Nelder-Mead",
                                lower = NULL, upper = NULL,
                                X = NULL, n = 1,
                                distance = "Eucl", radius = 1,
                                anisopars = NULL,
                                est.aniso = c(FALSE, FALSE),
                                sparse = FALSE,
                                B = 1000,
                                parallel = TRUE,
                                ncores = NULL,
                                progress = TRUE,
                                seed = NULL) {

  ## This implementation compares:
  ## H0: independent observations at distinct sites;
  ## H1: a spatially dependent model with nugget < 1.
  ##
  ## H0 is fitted with the marginal independence likelihood, whereas H1 is
  ## fitted with the full likelihood. Parametric bootstrap calibration is used
  ## because nugget = 1 is a boundary point and the correlation parameters are
  ## not identified under H0.

  ## ====== preserve global future/progress settings ======
  future_plan_original <- NULL
  if (requireNamespace("future", quietly = TRUE)) {
    future_plan_original <- future::plan()
    on.exit(
      try(future::plan(future_plan_original), silent = TRUE),
      add = TRUE
    )
  }

  old_handlers <- NULL
  if (requireNamespace("progressr", quietly = TRUE)) {
    old_handlers <- progressr::handlers()
    on.exit({
      if (length(old_handlers) > 0L) {
        progressr::handlers(old_handlers)
      } else {
        progressr::handlers("default")
      }
    }, add = TRUE)
  }

  ## ====== reproducibility without permanent RNG side effects ======
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("seed must be NULL or a single finite numeric value", call. = FALSE)
    }

    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv,
                           inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }

    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv,
                   inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)

    set.seed(as.integer(seed))
  }

  ## ====== input validation ======
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

  if (!is.null(X)) {
    X <- as.matrix(X)
    if (nrow(X) != length(data)) {
      stop("nrow(X) must match data length", call. = FALSE)
    }
    if (!is.numeric(X) || any(!is.finite(X))) {
      stop("X must contain only finite numeric values", call. = FALSE)
    }
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
    stop("start and fixed cannot contain duplicated parameter names",
         call. = FALSE)
  }

  overlap <- intersect(names(start), names(fixed))
  if (length(overlap) > 0L) {
    stop(
      "Parameters cannot be present in both start and fixed: ",
      paste(overlap, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.character(corrmodel) || length(corrmodel) != 1L) {
    stop("corrmodel must be a single character string", call. = FALSE)
  }
  if (!is.character(model) || length(model) != 1L) {
    stop("model must be a single character string", call. = FALSE)
  }
  if (!is.character(optimizer) || length(optimizer) != 1L) {
    stop("optimizer must be a single character string", call. = FALSE)
  }
  if (!is.character(optimizer_ind) || length(optimizer_ind) != 1L) {
    stop("optimizer_ind must be a single character string", call. = FALSE)
  }

  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 1) {
    stop("B must be a positive integer", call. = FALSE)
  }
  B <- as.integer(B)

  if (!is.logical(parallel) || length(parallel) != 1L) {
    stop("parallel must be logical (TRUE or FALSE)", call. = FALSE)
  }
  if (!is.logical(progress) || length(progress) != 1L) {
    stop("progress must be logical (TRUE or FALSE)", call. = FALSE)
  }
  if (!is.logical(sparse) || length(sparse) != 1L) {
    stop("sparse must be logical (TRUE or FALSE)", call. = FALSE)
  }
  if (!is.numeric(radius) || length(radius) != 1L ||
      !is.finite(radius) || radius <= 0) {
    stop("radius must be a positive finite number", call. = FALSE)
  }
  if (!is.logical(est.aniso) || length(est.aniso) != 2L) {
    stop("est.aniso must be a logical vector of length two", call. = FALSE)
  }
  if (!is.null(anisopars) && !is.list(anisopars)) {
    stop("anisopars must be NULL or a list", call. = FALSE)
  }

  if (!is.null(lower) &&
      (!is.list(lower) || is.null(names(lower)) || any(names(lower) == ""))) {
    stop("lower must be NULL or a named list", call. = FALSE)
  }
  if (!is.null(upper) &&
      (!is.list(upper) || is.null(names(upper)) || any(names(upper) == ""))) {
    stop("upper must be NULL or a named list", call. = FALSE)
  }

  ## ====== progress handlers ======
  use_progressr <- isTRUE(progress) &&
    requireNamespace("progressr", quietly = TRUE)

  if (isTRUE(progress) && !use_progressr) {
    warning(
      "progress=TRUE but progressr is not available; progress bars disabled.",
      call. = FALSE
    )
  }

  if (use_progressr) {
    progressr::handlers(global = TRUE)
    progressr::handlers(progressr::handler_txtprogressbar(clear = TRUE))
  } else if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::handlers("void")
  }

  ## ====== core handling ======
  coremax <- parallel::detectCores()
  if (is.na(coremax) || coremax <= 1L) {
    parallel <- FALSE
    ncores <- 1L
  } else {
    if (is.null(ncores)) ncores <- min(coremax - 1L, B)
    if (!is.numeric(ncores) || length(ncores) != 1L ||
        !is.finite(ncores) || ncores < 1) {
      stop("ncores must be NULL or a positive integer", call. = FALSE)
    }
    ncores <- max(1L, min(as.integer(ncores), B, coremax))
    if (ncores == 1L) parallel <- FALSE
  }

  ## ====== minimum distance, used only for exact null simulation ======
  ## Use the same metric and radius as the fitted spatial model. Full
  ## likelihood fitting already requires dense calculations, so constructing
  ## this distance matrix does not change the intended computational regime.
  dmat <- GeoModels::GeoDistances(
    coordx,
    distance = distance,
    radius = radius
  )
  diag(dmat) <- Inf
  d_min <- min(dmat)
  rm(dmat)

  if (!is.finite(d_min) || d_min <= 0) {
    stop(
      "The minimum inter-site distance is not positive; check duplicated coordinates.",
      call. = FALSE
    )
  }

  ## ====== construct H1 parameterization ======
  corr_names <- CorrParam(corrmodel)
  if (!length(corr_names)) {
    stop("Unable to determine the correlation parameters", call. = FALSE)
  }

  start1 <- start
  fixed1 <- fixed

  ## nugget must be estimated under H1. If supplied in fixed, promote it to
  ## start and move it away from the lower boundary when necessary.
  if (!("nugget" %in% names(start1))) {
    if (!("nugget" %in% names(fixed1))) {
      stop(
        "nugget must be supplied in start or fixed; it is estimated under H1.",
        call. = FALSE
      )
    }

    nugget_start <- as.numeric(fixed1$nugget)
    fixed1$nugget <- NULL

    if (!is.finite(nugget_start)) nugget_start <- 0.1
    nugget_start <- min(max(nugget_start, 0.05), 0.95)
    start1$nugget <- nugget_start
  }

  ## Current implementation is for univariate models with one nugget
  ## proportion. More complicated nugget parameterizations require a separate
  ## definition of the independence boundary.
  other_nuggets <- grep("^nugget", c(names(start1), names(fixed1)),
                         value = TRUE)
  other_nuggets <- setdiff(unique(other_nuggets), "nugget")
  if (length(other_nuggets) > 0L) {
    stop(
      "GeoTestIndependence currently supports models with a single parameter named 'nugget'.",
      call. = FALSE
    )
  }

  ## Validate that user-supplied bounds refer only to estimated H1 parameters.
  if (!is.null(lower)) {
    unknown <- setdiff(names(lower), names(start1))
    if (length(unknown) > 0L) {
      stop(
        "lower contains names not estimated under H1: ",
        paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
  }
  if (!is.null(upper)) {
    unknown <- setdiff(names(upper), names(start1))
    if (length(unknown) > 0L) {
      stop(
        "upper contains names not estimated under H1: ",
        paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
  }

  ## Bounded optimizers used by GeoFit require finite boxes.  In particular,
  ## bobyqa cannot work with -Inf lower bounds.  Preserve all finite bounds
  ## supplied by the user and replace only non-finite endpoints by broad,
  ## data-adaptive finite values.
  dat_sd <- suppressWarnings(stats::sd(data))
  if (!is.finite(dat_sd) || dat_sd <= 0) dat_sd <- 1
  dat_max <- suppressWarnings(max(abs(data)))
  if (!is.finite(dat_max)) dat_max <- 1
  BIG_M <- max(1e3, 100 * dat_sd, 100 * dat_max)
  POS_EPS <- 1e-8

  complete_bounds <- function(b, parameter_names, default) {
    out <- setNames(rep(default, length(parameter_names)), parameter_names)
    if (!is.null(b) && length(b) > 0L) {
      out[names(b)] <- unlist(b, use.names = FALSE)
    }
    as.list(out)
  }

  make_finite_bounds <- function(lows, ups, starts) {
    positive_names <- intersect(
      c("sill", "scale", "scale_s", "scale_t"),
      names(starts)
    )

    for (nm in names(starts)) {
      lo <- as.numeric(lows[[nm]])
      up <- as.numeric(ups[[nm]])

      if (length(lo) != 1L || is.na(lo)) lo <- -Inf
      if (length(up) != 1L || is.na(up)) up <- Inf

      if (!is.finite(lo)) {
        lo <- if (nm %in% positive_names) POS_EPS else -BIG_M
      }
      if (!is.finite(up)) up <- BIG_M

      if (nm %in% positive_names && lo <= 0) lo <- POS_EPS

      if (lo >= up) {
        stop(
          "Invalid bounds for parameter '", nm,
          "': lower must be strictly smaller than upper",
          call. = FALSE
        )
      }

      lows[[nm]] <- lo
      ups[[nm]] <- up
    }

    list(lower = lows, upper = ups)
  }

  lower1 <- complete_bounds(lower, names(start1), -Inf)
  upper1 <- complete_bounds(upper, names(start1), Inf)

  nugget_eps <- 1e-8
  lower_nugget_user <- if (!is.null(lower1$nugget)) {
    as.numeric(lower1$nugget)
  } else {
    0
  }
  upper_nugget_user <- if (!is.null(upper1$nugget)) {
    as.numeric(upper1$nugget)
  } else {
    Inf
  }

  if (!is.finite(lower_nugget_user)) lower_nugget_user <- 0
  if (lower_nugget_user < 0 || lower_nugget_user >= 1) {
    stop("The lower bound for nugget must lie in [0, 1)", call. = FALSE)
  }

  ## The alternative must approach nugget = 1; otherwise H0 is not in the
  ## closure of H1. GeoModels requires nugget < 1, hence 1 - nugget_eps.
  if (is.finite(upper_nugget_user) &&
      upper_nugget_user < 1 - 1e-6) {
    stop(
      "The upper bound for nugget must be 1 (or omitted) so that H0 is in the closure of H1.",
      call. = FALSE
    )
  }

  lower1$nugget <- lower_nugget_user
  upper1$nugget <- 1 - nugget_eps

  fb1 <- make_finite_bounds(lower1, upper1, start1)
  lower1 <- fb1$lower
  upper1 <- fb1$upper

  clamp_value <- function(x, lo = -Inf, up = Inf) {
    x <- as.numeric(x)
    lo <- as.numeric(lo)
    up <- as.numeric(up)

    if (!is.finite(x)) {
      if (is.finite(lo) && is.finite(up)) {
        x <- (lo + up) / 2
      } else if (is.finite(lo)) {
        x <- lo + max(1, abs(lo))
      } else if (is.finite(up)) {
        x <- up - max(1, abs(up))
      } else {
        x <- 0
      }
    }

    if (is.finite(lo) && x <= lo) {
      x <- lo + max(1e-10, 1e-8 * max(1, abs(lo)))
    }
    if (is.finite(up) && x >= up) {
      x <- up - max(1e-10, 1e-8 * max(1, abs(up)))
    }
    x
  }

  clamp_start <- function(starts, lows, ups) {
    for (nm in names(starts)) {
      lo <- if (!is.null(lows[[nm]])) lows[[nm]] else -Inf
      up <- if (!is.null(ups[[nm]])) ups[[nm]] else Inf
      starts[[nm]] <- clamp_value(starts[[nm]], lo, up)
    }
    starts
  }

  start1 <- clamp_start(start1, lower1, upper1)

  ## ====== construct H0 parameterization ======
  dependence_names <- unique(c(corr_names, "nugget"))

  start0 <- start1[setdiff(names(start1), dependence_names)]
  fixed0 <- fixed1[setdiff(names(fixed1), dependence_names)]

  if (length(start0) < 1L) {
    stop(
      "At least one marginal parameter must be estimated under independence.",
      call. = FALSE
    )
  }

  lower0 <- lower1[names(start0)]
  upper0 <- upper1[names(start0)]

  fixed0_arg <- if (length(fixed0)) fixed0 else NULL
  fixed1_arg <- if (length(fixed1)) fixed1 else NULL

  ## ====== fit helpers ======
  do_fit0 <- function(data_, starts = start0,
                      optimizer0 = optimizer_ind) {
    ## The independence likelihood is a low-dimensional marginal fit.
    ## In particular, the current CompIndLik2 implementation can return a
    ## non-finite objective with bounded optimizers such as bobyqa even when
    ## the same fit is regular with Nelder-Mead.  Therefore H0 has its own
    ## optimizer, with Nelder-Mead as the stable default.
    args0 <- list(
      data = data_,
      coordx = coordx,
      corrmodel = NULL,
      model = model,
      start = starts,
      fixed = fixed0_arg,
      optimizer = optimizer0,
      likelihood = "Marginal",
      type = "Independence",
      X = X,
      n = n,
      distance = distance,
      radius = radius
    )

    bounded0 <- optimizer0 %in% c(
      "L-BFGS-B", "nlminb", "nmkb", "multinlminb",
      "multiNelder-Mead", "bobyqa", "sbplx"
    )
    if (bounded0) {
      args0$lower <- lower0
      args0$upper <- upper0
    }

    do.call(GeoFit, args0)
  }

  fit_independence <- function(data_, starts = start0) {
    optimizers0 <- unique(c(optimizer_ind, "Nelder-Mead"))
    candidates <- list()
    used <- character(0)

    for (opt0 in optimizers0) {
      ff <- suppressWarnings(
        try(do_fit0(data_, starts = starts, optimizer0 = opt0),
            silent = TRUE)
      )
      ll <- loglik_value(ff)
      if (is.finite(ll)) {
        candidates[[length(candidates) + 1L]] <- ff
        used <- c(used, opt0)
      }
    }

    if (!length(candidates)) {
      return(list(fit = NULL, loglik = NA_real_, optimizer = NA_character_))
    }

    vals <- vapply(candidates, loglik_value, numeric(1))
    ii <- which.max(vals)
    list(fit = candidates[[ii]], loglik = vals[ii], optimizer = used[ii])
  }

  do_fit1 <- function(data_, starts = start1) {
    GeoFit(
      data = data_,
      coordx = coordx,
      corrmodel = corrmodel,
      model = model,
      start = starts,
      fixed = fixed1_arg,
      optimizer = optimizer,
      lower = lower1,
      upper = upper1,
      likelihood = "Full",
      type = "Standard",
      X = X,
      n = n,
      distance = distance,
      radius = radius,
      anisopars = anisopars,
      est.aniso = est.aniso,
      sparse = sparse
    )
  }

  loglik_value <- function(fit_obj) {
    if (inherits(fit_obj, "try-error") || is.null(fit_obj)) {
      return(NA_real_)
    }
    val <- if (!is.null(fit_obj$logCompLik)) {
      fit_obj$logCompLik
    } else {
      fit_obj$logLik
    }
    val <- as.numeric(val)
    if (length(val) != 1L || !is.finite(val)) NA_real_ else val
  }

  update_marginal_start <- function(template, fit0_obj) {
    out <- template
    est <- unlist(fit0_obj$param, use.names = TRUE)
    common <- intersect(names(out), names(est))
    out[common] <- as.list(est[common])
    clamp_start(out, lower1, upper1)
  }

  make_h1_candidates <- function(primary, fit0_obj = NULL,
                                 include_fallbacks = FALSE) {
    candidates <- list(clamp_start(primary, lower1, upper1))

    if (!include_fallbacks) return(candidates)

    base <- primary
    if (!is.null(fit0_obj)) {
      base <- update_marginal_start(base, fit0_obj)
    }

    nugget_grid <- c(0.95, 0.75, 0.50, 0.25)
    nugget_grid <- nugget_grid[
      nugget_grid > as.numeric(lower1$nugget) &
        nugget_grid < as.numeric(upper1$nugget)
    ]

    for (ng in nugget_grid) {
      ss <- base
      ss$nugget <- ng
      candidates[[length(candidates) + 1L]] <-
        clamp_start(ss, lower1, upper1)
    }

    ## Remove exact duplicate starts.
    keys <- vapply(candidates, function(x) {
      paste(names(x), signif(unlist(x, use.names = FALSE), 12),
            sep = "=", collapse = "|")
    }, character(1))
    candidates[!duplicated(keys)]
  }

  fit_spatial <- function(data_, start_primary, fit0_ref, ll0_ref) {
    fit_candidates <- list()

    first <- suppressWarnings(
      try(do_fit1(data_, start_primary), silent = TRUE)
    )
    if (is.finite(loglik_value(first))) {
      fit_candidates[[1L]] <- first
    }

    best_ll <- if (length(fit_candidates)) {
      max(vapply(fit_candidates, loglik_value, numeric(1)))
    } else {
      -Inf
    }

    tol <- 1e-8 * max(1, abs(ll0_ref), abs(best_ll))
    need_fallback <- !is.finite(best_ll) || best_ll < ll0_ref - tol
    fallback_used <- FALSE

    if (need_fallback) {
      fallback_used <- TRUE
      starts <- make_h1_candidates(
        start_primary,
        fit0_obj = fit0_ref,
        include_fallbacks = TRUE
      )

      ## The first start has already been attempted.
      if (length(starts) > 1L) {
        for (ii in 2:length(starts)) {
          ff <- suppressWarnings(
            try(do_fit1(data_, starts[[ii]]), silent = TRUE)
          )
          if (is.finite(loglik_value(ff))) {
            fit_candidates[[length(fit_candidates) + 1L]] <- ff
          }
        }
      }
    }

    if (!length(fit_candidates)) {
      return(list(
        fit = NULL,
        loglik = NA_real_,
        fallback = fallback_used,
        boundary = FALSE,
        valid = FALSE
      ))
    }

    vals <- vapply(fit_candidates, loglik_value, numeric(1))
    idx <- which.max(vals)
    best_fit <- fit_candidates[[idx]]
    best_ll <- vals[idx]
    tol <- 1e-8 * max(1, abs(ll0_ref), abs(best_ll))

    ## Since H0 belongs to the closure of H1, the effective unrestricted
    ## supremum cannot be below ll0. If the best spatial fit does not improve
    ## ll0 after the adaptive retries, the MLE is treated as a boundary
    ## solution and the LR statistic is zero.
    boundary <- best_ll <= ll0_ref + tol

    list(
      fit = best_fit,
      loglik = best_ll,
      fallback = fallback_used,
      boundary = boundary,
      valid = TRUE
    )
  }

  ## ====== observed fits ======
  message("Testing H0: absence of spatial dependence (independence)")

  ind_obs <- fit_independence(data, starts = start0)
  fit0 <- ind_obs$fit
  ll0 <- ind_obs$loglik
  optimizer_ind_used <- ind_obs$optimizer

  if (!is.finite(ll0) || is.null(fit0)) {
    stop(
      "The observed fit under independence failed with both optimizer_ind='",
      optimizer_ind, "' and the Nelder-Mead fallback.",
      call. = FALSE
    )
  }

  spatial_obs <- fit_spatial(
    data_ = data,
    start_primary = start1,
    fit0_ref = fit0,
    ll0_ref = ll0
  )

  if (!isTRUE(spatial_obs$valid) || is.null(spatial_obs$fit)) {
    stop("The observed spatial fit failed", call. = FALSE)
  }

  fit1 <- spatial_obs$fit
  ll1 <- spatial_obs$loglik

  if (isTRUE(spatial_obs$boundary)) {
    Lambda_obs <- 0
  } else {
    Lambda_obs <- 2 * (ll1 - ll0)
  }

  if (!is.finite(Lambda_obs) || Lambda_obs < 0) {
    stop("Invalid observed likelihood-ratio statistic", call. = FALSE)
  }

  par1 <- unlist(c(fit1$param, fit1$fixed), use.names = TRUE)
  nugget_hat <- if ("nugget" %in% names(par1)) {
    as.numeric(par1["nugget"])
  } else {
    NA_real_
  }

  ## A zero LR statistic implies a bootstrap p-value equal to one because all
  ## bootstrap LR statistics are nonnegative.
  if (Lambda_obs == 0) {
    cat("The unrestricted fit is on the independence boundary; bootstrap not needed.\n")

    return(invisible(list(
      lambda_obs = 0,
      pvalue = 1,
      nugget_hat = nugget_hat,
      boundary_observed = TRUE,
      optimizer_independence = optimizer_ind_used,
      optimizer_spatial = optimizer,
      seed = seed,
      B_requested = B,
      bootstrap_successful = 0L,
      bootstrap_failed = 0L,
      bootstrap_success_rate = NA_real_,
      fallback_observed = spatial_obs$fallback,
      fallback_bootstrap = 0L,
      boundary_bootstrap = 0L,
      B_rep = numeric(0),
      fit_H0 = fit0,
      fit_H1 = fit1
    )))
  }

  ## ====== parameters for exact independent null simulation ======
  ## GeoSim does not accept nugget = 1. Exact independence at the observed
  ## sites is therefore generated with an auxiliary compactly supported
  ## GenWend correlation whose support is strictly below d_min. This changes
  ## only the simulation device, not the null marginal model.
  par0_all <- unlist(c(fit0$param, fit0$fixed), use.names = TRUE)
  marginal_names <- unique(c(names(start0), names(fixed0)))
  missing_marginal <- setdiff(marginal_names, names(par0_all))
  if (length(missing_marginal) > 0L) {
    stop(
      "Unable to recover fitted marginal parameters under H0: ",
      paste(missing_marginal, collapse = ", "),
      call. = FALSE
    )
  }

  param_H0_marginal <- as.list(par0_all[marginal_names])
  sim_scale <- 0.5 * d_min
  param_H0_sim <- c(
    param_H0_marginal,
    list(
      nugget = 0,
      scale = sim_scale,
      smooth = 0.5,
      power2 = 4
    )
  )

  message("Running parametric bootstrap under independence (B = ", B, ") ...")

  data_sim <- GeoModels::GeoSim(
    coordx = coordx,
    corrmodel = "GenWend",
    distance = distance,
    model = model,
    n = n,
    param = param_H0_sim,
    radius = radius,
    sparse = TRUE,
    X = X,
    nrep = B,
    progress = FALSE
  )

  get_rep <- function(sim, b) {
    xx <- sim$data
    if (is.list(xx)) return(xx[[b]])
    if (is.null(dim(xx))) return(xx)
    xx[, b]
  }

  ## Use observed estimates as bootstrap starting values.
  start0_boot <- as.list(fit0$param)
  start1_boot <- as.list(fit1$param)
  start1_boot <- clamp_start(start1_boot, lower1, upper1)

  estimate_fun <- function(Zb) {
    ind_b <- fit_independence(Zb, starts = start0_boot)
    f0 <- ind_b$fit
    ll0b <- ind_b$loglik

    if (!is.finite(ll0b) || is.null(f0)) {
      return(c(lambda = NA_real_, fallback = 0, boundary = 0))
    }

    spatial <- fit_spatial(
      data_ = Zb,
      start_primary = start1_boot,
      fit0_ref = f0,
      ll0_ref = ll0b
    )

    if (!isTRUE(spatial$valid) || is.null(spatial$fit)) {
      return(c(
        lambda = NA_real_,
        fallback = as.numeric(spatial$fallback),
        boundary = 0
      ))
    }

    if (isTRUE(spatial$boundary)) {
      lambda <- 0
    } else {
      lambda <- 2 * (spatial$loglik - ll0b)
    }

    if (!is.finite(lambda) || lambda < 0) {
      return(c(
        lambda = NA_real_,
        fallback = as.numeric(spatial$fallback),
        boundary = as.numeric(spatial$boundary)
      ))
    }

    c(
      lambda = lambda,
      fallback = as.numeric(spatial$fallback),
      boundary = as.numeric(spatial$boundary)
    )
  }

  ## ====== sequential or parallel bootstrap ======
  if (!parallel) {
    message(sprintf(
      "Performing %d fits under H0 and H1 sequentially...", B
    ))

    if (use_progressr) {
      result_list <- progressr::with_progress({
        pb <- progressr::progressor(along = seq_len(B))
        out <- vector("list", B)
        for (b in seq_len(B)) {
          out[[b]] <- estimate_fun(get_rep(data_sim, b))
          pb(sprintf("b=%d", b))
        }
        out
      })
    } else {
      result_list <- vector("list", B)
      for (b in seq_len(B)) {
        result_list[[b]] <- estimate_fun(get_rep(data_sim, b))
      }
    }
  } else {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      warning(
        "Parallel bootstrap requested but future/future.apply is unavailable; using sequential evaluation.",
        call. = FALSE
      )

      result_list <- vector("list", B)
      for (b in seq_len(B)) {
        result_list[[b]] <- estimate_fun(get_rep(data_sim, b))
      }
    } else {
      message(sprintf(
        "Performing %d fits under H0 and H1 using %d cores...",
        B, ncores
      ))

      old_plan <- future::plan()
      on.exit(try(future::plan(old_plan), silent = TRUE), add = TRUE)
      future::plan(future::multisession, workers = ncores)

      ## Save one simulated dataset per file so each worker reads only its own
      ## replicate instead of receiving the complete simulation object.
      temp_dir <- tempdir()
      tag <- paste0(
        format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()
      )
      data_files <- vapply(seq_len(B), function(b) {
        fp <- file.path(
          temp_dir,
          sprintf("independence_%s_%05d.rds", tag, b)
        )
        saveRDS(get_rep(data_sim, b), fp, compress = FALSE)
        fp
      }, character(1))
      on.exit(unlink(data_files), add = TRUE)

      rm(data_sim)
      gc(verbose = FALSE, full = TRUE)

      old_max_size <- getOption("future.globals.maxSize")
      on.exit(options(future.globals.maxSize = old_max_size), add = TRUE)
      required_size <- as.numeric(utils::object.size(coordx)) +
        as.numeric(utils::object.size(estimate_fun))
      options(
        future.globals.maxSize = max(
          500 * 1024^2,
          1.5 * required_size
        )
      )

      if (use_progressr) {
        result_list <- progressr::with_progress({
          pb <- progressr::progressor(along = seq_len(B))
          future.apply::future_lapply(
            seq_len(B),
            function(b) {
              ans <- estimate_fun(readRDS(data_files[[b]]))
              pb(sprintf("b=%d", b))
              ans
            },
            future.seed = TRUE
          )
        })
      } else {
        result_list <- future.apply::future_lapply(
          seq_len(B),
          function(b) estimate_fun(readRDS(data_files[[b]])),
          future.seed = TRUE
        )
      }
    }
  }

  result_mat <- do.call(rbind, result_list)
  Lambda_all <- result_mat[, "lambda"]
  fallback_all <- result_mat[, "fallback"]
  boundary_all <- result_mat[, "boundary"]

  valid <- is.finite(Lambda_all)
  Lambda_boot <- Lambda_all[valid]

  n_valid <- length(Lambda_boot)
  failed <- B - n_valid
  success_rate <- n_valid / B
  fallback_boot <- sum(fallback_all > 0, na.rm = TRUE)
  boundary_boot <- sum(boundary_all[valid] > 0, na.rm = TRUE)

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
    warning(sprintf(
      "Only %.1f%% of bootstrap replications were successful; results may be unreliable.",
      100 * success_rate
    ), call. = FALSE)
  }

  if (B < 99L) {
    warning(
      "B < 99 gives a coarse bootstrap p-value; use a larger B for final analyses.",
      call. = FALSE
    )
  }

  pval <- (1 + sum(Lambda_boot >= Lambda_obs)) / (n_valid + 1)

  res <- list(
    lambda_obs = Lambda_obs,
    pvalue = pval,
    nugget_hat = nugget_hat,
    boundary_observed = spatial_obs$boundary,
    optimizer_independence = optimizer_ind_used,
    optimizer_spatial = optimizer,
    seed = seed,
    B_requested = B,
    bootstrap_successful = n_valid,
    bootstrap_failed = failed,
    bootstrap_success_rate = success_rate,
    fallback_observed = spatial_obs$fallback,
    fallback_bootstrap = fallback_boot,
    boundary_bootstrap = boundary_boot,
    simulation_support = sim_scale,
    B_rep = Lambda_boot,
    fit_H0 = fit0,
    fit_H1 = fit1
  )

  cat("\n=== FINAL RESULTS ===\n")
  cat("Null hypothesis                    = spatial independence\n")
  cat("Observed LR statistic              =", round(Lambda_obs, 5), "\n")
  cat("Optimizer under H0                 =", optimizer_ind_used, "\n")
  cat("Optimizer under H1                 =", optimizer, "\n")
  cat("Estimated nugget under H1          =", round(nugget_hat, 5), "\n")
  cat("Valid bootstrap reps               =", n_valid, "/", B, "\n")
  cat("Boundary bootstrap solutions       =", boundary_boot, "\n")
  cat("Adaptive H1 fallbacks              =", fallback_boot,
      "bootstrap reps\n")
  cat("p-value                            =", round(pval, 5), "\n\n")

  conclusion <- if (pval < 0.05) {
    "Reject H0 => significant spatial dependence detected"
  } else {
    "Do not reject H0 => insufficient evidence of spatial dependence"
  }
  cat("Conclusion:", conclusion, "\n\n")

  invisible(res)
}