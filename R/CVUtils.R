## Internal utilities for GeoCV.
## These functions centralize cross-validation plumbing; spatial, space-time
## and bivariate fold construction remain in GeoCV itself.

.GeoCVMetricNames <- c(
  "rmse", "mae", "mad", "lscore", "brie",
  "crps", "pit_mean", "intscore", "coverage",
  "n_requested", "n_valid", "n_failed", "n_invalid_mse", "n_valid_prob"
)

.GeoCVSummarizeScores <- function(obs, pred, mse, gaussian_scores,
                                   strict = getOption("GeoModels.cv_strict", FALSE)) {
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict))
    stop("GeoModels.cv_strict must be TRUE or FALSE.", call. = FALSE)
  obs <- as.numeric(c(obs))
  pred <- as.numeric(c(pred))
  mse <- as.numeric(c(mse))

  if (length(obs) != length(pred))
    stop("length of observations and predictions does not match", call. = FALSE)
  if (!gaussian_scores && length(mse) == 0L)
    mse <- rep(NA_real_, length(pred))
  if (length(mse) != length(pred))
    stop("length of predictive MSE and predictions does not match", call. = FALSE)

  ok_err <- is.finite(obs) & is.finite(pred)
  ok_mse <- is.finite(mse) & mse > 0
  n_invalid_mse <- if (gaussian_scores) sum(ok_err & !ok_mse) else NA_integer_
  if (strict && (any(!ok_err) || (gaussian_scores && any(ok_err & !ok_mse))))
    stop("Cross-validation excluded invalid predictions or predictive MSE (strict mode).",
         call. = FALSE)
  err <- obs[ok_err] - pred[ok_err]
  ans <- c(
    rmse = if (length(err)) sqrt(mean(err^2)) else NA_real_,
    mae = if (length(err)) mean(abs(err)) else NA_real_,
    mad = if (length(err)) stats::median(abs(err)) else NA_real_,
    lscore = NA_real_, brie = NA_real_, crps = NA_real_,
    pit_mean = NA_real_, intscore = NA_real_, coverage = NA_real_,
    n_requested = length(obs), n_valid = sum(ok_err), n_failed = sum(!ok_err),
    n_invalid_mse = n_invalid_mse,
    n_valid_prob = if (gaussian_scores) sum(ok_err & ok_mse) else NA_integer_
  )

  if (gaussian_scores) {
    ok_prob <- is.finite(obs) & is.finite(pred) & is.finite(mse) & mse > 0
    if (any(ok_prob)) {
      pp <- GeoScores(
        obs[ok_prob], pred = pred[ok_prob], mse = mse[ok_prob],
        score = c("brie", "crps", "lscore", "pit", "intscore", "coverage")
      )
      required <- c("LogScore", "Brier", "CRPS", "PIT", "IS95", "Cvg95")
      missing_scores <- setdiff(required, names(pp))
      if (length(missing_scores)) {
        stop(
          "GeoScores did not return the required components: ",
          paste(missing_scores, collapse = ", "), call. = FALSE
        )
      }
      pit_value <- mean(as.numeric(pp$PIT), na.rm = TRUE)
      if (is.nan(pit_value)) pit_value <- NA_real_
      ans[c("lscore", "brie", "crps", "pit_mean", "intscore", "coverage")] <- c(
        as.numeric(pp$LogScore), as.numeric(pp$Brier), as.numeric(pp$CRPS),
        pit_value, as.numeric(pp$IS95), as.numeric(pp$Cvg95)
      )
    }
  }

  ans[.GeoCVMetricNames]
}

.GeoCVBindScores <- function(x, metric_names = .GeoCVMetricNames) {
  if (!is.list(x) || !length(x))
    stop("No cross-validation results were returned", call. = FALSE)
  sizes <- vapply(x, length, integer(1))
  bad <- which(sizes != length(metric_names))
  if (length(bad)) {
    stop(
      "Invalid score vector returned by cross-validation iteration(s): ",
      paste(bad, collapse = ", "), call. = FALSE
    )
  }
  ans <- do.call(rbind, x)
  if (is.null(dim(ans)) || ncol(ans) != length(metric_names))
    stop("Cross-validation score results have invalid dimensions", call. = FALSE)
  colnames(ans) <- metric_names
  warning_messages <- unlist(lapply(x, attr, which = "GeoModels.task_warnings"),
                             use.names = FALSE)
  if (length(warning_messages)) {
    tab <- sort(table(warning_messages), decreasing = TRUE)
    shown <- head(tab, 10L)
    warning("Cross-validation warnings: ",
            paste(sprintf("%s [n=%d]", names(shown), as.integer(shown)), collapse = "; "),
            if (length(tab) > 10L) "; additional warning types are stored in the result.",
            call. = FALSE)
    attr(ans, "task_warnings") <- tab
  }
  ans
}

.GeoCVRunIterations <- function(iteration_fun, K, parallel, ncores, progress,
                                metric_names = .GeoCVMetricNames) {
  strict <- getOption("GeoModels.cv_strict", FALSE)
  warn_policy <- getOption("warn", 0L)
  run_iteration <- function(i) {
    old_options <- options(GeoModels.cv_strict = strict)
    on.exit(options(old_options), add = TRUE)
    messages <- character()
    ans <- withCallingHandlers(iteration_fun(i), warning = function(w) {
      ## Do not turn warn=2 into a silently accepted result.
      if (warn_policy >= 2L) stop(conditionMessage(w), call. = FALSE)
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    attr(ans, "GeoModels.task_warnings") <- messages
    ans
  }
  n.cores <- .GeoResolveWorkers(parallel, ncores, n_jobs = K)
  if (n.cores <= 1L) {
    out <- vector("list", K)
    if (isTRUE(progress)) {
      cat("Performing", K, "cross-validations...\n")
      pb <- txtProgressBar(min = 0, max = K, style = 3)
      on.exit(close(pb), add = TRUE)
    }
    for (i in seq_len(K)) {
      out[[i]] <- run_iteration(i)
      if (isTRUE(progress)) setTxtProgressBar(pb, i)
    }
    return(.GeoCVBindScores(out, metric_names = metric_names))
  }

  old_plan <- future::plan()
  on.exit(try(future::plan(old_plan), silent = TRUE), add = TRUE)
  future::plan(future::multisession, workers = n.cores)

  if (isTRUE(progress)) {
    out <- progressr::with_progress({
      p <- progressr::progressor(along = seq_len(K))
      future.apply::future_lapply(
        seq_len(K),
        function(i) {
          z <- run_iteration(i)
          p()
          z
        },
        future.seed = TRUE, future.stdout = FALSE, future.conditions = "condition"
      )
    })
  } else {
    out <- future.apply::future_lapply(
      seq_len(K), run_iteration,
      future.seed = TRUE, future.stdout = FALSE, future.conditions = "condition"
    )
  }
  .GeoCVBindScores(out, metric_names = metric_names)
}

.GeoCVMakeTestSamples <- function(n, K, n.fold, groups = NULL) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 2L)
    stop("At least two observations are required for cross-validation.", call. = FALSE)

  ntest <- max(1L, min(n - 1L, as.integer(round(n * n.fold))))
  if (is.null(groups))
    return(lapply(seq_len(K), function(i) sort(sample.int(n, ntest))))

  if (length(groups) != n)
    stop("Internal error: fold groups are not aligned with observations.", call. = FALSE)
  groups <- as.integer(factor(groups, levels = unique(groups)))
  group_sizes <- tabulate(groups)
  max_test <- sum(pmax(group_sizes - 1L, 0L))
  if (max_test < 1L) {
    stop(
      "Space-time cross-validation cannot remove observations while retaining at least one training observation at every time.",
      call. = FALSE
    )
  }
  ntest <- min(ntest, max_test)

  lapply(seq_len(K), function(i) {
    remaining <- group_sizes
    selected <- integer(ntest)
    n_selected <- 0L
    for (idx in sample.int(n)) {
      g <- groups[idx]
      if (remaining[g] > 1L) {
        n_selected <- n_selected + 1L
        selected[n_selected] <- idx
        remaining[g] <- remaining[g] - 1L
        if (n_selected == ntest) break
      }
    }
    if (n_selected != ntest)
      stop("Unable to construct a valid space-time cross-validation fold.", call. = FALSE)
    sort(selected)
  })
}

.GeoCVFlattenDesign <- function(X, ns, context) {
  if (is.null(X)) return(NULL)
  total <- sum(ns)
  if (is.list(X) && !is.data.frame(X) && !is.matrix(X)) {
    if (length(X) != length(ns))
      stop(context, ": X must contain one design matrix per block.", call. = FALSE)
    blocks <- lapply(seq_along(ns), function(i) {
      z <- as.matrix(X[[i]])
      if (!is.numeric(z) || nrow(z) != ns[i] || any(!is.finite(z)))
        stop(context, ": X block ", i, " is not aligned with the data.", call. = FALSE)
      z
    })
    p <- unique(vapply(blocks, ncol, integer(1)))
    if (length(p) != 1L)
      stop(context, ": all X blocks must have the same number of columns.", call. = FALSE)
    return(do.call(rbind, blocks))
  }
  X <- as.matrix(X)
  if (!is.numeric(X) || nrow(X) != total || any(!is.finite(X)))
    stop(context, ": X must have one finite numeric row per observation.", call. = FALSE)
  X
}

.GeoCVFlattenObservationVector <- function(x, ns, name) {
  if (is.null(x)) stop(name, " is missing", call. = FALSE)
  if (length(x) == 1L && !is.list(x)) {
    if (!is.numeric(x) || !is.finite(x))
      stop(name, " must be finite and numeric", call. = FALSE)
    return(as.numeric(x))
  }
  z <- if (is.list(x)) unlist(x, use.names = FALSE) else as.numeric(x)
  if (!is.numeric(z) || length(z) != sum(ns) || any(!is.finite(z)))
    stop(name, " must be scalar or contain one finite value per observation", call. = FALSE)
  as.numeric(z)
}

.GeoCVSubsetObservationVector <- function(x, idx) {
  if (is.null(x) || length(x) == 1L) x else x[idx]
}

.GeoCVNormalizeCountSize <- function(x, model, name = "fit$n") {
  if (is.null(x)) return(x)
  model <- .GeoPredictionCanonicalModel(as.character(model)[1L])
  z <- if (is.list(x)) unlist(x, use.names = FALSE) else as.numeric(x)
  if (!length(z) || any(!is.finite(z)))
    stop(name, " must contain finite numeric values", call. = FALSE)

  if (model %in% c("BinomialNeg", "BinomialNegZINB")) {
    if (any(z < 1) || any(abs(z - round(z)) > sqrt(.Machine$double.eps)))
      stop(name, " must contain positive integers for direct Negative-Binomial fields.", call. = FALSE)
    if (length(z) > 1L && any(z != z[1L]))
      stop("GeoCV: for direct Negative-Binomial fields n is the common number r of successes and cannot vary by observation.", call. = FALSE)
    return(as.integer(round(z[1L])))
  }

  if (model %in% c("Binary", "Bernoulli", "Geom", "Geometric")) {
    if (any(abs(z - 1) > sqrt(.Machine$double.eps)))
      stop("GeoCV: Binary/Bernoulli and Geometric fields require the common size parameter n=1.", call. = FALSE)
    return(1L)
  }
  x
}

.GeoCVVaryingFixedMean <- function(fixed, name, expected) {
  value <- fixed[[name]]
  if (is.null(value) || length(value) <= 1L) return(NULL)
  value <- as.numeric(unlist(value, use.names = FALSE))
  if (length(value) != expected || any(!is.finite(value)))
    stop("fit$fixed$", name, " is not aligned with the observations.", call. = FALSE)
  value
}

.GeoCVExplicitCoordinates <- function(fit, context) {
  .GeoKrigCoordinateMatrix(
    coordx = fit$coordx, coordy = fit$coordy, coordz = fit$coordz,
    grid = isTRUE(fit$grid), context = context
  )
}

.GeoCVDropMeanParameters <- function(z, prefix) {
  if (is.null(z) || is.null(names(z))) return(z)
  z[!grepl(paste0("^", prefix, "[0-9]*$"), names(z))]
}

.GeoCVRefit <- function(fit, data, X, n, fixed, start, model,
                        optimizer, lower, upper, sparse, method, varest,
                        coordx = NULL, coordx_dyn = NULL, coordt = NULL,
                        universal = FALSE) {
  args <- list(
    data = data, coordx = coordx, coordx_dyn = coordx_dyn,
    corrmodel = fit$corrmodel, X = X,
    likelihood = fit$likelihood, type = fit$type, grid = FALSE,
    copula = fit$copula, anisopars = fit$anisopars,
    est.aniso = fit$est.aniso, model = model,
    radius = fit$radius, n = n,
    maxdist = fit$maxdist, neighb = fit$neighb,
    p_neighb = fit$p_neighb, distance = fit$distance,
    optimizer = optimizer, lower = lower, upper = upper,
    start = start, fixed = fixed,
    weighted = isTRUE(fit$weighted),
    thin_method = if (!is.null(fit$thin_method)) fit$thin_method else "bernoulli",
    memdist = TRUE, method = method, sparse = sparse, varest = varest
  )
  if (!is.null(coordt)) {
    args$coordt <- coordt
    args$maxtime <- fit$maxtime
  }
  ans <- suppressWarnings(.GeoWithDuplicateCheckSuppressed(do.call(GeoFit, args)))
  if (isTRUE(universal) && is.null(ans$varcov)) {
    stop("The refitted model did not return varcov required for Universal kriging.", call. = FALSE)
  }
  ans
}

.GeoCVPredict <- function(fit, loc, Xloc, Mloc, nloc,
                          local, neighb, maxdist, maxtime,
                          type_krig, sparse, method,
                          time = NULL, which = NULL) {
  args <- list(
    fit, loc = loc, mse = TRUE,
    Xloc = Xloc, Mloc = Mloc, nloc = nloc,
    type_krig = type_krig, sparse = sparse, method = method
  )
  if (!is.null(time)) args$time <- time
  if (!is.null(which)) args$which <- which

  if (!local) return(.GeoWithDuplicateCheckSuppressed(do.call(GeoKrig, args)))

  args$neighb <- neighb
  args$maxdist <- maxdist
  if (!is.null(time)) args$maxtime <- maxtime
  args$progress <- FALSE
  .GeoWithDuplicateCheckSuppressed(do.call(GeoKrigloc, args))
}


.GeoCVConditionalPredict <- function(
    fit, loc, Xloc, Mloc, nloc,
    nrep, n_iter, mcmc_thin
) {
  model <- .GeoPredictionCanonicalModel(as.character(fit$model)[1L])
  if (!(model %in% c("Binomial", "BinomialNeg"))) {
    stop(
      "GeoCV conditional prediction is currently available only for direct Binomial and BinomialNeg random fields.",
      call. = FALSE
    )
  }
  if (!is.null(fit$copula)) {
    stop(
      "GeoCV conditional prediction is currently available only for the direct count random fields, not copula models.",
      call. = FALSE
    )
  }

  ans <- .GeoWithDuplicateCheckSuppressed(do.call(
    GeoSimcond,
    list(
      estobj = fit,
      loc = loc,
      Xloc = Xloc,
      Mloc = Mloc,
      nloc = nloc,
      nrep = nrep,
      n_iter = n_iter,
      mcmc_thin = mcmc_thin,
      method = "Cholesky",
      local = FALSE,
      sparse = FALSE,
      parallel = FALSE,
      progress = FALSE
    )
  ))

  pred <- as.numeric(ans$cond_mean)
  if (length(pred) != nrow(as.matrix(loc)) || any(!is.finite(pred))) {
    stop(
      "GeoSimcond returned an invalid conditional-mean prediction inside GeoCV.",
      call. = FALSE
    )
  }

  list(pred = pred, cond_var = as.numeric(ans$cond_var))
}
