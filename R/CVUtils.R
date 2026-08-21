## Internal utilities for GeoCV.
## These functions centralize cross-validation plumbing; spatial, space-time
## and bivariate fold construction remain in GeoCV itself.

.GeoCVMetricNames <- c(
  "rmse", "mae", "mad", "lscore", "brie",
  "crps", "pit_mean", "intscore", "coverage"
)

.GeoCVSummarizeScores <- function(obs, pred, mse, gaussian_scores) {
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
  if (!any(ok_err))
    stop("no finite observation/prediction pairs are available", call. = FALSE)

  err <- obs[ok_err] - pred[ok_err]
  ans <- c(
    rmse = sqrt(mean(err^2)),
    mae = mean(abs(err)),
    mad = stats::median(abs(err)),
    lscore = NA_real_, brie = NA_real_, crps = NA_real_,
    pit_mean = NA_real_, intscore = NA_real_, coverage = NA_real_
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

.GeoCVBindScores <- function(x) {
  if (!is.list(x) || !length(x))
    stop("No cross-validation results were returned", call. = FALSE)
  sizes <- vapply(x, length, integer(1))
  bad <- which(sizes != length(.GeoCVMetricNames))
  if (length(bad)) {
    stop(
      "Invalid score vector returned by cross-validation iteration(s): ",
      paste(bad, collapse = ", "), call. = FALSE
    )
  }
  ans <- do.call(rbind, x)
  if (is.null(dim(ans)) || ncol(ans) != length(.GeoCVMetricNames))
    stop("Cross-validation score results have invalid dimensions", call. = FALSE)
  colnames(ans) <- .GeoCVMetricNames
  ans
}

.GeoCVRunIterations <- function(iteration_fun, K, parallel, ncores, progress) {
  n.cores <- .GeoResolveWorkers(parallel, ncores, n_jobs = K)
  if (n.cores <= 1L) {
    out <- vector("list", K)
    if (isTRUE(progress)) {
      cat("Performing", K, "cross-validations...\n")
      pb <- txtProgressBar(min = 0, max = K, style = 3)
      on.exit(close(pb), add = TRUE)
    }
    for (i in seq_len(K)) {
      out[[i]] <- iteration_fun(i)
      if (isTRUE(progress)) setTxtProgressBar(pb, i)
    }
    return(.GeoCVBindScores(out))
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
          z <- iteration_fun(i)
          p()
          z
        },
        future.seed = TRUE, future.stdout = FALSE, future.conditions = "none"
      )
    })
  } else {
    out <- future.apply::future_lapply(
      seq_len(K), iteration_fun,
      future.seed = TRUE, future.stdout = FALSE, future.conditions = "none"
    )
  }
  .GeoCVBindScores(out)
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
    selected <- integer(0)
    for (idx in sample.int(n)) {
      g <- groups[idx]
      if (remaining[g] > 1L) {
        selected <- c(selected, idx)
        remaining[g] <- remaining[g] - 1L
        if (length(selected) == ntest) break
      }
    }
    if (length(selected) != ntest)
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
