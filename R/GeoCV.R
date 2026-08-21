GeoCV <- function(fit, K = 100, estimation = TRUE,
                  optimizer = NULL, lower = NULL, upper = NULL,
                  n.fold = 0.05, local = FALSE, neighb = NULL, maxdist = NULL,
                  maxtime = NULL, sparse = FALSE, type_krig = "Simple",
                  which = 1, parallel = FALSE, ncores = 6,
                  progress = TRUE, seed = NULL) {
  .dup_old <- .GeoSuppressDuplicateChecksForNestedCalls()
  on.exit(.GeoRestoreDuplicateCheckOption(.dup_old), add = TRUE)

  ## ---------- sanity checks ----------
  if (!inherits(fit, "GeoFit"))
    stop("fit must be an object of class 'GeoFit'", call. = FALSE)
  if (!is.numeric(n.fold) || length(n.fold) != 1L || !is.finite(n.fold) ||
      n.fold < 0.01 || n.fold > 0.99)
    stop("n.fold must be between 0.01 and 0.99", call. = FALSE)
  if (!is.logical(estimation) || length(estimation) != 1L || is.na(estimation))
    stop("estimation must be TRUE or FALSE", call. = FALSE)
  if (!is.logical(local) || length(local) != 1L || is.na(local))
    stop("local must be TRUE or FALSE", call. = FALSE)
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel))
    stop("parallel must be TRUE or FALSE", call. = FALSE)
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress))
    stop("progress must be TRUE or FALSE", call. = FALSE)
  if (!is.logical(sparse) || length(sparse) != 1L || is.na(sparse))
    stop("sparse must be TRUE or FALSE", call. = FALSE)

  fit_copula_phase1 <- if (is.null(fit$copula)) NULL else as.character(fit$copula)[1L]
  fit_model_code_phase1 <- suppressWarnings(as.integer(CkModel(fit$model)))
  phase1_count_copula <- !is.null(fit_copula_phase1) &&
    fit_copula_phase1 %in% c("Clayton", "SkewGaussian") &&
    fit_model_code_phase1 %in% .GeoGaussianCopulaDiscreteModels
  if (isTRUE(phase1_count_copula) && isTRUE(local))
    stop("GeoCV: Clayton/SkewGaussian copula with Poisson/Binomial/BinomialNeg is currently implemented for global GeoKrig (local=FALSE) only.", call. = FALSE)
  if (identical(fit_copula_phase1, "Clayton") &&
      fit_model_code_phase1 %in% .GeoGaussianCopulaDiscreteModels &&
      isTRUE(estimation))
    stop("GeoCV: fold refitting for discrete Clayton models still requires the Clayton phase-2 GeoFit likelihood. Use estimation=FALSE.", call. = FALSE)

  type_key <- tolower(gsub("[[:blank:]]", "", type_krig))
  if (!(type_key %in% c("simple", "universal")))
    stop("type_krig must be 'Simple' or 'Universal'", call. = FALSE)
  type_krig <- if (identical(type_key, "universal")) "Universal" else "Simple"

  tol_int <- sqrt(.Machine$double.eps)
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) || K < 2 ||
      abs(K - round(K)) > tol_int)
    stop("K must be an integer >= 2", call. = FALSE)
  K <- as.integer(K)
  if (K > 1000L) stop("K is too large", call. = FALSE)

  ncores <- .GeoValidateNcores(ncores)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        abs(seed - round(seed)) > tol_int || abs(seed) > .Machine$integer.max)
      stop("seed must be a finite integer scalar", call. = FALSE)
    seed <- as.integer(seed)
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
          rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  if (parallel &&
      (!requireNamespace("future", quietly = TRUE) ||
       !requireNamespace("future.apply", quietly = TRUE))) {
    warning(
      "Parallel execution requested but future/future.apply is unavailable; using sequential execution.",
      call. = FALSE
    )
    parallel <- FALSE
  }

  has_progressr <- requireNamespace("progressr", quietly = TRUE)
  if (parallel && progress && !has_progressr) {
    warning(
      "Package 'progressr' is required for progress bars in parallel mode; disabling progress.",
      call. = FALSE
    )
    progress <- FALSE
  }

  message("Cross-validation kriging can be time consuming ...")

  ## ---------- model classification ----------
  corr_code <- CkCorrModel(fit$corrmodel)
  spacetime <- CheckST(corr_code)
  bivariate <- CheckBiv(corr_code)
  space <- !spacetime && !bivariate

  if (bivariate && !(which %in% c(1, 2)))
    stop("which must be 1 or 2", call. = FALSE)
  if (local && is.null(maxdist) && is.null(neighb))
    stop("maxdist or neighb is required for local kriging", call. = FALSE)

  fixed_base <- if (is.null(fit$fixed)) {
    list()
  } else if (is.list(fit$fixed)) {
    fit$fixed
  } else {
    as.list(fit$fixed)
  }
  target_mean_name <- if (bivariate) paste0("mean_", which) else "mean"
  target_external_mean <- fixed_base[[target_mean_name]]
  target_external_mean <- !is.null(target_external_mean) &&
    length(target_external_mean) > 1L
  if (identical(type_krig, "Universal") && target_external_mean) {
    warning(
      "The target mean is supplied externally and is treated as known; ",
      "GeoCV uses Simple kriging for these folds.",
      call. = FALSE
    )
    type_krig <- "Simple"
  }
  mean_pattern <- if (bivariate) "^mean_[12][0-9]*$" else "^mean[0-9]*$"
  estimated_mean_names <- if (is.null(names(fit$param))) character(0) else
    grep(mean_pattern, names(fit$param), value = TRUE)
  if (identical(type_krig, "Universal") && !length(estimated_mean_names)) {
    warning(
      "No mean coefficients are estimated in fit; Universal and Simple kriging coincide. ",
      "GeoCV uses Simple kriging.",
      call. = FALSE
    )
    type_krig <- "Simple"
  }

  full_likelihood <- identical(fit$likelihood, "Full")
  if (identical(type_krig, "Universal") && isTRUE(estimation) && !full_likelihood) {
    stop(
      "Universal cross-validation with estimation=TRUE is not available for composite-likelihood refits, because each fold would require its own Godambe covariance matrix. Use type_krig='Simple', set estimation=FALSE with a fit already updated by GeoVarest(), or compute fold-specific Godambe matrices outside GeoCV.",
      call. = FALSE
    )
  }
  if (identical(type_krig, "Universal") && !isTRUE(estimation) && is.null(fit$varcov)) {
    stop(
      "Universal cross-validation with estimation=FALSE requires fit$varcov. For a composite-likelihood fit, update the fitted object with GeoVarest() first.",
      call. = FALSE
    )
  }

  ## Keep the data-generating model (fit$model) distinct from the criterion
  ## used for refitting.  Current GeoFit objects store estimation_model
  ## explicitly; the fallback is only for older misspecified fits.
  model1 <- fit$estimation_model
  if (is.null(model1) || length(model1) != 1L || is.na(model1) || !nzchar(model1)) {
    model1 <- as.character(fit$model)[1L]
    if (isTRUE(fit$missp)) {
      fallback <- c(
        StudentT = "Gaussian_misp_StudentT",
        Poisson = "Gaussian_misp_Poisson",
        PoissonZIP = "Gaussian_misp_PoissonZIP",
        PoissonGamma = "Gaussian_misp_PoissonGamma",
        Binomial = "Gaussian_misp_Binomial",
        SkewStudentT = "Gaussian_misp_SkewStudentT",
        Tukeygh = "Gaussian_misp_Tukeygh"
      )
      if (identical(model1, "BinomialNeg")) {
        stop(
          "GeoCV: this older misspecified BinomialNeg fit does not record the working likelihood. Refit it with the current GeoFit implementation before cross-validation.",
          call. = FALSE
        )
      }
      if (!(model1 %in% names(fallback))) {
        stop(
          "GeoCV: cannot recover the misspecified estimation model from this older GeoFit object. Refit it with the current GeoFit implementation.",
          call. = FALSE
        )
      }
      model1 <- unname(fallback[[model1]])
    }
  } else {
    model1 <- as.character(model1)[1L]
  }

  ## Exact probability scores below use a Gaussian predictive distribution.
  gaussian_scores <- identical(gsub("[[:blank:]]", "", fit$model), "Gaussian")
  if (!gaussian_scores) {
    warning(
      "For non-Gaussian fitted models GeoCV reports RMSE, MAE and MAD only; ",
      "Gaussian probability scores (log-score, Brier, CRPS, PIT, interval score and coverage) are returned as NA.",
      call. = FALSE
    )
  }

  if (is.null(optimizer)) optimizer <- fit$optimizer
  if (is.null(lower)) lower <- fit$lower
  if (is.null(upper)) upper <- fit$upper

  fit_method <- "cholesky"
  if (!is.null(fit$method)) {
    fit_method <- fit$method
  } else if (!is.null(fit$call)) {
    method_expr <- as.list(fit$call)$method
    if (!is.null(method_expr)) {
      method_value <- tryCatch(
        eval(method_expr, envir = parent.frame()),
        error = function(e) NULL
      )
      if (is.character(method_value) && length(method_value) == 1L)
        fit_method <- method_value
    }
  }
  fit_method <- tolower(gsub("[[:blank:]]", "", fit_method))
  if (!(fit_method %in% c("cholesky", "svd"))) fit_method <- "cholesky"

  refit_varest <- identical(type_krig, "Universal") && full_likelihood
  maxtime_use <- if (is.null(maxtime)) fit$maxtime else maxtime

  ## angle/ratio are represented once through anisopars.  GeoFit will append
  ## them to the estimated/fixed parameter blocks according to est.aniso.
  refit_start <- fit$param
  refit_fixed_base <- fixed_base
  if (!is.null(fit$anisopars)) {
    anames <- intersect(c("angle", "ratio"), names(unlist(fit$anisopars, use.names = TRUE)))
    if (length(anames)) {
      if (!is.null(refit_start)) refit_start[intersect(names(refit_start), anames)] <- NULL
      if (!is.null(refit_fixed_base)) refit_fixed_base[intersect(names(refit_fixed_base), anames)] <- NULL
      if (!is.null(lower) && !is.null(names(lower))) lower[intersect(names(lower), anames)] <- NULL
      if (!is.null(upper) && !is.null(names(upper))) upper[intersect(names(upper), anames)] <- NULL
    }
  }

  ## Normalize count-size information once; the remaining fold utilities are
  ## shared in CVUtils.R rather than recreated inside every GeoCV call.
  n_cv <- .GeoCVNormalizeCountSize(fit$n, fit$model)

  ## ---------- spatial case ----------
  results <- NULL
  if (space) {
    coords <- .GeoCVExplicitCoordinates(fit, "GeoCV spatial coordinates")
    N <- nrow(coords)
    data <- .GeoKrigNormalizeData(
      fit$data, nspace = N, coordt = NULL, grid = isTRUE(fit$grid),
      context = "GeoCV spatial data"
    )
    Xfull <- .GeoCVFlattenDesign(fit$X, N, "GeoCV spatial")
    Mfull <- .GeoCVVaryingFixedMean(fixed_base, "mean", N)
    nfull <- .GeoCVFlattenObservationVector(n_cv, N, "fit$n")
    test_samples <- .GeoCVMakeTestSamples(N, K, n.fold)

    cv_iteration <- function(i) {
      test_idx <- test_samples[[i]]
      train_idx <- setdiff(seq_len(N), test_idx)

      data_est <- data[train_idx]
      data_pred <- data[test_idx]
      coords_est <- coords[train_idx, , drop = FALSE]
      coords_pred <- coords[test_idx, , drop = FALSE]
      X_est <- if (!is.null(Xfull)) Xfull[train_idx, , drop = FALSE] else NULL
      X_pred <- if (!is.null(Xfull)) Xfull[test_idx, , drop = FALSE] else NULL
      M_pred <- if (!is.null(Mfull)) Mfull[test_idx] else NULL
      n_est <- .GeoCVSubsetObservationVector(nfull, train_idx)
      n_pred <- .GeoCVSubsetObservationVector(nfull, test_idx)

      fixed_for_est <- refit_fixed_base
      if (!is.null(Mfull)) fixed_for_est$mean <- Mfull[train_idx]

      if (estimation) {
        fit_s <- .GeoCVRefit(
          fit = fit, data = data_est, X = X_est, n = n_est,
          fixed = fixed_for_est, start = refit_start, model = model1,
          optimizer = optimizer, lower = lower, upper = upper,
          sparse = sparse, method = fit_method, varest = refit_varest,
          coordx = coords_est,
          universal = identical(type_krig, "Universal")
        )
      } else {
        fit_s <- fit
        fit_s$data <- data_est
        fit_s$coordx <- coords_est[, 1L]
        fit_s$coordy <- coords_est[, 2L]
        fit_s$coordz <- if (ncol(coords_est) == 3L) coords_est[, 3L] else NULL
        fit_s$coordx_dyn <- NULL
        fit_s$grid <- FALSE
        fit_s$ns <- length(train_idx)
        fit_s$X <- X_est
        fit_s$n <- n_est
        fit_s$fixed <- fixed_for_est
      }

      pr <- .GeoCVPredict(
        fit_s, loc = coords_pred, Xloc = X_pred, Mloc = M_pred, nloc = n_pred,
        local = local, neighb = neighb, maxdist = maxdist, maxtime = maxtime_use,
        type_krig = type_krig, sparse = sparse, method = fit_method
      )

      .GeoCVSummarizeScores(data_pred, pr$pred, pr$mse, gaussian_scores)
    }

    results <- .GeoCVRunIterations(cv_iteration, K, parallel, ncores, progress)
  }

  ## ---------- space-time case ----------
  if (spacetime) {
    coordt <- as.numeric(fit$coordt)
    if (!length(coordt) || any(!is.finite(coordt)))
      stop("fit$coordt is invalid", call. = FALSE)
    T <- length(coordt)

    if (!is.null(fit$coordx_dyn)) {
      coord_blocks <- .GeoNormalizeDynamicCoordinates(
        fit$coordx_dyn, coordt = coordt, context = "GeoCV space-time"
      )
      ns <- vapply(coord_blocks, nrow, integer(1))
      if (is.list(fit$data)) {
        if (length(fit$data) != T)
          stop("GeoCV: fit$data is not aligned with fit$coordt.", call. = FALSE)
        data_blocks <- lapply(seq_len(T), function(k) {
          z <- as.numeric(fit$data[[k]])
          if (length(z) != ns[k] || any(!is.finite(z)))
            stop("GeoCV: dynamic data block ", k, " is invalid.", call. = FALSE)
          z
        })
      } else {
        data_flat0 <- as.numeric(fit$data)
        if (length(data_flat0) != sum(ns) || any(!is.finite(data_flat0)))
          stop("GeoCV: dynamic data are not aligned with coordx_dyn.", call. = FALSE)
        data_blocks <- split(data_flat0, rep(seq_len(T), ns))
      }
    } else {
      coords_fixed <- .GeoCVExplicitCoordinates(fit, "GeoCV fixed space-time coordinates")
      ns <- rep(nrow(coords_fixed), T)
      data_flat0 <- .GeoKrigNormalizeData(
        fit$data, nspace = nrow(coords_fixed), coordt = coordt,
        grid = isTRUE(fit$grid), context = "GeoCV fixed space-time data"
      )
      data_blocks <- split(data_flat0, rep(seq_len(T), ns))
      coord_blocks <- lapply(seq_len(T), function(k) coords_fixed)
    }

    NT <- sum(ns)
    time_group <- rep(seq_len(T), ns)
    coords_flat <- do.call(rbind, coord_blocks)
    data_flat <- as.numeric(unlist(data_blocks, use.names = FALSE))
    Xflat <- .GeoCVFlattenDesign(fit$X, ns, "GeoCV space-time")
    Mflat <- .GeoCVVaryingFixedMean(fixed_base, "mean", NT)
    if (!is.null(Mflat) && !is.null(Xflat))
      stop("GeoCV: covariates and a fixed varying mean are not compatible.", call. = FALSE)
    nflat <- .GeoCVFlattenObservationVector(n_cv, ns, "fit$n")
    test_samples <- .GeoCVMakeTestSamples(NT, K, n.fold, groups = time_group)

    cv_iteration_st <- function(i) {
      test_idx <- test_samples[[i]]
      train_idx <- setdiff(seq_len(NT), test_idx)

      train_by_time <- lapply(seq_len(T), function(k)
        train_idx[time_group[train_idx] == k])
      test_by_time <- lapply(seq_len(T), function(k)
        test_idx[time_group[test_idx] == k])

      coord_train <- lapply(train_by_time, function(idx)
        coords_flat[idx, , drop = FALSE])
      data_train <- lapply(train_by_time, function(idx) data_flat[idx])
      X_train <- if (!is.null(Xflat)) {
        lapply(train_by_time, function(idx) Xflat[idx, , drop = FALSE])
      } else NULL
      M_train <- if (!is.null(Mflat)) Mflat[unlist(train_by_time)] else NULL
      n_train <- .GeoCVSubsetObservationVector(nflat, unlist(train_by_time))

      fixed_for_est <- refit_fixed_base
      if (!is.null(Mflat)) fixed_for_est$mean <- M_train

      if (estimation) {
        fit_s <- .GeoCVRefit(
          fit = fit, data = data_train, X = X_train, n = n_train,
          fixed = fixed_for_est, start = refit_start, model = model1,
          optimizer = optimizer, lower = lower, upper = upper,
          sparse = sparse, method = fit_method, varest = refit_varest,
          coordx_dyn = coord_train, coordt = coordt,
          universal = identical(type_krig, "Universal")
        )
      } else {
        fit_s <- fit
        fit_s$data <- data_train
        fit_s$coordx_dyn <- coord_train
        fit_s$coordt <- coordt
        fit_s$coordx <- unlist(lapply(coord_train, function(z) z[, 1L]))
        fit_s$coordy <- unlist(lapply(coord_train, function(z) z[, 2L]))
        fit_s$coordz <- if (ncol(coords_flat) == 3L) {
          unlist(lapply(coord_train, function(z) z[, 3L]))
        } else NULL
        fit_s$grid <- FALSE
        fit_s$ns <- vapply(coord_train, nrow, integer(1))
        fit_s$X <- X_train
        fit_s$n <- n_train
        fit_s$fixed <- fixed_for_est
      }

      pred_values <- pred_mse <- obs_values <- list()
      out_pos <- 0L
      for (k in seq_len(T)) {
        idx <- test_by_time[[k]]
        if (!length(idx)) next
        out_pos <- out_pos + 1L
        loc_k <- coords_flat[idx, , drop = FALSE]
        Xloc_k <- if (!is.null(Xflat)) Xflat[idx, , drop = FALSE] else NULL
        Mloc_k <- if (!is.null(Mflat)) Mflat[idx] else NULL
        nloc_k <- .GeoCVSubsetObservationVector(nflat, idx)

        pr <- .GeoCVPredict(
          fit_s, loc = loc_k, Xloc = Xloc_k, Mloc = Mloc_k, nloc = nloc_k,
          local = local, neighb = neighb, maxdist = maxdist, maxtime = maxtime_use,
          type_krig = type_krig, sparse = sparse, method = fit_method,
          time = coordt[k]
        )
        pred_values[[out_pos]] <- pr$pred
        pred_mse[[out_pos]] <- if (is.null(pr$mse)) {
          rep(NA_real_, length(pr$pred))
        } else pr$mse
        obs_values[[out_pos]] <- data_flat[idx]
      }

      .GeoCVSummarizeScores(
        unlist(obs_values, use.names = FALSE),
        unlist(pred_values, use.names = FALSE),
        unlist(pred_mse, use.names = FALSE),
        gaussian_scores
      )
    }

    results <- .GeoCVRunIterations(cv_iteration_st, K, parallel, ncores, progress)
  }

  ## ---------- bivariate case ----------
  if (bivariate) {
    ns <- as.integer(fit$ns)
    if (length(ns) != 2L || any(!is.finite(ns)) || any(ns < 1L))
      stop("GeoCV: invalid bivariate block sizes.", call. = FALSE)
    if (ns[which] < 2L)
      stop("GeoCV requires at least two observations in the target bivariate component.", call. = FALSE)

    if (!is.null(fit$coordx_dyn)) {
      if (length(fit$coordx_dyn) != 2L)
        stop("GeoCV: bivariate coordx_dyn must contain two coordinate matrices.", call. = FALSE)
      coords1 <- as.matrix(fit$coordx_dyn[[1L]])
      coords2 <- as.matrix(fit$coordx_dyn[[2L]])
    } else {
      coords_common <- .GeoCVExplicitCoordinates(fit, "GeoCV bivariate coordinates")
      coords1 <- coords2 <- coords_common
    }
    if (nrow(coords1) != ns[1L] || nrow(coords2) != ns[2L])
      stop("GeoCV: bivariate coordinates are not aligned with fit$ns.", call. = FALSE)

    biv_data <- .GeoBivariateDataVector(
      fit$data, ns, context = "GeoCV bivariate"
    )
    data1 <- biv_data[seq_len(ns[1L])]
    data2 <- biv_data[ns[1L] + seq_len(ns[2L])]

    Xflat <- if (!is.null(fit$X)) {
      .GeoBivariateDesign(fit$X, ns, context = "GeoCV bivariate")
    } else NULL
    X1 <- if (!is.null(Xflat)) Xflat[seq_len(ns[1L]), , drop = FALSE] else NULL
    X2 <- if (!is.null(Xflat)) Xflat[ns[1L] + seq_len(ns[2L]), , drop = FALSE] else NULL

    M1 <- .GeoCVVaryingFixedMean(fixed_base, "mean_1", ns[1L])
    M2 <- .GeoCVVaryingFixedMean(fixed_base, "mean_2", ns[2L])
    if ((!is.null(M1) || !is.null(M2)) && !is.null(Xflat))
      stop("GeoCV: bivariate covariates and fixed varying means are not compatible.", call. = FALSE)

    ## GeoKrig's bivariate path uses coefficient-based observed means.  Legacy
    ## fitted objects may nevertheless contain site-specific mean_1 and/or
    ## mean_2 vectors.  Handle each component independently by centering that
    ## component and fixing its coefficient mean to zero; Mloc adds the known
    ## target mean back at prediction sites.
    data1_model <- if (!is.null(M1)) data1 - M1 else data1
    data2_model <- if (!is.null(M2)) data2 - M2 else data2
    fixed_biv <- refit_fixed_base
    start_biv <- refit_start
    if (!is.null(M1)) {
      fixed_biv$mean_1 <- 0
      start_biv <- .GeoCVDropMeanParameters(start_biv, "mean_1")
    }
    if (!is.null(M2)) {
      fixed_biv$mean_2 <- 0
      start_biv <- .GeoCVDropMeanParameters(start_biv, "mean_2")
    }

    nflat <- .GeoCVFlattenObservationVector(n_cv, ns, "fit$n")
    n1 <- if (length(nflat) == 1L) nflat else nflat[seq_len(ns[1L])]
    n2 <- if (length(nflat) == 1L) nflat else nflat[ns[1L] + seq_len(ns[2L])]

    target_n <- ns[which]
    test_samples <- .GeoCVMakeTestSamples(target_n, K, n.fold)

    cv_iteration_biv <- function(i) {
      test_idx <- test_samples[[i]]
      train_idx <- setdiff(seq_len(target_n), test_idx)
      current_fixed <- fixed_biv

      if (which == 1L) {
        data_est <- list(data1_model[train_idx], data2_model)
        data_pred <- data1[test_idx]
        coords_est <- list(coords1[train_idx, , drop = FALSE], coords2)
        coords_pred <- coords1[test_idx, , drop = FALSE]
        X_est <- if (!is.null(Xflat)) {
          rbind(X1[train_idx, , drop = FALSE], X2)
        } else NULL
        X_pred <- if (!is.null(X1)) X1[test_idx, , drop = FALSE] else NULL
        M_pred <- if (!is.null(M1)) M1[test_idx] else NULL
        n_est <- if (length(nflat) == 1L) nflat else c(n1[train_idx], n2)
        n_pred <- if (length(nflat) == 1L) nflat else n1[test_idx]
      } else {
        data_est <- list(data1_model, data2_model[train_idx])
        data_pred <- data2[test_idx]
        coords_est <- list(coords1, coords2[train_idx, , drop = FALSE])
        coords_pred <- coords2[test_idx, , drop = FALSE]
        X_est <- if (!is.null(Xflat)) {
          rbind(X1, X2[train_idx, , drop = FALSE])
        } else NULL
        X_pred <- if (!is.null(X2)) X2[test_idx, , drop = FALSE] else NULL
        M_pred <- if (!is.null(M2)) M2[test_idx] else NULL
        n_est <- if (length(nflat) == 1L) nflat else c(n1, n2[train_idx])
        n_pred <- if (length(nflat) == 1L) nflat else n2[test_idx]
      }

      if (estimation) {
        fit_s <- .GeoCVRefit(
          fit = fit, data = data_est, X = X_est, n = n_est,
          fixed = current_fixed, start = start_biv, model = model1,
          optimizer = optimizer, lower = lower, upper = upper,
          sparse = sparse, method = fit_method, varest = refit_varest,
          coordx_dyn = coords_est,
          universal = identical(type_krig, "Universal")
        )
      } else {
        fit_s <- fit
        fit_s$data <- data_est
        fit_s$coordx_dyn <- coords_est
        fit_s$coordx <- fit_s$coordy <- fit_s$coordz <- NULL
        fit_s$grid <- FALSE
        fit_s$ns <- vapply(coords_est, nrow, integer(1))
        fit_s$X <- X_est
        fit_s$n <- n_est
        fit_s$fixed <- current_fixed
      }

      pr <- .GeoCVPredict(
        fit_s, loc = coords_pred, Xloc = X_pred, Mloc = M_pred, nloc = n_pred,
        local = local, neighb = neighb, maxdist = maxdist, maxtime = maxtime_use,
        type_krig = type_krig, sparse = sparse, method = fit_method,
        which = which
      )

      .GeoCVSummarizeScores(data_pred, pr$pred, pr$mse, gaussian_scores)
    }

    results <- .GeoCVRunIterations(cv_iteration_biv, K, parallel, ncores, progress)
  }

  if (is.null(results))
    stop("Unsupported GeoFit configuration for GeoCV.", call. = FALSE)

  list(
    rmse = results[, "rmse"],
    mae = results[, "mae"],
    mad = results[, "mad"],
    brie = results[, "brie"],
    crps = results[, "crps"],
    lscore = results[, "lscore"],
    pit_mean = results[, "pit_mean"],
    pit = results[, "pit_mean"],
    intscore = results[, "intscore"],
    coverage = results[, "coverage"],
    seed = seed,
    estimation_model = model1
  )
}
