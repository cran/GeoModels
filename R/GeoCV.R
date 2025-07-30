GeoCV <- function(fit, K = 100, estimation = TRUE,
                  optimizer = NULL, lower = NULL, upper = NULL,
                  n.fold = 0.05, local = FALSE, neighb = NULL, maxdist = NULL,
                  maxtime = NULL, sparse = FALSE, type_krig = "Simple",
                  which = 1, parallel = TRUE, ncores = NULL,
                  progress = TRUE) {

  ## ---------- sanity checks ----------
  if (!inherits(fit, "GeoFit"))
    stop("fit must be an object of class 'GeoFit'")
  if (n.fold > 0.99 || n.fold < 0.01)
    stop("n.fold must be between 0.01 and 0.99")
  if (!is.logical(parallel))
    stop("parallel must be TRUE or FALSE")
  if (!is.logical(progress))
    stop("progress must be TRUE or FALSE")

  ## Se l'utente non ha progressr installato, disattiviamo il progresso
  has_progressr <- requireNamespace("progressr", quietly = TRUE)
  if (parallel && progress && !has_progressr) {
    warning("Package 'progressr' is required for progress bars in parallel mode; disabling progress")
    progress <- FALSE
  }

  message("Cross-validation kriging can be time consuming ...")

  ## ---------- utilities ----------
  choose_cores <- function(nc) {
    coremax <- parallel::detectCores()
    if (is.na(coremax) || coremax == 1) return(1L)
    if (is.null(nc)) return(max(1L, coremax - 1L))
    if (is.numeric(nc) && nc >= 1 && nc <= coremax) return(as.integer(nc))
    stop("ncores not valid")
  }

  ## ---------- setup ----------
  if (is.null(fit$X)) {
    X <- Xloc <- NULL
    tempX <- NULL
  }
  mae <- rmse <- lscore <- crps <- mad <- brie <- NULL
  space_dyn <- FALSE
  if (is.null(optimizer)) {
    optimizer <- fit$optimizer
    lower     <- fit$lower
    upper     <- fit$upper
  }

  if (is.list(fit$data)) space_dyn <- TRUE
  spacetime <- CheckST(CkCorrModel(fit$corrmodel))
  bivariate <- CheckBiv(CkCorrModel(fit$corrmodel))
  K <- round(K)
  if (K < 2) stop("K must be >= 2")
  if (K > 1000) stop("K is too large")
  if (bivariate && !(which == 1 || which == 2)) stop("which must be 1 or 2")
  if (local && is.null(maxdist) && is.null(neighb))
    stop("maxdist or neighb required for local kriging")

  model1 <- fit$model
  if (fit$missp) {
    model1 <- switch(
      fit$model,
      StudentT     = "Gaussian_misp_StudentT",
      Poisson      = "Gaussian_misp_Poisson",
      PoissonZIP   = "Gaussian_misp_PoissonZIP",
      SkewStudentT = "Gaussian_misp_SkewStudenT",
      Tukeygh      = "Gaussian_misp_Tukeygh"
    )
  }

  space <- !spacetime && !bivariate

  ################################################################
  ### Spatial case  ##############################################
  ################################################################
  if (space) {
    N      <- length(fit$data)
    coords <- cbind(fit$coordx, fit$coordy, fit$coordz)
    data   <- fit$data
    if (length(fit$fixed$mean) > 1) tempM <- fit$fixed$mean
    if (!is.null(fit$X)) tempX <- fit$X

    sel_list <- replicate(K, sample(1:N, round(N * (1 - n.fold))), simplify = FALSE)

    cv_iteration <- function(i) {
      sel_data <- sel_list[[i]]
      X <- Xloc <- Mloc <- NULL
      if (!is.null(fit$X)) {
        X    <- tempX[sel_data, , drop = FALSE]
        Xloc <- tempX[-sel_data, , drop = FALSE]
      }
      if (length(fit$fixed$mean) > 1) {
        fit$fixed$mean <- tempM[sel_data]
        Mloc <- tempM[-sel_data]
      }

      data_to_pred <- data[-sel_data]
      data_to_est  <- data[sel_data]
      coords_est   <- coords[sel_data, , drop = FALSE]
      coords_pred  <- coords[-sel_data, , drop = FALSE]

      param <- c(fit$param, fit$fixed)
      if (estimation) {
        fit_s <- suppressWarnings(
          GeoFit(
            data = data_to_est, coordx = coords_est, corrmodel = fit$corrmodel,
            X = X, likelihood = fit$likelihood, type = fit$type, grid = fit$grid,
            copula = fit$copula, anisopars = fit$anisopars, est.aniso = fit$est.aniso,
            model = model1, radius = fit$radius, n = fit$n,
            maxdist = fit$maxdist, neighb = fit$neighb, distance = fit$distance,
            optimizer = optimizer, lower = lower, upper = upper,
            start = fit$param, fixed = fit$fixed
          )
        )
        param <- c(fit_s$param, fit_s$fixed)
      } else {
        fit_s <- fit
      }

      pr <- if (!local) {
        GeoKrig(fit_s, loc = coords_pred, mse = TRUE,
                param = param, Xloc = Xloc, Mloc = Mloc)
      } else {
        GeoKrigloc(fit_s, loc = coords_pred, mse = TRUE,
                   param = param, Xloc = Xloc, Mloc = Mloc,
                   neighb = neighb, maxdist = maxdist, progress = FALSE)
      }

      pp <- GeoScores(data_to_pred, pred = pr$pred, mse = pr$mse,
                      score = c("brie", "crps", "lscore", "pe"))
      c(pp$rmse, pp$mae, pp$mad, pp$lscore, pp$brie, pp$crps)
    }

    ## ---- esecuzione effettiva ---------------------------------
    if (!parallel) {
      rmse <- crps <- mae <- mad <- lscore <- brie <- double(K)
      if (isTRUE(progress)) cat("Performing", K, "cross-validations...\n")
      if (isTRUE(progress)) {
        pb <- txtProgressBar(min = 0, max = K, style = 3)
        on.exit(close(pb), add = TRUE)
      }
      for (i in 1:K) {
        res <- cv_iteration(i)
        rmse[i] <- res[1]; mae[i] <- res[2]; mad[i] <- res[3]
        lscore[i] <- res[4]; brie[i] <- res[5]; crps[i] <- res[6]
        if (isTRUE(progress)) setTxtProgressBar(pb, i)
      }
    } else {
      n.cores <- choose_cores(ncores)
      future::plan(future::multisession, workers = n.cores)
      on.exit(future::plan(future::sequential), add = TRUE)

      if (isTRUE(progress)) {
        p <- progressr::progressor(along = 1:K)
        res <- future.apply::future_lapply(
          1:K,
          function(i) { out <- cv_iteration(i); p(); out },
          future.seed = TRUE,
          future.stdout = FALSE,
          future.conditions = "none"
        )
      } else {
        res <- future.apply::future_lapply(
          1:K, cv_iteration,
          future.seed = TRUE,
          future.stdout = FALSE,
          future.conditions = "none"
        )
      }
      res <- do.call(rbind, res)
      rmse <- res[, 1]; mae <- res[, 2]; mad <- res[, 3]
      lscore <- res[, 4]; brie <- res[, 5]; crps <- res[, 6]
    }
  }

  ################################################################
  ### Spatio-temporal case  ######################################
  ################################################################
  if (spacetime) {
    coords <- cbind(fit$coordx, fit$coordy, fit$coordz)
    ns     <- fit$ns
    coordt <- fit$coordt
    T      <- length(coordt)
    NT     <- sum(ns)

    X <- if (!is.null(fit$X)) fit$X else NULL
    NS <- cumsum(ns)
    NS <- c(0, NS[-length(NS)])

    datos <- if (is.list(fit$data)) do.call(c, fit$data) else fit$data

    if (!space_dyn) {
      data_tot <- NULL
      for (k in 1:T) {
        data_tot <- rbind(data_tot, cbind(rep(coordt[k], ns[k]), coords, datos[k, ]))
      }
      if (!is.null(X)) data_tot <- cbind(data_tot, X)
    } else {
      ct        <- rep(coordt, times = ns)
      data_tot  <- cbind(ct, fit$coordx, fit$coordy, fit$coordz, datos)
      if (!is.null(X)) data_tot <- cbind(data_tot, do.call(rbind, X))
    }

    MM <- length(fit$fixed$mean) > 1
    if (MM && !is.null(X)) stop("covariates and fixed varying mean are not compatible")
    if (MM) {
      tempM <- if (!space_dyn) fit$fixed$mean else do.call(rbind, fit$fixed$mean)
      data_tot <- cbind(data_tot, tempM)
    }

    folds <- split(sample(1:NT), rep(1:K, length.out = NT))

    cv_iteration_st <- function(i) {
      test_idx  <- folds[[i]]
      train_idx <- setdiff(1:NT, test_idx)
      data_sel  <- as.matrix(data_tot[train_idx, , drop = FALSE])
      data_pred <- as.matrix(data_tot[test_idx, , drop = FALSE])

      data_sel_ord  <- data_sel[order(data_sel[, 1]), ]
      data_pred_ord <- data_pred[order(data_pred[, 1]), ]

      DD <- ncol(data_sel_ord)

      utt   <- unique(data_sel_ord[, 1])
      utt_1 <- unique(data_pred_ord[, 1])

      coordx_dynnew     <- Xnew <- datanew <- list()
      coordx_dynnew_loc <- Xnew_loc <- Mnew_loc <- list()

      for (k in seq_along(utt_1)) {
        ll <- data_pred_ord[data_pred_ord[, 1] == utt_1[k], , drop = FALSE]
        if (!is.null(fit$coordz)) {
          coordx_dynnew_loc[[k]] <- matrix(ll[, 2:4], ncol = 3)
          if (!MM && !is.null(X)) Xnew_loc[[k]] <- matrix(ll[, 6:DD], ncol = DD - 5)
          if (MM) Mnew_loc[[k]] <- ll[, 6]
        } else {
          coordx_dynnew_loc[[k]] <- matrix(ll[, 2:3], ncol = 2)
          if (!MM && !is.null(X)) Xnew_loc[[k]] <- matrix(ll[, 5:DD], ncol = DD - 4)
          if (MM) Mnew_loc[[k]] <- ll[, 5 + !is.null(X)]
        }
      }

      for (k in seq_along(utt)) {
        ss <- data_sel_ord[data_sel_ord[, 1] == utt[k], , drop = FALSE]
        if (!is.null(fit$coordz)) {
          coordx_dynnew[[k]] <- matrix(ss[, 2:4], ncol = 3)
          datanew[[k]]       <- as.vector(ss[, 5])
          if (!is.null(X)) Xnew[[k]] <- matrix(ss[, 6:DD], ncol = DD - 5)
        } else {
          coordx_dynnew[[k]] <- matrix(ss[, 2:3], ncol = 2)
          datanew[[k]]       <- as.vector(ss[, 4])
          if (!is.null(X)) Xnew[[k]] <- matrix(ss[, 5:DD], ncol = DD - 4)
        }
      }

      param <- c(fit$param, fit$fixed)
      if (estimation) {
        fit_s <- suppressWarnings(
          GeoFit(
            data = datanew, coordx_dyn = coordx_dynnew, coordt = utt,
            corrmodel = fit$corrmodel, X = Xnew, likelihood = fit$likelihood,
            type = fit$type, grid = fit$grid, copula = fit$copula,
            anisopars = fit$anisopars, est.aniso = fit$est.aniso,
            model = model1, radius = fit$radius, n = fit$n,
            maxdist = fit$maxdist, neighb = fit$neighb, maxtime = fit$maxtime,
            distance = fit$distance, optimizer = optimizer,
            lower = lower, upper = upper,
            start = fit$param, fixed = fit$fixed
          )
        )
        param <- c(fit_s$param, fit_s$fixed)
      } else {
        fit_s <- fit
      }

      pr_st  <- list()
      pr_mse <- list()
      for (j in seq_along(utt_1)) {
        pr <- if (!local) {
          GeoKrig(fit_s, loc = coordx_dynnew_loc[[j]], time = utt_1[j], mse = TRUE,
                  param = param, Xloc = Xnew_loc[[j]], Mloc = Mnew_loc[[j]])
        } else {
          GeoKrigloc(fit_s, loc = coordx_dynnew_loc[[j]], time = utt_1[j], mse = TRUE,
                     param = param, Xloc = Xnew_loc[[j]], Mloc = Mnew_loc[[j]],
                     neighb = neighb, maxdist = maxdist, maxtime = maxtime, progress = FALSE)
        }
        pr_st[[j]]  <- pr$pred
        pr_mse[[j]] <- pr$mse
      }

      pp <- GeoScores(data_pred[, 4], pred = unlist(pr_st), mse = unlist(pr_mse),
                      score = c("brie", "crps", "lscore", "pe"))
      c(pp$rmse, pp$mae, pp$mad, pp$lscore, pp$brie, pp$crps)
    }

    ## ---- esecuzione effettiva ---------------------------------
    if (!parallel) {
      rmse <- crps <- mae <- mad <- lscore <- brie <- double(K)
      if (isTRUE(progress)) cat("Performing", K, "cross-validations...\n")
      if (isTRUE(progress)) {
        pb <- txtProgressBar(min = 0, max = K, style = 3)
        on.exit(close(pb), add = TRUE)
      }
      for (i in 1:K) {
        res <- cv_iteration_st(i)
        rmse[i] <- res[1]; mae[i] <- res[2]; mad[i] <- res[3]
        lscore[i] <- res[4]; brie[i] <- res[5]; crps[i] <- res[6]
        if (isTRUE(progress)) setTxtProgressBar(pb, i)
      }
    } else {
      n.cores <- choose_cores(ncores)
      future::plan(future::multisession, workers = n.cores)
      on.exit(future::plan(future::sequential), add = TRUE)

      if (isTRUE(progress)) {
        p <- progressr::progressor(along = 1:K)
        res <- future.apply::future_lapply(
          1:K,
          function(i) { out <- cv_iteration_st(i); p(); out },
          future.seed = TRUE,
          future.stdout = FALSE,
          future.conditions = "none"
        )
      } else {
        res <- future.apply::future_lapply(
          1:K, cv_iteration_st,
          future.seed = TRUE,
          future.stdout = FALSE,
          future.conditions = "none"
        )
      }
      res <- do.call(rbind, res)
      rmse <- res[, 1]; mae <- res[, 2]; mad <- res[, 3]
      lscore <- res[, 4]; brie <- res[, 5]; crps <- res[, 6]
    }
  }

  ################################################################
  ### Bivariate case  ############################################
  ################################################################
  if (bivariate) {
    fixmeans <- length(fit$fixed$mean_1) > 1 && length(fit$fixed$mean_2) > 1
    if (fixmeans) {
      tempM1 <- fit$fixed$mean_1
      tempM2 <- fit$fixed$mean_2
    }

    ns <- fit$ns
    if (space_dyn) {
      data1  <- fit$data[[1]]; data2  <- fit$data[[2]]
      coords1 <- fit$coordx_dyn[[1]]; coords2 <- fit$coordx_dyn[[2]]
    } else {
      data1  <- fit$data[1, ]; data2  <- fit$data[2, ]
      coords <- cbind(fit$coordx, fit$coordy, fit$coordz)
      coords1 <- coords2 <- coords
    }

    X1 <- X2 <- NULL
    if (!is.null(fit$X)) {
      if (!space_dyn && !is.list(fit$X)) {
        X1 <- fit$X[1:ns[1], , drop = FALSE]
        X2 <- fit$X[(ns[1] + 1):(ns[1] + ns[2]), , drop = FALSE]
      } else if (space_dyn && is.list(fit$X)) {
        X1 <- fit$X[[1]]; X2 <- fit$X[[2]]
      }
    }

    samples <- lapply(1:K, function(i) {
      if (which == 1) sample(1:ns[1], round(ns[1] * (1 - n.fold)))
      else sample(1:ns[2], round(ns[2] * (1 - n.fold)))
    })

    cv_iteration_biv <- function(i) {
      sel_data <- samples[[i]]
      if (which == 1) {
        data_to_est <- list(data1[sel_data], data2)
        data_to_pred <- data1[-sel_data]
        coords_est   <- list(coords1[sel_data, , drop = FALSE], coords2)
        coords_pred  <- coords1[-sel_data, , drop = FALSE]
        X_use   <- if (!is.null(X1)) rbind(X1[sel_data, , drop = FALSE], X2) else NULL
        Xloc_use <- if (!is.null(X1)) rbind(X1[-sel_data, , drop = FALSE], X2) else NULL
        if (fixmeans) {
          current_fixed <- fit$fixed
          current_fixed$mean_1 <- tempM1[sel_data]
          Mloc <- list(tempM1[-sel_data], tempM2)
        } else {
          current_fixed <- fit$fixed
          Mloc <- NULL
        }
      } else {
        data_to_est <- list(data1, data2[sel_data])
        data_to_pred <- data2[-sel_data]
        coords_est   <- list(coords1, coords2[sel_data, , drop = FALSE])
        coords_pred  <- coords2[-sel_data, , drop = FALSE]
        X_use   <- if (!is.null(X2)) rbind(X1, X2[sel_data, , drop = FALSE]) else NULL
        Xloc_use <- if (!is.null(X2)) rbind(X1, X2[-sel_data, , drop = FALSE]) else NULL
        if (fixmeans) {
          current_fixed <- fit$fixed
          current_fixed$mean_2 <- tempM2[sel_data]
          Mloc <- list(tempM1, tempM2[-sel_data])
        } else {
          current_fixed <- fit$fixed
          Mloc <- NULL
        }
      }

      param <- if (estimation) {
        fit_s <- suppressWarnings(
          GeoFit(
            data = data_to_est, coordx = NULL, coordx_dyn = coords_est,
            corrmodel = fit$corrmodel, X = X_use,
            likelihood = fit$likelihood, type = fit$type, grid = fit$grid,
            model = "Gaussian", radius = fit$radius, n = fit$n,
            copula = fit$copula,
            maxdist = fit$maxdist, neighb = fit$neighb, distance = fit$distance,
            optimizer = optimizer, lower = lower, upper = upper,
            start = fit$param, fixed = current_fixed
          )
        )
        append(fit_s$param, fit_s$fixed)
      } else {
        append(fit$param, current_fixed)
      }

      pr <- if (!local) {
        GeoKrig(fit_s, loc = coords_pred, mse = TRUE,
                param = param, Xloc = Xloc_use, Mloc = Mloc, which = which)
      } else {
        GeoKrigloc(fit_s, loc = coords_pred, mse = TRUE,
                   neighb = neighb, maxdist = maxdist,
                   param = param, Xloc = Xloc_use, Mloc = Mloc, which = which,
                   progress = FALSE)
      }

      pp <- GeoScores(data_to_pred, pred = pr$pred, mse = pr$mse,
                      score = c("brie", "crps", "lscore", "pe"))
      c(pp$rmse, pp$mae, pp$mad, pp$lscore, pp$brie, pp$crps)
    }

    ## ---- esecuzione effettiva ---------------------------------
    if (!parallel) {
      rmse <- crps <- mae <- mad <- lscore <- brie <- double(K)
      if (isTRUE(progress)) cat("Performing", K, "cross-validations...\n")
      if (isTRUE(progress)) {
        pb <- txtProgressBar(min = 0, max = K, style = 3)
        on.exit(close(pb), add = TRUE)
      }
      for (i in 1:K) {
        res <- cv_iteration_biv(i)
        rmse[i] <- res[1]; mae[i] <- res[2]; mad[i] <- res[3]
        lscore[i] <- res[4]; brie[i] <- res[5]; crps[i] <- res[6]
        if (isTRUE(progress)) setTxtProgressBar(pb, i)
      }
    } else {
      n.cores <- choose_cores(ncores)
      future::plan(future::multisession, workers = n.cores)
      on.exit(future::plan(future::sequential), add = TRUE)

      if (isTRUE(progress)) {
        p <- progressr::progressor(along = 1:K)
        res <- future.apply::future_lapply(
          1:K,
          function(i) { out <- cv_iteration_biv(i); p(); out },
          future.seed = TRUE,
          future.stdout = FALSE,
          future.conditions = "none"
        )
      } else {
        res <- future.apply::future_lapply(
          1:K, cv_iteration_biv,
          future.seed = TRUE,
          future.stdout = FALSE,
          future.conditions = "none"
        )
      }
      res <- do.call(rbind, res)
      rmse <- res[, 1]; mae <- res[, 2]; mad <- res[, 3]
      lscore <- res[, 4]; brie <- res[, 5]; crps <- res[, 6]
    }
  }

  list(rmse = rmse, mae = mae, mad = mad, brie = brie, crps = crps, lscore = lscore)
}