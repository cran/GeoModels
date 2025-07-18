GeoVarestbootstrap <- function(fit, K = 100, sparse = FALSE, optimizer = NULL,
                               lower = NULL, upper = NULL, method = "cholesky",
                               alpha = 0.95, L = 1000, parallel = TRUE,
                               ncores = NULL) {

  ## -------- 0. VALIDAZIONI INIZIALI (inalterate) --------
  if (length(fit$coordt) == 1) fit$coordt <- NULL
  if (is.null(fit$sensmat))
    stop("Sensitivity matrix is missing: use sensitivity=TRUE in GeoFit")
  if (!(method %in% c("cholesky", "TB", "CE")))
    stop("The method of simulation is not correct")
  if (is.numeric(alpha) && !(alpha > 0 && alpha < 1))
    stop("alpha must be numeric between 0 and 1")
  if (!is.numeric(K) || K < 1) stop("K must be a positive integer")

  if (is.null(optimizer)) {
    optimizer <- fit$optimizer
    lower     <- fit$lower
    upper     <- fit$upper
  }

  future.seed = TRUE

  model <- fit$model
  cat("Parametric bootstrap can be time consuming ...\n")

  if (fit$missp) {
    model_map <- c(
      StudentT = "Gaussian_misp_StudentT",
      Poisson  = "Gaussian_misp_Poisson",
      PoissonZIP = "Gaussian_misp_PoissonZIP",
      SkewStudentT = "Gaussian_misp_SkewStudentT",
      Tukeygh = "Gaussian_misp_Tukeygh"
    )
    if (model %in% names(model_map)) model <- model_map[[model]]
  }

  dimat <- fit$numtime * fit$numcoord
  X_use <- if (!is.null(fit$X) &&
               sum(fit$X[1:dimat] == 1) == dimat &&
               ncol(fit$X) == 1) NULL else fit$X

  coords <- cbind(fit$coordx, fit$coordy)
  if (fit$bivariate && is.null(fit$coordx_dyn))
    coords <- coords[seq_len(nrow(coords) / 2), , drop = FALSE]

  ## -------- 1. SIMULAZIONE DATI --------
  sim_args <- list(
    coordx = coords, coordt = fit$coordt, coordx_dyn = fit$coordx_dyn,
    anisopars = fit$anisopars, corrmodel = fit$corrmodel, model = model,
    param = append(fit$param, fit$fixed), grid = fit$grid, X = fit$X,
    n = fit$n, distance = fit$distance, radius = fit$radius, nrep = K,
    progress = FALSE
  )

  if (is.null(fit$copula)) {
    if (method == "cholesky") {
      data_sim <- do.call(GeoSim, c(sim_args, list(sparse = sparse)))
    } else if (method %in% c("TB", "CE")) {
      data_sim <- do.call(GeoSimapprox,
                          c(sim_args,
                            list(method = method, parallel = parallel,
                                 ncores = ncores, L = L)))
    } else {
      stop("Unsupported method for simulation")
    }
  } else {
    if (method == "cholesky") {
      data_sim <- do.call(
        GeoSimCopula,
        c(sim_args,
          list(copula = fit$copula, sparse = sparse))
      )
    } else {
      stop("Unsupported method for copula simulation")
    }
  }

  ## -------- 2. IMPOSTAZIONI PARALLELE --------
  coremax <- parallel::detectCores()
  if (is.na(coremax) || coremax == 1) parallel <- FALSE
  if (parallel && is.null(ncores))
    ncores <- max(1, min(coremax - 1, K))
  ncores <- max(1, min(ncores, K))

  ## -------- 3. FUNZIONE DI STIMA --------
  estimate_fun <- function(k) {
    res <- tryCatch(
      GeoFit(
        data = data_sim$data[[k]], start = fit$param, fixed = fit$fixed,
        coordx = coords, coordt = fit$coordt, coordx_dyn = fit$coordx_dyn,
        copula = fit$copula, anisopars = fit$anisopars,
        est.aniso = fit$est.aniso, lower = lower, upper = upper,
        neighb = fit$neighb, corrmodel = fit$corrmodel, model = model,
        sparse = FALSE, n = fit$n, maxdist = fit$maxdist,
        maxtime = fit$maxtime, optimizer = optimizer,
        grid = fit$grid, likelihood = fit$likelihood,
        type = fit$type, X = X_use, distance = fit$distance,
        radius = fit$radius
      ),
      error = function(e) NULL
    )

    if (!is.null(res) &&
        identical(res$convergence, "Successful") &&
        is.finite(res$logCompLik) && res$logCompLik < 1e8) {
      unlist(res$param)
    } else {
      NULL
    }
  }

  ## -------- 4. ESECUZIONE con barra “intelligente” --------
  if (!parallel) {
    cat("Performing", K, "estimations sequentially...\n")
    res_list <- lapply(seq_len(K), estimate_fun)
  } else {
    future::plan(future::multisession, workers = ncores)
    on.exit(future::plan(future::sequential), add = TRUE)

    # stima grossolana per decidere se mostrare la barra
    time_one <- system.time(estimate_fun(1L))[["elapsed"]]
    use_bar  <- time_one > 5 && K >= 100

    if (use_bar && requireNamespace("progressr", quietly = TRUE)) {
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(K)
      res_list <- future.apply::future_lapply(
        seq_len(K),
        function(k) { pb(); estimate_fun(k) },
        future.seed   = future.seed,
        future.stdout = FALSE,
        future.conditions = "none"
      )
    } else {
      res_list <- future.apply::future_lapply(
        seq_len(K), estimate_fun,
        future.seed   = future.seed,
        future.stdout = FALSE,
        future.conditions = "none"
      )
    }
  }
  res <- do.call(rbind, Filter(Negate(is.null), res_list))

  if (nrow(res) < 2)
    stop("Insufficient successful bootstrap iterations for variance estimation")
  cat("Successful bootstrap iterations:", nrow(res), "out of", K, "\n")

  ## -------- 5. POST-PROCESSING (invariato) --------
  numparam <- length(fit$param)
  invG <- var(res)
  G <- tryCatch(solve(invG),
                error = function(e) {
                  warning("Bootstrap estimated Godambe matrix is singular")
                  NULL
                })

  stderr <- sqrt(pmax(0, diag(invG)))
  names(stderr) <- names(fit$param)

  if ((fit$likelihood == "Marginal" && fit$type %in% c("Independence", "Pairwise")) ||
      (fit$likelihood == "Conditional" && fit$type == "Pairwise")) {
    H <- fit$sensmat
    penalty <- sum(diag(H %*% invG))
    claic <- -2 * fit$logCompLik + 2 * penalty
    clbic <- -2 * fit$logCompLik + log(dimat) * penalty
    fit$varimat <- H %*% invG %*% H
  } else if (fit$likelihood == "Full" && fit$type == "Standard") {
    claic <- -2 * fit$logCompLik + 2 * numparam
    clbic <- -2 * fit$logCompLik + log(dimat) * numparam
  } else {
    claic <- clbic <- NA_real_
  }

  fit$claic      <- claic
  fit$clbic      <- clbic
  fit$stderr     <- stderr
  fit$varcov     <- invG
  fit$estimates  <- apply(res, 2, as.numeric)

  if (all(stderr > 0)) {
    aa <- qnorm(1 - (1 - alpha) / 2) * stderr
    pp <- as.numeric(fit$param)
    names(pp) <- names(fit$param)
    fit$CI <- data.frame(
      Lower = pp - aa,
      Upper = pp + aa,
      row.names = names(fit$param)
    )
    fit$pvalues <- 2 * pnorm(-abs(pp / stderr))
  } else {
    fit$CI <- fit$pvalues <- NULL
  }

  fit
}