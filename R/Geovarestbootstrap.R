GeoVarestbootstrap <- function(fit,K = 100,sparse = FALSE,
                               optimizer = NULL,lower = NULL,upper = NULL,method = "cholesky",alpha = 0.95, 
                               L = 1000,parallel = TRUE, ncores = NULL, progress = TRUE) {

  ## ------------------------------------------------------------------
  ## ------------------------------------------------------------------

  old_handlers <- progressr::handlers()
  if (!progress) {
    progressr::handlers("void")
  }
  on.exit(progressr::handlers(old_handlers), add = TRUE)
  ## ------------------------------------------------------------------
  ## 1. Check  input
  ## ------------------------------------------------------------------
  if (length(fit$coordt) == 1) fit$coordt <- NULL
  if (is.null(fit$sensmat)) stop("Sensitivity matrix is missing: use sensitivity = TRUE in GeoFit")
  if (!(method %in% c("cholesky", "TB", "CE"))) stop("The method of simulation is not correct")
  if (!is.numeric(K) || K < 1) stop("K must be a positive integer")
  if (is.numeric(alpha) && !(alpha > 0 && alpha < 1)) stop("alpha must be numeric between 0 and 1")
  if (!is.logical(progress)) stop("progress must be logical (TRUE or FALSE)")

  if (is.null(optimizer)) {
    optimizer <- fit$optimizer
    lower     <- fit$lower
    upper     <- fit$upper
  }

  model <- fit$model
  cat("Parametric bootstrap can be time consuming ...\n")

  ## ------------------------------------------------------------------
  ## 2. Misspecification mapping
  ## ------------------------------------------------------------------
  if (fit$missp) {
    model_map <- c(
      StudentT     = "Gaussian_misp_StudentT",
      Poisson      = "Gaussian_misp_Poisson",
      PoissonZIP   = "Gaussian_misp_PoissonZIP",
      SkewStudentT = "Gaussian_misp_SkewStudentT",
      Tukeygh      = "Gaussian_misp_Tukeygh"
    )
    if (model %in% names(model_map)) model <- model_map[[model]]
  }

  ## ------------------------------------------------------------------
  ## 3. X matrix computing
  ## ------------------------------------------------------------------
  dimat <- fit$numtime * fit$numcoord
  X_use <- if (!is.null(fit$X) && ncol(fit$X) == 1 && all(fit$X[1:dimat] == 1)) NULL else fit$X

  coords <- cbind(fit$coordx, fit$coordy)
  if (fit$bivariate && is.null(fit$coordx_dyn)) {
    coords <- coords[seq_len(nrow(coords) / 2), , drop = FALSE]
  }

  ## ------------------------------------------------------------------
  ## 4. Simulation
  ## ------------------------------------------------------------------
  sim_args <- list(
    coordx     = coords,
    coordt     = fit$coordt,
    coordx_dyn = fit$coordx_dyn,
    anisopars  = fit$anisopars,
    corrmodel  = fit$corrmodel,
    model      = fit$model,
    param      = append(fit$param, fit$fixed),
    grid       = fit$grid,
    X          = fit$X,
    n          = fit$n,
    distance   = fit$distance,
    radius     = fit$radius,
    nrep       = K,
    progress=progress
  )

  if (is.null(fit$copula)) {
    if (method == "cholesky") {
      #cat("Performing", K, "simulations using Cholesky method...\n")
      data_sim <- do.call(GeoSim, c(sim_args, list(sparse = sparse, method = method)))
    } else if (method %in% c("TB", "CE")) {
     # cat("Performing", K, "simulation using", method ," method"," ...\n")
      data_sim <- do.call(GeoSimapprox,
                          c(sim_args, list(method = method, L = L, parallel = parallel)))
    } else {
      stop("Unsupported method for simulation")
    }
  } else {
    if (method == "cholesky") {
      data_sim <- do.call(GeoSimCopula,
                          c(sim_args, list(copula = fit$copula, sparse = sparse, method = method)))
    } else {
      stop("Unsupported method for copula simulation")
    }
  }
  
  ## ------------------------------------------------------------------
  ## 5. setting parallelization
  ## ------------------------------------------------------------------
  coremax <- parallel::detectCores()
  if (is.na(coremax) || coremax <= 1) {
    parallel <- FALSE
    ncores   <- 1
  } else {
    if (is.null(parallel)) parallel <- TRUE
    ncores <- max(1, min(if (is.null(ncores)) coremax - 1 else ncores, K))
  }

  ## ------------------------------------------------------------------
  ## 6. estimation
  ## ------------------------------------------------------------------
estimate_fun <- function(k) {
  tryCatch(
    {
      # Cattura l'output ma salva il risultato della funzione
      captured_output <- capture.output({
        result <- GeoFit(
          data       = data_sim$data[[k]],
          start      = fit$param,
          fixed      = fit$fixed, 
          coordx     = coords,
          coordt     = fit$coordt, 
          coordx_dyn = fit$coordx_dyn, 
          copula     = fit$copula, 
          anisopars  = fit$anisopars,
          est.aniso  = fit$est.aniso, 
          lower      = lower, 
          upper      = upper, 
          neighb     = fit$neighb,
          corrmodel  = fit$corrmodel, 
          model      = model, 
          sparse     = FALSE, 
          n          = fit$n,
          maxdist    = fit$maxdist, 
          maxtime    = fit$maxtime, 
          optimizer  = optimizer, 
          grid       = fit$grid,
          likelihood = fit$likelihood, 
          type       = fit$type, 
          X          = X_use, 
          distance   = fit$distance,
          radius     = fit$radius
        )
      }, file = nullfile())
      
      # Restituisce il risultato della funzione GeoFit
      result
    },
    error = function(e) NULL
  )
}
#estimate_fun <- function(k) {
#  tryCatch(
#    withr::with_output_sink(
#      tempfile(),  # file temporaneo che sarà cancellato
#      GeoFit(
#        data       = data_sim$data[[k]],start      = fit$param,fixed      = fit$fixed, coordx     = coords,
#        coordt     = fit$coordt, coordx_dyn = fit$coordx_dyn, copula     = fit$copula, anisopars  = fit$anisopars,
#        est.aniso  = fit$est.aniso, lower      = lower, upper      = upper, neighb     = fit$neighb,
#        corrmodel  = fit$corrmodel, model      = model, sparse     = FALSE, n          = fit$n,
#        maxdist    = fit$maxdist, maxtime    = fit$maxtime, optimizer  = optimizer, grid       = fit$grid,
#        likelihood = fit$likelihood, type       = fit$type, X          = X_use, distance   = fit$distance,radius     = fit$radius
#      )
#    ),
#    error = function(e) NULL
#  )
#}
  ## ------------------------------------------------------------------
  ## 7.  bootstrap
  ## ------------------------------------------------------------------
  if (!parallel) {
    cat("Performing", K, "estimations sequentially...\n")
    if (progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
    }
    pb <- progressr::progressor(along = seq_len(K))

    res_list <- vector("list", K)
    for (k in seq_len(K)) {
      if (progress) {
        pb(sprintf("k=%d", k))
      }
      res_est <- suppressWarnings(estimate_fun(k))
      if (!is.null(res_est) &&
          res_est$convergence == "Successful" &&
          is.finite(res_est$logCompLik) &&
          res_est$logCompLik < 1.0e8) {
        res_list[[k]] <- unlist(res_est$param)
      } else {
        res_list[[k]] <- NULL
      }
    }
    res <- do.call(rbind, Filter(Negate(is.null), res_list))
  } else {
    cat("Performing", K, "estimations using", ncores, "cores...\n")

    future::plan(future::multisession, workers = ncores)
    on.exit(future::plan(future::sequential), add = TRUE)

    if (progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = seq_len(K))
    }

    num_params <- length(fit$param)
    xx <- foreach::foreach(
      k = seq_len(K),
      .combine = rbind,
      .options.future = list(
        packages = "GeoModels",
        seed     = TRUE,
        stdout   = FALSE,
        conditions = "none"
      )
    ) %dofuture% {
      if (progress) {
        pb(sprintf("k=%d", k))
      }
      res_est <- suppressWarnings(estimate_fun(k))

      if (is.null(res_est) ||
          res_est$convergence != "Successful" ||
          !is.finite(res_est$logCompLik) ||
          res_est$logCompLik >= 1.0e8) {
        c(rep(NA_real_, num_params), convergence = FALSE, logCompLik = NA_real_)
      } else {
        c(unlist(res_est$param), convergence = TRUE, logCompLik = res_est$logCompLik)
      }
    }

    valid <- xx[, "convergence"] == TRUE & !is.na(xx[, "logCompLik"]) & xx[, "logCompLik"] < 1.0e8
    res   <- xx[valid, seq_len(num_params), drop = FALSE]
  }

  ## ------------------------------------------------------------------
  ## 8. Post-processing
  ## ------------------------------------------------------------------
  n_successful <- nrow(res)
  if (n_successful < 2) stop("Insufficient successful bootstrap iterations for variance estimation")
  cat("Successful bootstrap iterations:", n_successful, "out of", K, "\n")

  invG   <- var(res)
  G      <- try(solve(invG), silent = TRUE)
  if (!is.matrix(G)) warning("Bootstrap estimated Godambe matrix is singular")

  stderr <- sqrt(pmax(0, diag(invG)))
  names(stderr) <- names(fit$param)

  numparam <- length(fit$param)
  if ((fit$likelihood == "Marginal" && fit$type %in% c("Independence", "Pairwise")) ||
      (fit$likelihood == "Conditional" && fit$type == "Pairwise")) {
    H       <- fit$sensmat
    penalty <- sum(diag(H %*% invG))
    claic   <- -2 * fit$logCompLik + 2 * penalty
    clbic   <- -2 * fit$logCompLik + log(dimat) * penalty
    fit$varimat <- H %*% invG %*% H
  } else if (fit$likelihood == "Full" && fit$type == "Standard") {
    claic <- -2 * fit$logCompLik + 2 * numparam
    clbic <- -2 * fit$logCompLik + log(dimat) * 2 * numparam
  } else {
    claic <- clbic <- NA_real_
  }

  ## ------------------------------------------------------------------
  ## 9. final output
  ## ------------------------------------------------------------------
  fit$claic      <- claic
  fit$clbic      <- clbic
  fit$stderr     <- stderr
  fit$varcov     <- invG
  fit$estimates  <- res

  z_alpha <- qnorm(1 - (1 - alpha) / 2)
  aa      <- z_alpha * stderr
  pp      <- as.numeric(fit$param)
  fit$conf.int <- rbind(pp - aa, pp + aa)
  colnames(fit$conf.int) <- names(fit$param)
  rownames(fit$conf.int) <- c("Lower", "Upper")

  fit$pvalues <- 2 * pnorm(-abs(pp / stderr))
  names(fit$pvalues) <- names(fit$param)

  return(fit)
}