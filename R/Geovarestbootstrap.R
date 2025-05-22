GeoVarestbootstrap <- function(fit, K = 100, sparse = FALSE, GPU = NULL, local = c(1, 1), optimizer = NULL,
                              lower = NULL, upper = NULL, method = "cholesky", alpha = 0.95, L = 1000,
                              parallel = NULL, ncores = NULL) {

  if (length(fit$coordt) == 1) fit$coordt <- NULL
  if (is.null(fit$sensmat)) stop("Sensitivity matrix is missing: use sensitivity=TRUE in GeoFit")
  if (!(method %in% c("cholesky", "TB", "CE"))) stop("The method of simulation is not correct")
  if (is.null(optimizer)) {
    optimizer <- fit$optimizer
    lower <- fit$lower
    upper <- fit$upper
  }
  if (is.numeric(alpha) && !(alpha > 0 && alpha < 1)) stop("alpha must be numeric between 0 and 1")

  model <- fit$model

  cat("Parametric bootstrap can be time consuming ...\n")

  if (fit$missp) {  # misspecification corrections
    model_map <- c(
      StudentT = "Gaussian_misp_StudentT",
      Poisson = "Gaussian_misp_Poisson",
      PoissonZIP = "Gaussian_misp_PoissonZIP",
      SkewStudentT = "Gaussian_misp_SkewStudentT",
      Tukeygh = "Gaussian_misp_Tukeygh"
    )
    if (model %in% names(model_map)) model <- model_map[[model]]
  }

  dimat <- fit$numtime * fit$numcoord
  # Ottimizzazione punto 1: evita modifiche ripetute a fit$X
  if (!is.null(fit$X) && sum(fit$X[1:dimat] == 1) == dimat && ncol(fit$X) == 1) {
    X_use <- NULL
  } else {
    X_use <- fit$X
  }

  coords <- cbind(fit$coordx, fit$coordy)
  if (fit$bivariate && is.null(fit$coordx_dyn)) coords <- coords[1:(length(fit$coordx) / 2), ]
  N <- nrow(coords)

  #cat("Performing", K, "simulations....\n")

  # Simulazione dati
  if (is.null(fit$copula)) {  # non copula models
    if (method == "cholesky") {
      data_sim <- GeoSim(
        coordx = coords, coordt = fit$coordt, coordx_dyn = fit$coordx_dyn, anisopars = fit$anisopars,
        corrmodel = fit$corrmodel, model = fit$model, param = append(fit$param, fit$fixed),
        sparse = sparse, grid = fit$grid, X = fit$X, n = fit$n, method = method,
        distance = fit$distance, radius = fit$radius, nrep = K
      )
    } else if (method %in% c("TB", "CE")) {
      data_sim <- GeoSimapprox(
        coordx = coords, coordt = fit$coordt, coordx_dyn = fit$coordx_dyn, anisopars = fit$anisopars,
        corrmodel = fit$corrmodel, model = fit$model, param = append(fit$param, fit$fixed),
        grid = fit$grid, X = fit$X, n = fit$n, method = method,
        parallel = parallel, ncores = ncores, L = L,
        distance = fit$distance, radius = fit$radius, nrep = K
      )
    } else {
      stop("Unsupported method for simulation")
    }
  } else {  # copula models
    if (method == "cholesky") {
      data_sim <- GeoSimCopula(
        coordx = coords, coordt = fit$coordt, coordx_dyn = fit$coordx_dyn, anisopars = fit$anisopars,
        corrmodel = fit$corrmodel, model = fit$model, copula = fit$copula,
        param = append(fit$param, fit$fixed),  sparse = sparse,
        grid = fit$grid, X = fit$X, n = fit$n, method = method,
        distance = fit$distance, radius = fit$radius, nrep = K
      )
    } else {
      stop("Unsupported method for copula simulation")
    }
  }

  # Ottimizzazione punto 2: gestione automatica parallelizzazione e numero core
  coremax <- parallel::detectCores()
  if (is.na(coremax) || coremax == 1) {
    parallel <- FALSE
  } else {
    if (is.null(parallel)) parallel <- TRUE
    if (is.null(ncores)) ncores <- max(1, coremax - 1)
  }

  # Funzione di stima bootstrap
  estimate_fun <- function(k) {
    GeoFit(
      data = data_sim$data[[k]], start = fit$param, fixed = fit$fixed,
      coordx = coords, coordt = fit$coordt, coordx_dyn = fit$coordx_dyn,
      copula = fit$copula, anisopars = fit$anisopars, est.aniso = fit$est.aniso,
      lower = lower, upper = upper, neighb = fit$neighb,
      corrmodel = fit$corrmodel, model = model, sparse = FALSE, n = fit$n,
      maxdist = fit$maxdist, maxtime = fit$maxtime,
      optimizer = optimizer, grid = fit$grid, likelihood = fit$likelihood,
      type = fit$type, X = X_use, distance = fit$distance, radius = fit$radius
    )
  }

  if (!parallel) {
    cat("Performing", K, "estimations sequentially...\n")
    # Disabilito progress bar per velocità (opzionale)
    res_list <- vector("list", K)
    for (k in seq_len(K)) {
      res_est <- estimate_fun(k)
      if (res_est$convergence == "Successful" && res_est$logCompLik < 1.0e8) {
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

    progressr::handlers(global = TRUE)
    progressr::handlers("txtprogressbar")
    pb <- progressr::progressor(along = 1:K)

    xx <- foreach::foreach(k = seq_len(K), .combine = rbind,
                          .options.future = list(seed = TRUE)) %dofuture% {
      pb(sprintf("k=%d", k))
      res_est <- estimate_fun(k)
      c(unlist(res_est$param), convergence = res_est$convergence, logCompLik = res_est$logCompLik)
    }

    # Filtra risultati convergenti e con logCompLik valido
    conv_idx <- xx[, "convergence"] == "Successful" & xx[, "logCompLik"] < 1.0e8
    res <- xx[conv_idx, colnames(xx) != "convergence" & colnames(xx) != "logCompLik", drop = FALSE]

    
  }

  # Calcolo varianza e Godambe
  numparam <- length(fit$param)
  invG <- var(res)
  G <- try(solve(invG), silent = TRUE)
  if (!is.matrix(G)) warning("Bootstrap estimated Godambe matrix is singular")

  stderr <- sqrt(diag(invG))

  # Calcolo criteri di informazione
  if ((fit$likelihood == "Marginal" && fit$type %in% c("Independence", "Pairwise")) ||
      (fit$likelihood == "Conditional" && fit$type == "Pairwise")) {

    H <- fit$sensmat
    penalty <- sum(diag(H %*% invG))
    claic <- -2 * fit$logCompLik + 2 * penalty
    clbic <- -2 * fit$logCompLik + log(dimat) * penalty
    fit$varimat <- H %*% invG %*% H
  } else if (fit$likelihood == "Full" && fit$type == "Standard") {
    claic <- -2 * fit$logCompLik + 2 * numparam
    clbic <- -2 * fit$logCompLik + log(dimat) * 2 * numparam
  } else {
    claic <- NA_real_
    clbic <- NA_real_
  }

  # Aggiornamento risultati
  fit$claic <- claic
  fit$clbic <- clbic
  fit$stderr <- stderr
  fit$varcov <- invG
  fit$estimates <- res

  # Intervalli di confidenza e p-values
  aa <- qnorm(1 - (1 - alpha) / 2) * stderr
  pp <- as.numeric(fit$param)
  fit$conf.int <- rbind(pp - aa, pp + aa)
  fit$pvalues <- 2 * pnorm(-abs(pp / stderr))

  return(fit)
}
