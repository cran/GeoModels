GeoVarestbootstrap <- function(fit, K = 100, sparse = FALSE, GPU = NULL, local = c(1, 1), optimizer = NULL,
                              lower = NULL, upper = NULL, method = "cholesky", alpha = 0.95, L = 1000,
                              parallel = TRUE, ncores = NULL) {

  # Controlli di input
  if (length(fit$coordt) == 1) fit$coordt <- NULL
  if (is.null(fit$sensmat)) stop("Sensitivity matrix is missing: use sensitivity=TRUE in GeoFit")
  if (!(method %in% c("cholesky", "TB", "CE"))) stop("The method of simulation is not correct")
  if (is.null(optimizer)) {
    optimizer <- fit$optimizer
    lower <- fit$lower
    upper <- fit$upper
  }
  if (is.numeric(alpha) && !(alpha > 0 && alpha < 1)) stop("alpha must be numeric between 0 and 1")
  if (!is.numeric(K) || K < 1) stop("K must be a positive integer")

  model <- fit$model

  cat("Parametric bootstrap can be time consuming ...\n")

  # Gestione modelli con misspecification
  if (fit$missp) {
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
  
  # Ottimizzazione gestione matrice X
  if (!is.null(fit$X) && sum(fit$X[1:dimat] == 1) == dimat && ncol(fit$X) == 1) {
    X_use <- NULL
  } else {
    X_use <- fit$X
  }

  # Preparazione coordinate
  coords <- cbind(fit$coordx, fit$coordy)
  if (fit$bivariate && is.null(fit$coordx_dyn)) {
    coords <- coords[1:(length(fit$coordx) / 2), , drop = FALSE]
  }
  N <- nrow(coords)

  # Simulazione dati
  if (is.null(fit$copula)) {  # modelli non-copula
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
  } else {  # modeli copula
    if (method == "cholesky") {
      data_sim <- GeoSimCopula(
        coordx = coords, coordt = fit$coordt, coordx_dyn = fit$coordx_dyn, anisopars = fit$anisopars,
        corrmodel = fit$corrmodel, model = fit$model, copula = fit$copula,
        param = append(fit$param, fit$fixed), sparse = sparse,
        grid = fit$grid, X = fit$X, n = fit$n, method = method,
        distance = fit$distance, radius = fit$radius, nrep = K
      )
    } else {
      stop("Unsupported method for copula simulation")
    }
  }

  # Gestione parallelizzazione
  coremax <- parallel::detectCores()
  if (is.na(coremax) || coremax == 1) {
    parallel <- FALSE
  } else {
    if (is.null(ncores)) ncores <- max(1, min(coremax - 1, K))
    ncores <- max(1, min(ncores, K))  # Assicura che ncores non superi K
  }

  # Funzione di stima bootstrap
  estimate_fun <- function(k) {
    tryCatch({
      suppressWarnings({
        result <- GeoFit(
          data = data_sim$data[[k]], start = fit$param, fixed = fit$fixed,
          coordx = coords, coordt = fit$coordt, coordx_dyn = fit$coordx_dyn,
          copula = fit$copula, anisopars = fit$anisopars, est.aniso = fit$est.aniso,
          lower = lower, upper = upper, neighb = fit$neighb,
          corrmodel = fit$corrmodel, model = model, sparse = FALSE, n = fit$n,
          maxdist = fit$maxdist, maxtime = fit$maxtime,
          optimizer = optimizer, grid = fit$grid, likelihood = fit$likelihood,
          type = fit$type, X = X_use, distance = fit$distance, radius = fit$radius
        )
        return(result)
      })
    }, error = function(e) {
      return(list(convergence = "Failed", logCompLik = Inf, param = NULL))
    })
  }

  # Esecuzione bootstrap
  if (!parallel) {
    cat("Performing", K, "estimations sequentially...\n")
    res_list <- vector("list", K)
    
    for (k in seq_len(K)) {
      if (k %% max(1, K %/% 10) == 0) {
        cat("Progress:", round(100 * k / K, 1), "%\n")
      }
      
      res_est <- estimate_fun(k)
      if (!is.null(res_est$convergence) && res_est$convergence == "Successful" && 
          !is.null(res_est$logCompLik) && res_est$logCompLik < 1.0e8) {
        res_list[[k]] <- unlist(res_est$param)
      } else {
        res_list[[k]] <- NULL
      }
    }
    
    valid_results <- Filter(Negate(is.null), res_list)
    if (length(valid_results) == 0) {
      stop("No successful bootstrap iterations")
    }
    res <- do.call(rbind, valid_results)
    
  } else {
    cat("Performing", K, "estimations using", ncores, "cores...\n")
    
    # Controllo disponibilità pacchetti
    if (!requireNamespace("future", quietly = TRUE) || 
        !requireNamespace("foreach", quietly = TRUE) ||
        !requireNamespace("doFuture", quietly = TRUE) ||
        !requireNamespace("progressr", quietly = TRUE)) {
      warning("Required packages for parallel execution not available. Running sequentially.")
      parallel <- FALSE
    }
    
    if (parallel) {
      future::plan(future::multisession, workers = ncores)
      on.exit(future::plan(future::sequential), add = TRUE)
      
      tryCatch({
        progressr::handlers(global = TRUE)
        progressr::handlers("txtprogressbar")
        pb <- progressr::progressor(along = 1:K)
        
        xx <- foreach::foreach(k = seq_len(K), .combine = rbind,
                              .options.future = list(seed = TRUE, stdout = NA,
                                                   conditions = character(0))) %dofuture% {
          pb(sprintf("k=%d", k))
          res_est <- estimate_fun(k)
          c(unlist(res_est$param), convergence = res_est$convergence, 
            logCompLik = res_est$logCompLik)
        }
        
        # Filtraggio risultati validi
        conv_idx <- !is.na(xx[, "convergence"]) & 
                   xx[, "convergence"] == "Successful" & 
                   !is.na(xx[, "logCompLik"]) & 
                   xx[, "logCompLik"] < 1.0e8
        
        if (sum(conv_idx) == 0) {
          stop("No successful bootstrap iterations")
        }
        
        res <- xx[conv_idx, !colnames(xx) %in% c("convergence", "logCompLik"), drop = FALSE]
        
      }, error = function(e) {
        warning("Parallel execution failed, falling back to sequential: ", e$message)
        parallel <<- FALSE
      })
    }
    
    # Fallback sequenziale se parallel fallisce
    if (!parallel) {
      # Ripeti il codice sequenziale
      res_list <- vector("list", K)
      for (k in seq_len(K)) {
        if (k %% max(1, K %/% 10) == 0) {
          cat("Progress:", round(100 * k / K, 1), "%\n")
        }
        res_est <- estimate_fun(k)
        if (!is.null(res_est$convergence) && res_est$convergence == "Successful" && 
            !is.null(res_est$logCompLik) && res_est$logCompLik < 1.0e8) {
          res_list[[k]] <- unlist(res_est$param)
        } else {
          res_list[[k]] <- NULL
        }
      }
      valid_results <- Filter(Negate(is.null), res_list)
      if (length(valid_results) == 0) {
        stop("No successful bootstrap iterations")
      }
      res <- do.call(rbind, valid_results)
    }
  }

  # Controllo risultati
  if (nrow(res) < 2) {
    stop("Insufficient successful bootstrap iterations for variance estimation")
  }
  
  cat("Successful bootstrap iterations:", nrow(res), "out of", K, "\n")

  # Calcolo varianza e matrice di Godambe
  numparam <- length(fit$param)
  invG <- var(res)
  
  # Controllo singolarità
  G <- tryCatch({
    solve(invG)
  }, error = function(e) {
    warning("Bootstrap estimated Godambe matrix is singular")
    return(NULL)
  })
  
  if (is.null(G) || !is.matrix(G)) {
    warning("Bootstrap estimated Godambe matrix is singular - cannot compute inverse")
    G <- NULL
  }

  stderr <- sqrt(pmax(0, diag(invG)))  # Evita radici negative
  names(stderr) <- names(fit$param)  # Assegna i nomi dei parametri

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
    clbic <- -2 * fit$logCompLik + log(dimat) * numparam
  } else {
    claic <- NA_real_
    clbic <- NA_real_
  }

  # Aggiornamento risultati nel fit object
  fit$claic <- claic
  fit$clbic <- clbic
  fit$stderr <- stderr
  fit$varcov <- invG
  fit$estimates <- apply(res, 2, function(x) as.numeric(gsub('"', '', x)))

  # Intervalli di confidenza e p-values
  if (all(stderr > 0)) {
    aa <- qnorm(1 - (1 - alpha) / 2) * stderr
    pp <- as.numeric(fit$param)
    names(pp) <- names(fit$param)  # Mantieni i nomi
    fit$conf.int <- rbind(pp - aa, pp + aa)
    colnames(fit$conf.int) <- names(fit$param)  # Nomi per gli intervalli
    rownames(fit$conf.int) <- c("Lower", "Upper")
    fit$pvalues <- 2 * pnorm(-abs(pp / stderr))
    names(fit$pvalues) <- names(fit$param)  # Nomi per i p-values
  } else {
    warning("Some standard errors are zero - confidence intervals may be unreliable")
    fit$conf.int <- matrix(NA, nrow = 2, ncol = length(fit$param))
    colnames(fit$conf.int) <- names(fit$param)
    rownames(fit$conf.int) <- c("Lower", "Upper")
    fit$pvalues <- rep(NA, length(fit$param))
    names(fit$pvalues) <- names(fit$param)
  }

  return(fit)
}