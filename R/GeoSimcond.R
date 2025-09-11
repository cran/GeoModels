GeoSimcond <- function(estobj = NULL, data, coordx, coordy = NULL, coordz = NULL, coordt = NULL,
                       coordx_dyn = NULL, corrmodel, distance = "Eucl",
                       grid = FALSE, loc, maxdist = NULL, maxtime = NULL, method = "Cholesky",
                       model = "Gaussian", n = 1, nrep = 1, local = FALSE, L = 1000,
                       neighb = NULL, param, anisopars = NULL, radius = 1, sparse = FALSE,
                       time = NULL, copula = NULL, X = NULL, Xloc = NULL, Mloc = NULL,
                       parallel = FALSE, ncores = NULL, progress=FALSE) {



################################################################################# 
### Helper function for conditional skew gaussian simulation 
################################################################################# 
SkewGaussianSimcond <- function(coord_obs, loc, data, param, corrmodel,
                                nrep = 1, method = "Cholesky", local=FALSE, neighb=NULL, L=NULL, 
                                parallel = FALSE, ncores = 1, progress = FALSE) {
  mm          <- param$mean
  sigma_sqrt  <- sqrt(param$sill)
  eta         <- param$skew
  keep        <- CorrParam(corrmodel)
  corr_param  <- param[keep]
  
  cat("Computing Kriging weights ...\n")
  
  if(!local)
    GeoW        <- GeoKrigWeights(
      coordx   = coord_obs, corrmodel = corrmodel,
      loc      = loc, model = "Gaussian",
      param    = c(sill = 1, corr_param))
  else stop("Local option cannot be used for gibbs sampling ...\n")

  # Pre-calcoli per ottimizzazione
  Sigma.mat.inv <- MatInv(mtx = GeoW$covmatrix)
  M0 <- nrow(coord_obs)
  N0 <- nrow(loc)
  Wglob <- GeoW$weights
  mm_loc <- mm
  s_sqrt <- sigma_sqrt
  
  # Pre-calcoli costanti per evitare ricalcoli
  sign_eta <- sign(eta)
  abs_eta <- abs(eta)
  eta_over_s_sqrt <- eta / s_sqrt
  sqrt_param_sill <- sqrt(param$sill)
  Z_trans_factor <- sign_eta * sigma_sqrt  # Pre-calcolo per Z.trans
  
  # Gaussian simulations
  if (method == "TB" || method == "CE") {
    Gauss.all <- GeoSimapprox(
      coordx   = rbind(coord_obs, loc),
      corrmodel = corrmodel, method = method, model = "Gaussian",
      param    = c(corr_param, mean = 0, nugget = 0), L=L,
      nrep     = 2 * nrep, progress = FALSE,
      parallel = FALSE
    )
  } else {
    Gauss.all <- GeoSim(
      coordx   = rbind(coord_obs, loc),
      corrmodel = corrmodel, model = "Gaussian",
      param    = c(corr_param, mean = 0, nugget = 0),
      nrep     = 2 * nrep, progress = FALSE
    )
  }
  
  pair_idx <- split(seq_len(2 * nrep), ceiling(seq_len(2 * nrep) / 2))

  ####### OTTIMIZZAZIONE: INI VAL SAMPLER FUNCTION VETTORIZZATA #######
  Ini.val.fun.vectorized <- function(x_vec, A) {
    if(length(A) != 2) stop("A deve essere un vettore di lunghezza 2")
    
    n <- length(x_vec)
    result_matrix <- matrix(0, nrow = 2, ncol = n)
    
    # Chiamata vettorizzata più efficiente
    for(i in seq_len(n)) {
      vvv <- dotCall64::.C64("rnorm_constraint_simple",
                             SIGNATURE = c("double", "double", "double", "double", "double"), 
                             A = as.double(A),
                             b = as.double(x_vec[i]),
                             mu = as.double(c(0,0)),
                             sigma = as.double(1),
                             result = double(2),
                             INTENT = c("r", "r", "r", "r", "w"))$result
      result_matrix[, i] <- vvv
    }
    return(result_matrix)
  }
  ####### END INI VAL SAMPLER FUNCTION #######
  
  # Funzione helper per Gibbs sampling - evita duplicazione
  run_skew_gibbs_sampler <- function(Z_trans, M0, Sigma.mat.inv, sigma_sqrt, abs_eta, sqrt_param_sill, eta) {
    # Genera valori iniziali in modo ottimizzato
    A_vec <- c(sqrt_param_sill, eta)
    ini.vals <- Ini.val.fun.vectorized(Z_trans, A_vec)
    
    if(!is.matrix(ini.vals) || nrow(ini.vals) != 2 || ncol(ini.vals) != M0) {
      stop("Errore nella generazione dei valori iniziali")
    }
    
    # Gibbs sampler
    gibbs.result <- dotCall64::.C64(
      "skew_gaussian_gibbs_sampler",
      SIGNATURE = c("double", "integer", "double", "double", "integer", "double", "double"), 
      data_obs     = as.double(Z_trans),
      n            = as.integer(M0),
      Sigma_mat_inv= as.double(Sigma.mat.inv),
      eta          = as.double(c(sigma_sqrt, abs_eta)),
      n_iter       = as.integer(1000),
      data_x       = as.double(ini.vals[1,]),
      data_y       = as.double(ini.vals[2,]),
      INTENT = c("r", "r", "r", "r", "r", "rw", "rw"), 
      VERBOSE = 0, 
      NAOK = TRUE, 
      PACKAGE = "GeoModels"
    )
    
    if(is.null(gibbs.result$data_x) || is.null(gibbs.result$data_y)) {
      return(NULL)
    }
    
    return(list(data_x = gibbs.result$data_x, data_y = gibbs.result$data_y))
  }
  
  # Funzione helper per processamento simulazione ottimizzata
  process_skew_simulation <- function(gibbs_result, idx_xy, Gauss_data, M0, N0, Wglob, 
                                      sign_eta, eta_over_s_sqrt, mm_loc) {
    simX <- Gauss_data[[idx_xy[1]]]
    simY <- Gauss_data[[idx_xy[2]]]
    
    # Processamento X
    sim.train.x <- simX[seq_len(M0)]
    sim.valid.x <- simX[(M0 + 1):(N0 + M0)]
    sk.pred.x <- as.numeric(crossprod(Wglob, gibbs_result$data_x - sim.train.x))
    geosim.x <- sim.valid.x + sk.pred.x
    
    # Processamento Y
    sim.train.y <- simY[seq_len(M0)]
    sim.valid.y <- simY[(M0 + 1):(N0 + M0)]
    sk.pred.y <- as.numeric(crossprod(Wglob, gibbs_result$data_y - sim.train.y))
    geosim.y <- sim.valid.y + sk.pred.y
    
    # Calcolo finale ottimizzato - operazioni vettorizzate
    res_i <- sign_eta * geosim.x + eta_over_s_sqrt * abs(geosim.y)
    return(mm_loc + res_i)
  }
  
  if (!parallel) {
    cat("Performing", nrep, "conditional simulations ...\n")
    
    if (progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = seq_len(nrep))
    }
    
    res <- vector("list", nrep)
    
    # Pre-calcolo Z.trans una volta sola per tutte le simulazioni
    Z_trans <- Z_trans_factor * data
    
    for (i in seq_len(nrep)) {
      # Gibbs sampling usando funzione helper
      gibbs_result <- run_skew_gibbs_sampler(Z_trans, M0, Sigma.mat.inv, 
                                             sigma_sqrt, abs_eta, sqrt_param_sill, eta)
      
      if(is.null(gibbs_result)) {
        warning(paste("Gibbs sampler fallito alla simulazione", i))
        next
      }
      
      # Processamento usando funzione helper ottimizzata
      idx_xy <- pair_idx[[i]]
      res[[i]] <- process_skew_simulation(gibbs_result, idx_xy, Gauss.all$data, 
                                          M0, N0, Wglob, sign_eta, eta_over_s_sqrt, mm_loc)
      
      if (progress) {
        pb(sprintf("Simulation %d/%d completed", i, nrep))
      }
    }
    
  } else {
    cat("Performing", nrep, "conditional simulations using", ncores, "cores ...\n")
    
    if (progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = seq_len(nrep))
    }
    
    # Ottimizzazione memoria: limite più conservativo
    old_limit <- options(future.globals.maxSize = 1500 * 1024^2)  # 1.5 GB
    on.exit(options(old_limit), add = TRUE)
    
    # Variabili essenziali ottimizzate - solo pre-calcoli
    essential_vars <- list(
      eta = eta,
      abs_eta = abs_eta,
      sign_eta = sign_eta,
      eta_over_s_sqrt = eta_over_s_sqrt,
      sigma_sqrt = sigma_sqrt,
      sqrt_param_sill = sqrt_param_sill,
      Z_trans = Z_trans_factor * data,  # Pre-calcolato
      M0 = M0,
      N0 = N0,
      Sigma.mat.inv = Sigma.mat.inv,
      Wglob = Wglob,
      mm_loc = mm_loc,
      pair_idx = pair_idx,
      Gauss_data = Gauss.all$data
    )
    
    # Worker function ottimizzata
    run_simulation_worker <- function(i, vars) {
      # Estrai variabili pre-calcolate
      eta <- vars$eta
      abs_eta <- vars$abs_eta
      sign_eta <- vars$sign_eta
      eta_over_s_sqrt <- vars$eta_over_s_sqrt
      sigma_sqrt <- vars$sigma_sqrt
      sqrt_param_sill <- vars$sqrt_param_sill
      Z_trans <- vars$Z_trans
      M0 <- vars$M0
      N0 <- vars$N0
      Sigma.mat.inv <- vars$Sigma.mat.inv
      Wglob <- vars$Wglob
      mm_loc <- vars$mm_loc
      pair_idx <- vars$pair_idx
      Gauss_data <- vars$Gauss_data
      
      # Gibbs sampling
      gibbs_result <- run_skew_gibbs_sampler(Z_trans, M0, Sigma.mat.inv, 
                                             sigma_sqrt, abs_eta, sqrt_param_sill, eta)
      
      if(is.null(gibbs_result)) {
        warning(paste("Gibbs sampler fallito alla simulazione", i))
        return(NULL)
      }
      
      # Processamento ottimizzato
      idx_xy <- pair_idx[[i]]
      return(process_skew_simulation(gibbs_result, idx_xy, Gauss_data, 
                                     M0, N0, Wglob, sign_eta, eta_over_s_sqrt, mm_loc))
    }
    
    # Parallelizzazione
    future::plan(future::multisession, workers = ncores)
    on.exit(future::plan(future::sequential), add = TRUE)
    
    res <- future.apply::future_lapply(seq_len(nrep), function(i) {
      result <- run_simulation_worker(i, essential_vars)
      if (progress) {
        pb(sprintf("Parallel simulation %d/%d completed", i, nrep))
      }
      return(result)
    }, future.seed = TRUE, future.globals = FALSE)
  
    # Rimuove risultati nulli e avvisa
    res <- res[!sapply(res, is.null)]
    failed_sims <- nrep - length(res)
    if (failed_sims > 0) {
      warning(paste(failed_sims, "simulazioni su", nrep, "sono fallite"))
    }
  }
  
  return(res)
} 
################################################################################# 
### Helper function for conditional Gamma simulation 
################################################################################# 
GammaSimcond <- function(coord_obs, loc, data, param, corrmodel, 
                         nrep = 1, method = "Cholesky", local = FALSE, neighb = NULL, 
                         L = NULL, parallel = FALSE, ncores = 1, progress = FALSE) {
  mm <- param$mean
  shape <- round(param$shape)
  keep <- CorrParam(corrmodel)
  corr_param <- param[keep]
  
  cat("Computing Kriging weights ...\n")
  if (!local) 
    GeoW <- GeoKrigWeights(coordx = coord_obs, corrmodel = corrmodel, 
                           loc = loc, model = "Gaussian", param = c(sill = 1, 
                                                                    corr_param))
  else stop("Local option cannot be used for gibbs sampling ...\n")
  
  # Pre-calcoli per ottimizzazione
  Sigma.mat.inv <- MatInv(mtx = GeoW$covmatrix)
  M0 <- nrow(coord_obs)
  N0 <- nrow(loc)
  Wglob <- GeoW$weights
  mm_loc <- mm
  data_transformed <- exp(-mm_loc) * data  # Pre-calcolo
  sqrt_factor <- sqrt(2 / shape)  # Pre-calcolo costante
  exp_mm_half <- exp(mm_loc) * 0.5  # Pre-calcolo per output finale

  # Generazione simulazioni Gaussiane
  if (method == "TB" || method == "CE") {
    Gauss.all <- GeoSimapprox(coordx = rbind(coord_obs, loc), 
                              corrmodel = corrmodel, method = method, 
                              model = "Gaussian", param = c(corr_param, mean = 0, 
                                                            nugget = 0, sill = 1), 
                              L = L, nrep = shape * nrep, progress = FALSE, parallel = FALSE)
  } else {
    Gauss.all <- GeoSim(coordx = rbind(coord_obs, loc), 
                        corrmodel = corrmodel, model = "Gaussian", 
                        param = c(corr_param, mean = 0, nugget = 0, sill = 1), 
                        nrep = shape * nrep, progress = FALSE)
  }
  
  pair_idx <- split(seq_len(shape * nrep), ceiling(seq_len(shape * nrep)/shape))
  
  # Funzione helper per Gibbs sampler - evita duplicazione codice
  run_gibbs_sampler <- function(data_transformed, M0, shape, Sigma.mat.inv, sqrt_factor) {
    # Genera segni casuali e valori iniziali
    random_sign <- matrix(ifelse(runif(M0 * shape) <= 0.5, 1, -1), nrow = M0, ncol = shape)
    ini_vals <- random_sign * matrix(rep(sqrt(data_transformed), shape) * sqrt_factor, 
                                     nrow = M0, ncol = shape)
  gibbs.result <- dotCall64::.C64("gamma_gibbs_sampler",
                                      SIGNATURE = c("double","double","int","int","double","int","int"),
                                      SigmaInv = as.double(Sigma.mat.inv),
                                      y        = as.double(data_transformed),
                                      n        = as.integer(M0),
                                      v        = as.integer(shape),
                                      U        = as.double(ini_vals),
                                      nIte     = as.integer(1000),
                                      nRep     = as.integer(150),  # Unificato a 150
                                      INTENT   = c("r","r","r","r","rw","r","r"), 
                                      VERBOSE = 0, 
                                      NAOK = FALSE, 
                                      PACKAGE = "GeoModels")
    
    if (is.null(gibbs.result$U)) return(NULL)
    
    return(matrix(gibbs.result$U, nrow = M0, ncol = shape))
  }
  
  # Funzione helper per processamento vettorizzato
  process_simulation_vectorized <- function(gibbs_U, idx_xy, Gauss_data, M0, N0, shape, Wglob) {
    # Pre-alloca matrice risultato
    res_i <- matrix(0, nrow = N0, ncol = shape)
    gauss_matrices <- lapply(idx_xy, function(k) Gauss_data[[k]])
    for (kk in seq_len(shape)) {
      simX <- gauss_matrices[[kk]]
      sim_train_x <- simX[seq_len(M0)]
      sim_valid_x <- simX[(M0 + 1):(N0 + M0)]
      sk_pred_x <- as.numeric(crossprod(Wglob, gibbs_U[, kk] - sim_train_x))
      res_i[, kk] <- sim_valid_x + sk_pred_x
    }
    
    return(res_i)
  }
  
  if (!parallel) {
    cat("Performing", nrep, "conditional simulations ...\n")
    if (progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = seq_len(nrep))
    }
    
    res <- vector("list", nrep)
    
    for (i in seq_len(nrep)) {
      # Gibbs sampling
      gibbs_U <- run_gibbs_sampler(data_transformed, M0, shape, Sigma.mat.inv, sqrt_factor)
      
      if (is.null(gibbs_U)) {
        warning(paste("Gibbs sampler failed at simulation", i))
        next
      }
      
      # Processamento vettorizzato
      idx_xy <- pair_idx[[i]]
      res_i <- process_simulation_vectorized(gibbs_U, idx_xy, Gauss.all$data, 
                                             M0, N0, shape, Wglob)
      
      # Calcolo finale ottimizzato
      res[[i]] <- exp_mm_half * rowSums(res_i^2)
      
      if (progress) {
        pb(sprintf("Simulation %d/%d completed", i, nrep))
      }
    }
  } else {
    cat("Performing", nrep, "conditional simulations using", ncores, "cores ...\n")
    if (progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = seq_len(nrep))
    }
    
    # Ottimizzazione memoria: reduce global size limit più conservativo
    old_limit <- options(future.globals.maxSize = 1000 * 1024^2)  # Ridotto a 1GB
    on.exit(options(old_limit), add = TRUE)
    
    # Variabili essenziali ottimizzate - solo quello che serve
    essential_vars <- list(
      shape = shape,
      data_transformed = data_transformed,
      M0 = M0,
      N0 = N0,
      Sigma.mat.inv = Sigma.mat.inv,
      Wglob = Wglob,
      exp_mm_half = exp_mm_half,
      sqrt_factor = sqrt_factor,
      pair_idx = pair_idx,
      Gauss_data = Gauss.all$data
    )
    
    # Worker function ottimizzata
    run_simulation_worker <- function(i, vars) {
      # Estrae variabili
      shape <- vars$shape
      data_transformed <- vars$data_transformed
      M0 <- vars$M0
      N0 <- vars$N0
      Sigma.mat.inv <- vars$Sigma.mat.inv
      Wglob <- vars$Wglob
      exp_mm_half <- vars$exp_mm_half
      sqrt_factor <- vars$sqrt_factor
      pair_idx <- vars$pair_idx
      Gauss_data <- vars$Gauss_data
      
      # Gibbs sampling usando la funzione helper
      gibbs_U <- run_gibbs_sampler(data_transformed, M0, shape, Sigma.mat.inv, sqrt_factor)
      
      if (is.null(gibbs_U)) {
        warning(paste("Gibbs sampler failed at simulation ", i))
        return(NULL)
      }
      
      # Processamento vettorizzato
      idx_xy <- pair_idx[[i]]
      res_i <- process_simulation_vectorized(gibbs_U, idx_xy, Gauss_data, 
                                             M0, N0, shape, Wglob)
      
      # Risultato finale
      return(exp_mm_half * rowSums(res_i^2))
    }
    
    # Parallelizzazione
    future::plan(future::multisession, workers = ncores)
    on.exit(future::plan(future::sequential), add = TRUE)
    
    res <- future.apply::future_lapply(seq_len(nrep), 
                                       function(i) {
                                         result <- run_simulation_worker(i, essential_vars)
                                         if (progress) {
                                           pb(sprintf("Parallel simulation %d/%d completed", i, nrep))
                                         }
                                         return(result)
                                       }, 
                                       future.seed = TRUE, 
                                       future.globals = FALSE)
    
    # Rimuove risultati nulli
    res <- res[!sapply(res, is.null)]
    

  }
  
  return(res)
}


################################################################################# 
### internal function: conditional simulation of Gaussian RF 
################################################################################# 
compute_cond_sim <- function(i, sim_data_list, krig_pred_vec, local_flag, 
                             weights_obj, n_obs, n_loc) {
  sim_i <- sim_data_list[[i]]
  sim_nc_obs_data <- sim_i[1L:n_obs]
  sim_nc_loc_data <- sim_i[(n_obs + 1L):(n_obs + n_loc)]
  if (local_flag) {
    sim_nc_obs_pred <- numeric(n_loc)
    
    for (j in seq_len(n_loc)) {
      wj <- weights_obj[[j]]
      if (!length(wj$neighbor_indices)) next
      
      idx <- wj$neighbor_indices
      Wj  <- wj$weights
      len <- min(length(idx), length(Wj))
      
      if (len > 0L) {
        sim_nc_obs_pred[j] <- sum(Wj[1:len] * sim_nc_obs_data[idx[1:len]])
      }
    }
  } else {
    W <- weights_obj
    xvec <- as.vector(sim_nc_obs_data)  
    
    if (ncol(W) != length(xvec)) {
      sim_nc_obs_pred <- as.vector(t(W) %*% xvec)
    } else {
      sim_nc_obs_pred <- as.vector(W %*% xvec)
    }
  }
  
  krig_pred_vec + sim_nc_loc_data - sim_nc_obs_pred
}

################################################################################# 
### main function: conditional Gaussian simulations
################################################################################# 
Gauss_cd <- function(data, corrmodel, nrep, method, L,
                     param,
                     coord_obs, loc, coordt_use, time, X, Xloc, Mloc, distance, radius,
                     local, neighb, maxdist, maxtime,
                     space, spacetime, bivariate, parallel, ncores, progress) {
  #############################################
  # Build combined coordinates
  #############################################
  coord_sim <- rbind(coord_obs, loc)
  time_sim  <- c(coordt_use, time)
  X_sim     <- rbind(X, Xloc)
  n_obs <- nrow(coord_obs)
  n_loc <- nrow(loc)

  # Unconditional simulations
  cat("Performing", nrep, "unconditional simulations ...\n")
  if (method == "Cholesky") {
    sim_args <- list(coordx = coord_sim, coordt = time_sim, corrmodel = corrmodel, progress = FALSE,
                     X = X_sim, nrep = nrep, distance = distance, radius = radius)
    sim_args <- c(sim_args, list(param = param))
    sim_nc <- do.call(GeoSim, sim_args)
  }
  if (method == "TB" || method == "CE") {
    sim_args_approx <- list(coordx = coord_sim, coordt = time_sim, corrmodel = corrmodel, progress = FALSE,
                            method = method, L = L, parallel = parallel,
                            X = X_sim, nrep = nrep, distance = distance, radius = radius)
    sim_args_approx <- c(sim_args_approx, list(param = param))
    sim_nc <- do.call(GeoSimapprox, sim_args_approx)
  }

  #############################################
  # Kriging predictions
  #############################################
  krig_sim_args <- list(coordx = coord_obs, coordt = coordt_use,
                        data = data, corrmodel = corrmodel,
                        loc = loc, X = X, Xloc = Xloc, Mloc = Mloc, time = time,
                        distance = distance, radius = radius)
  krig_sim_args <- c(krig_sim_args, list(param = param))

  if (local) {
    if (!is.null(neighb))   krig_sim_args$neighb  <- neighb
    if (!is.null(maxdist))  krig_sim_args$maxdist <- maxdist
    if (!is.null(maxtime))  krig_sim_args$maxtime <- maxtime
    krig_sim_args$progress <- FALSE
    krig_sim <- do.call(GeoKrigloc, krig_sim_args)
  } else {
    krig_sim <- do.call(GeoKrig, krig_sim_args)
  }

  #############################################
  # Pre-compute kriging weights
  #############################################
  if(!local) cat("Computing kriging weights ...\n")
   else  cat("Computing local kriging weights ...\n")
  weights_args <- list(coordx = coord_obs, coordt = coordt_use,
                       corrmodel = corrmodel, loc = loc, 
                       X = X, Xloc = Xloc, Mloc = Mloc, time = time,
                       distance = distance, radius = radius)
  weights_args <- c(weights_args, list(param = param))

  if (local) {
    if (!is.null(neighb))   weights_args$neighb  <- neighb
    if (!is.null(maxdist))  weights_args$maxdist <- maxdist
    if (!is.null(maxtime))  weights_args$maxtime <- maxtime
    if(parallel) weights_args$parallel <- parallel
    krig_weights_obj <- do.call(GeoKriglocWeights, weights_args)
    weights_list <- krig_weights_obj$weights
  } else {
    krig_weights_obj <- do.call(GeoKrigWeights, weights_args)
    W <- krig_weights_obj$weights
  }

  #############################################
  # Conditional simulations with progress bar
  #############################################
  cat("Performing", nrep, "conditional simulations ...\n")

  # Setup progressr
  if (progress) {
    progressr::handlers(global = TRUE)
    progressr::handlers("txtprogressbar")
    p <- progressr::progressor(along = seq_len(nrep))
  } else {
    p <- function(...) NULL
  }

  if (space) {
    if (parallel && !is.null(ncores) && nrep > 1) {
      memory_needed <- object.size(sim_nc$data) + 
                       object.size(if(local) weights_list else W) + 
                       object.size(krig_sim$pred)
      memory_gb <- as.numeric(memory_needed) / 1024^3
      
      if (memory_gb > 10) {
        # Chunked parallel
        chunk_size <- max(1, min(20, floor(nrep/4)))
        chunks <- split(seq_len(nrep), ceiling(seq_len(nrep) / chunk_size))
        sim_cond <- vector("list", nrep)

        for (chunk_idx in seq_along(chunks)) {
          chunk_indices <- chunks[[chunk_idx]]
          cl <- parallel::makeCluster(min(ncores, length(chunk_indices)))
          parallel::clusterExport(cl, c("compute_cond_sim"), envir = environment())

          chunk_sim_data_only <- sim_nc$data[chunk_indices]
          krig_pred_only <- krig_sim$pred
          weights_only <- if(local) weights_list else W

          chunk_results <- parallel::parLapply(cl, seq_along(chunk_indices), function(local_i) {
            res <- compute_cond_sim(local_i, chunk_sim_data_only, krig_pred_only, local,
                                    weights_only, n_obs, n_loc)
            p()  # progress update
            res
          })
          parallel::stopCluster(cl)

          for (local_i in seq_along(chunk_indices)) {
            sim_cond[[chunk_indices[local_i]]] <- chunk_results[[local_i]]
          }
          rm(chunk_sim_data_only, krig_pred_only, weights_only, chunk_results)
          gc()
        }
        
      } else {
        # Standard parallel
        if (is.null(ncores)) {
          coremax <- parallel::detectCores()
          ncores  <- if (is.na(coremax) || coremax < 2) 1 else max(1, coremax - 1)
        }
        cl <- parallel::makeCluster(ncores)
        on.exit(parallel::stopCluster(cl), add = TRUE)

        parallel::clusterExport(cl, c("compute_cond_sim"), envir = environment())

        sim_data_list <- sim_nc$data
        krig_pred_vec <- krig_sim$pred
        weights_obj <- if(local) weights_list else W

        sim_cond <- parallel::parLapply(cl, seq_len(nrep), function(i) {
          res <- compute_cond_sim(i, sim_data_list, krig_pred_vec, local,
                                  weights_obj, n_obs, n_loc)
          p()  # progress update
          res
        })
      }
    } else {
      # Sequential
      sim_cond <- vector("list", nrep)
      for (i in seq_len(nrep)) {
        sim_cond[[i]] <- compute_cond_sim(i, sim_nc$data, krig_sim$pred, local,
                                         if(local) weights_list else W, n_obs, n_loc)
        p()  # progress update
      }
    }
  } else {
    stop("Gauss_cd currently not implemented for spacetime models ")
  }
  return(sim_cond)
}
###############################################################
############### end helper functions ##########################
###############################################################

###############################################################
############################ MAIN FUNCTION ####################
###############################################################
  call <- match.call()

  ###### Check GeoFit object ######
  if (!is.null(estobj)) {
    if (!inherits(estobj, "GeoFit"))
      stop("need a 'GeoFit' object as input\n")

    data <- estobj$data
    if (!estobj$grid) {
      if (!estobj$bivariate) {
        if (is.null(estobj$coordx_dyn)) coordx <- cbind(estobj$coordx, estobj$coordy, estobj$coordz)
        else coordx <- NULL
      } else {
        if (is.null(estobj$coordx_dyn)) {
          coordx <- estobj$coordx[1:estobj$ns[1]]
          coordy <- estobj$coordy[1:estobj$ns[2]]
        } else {
          coordx <- NULL
        }
      }
    } else {
      coordx <- estobj$coordx; coordy <- estobj$coordy; coordz <- estobj$coordz
    }
    if (length(estobj$coordt) == 1) coordt <- NULL else coordt <- estobj$coordt
    coordx_dyn <- estobj$coordx_dyn
    corrmodel  <- estobj$corrmodel
    model      <- estobj$model
    distance   <- estobj$distance
    grid       <- estobj$grid
    n          <- estobj$n
    param      <- append(estobj$param, estobj$fixed)
    radius     <- estobj$radius
    copula     <- estobj$copula
    anisopars  <- estobj$anisopars
    X          <- estobj$X
  }
  #### end check

  ### Additional checks
  if (is.null(CkModel(model))) stop("The name of the model is not correct\n")
  if (!is.character(corrmodel) || is.null(CkCorrModel(corrmodel)))
    stop("The name of the correlation model is wrong\n")
  if (!(method %in% c("Cholesky", "TB", "CE"))) {
    stop("The method of unconditional simulation is not correct\n")
  }
  if (local && is.null(neighb)) stop("For local kriging you need to specify neighb \n")
  corrmodel <- gsub("[[:blank:]]", "", corrmodel)
  model     <- gsub("[[:blank:]]", "", model)
  distance  <- gsub("[[:blank:]]", "", distance)

  #### checking cores ######
  coremax <- parallel::detectCores()
  if (is.na(coremax) || coremax == 1) {
    parallel <- FALSE
  } else {
    if (is.null(parallel)) parallel <- TRUE
    if (is.null(ncores)) ncores <- max(1, coremax - 1)
  }

  # Determine model type (space / spacetime / bivariate)
  model_ck  <- CkCorrModel(corrmodel)
  bivariate <- CheckBiv(model_ck)
  spacetime <- CheckST(model_ck)
  space     <- !spacetime && !bivariate

  # Observed coordinates matrix
  if (is.matrix(coordx)) {
    coord_obs <- coordx
  } else {
    coord_obs <- data.frame(x = coordx)
    if (!is.null(coordy)) coord_obs$y <- coordy
    if (!is.null(coordz)) coord_obs$z <- coordz
    coord_obs <- as.matrix(coord_obs)
  }

  # Covariate checks and conversions
  if (!is.null(Xloc)) {
    if (is.vector(Xloc)) {
      Xloc <- matrix(Xloc, nrow = 1)
    } else if (!is.matrix(Xloc)) {
      Xloc <- as.matrix(Xloc)
    }
  }
  if (!is.null(X) && !is.matrix(X)) X <- as.matrix(X)
  if (is.matrix(X) && is.null(Xloc)) stop("Covariates for prediction locations are missing \n")
  if (is.null(X) && is.matrix(Xloc)) stop("Covariates are missing \n")
  if (CheckST(CkCorrModel(corrmodel))) {
    if (is.null(time)) stop("At least one temporal instant is needed for space-time kriging\n ")
  }
  if (!is.null(Mloc) && !is.null(Xloc)) stop("Either Mloc or Xloc must be fixed\n")
  if ((length(param$mean) > 1) && is.null(Mloc)) stop("Mloc must be fixed \n")
  loc_orig <- loc
  if (!is.null(anisopars)) {
    loc <- GeoAniso(loc, c(anisopars$angle, anisopars$ratio))
  }
  # Temporal coordinates to use
  coordt_use <- if (spacetime || bivariate) coordt else NULL

############################################################################
######################### not copula models ################################
############################################################################
if (is.null(copula)) {
  if (model == "Gaussian") {
    res <- Gauss_cd(data, corrmodel, nrep, method, L,
                    param,
                    coord_obs, loc, coordt_use, time, X, Xloc, Mloc, distance, radius,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores,progress)
  }

  ################################################################################
  ############### monotone transformations of one Gaussian RF ####################
  ################################################################################
  if (model == "LogGaussian") {
    mm <- exp(param$mean)
    vv <- param$sill
    datanorm <- (log(data/mm) + vv/2)/sqrt(vv)   # gaussian scale
    param$mean <- 0; param$sill <- 1
    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores,progress)
    res <- lapply(res, function(r) mm * exp(sqrt(vv) * r - vv/2)) # back-transformation
  }

  if (model == "Tukeyh") {
    inverse_lamb <- function(x, tail) {
      value <- sqrt(VGAM::lambertW(tail * x * x)/tail)
      return(sign(x) * value)
    }
    mm <- param$mean
    vv <- param$sill
    tail <- param$tail
    datanorm <- inverse_lamb((data - mm)/sqrt(vv), tail)
    param$mean <- 0; param$sill <- 1; param$tail <- NULL
    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores,progress)
    res <- lapply(res, function(r) mm + sqrt(vv) * r * exp(tail * r^2/2))
  }
##################################################################################################
  if (model == "Tukeyh2") {
    inverse_lamb <- function(x, tail) {
      value <- sqrt(VGAM::lambertW(tail * x * x)/tail)
      return(sign(x) * value)
    }
    mm <- param$mean
    vv <- param$sill
    tail1 <- param$tail1
    tail2 <- param$tail2
    inverse_lamb2 <- function(x, tail1, tail2) {
      aa <- 1:length(x)
      sel1 <- I(x >= 0) * aa
      sel2 <- I(x < 0) * aa
      a1 <- inverse_lamb(x[sel1], tail1)
      b1 <- inverse_lamb(x[sel2], tail2)
      sel1[sel1 > 0] <- a1
      sel2[sel2 > 0] <- b1
      return(c(sel1 + sel2))
    }
    datanorm <- inverse_lamb2((data - mm)/sqrt(vv), tail1, tail2)
    param$mean <- 0; param$sill <- 1; param$tail1 <- NULL; param$tail2 <- NULL
    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores,progress)
    ff <- function(r, mm, vv, tail1, tail2) {
      aa <- 1:length(r)
      sel1 <- I(r >= 0) * aa
      sel2 <- I(r < 0) * aa
      res1 <- mm + sqrt(vv) * r[sel1] * exp(tail1 * r[sel1]^2/2)
      res2 <- mm + sqrt(vv) * r[sel2] * exp(tail2 * r[sel2]^2/2)
      sel1[sel1 > 0] <- res1
      sel2[sel2 > 0] <- res2
      return(c(sel1 + sel2))
    }
    res <- lapply(res, ff, mm = mm, vv = vv, tail1 = tail1, tail2 = tail2)
  }
##################################################################################################
  if (model == "SinhAsinh") {
    mm <-param$mean
    vv <- param$sill
    skew <- param$skew
    tail <- param$tail
  datanorm <- sinh((asinh((data - mm)/sqrt(vv)) / tail) - skew)
    param$mean <- 0; param$sill <- 1; param$skew <- NULL; param$tail <- NULL
    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores,progress)
    res <- lapply(res, function(r) mm + sqrt(vv) * sinh(tail * (asinh(r) + skew)))
  }
  ################################################################################
  ############### non monotone transformations of more thhan one Gaussian RF #####
  ################################################################################
  if(model=="SkewGaussian") 
  {    res <- SkewGaussianSimcond(coord_obs, loc, data, param, corrmodel,local=local,neighb=neighb,L=L,
                       nrep = nrep, method = method, parallel = parallel, ncores = ncores,progress=progress)
  }

##################################################################################################
  if (model == "Gamma")        { 
    #Made by David Rivas
    res <- GammaSimcond(coord_obs, loc, data, param, corrmodel,local=local,neighb=neighb,L=L,
                       nrep = nrep, method = method, parallel = parallel, ncores = ncores,progress=progress)
  }
  if (model == "Weibull")      {
    #Made by David Rivas
    mm <- param$mean
    param.weibull <- param
    param.weibull$shape <- 2
    param.weibull$mean <- 0
    ck <- gamma( (param$shape + 1)/param$shape )
    data.t <- (exp(-mm)*data*ck)^param$shape
    res <- GammaSimcond(coord_obs, loc, data.t, param.weibull, corrmodel,local=local,neighb=neighb,L=L,
                        nrep = nrep, method = method, parallel = parallel, ncores = ncores,progress=progress)
    res <- lapply(res, function(x) exp(mm)*(x^(1/param$shape))/ck )
  }
  if (model == "Poisson")      { stop("Poisson  not implemented...") }
  if (model == "Binomial")     { stop("Binomial  not implemented...") }

} else { 
###########################################################################
################## copula models ###########################################
###########################################################################

####################### gaussian copula ####################################################
if (copula == "Gaussian") {

  if (model == "Weibull") {
    mm <- exp(param$mean)
    ss <- param$shape
    datanorm1 <- pweibull(data, shape = ss, scale = mm/(gamma(1 + 1/ss)))
    datanorm  <- qnorm(datanorm1)
    param$mean <- 0; param$sill <- 1; param$shape <- NULL
    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores,progress)
    res <- lapply(res, function(r) pnorm(r))
    res <- lapply(res, function(r) qweibull(r, shape = ss, scale = mm/(gamma(1 + 1/ss))))
  }

  if (model == "LogGaussian") {
    mm <- exp(param$mean)
    vv <- param$sill
    datanorm1 <- plnorm(data, meanlog = mm - vv/2, sdlog = sqrt(vv), log.p = FALSE)
    datanorm  <- qnorm(datanorm1)
    param$mean <- 0; param$sill <- 1
    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores,progress)
    res <- lapply(res, function(r) pnorm(r))
    res <- lapply(res, function(r) qlnorm(r, meanlog = mm - vv/2, sdlog = sqrt(vv), log.p = FALSE))
  }

  if (model == "Gamma") {
    mm <- exp(param$mean)
    ss <- param$shape
    datanorm1 <- pgamma(data, shape = ss/2, rate = ss/(2 * mm))
    datanorm  <- qnorm(datanorm1)
    param$mean <- 0; param$sill <- 1; param$shape <- NULL
    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores,progress)
    res <- lapply(res, function(r) pnorm(r))
    res <- lapply(res, function(r) qgamma(r, shape = ss/2, rate = ss/(2 * mm)))
  }

  if (model == "Beta2") {
    mm <- 1/(1 + exp(-param$mean))
    ss <- param$shape
    pmin <- param$min
    pmax <- param$max
    datanorm1 <- pbeta((data - pmin)/(pmax - pmin), shape1 = mm * ss, shape2 = (1 - mm) * ss)
    datanorm  <- qnorm(datanorm1)
    param$mean <- 0; param$sill <- 1; param$shape <- NULL; param$min <- NULL; param$max <- NULL
    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores,progress)
    res <- lapply(res, function(r) pnorm(r))
    res <- lapply(res, function(r) pmin + (pmax - pmin) * qbeta(r, shape1 = mm * ss, shape2 = (1 - mm) * ss))
  }

  if (model == "StudentT") {
    mm <- param$mean
    vv <- param$sill
    df <- param$df
    datanorm1 <- pt((data - mm)/sqrt(vv), df = 1/df)
    datanorm  <- qnorm(datanorm1)
    param$mean <- 0; param$sill <- 1; param$df <- NULL
    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores,progress)
    res <- lapply(res, function(r) pnorm(r))
    res <- lapply(res, function(r) mm + sqrt(vv) * qt(r, df = 1/df))
  }

  if (model == "Gaussian") {
    mm <- param$mean
    vv <- param$sill
    datanorm1 <- pnorm(data, mm, sqrt(vv))
    datanorm  <- qnorm(datanorm1)
    param$mean <- 0; param$sill <- 1
    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores,progress)
    res <- lapply(res, function(r) pnorm(r))
    res <- lapply(res, function(r) qnorm(r, mm, sqrt(vv)))
  }

  if (model == "Binomial") {
    mm <- param$mean; prob <- pnorm(mm)
    datanorm1 <- (pbinom(data - 1, size = n, prob = prob) + pbinom(data, size = n, prob = prob))/2
    datanorm  <- qnorm(datanorm1)
    param$mean <- 0; param$sill <- 1
    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores,progress)
    res <- lapply(res, function(r) pnorm(r))
    res <- lapply(res, function(r) qbinom(p = r, size = n, prob = prob))
  }

  if (model == "BinomialNeg") {
    mm <- param$mean; prob <- pnorm(mm); size <- n
    datanorm1 <- (pnbinom(data - 1, size = size, prob = prob) + pnbinom(data, size = size, prob = prob))/2
    datanorm  <- qnorm(datanorm1)
    param$mean <- 0; param$sill <- 1
    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores,progress)
    res <- lapply(res, function(r) pnorm(r))
    res <- lapply(res, function(r) qnbinom(p = r, size = size, prob = prob))
  }

  if (model == "Poisson") {
    mu <- exp(param$mean)
    datanorm1 <- (ppois(data - 1, lambda = mu) + ppois(data, lambda = mu))/2
    datanorm  <- qnorm(datanorm1)
    param$mean <- 0; param$sill <- 1
    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores,progress)
    res <- lapply(res, function(r) pnorm(r))
    res <- lapply(res, function(r) qpois(r, lambda = mu))
  }

}
####################### skewgaussian copula ####################################################
if (copula == "SkewGaussian") {

    if (model == "Weibull") {

    param$skew=param$nu
    mm <- exp(param$mean)
    ss <- param$shape
    datanorm1 <- pweibull(data, shape = ss, scale = mm/(gamma(1 + 1/ss))) 
    param$mean <- 0; param$sill <- 1; param$shape <- NULL

    omega=as.numeric(sqrt((param$nu^2 + param$sill)/param$sill))
    alpha=as.numeric(param$nu/param$sill^0.5)
    datanorm  <- sn::qsn(datanorm1, xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))

    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param, corrmodel,local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, parallel = parallel, ncores = ncores,progress=progress)

    res=lapply(res, function(r) sn::psn(r ,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha)))
    res <- lapply(res, function(r) qweibull(r, shape = ss, scale = mm/(gamma(1 + 1/ss))))
  }

    if (model == "LogGaussian") {
    param$skew=param$nu  
    mm <- exp(param$mean)
    vv <- param$sill
    datanorm1 <- plnorm(data, meanlog = mm - vv/2, sdlog = sqrt(vv), log.p = FALSE)
    param$mean <- 0;param$sill<-1
    omega=as.numeric(sqrt((param$nu^2 + param$sill)/param$sill))
    alpha=as.numeric(param$nu/param$sill^0.5)
    datanorm  <- sn::qsn(datanorm1, xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))

    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param, corrmodel,local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, parallel = parallel, ncores = ncores,progress=progress)
    res=lapply(res, function(r) sn::psn(r ,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha)))
    res <- lapply(res, function(r) qlnorm(r, meanlog = mm - vv/2, sdlog = sqrt(vv), log.p = FALSE))
  }

  if (model == "Gamma") {
    param$skew=param$nu 
    mm <- exp(param$mean)
    ss <- param$shape
    datanorm1 <- pgamma(data, shape = ss/2, rate = ss/(2 * mm))
    param$mean <- 0; param$sill <- 1; param$shape <- NULL
    omega=as.numeric(sqrt((param$nu^2 + param$sill)/param$sill))
    alpha=as.numeric(param$nu/param$sill^0.5)
    datanorm  <- sn::qsn(datanorm1, xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))
    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param, corrmodel,local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, parallel = parallel, ncores = ncores,progress=progress)
    res=lapply(res, function(r) sn::psn(r ,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha)))
    res <- lapply(res, function(r) qgamma(r, shape = ss/2, rate = ss/(2 * mm)))
  }

    if (model == "Beta2") {
    param$skew=param$nu 
    mm <- 1/(1 + exp(-param$mean))
    ss <- param$shape
    pmin <- param$min
    pmax <- param$max
    datanorm1 <- pbeta((data - pmin)/(pmax - pmin), shape1 = mm * ss, shape2 = (1 - mm) * ss)
    param$mean <- 0; param$sill <- 1;
    omega=as.numeric(sqrt((param$nu^2 + param$sill)/param$sill))
    alpha=as.numeric(param$nu/param$sill^0.5)
    datanorm  <- sn::qsn(datanorm1, xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))
    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param, corrmodel, local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, parallel = parallel, ncores = ncores,progress=progress)
    res=lapply(res, function(r) sn::psn(r ,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha)))
    res <- lapply(res, function(r) pmin + (pmax - pmin) * qbeta(r, shape1 = mm * ss, shape2 = (1 - mm) * ss))
  }

    if (model == "StudentT") {
    param$skew=param$nu 
    mm <- param$mean
    vv <- param$sill
    df <- param$df
    datanorm1 <- pt((data - mm)/sqrt(vv), df = 1/df)
    param$mean <- 0; param$sill <- 1; param$df <- NULL
    omega=as.numeric(sqrt((param$nu^2 + param$sill)/param$sill))
    alpha=as.numeric(param$nu/param$sill^0.5)
    datanorm  <- sn::qsn(datanorm1, xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))
    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param, corrmodel, local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, parallel = parallel, ncores = ncores,progress=progress)
    res=lapply(res, function(r) sn::psn(r ,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha)))
    res <- lapply(res, function(r) mm + sqrt(vv) * qt(r, df = 1/df))
  }

    if (model == "Gaussian") {
    param$skew=param$nu 
    mm <- param$mean
    vv <- param$sill
    datanorm1 <- pnorm(data, mm, sqrt(vv))
    param$mean <- 0; param$sill <- 1
    omega=as.numeric(sqrt((param$nu^2 + param$sill)/param$sill))
    alpha=as.numeric(param$nu/param$sill^0.5)
    datanorm  <- sn::qsn(datanorm1, xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))
    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param, corrmodel, local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, parallel = parallel, ncores = ncores,progress=progress)
    res=lapply(res, function(r) sn::psn(r ,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha)))
    res <- lapply(res, function(r) qnorm(r, mm, sqrt(vv)))
  }

    if (model == "Binomial") {
    param$skew=param$nu   
    mm <- param$mean; prob <- pnorm(mm)
    datanorm1 <- (pbinom(data - 1, size = n, prob = prob) + pbinom(data, size = n, prob = prob))/2
    param$mean <- 0; param$sill <- 1
    omega=as.numeric(sqrt((param$nu^2 + param$sill)/param$sill))
    alpha=as.numeric(param$nu/param$sill^0.5)
    datanorm  <- sn::qsn(datanorm1, xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))
    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param, corrmodel, local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, parallel = parallel, ncores = ncores,progress=progress)
    res=lapply(res, function(r) sn::psn(r ,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha)))
    res <- lapply(res, function(r) qbinom(p = r, size = n, prob = prob))
  }

  if (model == "BinomialNeg") {
    param$skew=param$nu  
    mm <- param$mean; prob <- pnorm(mm); size <- n
    datanorm1 <- (pnbinom(data - 1, size = size, prob = prob) + pnbinom(data, size = size, prob = prob))/2
    param$mean <- 0; param$sill <- 1
    omega=as.numeric(sqrt((param$nu^2 + param$sill)/param$sill))
    alpha=as.numeric(param$nu/param$sill^0.5)
    datanorm  <- sn::qsn(datanorm1, xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))
    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param, corrmodel, local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, parallel = parallel, ncores = ncores,progress=progress)
    res=lapply(res, function(r) sn::psn(r ,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha)))
    res <- lapply(res, function(r) qnbinom(p = r, size = size, prob = prob))
  }

  if (model == "Poisson") {
    param$skew=param$nu  
    mu <- exp(param$mean)
    datanorm1 <- (ppois(data - 1, lambda = mu) + ppois(data, lambda = mu))/2
    param$mean <- 0; param$sill <- 1
    omega=as.numeric(sqrt((param$nu^2 + param$sill)/param$sill))
    alpha=as.numeric(param$nu/param$sill^0.5)
    datanorm  <- sn::qsn(datanorm1, xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))
    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param, corrmodel, local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, parallel = parallel, ncores = ncores,progress=progress)
    res=lapply(res, function(r) sn::psn(r ,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha)))
    res <- lapply(res, function(r) qpois(r, lambda = mu))
  }

}
###########################################################################
}


sim_result <- do.call(rbind, res)
cond_mean=colMeans(sim_result)  # conditional mean
cond_var =apply(sim_result,2, var) # conditional var
############################################################################
############################################################################
############################ OUTPUT ########################################
GeoSimcond_out <- list(
    bivariate = bivariate,
    coordx = coordx,
    coordy = coordy,
    coordz = coordz,
    coordt = coordt,
    corrmodel = corrmodel,
    data = data,
    distance = distance,
    grid = grid,
    loc = loc,
    copula = copula,
    numcoord = nrow(coord_obs),
    numloc = nrow(loc),
    numtime = length(coordt),
    numt = length(time),
    maxdist = maxdist,
    maxtime = maxtime,
    model = model,
    n = n,
    param = param,
    condsim = res,            ##### results!
    cond_mean=cond_mean,     ###  conditional mean
    cond_var=cond_var,       ### conditioonal variance
    radius = radius,
    spacetime = spacetime,
    time = time
  )
  structure(c(GeoSimcond_out, call = call), class = c("GeoSimcond"))
}