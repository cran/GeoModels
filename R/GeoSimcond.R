GeoSimcond <- function(estobj = NULL, data, coordx, coordy = NULL, coordz = NULL, coordt = NULL,
                       coordx_dyn = NULL, corrmodel, distance = "Eucl",
                       grid = FALSE, loc, maxdist = NULL, maxtime = NULL, method = "Cholesky",
                       model = "Gaussian", n = 1, nrep = 1, local = FALSE, L = 1000,
                       neighb = NULL, param, anisopars = NULL, radius = 1, sparse = FALSE,
                       time = NULL, copula = NULL, X = NULL, Xloc = NULL, Mloc = NULL,
                       parallel = FALSE, ncores = NULL, progress = FALSE) {

################################################################################# 
### Helper function for conditional skew gaussian simulation 
################################################################################# 
SkewGaussianSimcond <- function(
    coord_obs, loc, data, param, corrmodel,
    nrep = 1, method = "Cholesky", local = FALSE,
    neighb = NULL, L = NULL,
    distance = "Eucl", radius = 1, anisopars = NULL,
    parallel = FALSE, ncores = 1, progress = FALSE,
    n_iter = 1000L, batch_size = NULL, max_batch_mb = 256
) {

  ##########################################################################
  ## Parameters and checks
  ##########################################################################

  mm         <- as.numeric(param$mean)
  sill       <- as.numeric(param$sill)
  skew       <- as.numeric(param$skew)
  sigma_sqrt <- sqrt(sill)

  if (length(mm) != 1L || !is.finite(mm)) {
    stop("SkewGaussianSimcond currently requires a scalar mean.")
  }

  if (length(sill) != 1L || !is.finite(sill) || sill <= 0) {
    stop("param$sill must be strictly positive.")
  }

  if (length(skew) != 1L || !is.finite(skew)) {
    stop("param$skew must be finite.")
  }

  if (length(nrep) != 1L || !is.finite(nrep) || nrep < 1) {
    stop("nrep must be at least 1.")
  }

  if (length(n_iter) != 1L || !is.finite(n_iter) || n_iter < 1) {
    stop("n_iter must be at least 1.")
  }

  if (local) {
    stop(
      "Local prediction cannot be used with the ",
      "skew-Gaussian Gibbs sampler."
    )
  }

  if (!(method %in% c("Cholesky", "TB", "CE"))) {
    stop("method must be one of 'Cholesky', 'TB', or 'CE'.")
  }

  nrep  <- as.integer(nrep)
  n_iter <- as.integer(n_iter)

  coord_obs <- as.matrix(coord_obs)
  loc       <- as.matrix(loc)
  data      <- as.numeric(data)

  M0 <- nrow(coord_obs)
  N0 <- nrow(loc)

  if (length(data) != M0) {
    stop(
      "The length of data must match the number ",
      "of observed locations."
    )
  }

  if (N0 < 1L) {
    stop("loc must contain at least one prediction location.")
  }

  ##########################################################################
  ## Standardized latent representation
  ##
  ## Y(s) = mm + skew |G1(s)| + sqrt(sill) G2(s)
  ##
  ## Z(s) = gamma |G1(s)| + G2(s),
  ## gamma = skew / sqrt(sill)
  ##########################################################################

  gamma <- skew / sigma_sqrt
  z_obs <- (data - mm) / sigma_sqrt

  if (any(!is.finite(z_obs))) {
    stop("The standardized observations contain non-finite values.")
  }

  gaussian_case <- isTRUE(gamma == 0)

  keep       <- CorrParam(corrmodel)
  corr_param <- param[keep]

  latent_param <- c(
    corr_param,
    mean = 0,
    nugget = 0,
    sill = 1
  )

  ##########################################################################
  ## Gaussian covariance matrix and kriging weights
  ##########################################################################

  cat("Computing Kriging weights ...\n")

  GeoW <- GeoKrigWeights(
    coordx = coord_obs,
    corrmodel = corrmodel,
    loc = loc,
    model = "Gaussian",
    distance = distance,
    radius = radius,
    anisopars = anisopars,
    param = latent_param
  )

  Wglob <- GeoW$weights

  if (is.null(Wglob)) {
    stop("GeoKrigWeights did not return the kriging weights.")
  }

  if (is.null(dim(Wglob))) {
    Wglob <- matrix(Wglob, nrow = M0)
  }

  if (nrow(Wglob) == M0) {
    weights_orientation <- 1L
  } else if (ncol(Wglob) == M0) {
    weights_orientation <- 2L
  } else {
    stop(
      "The dimensions of the kriging weights ",
      "are not compatible with the data."
    )
  }

  # The precision matrix is unnecessary when gamma = 0.
  if (gaussian_case) {

    Sigma.mat.inv <- numeric(0L)

  } else {

    if (is.null(GeoW$covmatrix)) {
      stop(
        "GeoKrigWeights did not return the ",
        "latent covariance matrix."
      )
    }

    Sigma.mat.inv <- MatInv(mtx = GeoW$covmatrix)

    if (any(!is.finite(Sigma.mat.inv))) {
      stop(
        "The inverse latent covariance matrix ",
        "contains non-finite values."
      )
    }
  }

  all_coord <- rbind(coord_obs, loc)

  ##########################################################################
  ## Batch size
  ##########################################################################

  fields_per_rep <- if (gaussian_case) 1L else 2L

  if (is.null(batch_size)) {

    if (length(max_batch_mb) != 1L ||
        !is.finite(max_batch_mb) ||
        max_batch_mb <= 0) {
      stop("max_batch_mb must be strictly positive.")
    }

    bytes_per_rep <- fields_per_rep * (M0 + N0) * 8

    batch_size <- floor(
      max_batch_mb * 1024^2 / max(bytes_per_rep, 1)
    )

    batch_size <- max(
      1L,
      min(nrep, as.integer(batch_size))
    )

  } else {

    if (length(batch_size) != 1L ||
        !is.finite(batch_size) ||
        batch_size < 1) {
      stop("batch_size must be at least 1.")
    }

    batch_size <- max(
      1L,
      min(nrep, as.integer(batch_size))
    )
  }

  ##########################################################################
  ## Normalize the output of GeoSim or GeoSimapprox
  ##########################################################################

  normalize_simulations <- function(x, expected) {

    if (is.list(x)) {

      if (length(x) != expected) {
        stop(
          "Unexpected number of unconditional ",
          "Gaussian simulations."
        )
      }

      return(x)
    }

    if (is.matrix(x)) {

      if (ncol(x) == expected) {
        return(
          lapply(
            seq_len(expected),
            function(j) x[, j]
          )
        )
      }

      if (nrow(x) == expected) {
        return(
          lapply(
            seq_len(expected),
            function(j) x[j, ]
          )
        )
      }
    }

    if (expected == 1L) {
      return(list(as.numeric(x)))
    }

    stop(
      "Unexpected format of the unconditional ",
      "Gaussian simulations."
    )
  }

  ##########################################################################
  ## Generate one batch of unconditional Gaussian simulations
  ##########################################################################

  generate_batch <- function(current_nrep) {

    current_fields <- fields_per_rep * current_nrep

    if (method %in% c("TB", "CE")) {

      sim_obj <- GeoSimapprox(
        coordx = all_coord,
        corrmodel = corrmodel,
        method = method,
        model = "Gaussian",
        param = latent_param,
        L = L,
        nrep = current_fields,
        progress = FALSE,
        parallel = FALSE,
        distance = distance,
        radius = radius,
        anisopars = anisopars
      )

    } else {

      sim_obj <- GeoSim(
        coordx = all_coord,
        corrmodel = corrmodel,
        model = "Gaussian",
        param = latent_param,
        nrep = current_fields,
        progress = FALSE,
        distance = distance,
        radius = radius,
        anisopars = anisopars
      )
    }

    sim_data <- normalize_simulations(
      sim_obj$data,
      current_fields
    )

    if (gaussian_case) {

      tasks <- lapply(
        seq_len(current_nrep),
        function(j) {
          list(
            g2 = as.numeric(sim_data[[j]])
          )
        }
      )

    } else {

      tasks <- lapply(
        seq_len(current_nrep),
        function(j) {
          list(
            g1 = as.numeric(sim_data[[2L * j - 1L]]),
            g2 = as.numeric(sim_data[[2L * j]])
          )
        }
      )
    }

    rm(sim_obj, sim_data)

    tasks
  }

  ##########################################################################
  ## Self-contained worker
  ##
  ## Only the Gaussian simulations required by the current task are passed
  ## to the worker. The complete batch is not captured as a global object.
  ##########################################################################

  worker_fun <- function(
      task, z_obs, M0, N0, gamma, mm, skew, sigma_sqrt,
      Sigma_mat_inv, Wglob, weights_orientation,
      n_iter, gaussian_case
  ) {

    tryCatch({

      apply_weights <- function(delta) {

        if (weights_orientation == 1L) {
          as.numeric(crossprod(Wglob, delta))
        } else {
          as.numeric(Wglob %*% delta)
        }
      }

      conditional_field <- function(
          simulation,
          latent_obs
      ) {

        if (length(simulation) != M0 + N0) {
          stop(
            "An unconditional Gaussian simulation ",
            "has an incorrect length."
          )
        }

        sim_obs <- simulation[seq_len(M0)]
        sim_loc <- simulation[M0 + seq_len(N0)]

        sim_loc + apply_weights(
          latent_obs - sim_obs
        )
      }

      ######################################################################
      ## Gaussian special case
      ######################################################################

      if (gaussian_case) {

        g2_loc <- conditional_field(
          task$g2,
          z_obs
        )

        return(
          list(
            ok = TRUE,
            value = mm + sigma_sqrt * g2_loc
          )
        )
      }

      ######################################################################
      ## Initial state satisfying
      ##
      ## z_i = gamma |G1_i| + G2_i
      ######################################################################
      g1_init <- stats::rnorm(M0)
      g2_init <- z_obs - gamma * abs(g1_init)
      ######################################################################
      ## Exact Gibbs kernel implemented in C
      ######################################################################
      gibbs_result <- dotCall64::.C64(
        "skew_gaussian_gibbs_sampler",
        SIGNATURE = c(
          "double", "integer", "double", "double",
          "integer", "double", "double"
        ),
        data_obs = as.double(z_obs),
        n = as.integer(M0),
        Sigma_mat_inv = as.double(Sigma_mat_inv),
        eta = as.double(c(gamma, 0)),
        n_iter = as.integer(n_iter),
        data_x = as.double(g1_init),
        data_y = as.double(g2_init),
        INTENT = c(
          "r", "r", "r", "r",
          "r", "rw", "rw"
        ),
        VERBOSE = 0,
        NAOK = FALSE,
        PACKAGE = "GeoModels"
      )

      g1_obs <- as.numeric(gibbs_result$data_x)
      g2_obs <- as.numeric(gibbs_result$data_y)

      if (length(g1_obs) != M0 ||
          length(g2_obs) != M0 ||
          any(!is.finite(g1_obs)) ||
          any(!is.finite(g2_obs))) {
        stop(
          "The Gibbs sampler returned invalid ",
          "latent values."
        )
      }

      ######################################################################
      ## Check the nonlinear restriction
      ######################################################################
      constraint_error <- max(
        abs(
          z_obs -
            gamma * abs(g1_obs) -
            g2_obs
        )
      )
      constraint_tolerance <- 1e-7 * (
        1 + max(abs(z_obs))
      )

      if (!is.finite(constraint_error) ||
          constraint_error > constraint_tolerance) {
        stop(
          "The Gibbs output does not satisfy the ",
          "skew-Gaussian constraint."
        )
      }
      ######################################################################
      ## Gaussian substitution conditional simulations
      ######################################################################

      g1_loc <- conditional_field(
        task$g1,
        g1_obs
      )

      g2_loc <- conditional_field(
        task$g2,
        g2_obs
      )

      ######################################################################
      ## Back transformation
      ##
      ## Y(s) = mm + skew |G1(s)| + sqrt(sill) G2(s)
      ######################################################################

      list(
        ok = TRUE,
        value =
          mm +
          skew * abs(g1_loc) +
          sigma_sqrt * g2_loc
      )

    }, error = function(e) {

      list(
        ok = FALSE,
        error = conditionMessage(e)
      )
    })
  }

  ##########################################################################
  ## Parallel configuration
  ##########################################################################

  if (parallel && nrep > 1L) {

    if (is.null(ncores) ||
        length(ncores) != 1L ||
        !is.finite(ncores)) {
      coremax <- parallel::detectCores()
      ncores <- if (is.na(coremax) || coremax < 2L) {
        1L
      } else {
        coremax - 1L
      }
    }
    ncores <- max(
      1L,
      min(as.integer(ncores), nrep)
    )
    parallel <- ncores > 1L
  } else {
    parallel <- FALSE
    ncores <- 1L
  }
  if (parallel) {
    old_plan <- future::plan()

    on.exit(
      future::plan(old_plan),
      add = TRUE
    )

    old_limit <- getOption(
      "future.globals.maxSize"
    )

    static_size <- sum(
      as.numeric(object.size(Sigma.mat.inv)),
      as.numeric(object.size(Wglob)),
      as.numeric(object.size(z_obs))
    )

    required_limit <- max(
      500 * 1024^2,
      2 * static_size + max_batch_mb * 1024^2
    )

    if (!is.null(old_limit)) {
      required_limit <- max(
        required_limit,
        old_limit
      )
    }

    options(
      future.globals.maxSize = required_limit
    )
    on.exit(
      options(
        future.globals.maxSize = old_limit
      ),
      add = TRUE
    )
    future::plan(
      future::multisession,
      workers = ncores
    )
  }
  ##########################################################################
  ## Batch execution
  ##########################################################################

  run_batches <- function(progressor) {

    res <- vector("list", nrep)
    failure_messages <- character(0L)
    batch_starts <- seq.int(
      1L,
      nrep,
      by = batch_size
    )
    for (batch_start in batch_starts) {

      batch_end <- min(
        nrep,
        batch_start + batch_size - 1L
      )
      current_indices <- batch_start:batch_end
      current_nrep <- length(current_indices)
      tasks <- generate_batch(current_nrep)
      if (parallel && current_nrep > 1L) {

        batch_result <- future.apply::future_lapply(
          tasks,
          worker_fun,
          z_obs = z_obs,
          M0 = M0,
          N0 = N0,
          gamma = gamma,
          mm = mm,
          skew = skew,
          sigma_sqrt = sigma_sqrt,
          Sigma_mat_inv = Sigma.mat.inv,
          Wglob = Wglob,
          weights_orientation = weights_orientation,
          n_iter = n_iter,
          gaussian_case = gaussian_case,
          future.seed = TRUE,

          # One future chunk per worker. This avoids repeatedly
          # serializing Sigma.mat.inv and Wglob for every simulation.
          future.scheduling = 1,

          # worker_fun is self-contained and all data are passed
          # explicitly.
          future.globals = FALSE,

          future.packages = c(
            "GeoModels",
            "dotCall64"
          )
        )

      } else {

        batch_result <- lapply(
          tasks,
          worker_fun,
          z_obs = z_obs,
          M0 = M0,
          N0 = N0,
          gamma = gamma,
          mm = mm,
          skew = skew,
          sigma_sqrt = sigma_sqrt,
          Sigma_mat_inv = Sigma.mat.inv,
          Wglob = Wglob,
          weights_orientation = weights_orientation,
          n_iter = n_iter,
          gaussian_case = gaussian_case
        )
      }

      for (j in seq_along(batch_result)) {

        global_index <- current_indices[j]

        if (isTRUE(batch_result[[j]]$ok)) {

          res[[global_index]] <-
            batch_result[[j]]$value

        } else {

          failure_messages <- c(
            failure_messages,
            sprintf(
              "simulation %d: %s",
              global_index,
              batch_result[[j]]$error
            )
          )
        }

        progressor()
      }

      rm(tasks, batch_result)
      gc(verbose = FALSE)
    }

    ########################################################################
    ## Remove failed simulations
    ########################################################################

    failed <- vapply(
      res,
      is.null,
      logical(1L)
    )
    if (any(failed)) {
      details <- if (length(failure_messages)) {
        paste(
          utils::head(failure_messages, 3L),
          collapse = "; "
        )
      } else {
        "unknown error"
      }
      warning(
        sprintf(
          "%d simulations out of %d failed. %s",
          sum(failed),
          nrep,
          details
        ),
        call. = FALSE
      )
      res <- res[!failed]
    }
    if (!length(res)) {
      stop("All conditional simulations failed.")
    }

    res
  }
  ##########################################################################
  ## Run
  ##########################################################################
  cat(
    "Performing",
    nrep,
    "conditional simulations",
    if (parallel) {
      paste0(" using ", ncores, " cores")
    } else {
      ""
    },
    "...\n"
  )
  if (progress) {
    progressr::handlers(global = TRUE)
    progressr::handlers("txtprogressbar")
    res <- progressr::with_progress({
      p <- progressr::progressor(
        steps = nrep
      )
      run_batches(p)
    })
  } else {
    res <- run_batches(
      function(...) invisible(NULL)
    )
  }
  res
}
################################################################################# 
### Helper function for conditional Gamma simulation 
################################################################################# 
################################################################################# 
GammaSimcond <- function(coord_obs, loc, data, param, corrmodel,
                         nrep = 1, method = "Cholesky", local = FALSE, neighb = NULL,
                         L = NULL, distance = "Eucl", radius = 1, anisopars = NULL,
                         parallel = FALSE, ncores = 1, progress = FALSE) {
  mm <- param$mean
  shape <- round(param$shape)
  keep <- CorrParam(corrmodel)
  corr_param <- param[keep]
  
  cat("Computing Kriging weights ...\n")
  if (!local) 
    GeoW <- GeoKrigWeights(coordx = coord_obs, corrmodel = corrmodel,
                           loc = loc, model = "Gaussian",
                           distance = distance, radius = radius,
                           anisopars = anisopars,
                           param = c(sill = 1, corr_param))
  else stop("Local option cannot be used for gibbs sampling ...\n")
  
  # Pre-calcoli per ottimizzazione
  Sigma.mat.inv <- MatInv(mtx = GeoW$covmatrix)
  M0 <- nrow(coord_obs)
  N0 <- nrow(loc)
  Wglob <- GeoW$weights
  mm_loc <- mm
  data_transformed <- exp(-mm_loc) * data
  sqrt_factor <- sqrt(2 / shape)
  exp_mm_half <- exp(mm_loc) * 0.5

  # STRATEGIA MEMORY-EFFICIENT: Genera simulazioni in batch
  if (parallel && nrep > 1) {
    
    # Calcola memoria necessaria per un batch
    coord_sim <- rbind(coord_obs, loc)
    n_total <- nrow(coord_sim)
    
    # Stima: ogni simulazione ~ 8 bytes * n_total
    memory_per_sim_mb <- (8 * n_total) / (1024^2)
    memory_per_batch_mb <- memory_per_sim_mb * shape
    
    # Determina batch size: massimo 200 MB per batch
    max_batch_memory_mb <- 200
    sims_per_batch <- max(1, floor(max_batch_memory_mb / memory_per_batch_mb))
    
    #cat(sprintf("Memory-efficient mode: processing %d simulations per batch\n", sims_per_batch))
    #cat(sprintf("Total batches needed: %d\n", ceiling(nrep / sims_per_batch)))
    
    # Setup progress
    if (progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = seq_len(nrep))
    }
    
    # Parallelizzazione per batch
    res <- vector("list", nrep)
    n_batches <- ceiling(nrep / sims_per_batch)
    
    for (batch_idx in seq_len(n_batches)) {
      start_idx <- (batch_idx - 1) * sims_per_batch + 1
      end_idx <- min(batch_idx * sims_per_batch, nrep)
      batch_nrep <- end_idx - start_idx + 1
      
      cat(sprintf("\n=== Batch %d/%d: simulations %d to %d ===\n", 
                  batch_idx, n_batches, start_idx, end_idx))
      
      # Genera simulazioni Gaussiane SOLO per questo batch
      if (method == "TB" || method == "CE") {
        Gauss_batch <- GeoSimapprox(
          coordx = coord_sim, 
          corrmodel = corrmodel, 
          method = method, 
          model = "Gaussian", 
          param = c(corr_param, mean = 0, nugget = 0, sill = 1), 
          L = L, 
          nrep = shape * batch_nrep, 
          progress = FALSE, 
          parallel = FALSE
        )
      } else {
        Gauss_batch <- GeoSim(
          coordx = coord_sim, 
          corrmodel = corrmodel, 
          model = "Gaussian", 
          param = c(corr_param, mean = 0, nugget = 0, sill = 1), 
          nrep = shape * batch_nrep, 
          progress = FALSE
        )
      }
      
      pair_idx_batch <- split(seq_len(shape * batch_nrep), 
                             ceiling(seq_len(shape * batch_nrep) / shape))
      
      # Estrai solo i dati necessari
      Gauss_data_batch <- Gauss_batch$data
      
      # Variabili essenziali per questo batch
      essential_vars <- list(
        shape = shape,
        data_transformed = data_transformed,
        M0 = M0,
        N0 = N0,
        Sigma.mat.inv = Sigma.mat.inv,
        Wglob = Wglob,
        exp_mm_half = exp_mm_half,
        sqrt_factor = sqrt_factor,
        pair_idx = pair_idx_batch
      )
      
      # Worker function ottimizzata
      run_simulation_worker <- function(local_i, vars, Gauss_data) {
        shape <- vars$shape
        data_transformed <- vars$data_transformed
        M0 <- vars$M0
        N0 <- vars$N0
        Sigma.mat.inv <- vars$Sigma.mat.inv
        Wglob <- vars$Wglob
        exp_mm_half <- vars$exp_mm_half
        sqrt_factor <- vars$sqrt_factor
        pair_idx <- vars$pair_idx
        
        # Gibbs sampling
        random_sign <- matrix(ifelse(runif(M0 * shape) <= 0.5, 1, -1), 
                             nrow = M0, ncol = shape)
        ini_vals <- random_sign * matrix(rep(sqrt(data_transformed), shape) * sqrt_factor, 
                                         nrow = M0, ncol = shape)
        
        gibbs.result <- tryCatch({
          dotCall64::.C64("gamma_gibbs_sampler",
                          SIGNATURE = c("double", "double", "int", "int", 
                                       "double", "int", "int"),
                          SigmaInv = as.double(Sigma.mat.inv),
                          y        = as.double(data_transformed),
                          n        = as.integer(M0),
                          v        = as.integer(shape),
                          U        = as.double(ini_vals),
                          nIte     = as.integer(1000),
                          nRep     = as.integer(150),
                          INTENT   = c("r", "r", "r", "r", "rw", "r", "r"), 
                          VERBOSE = 0, 
                          NAOK = FALSE, 
                          PACKAGE = "GeoModels")
        }, error = function(e) {
          return(list(U = NULL))
        })
        
        if (is.null(gibbs.result$U)) return(NULL)
        
        gibbs_U <- matrix(gibbs.result$U, nrow = M0, ncol = shape)
        
        # Processamento
        idx_xy <- pair_idx[[local_i]]
        res_i <- matrix(0, nrow = N0, ncol = shape)
        
        for (kk in seq_len(shape)) {
          simX <- Gauss_data[[idx_xy[kk]]]
          sim_train_x <- simX[seq_len(M0)]
          sim_valid_x <- simX[(M0 + 1):(N0 + M0)]
          sk_pred_x <- as.numeric(crossprod(Wglob, gibbs_U[, kk] - sim_train_x))
          res_i[, kk] <- sim_valid_x + sk_pred_x
        }
        
        return(exp_mm_half * rowSums(res_i^2))
      }
      
      # Parallelizzazione del batch corrente
      old_limit <- options(future.globals.maxSize = 2000 * 1024^2)
      on.exit(options(old_limit), add = TRUE)
      
      future::plan(future::multisession, workers = ncores)
      
      batch_results <- future.apply::future_lapply(
        seq_len(batch_nrep),
        function(local_i) {
          result <- run_simulation_worker(local_i, essential_vars, Gauss_data_batch)
          return(result)
        },
        future.seed = TRUE,
        future.packages = c("dotCall64", "GeoModels"),
        future.globals = structure(TRUE, add = c("essential_vars", "Gauss_data_batch"))
      )
      
      # Salva risultati del batch
      for (local_i in seq_len(batch_nrep)) {
        res[[start_idx + local_i - 1]] <- batch_results[[local_i]]
        if (progress) {
          pb(sprintf("Simulation %d/%d completed", start_idx + local_i - 1, nrep))
        }
      }
      
      # Pulizia memoria dopo ogni batch
      rm(Gauss_batch, Gauss_data_batch, batch_results, essential_vars)
      gc()
      
      future::plan(future::sequential)
    }
    
  } else {
    # Modalità sequenziale (codice originale)
    if (method == "TB" || method == "CE") {
      Gauss.all <- GeoSimapprox(
        coordx = rbind(coord_obs, loc), 
        corrmodel = corrmodel, 
        method = method, 
        model = "Gaussian", 
        param = c(corr_param, mean = 0, nugget = 0, sill = 1), 
        L = L, 
        nrep = shape * nrep, 
        progress = FALSE, 
        parallel = FALSE
      )
    } else {
      Gauss.all <- GeoSim(
        coordx = rbind(coord_obs, loc), 
        corrmodel = corrmodel, 
        model = "Gaussian", 
        param = c(corr_param, mean = 0, nugget = 0, sill = 1), 
        nrep = shape * nrep, 
        progress = FALSE
      )
    }
    
    pair_idx <- split(seq_len(shape * nrep), ceiling(seq_len(shape * nrep) / shape))
    
    cat("Performing", nrep, "conditional simulations ...\n")
    if (progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = seq_len(nrep))
    }
    
    res <- vector("list", nrep)
    
    for (i in seq_len(nrep)) {
      random_sign <- matrix(ifelse(runif(M0 * shape) <= 0.5, 1, -1), 
                           nrow = M0, ncol = shape)
      ini_vals <- random_sign * matrix(rep(sqrt(data_transformed), shape) * sqrt_factor, 
                                       nrow = M0, ncol = shape)
      
      gibbs.result <- dotCall64::.C64("gamma_gibbs_sampler",
                                      SIGNATURE = c("double", "double", "int", "int", 
                                                   "double", "int", "int"),
                                      SigmaInv = as.double(Sigma.mat.inv),
                                      y        = as.double(data_transformed),
                                      n        = as.integer(M0),
                                      v        = as.integer(shape),
                                      U        = as.double(ini_vals),
                                      nIte     = as.integer(1000),
                                      nRep     = as.integer(150),
                                      INTENT   = c("r", "r", "r", "r", "rw", "r", "r"), 
                                      VERBOSE = 0, 
                                      NAOK = FALSE, 
                                      PACKAGE = "GeoModels")
      
      if (is.null(gibbs.result$U)) {
        warning(paste("Gibbs sampler failed at simulation", i))
        next
      }
      
      gibbs_U <- matrix(gibbs.result$U, nrow = M0, ncol = shape)
      
      idx_xy <- pair_idx[[i]]
      res_i <- matrix(0, nrow = N0, ncol = shape)
      
      for (kk in seq_len(shape)) {
        simX <- Gauss.all$data[[idx_xy[kk]]]
        sim_train_x <- simX[seq_len(M0)]
        sim_valid_x <- simX[(M0 + 1):(N0 + M0)]
        sk_pred_x <- as.numeric(crossprod(Wglob, gibbs_U[, kk] - sim_train_x))
        res_i[, kk] <- sim_valid_x + sk_pred_x
      }
      
      res[[i]] <- exp_mm_half * rowSums(res_i^2)
      
      if (progress) {
        pb(sprintf("Simulation %d/%d completed", i, nrep))
      }
    }
  }
  
  # Rimuove risultati nulli
  res <- res[!sapply(res, is.null)]
  failed_sims <- nrep - length(res)
  if (failed_sims > 0) {
    warning(paste(failed_sims, "simulazioni su", nrep, "sono fallite"))
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
                     coord_obs, loc, coordt_use, time, X, Xloc, Mloc, distance, radius, anisopars,
                     local, neighb, maxdist, maxtime,
                     space, spacetime, bivariate, parallel, ncores, progress) {
  
  #############################################
  # Validazione input
  #############################################
  if (is.null(loc) || nrow(loc) == 0) {
    stop("loc must be a non-empty matrix")
  }
  if (nrep < 1) {
    stop("nrep must be at least 1")
  }
  
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
                     X = X_sim, nrep = nrep, distance = distance, radius = radius, anisopars = anisopars)
    sim_args <- c(sim_args, list(param = param))
    sim_nc <- do.call(GeoSim, sim_args)
  }
  if (method == "TB" || method == "CE") {
    sim_args_approx <- list(coordx = coord_sim, coordt = time_sim, corrmodel = corrmodel, progress = FALSE,
                            method = method, L = L, parallel = parallel,
                            X = X_sim, nrep = nrep, distance = distance, radius = radius, anisopars = anisopars)
    sim_args_approx <- c(sim_args_approx, list(param = param))
    sim_nc <- do.call(GeoSimapprox, sim_args_approx)
  }

  #############################################
  # Kriging predictions
  #############################################
  krig_sim_args <- list(coordx = coord_obs, coordt = coordt_use,
                        data = data, corrmodel = corrmodel,
                        loc = loc, X = X, Xloc = Xloc, Mloc = Mloc, time = time,
                        distance = distance, radius = radius, anisopars = anisopars)
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
  if (!local) cat("Computing kriging weights ...\n")
  else cat("Computing local kriging weights ...\n")
  
  weights_args <- list(coordx = coord_obs, coordt = coordt_use,
                       corrmodel = corrmodel, loc = loc, 
                       X = X, Xloc = Xloc, Mloc = Mloc, time = time,
                       distance = distance, radius = radius, anisopars = anisopars)
  weights_args <- c(weights_args, list(param = param))

  if (local) {
    if (!is.null(neighb))   weights_args$neighb  <- neighb
    if (!is.null(maxdist))  weights_args$maxdist <- maxdist
    if (!is.null(maxtime))  weights_args$maxtime <- maxtime
    if (parallel) weights_args$parallel <- parallel
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
                       object.size(if (local) weights_list else W) + 
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
          on.exit(parallel::stopCluster(cl), add = TRUE)
          parallel::clusterExport(cl, c("compute_cond_sim"), envir = environment())

          chunk_sim_data_only <- sim_nc$data[chunk_indices]
          krig_pred_only <- krig_sim$pred
          weights_only <- if (local) weights_list else W

          chunk_results <- parallel::parLapply(cl, seq_along(chunk_indices), function(local_i) {
            res <- compute_cond_sim(local_i, chunk_sim_data_only, krig_pred_only, local,
                                    weights_only, n_obs, n_loc)
            res
          })
          parallel::stopCluster(cl)

          for (local_i in seq_along(chunk_indices)) {
            sim_cond[[chunk_indices[local_i]]] <- chunk_results[[local_i]]
            p()  # Progress update
          }
          rm(chunk_sim_data_only, krig_pred_only, weights_only, chunk_results)
          if (chunk_idx %% 5 == 0) gc()  # gc() solo ogni 5 chunk
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
        weights_obj <- if (local) weights_list else W

        sim_cond <- parallel::parLapply(cl, seq_len(nrep), function(i) {
          res <- compute_cond_sim(i, sim_data_list, krig_pred_vec, local,
                                  weights_obj, n_obs, n_loc)
          res
        })
        
        # Progress update dopo il completamento parallelo
        for (i in seq_len(nrep)) {
          p()
        }
      }
    } else {
      # Sequential
      sim_cond <- vector("list", nrep)
      for (i in seq_len(nrep)) {
        sim_cond[[i]] <- tryCatch({
          compute_cond_sim(i, sim_nc$data, krig_sim$pred, local,
                          if (local) weights_list else W, n_obs, n_loc)
        }, error = function(e) {
          warning(paste("Simulation", i, "failed:", e$message))
          return(NULL)
        })
        p()  # Progress update
      }
      
      # Rimuove risultati nulli
      sim_cond <- sim_cond[!sapply(sim_cond, is.null)]
      failed_sims <- nrep - length(sim_cond)
      if (failed_sims > 0) {
        warning(paste(failed_sims, "simulazioni su", nrep, "sono fallite"))
      }
    }
  } else {
    stop("Gauss_cd currently not implemented for spacetime models")
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
  # Do not anisotropize loc here.  Anisotropy is propagated to
  # GeoSim/GeoKrig/GeoKrigWeights/GeoKriglocWeights, which handle it consistently.
  loc_orig <- loc
  # Temporal coordinates to use
  coordt_use <- if (spacetime || bivariate) coordt else NULL

  # Helper for multiplicative models: eta(s) = log(mu(s)) = X(s)' beta.
  # If X/Xloc are supplied, they are used; otherwise a constant mean is assumed.
  # Mloc, when supplied, is interpreted on the eta = log(mu) scale.
  get_eta_obs_loc <- function(param, X, Xloc, Mloc, n_obs, n_loc) {
    sel <- substr(names(param), 1, 4) == "mean"
    beta <- as.numeric(unlist(param[sel]))
    if (length(beta) == 0L) stop("Mean parameter(s) are missing\n")

    if (!is.null(X)) {
      Xobs <- as.matrix(X)
      if (ncol(Xobs) != length(beta)) {
        stop("ncol(X) must match the number of mean parameters\n")
      }
      eta_obs <- as.numeric(Xobs %*% beta)
    } else {
      if (length(param$mean) == n_obs) eta_obs <- as.numeric(param$mean)
      else eta_obs <- rep(as.numeric(param$mean)[1], n_obs)
    }

    if (!is.null(Mloc)) {
      eta_loc <- as.numeric(Mloc)
    } else if (!is.null(Xloc)) {
      Xloc_mat <- as.matrix(Xloc)
      if (ncol(Xloc_mat) != length(beta)) {
        stop("ncol(Xloc) must match the number of mean parameters\n")
      }
      eta_loc <- as.numeric(Xloc_mat %*% beta)
    } else {
      eta_loc <- rep(as.numeric(param$mean)[1], n_loc)
    }

    list(obs = eta_obs, loc = eta_loc)
  }

############################################################################
######################### not copula models ################################
############################################################################
if (is.null(copula)) {
  if (model == "Gaussian") {
    res <- Gauss_cd(data, corrmodel, nrep, method, L,
                    param,
                    coord_obs, loc, coordt_use, time, X, Xloc, Mloc, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)
  }

  ################################################################################
  ############### monotone transformations of one Gaussian RF ####################
  ################################################################################
  if (model == "LogGaussian") {
    eta <- get_eta_obs_loc(param, X, Xloc, Mloc, length(data), nrow(loc))
    vv <- param$sill

    # Multiplicative parametrization:
    # Y(s) = exp(eta(s)) * exp(sqrt(vv) Z(s) - vv/2), hence E[Y(s)] = exp(eta(s)).
    datanorm <- (log(data) - eta$obs + vv/2)/sqrt(vv)   # gaussian scale

    param_gauss <- param
    param_gauss$mean <- 0
    param_gauss$sill <- 1

    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param_gauss,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)

    res <- lapply(res, function(r) exp(eta$loc) * exp(sqrt(vv) * r - vv/2)) # back-transformation
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
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)
    res <- lapply(res, function(r) mm + sqrt(vv) * r * exp(tail * r^2/2))
  }

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
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)
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

  if (model == "SinhAsinh") {
    mm <- param$mean
    vv <- param$sill
    skew <- param$skew
    tail <- param$tail
    datanorm <- sinh((asinh((data - mm)/sqrt(vv)) / tail) - skew)
    param$mean <- 0; param$sill <- 1; param$skew <- NULL; param$tail <- NULL
    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)
    res <- lapply(res, function(r) mm + sqrt(vv) * sinh(tail * (asinh(r) + skew)))
  }
  
  ################################################################################
  ############### non monotone transformations of more than one Gaussian RF ######
  ################################################################################
  if (model == "SkewGaussian") {    
    res <- SkewGaussianSimcond(coord_obs, loc, data, param, corrmodel, local = local, 
                               neighb = neighb, L = L, nrep = nrep, method = method,
                               distance = distance, radius = radius, anisopars = anisopars,
                               parallel = parallel, ncores = ncores, progress = progress)
  }

  if (model == "Gamma") { 
    res <- GammaSimcond(coord_obs, loc, data, param, corrmodel, local = local, 
                        neighb = neighb, L = L, nrep = nrep, method = method,
                        distance = distance, radius = radius, anisopars = anisopars,
                        parallel = parallel, ncores = ncores, progress = progress)
  }
  
  if (model == "Weibull") {
    mm <- param$mean
    param.weibull <- param
    param.weibull$shape <- 2
    param.weibull$mean <- 0
    ck <- gamma((param$shape + 1)/param$shape)
    data.t <- (exp(-mm) * data * ck)^param$shape
    res <- GammaSimcond(coord_obs, loc, data.t, param.weibull, corrmodel, local = local, 
                        neighb = neighb, L = L, nrep = nrep, method = method,
                        distance = distance, radius = radius, anisopars = anisopars,
                        parallel = parallel, ncores = ncores, progress = progress)
    res <- lapply(res, function(x) exp(mm) * (x^(1/param$shape))/ck)
  }
  
  if (model == "Poisson")  { stop("Poisson not implemented...") }
  if (model == "Binomial") { stop("Binomial not implemented...") }

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
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)
    res <- lapply(res, function(r) pnorm(r))
    res <- lapply(res, function(r) qweibull(r, shape = ss, scale = mm/(gamma(1 + 1/ss))))
  }

  if (model == "LogGaussian") {
    eta <- get_eta_obs_loc(param, X, Xloc, Mloc, length(data), nrow(loc))
    vv  <- param$sill
    eps <- .Machine$double.eps

    # Multiplicative parametrization with mean exp(eta):
    # log Y(s) ~ N(eta(s) - vv/2, vv).
    datanorm1 <- plnorm(data, meanlog = eta$obs - vv/2, sdlog = sqrt(vv), log.p = FALSE)
    datanorm1 <- pmin(pmax(datanorm1, eps), 1 - eps)
    datanorm  <- qnorm(datanorm1)

    param_gauss <- param
    param_gauss$mean <- 0
    param_gauss$sill <- 1

    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param_gauss,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)

    res <- lapply(res, function(r) pmin(pmax(pnorm(r), eps), 1 - eps))
    res <- lapply(res, function(u) qlnorm(u, meanlog = eta$loc - vv/2, sdlog = sqrt(vv), log.p = FALSE))
  }

  if (model == "Gamma") {
    mm <- exp(param$mean)
    ss <- param$shape
    datanorm1 <- pgamma(data, shape = ss/2, rate = ss/(2 * mm))
    datanorm  <- qnorm(datanorm1)
    param$mean <- 0; param$sill <- 1; param$shape <- NULL
    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)
    res <- lapply(res, function(r) pnorm(r))
    res <- lapply(res, function(r) qgamma(r, shape = ss/2, rate = ss/(2 * mm)))
  }

  if (model == "SkewLaplace") {
    pSkewLaplace <- function(x, sk, mm, vv) {
      s <- sqrt(vv)
      F <- numeric(length(x))
      sel1 <- x < mm
      F[sel1] <- sk * exp((1 - sk) * (x[sel1] - mm) / s)
      sel2 <- !sel1
      F[sel2] <- 1 - (1 - sk) * exp(-sk * (x[sel2] - mm) / s)
      return(F)
    }

    qtpSkewLaplace22 <- function(u, sk, mm, vv) {
      s  <- sqrt(vv)
      res <- rep(NA_real_, length(u))
      hi <- u >= sk
      lo <- u <  sk
      res[hi] <- mm - (s / sk) * (log1p(-u[hi]) - log1p(-sk))
      res[lo] <- mm + (s / (1 - sk)) * (log(u[lo]) - log(sk))
      res
    }
    
    mm <- param$mean
    skew <- param$skew
    vv <- param$sill
    datanorm1 <- pSkewLaplace(data, skew, mm, vv)
    datanorm  <- qnorm(datanorm1)
    param$mean <- 0; param$sill <- 1; param$skew <- NULL
    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)
    res <- lapply(res, function(r) pnorm(r))
    res <- lapply(res, function(r) qtpSkewLaplace22(r, skew, mm, vv))
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
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)
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
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)
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
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)
    res <- lapply(res, function(r) pnorm(r))
    res <- lapply(res, function(r) qnorm(r, mm, sqrt(vv)))
  }

  if (model == "Binomial") {
    eps <- .Machine$double.eps
    eta <- get_eta_obs_loc(param, X, Xloc, Mloc, length(data), nrow(loc))
    prob_obs <- pmin(pmax(pnorm(eta$obs), eps), 1 - eps)
    prob_loc <- pmin(pmax(pnorm(eta$loc), eps), 1 - eps)

    # Mid-PIT approximation for discrete observed data.
    datanorm1 <- (
      pbinom(data - 1, size = n, prob = prob_obs) +
      pbinom(data,     size = n, prob = prob_obs)
    ) / 2
    datanorm1 <- pmin(pmax(datanorm1, eps), 1 - eps)
    datanorm  <- qnorm(datanorm1)

    param_gauss <- param
    param_gauss$mean <- 0
    param_gauss$sill <- 1

    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param_gauss,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)
    res <- lapply(res, function(r) pmin(pmax(pnorm(r), eps), 1 - eps))
    res <- lapply(res, function(u) qbinom(p = u, size = n, prob = prob_loc))
  }

  if (model == "BinomialNeg") {
    eps <- .Machine$double.eps
    eta <- get_eta_obs_loc(param, X, Xloc, Mloc, length(data), nrow(loc))
    prob_obs <- pmin(pmax(pnorm(eta$obs), eps), 1 - eps)
    prob_loc <- pmin(pmax(pnorm(eta$loc), eps), 1 - eps)
    size <- n

    # Mid-PIT approximation for discrete observed data.
    datanorm1 <- (
      pnbinom(data - 1, size = size, prob = prob_obs) +
      pnbinom(data,     size = size, prob = prob_obs)
    ) / 2
    datanorm1 <- pmin(pmax(datanorm1, eps), 1 - eps)
    datanorm  <- qnorm(datanorm1)

    param_gauss <- param
    param_gauss$mean <- 0
    param_gauss$sill <- 1

    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param_gauss,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)
    res <- lapply(res, function(r) pmin(pmax(pnorm(r), eps), 1 - eps))
    res <- lapply(res, function(u) qnbinom(p = u, size = size, prob = prob_loc))
  }

  if (model == "Poisson") {
    eps <- .Machine$double.eps
    eta <- get_eta_obs_loc(param, X, Xloc, Mloc, length(data), nrow(loc))
    mu_obs <- exp(eta$obs)
    mu_loc <- exp(eta$loc)

    # Mid-PIT approximation for discrete observed data.
    datanorm1 <- (
      ppois(data - 1, lambda = mu_obs) +
      ppois(data,     lambda = mu_obs)
    ) / 2
    datanorm1 <- pmin(pmax(datanorm1, eps), 1 - eps)
    datanorm  <- qnorm(datanorm1)

    param_gauss <- param
    param_gauss$mean <- 0
    param_gauss$sill <- 1

    res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param_gauss,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius, anisopars,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores, progress)
    res <- lapply(res, function(r) pmin(pmax(pnorm(r), eps), 1 - eps))
    res <- lapply(res, function(u) qpois(u, lambda = mu_loc))
  }
}

####################### skewgaussian copula ####################################################
if (copula == "SkewGaussian") {

### reparametrization ## from (-1,1) to (-Inf,+Inf)
if (abs(param$nu) >= 1) stop("nu parameter must be between -1 and 1", call. = FALSE)
param$skew <- param$nu / sqrt(1 - param$nu * param$nu)


    if (model == "Weibull") {

    mm <- exp(param$mean)
    ss <- param$shape
    datanorm1 <- pweibull(data, shape = ss, scale = mm/(gamma(1 + 1/ss))) 
    param$mean <- 0; param$sill <- 1; param$shape <- NULL

    omega <- sqrt(1 + param$skew^2)
    alpha <- as.numeric(param$skew)
    datanorm  <- sn::qsn(datanorm1, xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))

    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param, corrmodel,local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, distance = distance, radius = radius, anisopars = anisopars, parallel = parallel, ncores = ncores,progress=progress)

    res=lapply(res, function(r) sn::psn(r ,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha)))
    res <- lapply(res, function(r) qweibull(r, shape = ss, scale = mm/(gamma(1 + 1/ss))))
  }

    if (model == "LogGaussian") {
    eta <- get_eta_obs_loc(param, X, Xloc, Mloc, length(data), nrow(loc))
    vv  <- param$sill
    eps <- .Machine$double.eps

    # Multiplicative parametrization with mean exp(eta):
    # log Y(s) ~ N(eta(s) - vv/2, vv).
    datanorm1 <- plnorm(data, meanlog = eta$obs - vv/2, sdlog = sqrt(vv), log.p = FALSE)
    datanorm1 <- pmin(pmax(datanorm1, eps), 1 - eps)

    param_skew <- param
    param_skew$mean <- 0
    param_skew$sill <- 1

    omega <- sqrt(1 + param_skew$skew^2)
    alpha <- as.numeric(param_skew$skew)
    datanorm  <- sn::qsn(datanorm1, xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))

    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param_skew, corrmodel,local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, distance = distance, radius = radius, anisopars = anisopars, parallel = parallel, ncores = ncores,progress=progress)
    res=lapply(res, function(r) pmin(pmax(sn::psn(r ,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha)), eps), 1 - eps))
    res <- lapply(res, function(u) qlnorm(u, meanlog = eta$loc - vv/2, sdlog = sqrt(vv), log.p = FALSE))
  }

  if (model == "Gamma") {

    mm <- exp(param$mean)
    ss <- param$shape
    datanorm1 <- pgamma(data, shape = ss/2, rate = ss/(2 * mm))
    param$mean <- 0; param$sill <- 1; param$shape <- NULL
    omega <- sqrt(1 + param$skew^2)
    alpha <- as.numeric(param$skew)
    datanorm  <- sn::qsn(datanorm1, xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))
    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param, corrmodel,local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, distance = distance, radius = radius, anisopars = anisopars, parallel = parallel, ncores = ncores,progress=progress)
    res=lapply(res, function(r) sn::psn(r ,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha)))
    res <- lapply(res, function(r) qgamma(r, shape = ss/2, rate = ss/(2 * mm)))
  }

    if (model == "Beta2") {

    mm <- 1/(1 + exp(-param$mean))
    ss <- param$shape
    pmin <- param$min
    pmax <- param$max
    datanorm1 <- pbeta((data - pmin)/(pmax - pmin), shape1 = mm * ss, shape2 = (1 - mm) * ss)
    param$mean <- 0; param$sill <- 1;
    omega <- sqrt(1 + param$skew^2)
    alpha <- as.numeric(param$skew)
    datanorm  <- sn::qsn(datanorm1, xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))
    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param, corrmodel, local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, distance = distance, radius = radius, anisopars = anisopars, parallel = parallel, ncores = ncores,progress=progress)
    res=lapply(res, function(r) sn::psn(r ,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha)))
    res <- lapply(res, function(r) pmin + (pmax - pmin) * qbeta(r, shape1 = mm * ss, shape2 = (1 - mm) * ss))
  }

    if (model == "StudentT") {
    mm <- param$mean
    vv <- param$sill
    df <- param$df
    datanorm1 <- pt((data - mm)/sqrt(vv), df = 1/df)
    param$mean <- 0; param$sill <- 1; param$df <- NULL
    omega <- sqrt(1 + param$skew^2)
    alpha <- as.numeric(param$skew)
    datanorm  <- sn::qsn(datanorm1, xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))
    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param, corrmodel, local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, distance = distance, radius = radius, anisopars = anisopars, parallel = parallel, ncores = ncores,progress=progress)
    res=lapply(res, function(r) sn::psn(r ,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha)))
    res <- lapply(res, function(r) mm + sqrt(vv) * qt(r, df = 1/df))
  }

    if (model == "Gaussian") {
    mm <- param$mean
    vv <- param$sill
    datanorm1 <- pnorm(data, mm, sqrt(vv))
    param$mean <- 0; param$sill <- 1
    omega <- sqrt(1 + param$skew^2)
    alpha <- as.numeric(param$skew)
    datanorm  <- sn::qsn(datanorm1, xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))
    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param, corrmodel, local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, distance = distance, radius = radius, anisopars = anisopars, parallel = parallel, ncores = ncores,progress=progress)
    res=lapply(res, function(r) sn::psn(r ,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha)))
    res <- lapply(res, function(r) qnorm(r, mm, sqrt(vv)))
  }

    if (model == "Binomial") { 
    eps <- .Machine$double.eps
    eta <- get_eta_obs_loc(param, X, Xloc, Mloc, length(data), nrow(loc))
    prob_obs <- pmin(pmax(pnorm(eta$obs), eps), 1 - eps)
    prob_loc <- pmin(pmax(pnorm(eta$loc), eps), 1 - eps)

    # Mid-PIT approximation for discrete observed data.
    datanorm1 <- (
      pbinom(data - 1, size = n, prob = prob_obs) +
      pbinom(data,     size = n, prob = prob_obs)
    ) / 2
    datanorm1 <- pmin(pmax(datanorm1, eps), 1 - eps)

    param_skew <- param
    param_skew$mean <- 0
    param_skew$sill <- 1
    omega <- sqrt(1 + param_skew$skew^2)
    alpha <- as.numeric(param_skew$skew)

    datanorm <- sn::qsn(datanorm1, xi = 0, omega = as.numeric(omega), alpha = as.numeric(alpha))
    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param_skew, corrmodel, local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, distance = distance, radius = radius, anisopars = anisopars, parallel = parallel, ncores = ncores,progress=progress)
    res <- lapply(res, function(r) pmin(pmax(sn::psn(r, xi = 0, omega = as.numeric(omega), alpha = as.numeric(alpha)), eps), 1 - eps))
    res <- lapply(res, function(u) qbinom(p = u, size = n, prob = prob_loc))
  }

  if (model == "BinomialNeg") {
    eps <- .Machine$double.eps
    eta <- get_eta_obs_loc(param, X, Xloc, Mloc, length(data), nrow(loc))
    prob_obs <- pmin(pmax(pnorm(eta$obs), eps), 1 - eps)
    prob_loc <- pmin(pmax(pnorm(eta$loc), eps), 1 - eps)
    size <- n

    # Mid-PIT approximation for discrete observed data.
    datanorm1 <- (
      pnbinom(data - 1, size = size, prob = prob_obs) +
      pnbinom(data,     size = size, prob = prob_obs)
    ) / 2
    datanorm1 <- pmin(pmax(datanorm1, eps), 1 - eps)

    param_skew <- param
    param_skew$mean <- 0
    param_skew$sill <- 1
    omega <- sqrt(1 + param_skew$skew^2)
    alpha <- as.numeric(param_skew$skew)

    datanorm <- sn::qsn(datanorm1, xi = 0, omega = as.numeric(omega), alpha = as.numeric(alpha))
    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param_skew, corrmodel, local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, distance = distance, radius = radius, anisopars = anisopars, parallel = parallel, ncores = ncores,progress=progress)
    res <- lapply(res, function(r) pmin(pmax(sn::psn(r, xi = 0, omega = as.numeric(omega), alpha = as.numeric(alpha)), eps), 1 - eps))
    res <- lapply(res, function(u) qnbinom(p = u, size = size, prob = prob_loc))
  }

  if (model == "Poisson") {
    eps <- .Machine$double.eps
    eta <- get_eta_obs_loc(param, X, Xloc, Mloc, length(data), nrow(loc))
    mu_obs <- exp(eta$obs)
    mu_loc <- exp(eta$loc)

    # Mid-PIT approximation for discrete observed data.
    datanorm1 <- (
      ppois(data - 1, lambda = mu_obs) +
      ppois(data,     lambda = mu_obs)
    ) / 2
    datanorm1 <- pmin(pmax(datanorm1, eps), 1 - eps)

    param_skew <- param
    param_skew$mean <- 0
    param_skew$sill <- 1
    omega <- sqrt(1 + param_skew$skew^2)
    alpha <- as.numeric(param_skew$skew)

    datanorm <- sn::qsn(datanorm1, xi = 0, omega = as.numeric(omega), alpha = as.numeric(alpha))
    res <- SkewGaussianSimcond(coord_obs, loc, datanorm, param_skew, corrmodel, local=local,neighb=neighb,L=L,
                            nrep = nrep, method = method, distance = distance, radius = radius, anisopars = anisopars, parallel = parallel, ncores = ncores,progress=progress)
    res <- lapply(res, function(r) pmin(pmax(sn::psn(r, xi = 0, omega = as.numeric(omega), alpha = as.numeric(alpha)), eps), 1 - eps))
    res <- lapply(res, function(u) qpois(u, lambda = mu_loc))
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