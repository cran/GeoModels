GeoNeighbSelect <- function(data, coordx, coordy=NULL, coordz=NULL, coordt=NULL, coordx_dyn=NULL,
  copula=NULL, corrmodel=NULL, distance="Eucl", fixed=NULL, anisopars=NULL,
  est.aniso=c(FALSE,FALSE), grid=FALSE, likelihood='Marginal', 
  lower=NULL, neighb=c(1,2,3,4,5),
  maxtime=Inf, memdist=TRUE, model='Gaussian', n=1, ncores=NULL,
  optimizer='Nelder-Mead', parallel=FALSE, bivariate=FALSE, radius=1, start=NULL, 
  type='Pairwise', upper=NULL, weighted=FALSE,
  X=NULL, nosym=FALSE, spobj=NULL, spdata=NULL, vario=NULL, progress=TRUE)
{
  # Input validation for neighb parameter
  if(!is.numeric(neighb)) stop("neighb must be a numeric vector")
  if(sum(neighb - floor(neighb))) stop("neighb must be a positive integer numeric vector")

  estimates = best_T = best_K = NULL

  # Validate correlation model
  corrmodel_valid <- try(CkCorrModel(corrmodel), silent = TRUE)
  if(inherits(corrmodel_valid, "try-error")) {
    stop("correlation model is not valid\n")
  }
  
  bivariate <- CheckBiv(corrmodel_valid)
  spacetime <- CheckST(corrmodel_valid)
  space <- !spacetime && !bivariate

  # Validate variogram object
  if(!is.null(vario)) {
    if(!inherits(vario, "GeoVariogram")) {
      stop("A GeoVariogram object is needed as input for vario\n")
    }
    if(!((vario$bivariate && bivariate) || (!is.null(vario$bint) && spacetime) || 
         (is.null(vario$bint) && !bivariate))) {
      stop("The GeoVariogram object is not of the same type of the correlation model\n")
    }
    semiv <- vario
  } else {
    stop("A GeoVariogram object is needed as input for vario\n")
  }

  # Validate start parameter
  if(is.null(start) || length(start) == 0) {
    stop("start parameter must be provided\n")
  }

  # Pre-compute common values
  coremax <- parallel::detectCores()
  if(is.na(coremax) || coremax == 1) parallel <- FALSE
  K <- length(neighb)
  P <- NULL
  
  # Optimized pre-allocation
  if(spacetime) {
    if(!is.numeric(maxtime)) stop("maxtime must be a numeric vector")
    P <- length(maxtime)
    estimates <- matrix(NA_real_, nrow = K * P, ncol = length(start))
    res <- rep(NA_real_, K * P)
  } else {
    estimates <- matrix(NA_real_, nrow = K, ncol = length(start))
    res <- rep(NA_real_, K)
  }

  # Pre-compute common fixed parameters
  common_params <- list(
    data = data, coordx = coordx, coordy = coordy, coordz = coordz, coordt = coordt, 
    coordx_dyn = coordx_dyn, copula = copula, corrmodel = corrmodel, distance = distance,
    fixed = fixed, anisopars = anisopars, est.aniso = est.aniso, grid = grid, 
    likelihood = likelihood, lower = lower, memdist = memdist, model = model, n = n, 
    optimizer = optimizer, radius = radius, start = start, type = type, upper = upper, 
    weighted = weighted, X = X, nosym = nosym, spobj = spobj, spdata = spdata
  )

  # Core fitting function - computes fit for given neighb and maxtime values
  compute_fit <- function(neighb_val, maxtime_val = Inf) {
    params <- common_params
    params$neighb <- neighb_val
    params$maxtime <- maxtime_val
    
    # Fit the model
    aa <- suppressWarnings(do.call(GeoFit, params))
    estimates_val <- unlist(aa$param)
    
    # Compute residuals if fit is valid
    if(!(aa$logCompLik == 0 || aa$logCompLik > 1e+13)) {
      cc <- suppressWarnings(
        GeoCorrFct(semiv$centers, t = semiv$centert, corrmodel = corrmodel, 
                   model = model, distance = distance, 
                   param = c(aa$param, aa$fixed), radius = radius, n = n, 
                   covariance = TRUE, variogram = TRUE)$corr
      )
      res_val <- sum((cc - semiv$variograms)^2)
    } else {
      res_val <- Inf
    }
    
    list(estimates = estimates_val, res = res_val)
  }

  #################### SPATIAL CASE ####################
  if(space || bivariate) {
    if(!parallel) {
      # Sequential processing with progress bar
      if(progress) {
        progressr::handlers(global = TRUE)
        progressr::handlers("txtprogressbar")
        pb <- progressr::progressor(along = 1:K)
      }
      
      for(M in seq_len(K)) {
        if(progress) pb(sprintf("k=%g", M)) 
        
        result <- compute_fit(neighb[M])
        estimates[M, ] <- result$estimates
        res[M] <- result$res
      }
    } else {
      # Optimized parallel configuration
      n.cores <- if(is.null(ncores)) {
        min(coremax - 1L, K)  # Don't use more cores than necessary
      } else {
        if(!is.numeric(ncores) || ncores > coremax || ncores < 1) {
          stop("number of cores not valid\n")
        }
        min(ncores, K)
      }
      
      cat("Performing", K, "estimations using", n.cores, "cores...\n")
      future::plan(future::multisession, workers = n.cores)
      
      if(progress) {
        progressr::handlers(global = TRUE)
        progressr::handlers("txtprogressbar")
        pb <- progressr::progressor(along = 1:K)
      }
      
      # Single fit function for parallel execution
      single_fit <- function(M) {
        result <- compute_fit(neighb[M])
        if(progress) pb(sprintf("k=%g", M))
        c(result$estimates, res = result$res)
      }
      
      # Execute parallel computations
      results <- future.apply::future_lapply(
        seq_len(K), single_fit, 
        future.seed = TRUE, 
        future.stdout = FALSE,
        future.conditions = "none"
      )
      
      # Optimized results conversion
      n_params <- length(start)
      estimates_list <- lapply(results, function(x) x[seq_len(n_params)])
      estimates <- do.call(rbind, estimates_list)
      res <- vapply(results, function(x) x[n_params + 1L], numeric(1))
      
      future::plan(sequential)
    }
  }

  #################### SPATIO-TEMPORAL CASE ####################
  if(spacetime) {
    total_iter <- K * P
    
    if(!parallel) {
      # Sequential processing with progress bar
      if(progress) {
        progressr::handlers(global = TRUE)
        progressr::handlers("txtprogressbar")
        pb <- progressr::progressor(along = 1:total_iter)
      }
      
      idx <- 1L
      for(L in seq_len(P)) { 
        for(M in seq_len(K)) {
          if(progress) pb(sprintf("k=%g", idx)) 
          
          result <- compute_fit(neighb[M], maxtime[L])
          estimates[idx, ] <- result$estimates
          res[idx] <- result$res
          idx <- idx + 1L
        }
      }
    } else {
      # Optimized parallel configuration
      n.cores <- if(is.null(ncores)) {
        min(coremax - 1L, total_iter)
      } else {
        if(!is.numeric(ncores) || ncores > coremax || ncores < 1) {
          stop("number of cores not valid\n")
        }
        min(ncores, total_iter)
      }
      
      cat("Performing", total_iter, "estimations using", n.cores, "cores...\n")
      future::plan(future::multisession, workers = n.cores)
      
      if(progress) {
        progressr::handlers(global = TRUE)
        progressr::handlers("txtprogressbar")
        pb <- progressr::progressor(along = 1:total_iter)
      }
      
      # Pre-compute combinations for parallel execution
      combinations <- expand.grid(L = seq_len(P), M = seq_len(K))
      
      # Optimized spatio-temporal parallel function
      single_fit_st <- function(i) {
        L <- combinations$L[i]
        M <- combinations$M[i]
        
        result <- compute_fit(neighb[M], maxtime[L])
        if(progress) pb(sprintf("k=%g", i))
        c(result$estimates, res = result$res)
      }
      
      # Execute parallel computations
      results <- future.apply::future_lapply(
        seq_len(nrow(combinations)), single_fit_st, 
        future.seed = TRUE, 
        future.stdout = FALSE,
        future.conditions = "none"
      )
      
      # Optimized results conversion
      n_params <- length(start)
      estimates_list <- lapply(results, function(x) x[seq_len(n_params)])
      estimates <- do.call(rbind, estimates_list)
      res <- vapply(results, function(x) x[n_params + 1L], numeric(1))
      
      future::plan(sequential)
    }
  }

  # Optimized final results computation
  indexmin <- which.min(res)
  
  if(space || bivariate) {
    bestK <- neighb[indexmin]
  }
  
  if(spacetime) {
    grid_result <- expand.grid(neighb, maxtime)
    bestK <- grid_result[indexmin, 1L]
    best_T <- grid_result[indexmin, 2L]
  }

  # Compute optimized suggestions
  sugg_neighb <- if(bestK <= 20) {
    bestK
  } else {
    sorted_neighb <- sort(neighb)
    median_neighb <- as.integer(quantile(sorted_neighb, 0.5))
    min(median_neighb, bestK)
  }

  sugg_time <- if(spacetime) {
    if(best_T <= 3) {
      best_T
    } else {
      sorted_maxtime <- sort(maxtime)
      median_time <- as.integer(quantile(sorted_maxtime, 0.5))
      min(median_time, best_T)
    }
  } else NULL

  # Return results list
  list(
    best_neighb = as.numeric(bestK),
    best_maxtime = best_T,
    res = unname(res),
    estimates = estimates,
    sugg_neighb = as.numeric(sugg_neighb),
    sugg_time = sugg_time
  )
}