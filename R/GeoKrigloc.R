####################################################
### File name: GeoKrigloc_optimized.r
### Optimized version with reduced code duplication
### and improved efficiency
####################################################

GeoKrigloc = function(estobj=NULL, data, coordx, coordy=NULL, coordz=NULL, coordt=NULL, coordx_dyn=NULL,
                      corrmodel, distance="Eucl", grid=FALSE, loc, neighb=NULL,
                      maxdist=NULL, maxtime=NULL, method="cholesky", model="Gaussian",
                      n=1, nloc=NULL, mse=FALSE, param, anisopars=NULL,
                      radius=1, sparse=FALSE, time=NULL, type="Standard",
                      type_krig="Simple", weigthed=TRUE, which=1,
                      copula=NULL, X=NULL, Xloc=NULL, Mloc=NULL, varcov=NULL,
                      spobj=NULL, spdata=NULL, parallel=FALSE, ncores=NULL, progress=TRUE)
{
  ###########################################################################
  ##  0. Helper functions (optimized)                                      ##
  ###########################################################################
  
  setup_parallel <- function(ncores) {
    n.cores <- if (is.null(ncores)) {
      max(1, min(4, parallel::detectCores() - 1))
    } else {
      if (!is.numeric(ncores) || ncores < 1 || ncores > parallel::detectCores())
        stop("Invalid number of cores specified")
      as.integer(ncores)
    }
    
    old_opts <- options(
      warn = -1,
      verbose = FALSE,
      future.globals.maxSize = 4 * 1024^3
    )
    
    future::plan(future::multisession, workers = n.cores)
    
    list(n.cores = n.cores, old_opts = old_opts)
  }

  cleanup_parallel <- function(cfg) {
    options(cfg$old_opts)
    future::plan(future::sequential)
    invisible(gc(verbose = FALSE))
  }

  estimate_task_ram <- function(neigh, k) {
    n <- nrow(neigh$coordx[[k]])
    p <- if(!is.null(neigh$X[[k]])) ncol(neigh$X[[k]]) else 0L
    # Estimate: covariance matrix + X matrix + overhead
    8 * (n * n + n * p) * 3
  }

  choose_workers_safe <- function(neigh, total_tasks, max_cores) {
    free_ram <- tryCatch({
      if(requireNamespace("memuse", quietly = TRUE)) {
        mu <- memuse::Sys.meminfo()
        as.numeric(mu["freeram"], units = "b")
      } else {
        2 * 1024^3  # Default: 2GB
      }
    }, error = function(e) 2 * 1024^3)

    # Sample up to 100 tasks to estimate worst-case RAM
    sample_size <- min(100, total_tasks)
    worst_ram <- max(vapply(seq_len(sample_size), 
                            function(k) estimate_task_ram(neigh, k),
                            numeric(1)))
    
    needed_per_task <- worst_ram * 1.5  # Safety margin
    safe_cores <- max(1L, as.integer(floor(free_ram / needed_per_task)))
    
    min(max_cores, safe_cores)
  }

  ###########################################################################
  ##  1. Input processing (unchanged)                                      ##
  ###########################################################################
  
  call = match.call()

  if(!is.null(estobj)){
    if(!inherits(estobj,"GeoFit")) stop("need a 'GeoFit' object as input\n")
    data = estobj$data
    if(!estobj$grid){
      if(!estobj$bivariate){
        if(is.null(estobj$coordx_dyn))
          coordx = cbind(estobj$coordx, estobj$coordy, estobj$coordz)
        else
          coordx_dyn = estobj$coordx_dyn
      } else {
        if(is.null(estobj$coordx_dyn)){
          coordx = estobj$coordx[1:estobj$ns[1]]
          coordy = estobj$coordy[1:estobj$ns[2]]
        } else {
          coordx_dyn = estobj$coordx_dyn
        }
      }
    } else {
      coordx = estobj$coordx
      coordy = estobj$coordy
      coordz = estobj$coordz
    }
    if(length(estobj$coordt) == 1) coordt = NULL else coordt = estobj$coordt
    coordx_dyn = estobj$coordx_dyn
    corrmodel  = estobj$corrmodel
    model      = estobj$model
    distance   = estobj$distance
    grid       = estobj$grid
    n          = estobj$n
    param      = append(estobj$param, estobj$fixed)
    radius     = estobj$radius
    copula     = estobj$copula
    anisopars  = estobj$anisopars
    X          = estobj$X
    varcov     = estobj$varcov
  }

  bivariate  = CheckBiv(CkCorrModel(corrmodel))
  spacetime  = CheckST(CkCorrModel(corrmodel))
  space      = !spacetime && !bivariate

  if(!is.null(spobj)) {
    a = sp2Geo(spobj, spdata); coordx = a$coords
    if(!a$pj && distance != "Chor") distance = "Geod"
    if(spacetime) coordt = a$coordt
    if(!is.null(a$Y) && !is.null(a$X)) { data = a$Y; X = a$X }
  }

  if(is.null(coordx_dyn)){
    coords = coordx
    if(!is.null(coordy)){
      if(!grid) coords = cbind(coordx, coordy, coordz)
      else {
        if(!is.null(coordz))
          coords = as.matrix(expand.grid(coordx, coordy, coordz))
        else
          coords = as.matrix(expand.grid(coordx, coordy))
      }
    }
  } else {
    coordx = coordy = coordz = coords = NULL
  }

  loc  = as.matrix(loc)
  Nloc = nrow(loc)
  if(is.null(Nloc)) Nloc = 1
  Tloc = length(time)
  if(bivariate) Tloc = 1
  if(length(param$mean) > 1) M = param$mean else M = NULL

  ###########################################################################
  ##  2. Unified kriging function                                          ##
  ###########################################################################
  
  # Pre-compute task indices for spacetime/bivariate
  task_indices <- if(spacetime || bivariate) {
    # Pre-allocate matrix: each row is (i, j) for task k
    idx <- matrix(0L, nrow = Nloc * Tloc, ncol = 2L)
    k <- 1L
    for(i in seq_len(Nloc)) {
      for(j in seq_len(Tloc)) {
        idx[k, 1L] <- i
        idx[k, 2L] <- j
        k <- k + 1L
      }
    }
    idx
  } else {
    NULL
  }
  
  # Unified kriging function for all cases
  perform_kriging <- function(k, neigh, param_base, is_space) {
    # Create local parameter copy
    local_param <- param_base
    if(!is.null(M)) local_param$mean <- neigh$M[[k]]
    
    if(is_space) {
      # Space only
      pr <- GeoKrig(
        loc = loc[k, ], 
        data = neigh$data[[k]],
        coordx = neigh$coordx[[k]], 
        corrmodel = corrmodel,
        distance = distance, 
        n = n, 
        X = neigh$X[[k]],
        Xloc = Xloc[k, ], 
        Mloc = Mloc[k], 
        type_krig = type_krig,
        sparse = sparse, 
        model = model, 
        param = local_param,
        anisopars = anisopars, 
        radius = radius, 
        mse = mse,
        copula = copula, 
        varcov = varcov
      )
    } else {
      # Spacetime or bivariate
      i <- task_indices[k, 1L]
      j <- task_indices[k, 2L]
      
      if(spacetime) {
        pr <- GeoKrig(
          data = neigh$data[[k]], 
          coordx = neigh$coordx[[k]],
          coordt = neigh$coordt[[k]], 
          loc = loc[i, ],
          time = time[j], 
          X = neigh$X[[k]],
          Mloc = Mloc[k],
          Xloc = Xloc[k, ],
          type_krig = type_krig,
          sparse = sparse, 
          corrmodel = corrmodel,
          distance = distance, 
          model = model, 
          param = local_param,
          anisopars = anisopars, 
          radius = radius, 
          mse = mse,
          varcov = varcov,
          copula = copula, 
          n = n
        )
      } else {  # bivariate
        pr <- GeoKrig(
          data = neigh$data[[k]], 
          coordx = neigh$coordx[[k]],
          loc = loc[i, ], 
          X = neigh$X[[k]],
          Mloc = Mloc[k],
          Xloc = Xloc[k, ],
          type_krig = type_krig,
          sparse = sparse, 
          corrmodel = corrmodel,
          distance = distance, 
          model = model, 
          param = local_param,
          varcov = varcov,
          anisopars = anisopars, 
          radius = radius, 
          mse = mse,
          copula = copula, 
          n = n
        )
      }
    }
    
    list(pred = pr$pred, mse = if(mse) pr$mse else NULL)
  }

  ###########################################################################
  ##  3. Main kriging execution                                            ##
  ###########################################################################
  
  # Setup neighborhood
  if(space) {
    neigh = GeoNeighborhood(data, coordx = coords, distance = distance, loc = loc,
                            neighb = neighb, maxdist = maxdist, X = X, M = M,
                            parallel = FALSE, ncores = ncores)
    total_tasks <- Nloc
  } else if(spacetime) {
    neigh = GeoNeighborhood(data, coordx = coords, coordt = coordt,
                            distance = distance, neighb = neighb,
                            coordx_dyn = coordx_dyn, loc = loc, time = time,
                            maxdist = maxdist, maxtime = maxtime, X = X, M = M,
                            parallel = FALSE, ncores = ncores)
    total_tasks <- Nloc * Tloc
  } else {  # bivariate
    neigh = GeoNeighborhood(data, coordx = coords, distance = distance,
                            neighb = neighb, coordx_dyn = coordx_dyn,
                            loc = loc, maxdist = maxdist, X = X, M = M,
                            parallel = FALSE, ncores = ncores)
    total_tasks <- Nloc * Tloc
  }

  # Pre-allocate results
  res1 <- numeric(total_tasks)
  res2 <- if(mse) numeric(total_tasks) else NULL
  
  # Setup progress handlers
  if (!is.logical(progress)) stop("progress must be logical (TRUE/FALSE)")
  
  if(progress) {
    progressr::handlers(global = TRUE)
    progressr::handlers("txtprogressbar")
    pb <- progressr::progressor(along = seq_len(total_tasks))
  } else {
    pb <- function(...) {}
  }

  # Execute kriging (sequential or parallel)
  if(!parallel) {
    # Sequential execution
    for(k in seq_len(total_tasks)) {
      pb(sprintf("k=%d", k))
      pr <- perform_kriging(k, neigh, param, space)
      res1[k] <- pr$pred
      if(mse) res2[k] <- pr$mse
    }
    
  } else {
    # Parallel execution
    config <- setup_parallel(ncores)
    on.exit(cleanup_parallel(config), add = TRUE)
    
    max_cores <- config$n.cores
    safe_cores <- choose_workers_safe(neigh, total_tasks, max_cores)
    
    if(safe_cores < max_cores) {
      warning("Reducing cores from ", max_cores, " to ", safe_cores,
              " to avoid memory pressure.")
      future::plan(future::multisession, workers = safe_cores)
    }
    
    cat("Performing local kriging using ", safe_cores, " cores...\n")
    
    # Execute with progress inside workers
    res <- future.apply::future_lapply(
      seq_len(total_tasks),
      function(k) {
        suppressPackageStartupMessages({
          result <- perform_kriging(k, neigh, param, space)
        })
        if(progress) pb(sprintf("k=%d", k))
        result
      },
      future.seed = TRUE,
      future.stdout = FALSE,
      future.conditions = "none"
    )
    
    # Extract results efficiently
    res1 <- vapply(res, `[[`, numeric(1), "pred")
    if(mse) res2 <- vapply(res, `[[`, numeric(1), "mse")
  }

  ###########################################################################
  ##  4. Format output                                                     ##
  ###########################################################################
  
###########################################################################
##  4. Format output                                                     ##
###########################################################################

varpred <- NULL

if (spacetime || bivariate) {
  pred <- matrix(res1, nrow = Tloc, ncol = Nloc)   # byrow = FALSE (default)

  if (mse) varpred <- matrix(res2, nrow = Tloc, ncol = Nloc)
} else {
  pred <- res1
  if (mse) varpred <- res2
}

if (Tloc == 1) {
  pred <- c(pred)
  varpred <- c(varpred)
}

  GeoKrigloc_out <- list(
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
    numcoord = if(exists("coords", inherits = FALSE)) nrow(coords) else NA,
    numloc = Nloc,
    numtime = length(coordt),
    numt = Tloc,
    maxdist = maxdist,
    maxtime = maxtime,
    model = model,
    n = n,
    param = param,
    pred = pred,
    radius = radius,
    spacetime = spacetime,
    time = time,
    type_krig = type_krig,
    mse = varpred,
    varcov = varcov
  )
  
  structure(c(GeoKrigloc_out, call = call), class = c("GeoKrigloc"))
}