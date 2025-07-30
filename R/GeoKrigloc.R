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
  ##  0.  helpers (memory & globals)                                       ##
  ###########################################################################
  # quiet parallel setup
  setup_parallel <- function(ncores, total_tasks) {
    if (is.null(ncores)){
      n.cores <- max(1, min(4, parallel::detectCores() - 1))
    } else {
      if (!is.numeric(ncores) || ncores < 1 || ncores > parallel::detectCores())
        stop("Invalid number of cores specified")
      n.cores <- ncores
    }
    old_warn  <- options(warn = -1)
    old_verb  <- options(verbose = FALSE)
    future::plan(future::multisession, workers = n.cores)
    oopts <- options(future.globals.maxSize = 4 * 1024^3)  # raise cap
    list(n.cores = n.cores, oopts = oopts,
         old_warn = old_warn, old_verb = old_verb)
  }
  cleanup_parallel <- function(cfg) {
    tryCatch({
      options(cfg$oopts)
      options(cfg$old_warn)
      options(cfg$old_verb)
      future::plan(future::sequential)
      invisible(gc(verbose = FALSE))
    }, error = function(e) warning("Cleanup error: ", e$message))
  }

  # crude RAM estimator
  estimate_task_ram <- function(neigh, k) {
    tryCatch({
      n <- nrow(neigh$coordx[[k]])
      p <- if(!is.null(neigh$X[[k]])) ncol(neigh$X[[k]]) else 0L
      8 * (n * n + n * p) * 3   # 3 copies, double
    }, error = function(e) 256 * 1024^2)   # fallback
  }
  choose_workers_safe <- function(neigh, total_tasks, max_cores) {
    free <- tryCatch({
      if(requireNamespace("memuse", quietly = TRUE)) {
        mu <- memuse::Sys.meminfo()
        as.numeric(mu["freeram"], units = "b")
      } else 2 * 1024^3
    }, error = function(e) 2 * 1024^3)
    worst  <- max(sapply(seq_len(min(100, total_tasks)),
                         function(k) estimate_task_ram(neigh, k)))
    needed <- worst * 1.5
    safe   <- max(1, floor(free / needed))
    min(max_cores, safe)
  }

  ###########################################################################
  ##  1.  original body (unchanged)                                        ##
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
    corrmodel = estobj$corrmodel
    model = estobj$model
    distance = estobj$distance
    grid = estobj$grid
    n = estobj$n
    param = append(estobj$param, estobj$fixed)
    radius = estobj$radius
    copula = estobj$copula
    anisopars = estobj$anisopars
    X = estobj$X
    varcov=estobj$varcov
  }

  bivariate = CheckBiv(CkCorrModel(corrmodel))
  spacetime = CheckST(CkCorrModel(corrmodel))
  space = !spacetime && !bivariate

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

  loc = as.matrix(loc)
  Nloc = nrow(loc)
  if(is.null(Nloc)) Nloc = 1
  Tloc = length(time)
  if(bivariate) Tloc = 1
  if(length(param$mean) > 1) M = param$mean else M = NULL

  ###########################################################################
  ##  2.  kriging loops – memory-safe & globals-light                      ##
  ###########################################################################
  ## -------- space --------------------------------------------------------
  if(space){
    neigh = GeoNeighborhood(data, coordx = coords, distance = distance, loc = loc,
                            neighb = neighb, maxdist = maxdist, X = X, M = M,
                            parallel = FALSE, ncores = ncores)

    res1 = numeric(Nloc)
    res2 = if(mse) numeric(Nloc) else NULL
    total_tasks = Nloc

    if(progress){
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = seq_len(total_tasks))
    } else pb <- function(...) {}

    if(!parallel){
      for(i in seq_len(total_tasks)){
        pb(sprintf("i=%d", i))
        if(!is.null(M)) param$mean = neigh$M[[i]]
        pr = GeoKrig(loc = loc[i, ], data = neigh$data[[i]],
                     coordx = neigh$coordx[[i]], corrmodel = corrmodel,
                     distance = distance, n = n, X = neigh$X[[i]],
                     Xloc = Xloc[i, ], Mloc = Mloc[i], type_krig = type_krig,
                     sparse = sparse, model = model, param = param,
                     anisopars = anisopars, radius = radius, mse = mse,
                     copula = copula,varcov=varcov)
        res1[i] = pr$pred
        if(mse) res2[i] = pr$mse
      }
    } else {
      config = setup_parallel(ncores, total_tasks)
      max_cores = config$n.cores
      safe_cores = choose_workers_safe(neigh, total_tasks, max_cores)
      if(safe_cores < max_cores){
        warning("Reducing cores from ", max_cores, " to ", safe_cores,
                " to avoid memory pressure.")
        future::plan(future::multisession, workers = safe_cores)
      }
      cat("Performing local kriging using ", safe_cores, " cores...\n")


res = future.apply::future_lapply(seq_len(total_tasks), function(i){
  suppressPackageStartupMessages({
    local_param = param
    if(!is.null(M)) local_param$mean = neigh$M[[i]]
    pr = GeoKrig(loc = loc[i, ], data = neigh$data[[i]],
                 coordx = neigh$coordx[[i]], corrmodel = corrmodel,
                 distance = distance, n = n, X = neigh$X[[i]],
                 Xloc = Xloc[i, ], Mloc = Mloc[i],
                 type_krig = type_krig, sparse = sparse,
                 model = model, param = local_param,
                 anisopars = anisopars, radius = radius,
                 mse = mse, copula = copula,varcov)
  })
  if(progress) pb(sprintf("i=%d", i))
  list(pred = pr$pred, mse = if(mse) pr$mse else NULL)
}, future.seed = TRUE, future.stdout = FALSE,
        future.conditions = "none")

      res1 = sapply(res, `[[`, "pred")
      if(mse) res2 = sapply(res, `[[`, "mse")

      cleanup_parallel(config)
    }
  }

  ## -------- spacetime -----------------------------------------------------
  if(spacetime){
    neigh = GeoNeighborhood(data, coordx = coords, coordt = coordt,
                            distance = distance, neighb = neighb,
                            coordx_dyn = coordx_dyn, loc = loc, time = time,
                            maxdist = maxdist, maxtime = maxtime, X = X, M = M,
                            parallel = FALSE, ncores = ncores)

    total_tasks = Nloc * Tloc
    res1 = res2 = double(total_tasks)

    if(progress){
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = 1:total_tasks)
    } else pb <- function(...) {}

    if(!parallel){
      k = 1
      for(i in 1:Nloc){
        for(j in 1:Tloc){
          pb(sprintf("k=%d", k))
          if(!is.null(M)) param$mean = neigh$M[[k]]
          pr = GeoKrig(data = neigh$data[[k]], coordx = neigh$coordx[[k]],
                       coordt = neigh$coordt[[k]], loc = loc[i, ],
                       time = time[j], X = neigh$X[[k]],
                       Mloc = Mloc[i + Nloc*(j-1)],
                       Xloc = Xloc[i + Nloc*(j-1), ], type_krig = type_krig,
                       sparse = sparse, corrmodel = corrmodel,
                       distance = distance, model = model, param = param,
                       anisopars = anisopars, radius = radius, mse = mse,
                       varcov=varcov,
                       copula = copula, n = n)
          res1[k] = pr$pred
          if(mse) res2[k] = pr$mse
          k = k + 1
        }
      }
    } else {
      config = setup_parallel(ncores, total_tasks)
      max_cores = config$n.cores
      safe_cores = choose_workers_safe(neigh, total_tasks, max_cores)
      if(safe_cores < max_cores){
        warning("Reducing cores from ", max_cores, " to ", safe_cores,
                " to avoid memory pressure.")
        future::plan(future::multisession, workers = safe_cores)
      }
       cat("Performing local kriging using ", safe_cores, " cores...\n")

      res = future.apply::future_lapply(seq_len(total_tasks), function(k){
        suppressMessages({
          local_param = param
          if(!is.null(M)) local_param$mean = neigh$M[[k]]
          j = ((k-1) %% Tloc) + 1
          i = ((k-1) %/% Tloc) + 1
          pr = GeoKrig(data = neigh$data[[k]], coordx = neigh$coordx[[k]],
                       coordt = neigh$coordt[[k]], loc = loc[i, ],
                       time = time[j], X = neigh$X[[k]],
                       Mloc = Mloc[k],
                       Xloc = Xloc[k, ], type_krig = type_krig,
                       sparse = sparse, corrmodel = corrmodel,
                       distance = distance, model = model, param = local_param,
                       anisopars = anisopars, radius = radius, mse = mse,varcov=varcov,
                       copula = copula, n = n)
        })
        if(progress) pb(sprintf("k=%d", k))
        list(pred = pr$pred, mse = if(mse) pr$mse else NULL)
      }, future.seed = TRUE, future.stdout = FALSE,
        future.conditions = "none")

 

      res1 = sapply(res, `[[`, "pred")
      if(mse) res2 = sapply(res, `[[`, "mse")

      cleanup_parallel(config)
    }
  }

  ## -------- bivariate -----------------------------------------------------
  if(bivariate){
    neigh = GeoNeighborhood(data, coordx = coords, distance = distance,
                            neighb = neighb, coordx_dyn = coordx_dyn,
                            loc = loc, maxdist = maxdist, X = X, M = M,
                            parallel = FALSE, ncores = ncores)

    total_tasks = Nloc * Tloc
    res1 = res2 = double(total_tasks)

    if(progress){
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = 1:total_tasks)
    } else pb <- function(...) {}

    if(!parallel){
      k = 1
      for(i in 1:Nloc){
        for(j in 1:Tloc){
          pb(sprintf("k=%d", k))
          if(!is.null(M)) param$mean = neigh$M[[k]]
          pr = GeoKrig(data = neigh$data[[k]], coordx = neigh$coordx[[k]],
                       loc = loc[i, ], X = neigh$X[[k]],
                       Mloc = Mloc[i + Nloc*(j-1)],
                       Xloc = Xloc[i + Nloc*(j-1), ], type_krig = type_krig,
                       sparse = sparse, corrmodel = corrmodel,
                       distance = distance, model = model, param = param,varcov=varcov,
                       anisopars = anisopars, radius = radius, mse = mse,
                       copula = copula, n = n)
          res1[k] = pr$pred
          if(mse) res2[k] = pr$mse
          k = k + 1
        }
      }
    } else {
      config = setup_parallel(ncores, total_tasks)
      max_cores = config$n.cores
      safe_cores = choose_workers_safe(neigh, total_tasks, max_cores)
      if(safe_cores < max_cores){
        warning("Reducing cores from ", max_cores, " to ", safe_cores,
                " to avoid memory pressure.")
        future::plan(future::multisession, workers = safe_cores)
      }
       cat("Performing local kriging using ", safe_cores, " cores...\n")

      res = future.apply::future_lapply(seq_len(total_tasks), function(k){
        suppressMessages({
          local_param = param
          if(!is.null(M)) local_param$mean = neigh$M[[k]]
          j = ((k-1) %% Tloc) + 1
          i = ((k-1) %/% Tloc) + 1
          pr = GeoKrig(data = neigh$data[[k]], coordx = neigh$coordx[[k]],
                       loc = loc[i, ], X = neigh$X[[k]],
                       Mloc = Mloc[k],
                       Xloc = Xloc[k, ], type_krig = type_krig,
                       sparse = sparse, corrmodel = corrmodel,
                       distance = distance, model = model, param = local_param,varcov=varcov,
                       anisopars = anisopars, radius = radius, mse = mse,
                       copula = copula, n = n)
        })
        if(progress) pb(sprintf("k=%d", k))
        list(pred = pr$pred, mse = if(mse) pr$mse else NULL)
      }, future.seed = TRUE, future.stdout = FALSE,
        future.conditions = "none")

      res1 = sapply(res, `[[`, "pred")
      if(mse) res2 = sapply(res, `[[`, "mse")

      cleanup_parallel(config)
    }
  }

  ###########################################################################
  ##  3.  output (unchanged)                                              ##
  ###########################################################################
  varpred = NULL
  if(spacetime || bivariate){
    pred = matrix(t(res1), nrow = Tloc, ncol = Nloc)
    if(mse) varpred = matrix(c(res2), nrow = Tloc, ncol = Nloc)
  } else {
    pred = c(res1)
    if(mse) varpred = c(res2)
  }
  if(Tloc == 1){ pred = c(pred); varpred = c(varpred) }

  GeoKrigloc_out = list(
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
    varcov=varcov
  )
  structure(c(GeoKrigloc_out, call = call), class = c("GeoKrigloc"))
}