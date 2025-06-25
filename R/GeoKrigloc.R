GeoKrigloc = function(estobj=NULL, data, coordx, coordy=NULL, coordz=NULL, coordt=NULL, coordx_dyn=NULL, corrmodel, distance="Eucl", grid=FALSE, loc, neighb=NULL,
                      maxdist=NULL, maxtime=NULL, method="cholesky", model="Gaussian", n=1, nloc=NULL, mse=FALSE,  param, anisopars=NULL, 
                      radius=6371, sparse=FALSE, time=NULL, type="Standard",
                      type_mse=NULL, type_krig="Simple", weigthed=TRUE, which=1, copula=NULL, X=NULL, Xloc=NULL, Mloc=NULL,
                      spobj=NULL, spdata=NULL, parallel=FALSE, ncores=NULL, progress=TRUE) {

  call = match.call()    

  # --- Estrazione da GeoFit ---
  if(!is.null(estobj)){
    if(!inherits(estobj,"GeoFit")) stop("need a 'GeoFit' object as input\n")
    data = estobj$data
    if(!estobj$grid){
      if(!estobj$bivariate){
        if(is.null(estobj$coordx_dyn)) coordx = cbind(estobj$coordx, estobj$coordy, estobj$coordz)
        else coordx_dyn = estobj$coordx_dyn
      } else {
        if(is.null(estobj$coordx_dyn)) { 
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
  }

  # --- Coordinate e setup ---
  bivariate = CheckBiv(CkCorrModel(corrmodel))
  spacetime = CheckST(CkCorrModel(corrmodel))
  space = !spacetime && !bivariate

  if(!is.null(spobj)) {
    if(space || bivariate){
      a = sp2Geo(spobj, spdata); coordx = a$coords 
      if(!a$pj) {if(distance != "Chor") distance = "Geod"}
    }
    if(spacetime){
      a = sp2Geo(spobj, spdata); coordx = a$coords ; coordt = a$coordt 
      if(!a$pj) {if(distance != "Chor") distance = "Geod"}
    }
    if(!is.null(a$Y) && !is.null(a$X)) { data = a$Y ; X = a$X }
  }

  if(is.null(coordx_dyn)){
    coords = coordx
    if(!is.null(coordy)){
      if(!grid) coords = cbind(coordx, coordy, coordz)
      else {
        if(!is.null(coordz)) coords = as.matrix(expand.grid(coordx, coordy, coordz))
        else coords = as.matrix(expand.grid(coordx, coordy))
      }
    }
  } else {
    coordx = NULL; coordy = NULL; coordz = NULL; coords = NULL
  }

  Nloc = nrow(loc)
  if(is.null(Nloc)) Nloc = 1
  Tloc = length(time)
  if(bivariate) Tloc = 1
  if(length(param$mean) > 1) M = param$mean else M = NULL

  coremax = parallel::detectCores()
  if(is.na(coremax) || coremax == 1) parallel = FALSE

  #------------------- SPAZIO ----------------------
  if (space) {
    neigh <- GeoNeighborhood(
      data, coordx = coords, distance = distance, loc = loc, neighb = neighb,
      maxdist = maxdist, X = X, M = M, parallel = FALSE, ncores = ncores
    )

    res1 <- numeric(Nloc)
    res2 <- if (mse) numeric(Nloc) else NULL

    if(progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = seq_len(Nloc))
    } else {
      pb <- function(...) {}
    }

    if (!parallel) {
      for (i in seq_len(Nloc)) {
       pb(sprintf("i=%d", i))
        if (!is.null(M)) param$mean <- neigh$M[[i]]
        pr <- GeoKrig(  loc = loc[i, ], data = neigh$data[[i]], coordx = neigh$coordx[[i]],
          corrmodel = corrmodel, distance = distance, n = n,
          X = neigh$X[[i]], Xloc = Xloc[i, ], Mloc = Mloc[i], type_krig = type_krig,
          sparse = sparse, model = model, param = param, anisopars = anisopars,
          mse = mse, copula = copula
        )
        res1[i] <- pr$pred
        if (mse) res2[i] <- pr$mse
      }
    } else {
      if (is.null(ncores)) {
        n.cores <- max(1, future::availableCores() - 1)
      } else {
        if (!is.numeric(ncores) || ncores < 1 || ncores > future::availableCores()) {
          stop("Invalid number of cores specified")
        }
        n.cores <- ncores
      }
      future::plan(multisession, workers = n.cores)
      on.exit(future::plan(sequential), add = TRUE)
      cat("Performing local kriging using ", n.cores, " cores...\n")
      oopts <- options(future.globals.maxSize = 8L * 1024^3)
      on.exit(options(oopts), add = TRUE)

      xx <- foreach::foreach(
        i = seq_len(Nloc),
        .combine = rbind,
        .options.future = list(seed = TRUE, globals = structure(TRUE, add = "param"),stdout = NA,  
                        conditions = character(0))) %dofuture% {
        pb(sprintf("i=%d", i))
        if (!is.null(M)) param$mean <- neigh$M[[i]]
        pr <- GeoKrig(
           loc = loc[i, ], data = neigh$data[[i]], coordx = neigh$coordx[[i]],
          corrmodel = corrmodel, distance = distance, n = n,
          X = neigh$X[[i]], Xloc = Xloc[i, ], Mloc = Mloc[i], type_krig = type_krig,
          sparse = sparse, model = model, param = param, anisopars = anisopars,
          mse = mse, copula = copula
        )
        c(pr$pred, if (mse) pr$mse else NA_real_)
      }
      xx <- matrix(xx, nrow = Nloc)
      res1 <- as.numeric(xx[, 1])
      if (mse) res2 <- as.numeric(xx[, 2]) else res2 <- NULL
    }
  }

  #------------------- SPAZIO-TEMPO ----------------------
  if(spacetime) {  
    neigh = GeoNeighborhood(data, coordx=coords, coordt=coordt, distance=distance, neighb=neighb, coordx_dyn=coordx_dyn,
                           loc=loc, time=time, maxdist=maxdist, maxtime=maxtime, X=X, M=M, parallel=FALSE, ncores=ncores)
    res1 = res2 = double(Nloc * Tloc)
    k = 1

    if(progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = 1:(Nloc * Tloc))
    } else {
      pb <- function(...) {}
    }

    if(!parallel) {
      for(i in 1:Nloc){
        for(j in 1:Tloc){
          pb(sprintf("k=%d", k))
          if(!is.null(M)) param$mean=neigh$M[[k]]
          pr = GeoKrig( data=neigh$data[[k]], coordx=neigh$coordx[[k]], coordt=neigh$coordt[[k]], loc=loc[i,], time=time[j],
                       X=neigh$X[[k]], Mloc=Mloc[i+(Nloc)*(j-1)], Xloc=Xloc[i+(Nloc)*(j-1),], type_krig=type_krig, sparse=sparse,
                       corrmodel=corrmodel, distance=distance, model=model, param=param, anisopars=anisopars, mse=mse, copula=copula, n=n)
          res1[k] = pr$pred
          if(mse) res2[k] = pr$mse
          k = k+1
        }
      }
    } else {
      if (is.null(ncores)) {
        n.cores <- max(1, future::availableCores() - 1)
      } else {
        if (!is.numeric(ncores) || ncores < 1 || ncores > future::availableCores()) {
          stop("Invalid number of cores specified")
        }
        n.cores <- ncores
      }
      future::plan(multisession, workers = n.cores)
      on.exit(future::plan(sequential), add = TRUE)
      cat("Performing local spatio-temporal kriging using ", n.cores, " cores...\n")
      oopts <- options(future.globals.maxSize = 8L * 1024^3)
      on.exit(options(oopts), add = TRUE)

      xx <- foreach::foreach(
        k = seq_len(Nloc*Tloc),
        .combine = rbind,
        .options.future = list(seed = TRUE, globals = structure(TRUE, add = "param"),stdout = NA,  
                        conditions = character(0))) %dofuture% {
        pb(sprintf("k=%d", k))
        if(!is.null(M)) param$mean=neigh$M[[k]]
        pr = GeoKrig(data=neigh$data[[k]], coordx=neigh$coordx[[k]], coordt=neigh$coordt[[k]], loc=loc[ ((k-1) %% Nloc) + 1,], time=time[ ((k-1) %/% Nloc) + 1],
                     X=neigh$X[[k]], Mloc=Mloc[k], Xloc=Xloc[k,], type_krig=type_krig, sparse=sparse,
                     corrmodel=corrmodel, distance=distance, model=model, param=param, anisopars=anisopars, mse=mse, copula=copula, n=n)
        c(pr$pred, if (mse) pr$mse else NA_real_)
      }
      xx <- matrix(xx, nrow = Nloc*Tloc)
      res1 <- as.numeric(xx[, 1])
      if (mse) res2 <- as.numeric(xx[, 2]) else res2 <- NULL
    }
  }

  #------------------- BIVARIATO ----------------------
  if(bivariate) {
    neigh = GeoNeighborhood(data, coordx=coords, distance=distance, neighb=neighb, coordx_dyn=coordx_dyn,
                           loc=loc, maxdist=maxdist, X=X, M=M, parallel=FALSE, ncores=ncores)
    res1 = res2 = double(Nloc * Tloc)
    k = 1

    if(progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = 1:(Nloc * Tloc))
    } else {
      pb <- function(...) {}
    }

    if(!parallel) {
      for(i in 1:Nloc){
        for(j in 1:Tloc){
          pb(sprintf("k=%d", k))
          if(!is.null(M)) param$mean=neigh$M[[k]]
          pr = GeoKrig( data=neigh$data[[k]], coordx=neigh$coordx[[k]], loc=loc[i,], 
                       X=neigh$X[[k]], Mloc=Mloc[i+(Nloc)*(j-1)], Xloc=Xloc[i+(Nloc)*(j-1),], type_krig=type_krig, sparse=sparse,
                       corrmodel=corrmodel, distance=distance, model=model, param=param, anisopars=anisopars, mse=mse, copula=copula, n=n)
          res1[k] = pr$pred
          if(mse) res2[k] = pr$mse
          k = k+1
        }
      }
    } else {
      if (is.null(ncores)) {
        n.cores <- max(1, future::availableCores() - 1)
      } else {
        if (!is.numeric(ncores) || ncores < 1 || ncores > future::availableCores()) {
          stop("Invalid number of cores specified")
        }
        n.cores <- ncores
      }
      future::plan(multisession, workers = n.cores)
      on.exit(future::plan(sequential), add = TRUE)
      cat("Performing local bivariate kriging using ", n.cores, " cores...\n")
      oopts <- options(future.globals.maxSize = 8L * 1024^3)
      on.exit(options(oopts), add = TRUE)

      xx <- foreach::foreach(
        k = seq_len(Nloc*Tloc),
        .combine = rbind,
        .options.future = list(seed = TRUE, globals = structure(TRUE, add = "param"),stdout = NA,  
                        conditions = character(0))
      ) %dofuture% {
        pb(sprintf("k=%d", k))
        if(!is.null(M)) param$mean=neigh$M[[k]]
        pr = GeoKrig( data=neigh$data[[k]], coordx=neigh$coordx[[k]], loc=loc[ ((k-1) %% Nloc) + 1,], 
                     X=neigh$X[[k]], Mloc=Mloc[k], Xloc=Xloc[k,], type_krig=type_krig, sparse=sparse,
                     corrmodel=corrmodel, distance=distance, model=model, param=param, anisopars=anisopars, mse=mse, copula=copula, n=n)
        c(pr$pred, if (mse) pr$mse else NA_real_)
      }
      xx <- matrix(xx, nrow = Nloc*Tloc)
      res1 <- as.numeric(xx[, 1])
      if (mse) res2 <- as.numeric(xx[, 2]) else res2 <- NULL
    }
  }

  #------------------- OUTPUT ----------------------
  varpred=NULL
  if(spacetime||bivariate) {      
    pred=matrix(t(res1),nrow=Tloc,ncol=Nloc)
    if(mse) varpred=matrix(c(res2),nrow=Tloc,ncol=Nloc)
  } else {
    pred=c(res1)
    if(mse) varpred=c(res2)
  }
  if(Tloc==1)  {pred=c(pred);varpred=c(varpred)}

  GeoKrigloc_out = list(   
    bivariate=bivariate,
    coordx = coordx,
    coordy = coordy,
    coordz = coordz,
    coordt = coordt,
    corrmodel = corrmodel,
    data=data,
    distance = distance,
    grid=grid,
    loc=loc,
    copula=copula,
    numcoord = if(exists("coords")) nrow(coords) else NA,
    numloc= Nloc,
    numtime = length(coordt),
    numt = Tloc,
    maxdist=maxdist,
    maxtime=maxtime,
    model=model,
    n=n,
    param = param,
    pred=pred,
    radius=radius,
    spacetime = spacetime,
    time=time,
    type_krig=type_krig,
    mse=varpred
  )
  structure(c(GeoKrigloc_out, call = call), class = c("GeoKrigloc"))
}
