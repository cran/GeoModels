GeoKriglocWeights <- function(coordx, coordy=NULL, coordz=NULL, coordt=NULL, coordx_dyn=NULL,
                             corrmodel, distance="Eucl", grid=FALSE, loc, neighb=NULL,
                             maxdist=NULL, maxtime=NULL, method="cholesky", model="Gaussian",
                             n=1, nloc=NULL, param, anisopars=NULL,
                             radius=1, sparse=FALSE, time=NULL, which=1,
                             copula=NULL, X=NULL, Xloc=NULL, Mloc=NULL,
                             parallel=TRUE)  
{
  ###########################################################################
  ##  Wrapper per pesi salvati su file temporanei
  ###########################################################################
  KrigWeightRef <- function(file) {
    structure(list(file = file), class = "KrigWeightRef")
  }
  `$.KrigWeightRef` <- function(x, name) {
    obj <- readRDS(x$file)  
    obj[[name]]
  }
  `[[.KrigWeightRef` <- function(x, i, ...) {
    obj <- readRDS(x$file)
    obj[[i]]
  }
  print.KrigWeightRef <- function(x, ...) {
    cat("KrigWeightRef object (lazy-loaded from file):", x$file, "\n")
  }
  save_tmp <- function(obj) {
    f <- tempfile(fileext = ".rds")
    saveRDS(obj, f)
    f
  }

  ###########################################################################
  ##  Validazione input
  ###########################################################################
  call = match.call()
  
  if(is.null(CkModel(model))) stop("The name of the model is not correct\n")
  if(!is.character(corrmodel) || is.null(CkCorrModel(corrmodel))) stop("the name of the correlation model is wrong")
  
  corrmodel <- gsub("[[:blank:]]", "", corrmodel)
  model <- gsub("[[:blank:]]", "", model)
  distance <- gsub("[[:blank:]]", "", distance)
  method <- gsub("[[:blank:]]", "", method)

  bivariate = CheckBiv(CkCorrModel(corrmodel))
  spacetime = CheckST(CkCorrModel(corrmodel))
  space = !spacetime && !bivariate

  # coordinate
  if(is.null(coordx_dyn)){
    if(is.null(coordy)) {
      coords = as.matrix(coordx)
    } else {
      if(!grid) {
        coords = cbind(coordx, coordy)
        if(!is.null(coordz)) coords = cbind(coords, coordz)
      } else {
        if(!is.null(coordz)) {
          coords = as.matrix(expand.grid(coordx, coordy, coordz))
        } else {
          coords = as.matrix(expand.grid(coordx, coordy))
        }
      }
    }
  } else {
    coords = NULL
  }

  if(is.vector(loc)) loc = matrix(loc, nrow=1)
  loc = as.matrix(loc)
  Nloc = nrow(loc)
  
  if(is.null(time)) time = 0
  Tloc = length(time)
  if(bivariate) Tloc = 1
  
  M = NULL
  if(length(param$mean) > 1) M = param$mean

  ###########################################################################
  ##  Loop Kriging
  ###########################################################################
  
  process_loc <- function(i, neigh, coords, param, loc, X, Xloc, Mloc) {
    param_local = param
    if(!is.null(M) && !is.null(neigh$M)) param_local$mean = neigh$M[[i]]
    coords_neigh = neigh$coordx[[i]]
    if(is.null(coords_neigh) || nrow(coords_neigh) == 0) return(NULL)
    weights_result = GeoKrigWeights(
      coordx = coords_neigh,
      loc = matrix(loc[i, ], nrow=1), 
      corrmodel = corrmodel,
      distance = distance, 
      n = n, 
      X = if(!is.null(neigh$X)) neigh$X[[i]] else NULL,
      Xloc = if(!is.null(Xloc)) matrix(Xloc[i, ], nrow=1) else NULL,
      Mloc = if(!is.null(Mloc)) Mloc[i] else NULL,
      sparse = sparse, 
      model = model, 
      param = param_local,
      anisopars = anisopars, 
      radius = radius,
      copula = copula, 
      method = method,
      which = which
    )
    weights_result$neighbor_indices = neigh$indices[[i]]
    KrigWeightRef(save_tmp(weights_result))
  }

  ## -------- SPACE CASE --------------------------------------------------------
  if(space){
    neigh = GeoNeighborhood(coordx = coords, distance = distance, loc = loc,
                            neighb = neighb, maxdist = maxdist, X = X, M = M,
                            parallel = FALSE, ncores = NULL)
    if(parallel){

      plan(multisession, workers = max(1, parallel::detectCores() - 1))
      weights_list <- future.apply::future_lapply(seq_len(Nloc), function(i){
        process_loc(i, neigh, coords, param, loc, X, Xloc, Mloc)
      },future.seed=TRUE)
    } else {
      weights_list = vector("list", Nloc)
      for(i in seq_len(Nloc)){
        weights_list[[i]] <- process_loc(i, neigh, coords, param, loc, X, Xloc, Mloc)
      }
    }
  }

  ## -------- SPACETIME CASE -----------------------------------------------------
  if(spacetime){
    neigh = GeoNeighborhood(coordx = coords, coordt = coordt,
                            distance = distance, neighb = neighb,
                            coordx_dyn = coordx_dyn, loc = loc, time = time,
                            maxdist = maxdist, maxtime = maxtime, X = X, M = M,
                            parallel = FALSE, ncores = NULL)
    total_tasks = Nloc * Tloc
    weights_list = vector("list", total_tasks)
    k = 1
    for(i in 1:Nloc){
      for(j in 1:Tloc){
        param_local = param
        if(!is.null(M) && !is.null(neigh$M)) param_local$mean = neigh$M[[k]]
        coords_neigh = neigh$coordx[[k]]
        if(!is.null(coords_neigh) && length(coords_neigh) > 0){
          weights_result = GeoKrigWeights(
            coordx = coords_neigh, 
            coordt = neigh$coordt[[k]], 
            loc = matrix(loc[i, ], nrow=1),
            time = time[j], 
            X = if(!is.null(neigh$X)) neigh$X[[k]] else NULL,
            Mloc = if(!is.null(Mloc)) Mloc[i + Nloc*(j-1)] else NULL,
            Xloc = if(!is.null(Xloc)) matrix(Xloc[i + Nloc*(j-1), ], nrow=1) else NULL,
            sparse = sparse, 
            corrmodel = corrmodel,
            distance = distance, 
            model = model, 
            param = param_local,
            anisopars = anisopars, 
            radius = radius,
            copula = copula, 
            n = n, 
            method = method,
            which = which
          )
          weights_result$neighbor_indices = neigh$indices[[k]]
          weights_list[[k]] = KrigWeightRef(save_tmp(weights_result))
        }
        k = k + 1
      }
    }
  }

  ## -------- BIVARIATE CASE -----------------------------------------------------
  if(bivariate){
    neigh = GeoNeighborhood(coordx = coords, distance = distance,
                            neighb = neighb, coordx_dyn = coordx_dyn,
                            loc = loc, maxdist = maxdist, X = X, M = M,
                            parallel = FALSE, ncores = NULL)
    if(parallel){

      plan(multisession, workers = max(1, parallel::detectCores() - 1))
      weights_list <- future.apply::future_lapply(seq_len(Nloc), function(i){
        process_loc(i, neigh, coords, param, loc, X, Xloc, Mloc)
      },future.seed=TRUE)
    } else {
      weights_list = vector("list", Nloc)
      for(i in seq_len(Nloc)){
        weights_list[[i]] <- process_loc(i, neigh, coords, param, loc, X, Xloc, Mloc)
      }
    }
  }

  ###########################################################################
  ##  Output
  ###########################################################################
  GeoKriglocWeights_out = list(
    bivariate = bivariate,
    coordx = coordx,
    coordy = coordy,
    coordz = coordz,
    coordt = coordt,
    corrmodel = corrmodel,
    data = NULL,
    distance = distance,
    grid = grid,
    loc = loc,
    copula = copula,
    numcoord = if(!is.null(coords)) nrow(coords) else length(coordx),
    numloc = Nloc,
    numtime = if(!is.null(coordt)) length(coordt) else 0,
    numt = Tloc,
    maxdist = maxdist,
    maxtime = maxtime,
    model = model,
    n = n,
    param = param,
    weights = weights_list,
    radius = radius,
    spacetime = spacetime,
    time = time,
    neighb = neighb
  )
  
  structure(c(GeoKriglocWeights_out, call = call), class = c("GeoKriglocWeights"))
}
