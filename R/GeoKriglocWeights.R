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

  # coordinate osservate
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
  ##  SUGG. 1: anisotropia anche per la selezione dei vicini
  ###########################################################################
  coords_aniso <- coords
  loc_aniso    <- loc
  if (!is.null(anisopars)) {
    if (!is.null(coords)) coords_aniso <- GeoAniso(coords, c(anisopars$angle, anisopars$ratio))
    loc_aniso <- GeoAniso(loc, c(anisopars$angle, anisopars$ratio))
  }

  ###########################################################################
  ##  Gestione memoria + plan centralizzato (SUGG. 5)
  ###########################################################################
  if(parallel) {
    old_limit <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = 3000 * 1024^2)  # 3 GB
    on.exit(options(future.globals.maxSize = old_limit), add = TRUE)
  }
  .plan_set <- FALSE
  .ensure_plan <- function() {
    if (parallel && !.plan_set) {
      future::plan(future::multisession, workers = max(1, parallel::detectCores() - 1))
      on.exit(future::plan(future::sequential), add = TRUE)
      .plan_set <<- TRUE
    }
  }

  ###########################################################################
  ##  Loop Kriging - SPACE
  ###########################################################################
  if(space){
    neigh = GeoNeighborhood(coordx = coords_aniso, distance = distance, loc = loc_aniso,
                            neighb = neighb, maxdist = maxdist, X = X, M = M,
                            parallel = FALSE, ncores = NULL)
    
    if(parallel){
      .ensure_plan()
      # Variabili essenziali
      essential_vars <- list(
        neigh_coordx = neigh$coordx,
        neigh_indices = neigh$indices,
        neigh_X = neigh$X,
        neigh_M = neigh$M,
        loc = loc,
        Xloc = Xloc,
        Mloc = Mloc,
        param = param,
        M = M,
        corrmodel = corrmodel,
        distance = distance,
        n = n,
        sparse = sparse,
        model = model,
        anisopars = anisopars,
        radius = radius,
        copula = copula,
        method = method,
        which = which
      )
      
      worker_fn <- function(i, vars) {
        # estrai vicini
        coords_neigh <- vars$neigh_coordx[[i]]
        # SUGG. 4: check robusto
        if (is.null(coords_neigh) || (is.matrix(coords_neigh) && nrow(coords_neigh) == 0) || length(coords_neigh) == 0) {
          return(NULL)
        }
        # mean locale
        param_local <- vars$param
        if(!is.null(vars$M) && !is.null(vars$neigh_M)) {
          param_local$mean <- vars$neigh_M[[i]]
        }
        # SUGG. 3: Xloc riga robusta a vettore/matrice
        xloc_row <- if (!is.null(vars$Xloc)) {
          if (is.vector(vars$Xloc)) matrix(vars$Xloc[i], nrow = 1) else matrix(vars$Xloc[i, ], nrow = 1)
        } else NULL
        
        weights_result <- GeoKrigWeights(
          coordx = coords_neigh,
          loc = matrix(vars$loc[i, ], nrow=1), 
          corrmodel = vars$corrmodel,
          distance = vars$distance, 
          n = vars$n, 
          X = if(!is.null(vars$neigh_X)) vars$neigh_X[[i]] else NULL,
          Xloc = xloc_row,
          Mloc = if(!is.null(vars$Mloc)) vars$Mloc[i] else NULL,
          sparse = vars$sparse, 
          model = vars$model, 
          param = param_local,
          anisopars = vars$anisopars, 
          radius = vars$radius,
          copula = vars$copula, 
          method = vars$method,
          which = vars$which
        )
        weights_result$neighbor_indices <- vars$neigh_indices[[i]]
        return(KrigWeightRef(save_tmp(weights_result)))
      }
      
      weights_list <- future.apply::future_lapply(
        seq_len(Nloc), 
        worker_fn,
        vars = essential_vars,
        future.seed = TRUE,
        future.globals = FALSE
      )
      
    } else {
      weights_list = vector("list", Nloc)
      for(i in seq_len(Nloc)){
        param_local = param
        if(!is.null(M) && !is.null(neigh$M)) param_local$mean = neigh$M[[i]]
        coords_neigh = neigh$coordx[[i]]
        # SUGG. 4: check robusto
        if (is.null(coords_neigh) || (is.matrix(coords_neigh) && nrow(coords_neigh) == 0) || length(coords_neigh) == 0) {
          weights_list[[i]] <- NULL
          next
        }
        # SUGG. 3: Xloc riga robusta
        xloc_row <- if (!is.null(Xloc)) {
          if (is.vector(Xloc)) matrix(Xloc[i], nrow = 1) else matrix(Xloc[i, ], nrow = 1)
        } else NULL
        
        weights_result = GeoKrigWeights(
          coordx = coords_neigh,
          loc = matrix(loc[i, ], nrow=1), 
          corrmodel = corrmodel,
          distance = distance, 
          n = n, 
          X = if(!is.null(neigh$X)) neigh$X[[i]] else NULL,
          Xloc = xloc_row,
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
        weights_list[[i]] = KrigWeightRef(save_tmp(weights_result))
      }
    }
  }

  ###########################################################################
  ##  Loop Kriging - SPACETIME
  ###########################################################################
  if(spacetime){
    neigh = GeoNeighborhood(coordx = coords_aniso, coordt = coordt,
                            distance = distance, neighb = neighb,
                            coordx_dyn = coordx_dyn, loc = loc_aniso, time = time,
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
        # SUGG. 4: check robusto
        if (!is.null(coords_neigh) && !(is.matrix(coords_neigh) && nrow(coords_neigh) == 0) && length(coords_neigh) > 0){
          # SUGG. 3: Xloc riga robusta con indice spacetime
          idx_st <- i + Nloc*(j-1)
          xloc_row <- if (!is.null(Xloc)) {
            if (is.vector(Xloc)) matrix(Xloc[idx_st], nrow = 1) else matrix(Xloc[idx_st, ], nrow = 1)
          } else NULL
          
          weights_result = GeoKrigWeights(
            coordx = coords_neigh, 
            coordt = neigh$coordt[[k]], 
            loc = matrix(loc[i, ], nrow=1),
            time = time[j], 
            X = if(!is.null(neigh$X)) neigh$X[[k]] else NULL,
            Mloc = if(!is.null(Mloc)) Mloc[idx_st] else NULL,
            Xloc = xloc_row,
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
        } else {
          weights_list[[k]] <- NULL
        }
        k = k + 1
      }
    }
  }

  ###########################################################################
  ##  Loop Kriging - BIVARIATE
  ###########################################################################
  if(bivariate){
    neigh = GeoNeighborhood(coordx = coords_aniso, distance = distance,
                            neighb = neighb, coordx_dyn = coordx_dyn,
                            loc = loc_aniso, maxdist = maxdist, X = X, M = M,
                            parallel = FALSE, ncores = NULL)
    
    if(parallel){
      .ensure_plan()
      essential_vars <- list(
        neigh_coordx = neigh$coordx,
        neigh_indices = neigh$indices,
        neigh_X = neigh$X,
        neigh_M = neigh$M,
        loc = loc,
        Xloc = Xloc,
        Mloc = Mloc,
        param = param,
        M = M,
        corrmodel = corrmodel,
        distance = distance,
        n = n,
        sparse = sparse,
        model = model,
        anisopars = anisopars,
        radius = radius,
        copula = copula,
        method = method,
        which = which
      )
      
      worker_fn <- function(i, vars) {
        coords_neigh <- vars$neigh_coordx[[i]]
        # SUGG. 4: check robusto
        if (is.null(coords_neigh) || (is.matrix(coords_neigh) && nrow(coords_neigh) == 0) || length(coords_neigh) == 0) {
          return(NULL)
        }
        param_local <- vars$param
        if(!is.null(vars$M) && !is.null(vars$neigh_M)) {
          param_local$mean <- vars$neigh_M[[i]]
        }
        # SUGG. 3: Xloc riga robusta
        xloc_row <- if (!is.null(vars$Xloc)) {
          if (is.vector(vars$Xloc)) matrix(vars$Xloc[i], nrow = 1) else matrix(vars$Xloc[i, ], nrow = 1)
        } else NULL
        
        weights_result <- GeoKrigWeights(
          coordx = coords_neigh,
          loc = matrix(vars$loc[i, ], nrow=1), 
          corrmodel = vars$corrmodel,
          distance = vars$distance, 
          n = vars$n, 
          X = if(!is.null(vars$neigh_X)) vars$neigh_X[[i]] else NULL,
          Xloc = xloc_row,
          Mloc = if(!is.null(vars$Mloc)) vars$Mloc[i] else NULL,
          sparse = vars$sparse, 
          model = vars$model, 
          param = param_local,
          anisopars = vars$anisopars, 
          radius = vars$radius,
          copula = vars$copula, 
          method = vars$method,
          which = vars$which
        )
        weights_result$neighbor_indices <- vars$neigh_indices[[i]]
        return(KrigWeightRef(save_tmp(weights_result)))
      }
      
      weights_list <- future.apply::future_lapply(
        seq_len(Nloc), 
        worker_fn,
        vars = essential_vars,
        future.seed = TRUE,
        future.globals = FALSE
      )
      
    } else {
      weights_list = vector("list", Nloc)
      for(i in seq_len(Nloc)){
        param_local = param
        if(!is.null(M) && !is.null(neigh$M)) param_local$mean = neigh$M[[i]]
        coords_neigh = neigh$coordx[[i]]
        # SUGG. 4: check robusto
        if (is.null(coords_neigh) || (is.matrix(coords_neigh) && nrow(coords_neigh) == 0) || length(coords_neigh) == 0) {
          weights_list[[i]] <- NULL
          next
        }
        # SUGG. 3: Xloc riga robusta
        xloc_row <- if (!is.null(Xloc)) {
          if (is.vector(Xloc)) matrix(Xloc[i], nrow = 1) else matrix(Xloc[i, ], nrow = 1)
        } else NULL
        
        weights_result = GeoKrigWeights(
          coordx = coords_neigh,
          loc = matrix(loc[i, ], nrow=1), 
          corrmodel = corrmodel,
          distance = distance, 
          n = n, 
          X = if(!is.null(neigh$X)) neigh$X[[i]] else NULL,
          Xloc = xloc_row,
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
        weights_list[[i]] = KrigWeightRef(save_tmp(weights_result))
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
