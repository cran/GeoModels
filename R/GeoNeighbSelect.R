GeoNeighbSelect <- function(data, coordx, coordy=NULL, coordz=NULL, coordt=NULL, coordx_dyn=NULL,
  copula=NULL, corrmodel=NULL, distance="Eucl", fixed=NULL, anisopars=NULL,
  est.aniso=c(FALSE,FALSE), grid=FALSE, likelihood='Marginal', 
  lower=NULL, neighb=c(1,2,3,4,5),
  maxtime=Inf, memdist=TRUE, model='Gaussian', n=1, ncores=NULL,
  optimizer='Nelder-Mead', parallel=FALSE, bivariate=FALSE,radius=1, start=NULL, type='Pairwise', upper=NULL, weighted=FALSE,
  X=NULL, nosym=FALSE, spobj=NULL, spdata=NULL, vario=NULL)
{

  if(!is.numeric(neighb))  stop("neighb must be a numeric vector")
  if(sum(neighb-floor(neighb))) stop("neighb must be a positive integer numeric vector")

  estimates = best_T = best_K = NULL

  if(!is.null(try(CkCorrModel(corrmodel), TRUE))) {
    bivariate <- CheckBiv(CkCorrModel(corrmodel))
    spacetime <- CheckST(CkCorrModel(corrmodel))
    space = !spacetime && !bivariate
  } else { stop("correlation model is not valid\n") }

  if(!is.null(vario)) {
    if(!inherits(vario,"GeoVariogram"))  stop("A GeoVariogram object is needed as input for vario\n")
    if( !((vario$bivariate&&bivariate)||(!is.null(vario$bint)&&spacetime)||(is.null(vario$bint)&&!bivariate)) )
      stop("The GeoVariogram object is not of the same type of the correlation model\n")
    semiv=vario
  } else stop("A GeoVariogram object is needed as input for vario\n")

  coremax = parallel::detectCores()
  if(is.na(coremax)||coremax==1) parallel=FALSE
  K = length(neighb)
  res = double(K)
  P = NULL
  if(spacetime){
    if(!is.numeric(maxtime))  stop("maxtime must be a numeric vector")
    P = length(maxtime)
    res = double(K*P)
  }

  #################### SPATIAL ####################
  if(space||bivariate) {
    estimates = matrix(NA, nrow=K, ncol=length(start))
    res = rep(NA, K)
    if(!parallel) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = 1:K)
      for(M in 1:K) {
        pb(sprintf("k=%g", M)) 
        aa = GeoFit(data=data, coordx=coordx, coordy=coordy, coordz=coordz, coordt=coordt, coordx_dyn=coordx_dyn, copula=copula, corrmodel=corrmodel, distance=distance,
                    fixed=fixed, anisopars=anisopars, est.aniso=est.aniso, grid=grid, likelihood=likelihood, 
                    lower=lower, neighb=neighb[M], maxtime=maxtime, memdist=memdist, model=model, n=n, 
                    optimizer=optimizer, radius=radius, start=start, type=type, upper=upper, weighted=weighted, X=X, nosym=nosym, spobj=spobj, spdata=spdata)
        clest = aa$param
        estimates[M, ] = unlist(clest)
        if(aa$convergence=="Successful") {
          # Metodo alternativo per la valutazione
          cc = GeoCorrFct(semiv$centers, t=semiv$centert, corrmodel=corrmodel, model=model, distance=distance, 
                          param=c(aa$param, aa$fixed), radius=radius, n=n, covariance=TRUE, variogram=TRUE)$corr
          res[M] = sum((cc - semiv$variograms)^2)
        } else { res[M] = Inf }
      }
    } else {
      if(is.null(ncores)){ n.cores <- coremax - 1 }
      else {
        if(!is.numeric(ncores)) stop("number of cores not valid\n")
        if(ncores>coremax||ncores<1) stop("number of cores not valid\n")
        n.cores = ncores
      }
      cat("Performing", K, "estimations using", n.cores, "cores...\n")
      future::plan(future::multisession, workers = n.cores)
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = 1:K)
      results <- foreach::foreach(M = 1:K, .combine = rbind, .options.future = list(seed = TRUE,stdout = NA,  
                        conditions = character(0))) %dofuture% {
        aa <- GeoFit(data = data, coordx = coordx, coordy = coordy, coordz=coordz, coordt = coordt, coordx_dyn = coordx_dyn,
                     copula = copula, corrmodel = corrmodel, distance = distance, fixed = fixed, anisopars = anisopars,
                     est.aniso = est.aniso, grid = grid, likelihood = likelihood, lower = lower, neighb = neighb[M],
                     maxtime = maxtime, memdist = memdist, model = model, n = n, optimizer = optimizer,
                     radius = radius, start = start, type = type, upper = upper, weighted = weighted,
                     X = X, nosym = nosym, spobj = spobj, spdata = spdata)
        estimates <- unlist(aa$param)
        if(aa$convergence == "Successful") {
          cc = GeoCorrFct(semiv$centers, t=semiv$centert, corrmodel=corrmodel, model=model, distance=distance, 
                          param=c(aa$param, aa$fixed), radius=radius, n=n, covariance=TRUE, variogram=TRUE)$corr
          res = sum((cc - semiv$variograms)^2)
        } else {
          res = Inf
        }
        pb(sprintf("k=%g", M))
        c(estimates, res = res)
      }
      estimates <- results[, -ncol(results), drop=FALSE]
      rownames(estimates) <- NULL
      res <- results[, ncol(results)]
      future::plan(sequential)
    }
  }

  #################### SPATIO-TEMPORAL ####################
  if(spacetime) {
    estimates = matrix(NA, nrow=K*P, ncol=length(start))
    res = rep(NA, K*P)
    if(!parallel) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = 1:(K*P))
      F=1
      for(L in 1:P) { 
        for(M in 1:K) {
          pb(sprintf("k=%g", F)) 
          aa = GeoFit(data=data, coordx=coordx, coordy=coordy, coordz=coordz, coordt=coordt, coordx_dyn=coordx_dyn, copula=copula, corrmodel=corrmodel, distance=distance,
                      fixed=fixed, anisopars=anisopars, est.aniso=est.aniso, grid=grid, likelihood=likelihood, 
                      lower=lower, neighb=neighb[M], maxtime=maxtime[L], memdist=memdist, model=model, n=n, 
                      optimizer=optimizer, radius=radius, start=start, type=type, upper=upper, weighted=weighted, X=X, nosym=nosym, spobj=spobj, spdata=spdata)
          estimates[F, ] = unlist(aa$param)
          if(aa$convergence=="Successful") {
            cc = GeoCorrFct(semiv$centers, t=semiv$centert, corrmodel=corrmodel, model=model, distance=distance, 
                            param=c(aa$param, aa$fixed), radius=radius, n=n, covariance=TRUE, variogram=TRUE)$corr
            res[F] = sum((cc - semiv$variograms)^2)
          } else { res[F] = Inf }
          F = F + 1
        }
      }
    } else {
      if(is.null(ncores)){ n.cores <- coremax - 1 }
      else {
        if(!is.numeric(ncores)) stop("number of cores not valid\n")
        if(ncores>coremax||ncores<1) stop("number of cores not valid\n")
        n.cores = ncores
      }
      cat("Performing", K, "estimations using", n.cores, "cores...\n")
      `%:%` <- foreach::`%:%`
      future::plan(future::multisession, workers = n.cores)
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = 1:(K*P))
      results <- foreach::foreach(L = 1:P, .combine = rbind) %:%
        foreach::foreach(M = 1:K, .combine = rbind, .options.future = list(seed = TRUE,stdout = NA,  
                        conditions = character(0))) %dofuture% {
          F = (L - 1) * K + M
          pb(sprintf("k=%g", F))
          aa = GeoFit(data = data, coordx = coordx, coordy = coordy, coordz = coordz, coordt = coordt, coordx_dyn = coordx_dyn, copula = copula, corrmodel = corrmodel, distance = distance,
                      fixed = fixed, anisopars = anisopars, est.aniso = est.aniso, grid = grid, likelihood = likelihood,
                      lower = lower, neighb = neighb[M], maxtime = maxtime[L], memdist = memdist, model = model, n = n,
                      optimizer = optimizer, radius = radius, start = start, type = type, upper = upper, weighted = weighted, X = X, nosym = nosym, spobj = spobj, spdata = spdata)
          estimates = unlist(aa$param)
          if(aa$convergence == "Successful") {
            cc = GeoCorrFct(semiv$centers, t=semiv$centert, corrmodel=corrmodel, model=model, distance=distance, 
                            param=c(aa$param, aa$fixed), radius=radius, n=n, covariance=TRUE, variogram=TRUE)$corr
            res = sum((cc - semiv$variograms)^2)
          } else {
            res = Inf
          }
          c(estimates, res = res)
        }
      estimates = results[, -ncol(results), drop=FALSE]
      rownames(estimates) <- NULL
      res = results[, ncol(results)]
      future::plan(sequential)
    }
  }

  if(space||bivariate) {
    indexmin = which.min(res)
    bestK = neighb[indexmin]
  }
  if(spacetime) {
    indexmin = which.min(res)
    bestKT = as.matrix(expand.grid(neighb, maxtime))[indexmin, ]
    bestK = as.numeric(bestKT[1])
    best_T = as.numeric(bestKT[2])
  }


# Aggiungere questo codice prima del return finale

sugg_neighb <- NULL
sugg_time <- NULL

if(space || bivariate) {
  # Criterio per i vicini spaziali
  if(bestK <= 20) {
    sugg_neighb <- bestK
  } else {

    sorted_neighb <- sort(neighb)
    median_neighb <- round(quantile(sorted_neighb,0.5))
    
    if(median_neighb < bestK) {
      sugg_neighb <- median_neighb
    } else {
      sugg_neighb <- bestK
    }
  }
}

if(spacetime) {
  # Criterio per i vicini spaziali nel caso spazio-temporale
  if(bestK <= 20) {
    sugg_neighb <- bestK
  } else {
    sorted_neighb <- sort(neighb)
    median_neighb <- round(median(quantile(sorted_neighb,0.5)))
    
    if(median_neighb < bestK) {
      sugg_neighb <- median_neighb
    } else {
      sugg_neighb <- bestK
    }
  }
  
  # Criterio per i tempi
  if(best_T <= 3) {
    sugg_time <- best_T
  } else {
    sorted_maxtime <- sort(maxtime)
    median_time <- round(quantile(sorted_maxtime,0.5))
    
    if(median_time < best_T) {
      sugg_time <- median_time
    } else {
      sugg_time <- best_T
    }
  }
}

  a = list(best_neighb = as.numeric(bestK), best_maxtime = best_T, res = unname(res), estimates = estimates, sugg_neighb = as.numeric(sugg_neighb), sugg_time = sugg_time)
  return(a)
}
