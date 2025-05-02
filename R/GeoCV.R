GeoCV <- function(fit, K = 100, estimation = TRUE, 
                  optimizer = NULL, lower = NULL, upper = NULL,
                  n.fold = 0.05, local = FALSE, neighb = NULL, maxdist = NULL, 
                  maxtime = NULL, sparse = FALSE, type_krig = "Simple", 
                  which = 1, parallel = FALSE, ncores = NULL) {

  if (!inherits(fit, "GeoFit")) stop("fit must be an object of class 'GeoFit'")
  if (n.fold > 0.99 || n.fold < 0.01) stop("n.fold must be between 0.01 and 0.99")
  if (!is.logical(parallel)) stop("parallel must be TRUE or FALSE")
  print("Cross-validation kriging can be time consuming ...")

  if (is.null(fit$X)) { X = Xloc = NULL; tempX = NULL }
  mae = rmse = lscore = crps = mad = brie = NULL
  space_dyn = FALSE
  if (is.null(optimizer)) { optimizer = fit$optimizer; lower = fit$lower; upper = fit$upper }

  if (is.list(fit$data)) space_dyn = TRUE
  spacetime <- CheckST(CkCorrModel(fit$corrmodel))
  bivariate <- CheckBiv(CkCorrModel(fit$corrmodel))
  K = round(K)
  if (K < 2) stop("K must be greater or equal to 2")
  if (K > 1000) stop("K is too large")
  if (bivariate) {
    if (!(which == 1 || which == 2)) stop("which must be 1 or 2")
  }
  if (local) if (is.null(maxdist) && is.null(neighb)) stop("maxdist or neighb are required for local kriging")

  cat("Starting iteration from 1 to", K, " ...\n")
  space = !spacetime && !bivariate


  model1 = fit$model
  if (fit$missp) {  # Misspecification
    if (fit$model == "StudentT") model1 = "Gaussian_misp_StudentT"
    if (fit$model == "Poisson") model1 = "Gaussian_misp_Poisson"
    if (fit$model == "PoissonZIP") model1 = "Gaussian_misp_PoissonZIP"
    if (fit$model == "SkewStudentT") model1 = "Gaussian_misp_SkewStudenT"
    if (fit$model == "Tukeygh") model1 = "Gaussian_misp_Tukeygh"
  }


    coremax <- parallel::detectCores()
    if (is.na(coremax) || coremax == 1) {parallel <- FALSE}

  ########################################################################################################################
  ########### spatial case ###############################################################################################
  ########################################################################################################################
  if (space) {
    N <- length(fit$data)
    coords <- cbind(fit$coordx, fit$coordy, fit$coordz)
    data <- fit$data
    if (length(fit$fixed$mean) > 1) tempM <- fit$fixed$mean
    if (!is.null(fit$X)) tempX <- fit$X

    # samples
    sel_list <- replicate(K, sample(1:N, round(N * (1 - n.fold))), simplify = FALSE)

    ###################################
    ####### non parallel version 
    ###################################
    if (!parallel) {

      rmse <- crps <- mae <- mad <- lscore <- brie <- double(K)

      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = 1:K)

      cat("Performing", K, "cross-validations...\n")

       X <- Xloc <-Mloc<- NULL
      i = 1
      while (i <= K) {
        sel_data <- sel_list[[i]]   # random sample

        if (!is.null(fit$X)) {
          X <- tempX[sel_data, ]
          Xloc <- tempX[-sel_data, ]
        } 

        if (length(fit$fixed$mean) > 1) {
          fit$fixed$mean <- tempM[sel_data]
          Mloc <- tempM[-sel_data]
        }

        data_to_pred <- fit$data[-sel_data]
        data_to_est <- fit$data[sel_data]
        coords_est <- coords[sel_data, ]
        coords_to_pred <- coords[-sel_data, ]

        param <- c(fit$param, fit$fixed)
        if (estimation) {
          fit_s <- GeoFit(
            data = data_to_est, coordx = coords_est, corrmodel = fit$corrmodel,
            X = X, likelihood = fit$likelihood, type = fit$type, grid = fit$grid,
            copula = fit$copula, anisopars = fit$anisopars, est.aniso = fit$est.aniso,
            model = model1, radius = fit$radius, n = fit$n,
            maxdist = fit$maxdist, neighb = fit$neighb, distance = fit$distance,
            optimizer = optimizer, lower = lower, upper = upper,
            start = fit$param, fixed = fit$fixed
          )
          if (!is.null(fit$anisopars)) {
            fit$param$angle <- NULL; fit$param$ratio <- NULL
            fit$fixed$angle <- NULL; fit$fixed$ratio <- NULL
          }
          #param <- append(fit_s$param, fit_s$fixed)
          param <- c(fit_s$param, fit_s$fixed)
        } else {
          fit_s <- fit
        }

        if (!local) {
          pr <- GeoKrig(fit_s, loc = coords_to_pred, mse = TRUE,
                        param = param, Xloc = Xloc, Mloc = Mloc)
        } else {
          pr <- GeoKrigloc(fit_s, loc = coords_to_pred, mse = TRUE,
                           param = param, Xloc = Xloc, Mloc = Mloc,
                           neighb = neighb, maxdist = maxdist)
        }

        pp <- GeoScores(data_to_pred, pred = pr$pred, mse = pr$mse,
                        score = c("brie", "crps", "lscore", "pe"))

        rmse[i] <- pp$rmse
        mae[i] <- pp$mae
        mad[i] <- pp$mad
        lscore[i] <- pp$lscore
        brie[i] <- pp$brie
        crps[i] <- pp$crps

        pb(sprintf("i=%g", i))
        i <- i + 1
      }
    } # end non parallel

    ######################################################################
    #######  parallel version 
    ######################################################################
    if (parallel) {

      if (is.null(ncores)) {n.cores <- coremax - 1} 
      else { if (!is.numeric(ncores) || ncores > coremax || ncores < 1)
          stop("number of cores not valid\n")
          n.cores <- ncores
      }
    
      future::plan(multisession, workers = n.cores)
      on.exit(future::plan(sequential), add = TRUE)
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = 1:K)
      
      cat("Performing", K, "cross-validations using", n.cores, "cores...\n")
      
################################################################
############ starting foreach ##################################
################################################################
      results <- foreach::foreach(i = 1:K, .combine = 'rbind',
                                   .options.future = list(seed = TRUE)) %dofuture% {
     X <- Xloc <-Mloc<- NULL
        sel_data <- sel_list[[i]]  # random sample


  if (!is.null(fit$X)) {
          X <- tempX[sel_data, ]
          Xloc <- tempX[-sel_data, ]
        } 

        if (length(fit$fixed$mean) > 1) {
          fit$fixed$mean <- tempM[sel_data]
          Mloc <- tempM[-sel_data]
        }
        data_to_pred <- data[-sel_data]
        data_to_est <- data[sel_data]
        coords_est <- coords[sel_data, ]
        coords_to_pred <- coords[-sel_data, ]

        param <- c(fit$param, fit$fixed)

        if (estimation) {
          fit_s <- GeoFit(
            data = data_to_est, coordx = coords_est, corrmodel = fit$corrmodel,
            X = X, likelihood = fit$likelihood, type = fit$type,
            grid = fit$grid, copula = fit$copula, anisopars = fit$anisopars,
            est.aniso = fit$est.aniso, model = model1, radius = fit$radius,
            n = fit$n, local = fit$local, GPU = fit$GPU,
            maxdist = fit$maxdist, neighb = fit$neighb, distance = fit$distance,
            optimizer = optimizer, lower = lower, upper = upper,
            start = fit$param, fixed = fit$fixed
          )

          if (!is.null(fit$anisopars)) {
            fit_s$param$angle <- NULL; fit_s$param$ratio <- NULL
            fit_s$fixed$angle <- NULL; fit_s$fixed$ratio <- NULL
          }
        #param <- append(fit_s$param, fit_s$fixed)
        param <- c(fit_s$param, fit_s$fixed)
        } 
        else {fit_s <- fit}

        if (!local) {
          pr <- GeoKrig(fit_s, loc = coords_to_pred, mse = TRUE,
                        param = param, Xloc = Xloc, Mloc = Mloc)
        } else {
          pr <- GeoKrigloc(fit_s, loc = coords_to_pred, mse = TRUE,
                           param = param, Xloc = Xloc, Mloc = Mloc,
                           neighb = neighb, maxdist = maxdist)
        }

        pp <- GeoScores(data_to_pred, pred = pr$pred, mse = pr$mse,
                        score = c("brie", "crps", "lscore", "pe"))

        pb(sprintf("i=%g", i))

        return(c(pp$rmse, pp$mae, pp$mad, pp$lscore, pp$brie, pp$crps))
      }

      rmse <- unname(results[, 1])
      mae <- unname(results[, 2])
      mad <- unname(results[, 3])
      lscore <- unname(results[, 4])
      brie <- unname(results[, 5])
      crps <- unname(results[, 6])

      #future::plan(sequential)
      #on.exit(future::plan(sequential), add = TRUE)
    }
  }









############################################################
############################################################
############################################################
########### spatio temporal case ###########################
############################################################
############################################################
############################################################

if (spacetime) {

  
  coords = cbind(fit$coordx, fit$coordy, fit$coordz)
  ns = fit$ns; coordt = fit$coordt
  T = length(coordt); NT = sum(ns)

  if (!is.null(fit$X)) X <- fit$X
  else  X = NULL

  NS = cumsum(ns); NS = c(c(0, NS)[-(length(ns) + 1)], NT)

  if (is.list(fit$data)) {datos = do.call(c, args = c(fit$data))} 
  else {datos = fit$data}
########
  if (!space_dyn) {
    data_tot =  NULL
    for (k in 1:T) {
      data_tot = rbind(data_tot, cbind(rep(coordt[k], ns[k]), coords, datos[k, ]))  
    }
    data_tot = cbind(data_tot, X) 
  }

########
  if (space_dyn) {
  ct = NULL
  for (k in 1:T) {ct = c(ct, rep(coordt[k], ns[k]))}
  
  if (!is.null(X)) {
    data_tot = cbind(ct, fit$coordx, fit$coordy, fit$coordz, datos, do.call(rbind, args = c(X)))
  } else {
    data_tot = cbind(ct, fit$coordx, fit$coordy, fit$coordz, datos)
  }
}
########
MM=length(fit$fixed$mean) > 1
if(MM&!is.null(X)) stop("covariates and fixed varying mean are not compatible")

if (MM) {
  if (!space_dyn) tempM=fit$fixed$mean 
  else            tempM= do.call(rbind, args = c(fit$fixed$mean))
  data_tot=cbind(data_tot,tempM)
    }  
########

  sel_data = sample(1:NT, round(NT * (1 - n.fold))) 
  folds = split(sample(1:NT), rep(1:K, length.out = NT))


###################### not parallel case ############################

if (!parallel) {

  progressr::handlers(global = TRUE)
  progressr::handlers("txtprogressbar")
  pb <- progressr::progressor(along = 1:K)

  rmse = crps = mae = mad = lscore = brie = double(K)

   cat("Performing", K, "cross-validations...\n")

  for (i in 1:K) {

    test_idx = folds[[i]]
    train_idx = setdiff(1:NT, test_idx)

    data_sel = data_tot[train_idx, ]
    data_to_pred = data_tot[test_idx, ]

    if (is.vector(data_to_pred)) data_to_pred = matrix(data_to_pred, nrow = 1)
    if (is.vector(data_sel)) data_sel = matrix(data_sel, nrow = 1)

    data_sel_ord = as.matrix(data_sel[order(data_sel[, 1]), ])
    data_to_pred_ord = as.matrix(data_to_pred[order(data_to_pred[, 1]), ])

    DD = ncol(data_sel_ord)
    coordx_dynnew = Xnew = Xnew_loc = datanew = coordx_dynnew_loc = list()

    utt = unique(data_sel_ord[, 1])
    utt_1 = unique(data_to_pred_ord[, 1])

    
    Xnew_loc=Mnew_loc=NULL
    for (k in 1:length(utt_1)) {
      ll = data_to_pred_ord[data_to_pred_ord[, 1] == utt_1[k], , drop = FALSE]
   
      if (ncol(ll) >= 3) {
        if (!is.null(fit$coordz))  # 3d
        {
              coordx_dynnew_loc[[k]] = matrix(ll[, 2:4], ncol = 3)
              if(!MM){if (!is.null(X)) Xnew_loc[[k]] = matrix(ll[, 6:DD], ncol = DD - 5)}
                else   {Mnew_loc[[k]] = ll[, 6]}
        } 
      else {                     #2d
          coordx_dynnew_loc[[k]] = matrix(ll[, 2:3], ncol = 2)
           if(!MM){ if (!is.null(X)) Xnew_loc[[k]] = matrix(ll[, 5:DD], ncol = DD - 4) }
           else   {Mnew_loc[[k]] = ll[, 6]}
        }
      } 
    }


 
     Xnew = NULL
    for (k in 1:length(utt)) {
      ss = data_sel_ord[data_sel_ord[, 1] == utt[k], , drop = FALSE]
      if (!is.null(fit$coordz)) {
        coordx_dynnew[[k]] = matrix(ss[, 2:4], ncol = 3)
        datanew[[k]] = as.vector(ss[, 5])
        if (!is.null(X)) Xnew[[k]] = matrix(ss[, 6:DD], ncol = DD - 5)
      } else {
        coordx_dynnew[[k]] = matrix(ss[, 2:3], ncol = 2)
        datanew[[k]] = as.vector(ss[, 4])
        if (!is.null(X)) Xnew[[k]] = matrix(ss[, 5:DD], ncol = DD - 4)
      }
    }

    param = c(fit$param, fit$fixed)

    if (estimation) {
      fit_s = GeoFit(data = datanew, coordx_dyn = coordx_dynnew, coordt = utt,
                     corrmodel = fit$corrmodel, X = Xnew, likelihood = fit$likelihood, 
                     type = fit$type, grid = fit$grid, copula = fit$copula, 
                     anisopars = fit$anisopars, est.aniso = fit$est.aniso,
                     model = model1, radius = fit$radius, n = fit$n,
                     maxdist = fit$maxdist, neighb = fit$neighb, maxtime = fit$maxtime, 
                     distance = fit$distance, optimizer = optimizer, 
                     lower = lower, upper = upper,
                     start = fit$param, fixed = fit$fixed)

      if (!is.null(fit$anisopars)) {
        fit_s$param$angle = NULL
        fit_s$param$ratio = NULL
        fit_s$fixed$angle = NULL
        fit_s$fixed$ratio = NULL
      }

      #param = append(fit_s$param, fit_s$fixed)
      param = c(fit_s$param, fit_s$fixed)
    }

    pr_st = pr_mse = list()
    for (j in 1:length(utt_1)) {
      if (is.null(Xnew)) Xnew_loc[j] = list(NULL)
      if (!MM) Mnew_loc[j] = list(NULL)

      if (!local) {
        pr = GeoKrig(fit_s, loc = coordx_dynnew_loc[[j]], 
                     param = param, time = utt_1[j], mse = TRUE,  Mloc=Mnew_loc[[j]],
                     Xloc = Xnew_loc[[j]])
      } else {
        pr = GeoKrigloc(fit_s, loc = coordx_dynnew_loc[[j]], 
                        param = param, time = utt_1[j], mse = TRUE,  Mloc=Mnew_loc[[j]],
                        neighb = neighb, maxdist = maxdist, maxtime = maxtime, 
                        Xloc = Xnew_loc[[j]])
      }
      pr_st[[j]] = pr$pred
      pr_mse[[j]] = pr$mse
    }

    pp = GeoScores(c(data_to_pred[, 4]), pred = as.numeric(unlist(pr_st)), 
                   mse = as.numeric(unlist(pr_mse)),
                   score = c("brie", "crps", "lscore", "pe"))

    rmse[i] = unname(pp$rmse)
    mae[i] = unname(pp$mae)
    mad[i] = unname(pp$mad)
    lscore[i] = unname(pp$lscore)
    brie[i] = unname(pp$brie)
    crps[i] = unname(pp$crps)

    rm(pr, pp, data_sel, data_to_pred)
    gc(full = TRUE)

    
  }
}### End Not Parallel

######################  parallel case ############################
if (parallel) {
  
  if (is.null(ncores)) {n.cores <- coremax - 1} 
  else { 
    if (!is.numeric(ncores) || ncores > coremax || ncores < 1)
      stop("number of cores not valid\n")
    n.cores <- ncores
  }
  

      future::plan(multisession, workers = n.cores)
      on.exit(future::plan(sequential), add = TRUE)
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = 1:K)
      
      cat("Performing", K, "cross-validations using", n.cores, "cores...\n")
  
      results <- foreach(i = 1:K,
                   .combine = 'rbind',  
                   .options.future = list(seed = TRUE)) %dofuture% {

    test_idx = folds[[i]]
    train_idx = setdiff(1:NT, test_idx)

    data_sel = data_tot[train_idx, ]
    data_to_pred = data_tot[test_idx, ]

    if (is.vector(data_to_pred)) data_to_pred = matrix(data_to_pred, nrow = 1)
    if (is.vector(data_sel)) data_sel = matrix(data_sel, nrow = 1)

    data_sel_ord = as.matrix(data_sel[order(data_sel[, 1]), ])
    data_to_pred_ord = as.matrix(data_to_pred[order(data_to_pred[, 1]), ])

    DD = ncol(data_sel_ord)
    coordx_dynnew = Xnew = Xnew_loc = datanew = coordx_dynnew_loc = list()

    utt = unique(data_sel_ord[, 1])
    utt_1 = unique(data_to_pred_ord[, 1])

  Xnew_loc=Mnew_loc=NULL
    for (k in 1:length(utt_1)) {
      ll = data_to_pred_ord[data_to_pred_ord[, 1] == utt_1[k], , drop = FALSE]
   
      if (ncol(ll) >= 3) {
        if (!is.null(fit$coordz))  # 3d
        {
              coordx_dynnew_loc[[k]] = matrix(ll[, 2:4], ncol = 3)
              if(!MM){if (!is.null(X)) Xnew_loc[[k]] = matrix(ll[, 6:DD], ncol = DD - 5)}
                else   {Mnew_loc[[k]] = ll[, 6]}
        } 
      else {                     #2d
          coordx_dynnew_loc[[k]] = matrix(ll[, 2:3], ncol = 2)
           if(!MM){ if (!is.null(X)) Xnew_loc[[k]] = matrix(ll[, 5:DD], ncol = DD - 4) }
           else   {Mnew_loc[[k]] = ll[, 6]}
        }
      } 
    }
 
     Xnew = NULL
    for (k in 1:length(utt)) {
      ss = data_sel_ord[data_sel_ord[, 1] == utt[k], , drop = FALSE]
      if (!is.null(fit$coordz)) {
        coordx_dynnew[[k]] = matrix(ss[, 2:4], ncol = 3)
        datanew[[k]] = as.vector(ss[, 5])
        if (!is.null(X)) Xnew[[k]] = matrix(ss[, 6:DD], ncol = DD - 5)
      } else {
        coordx_dynnew[[k]] = matrix(ss[, 2:3], ncol = 2)
        datanew[[k]] = as.vector(ss[, 4])
        if (!is.null(X)) Xnew[[k]] = matrix(ss[, 5:DD], ncol = DD - 4)
      }
    }

    #param = append(fit$param, fit$fixed)
    param = c(fit$param, fit$fixed)

    if (estimation) {
      fit_s = GeoFit(data = datanew, coordx_dyn = coordx_dynnew, coordt = utt,
                     corrmodel = fit$corrmodel, X = Xnew, likelihood = fit$likelihood, 
                     type = fit$type, grid = fit$grid, copula = fit$copula, 
                     anisopars = fit$anisopars, est.aniso = fit$est.aniso,
                     model = model1, radius = fit$radius, n = fit$n,
                     maxdist = fit$maxdist, neighb = fit$neighb, maxtime = fit$maxtime, 
                     distance = fit$distance, optimizer = optimizer, 
                     lower = lower, upper = upper,
                     start = fit$param, fixed = fit$fixed)

      if (!is.null(fit$anisopars)) {
        fit_s$param$angle = NULL
        fit_s$param$ratio = NULL
        fit_s$fixed$angle = NULL
        fit_s$fixed$ratio = NULL
      }

      #param = append(fit_s$param, fit_s$fixed)
      param = c(fit_s$param, fit_s$fixed)
    }

    pr_st = pr_mse = list()
    for (j in 1:length(utt_1)) {
      if (is.null(Xnew)) Xnew_loc[j] = list(NULL)

      if (!local) {
        pr = GeoKrig(fit_s, loc = coordx_dynnew_loc[[j]], 
                     param = param, time = utt_1[j], mse = TRUE, Mloc=Mnew_loc[[j]],
                     Xloc = Xnew_loc[[j]])
      } else {
        pr = GeoKrigloc(fit_s, loc = coordx_dynnew_loc[[j]], 
                        param = param, time = utt_1[j], mse = TRUE,  Mloc=Mnew_loc[[j]],
                        neighb = neighb, maxdist = maxdist, maxtime = maxtime, 
                        Xloc = Xnew_loc[[j]])
      }

      pr_st[[j]] = pr$pred
      pr_mse[[j]] = pr$mse
    }


  pp = GeoScores(c(data_to_pred[, 4]), pred = as.numeric(unlist(pr_st)), 
                   mse = as.numeric(unlist(pr_mse)),
                   score = c("brie", "crps", "lscore", "pe"))
    
   
    pb(sprintf("i=%g", i))
    # Restituisce un vettore numerico invece di una lista
    c(pp$rmse, pp$mae, pp$mad, pp$lscore, pp$brie, pp$crps)

  }
  
  # Estrazione dei risultati dalla matrice
  rmse = unname(results[, 1])
  mae = unname(results[, 2])
  mad = unname(results[, 3])
  lscore = unname(results[, 4])
  brie = unname(results[, 5])
  crps = unname(results[, 6])
}

} ## end spacetime









############################################################
############################################################
############################################################
########### spatial bivariate case #########################
############################################################
############################################################
############################################################
if (bivariate) {

  fixmeans <- length(fit$fixed$mean_1) > 1 && length(fit$fixed$mean_2) > 1
  if (fixmeans) {
    tempM1 <- fit$fixed$mean_1
    tempM2 <- fit$fixed$mean_2
  }

  ns <- fit$ns
  if (space_dyn) {
    data1 <- fit$data[[1]]; data2 <- fit$data[[2]]
    coords1 <- fit$coordx_dyn[[1]]; coords2 <- fit$coordx_dyn[[2]]
  } else {
    data1 <- fit$data[1, ]; data2 <- fit$data[2, ]
    coords <- cbind(fit$coordx, fit$coordy, fit$coordz)
    coords1 <- coords; coords2 <- coords
  }


  X1 <- X2 <- NULL
  if (!is.null(fit$X)) {
    if (!space_dyn && !is.list(fit$X)) {
      X1 <- fit$X[1:ns[1], ]; X2 <- fit$X[(ns[1]+1):(ns[1]+ns[2]), ]
    } else if (space_dyn && is.list(fit$X)) {
      X1 <- fit$X[[1]]; X2 <- fit$X[[2]]
    }
  }

# ==================== VERSIONE NON PARALLELA ====================
  if (!parallel) {     
    metrics <- matrix(0, nrow = K, ncol = 6)
    colnames(metrics) <- c("rmse", "mae", "mad", "lscore", "brie", "crps")
    progressr::handlers(global = TRUE)
    progressr::handlers("txtprogressbar")
    pb <- progressr::progressor(along = 1:K)
     cat("Performing", K, "cross-validations...\n")
    for (i in 1:K) {
    
      if (which == 1) {
        sel_data <- sample(1:ns[1], round(ns[1] * (1 - n.fold)))
        data_to_est <- list(data1[sel_data], data2)
        data_to_pred <- data1[-sel_data]
        coords_est <- list(coords1[sel_data, ], coords2)
        coords_pred <- coords1[-sel_data, ]
        
        X_use <- if (!is.null(X1)) rbind(X1[sel_data, ], X2) else NULL
        Xloc_use <- if (!is.null(X1)) rbind(X1[-sel_data, ], X2) else NULL
        
        if (fixmeans) {
          current_fixed <- fit$fixed
          current_fixed$mean_1 <- tempM1[sel_data]
          Mloc <- list(tempM1[-sel_data], tempM2)
        } else {
          current_fixed <- fit$fixed
          Mloc <- NULL
        }
      } 
      else  ## case which 2
      {
        sel_data <- sample(1:ns[2], round(ns[2] * (1 - n.fold)))
        data_to_est <- list(data1, data2[sel_data])
        data_to_pred <- data2[-sel_data]
        coords_est <- list(coords1, coords2[sel_data, ])
        coords_pred <- coords2[-sel_data, ]
        
        X_use <- if (!is.null(X2)) rbind(X1, X2[sel_data, ]) else NULL
        Xloc_use <- if (!is.null(X2)) rbind(X1, X2[-sel_data, ]) else NULL
        
        if (fixmeans) {
          current_fixed <- fit$fixed
          current_fixed$mean_2 <- tempM2[sel_data]
          Mloc <- list(tempM1, tempM2[-sel_data])
        } else {
          current_fixed <- fit$fixed
          Mloc <- NULL
        }
      }

      # Stima parametri
      param <- if (estimation) {
        fit_s <- GeoFit(
          data = data_to_est, coordx = NULL, coordx_dyn = coords_est,
          corrmodel = fit$corrmodel, X = X_use,
          likelihood = fit$likelihood, type = fit$type, grid = fit$grid,
          model = "Gaussian", radius = fit$radius, n = fit$n,
          copula = fit$copula, 
          maxdist = fit$maxdist, neighb = fit$neighb, distance = fit$distance,
          optimizer = optimizer, lower = lower, upper = upper,
          start = fit$param, fixed = current_fixed
        )
        if (!is.null(fit$anisopars)) {
          fit_s$param$angle <- fit_s$param$ratio <- NULL
          fit_s$fixed$angle <- fit_s$fixed$ratio <- NULL
        }
        append(fit_s$param, fit_s$fixed)
      } else {
        append(fit$param, current_fixed)
      }

      pr <- if (!local) {
        GeoKrig(fit_s,
           coordx = NULL, 
          loc = coords_pred, mse = TRUE,
          param = param, Xloc = Xloc_use, Mloc = Mloc, which = which)
      } else {
        GeoKrigloc(fit_s,
           coordx = NULL, 
          loc = coords_pred, mse = TRUE,neighb = neighb, maxdist = maxdist,
          param = param, Xloc = Xloc_use, Mloc = Mloc, which = which)
      }

      # Calcolo metriche
      scores <- GeoScores(
    data_to_pred, pred = pr$pred, mse = pr$mse,score = c("brie", "crps", "lscore", "pe"))
      
      metrics[i, ] <- c(scores$rmse, scores$mae, scores$mad, 
                       scores$lscore, scores$brie, scores$crps)
      pb(sprintf("i=%g", i))
    }
    
    # Assegnazione risultati
    rmse <- metrics[, "rmse"]
    mae <- metrics[, "mae"]
    mad <- metrics[, "mad"]
    lscore <- metrics[, "lscore"]
    brie <- metrics[, "brie"]
    crps <- metrics[, "crps"]
    
  } 

  # ==================== VERSIONE PARALLELA ====================
  if(parallel) {   

  
  if (is.null(ncores)) {n.cores <- coremax - 1} 
  else { 
    if (!is.numeric(ncores) || ncores > coremax || ncores < 1)
      stop("number of cores not valid\n")
    n.cores <- ncores
  }

    samples <- lapply(1:K, function(i) {
                      if (which == 1) { sample(1:ns[1], round(ns[1] * (1 - n.fold)))} 
                      else {sample(1:ns[2], round(ns[2] * (1 - n.fold)))}
                       })
    future::plan(future::multisession, workers = n.cores)
    on.exit(future::plan(sequential), add = TRUE)

    progressr::handlers(global = TRUE)
    progressr::handlers("txtprogressbar")
    pb <- progressr::progressor(along = 1:K)
       cat("Performing", K, "cross-validations using", n.cores, "cores...\n")
    
    ############################################ 
    ########### starting foreach ############### 
    ############################################
    results <- foreach::foreach(i = 1:K,.combine = rbind,.options.future = list(seed =TRUE)  ) %dofuture% {

      sel_data <- samples[[i]] ##### sample!
      
      if (which == 1) {
        data_to_est <- list(data1[sel_data], data2)
        data_to_pred <- data1[-sel_data]
        coords_est <- list(coords1[sel_data, ], coords2)
        coords_pred <- coords1[-sel_data, ]
        X_use <- if (!is.null(X1)) rbind(X1[sel_data, ], X2) else NULL
        Xloc_use <- if (!is.null(X1)) rbind(X1[-sel_data, ], X2) else NULL
        
        if (fixmeans) {
          current_fixed <- fit$fixed
          current_fixed$mean_1 <- tempM1[sel_data]
          Mloc <- list(tempM1[-sel_data], tempM2)
        } else {
          current_fixed <- fit$fixed
          Mloc <- NULL
        }
      } else {
        data_to_est <- list(data1, data2[sel_data])
        data_to_pred <- data2[-sel_data]
        coords_est <- list(coords1, coords2[sel_data, ])
        coords_pred <- coords2[-sel_data, ]
        
        X_use <- if (!is.null(X2)) rbind(X1, X2[sel_data, ]) else NULL
        Xloc_use <- if (!is.null(X2)) rbind(X1, X2[-sel_data, ]) else NULL
        
        if (fixmeans) {
          current_fixed <- fit$fixed
          current_fixed$mean_2 <- tempM2[sel_data]
          Mloc <- list(tempM1, tempM2[-sel_data])
        } else {
          current_fixed <- fit$fixed
          Mloc <- NULL
        }
      }

      param <- if (estimation) {
        fit_s <- GeoFit(
          data = data_to_est, coordx = NULL, coordx_dyn = coords_est,
          corrmodel = fit$corrmodel, X = X_use,
          likelihood = fit$likelihood, type = fit$type, grid = fit$grid,
          model = "Gaussian", radius = fit$radius, n = fit$n,
          copula = fit$copula, 
          maxdist = fit$maxdist, neighb = fit$neighb, distance = fit$distance,
          optimizer = optimizer, lower = lower, upper = upper,
          start = fit$param, fixed = current_fixed
        )
        if (!is.null(fit$anisopars)) {
          fit_s$param$angle <- fit_s$param$ratio <- NULL
          fit_s$fixed$angle <- fit_s$fixed$ratio <- NULL
        }
        append(fit_s$param, fit_s$fixed)
      } else {
        append(fit$param, current_fixed)
      }

      pr <- if (!local) {
        GeoKrig(fit_s,
          loc = coords_pred,mse = TRUE,
          param = param, Xloc = Xloc_use, Mloc = Mloc, which = which)
      } else {
        GeoKrigloc(fit_s,
          loc = coords_pred,mse = TRUE,
           neighb = neighb, maxdist = maxdist,
          param = param, Xloc = Xloc_use, Mloc = Mloc, which = which)
      }
      
      scores <- GeoScores(
        data_to_pred, pred = pr$pred, mse = pr$mse,
        score = c("brie", "crps", "lscore", "pe")
      )
      pb(sprintf("i=%g", i))
      c(scores$rmse, scores$mae, scores$mad, 
        scores$lscore, scores$brie, scores$crps)
    }
  

    rmse <- unname(results[, 1])
    mae <- unname(results[, 2])
    mad <- unname(results[, 3])
    lscore <- unname(results[, 4])
    brie <- unname(results[, 5])
    crps <- unname(results[, 6])
  }
}  # end bivariate

return(list(rmse=rmse,mae=mae,mad=mad,brie=brie,crps=crps,lscore=lscore))
}
