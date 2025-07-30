GeoSimcond <- function(estobj = NULL, data, coordx, coordy = NULL, coordz = NULL, coordt = NULL,
                       coordx_dyn = NULL, corrmodel, distance = "Eucl",
                       grid = FALSE, loc, maxdist = NULL, maxtime = NULL, method = "Cholesky",
                       model = "Gaussian", n = 1, nrep = 1, local = FALSE, L = 1000,
                       neighb = NULL, param, anisopars = NULL, radius = 1, sparse = FALSE,
                       time = NULL, copula = NULL, X = NULL, Xloc = NULL, Mloc = NULL,
                       parallel = TRUE, ncores = NULL) {

################################################################################# 
### internal function conditional simulation  of Gaussian RF
################################################################################# 
Gauss_cd <- function(data,corrmodel,nrep, method, L,
                     param,
                     coord_obs, loc, coordt_use, time, X, Xloc, Mloc, distance, radius,
                     local, neighb, maxdist, maxtime,
                     space, spacetime, bivariate, parallel, ncores) {
  #############################################
  # Unconditional simulation on combined coordinates (observed + prediction locations)
  #############################################
  coord_sim <- rbind(coord_obs, loc)
  time_sim <- c(coordt_use, time)
  X_sim <- rbind(X, Xloc)
  n_obs <- nrow(coord_obs)
  n_loc <- nrow(loc)

 cat("Performing", nrep, "unconditional simulations ...\n")
  if (method == "Cholesky") {
    sim_args <- list(coordx = coord_sim, coordt = time_sim, corrmodel = corrmodel,progress = FALSE,
                     X = X_sim, nrep = nrep, distance = distance, radius = radius)
    sim_args <- c(sim_args, list(param = param))
    sim_nc <- do.call(GeoSim, sim_args)
  }
  if (method == "TB" || method == "CE") {
    sim_args_approx <- list(coordx = coord_sim, coordt = time_sim, corrmodel = corrmodel,progress = FALSE,
                            method = method, L = L, parallel=parallel,
                            X = X_sim, nrep = nrep, distance = distance, radius = radius)
    sim_args_approx <- c(sim_args_approx, list(param = param))
    sim_nc <- do.call(GeoSimapprox, sim_args_approx)
  }


  #############################################
  # Kriging on prediction locations using observed data
  #############################################
  krig_sim_args <- list(coordx = coord_obs, coordt = coordt_use,
                        data = data, corrmodel = corrmodel,
                        loc = loc, X = X, Xloc = Xloc, Mloc = Mloc, time = time,
                        distance = distance, radius = radius)
  krig_sim_args <- c(krig_sim_args, list(param = param))

  if (local) {
    if (!is.null(neighb)) krig_sim_args$neighb <- neighb
    if (!is.null(maxdist)) krig_sim_args$maxdist <- maxdist
    if (!is.null(maxtime)) krig_sim_args$maxtime <- maxtime
    krig_sim_args$progress <- FALSE
    krig_sim <- do.call(GeoKrigloc, krig_sim_args)
  } else {
    krig_sim <- do.call(GeoKrig, krig_sim_args)
  }
  #############################################
  # Construct conditional simulations
  #############################################
  sim_cond <- vector("list", nrep) ### return a list

  if (space) {

    if (!parallel) {
      cat("Performing", nrep, "conditional simulations ...\n")

      # Setup progress bar
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = 1:nrep)

      for (i in seq_len(nrep)) {

         pb(sprintf("i=%d", i))
        sim_nc_obs_data <- sim_nc$data[[i]][1:n_obs]
        sim_nc_loc_data <- sim_nc$data[[i]][(n_obs + 1):(n_obs + n_loc)]

     #### to speed up here  one can compute the chol decomp only one time (only data are different)
        sim_nc_obs_args <- list(coordx = coord_obs, coordt = coordt_use, data = sim_nc_obs_data,
                               corrmodel = corrmodel, loc = loc, time = time,
                               distance = distance, radius = radius)
        sim_nc_obs_args <- c(sim_nc_obs_args, list(param = param))

        if (local) {
          if (!is.null(neighb)) sim_nc_obs_args$neighb <- neighb
          if (!is.null(maxdist)) sim_nc_obs_args$maxdist <- maxdist
          if (!is.null(maxtime)) sim_nc_obs_args$maxtime <- maxtime
          sim_nc_obs_args$progress <- FALSE
          sim_nc_obs <- do.call(GeoKrigloc, sim_nc_obs_args)
        } else {
          sim_nc_obs <- do.call(GeoKrig, sim_nc_obs_args)
        }
        sim_cond[[i]] <- krig_sim$pred + sim_nc_loc_data - sim_nc_obs$pred
      }

    } else { # parallel
      cat("Performing", nrep, "conditional simulations using", ncores, "cores ...\n")
      future::plan(future::multisession, workers = ncores)
      on.exit(future::plan(future::sequential), add = TRUE)
      # Setup progress bar
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = 1:nrep)
      sim_cond <- foreach::foreach(i = seq_len(nrep), .options.future = list(seed = TRUE,stdout = NA,  
                        conditions = character(0))) %dofuture% {
         pb(sprintf("i=%d", i))
        sim_nc_obs_data <- sim_nc$data[[i]][1:n_obs]
        sim_nc_loc_data <- sim_nc$data[[i]][(n_obs + 1):(n_obs + n_loc)]
        sim_nc_obs_args <- list(coordx = coord_obs, coordt = coordt_use, data = sim_nc_obs_data,
                               corrmodel = corrmodel, loc = loc, time = time,
                               distance = distance, radius = radius)
        sim_nc_obs_args <- c(sim_nc_obs_args, list(param = param))

        if (local) {
          if (!is.null(neighb)) sim_nc_obs_args$neighb <- neighb
          if (!is.null(maxdist)) sim_nc_obs_args$maxdist <- maxdist
          if (!is.null(maxtime)) sim_nc_obs_args$maxtime <- maxtime
          sim_nc_obs_args$progress <- FALSE
          sim_nc_obs <- do.call(GeoKrigloc, sim_nc_obs_args)
        } else {
          sim_nc_obs <- do.call(GeoKrig, sim_nc_obs_args)
        }

        krig_sim$pred + sim_nc_loc_data - sim_nc_obs$pred
      }
    }
  } # end space
  return(sim_cond)
}
###############################################################
###############################################################
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
    corrmodel <- estobj$corrmodel
    model <- estobj$model
    distance <- estobj$distance
    grid <- estobj$grid
    n <- estobj$n
    param <- append(estobj$param, estobj$fixed)
    radius <- estobj$radius
    copula <- estobj$copula
    anisopars <- estobj$anisopars
    X <- estobj$X
  }
  #### end check

  ### Additional checks
  if (is.null(CkModel(model))) stop("The name of the model is not correct\n")
  if (!is.character(corrmodel) || is.null(CkCorrModel(corrmodel)))
    stop("The name of the correlation model is wrong\n")
  if (!(method %in% c("Cholesky", "TB", "CE"))) {
    stop("The method of unconditional simulation is not correct\n")
  }
  if(local&&is.null(neighb)) stop("Fro local kriking you need to specify neighb \n")
  corrmodel <- gsub("[[:blank:]]", "", corrmodel)
  model <- gsub("[[:blank:]]", "", model)
  distance <- gsub("[[:blank:]]", "", distance)

#### checkin cores ######
    coremax <- parallel::detectCores()
  if (is.na(coremax) || coremax == 1) {
    parallel <- FALSE
  } else {
    if (is.null(parallel)) parallel <- TRUE
    if (is.null(ncores)) ncores <- max(1, coremax - 1)
  }

  # Determine model type ans setting (space spacetime bivariate)
  model_ck <- CkCorrModel(corrmodel)
  bivariate <- CheckBiv(model_ck)
  spacetime <- CheckST(model_ck)
  space <- !spacetime && !bivariate

  # Convert observed coordinates to matrix if needed
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
if(is.null(copula)){
  if(model=="Gaussian") {
                      res=Gauss_cd(data,corrmodel,nrep,method,L,
                              param,
                              coord_obs,loc,coordt_use,time,X,Xloc,Mloc,distance,radius,
                              local,neighb,maxdist,maxtime,
                              space,spacetime,bivariate,parallel,ncores)
                      }

  ################################################################################
  ############### monotone tranformations of one Gausssian RF ####################
  ################################################################################
  if(model=="LogGaussian")
  {    
                      mm=exp(param$mean)   ### what happened for mu=XB???
                      vv=param$sill
                      datanorm=(log(data/mm)+vv/2)/sqrt(vv) ##transformation in the gaussian scale 
                      param$mean=0;param$sill=1
                      res=Gauss_cd(datanorm,corrmodel,nrep,method,L,param,
                          coord_obs,loc,coordt_use,time,NULL,NULL,NULL,distance,radius,
                          local,neighb,maxdist,maxtime,
                          space,spacetime,bivariate,parallel,ncores)
                      res <- lapply(res, function(r) mm * exp(sqrt(vv) * r - vv/2)) # backtranformation

  }
   ################################################################################
  if(model=="Tukeyh") {
                      inverse_lamb=function(x,tail){
                      value = sqrt(VGAM::lambertW(tail*x*x)/tail);
                      return(sign(x)*value);
                      }
                      mm=param$mean
                      vv=param$sill
                      tail=param$tail
                      datanorm=inverse_lamb((data-mm)/sqrt(vv),tail)   ##transformation in the gaussian scale
                      param$mean=0;param$sill=1; param$tail=NULL
                      res=Gauss_cd(datanorm,corrmodel,nrep,method,L,param,
                          coord_obs,loc,coordt_use,time,NULL,NULL,NULL,distance,radius,
                          local,neighb,maxdist,maxtime,
                          space,spacetime,bivariate,parallel,ncores)
                      res <- lapply(res, function(r) mm+ sqrt(vv)*r*exp(tail*r^2/2)) # backtranformation
  }
  ################################################################################
  if(model=="Tukeyh2") {
     ################
    inverse_lamb=function(x,tail){
                      value = sqrt(VGAM::lambertW(tail*x*x)/tail);
                      return(sign(x)*value);
                      }
                mm=param$mean
                vv=param$sill
                tail1=param$tail1
                tail2=param$tail2
    inverse_lamb2=function(x,tail1,tail2){
                      aa=1:length(x)
                      sel1=I(x>=0)*aa
                      sel2=I(x<0)*aa
                      a1=inverse_lamb(x[sel1],tail1)  
                      b1=inverse_lamb(x[sel2],tail2)
                      sel1[sel1>0]<-a1
                      sel2[sel2>0]<-b1
                      return(c(sel1+sel2))
                      }
       ################
                     datanorm= inverse_lamb2((data-mm)/sqrt(vv),tail1,tail2)        
                     param$mean=0;param$sill=1; param$tail1=NULL;param$tail2=NULL
                     res=Gauss_cd(datanorm,corrmodel,nrep,method,L,param,
                          coord_obs,loc,coordt_use,time,NULL,NULL,NULL,distance,radius,
                          local,neighb,maxdist,maxtime,
                          space,spacetime,bivariate,parallel,ncores)

    ff<- function(r,mm,vv,tail1,tail2){
                      aa=1:length(r)
                      sel1=I(r>=0)*aa
                      sel2=I(r<0)*aa
                      res1=mm+ sqrt(vv)*r[sel1]*exp(tail1*r[sel1]^2/2)
                      res2=mm+ sqrt(vv)*r[sel2]*exp(tail2*r[sel2]^2/2)
                      sel1[sel1>0]<-res1
                      sel2[sel2>0]<-res2                    
                      return(c(sel1+sel2))}
                     
                     res <- lapply(res,ff,mm = mm, vv = vv, tail1 = tail1, tail2 = tail2) # backtranformation
  }
    ################################################################################
  if(model=="SinhAsinh") {
                      mm=exp(param$mean)  
                      vv=param$sill
                      skew=param$skew
                      tail=param$tail
                      datanorm=sinh(tail*asinh((data-mm)/sqrt(vv))-skew)      ##transformation in the gaussian scale
                      param$mean=0;param$sill=1; param$skew=NULL;param$tail=NULL
                      res=Gauss_cd(datanorm,corrmodel,nrep,method,L,param,
                          coord_obs,loc,coordt_use,time,NULL,NULL,NULL,distance,radius,
                          local,neighb,maxdist,maxtime,
                          space,spacetime,bivariate,parallel,ncores)
                      res <- lapply(res, function(r) mm+sqrt(vv)*sinh( (1/tail)*(asinh(r)+skew))) # backtranformation

  }
 ###############################################################################
 ############### tranformations of independent Gausssian RFs ###################
 ###############################################################################

if(model=="SkewGaussian") {
                     # mm=exp(param$mean)  
                     # vv=param$sill
                     # skew=param$skew
                    #...
                    }
if(model=="Gamma") {
                    #...

                    }
if(model=="Weibull") {
                    #...
                    }
if(model=="Poisson") {
                    #...
                    }


 }       ### end not copula models
 else{ 
###########################################################################
################## copula models ###########################################
############################################################################
if(copula=="Gaussian")
    {
         if(model=="Weibull")
         {
         mm=exp(param$mean)  
         ss=param$shape
         datanorm1=  pweibull(data,shape=ss,scale=mm/(gamma(1+1/ss)))
         datanorm=qnorm(datanorm1)
         param$mean=0;param$sill=1;param$shape=NULL;
         res=Gauss_cd(datanorm,corrmodel,nrep,method,L,param,
                          coord_obs,loc,coordt_use,time,NULL,NULL,NULL,distance,radius,
                          local,neighb,maxdist,maxtime,
                          space,spacetime,bivariate,parallel,ncores)
         res <- lapply(res, function(r) pnorm(r))
         res <- lapply(res, function(r) qweibull(r,shape=ss,scale=mm/(gamma(1+1/ss )))) # backtranformation
         }

        if(model=="LogGaussian")
        {                      
        mm=exp(param$mean)
        vv=param$sill
        datanorm1=plnorm(data, meanlog = mm - vv/2, sdlog = sqrt(vv), log.p = FALSE)
        datanorm=qnorm(datanorm1)
        param$mean=0;param$sill=1;
        res=Gauss_cd(datanorm,corrmodel,nrep,method,L,param,
                          coord_obs,loc,coordt_use,time,NULL,NULL,NULL,distance,radius,
                          local,neighb,maxdist,maxtime,
                          space,spacetime,bivariate,parallel,ncores)
        res <- lapply(res, function(r) pnorm(r))
        res <- lapply(res, function(r) qlnorm(r,meanlog=mm-vv/2, sdlog = sqrt(vv),log.p = FALSE)) # backtranformation
        #res <- lapply(res, function(r) pnorm((r-mm-vv/2)/sqrt(vv)))
        }



        if(model=="Gamma")
        {        
        mm=exp(param$mean)  
        ss=param$shape
        datanorm1= pgamma(data,shape=ss/2,rate=ss/(2*mm))
        datanorm=qnorm(datanorm1)
        param$mean=0;param$sill=1;param$shape=NULL;
        res=Gauss_cd(datanorm,corrmodel,nrep,method,L,param,
                          coord_obs,loc,coordt_use,time,NULL,NULL,NULL,distance,radius,
                          local,neighb,maxdist,maxtime,
                          space,spacetime,bivariate,parallel,ncores)
         res <- lapply(res, function(r) pnorm(r))
         res <- lapply(res, function(r) qgamma(r,shape=ss/2,rate=ss/(2*mm))) # backtranformation}
        }
         
        if(model=="Beta2")
        {
        mm=1/(1+exp(-param$mean))
        ss=param$shape
        pmin=param$min
        pmax=param$max
        datanorm1=pbeta((data-pmin)/(pmax-pmin),shape1=mm*ss,shape2=(1-mm)*ss)
        datanorm=qnorm(datanorm1)
        param$mean=0;param$sill=1;param$shape=NULL;param$min=NULL;param$max=NULL;
        res=Gauss_cd(datanorm,corrmodel,nrep,method,L,param,
                          coord_obs,loc,coordt_use,time,NULL,NULL,NULL,distance,radius,
                          local,neighb,maxdist,maxtime,
                          space,spacetime,bivariate,parallel,ncores)
        res <- lapply(res, function(r) pnorm(r))
        res <- lapply(res, function(r) pmin+(pmax-pmin)*qbeta(r,shape1=mm*ss,shape2=(1-mm)*ss))
        }

        if(model=="StudentT")
        {
        mm=param$mean
        vv=param$sill
        df=param$df
        datanorm1=pt((data-mm)/sqrt(vv),df=1/df)
        datanorm=qnorm(datanorm1)
        param$mean=0;param$sill=1;param$df=NULL;
        res=Gauss_cd(datanorm,corrmodel,nrep,method,L,param,
                          coord_obs,loc,coordt_use,time,NULL,NULL,NULL,distance,radius,
                          local,neighb,maxdist,maxtime,
                          space,spacetime,bivariate,parallel,ncores)
        res <- lapply(res, function(r) pnorm(r))
        res <- lapply(res, function(r) mm+sqrt(vv)*qt(r,df=1/df))
        }

        if(model=="Gaussian")
        {
        mm=param$mean
        vv=param$sill
        datanorm1=pnorm(data,mm,sqrt(vv))
        datanorm=qnorm(datanorm1)
        param$mean=0;param$sill=1;
        res=Gauss_cd(datanorm,corrmodel,nrep,method,L,param,
                          coord_obs,loc,coordt_use,time,NULL,NULL,NULL,distance,radius,
                          local,neighb,maxdist,maxtime,
                          space,spacetime,bivariate,parallel,ncores)
        res <- lapply(res, function(r) pnorm(r))
        res <- lapply(res, function(r) qnorm(r,mm,sqrt(vv)))
        }

         if(model=="Binomial")
        {
         mm <- param$mean; prob <- pnorm(mm)
         datanorm1 <- (pbinom(data - 1, size = n, prob = prob) + pbinom(data, size = n, prob = prob)) / 2 # Applicare correzione di continuità per la discrezione
         datanorm <- qnorm(datanorm1)
         param$mean <- 0;param$sill <- 1
         res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius,
                local, neighb, maxdist, maxtime,
                space, spacetime, bivariate, parallel, ncores)
        res <- lapply(res, function(r) pnorm(r))
        res <- lapply(res, function(r) qbinom(p = r, size = n, prob = prob))
        }
        if (model == "BinomialNeg")
         {
           mm <- param$mean; prob <- pnorm(mm);size=n
           datanorm1 <- (pnbinom(data - 1, size = size, prob = prob) + pnbinom(data,     size = size, prob = prob)) / 2 # Applicare correzione di continuità per la discrezione
           datanorm <- qnorm(datanorm1)
           param$mean <- 0; param$sill <- 1
           res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores)
           res <- lapply(res, function(r) pnorm(r))
           res <- lapply(res, function(r) qnbinom(p = r, size = size, prob = prob))
         }
         if (model == "Poisson")
         {
            mu <- exp(param$mean)
            datanorm1 <- (ppois(data - 1, lambda = mu) + ppois(data,     lambda = mu)) / 2    # Trasformazione: discreto → normale (con correzione di continuità)
            datanorm <- qnorm(datanorm1)
            # Standardizzazione per copula gaussiana
            param$mean <- 0; param$sill <- 1
            # Simulazione condizionata
            res <- Gauss_cd(datanorm, corrmodel, nrep, method, L, param,
                    coord_obs, loc, coordt_use, time, NULL, NULL, NULL, distance, radius,
                    local, neighb, maxdist, maxtime,
                    space, spacetime, bivariate, parallel, ncores)
            # Back-transform: N(0,1) → Poisson(μ)
            res <- lapply(res, function(r) pnorm(r))
            res <- lapply(res, function(r) qpois(r, lambda = mu))
        }


    }
}

############################################################################
############################################################################
############################################################################
############################
GeoSimcond_out = list(   
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
    #copula=copula,
    numcoord = nrow(coord_obs),
    numloc= nrow(loc),
    numtime = length(coordt),
    numt = length(time),
    maxdist=maxdist,
    maxtime=maxtime,
    model=model,
    n=n,
    param = param,
    condsim=res,            #####result!
    radius=radius,
    spacetime = spacetime,
    time=time
  )
  structure(c(GeoSimcond_out, call = call), class = c("GeoSimcond"))
}


