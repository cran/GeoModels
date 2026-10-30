####################################################
### File name: WeightedLeastSquare.r
####################################################


### Procedures are in alphabetical order.

print.GeoWLS <- function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    if(x$weighted)
      method <- 'Weighted Least Squares'
    else
      method <- 'Least Squares'

    cat('\n##############################################################')
    cat('\nResults:', method,'Fitting of Gaussian Random Fields.\n')
    cat('\nCovariance model:', x$corrmodel, '\n')
    cat('Number of spatial coordinates:', x$numcoord, '\n')
    cat('Number of dependent temporal realisations:', x$numtime, '\n')
    cat('Number of estimated parameters:', length(x$param), '\n')
    cat('The value of the', method, 'at the minimum:',-x$wls,'\n')
    cat('Number of spatial bins', length(x$bins),'\n')
    cat('Number of temporal bins', length(x$bint),'\n')
    #cat('Min and max spatial distances:', x$srange,'\n')
    #if(length(x$coordt)>1) cat('Min and max temporal interval:', x$trange,'\n')

    cat('\nEstimated parameters:\n')
    print.default(x$param, digits = digits, print.gap = 2,
                  quote = FALSE)

    cat('\n##############################################################\n')
    invisible(x)
  }



WlsStart <- function(coordx, coordy, coordz,coordt, coordx_dyn, corrmodel, data, distance, fcall, fixed, grid,
                    likelihood, maxdist, neighb,maxtime, model, n, param, 
                    parscale, paramrange, radius, start, taper, tapsep, type, varest, 
                    weighted, copula,X,memdist,nosym,p_neighb = 1, thin_method = "bernoulli")
  {
    if (!isTRUE(memdist)) {
      warning(
        "'memdist = FALSE' is deprecated and ignored; WlsStart now always precomputes the selected pair structure for GeoFit.",
        call. = FALSE
      )
    }
    memdist <- TRUE

    fit_context <- getOption("GeoModels.fit_context_cache", NULL)
    if (is.environment(fit_context) &&
        exists("wls_initparam", envir = fit_context, inherits = FALSE)) {
      initparam <- get("wls_initparam", envir = fit_context, inherits = FALSE)

      ## Reproduce the data layout prepared by StartParam on the first call.
      ## Only values change between bootstrap replications; dimensions and
      ## ordering remain fixed.  This preserves the exact input seen by the
      ## likelihood while avoiding the repeated static setup.
      data_template <- initparam$data
      if (!is.list(data_template)) {
        data <- as.vector(data)
        if (!is.null(dim(data_template)) &&
            length(data) == length(data_template)) {
          dim(data) <- dim(data_template)
        }
      }
      initparam$data <- data
      return(initparam)
    }

    # Determines the range of the parameters for a given correlation
    SetRangeParam <- function(namesparam, numparam)
    {
        low <- 0#.Machine$double.eps
        big<- Inf
        lower <- NULL
        upper <- NULL
        # Check for the param set:
        for(i in 1:numparam){ 
            if (grepl("^mean[0-9]*$", namesparam[i]) ||
                grepl("^mean_[12][0-9]*$", namesparam[i])) {
                lower <- c(lower, -big)
                upper <- c(upper, big)
            }
            if(namesparam[i]=='skew'){
                lower <- c(lower, -big)
                upper <- c(upper, big)}
            if(namesparam[i]=='tail'){
                lower <- c(lower, low)
                upper <- c(upper, if (model %in% c('Tukeyh','tukeyh')) 0.5 - sqrt(.Machine$double.eps) else big)}
            if(namesparam[i] %in% c('tail1','tail2')){
                lower <- c(lower, low)
                upper <- c(upper, 0.5 - sqrt(.Machine$double.eps))}
           if(namesparam[i]=='df'){
                lower <- c(lower, 0)
                upper <- c(upper, big)}
            if(namesparam[i]=='nugget'){
                lower <- c(lower, low)
                upper <- c(upper, big)}
             if(namesparam[i]=='nugget1'){
                lower <- c(lower, low)
                upper <- c(upper, big)}
          if(namesparam[i]=='nugget2'){
                lower <- c(lower, low)
                upper <- c(upper, big)}
            if(namesparam[i]=='power'){
                lower <- c(lower, low)
                upper <- c(upper, 2)}
            if(namesparam[i]=='power_s'){
                lower <- c(lower, low)
                upper <- c(upper, 2)}
            if(namesparam[i]=='power_t'){
                lower <- c(lower, low)
                upper <- c(upper, 2)}
            if(namesparam[i]=='power1'){
                lower <- c(lower, low)
                upper <- c(upper, 2)}
            if(namesparam[i]=='scale'||namesparam[i]=='scale_11'||namesparam[i]=='scale_12'||namesparam[i]=='scale_22'){
                lower <- c(lower, low)
                upper <- c(upper, big)}
                  if(namesparam[i]=='smooth_11'||namesparam[i]=='smooth_12'||namesparam[i]=='smooth_22'){
                lower <- c(lower, low)
                upper <- c(upper, big)}  
            if(namesparam[i]=='power_2'||namesparam[i]=='power2_11'||namesparam[i]=='power2_12'||namesparam[i]=='power2_22'){
                lower <- c(lower, low)
                upper <- c(upper, big)}       
            if(namesparam[i]=='scale_s'){
                lower <- c(lower, low)
                upper <- c(upper, big)}
            if(namesparam[i]=='scale_t'){
                lower <- c(lower, low)
                upper <- c(upper, big)}
            if(namesparam[i]=='sep'){
                lower <- c(lower, if (!is.na(.GeoSTCode(corrmodel))) -1 else low)
                upper <- c(upper, 1)}
            if(namesparam[i]=='pcol'){
                lower <- c(lower, -1)
                upper <- c(upper,  1)}
            if(namesparam[i]=='sill'||namesparam[i]=='sill_1'||namesparam[i]=='sill_2'||namesparam[i]=='nugget_1'||namesparam[i]=='nugget2_2'){
                lower <- c(lower, low)
                upper <- c(upper, big)}
            if(namesparam[i] %in% c('smooth_s', 'smooth_t') &&
               !is.na(.GeoSTCode(corrmodel))) {
                lo <- if (isTRUE(CkCorrModel(corrmodel) == 85L))
                  -0.5 + sqrt(.Machine$double.eps) else .Machine$double.eps
                lower <- c(lower, lo)
                upper <- c(upper, big)
            }
            if(namesparam[i] %in% c('power2_s', 'power2_t') &&
               isTRUE(CkCorrModel(corrmodel) == 85L)) {
                lower <- c(lower, .Machine$double.eps)
                upper <- c(upper, big)
            }
            if(namesparam[i]=='smooth'){
                lower <- c(lower, low)
                upper <- c(upper, big)}
              }
              
        return(list(lower=lower, upper=upper))
    }
    
    ### Initialization parameters:
    
    initparam <- StartParam(coordx, coordy, coordz,coordt, coordx_dyn, corrmodel, data, distance, fcall, fixed,
                           grid, likelihood, maxdist,neighb, maxtime, model, n, 
                           param, parscale, paramrange, radius,  start, taper, tapsep,
                           "GeoWLS", type,  weighted,copula, X, memdist, nosym,
                           p_neighb, thin_method)
  

    
    if(!is.null(initparam$error))     stop(initparam$error)
    if(length(coordt)>0&&is.list(X)) X=X[[1]]

    bivariate<-CheckBiv(CkCorrModel(corrmodel))
    

    if(!bivariate) {if(is.null(X))  {X=1;num_betas=1} 
                    else num_betas=ncol(X)  }

    if( bivariate) {if(is.null(X))  {X=1;num_betas=c(1,1)} 
                    else 
                       { if(is.list(X)) num_betas=c(ncol(X[[1]]),ncol(X[[2]])) 
                         else num_betas=c(ncol(X),ncol(X))
                       } 
                   } 

    
    ### Set the initial type of likelihood objects:
    initparam$type <- CkType(type)
   # if(substr(model,1,6)=='Binary'||substr(model,1,8)=='Binomial'||substr(model,1,11)=='BinomialNeg'||substr(model,1,4)=='Geom'||substr(model,1,4)=='Pois') return(initparam)
    start_user <- if(is.null(start)) numeric(0) else unlist(start)
    fixed_user <- if(is.null(fixed)) numeric(0) else unlist(fixed)

    .reconcile_mean_parameters <- function(initparam, expected,
                                           start_values, fixed_values) {
      start_names <- names(start_values)
      fixed_names <- names(fixed_values)
      if(is.null(start_names)) start_names <- character(0)
      if(is.null(fixed_names)) fixed_names <- character(0)

      overlap <- intersect(expected, intersect(start_names, fixed_names))
      if(length(overlap)) {
        stop("Mean coefficient(s) cannot be supplied in both start and fixed: ",
             paste(overlap, collapse=", "), call.=FALSE)
      }

      ## StartParam supplies data-driven defaults for the mean coefficients.
      ## They are placeholders only: explicit fixed values have precedence,
      ## followed by explicit starting values.  Missing coefficients remain
      ## estimated and use the StartParam defaults as their initial values.
      defaults <- c(initparam$param, initparam$fixed)
      par_now <- initparam$param
      fix_now <- initparam$fixed
      if(is.null(par_now)) par_now <- numeric(0)
      if(is.null(fix_now)) fix_now <- numeric(0)

      for(nm in expected) {
        par_now <- par_now[names(par_now) != nm]
        fix_now <- fix_now[names(fix_now) != nm]

        if(nm %in% fixed_names) {
          value <- as.numeric(fixed_values[[nm]])
          if(length(value) != 1L || !is.finite(value))
            stop("Fixed mean coefficient '", nm, "' must be a finite scalar.",
                 call.=FALSE)
          fix_now[nm] <- value
          initparam$flagnuis[nm] <- 0
        } else {
          if(nm %in% start_names) {
            value <- as.numeric(start_values[[nm]])
          } else if(nm %in% names(defaults)) {
            value <- as.numeric(defaults[[nm]])
          } else {
            value <- 0
          }
          if(length(value) != 1L || !is.finite(value))
            stop("Starting value for mean coefficient '", nm,
                 "' must be a finite scalar.", call.=FALSE)
          par_now[nm] <- value
          initparam$flagnuis[nm] <- 1
        }
      }

      if(length(par_now)) par_now <- par_now[sort(names(par_now))]
      initparam$param <- par_now
      initparam$fixed <- if(length(fix_now)) fix_now else NULL
      initparam$namesparam <- names(initparam$param)
      initparam$numparam <- length(initparam$param)
      initparam$namesfixed <- names(initparam$fixed)
      initparam$numfixed <- length(initparam$fixed)
      initparam
    }

    ### Parameters omitted from start use the conservative defaults created by
    ### StartParam().  Repair defaults that would otherwise lie on an invalid
    ### or numerically dangerous boundary.  Explicit user starts are never
    ### changed here.
    expected_mean <- if(!bivariate) {
      .GeoMean_names(num_betas)
    } else {
      c(.GeoMean_names(num_betas[1], "mean_1"),
        .GeoMean_names(num_betas[2], "mean_2"))
    }
    start_names_user <- names(start_user)
    if(is.null(start_names_user)) start_names_user <- character(0)

    .repair_automatic_starts <- function(initparam, start_names, mean_names,
                                         model_name, data_values) {
      par_now <- initparam$param
      if(is.null(par_now) || !length(par_now)) return(initparam)

      auto_names <- setdiff(names(par_now), c(start_names, mean_names))
      zz <- suppressWarnings(as.numeric(unlist(data_values, use.names = FALSE)))
      zz <- zz[is.finite(zz)]
      vv <- if(length(zz) > 1L) stats::var(zz) else NA_real_
      if(!is.finite(vv) || vv <= 0) vv <- 1
      rr <- if(length(zz)) range(zz) else c(0, 1)
      span <- diff(rr)
      if(!is.finite(span) || span <= 0) span <- max(1, abs(rr[1L]))
      pad <- max(sqrt(.Machine$double.eps), 0.01 * span)

      for(nm in auto_names) {
        cur <- suppressWarnings(as.numeric(par_now[[nm]]))
        if(length(cur) != 1L || !is.finite(cur)) cur <- NA_real_

        if(grepl("^nugget($|_)", nm)) cur <- 0.05
        else if(grepl("^sill($|_)", nm)) cur <- vv
        else if(grepl("^scale($|_|[st]$)", nm) || grepl("^scale_[0-9]", nm)) {
          if(!is.finite(cur) || cur <= 0) cur <- 1
        }
        else if(grepl("^smooth($|_)", nm)) cur <- 0.5
        else if(grepl("^power($|_)", nm) || grepl("^power[12]($|_)", nm)) cur <- 1
        else if(nm == "sep") cur <- 0.5
        else if(nm == "pcol") cur <- 0
        else if(grepl("^shape($|[12]$|_)", nm)) cur <- 1
        else if(grepl("^df($|_)", nm)) cur <- 0.2
        else if(grepl("^tail($|[12]$|_)", nm)) cur <- 0.05
        else if(grepl("^skew($|_)", nm)) {
          cur <- if(model_name == "SkewLaplace") 0.5 else 0
        }
        else if(nm == "pmu") cur <- 0
        else if(nm == "min") cur <- rr[1L] - pad
        else if(nm == "max") cur <- rr[2L] + pad
        else if(!is.finite(cur)) cur <- 1

        par_now[[nm]] <- cur
      }

      ## Ensure an automatically estimated bounded-support interval has
      ## strictly positive width even for constant observations.
      if(all(c("min", "max") %in% names(par_now)) &&
         !any(c("min", "max") %in% start_names)) {
        if(!is.finite(par_now[["min"]]) || !is.finite(par_now[["max"]]) ||
           par_now[["min"]] >= par_now[["max"]]) {
          par_now[["min"]] <- rr[1L] - pad
          par_now[["max"]] <- rr[2L] + pad
        }
      }

      initparam$param <- par_now
      initparam$namesparam <- names(par_now)
      initparam$numparam <- length(par_now)
      initparam
    }

    initparam <- .repair_automatic_starts(
      initparam, start_names_user, expected_mean, model, data
    )
    initparam <- .GeoSTRepairStarts(initparam, corrmodel, start_names_user)


    ## Wrapped Gaussian: use circular rather than Euclidean automatic starts.
    ## If the user did not provide them explicitly, initialize the intercept
    ## from the sample mean direction and the latent variance from the
    ## resultant length R, since for a wrapped normal R = exp(-sill/2).
    if(identical(model, "Wrapped") && !bivariate) {
      theta <- suppressWarnings(as.numeric(unlist(data, use.names = FALSE)))
      theta <- theta[is.finite(theta)]
      if(length(theta)) {
        Cbar <- mean(cos(theta))
        Sbar <- mean(sin(theta))
        mean_dir <- atan2(Sbar, Cbar) %% (2*pi)
        Rbar <- sqrt(Cbar^2 + Sbar^2)
        Rbar <- min(1 - sqrt(.Machine$double.eps),
                    max(sqrt(.Machine$double.eps), Rbar))
        eta0 <- tan((mean_dir - pi)/2)
        if(!is.finite(eta0)) eta0 <- if(mean_dir < pi) -100 else 100
        eta0 <- max(-100, min(100, eta0))
        sill0 <- -2 * log(Rbar)

        ## StartParam can temporarily store automatically generated mean
        ## coefficients in either `param` or `fixed`; WlsStart later reconciles
        ## those placeholders.  Replace the placeholder wherever it lives, but
        ## never overwrite a value explicitly supplied by the user.
        if(!("mean" %in% start_names_user) &&
           !("mean" %in% names(fixed_user))) {
          if("mean" %in% names(initparam$param)) initparam$param[["mean"]] <- eta0
          if("mean" %in% names(initparam$fixed)) initparam$fixed[["mean"]] <- eta0
        }

        slope_names <- setdiff(expected_mean, "mean")
        for(nm in slope_names) {
          if(!(nm %in% start_names_user) && !(nm %in% names(fixed_user))) {
            if(nm %in% names(initparam$param)) initparam$param[[nm]] <- 0
            if(nm %in% names(initparam$fixed)) initparam$fixed[[nm]] <- 0
          }
        }

        if("sill" %in% names(initparam$param) &&
           !("sill" %in% start_names_user) &&
           !("sill" %in% names(fixed_user)))
          initparam$param[["sill"]] <- sill0
      }
    }

    ###### continuous and discrete models
    if((model %in% c('Gaussian','Gauss','Chisq','LogLogistic','Logistic','SkewLaplace','Gamma','Gamma2','Beta','Beta2','LogGaussian','LogGauss','Binomial_TwoPieceGaussian','Binomial_TwoPieceGauss',
          'Tukeygh','Tukeyh','Tukeyh2','Kumaraswamy','Kumaraswamy2','Weibull','SkewGaussian','SkewGauss','SinhAsinh','StudentT','SkewStudentT',
          "Gaussian_misp_StudentT","Gaussian_misp_Poisson","Gaussian_misp_Tukeygh",
          "Binomial","Geometric","PoisBin","BinomialNeg","Binary_misp_BinomialNeg","Binomial2","PoisBinNeg","Poisson","BinomialNegZINB","BinomialLogistic","Gaussian_misp_Binomial",
          "Gaussian_misp_BinomialNeg","PoissonZIP1","BinomialNegZINB1",
          "Gaussian_misp_SkewStudentT","PoissonGamma","PoissonWeibull","Gaussian_misp_PoissonGamma","PoissonGammaZIP1",
          "TwoPieceStudentT",'Wrapped',"TwoPieceGaussian","TwoPieceGauss","TwoPieceTukeyh","TwoPieceBimodal")) &
          (type %in% c('Standard','Pairwise','Tapering','Tapering1','Independence'))) {

      if(!initparam$bivariate) {
        mean_names <- .GeoMean_names(num_betas)
        initparam <- .reconcile_mean_parameters(
          initparam, mean_names, start_user, fixed_user
        )
      } else {
        mean_names <- c(.GeoMean_names(num_betas[1], "mean_1"),
                        .GeoMean_names(num_betas[2], "mean_2"))
        initparam <- .reconcile_mean_parameters(
          initparam, mean_names, start_user, fixed_user
        )
      }

    }

    paramrange <- SetRangeParam(names(initparam$param), length(initparam$param))
    initparam$lower <- paramrange$lower
    initparam$upper <- paramrange$upper

    if (is.environment(fit_context)) {
      assign("wls_initparam", initparam, envir = fit_context)
    }
    return(initparam)
  }

####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
GeoWLS <- function(data, coordx, coordy=NULL,coordz=NULL, coordt=NULL,  coordx_dyn=NULL, corrmodel, distance="Eucl",
                         fixed=NULL,grid=FALSE, maxdist=NULL, neighb=NULL,maxtime=NULL,
                         optimizer='Nelder-Mead', numbins=NULL, radius=1,  start=NULL,
                         weighted=FALSE,optimization=TRUE,numbins_t=NULL)
  {
    call <- match.call()
    grid_requested <- isTRUE(grid)
  .st_old <- .GeoSTEnter(corrmodel,
    coordx = coordx, coordy = coordy,
    coordz = coordz, coordx_dyn = coordx_dyn,
    distance = distance, context = "GeoWLS")
  on.exit(.GeoSTRestore(.st_old), add = TRUE)

    ## GeoWLS is a Gaussian-only variogram least-squares procedure.  Keep the
    ## model label internal because shared setup code and GeoCovariogram use it,
    ## but do not expose a public `model` argument.
    model <- 'Gaussian'

    if (grid_requested) {
      if (!is.null(coordx_dyn)) stop("grid=TRUE cannot be combined with coordx_dyn")
      is_biv_wls <- isTRUE(CheckBiv(CkCorrModel(corrmodel)))
      gg <- .GeoVariogram_normalize_grid(
        data=data, coordx=coordx, coordy=coordy, coordz=coordz,
        coordt=coordt, bivariate=is_biv_wls
      )
      data <- gg$data
      coordx <- gg$coordx
      coordy <- gg$coordy
      coordz <- gg$coordz
      grid <- FALSE
    }

    if (!is.null(maxdist)) {
      if (!is.numeric(maxdist) || any(!is.finite(maxdist)) || any(maxdist <= 0))
        stop("maxdist must contain positive finite numeric values")
    }
    .GeoVariogram_validate_data(data, context="GeoWLS")

    ### Check the parameters given in input:
    checkinput <- CkInput(coordx, coordy, coordz,coordt, coordx_dyn, corrmodel, data, distance,"Fitting", fixed, grid, 'None',
                              maxdist, maxtime,  model,NULL,  optimizer, NULL, radius, start, NULL,
                             NULL, 'GeoWLS', FALSE, weighted, NULL,NULL)
    

    if(!is.null(checkinput$error))
      stop(checkinput$error)
    # Check the number of spatial breakpoints (legacy GeoModels meaning:
    # numbins breakpoints define numbins - 1 empirical spatial bins).
    if(is.null(numbins)) numbins <- 13L
    if(!is.numeric(numbins) || length(numbins)!=1L || !is.finite(numbins) ||
       numbins < 2 || numbins %% 1 != 0)
      stop('numbins must be an integer >= 2\n')
    numbins <- as.integer(numbins)
    if(!is.null(numbins_t)) {
      if(!is.numeric(numbins_t) || length(numbins_t)!=1L || !is.finite(numbins_t) ||
         numbins_t < 1 || numbins_t %% 1 != 0)
        stop('numbins_t must be a positive integer\n')
      numbins_t <- as.integer(numbins_t)
    }
    ### Define the object function for the weighted least squares method:
    WLsquare <- function(param,bins, bint, corrmodel, fixed, fun, lenbins, moments,
                         namescorr, namesnuis, numbins, numbint)
      {
        param <- c(param, fixed)#set the parameters set:
        paramcorr <- param[namescorr]#set the correlation parameters:
        nuisance <- param[namesnuis]#set the nuisance parameters:

     if(fun=='WLeastSquare_G')
        result <- .C('WLeastSquare_G', as.double(bins), as.double(bint), as.integer(corrmodel),
                     as.double(lenbins), as.double(moments), as.integer(numbins),
                     as.integer(numbint), as.double(nuisance), as.double(paramcorr),
                     res=double(1), PACKAGE='GeoModels', NAOK=TRUE)$res
     if(fun=='LeastSquare_G')
        result <- .C('LeastSquare_G', as.double(bins), as.double(bint), as.integer(corrmodel),
                     as.double(lenbins), as.double(moments), as.integer(numbins),
                     as.integer(numbint), as.double(nuisance), as.double(paramcorr),
                     res=double(1), PACKAGE='GeoModels', NAOK=TRUE)$res
       return(result)
      }
    ### Initializes global variables:
    GeoWLS <- NULL
    fname <- NULL
    variogramt <- NULL
    variogramst <- NULL
    ### Initializes the parameter values:
    parscale <- NULL
    initparam <- StartParam(coordx, coordy, coordz,coordt, coordx_dyn,corrmodel, data, distance, "Fitting", fixed, grid,
    'None', maxdist, neighb,maxtime,  model, NULL,  NULL,parscale, TRUE, radius, start,NULL,  NULL,
                           'GeoWLS', 'GeoWLS', FALSE, NULL, NULL, FALSE, FALSE,
                           1, "bernoulli")


    if(!is.null(initparam$error))
      stop(initparam$error)
     coordx=initparam$coordx
     coordy=initparam$coordy
     coordz=initparam$coordz
     coordt=initparam$coordt
     on.exit({
       try(.C('DeleteGlobalVar', PACKAGE='GeoModels', NAOK=TRUE), silent=TRUE)
     }, add=TRUE)
     ## Native variogram routines always index the third coordinate.  Use an
     ## explicit zero vector in 2D instead of passing NULL or a scalar zero.
     if(is.null(coordz)) coordz <- rep(0, length(coordx))
    ###### ----------- START Estimation of the empirical variogram ---------- #####
    numvario <- numbins-1
    bins <- double(numbins) # vector of spatial bins
    moments <- double(numvario) # vector of spatial moments
    lenbins <- integer(numvario) # vector of spatial bin sizes
    ### Checks the type of variogram:
    fname <- 'Binned_Variogram'


    if(initparam$spacetime){### Computes the spatial-temporal variogram:
      spacetime_dyn <- !is.null(coordx_dyn)
      ns <- as.integer(initparam$ns)
      NS <- as.integer(c(0L, cumsum(ns))[-(length(ns) + 1L)])

      ti <- .GeoVariogram_time_bins(
        coordt, maxtime,
        numbins_t = if (is.null(numbins_t)) numbins - 1L else numbins_t
      )
      time_breaks <- ti$breaks
      temporal_centers <- ti$centers
      ntbin <- length(time_breaks) - 1L

      bint_native <- as.double(time_breaks)
      momentt <- double(ntbin)
      lenbint_native <- integer(ntbin)
      momentst <- double(numvario * ntbin)
      lenbinst_native <- integer(numvario * ntbin)

      ## The native routines expect observations in time-block order,
      ## consistently with the coordinates expanded by StartParam().
      if(!spacetime_dyn) data_native <- c(t(data))
      else data_native <- unlist(data, use.names = FALSE)

      if(spacetime_dyn) {
        EV=.C("Binned_Variogram_st2_dyn", bins=bins, bint=bint_native,
             as.double(coordx),as.double(coordy),as.double(coordz),
             as.double(coordt),as.double(data_native),
             lenbins=lenbins,lenbinst=lenbinst_native,lenbint=lenbint_native,
             moments=moments, momentst=momentst, momentt=momentt,
             as.integer(numbins), as.integer(length(time_breaks)),
             as.integer(ns),as.integer(NS), PACKAGE='GeoModels', NAOK=TRUE)
      } else {
        EV=.C("Binned_Variogram_st2", bins=bins, bint=bint_native,
             as.double(coordx),as.double(coordy),as.double(coordz),
             as.double(coordt),as.double(data_native),
             lenbins=lenbins,lenbinst=lenbinst_native,lenbint=lenbint_native,
             moments=moments, momentst=momentst, momentt=momentt,
             as.integer(numbins), as.integer(length(time_breaks)),
             as.integer(ns),as.integer(NS), PACKAGE='GeoModels', NAOK=TRUE)
      }

      bins_native <- EV$bins
      spatial_count <- EV$lenbins
      spatial_sum <- EV$moments
      temporal_count <- EV$lenbint
      temporal_sum <- EV$momentt
      st_count <- matrix(EV$lenbinst, nrow=numvario, ncol=ntbin, byrow=TRUE)
      st_sum <- matrix(EV$momentst, nrow=numvario, ncol=ntbin, byrow=TRUE)

      centers <- bins_native[-length(bins_native)] + diff(bins_native)/2
      variograms <- rep(NA_real_, numvario)
      ii <- spatial_count > 0L
      variograms[ii] <- spatial_sum[ii] / spatial_count[ii]

      variogramt <- rep(NA_real_, ntbin)
      ii <- temporal_count > 0L
      variogramt[ii] <- temporal_sum[ii] / temporal_count[ii]

      st_mean <- matrix(NA_real_, nrow=numvario, ncol=ntbin)
      ii <- st_count > 0L
      st_mean[ii] <- st_sum[ii] / st_count[ii]
      variogramst <- as.numeric(as.vector(t(st_mean)))

      ## Build the rectangular WLS system without dropping empty cells.  The
      ## C objective now skips zero-count cells, preserving the lag geometry.
      moment_full <- rbind(temporal_sum, st_sum)
      moment_full <- cbind(c(0, spatial_sum), moment_full)
      count_full <- rbind(temporal_count, st_count)
      count_full <- cbind(c(1L, spatial_count), count_full)

      bins <- c(-bins_native[1L], bins_native)
      bint <- c(0, temporal_centers)
      moments <- moment_full
      lenbins <- count_full
      numbins <- length(bins)
      numbint <- length(bint)
    }
    else{### Computes the spatial variogram:
      numbint <- 1L
      bint <- 0
      EV=.C("Binned_Variogram2", bins=bins, as.double(coordx),as.double(coordy),
         as.double(coordz),as.double(coordt),as.double(initparam$data), lenbins=lenbins,
         moments=moments, as.integer(numbins),PACKAGE='GeoModels', NAOK=TRUE)
      bins <- EV$bins
      lenbins <- EV$lenbins
      moments <- EV$moments

      centers <- bins[-length(bins)] + diff(bins)/2
      variograms <- rep(NA_real_, length(lenbins))
      ii <- lenbins > 0L
      variograms[ii] <- moments[ii] / lenbins[ii]

      if(!is.null(initparam$param['scale']) && any(is.finite(variograms))) {
        initparam$param['scale'] <- centers[which.max(replace(variograms, !is.finite(variograms), -Inf))]
      }
    }
    ###### ----------- END Estimation of the empirical variogram ---------- #####



    ###### ---------- START model fitting by weighted least squares method ----------######
    #if(model=='Gaussian'||model=='Gauss') # Gaussian random field:
      if(weighted) fname <- 'WLeastSquare_G'
      else fname <- 'LeastSquare_G'


if(optimization){
    ### Computes estimates by the weighted least squares method:
    if(optimizer=='L-BFGS-B')
      fitted <- optim(initparam$param, WLsquare, bins=bins, bint=bint, corrmodel=initparam$corrmodel,
                      fixed=initparam$fixed, fun=fname, lenbins=lenbins, method=optimizer,
                      moments=moments, namescorr=initparam$namescorr, namesnuis=initparam$namesnuis,
                      numbins=numbins, numbint=numbint, control=list(fnscale=-1, factr=1, pgtol=1e-14,
                      maxit = 1e8), lower=initparam$lower, upper=initparam$upper, hessian=FALSE)
    else
      fitted <- optim(initparam$param, WLsquare, bins=bins, bint=bint, corrmodel=initparam$corrmodel,
                      fixed=initparam$fixed, fun=fname, lenbins=lenbins, method=optimizer,
                      moments=moments, namescorr=initparam$namescorr, namesnuis=initparam$namesnuis,
                      numbins=numbins, numbint=numbint, control=list(fnscale=-1, reltol=1e-14, maxit=1e8),
                      hessian=FALSE)
     }  
     else{ res=WLsquare(unlist(start),bins=bins, bint=bint, corrmodel=initparam$corrmodel,
                      fixed=initparam$fixed, fun=fname, lenbins=lenbins, 
                      moments=moments, namescorr=initparam$namescorr, namesnuis=initparam$namesnuis,
                      numbins=numbins, numbint=numbint)
              fitted=list(convergence=NULL,counts=NULL,message=NULL,par=unlist(start),value=res)

            }


    ###### ---------- END model fitting by weighted least squares method ----------######
    ### Set the output:
    ## Keep the structural metadata required by downstream methods such as
    ## GeoCovariogram().  These fields describe the fitted Gaussian field;
    ## they do not re-introduce a user-selectable model argument.
    ##
    ## StartParam() represents the default intercept-only design by the scalar
    ## value 1.  Preserve the same meaning but store an explicit n x 1 matrix
    ## in GeoWLS objects.  Mean utilities intentionally reject ambiguous
    ## vectors, so this canonical representation is important for downstream
    ## methods (notably GeoCovariogram).
    X_out <- initparam$X
    if(!initparam$bivariate) {
      nobs_out <- length(unlist(initparam$data, use.names=FALSE))
      if(is.null(X_out) ||
         (is.atomic(X_out) && is.null(dim(X_out)) && length(X_out)==1L &&
          length(initparam$num_betas)==1L && initparam$num_betas==1L)) {
        X_out <- matrix(1, nrow=nobs_out, ncol=1L)
      } else {
        X_out <- .GeoMean_design(X_out, n=nobs_out, p=initparam$num_betas,
                                 name="GeoWLS$X")
      }
    }

    GeoWLS <- list(      bins=bins,
                         bint=bint,
                         centers=centers,
                         bivariate=initparam$bivariate,
                         spacetime=initparam$spacetime,
                         coordx = initparam$coordx,
                         coordy = initparam$coordy,
                         coordz = initparam$coordz,
                         coordt = initparam$coordt,
                         coordx_dyn=coordx_dyn,
                         convergence = fitted$convergence,
                         corrmodel = corrmodel,
                         copula=NULL,
                         data = initparam$data,
                         distance=initparam$distance,
                         fixed = initparam$fixed,
                         grid = grid_requested,
                         iterations = fitted$counts,
                         maxdist =maxdist,
                         maxtime = maxtime,
                         numbins_t = numbins_t,
                         message = fitted$message,
                         model=model,
                         ## `n` is irrelevant for continuous Gaussian fields,
                         ## but shared covariance utilities pass it to the native
                         ## correlation engine.  Use the same harmless scalar
                         ## default as GeoCorrFct() instead of returning NULL.
                         n=1L,
                         numbetas=initparam$num_betas,
                         numcoord=initparam$numcoord,
                         numtime=initparam$numtime,
                         param = fitted$par,
                         radius=initparam$radius,
                         X=X_out,
                         variograms = variograms,
                         variogramt = variogramt,
                         variogramst = variogramst,
                         weighted = weighted,
                         wls = fitted$value)
    structure(c(GeoWLS, call = call), class = c("GeoWLS"))
  }
