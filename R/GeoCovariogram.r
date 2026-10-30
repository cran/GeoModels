####################################################
### File name: GeoCovVariogram.r
####################################################

### Compute and plot the (estimated) covariance function and the variogram
### from  GeoFit object

## Select the empirical temporal profile used by GeoCovariogram for a
## dynamic space-time variogram.  The true temporal margin is gamma(0,u).
## When it is unavailable and fix.lags requests h=0, use the first retained
## positive-distance ST bin with finite values.  This helper changes only the
## diagnostic profile; it never alters the GeoVariogram object itself.
.GeoCovariogramDynamicTemporalProfile <- function(vario, evario, lags,
                                                   fix.lags, dynamic=FALSE) {
  idx <- fix.lags
  value <- if(!is.null(idx) && length(lags) >= idx) lags[idx] else NA_real_
  out <- list(index=idx, value=value, fallback=FALSE,
              lower=NA_real_, upper=NA_real_, first.bin=FALSE)

  if(!isTRUE(dynamic) || is.null(idx) || idx != 1L ||
     !is.matrix(evario) || nrow(evario) <= 1L || ncol(evario) <= 1L)
    return(out)

  ## If collocated dynamic pairs exist, retain the exact h=0 profile.
  if(any(is.finite(evario[1L,-1L]))) return(out)

  positive_surface <- evario[-1L,-1L,drop=FALSE]
  candidates <- which(rowSums(is.finite(positive_surface)) > 0L)
  if(!length(candidates)) return(out)

  out$index <- candidates[1L] + 1L
  if(length(lags) >= out$index) out$value <- lags[out$index]
  out$fallback <- TRUE

  ## GeoVariogram may omit empty retained centers.  Recover the original
  ## bin interval by matching the retained center to the full midpoint grid.
  br <- as.numeric(vario$bins)
  if(length(br) >= 2L && all(is.finite(br)) && is.finite(out$value)) {
    mids <- br[-length(br)] + diff(br)/2
    if(length(mids)) {
      jj <- which.min(abs(mids - out$value))
      out$lower <- br[jj]
      out$upper <- br[jj + 1L]
      out$first.bin <- identical(as.integer(jj),1L)
    }
  }
  out
}

GeoCovariogram <- function(fitted, distance="Eucl", answer.cov=FALSE, answer.vario=FALSE,
                         answer.range=FALSE, fix.lags=NULL, fix.lagt=NULL,
                        show.cov=FALSE, show.vario=TRUE,  show.range=FALSE,
                        add.cov=FALSE, add.vario=FALSE, pract.range=95,
                        vario=NULL,invisible=FALSE,...)
  {
    result <- NULL

    CheckDistance<- function(distance)
    {
        CheckDistance <- NULL
        CheckDistance <- switch(distance,
                                eucl=0,
                                Eucl=0,
                                chor=1,  
                                Chor=1,
                                geod=2,
                                Geod=2,
                                proj=3,
                                Proj=3)
        return(CheckDistance)
    }
CorrelationFct <- function(bivariate, corr_model, lags, lag_t, num_lags, num_lag_t, mu, model, nuisance, param, N) {
  if (!bivariate) {
    result <- dotCall64::.C64(
      "VectCorrelation",
      SIGNATURE = c("double", "integer", "double", "integer", "integer", "double", "integer", "double", "double", "double", "integer"),
      INTENT = c("w", rep("r", 10)),
      corr = dotCall64::vector_dc("double", num_lags * num_lag_t),
      corr_model = corr_model,
      lags = lags,num_lags = num_lags,num_lag_t = num_lag_t,mu = mu,
      model = model,nuisance = nuisance,param = param,lag_t = lag_t,N = N,PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
    
    cc <- result$corr
  } else {
    result <- dotCall64::.C64(
      "VectCorrelation_biv",
      SIGNATURE = c("double", "double", "integer", "double", "integer", "integer", "double", "integer", "double", "double", "double", "integer"),
      INTENT = c("w", "w", rep("r", 10)),
      corr = dotCall64::vector_dc("double", num_lags * 4),vario = dotCall64::vector_dc("double", num_lags * 4),corr_model = corr_model,
      lags = lags,num_lags = num_lags,num_lag_t = num_lag_t,mu = mu,
      model = model,nuisance = nuisance,param = param,lag_t = lag_t,N = N,PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
    
    cc <- c(result$corr, result$vario)
  }
  cc
}


.GeoCovariogram_sqerr <- function(empirical, theoretical) {
  empirical <- as.numeric(empirical)
  theoretical <- as.numeric(theoretical)
  if(length(empirical) != length(theoretical))
    stop("Empirical and theoretical semivariograms have incompatible dimensions.",
         call. = FALSE)
  ok <- is.finite(empirical) & is.finite(theoretical)
  if(!any(ok)) return(NA_real_)
  sum((empirical[ok] - theoretical[ok])^2)
}

.GeoCovariogram_nearest_index <- function(grid, value) {
  grid <- as.numeric(grid)
  if(!length(grid) || !is.finite(value)) return(NA_integer_)
  which.min(abs(grid - value))
}

.GeoCovariogram_finite_max <- function(x, default = 0) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if(!length(x)) return(as.numeric(default)[1L])
  max(x)
}

.GeoCovariogram_practical_range <- function(lag, variogram, covariance, fraction) {
  lag <- as.numeric(lag)
  variogram <- as.numeric(variogram)
  covariance <- as.numeric(covariance)
  ok <- is.finite(lag) & is.finite(variogram) & is.finite(covariance)
  if(sum(ok) < 2L) return(NA_real_)
  lag <- lag[ok]; variogram <- variogram[ok]; covariance <- covariance[ok]
  oo <- order(lag)
  lag <- lag[oo]; variogram <- variogram[oo]; covariance <- covariance[oo]
  variance <- stats::median(variogram + covariance, na.rm = TRUE)
  if(!is.finite(variance) || variance <= 0) return(NA_real_)
  target <- fraction * variance
  jj <- which(variogram >= target)[1L]
  if(is.na(jj)) return(NA_real_)
  if(jj <= 1L) return(lag[1L])
  y <- variogram[c(jj - 1L, jj)]
  x <- lag[c(jj - 1L, jj)]
  if(!all(is.finite(y)) || diff(y) == 0) return(x[2L])
  stats::approx(x = y, y = x, xout = target, ties = "ordered")$y
}


if(show.range || answer.range)  {
    # Practical range in the Gaussian case:
    PracRangeNorm <- function(corrmodel, lags, lagt, nuisance, numlags, numlagt, param, pract.range)
    { 
        return(nuisance["sill"]*(1-nuisance["nugget"])*CorrelationFct(bivariate,corrmodel, lags, lagt, numlags, numlagt, mu,CkModel(fitted$model),nuisance, param,fitted_n)
            -(nuisance["sill"]*(1-pract.range)))
    
    }
}
###########################################################     
#############  starting ###################################
########################################################### 

    if(!inherits(fitted,"GeoFit") && !inherits(fitted,"GeoWLS"))
        stop("Enter an object obtained from GeoFit or WLeastSquare")
    if(inherits(fitted, "GeoFit") && identical(fitted$residual_type, "Pearson")) {
        stop(
            "GeoCovariogram() cannot compare a Pearson-residual variogram with the covariance model of the original discrete response. Use GeoCovariogram() on the original fit, or inspect the residual variogram without a fitted response-scale curve.",
            call. = FALSE
        )
    }

    dyn=FALSE
    OLS=NULL
    isvario <- !is.null(vario) # is empirical variogram is passed?
    bivariate <- fitted$bivariate
    if(bivariate) fitted$numtime=1
    ispatim <- fitted$spacetime
    dyn<- is.list(fitted$coordx_dyn)

    ## Continuous fitted objects do not need a trial/count parameter, but the
    ## shared native correlation interface requires one integer scalar.  Older
    ## GeoWLS objects may therefore have `n = NULL`; use the same default as
    ## GeoCorrFct() so GeoCovariogram remains backward compatible.
    fitted_n <- fitted$n
    if(is.null(fitted_n) || !length(fitted_n)) fitted_n <- 1L

    
    if(!isvario) stop("an object vario is needed\n")


if (fitted$model %in% c("Weibull", "Poisson","PoissonGamma", "Binomial", "Gamma", 
        "LogLogistic", "BinomialNeg", "Bernoulli", "Geometric", 
        "Gaussian_misp_Poisson", "PoissonZIP", "Gaussian_misp_PoissonZIP", 
        "Gaussian_misp_PoissonGamma","PoissonGammaZIP","PoissonGammaZIP1",
        "BinomialNegZINB", "PoissonZIP1", "Gaussian_misp_PoissonZIP1", 
        "BinomialNegZINB1", "Beta2", "Kumaraswamy2", "Beta", 
        "Kumaraswamy")) { fitted$param$sill=1}


do_plot <- isTRUE(show.cov) || isTRUE(show.vario)
if(do_plot) {
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  if(!ispatim && !bivariate && show.cov && show.vario) par(mfrow=c(1,2))
  if(show.vario && ispatim) par(mfrow=c(1,2))
  if(show.vario && ispatim && !is.null(fix.lags) && !is.null(fix.lagt)) par(mfrow=c(2,2))
  if(show.cov && ispatim && !show.vario) {
    nprof <- as.integer(!is.null(fix.lags)) + as.integer(!is.null(fix.lagt))
    if(nprof == 0L) par(mfrow=c(1,1))
    if(nprof == 1L) par(mfrow=c(1,2))
    if(nprof == 2L) par(mfrow=c(2,2))
  }
  if(bivariate && !dyn && (show.vario || show.cov)) par(mfrow=c(2,2))
  if(bivariate&&dyn) par(mfrow=c(1,2))
}

    ## Preserve a site-specific fixed mean before flattening the fitted
    ## parameter containers.  Otherwise unlist() converts fixed$mean into
    ## separate elements named mean1, mean2, ..., and the external mean can
    ## no longer be identified by .GeoMean_external().
    external_mean_fit <- NULL
    if(!bivariate) {
      nobs_mean_fit <- length(unlist(fitted$data, use.names = FALSE))
      external_mean_fit <- .GeoMean_external(
        fitted$fixed, n = nobs_mean_fit
      )
      if(!is.null(external_mean_fit)) {
        if(!is.list(fitted$fixed)) fitted$fixed <- as.list(fitted$fixed)
        ## Keep a scalar placeholder for the native/nuisance parameter paths;
        ## the original location-specific vector is retained separately above.
        fitted$fixed[["mean"]] <- 0
      }
    }

    fitted$param <- unlist(fitted$param)
    fitted$fixed <- unlist(fitted$fixed)

    ## Normalize the optional copula field.  Older GeoFit objects may encode
    ## "no copula" as NULL, NA, or an empty string.
    fit_copula <- fitted$copula
    if(!is.null(fit_copula)) {
      fit_copula <- as.character(fit_copula)[1L]
      if(is.na(fit_copula) || !nzchar(trimws(fit_copula))) fit_copula <- NULL
    }
    fitted$copula <- fit_copula
    
    # START ---- check input --- #

    if( !inherits(fitted,"GeoFit")  & !inherits(fitted,"GeoWLS") )
        stop("Enter an object obtained GeoFit or WLeastSquare")
  
    if(isvario & !inherits(vario,"GeoVariogram") )
        stop("Enter an object obtained from the function GeoVariogram")

    vario_bivariate <- isTRUE(vario$bivariate)
    vario_spacetime <- !is.null(vario$bint) || !is.null(vario$centert)
    if(vario_bivariate != isTRUE(bivariate))
        stop("The GeoVariogram and fitted objects disagree on bivariate structure.", call. = FALSE)
    if(vario_spacetime != isTRUE(ispatim))
        stop("The GeoVariogram and fitted objects disagree on space-time structure.", call. = FALSE)

    ## When omitted, inherit the metric stored by the fit.  The empirical lag
    ## values themselves come from `vario`, so explicit disagreement is a
    ## useful diagnostic even though the correlation function only sees lags.
    if(missing(distance) && !is.null(fitted$distance)) distance <- fitted$distance
    if(!is.null(fitted$distance) && length(fitted$distance) &&
       !identical(tolower(as.character(distance)[1L]),
                  tolower(as.character(fitted$distance)[1L])))
        warning("`distance` differs from the metric stored in the fitted object; ",
                "the GeoVariogram and fit should use the same distance metric.",
                call. = FALSE)

    if(!is.numeric(pract.range) || length(pract.range) != 1L ||
       is.na(pract.range) || !is.finite(pract.range))
        stop("Enter a finite numeric scalar for pract.range.", call. = FALSE)
    if(pract.range < 0 || pract.range > 100)
        stop("pract.range must lie between 0 and 100.", call. = FALSE)
    pract.range <- pract.range / 100
    # Set the type of process:
    model <- CkModel(fitted$model)
    gaussian <- model==1
    skewgausssian<- model==10
    skewlaplace<-model==59
    gamma<- model==21
    skewstudentT<- model==18||model==37
    studentT<- model==12||model==35
    weibull<- model==26
    twopieceT<- model==27
    twopieceTukeyh<- model==38 
    twopieceGauss<- model==29
    twopiecebimodal<- model==39
    loggauss<- model==22
    binary <- model==2
    binomial <- model==11
    binomial2 <- model==19
    geom <- model==14
    binomialneg<-model==16
    binomialnegZINB<-model==45
    tukeyh<- model ==34
    tukeyh2<- model ==40
    sas<- model ==20
    poisson<- model==30
    poissongamma<- model==46
    ## Zero-inflated count models.  Codes 53 and 58 are the one-nugget
    ## versions and must follow the same observed-scale covariance path.
    poissonZIP<- model %in% c(43,44,53)
    poissongammaZIP<- model %in% c(57,58)
    loglogistic <- model==24
    tukeygh<- model==9||model==41
    Gaussian_misp_Binomial<-model==51
    Gaussian_misp_Poisson<-model==36
    Gaussian_misp_PoissonGamma<- model==47
    Gaussian_misp_BinomialNeg<-model==52

    copula_cont <- !bivariate && !is.null(fit_copula) &&
      fit_copula %in% c("Gaussian", "Clayton", "SkewGaussian") &&
      as.integer(model) %in% .GeoGaussianCopulaContinuousModels

    ## Clayton-like and skew-Gaussian covariance is currently defined here
    ## for the audited continuous copula margins only.  Gaussian-copula
    ## discrete margins retain the legacy GeoCovariogram path below.
    if(!is.null(fit_copula) && fit_copula %in% c("Clayton", "SkewGaussian") &&
       !copula_cont)
      stop("GeoCovariogram supports copula='", fit_copula,
           "' for the continuous copula margins only.", call. = FALSE)

    zero <- 0;slow=1e-3;
    if(copula_cont || gaussian||skewgausssian||skewlaplace||gamma||loggauss||binomial||Gaussian_misp_Binomial||Gaussian_misp_BinomialNeg||Gaussian_misp_Poisson||binomialneg||binomialnegZINB||geom||tukeyh||tukeyh2||sas||twopiecebimodal||skewstudentT
            ||twopieceGauss||twopieceTukeyh||twopieceT) slow=1e-6
    else slow=1e-5 
    # lags associated to empirical variogram estimation
    if(isvario){
    lags <- c(0,vario$centers);numlags <- length(lags)
    if(ispatim) {
     if(is.null(vario$centert)) lagt <-c(0,vario$bint) 
     else  lagt <-c(0,vario$centert)
    }
    else {lagt=0}
    numlagt <- length(lagt)
    }

    # Check profile indices.  They index the empirical lag grids (including
    # the zero-lag row/column), while the fitted curve is matched by lag value.
    .GeoCovariogram_check_profile_index <- function(x, upper, label) {
      if(is.null(x)) return(NULL)
      if(!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
         x < 1 || x > upper || x != as.integer(x))
        stop(label, " must be a single integer in 1:", upper, ".", call. = FALSE)
      as.integer(x)
    }
    fix.lags <- .GeoCovariogram_check_profile_index(fix.lags, numlags, "fix.lags")
    fix.lagt <- .GeoCovariogram_check_profile_index(fix.lagt, numlagt, "fix.lagt")
    # END ---- check input --- #
    # set range (for practical range)
    if(tolower(fitted$corrmodel)=="matern"){
        lower <- 1e-10
        upper <- 1e20}
    else{
        lower <- 0  
        upper <- 1e100}

    num_betas <- fitted$numbetas
    if(is.null(num_betas) || !length(num_betas)) num_betas <- 1L
    mu=0;nuisance=0
    mm=c(mean=0)
    ## selecting nuisance, correlation parameters and a representative mean
      if(!bivariate){
        allpar <- c(fitted$fixed, fitted$param)
        param <- allpar[CorrelationPar(CkCorrModel(fitted$corrmodel))]
        nuisance <- allpar[NuisParam2(fitted$model,FALSE,num_betas, copula = fit_copula)]
        mean_names <- .GeoMean_names(num_betas)
        nuisance <- nuisance[!(names(nuisance) %in% mean_names)]

        nobs_mean <- length(unlist(fitted$data, use.names = FALSE))
        external_mean <- external_mean_fit
        if(!is.null(external_mean)) {
          mean_reference <- mean(external_mean)
        } else {
          beta <- .GeoMean_beta(allpar, num_betas)
          Xfit <- fitted$X
          ## Backward compatibility for GeoWLS objects created before fix16:
          ## their intercept-only design could be stored as the scalar 1.
          ## Treat that legacy sentinel as the default intercept design rather
          ## than passing an ambiguous vector to .GeoMean_design().
          if(inherits(fitted, "GeoWLS") &&
             is.atomic(Xfit) && is.null(dim(Xfit)) && length(Xfit)==1L &&
             length(num_betas)==1L && num_betas==1L) Xfit <- NULL
          Xobs <- .GeoMean_design(Xfit, n = nobs_mean, p = num_betas, name = "fitted$X")
          mean_reference <- mean(as.numeric(Xobs %*% beta))
        }
        mm <- c(mean = mean_reference)
        }
      if(bivariate){
        param <- c(fitted$fixed,fitted$param)[CorrelationPar(CkCorrModel(fitted$corrmodel))]
        nuisance <- c(fitted$fixed,fitted$param)[NuisParam2(fitted$model,FALSE,num_betas)]
    }
    
     #############################################
     # computing the spatio-temporal distances where to compute the fitted model
     #############################################
    if(isvario) {

             if(!invisible) {   #### standard case
                              lags_m <- seq(slow,max(vario$centers),length.out =150)
                              if(ispatim)  lagt_m <-seq(slow,max(vario$bint),length.out =150)
                              else         lagt_m<-0
                            }
             else           {   ### invisible case
                                if (ispatim) {lags_m=c(slow,vario$centers)
                                              lagt_m <-c(slow,vario$bint)}
                                else         {lags_m=vario$centers
                                             lagt_m<-0}
                            }    
     }
    else{
        mmm <- double(2)
        type_dist <- CheckDistance(distance)
        p=.C("Maxima_Minima_dist",mmm=as.double(mmm),as.double(fitted$coordx),as.double(fitted$coordy),as.double(fitted$coordz),
        as.integer(fitted$numcoord),as.integer(type_dist),as.double(fitted$radius),
        PACKAGE='GeoModels', NAOK=TRUE)
        lags_m <- seq(slow,p$mmm[2],length=150)
        if (ispatim) {
            tt <- double(2)
            lagt_m <- seq(slow,max(c(dist(fitted$coordt))) ,length=150)
        }
        else lagt_m<-0
    }
    numlags_m <- length(lags_m)
    numlagt_m <- length(lagt_m)
if(!bivariate) {
    ## A covariogram indexed only by lag cannot retain a location-varying
    ## mean.  Use the average fitted linear predictor as representative.
    mu <- rep(as.numeric(mm["mean"]), numlags_m * numlagt_m)
  } else{}  
   ###############
    .st_old <- .GeoSTEnter(fitted$corrmodel,
      coordx = fitted$coordx, coordy = fitted$coordy, coordz = fitted$coordz,
      coordx_dyn = fitted$coordx_dyn,
      distance = if (is.null(fitted$distance)) distance else fitted$distance,
      context = "GeoCovariogram")
    on.exit(.GeoSTRestore(.st_old), add = TRUE)

    corrmodel <- CkCorrModel(fitted$corrmodel)

     nui=nuisance
     nui['sill']=1;
if(!(binomial||geom||binomialneg||binomialnegZINB||Gaussian_misp_Binomial||
       Gaussian_misp_BinomialNeg||Gaussian_misp_Poisson||Gaussian_misp_PoissonGamma)) nui['nugget']=0
else                                     nui['nugget']=nuisance['nugget']

    if(fitted$model=="Gaussian_misp_Binomial") fitted$model="Binomial"
    if(fitted$model=="Gaussian_misp_BinomialNeg") fitted$model="BinomialNeg"
    if(fitted$model=="Gaussian_misp_Poisson") fitted$model="Poisson"
    if(fitted$model=="Gaussian_misp_PoissonGamma") fitted$model="PoissonGamma"
  
    ## Copula covariances are functions of the correlation of the underlying
    ## Gaussian RF, not of a margin-specific transformed correlation.  In
    ## particular, using the Beta/Kumaraswamy model id here would apply an
    ## additional transformation before the copula covariance engine and make
    ## GeoCovariogram inconsistent with GeoCovmatrix/GeoKrig.
    ## For the zero-inflated models GeoCovariogram needs the correlation of
    ## the underlying latent Gaussian field and constructs the observed-scale
    ## covariance below.  VectCorrelation has no dedicated latent-only branch
    ## for all ZIP model codes (notably 53, 57 and 58), so request the Gaussian
    ## correlation explicitly rather than passing an unsupported ZIP code.
    correlation_model <- if(copula_cont || poissonZIP || poissongammaZIP)
                           CkModel("Gaussian") else CkModel(fitted$model)
    correlation <- CorrelationFct(bivariate,corrmodel, lags_m, lagt_m, numlags_m, numlagt_m,mu,
                                     correlation_model, nui,param,fitted_n)
        #correlation <- CorrelationFct(bivariate,corrmodel, lags_m, lagt_m, numlags_m, numlagt_m,mu,
         #                            CkModel(fitted$model), nui,param,fitted$n)

iscopula <- identical(fit_copula, "Gaussian")

    ## For the supported continuous copula margins, compute the theoretical
    ## covariance directly on the observed scale with exactly the same engine
    ## used by GeoCovmatrix() and GeoKrig().  A lag-only curve cannot retain a
    ## location-varying marginal predictor, so eta_ref is the representative
    ## mean predictor computed above.
    copula_variance <- NULL
    if(copula_cont) {
      eta_ref <- as.numeric(mm["mean"])[1L]
      if(!is.finite(eta_ref)) eta_ref <- 0
      nugget_cop <- suppressWarnings(as.numeric(nuisance["nugget"]))
      if(!length(nugget_cop) || !is.finite(nugget_cop[1L])) nugget_cop <- 0
      nugget_cop <- nugget_cop[1L]
      if(nugget_cop < 0 || nugget_cop >= 1)
        stop("GeoCovariogram requires a copula nugget in [0,1).", call. = FALSE)

      rho_cop <- (1 - nugget_cop) * as.numeric(correlation)
      params_cop <- c(mm, nuisance)
      covariance <- .GeoCopulaCrossCovFast(
        fit_copula, rho_cop, model, params_cop,
        eta1 = rep(eta_ref, length(rho_cop)),
        eta2 = rep(eta_ref, length(rho_cop))
      )
      copula_variance <- as.numeric(
        .GeoCopulaMarginalMoments(model, params_cop, eta_ref)$variance
      )[1L]
      if(!is.finite(copula_variance) || copula_variance <= 0)
        stop("Invalid marginal variance in copula covariogram.", call. = FALSE)
      variogram <- copula_variance - covariance
    }

##########################################
########### starting legacy cases ########
##########################################
if(!copula_cont) {
    if(gaussian){
    llm=length(lags_m)
        if(bivariate){  covariance11 <- correlation[(0*llm+1):(1*llm)]
                        covariance12 <- correlation[(1*llm+1):(2*llm)]
                        covariance22 <- correlation[(3*llm+1):(4*llm)]
                        variogram11  <- correlation[(4*llm+1):(5*llm)]
                        variogram12  <- correlation[(5*llm+1):(6*llm)]
                        variogram22  <- correlation[(7*llm+1):(8*llm)]
                           }
        else { 

        #covariance <- nuisance["nugget"]+nuisance["sill"]*correlation
        covariance <- as.numeric(nuisance["sill"])*correlation*(1-as.numeric(nuisance["nugget"]))
        #variogram <- nuisance["nugget"]+nuisance["sill"]*(1-correlation)
        variogram <-as.numeric(nuisance["sill"])*(1-correlation*(1-as.numeric(nuisance["nugget"])))
       
        }
    }
##########################################
 if(skewgausssian) {    
            if(bivariate) {}
              else {
        if(!iscopula){
              correlation1=(1-as.numeric(nuisance['nugget']) )*correlation  
              vv=as.numeric(nuisance['sill']);sk=as.numeric(nuisance['skew']);sk2=sk^2;corr2=correlation^2;  
              cc=(2*sk2)*(sqrt(1-corr2) + correlation*asin(correlation)-1)/(pi*vv+sk2*(pi-2)) + (correlation1*vv)/(vv+sk2*(1-2/pi))
              vs=(vv+sk2*(1-2/pi))
              covariance=vs*cc;variogram=vs*(1-cc) }
              else
              {
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)   )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
              }
          }
    }
##########################################
 if(skewlaplace) {    
            if(bivariate) {}
              else {
        if(!iscopula){
              correlation1=(1-as.numeric(nuisance['nugget']) )*correlation  
              vv=as.numeric(nuisance['sill']);sk=as.numeric(nuisance['skew']);sk2=sk*sk;corr2=correlation1^2;  
              cc=corr2
              vs= vv*(1-2*sk+2*sk2)/(sk2*(1-sk)^2)
              covariance=vs*cc;variogram=vs*(1-cc) }
              else
              {
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)   )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
              }
          }
    }
##########################################
   if(twopieceT)        { if(bivariate) {}
                        else {
                        if(!iscopula){
                              correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))
                              nu=1/as.numeric(nuisance['df']); sk=as.numeric(nuisance['skew']);sill=as.numeric(nuisance['sill'])
                              sk2=sk^2
                              vs=sill* ((nu/(nu-2))*(1+3*sk2) - 4*sk2*(nu/pi)*(gamma(0.5*(nu-1))/gamma(0.5*nu))^2)
                              corr2=correlation1^2;sk2=sk^2
                                  a1=Re(hypergeo::hypergeo(0.5,0.5,nu/2,corr2))
                                  a2=correlation1 *asin(correlation1) + (1-corr2)^(0.5)
                                  ll=qnorm((1-sk)/2)
                                  p11=pbivnorm::pbivnorm(ll,ll, rho = correlation, recycle = TRUE)
                                  a3=3*sk2 + 2*sk + 4*p11 - 1
                                  KK=( nu*(nu-2)*gamma((nu-1)/2)^2) / (nu*pi*gamma(nu/2)^2*(3*sk2+1)-4*sk2*nu*(nu-2)*gamma((nu-1)/2)^2 )
                                  cc= KK*(a1*a2*a3-4*sk2);
                              ##
                              covariance=vs*cc;variogram=vs*(1-cc)  
                          }
                                else
              {
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)  )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
              }

                          }
                  } 
  ##########################################
   if(twopieceTukeyh)        { if(bivariate) {}
                        else {
                          if(!iscopula){     
                              correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))
                              tail=as.numeric(nuisance['tail']); sk=as.numeric(nuisance['skew']);
                              sill=as.numeric(nuisance['sill'])
                              corr2=correlation1^2;sk2=sk^2
                              gg2=(1-(1-corr2)*tail)^2
                              xx=corr2/gg2
                              A=(asin(sqrt(xx))*sqrt(xx)+sqrt(1-xx))/(1-xx)^(1.5)
                              ll=qnorm((1-sk)/2)
                              p11=pbivnorm::pbivnorm(ll,ll, rho = correlation, recycle = TRUE)
                              a3=3*sk2 + 2*sk + 4*p11 - 1
                              mm=8*sk2/(pi*(1-tail)^2); 
                              ff=(1+3*sk2)/(1-2*tail)^(1.5)
                              M=(2*(1-corr2)^(3/2))/(pi*gg2)
                              cc=  (M*A*a3-mm)/( ff- mm)
                              vs=sill*(ff- mm) 
                              covariance=vs*cc;variogram=vs*(1-cc) 
                                  }
                                else
              {
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)   )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
                   }

                          }
                  } 
##########################################
   if(twopieceGauss)        { 
                        if(bivariate) {}
                        else {     
             if(!iscopula){     
                        correlation1=correlation*(1-as.numeric(nuisance['nugget']))
                        corr2=sqrt(1-correlation1^(2))
                        sk=as.numeric(nuisance['skew']); sk2=sk^2
                        ll=qnorm((1-sk)/2)
                        p11=pbivnorm::pbivnorm(ll,ll, rho = correlation, recycle = TRUE)
                        KK=3*sk2+2*sk+ 4*p11 - 1
                        cc=(2*((corr2 + correlation1*asin(correlation1))*KK)- 8*sk2)/(3*pi*sk2  -  8*sk2   +pi   )
                        vs= as.numeric(nuisance['sill'])*(1+3*sk2-8*sk2/pi)
                        covariance=vs*cc;variogram=vs*(1-cc) 
                          }
                                else
              {
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)  )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
                        }
                  } 
              }
##########################################
 if(twopiecebimodal)        { if(bivariate) {}
                        else {      
                            if(!iscopula){                         
                                  correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))
                                  nu=as.numeric(nuisance['df']); sk=as.numeric(nuisance['skew'])
                                  delta=as.numeric(nuisance['shape'])
                                  alpha=2*(delta+1)/nu
                                  nn=2^(1-alpha/2)
                                  ll=qnorm((1-sk)/2)
                                  p11=pbivnorm::pbivnorm(ll,ll, rho = correlation, recycle = TRUE)
                                  corr2=correlation1^2;sk2=sk^2
                                  a1=Re(hypergeo::hypergeo(-1/alpha ,-1/alpha,nu/2,corr2))
                                  a3=3*sk2 + 2*sk + 4*p11 - 1
                                  MM=(2^(2/alpha)*(gamma(nu/2 + 1/alpha))^2) 
                                  vari=2^(2/alpha)*(gamma(nu/2 + 2/alpha))*gamma(nu/2)* (1+3*sk2) - sk2*2^(2/alpha+2)*gamma(nu/2+1/alpha)^2 
                                  cc= MM*(a1*a3-4*sk2)/vari
                                  vs=as.numeric(nuisance['sill'])*vari/(nn^(2/alpha)*gamma(nu/2)^2)
                                  covariance=vs*cc; variogram=vs*(1-cc)
                                         }
                                else
              {
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)   )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
                        }
                  } 
              }
##########################################

   if(studentT)        { if(bivariate) {}
                        else {
                                if(!iscopula){   
                              correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))
                              nu=1/as.numeric(nuisance['df']);sill=as.numeric(nuisance['sill'])
      
                              vs=sill/(1-2/nu)
                              cc=.GeoStudentCorrelation(correlation, correlation1, nu)
                              covariance=vs*cc;variogram=vs*(1-cc) 
                                            }
                                else
              {

                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)  )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
                        }
                  } 
              }
##########################################  
   if(skewstudentT)        { if(bivariate) {}
                        else {
                                if(!iscopula){   
                               correlation1=correlation*(1-nuisance['nugget'] )
                               nu=as.numeric(1/nuisance['df']); sk=as.numeric(nuisance['skew'])
                               sill=as.numeric(nuisance['sill'])
                               cc <- .GeoSkewStudentCorrelation(correlation, correlation1, nu, sk)
                               vs <- sill * .GeoSkewStudentMoments(nu, sk)$variance
                               covariance=vs*cc;variogram=vs*(1-cc)
                                            }
                                else
              {
               
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)  )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
                        }
                  } 
              }
##########################################
  if(tukeyh)        { if(bivariate) {}
                        else {
                                if(!iscopula){   
                              correlation=correlation*(1-nuisance['nugget'] )
                              h=as.numeric(nuisance['tail'])
                              sill=as.numeric(nuisance['sill'])
                              vs=  (1-2*h)^(-1.5)     
                              cc=correlation*(1-2*h)^(1.5)/ ((1-h)^2-(h*correlation)^2)^(1.5)
                              covariance=sill*vs*cc;variogram=sill*vs*(1-cc)  
                                         }
                                else
              {

                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)  )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
                        }
                  } 
              }

if (tukeyh2) { 
  if (bivariate) {
    # (bivariate case not handled here)
  } else {
    if (!iscopula) {
      nugget <- as.numeric(nuisance['nugget'])
      hr <- as.numeric(nuisance['tail1'])
      hl <- as.numeric(nuisance['tail2'])
      sill <- as.numeric(nuisance['sill'])
      mom_hh <- .GeoTukeyh2StandardMoments(hr, hl)
      rho_lat <- (1 - nugget) * as.numeric(correlation)
      cc <- .GeoTukeyh2Correlation(rho_lat, hr, hl)
      covariance <- sill * mom_hh$variance * cc
      variogram <- sill * mom_hh$variance * (1 - cc)
    } else {
      cova <- gaussian_copula_cov((1 - as.numeric(nuisance['nugget'])) * correlation,
                                  model, c(mm, nuisance))
      vv <- variance_disp(model,nuisance)
      covariance <- cova
      variogram <- vv - cova
    }
  }
}

##########################################
  if(tukeygh)        { if(bivariate) {}
                        else {
                     if(!iscopula){   
                              correlation=correlation*(1-nuisance['nugget'] )
                              rho=correlation
                              tail=as.numeric(nuisance['tail'])
                              sill=as.numeric(nuisance['sill'])
                              eta=as.numeric(nuisance['skew'])
                              rho2=rho*rho; eta2=eta*eta; tail2=tail*tail;

                        if(tail>1e-05&&abs(eta)>1e-05){
                               rho2=rho*rho; eta2=eta*eta; tail2=tail*tail;
                               u=1-tail; a=1+rho;
                               mu=(exp(eta2/(2*u))-1)/(eta*sqrt(u));
                               vs=(exp(2*eta2/(1-2*tail))-2*exp(eta2/(2*(1-2*tail)))+1)/(eta2*sqrt(1-2*tail))-mu*mu;
                               A1=exp(a*eta2/(1-tail*a));
                               A2=2*exp(0.5*eta2*  (1-tail*(1-rho2))  / (u*u- tail2*rho2)  );
                               A3=eta2*sqrt(u*u- rho2*tail*tail);
                               cc=((A1-A2+1)/A3-mu*mu)/vs;
                               covariance=sill*vs*cc;variogram=sill*vs*(1-cc)
                            } 
                        if(tail<=1e-05&&abs(eta)>1e-05){
                              vs=( -exp(eta^2)+exp(eta^2*2))*eta^(-2)
                              cc= (( -exp(eta^2)+exp(eta^2*(1+rho)))*eta^(-2))/vs
                              covariance=sill*vs*cc;variogram=sill*vs*(1-cc) 
                            } 
                        if(tail>1e-05&&abs(eta)<=1e-05){
                              vs=  (1-2*tail)^(-1.5)     
                              cc=(-rho/((1+h*(tail-1))*(-1+tail+tail*rho)*(1+tail*(-2+tail-tail*rho^2))^0.5))/vs
                              covariance=sill*vs*cc;variogram=sill*vs*(1-cc) 
                            } 
                        if(tail<=1e-05&&abs(eta)<=1e-05){
                            covariance=sill*rho;variogram=sill*(1-rho)
                            }  
                                
                     }
                                else
              {
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)  )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
                        }
                  } 
              }
######################################################
  if(sas) { if(bivariate) {}
                        else {
             if(!iscopula){    
nugget <- as.numeric(nuisance['nugget'])  # in [0,1]
d      <- as.numeric(nuisance['tail'])
e      <- as.numeric(nuisance['skew'])
sill   <- as.numeric(nuisance['sill'])
mom_sas <- .GeoSASStandardMoments(e, d)
rho_lat <- (1 - nugget) * as.numeric(correlation)
cc <- corrsas(rho_lat, e, d)
sigma2_marg <- sill * mom_sas$variance
covariance <- sigma2_marg * cc
variogram  <- sigma2_marg * (1 - cc)

                        }
                                else
              {
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)  )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
                        }
                  } 
              }
##########################################
 if(gamma)        { if(bivariate) {}
                        else {
                if(!iscopula){ 
                              correlation=correlation*(1-as.numeric(nuisance['nugget'] ))
                              vs=2*exp(mm)^2/as.numeric(nuisance['shape'])
                              cc=correlation^2
                              covariance=vs*cc;variogram=vs*(1-cc)  
                        }
                                else
              {
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)   )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
                        }
                  } 
              }
##########################################
 if(weibull)        { if(bivariate) {} 
                        else {
                if(!iscopula){        
                            ssh=as.numeric(nuisance['shape'])
                        correlation=correlation*(1-as.numeric(nuisance['nugget'] ))  
                        vs=exp(mm)^2*(gamma(1+2/ssh)/gamma(1+1/ssh)^2-1)
                        cc=.GeoWeibullCorr(correlation, ssh)
                        covariance=vs*cc;variogram=vs*(1-cc) 
                      }
                                else
              {
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)  )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
                        }
                  } 
              }
##########################################
  if(loglogistic)    { if(bivariate) {}  
                      else { 
              if(!iscopula){   
                     correlation=correlation*(1-as.numeric(nuisance['nugget'] ))
                     sh=as.numeric(nuisance["shape"])
                     vs=exp(mm)^2*(2*sh*sin(pi/sh)^2/(pi*sin(2*pi/sh))-1)
                     cc=((pi*sin(2*pi/sh))/(2*sh*(sin(pi/sh))^2-pi*sin(2*pi/sh)))*
                                    (Re(hypergeo::hypergeo(-1/sh, -1/sh, 1,correlation^2))*
                                     Re(hypergeo::hypergeo( 1/sh,  1/sh, 1,correlation^2)) -1)
                      covariance=vs*cc;variogram=vs*(1-cc)  

                }
                                else
              {
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)  )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
                        }
                  } 
              }
##########################################
  if(loggauss)    { if(bivariate) {}  
                      else {  
             if(!iscopula){    
                    vvar=as.numeric(nuisance["sill"])
                    yy=vvar*correlation*(1-as.numeric(nuisance["nugget"]))
                    cc<-(exp(yy)-1)/(exp(vvar)-1)
                    eta_ref <- as.numeric(mm["mean"])[1L]
                    vs<-exp(2*eta_ref)*(exp(vvar)-1)
                    covariance=vs*cc;variogram=vs*(1-cc)   
                    
                     }
                                else
              {
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)   )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
                        }
                  } 
              }
   #                   
   if(binary||binomial||binomial2||geom||binomialneg||binomialnegZINB||Gaussian_misp_Binomial||Gaussian_misp_BinomialNeg) {
                    if(bivariate) {}
                    else {      
                    if(!iscopula){      
                           pp=pnorm(mu)
                           if(binary||binomial||binomial2) vv=min(fitted$n)*pp*(1-pp)
                           if(geom)             vv=(1-pp)/pp^2;
                           if(binomialneg)      vv=(fitted$n)*(1-pp)/pp^2;
                           if(binomialnegZINB) { 
                                     pg=pnorm(as.numeric(nuisance['pmu']))
                                     vv=fitted$n*(1-pp)*(1-pg)*(1+fitted$n*pg*(1-pp)) /pp^2
                                               }
                           if(Gaussian_misp_Binomial||Gaussian_misp_BinomialNeg) vv=1                    
                           covariance=vv*correlation
                           variogram=vv*(1-correlation)
                          }
                                else
              {
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)  )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
                        }
                  } 
              }
###########################################################################              
     if(poisson||Gaussian_misp_Poisson) {
                    if(bivariate) {}
                    else {   
                if(!iscopula){  
                           correlation=(1-as.numeric(nuisance['nugget']))*correlation   
                           corr2=correlation^2    
                           vv=exp(mu);
                           z=2*vv/(1-corr2)
                           cc=corr2*(1-(besselI(z,0,expon.scaled = TRUE)+besselI(z,1,expon.scaled = TRUE)))
                           if(Gaussian_misp_Poisson) vv=1   
                           covariance=vv*cc
                           variogram=vv*(1-cc)
                    }
                                else
              {
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)   )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
                        }
                  } 
              }
########################################################################### 
   if(poissongamma||Gaussian_misp_PoissonGamma) {

                    if(bivariate) {}
                    else {
                             if(!iscopula){
                       ## Correct stationary PEM / Poisson-Gamma correlation
                       ## (Corollary 5.1).  Keep GeoCovariogram on the same
                       ## compiled kernel used by GeoCorrFct/GeoCovmatrix,
                       ## rather than the obsolete R-level hypergeometric
                       ## representation previously used here.
                       nug=as.numeric(nuisance['nugget'])
                       a=as.numeric(nuisance['shape'])
                       mean_pg=exp(mu)

                       ## VectCorrelation() returns the latent correlation for
                       ## model 46/47.  The nugget acts on that latent field;
                       ## at the same space-time point the response correlation
                       ## is exactly one.
                       if(ispatim) {
                         zero_lag=rep(lags_m,times=numlagt_m)==0 &
                                  rep(lagt_m,each=numlags_m)==0
                       }
                       else zero_lag=rep(lags_m==0,times=numlagt_m)
                       latent_corr=(1-nug)*correlation + nug*zero_lag

                       cc=.Call("GeoCorrPGsVecCall",
                         as.double(latent_corr),as.double(mean_pg[1L]),
                         as.double(a),PACKAGE="GeoModels")
                       if(any(zero_lag)) cc[zero_lag]=1

                       ## Var{M(s)} = mu * (1 + mu/kappa).
                       vv=mean_pg*(1+mean_pg/a)
                       covariance=vv*cc
                       variogram=vv*(1-cc)
                           }
                                else
              {
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)   )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
                        }
                  } 
              }
########################################################################### 
     if(poissonZIP) {
                    if(bivariate) {}
               else {
                  if(!iscopula){
                       ## Y(s)=I(s)X(s), where X is the renewal-Poisson count
                       ## field and I is the non-structural-zero indicator.
                       ## Hence
                       ## Cov(Y_i,Y_j)=p11 Cov(X_i,X_j)
                       ##              + mu_i mu_j {p11-(1-p)^2}.
                       ## The previous code used the ZIP variance, rather than
                       ## the Poisson mean, inside the Bessel correlation and
                       ## did not assign `covariance` at all.
                       p=pnorm(as.numeric(nuisance['pmu']))
                       q=1-p
                       MM=exp(mu)

                       if(ispatim) {
                         zero_lag=rep(lags_m,times=numlagt_m)==0 &
                                  rep(lagt_m,each=numlags_m)==0
                       } else zero_lag=rep(lags_m==0,times=numlagt_m)

                       ## PoissonZIP1 (model 53) uses one shared nugget; the
                       ## standard PoissonZIP/misspecified version uses two.
                       if(model==53) {
                         nug_count=as.numeric(nuisance['nugget'])
                         nug_zero=nug_count
                       } else {
                         nug_count=as.numeric(nuisance['nugget1'])
                         nug_zero=as.numeric(nuisance['nugget2'])
                       }
                       rho_count=(1-nug_count)*correlation+nug_count*zero_lag
                       rho_zero=(1-nug_zero)*correlation+nug_zero*zero_lag

                       p11=numeric(length(correlation))
                       p11[zero_lag]=q
                       jj=!zero_lag
                       if(any(jj))
                         p11[jj]=1-2*p+pbivnorm::pbivnorm(
                           as.numeric(nuisance['pmu']),
                           as.numeric(nuisance['pmu']),
                           rho=rho_zero[jj],recycle=TRUE)

                       ## Stationary underlying Poisson correlation.  At zero
                       ## lag impose the exact endpoint and avoid 1-rho^2=0.
                       cc_count=numeric(length(correlation))
                       cc_count[zero_lag]=1
                       if(any(jj)) {
                         corr2=rho_count[jj]^2
                         z=2*MM[jj]/(1-corr2)
                         cc_count[jj]=corr2*(1-(
                           besselI(z,0,expon.scaled=TRUE)+
                           besselI(z,1,expon.scaled=TRUE)))
                       }
                       cov_count=MM*cc_count
                       covariance=p11*cov_count+MM^2*(p11-q^2)
                       vv=q*MM*(1+p*MM)
                       variogram=vv-covariance
            }
                                else
              {
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)   )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
                        }
                  }
              }
       if(poissongammaZIP) {
                    if(bivariate) {}
                    else {
                    if(!iscopula){
                       ## Zero-inflated PEM / Poisson-Gamma covariance.  Use
                       ## the corrected stationary PEM kernel for the count
                       ## component, then apply the exact zero-inflation
                       ## covariance identity.  This removes the obsolete
                       ## R-level hypergeometric expression from this branch.
                       p=pnorm(as.numeric(nuisance['pmu']))
                       q=1-p
                       MM=exp(mu)
                       a=as.numeric(nuisance['shape'])

                       if(ispatim) {
                         zero_lag=rep(lags_m,times=numlagt_m)==0 &
                                  rep(lagt_m,each=numlags_m)==0
                       } else zero_lag=rep(lags_m==0,times=numlagt_m)

                       ## PoissonGammaZIP1 (model 58) has a shared nugget;
                       ## model 57 has separate count/zero-process nuggets.
                       if(model==58) {
                         nug_count=as.numeric(nuisance['nugget'])
                         nug_zero=nug_count
                       } else {
                         nug_count=as.numeric(nuisance['nugget1'])
                         nug_zero=as.numeric(nuisance['nugget2'])
                       }
                       rho_count=(1-nug_count)*correlation+nug_count*zero_lag
                       rho_zero=(1-nug_zero)*correlation+nug_zero*zero_lag

                       p11=numeric(length(correlation))
                       p11[zero_lag]=q
                       jj=!zero_lag
                       if(any(jj))
                         p11[jj]=1-2*p+pbivnorm::pbivnorm(
                           as.numeric(nuisance['pmu']),
                           as.numeric(nuisance['pmu']),
                           rho=rho_zero[jj],recycle=TRUE)

                       cc_count=.Call("GeoCorrPGsVecCall",
                         as.double(rho_count),as.double(MM[1L]),
                         as.double(a),PACKAGE="GeoModels")
                       if(any(zero_lag)) cc_count[zero_lag]=1
                       var_count=MM*(1+MM/a)
                       cov_count=var_count*cc_count

                       covariance=p11*cov_count+MM^2*(p11-q^2)
                       vv=q*MM*(1+MM/a+p*MM)
                       variogram=vv-covariance
                 }
                                else
              {
                cova=gaussian_copula_cov((1-as.numeric(nuisance['nugget']))*correlation, model,c(mm,nuisance)   )
                vv=variance_disp(model,nuisance)
                covariance=cova
                variogram=vv-cova
                        }
                  }
              }



}

##########################################
############ starting graphics############
##########################################
      vario.main <- "Spatial semivariogram"
      vario.ylab <- "SemiVariogram"
        if(ispatim){
            dim(covariance) <- c(numlags_m, numlagt_m)
            dim(variogram) <- c(numlags_m, numlagt_m)
            vario.main <- "Space-time semivariogram"
            vario.zlab <- "SemiVariogram" }
        # Compute the practical range on the observed marginal scale.  For
        # Gaussian models retain the exact correlation-root calculation.  For
        # non-Gaussian models use the fitted covariance/variogram themselves,
        # so the definition is consistent with GeoCovmatrix/GeoKrig.
        Range <- NA_real_
        if(show.range || answer.range){
            if(bivariate) {
              Range <- NA_real_
            } else if(gaussian && !copula_cont) {
              Range <- tryCatch(
                uniroot(PracRangeNorm, c(lower, upper), corrmodel=corrmodel,
                        nuisance=nuisance, numlags=1, numlagt=1, lagt=0,
                        param=param, pract.range=pract.range)$root,
                error = function(e) NA_real_
              )
            } else {
              if(ispatim) {
                jt0 <- .GeoCovariogram_nearest_index(lagt_m, 0)
                vcurve <- variogram[, jt0]
                ccurve <- covariance[, jt0]
              } else {
                vcurve <- variogram
                ccurve <- covariance
              }
              Range <- .GeoCovariogram_practical_range(
                lags_m, vcurve, ccurve, pract.range
              )
            }
        }

        ## Sum-of-squares diagnostic.  Compute it independently of plotting so
        ## invisible=TRUE is useful with show.vario=FALSE as documented.
        OLS <- NULL
        if(isTRUE(invisible)) {
          if(bivariate) {
            OLS <- .GeoCovariogram_sqerr(vario$variograms[1,], variogram11) +
                   .GeoCovariogram_sqerr(vario$variograms[2,], variogram22)
            if(!dyn)
              OLS <- OLS + .GeoCovariogram_sqerr(vario$variogramst, variogram12)
          } else if(ispatim) {
            nbint_ols <- if(is.null(vario$centert)) length(vario$bint) else length(vario$centert)
            empirical_st <- matrix(vario$variogramst,
                                   nrow=length(vario$centers),
                                   ncol=nbint_ols, byrow=TRUE)
            empirical_st <- rbind(c(0,vario$variogramt),
                                  cbind(vario$variograms, empirical_st))
            OLS <- .GeoCovariogram_sqerr(empirical_st, variogram)
          } else {
            OLS <- .GeoCovariogram_sqerr(vario$variograms, variogram)
          }
        }



if(show.cov){

#####################################
#### bivariate case covariance ######
##################################### 
        if(bivariate&&!dyn){
            #par(mfrow=c(2,2))
       plot(lags_m, covariance11, type='l', ylim=c(min(covariance11),
                     max(covariance11)), main="First covariance",
                     xlab="Distance", ylab="Covariance",...)
       plot(lags_m, covariance12, type='l', ylim=c(min(covariance12),
                     max(covariance12)), main="Cross covariance",
                     xlab="Distance", ylab="Covariance",...)      
       plot(lags_m, covariance12, type='l', ylim=c(min(covariance12),
                     max(covariance12)), main="Cross covariance",
                     xlab="Distance", ylab="Covariance",...)   
       plot(lags_m, covariance22, type='l', ylim=c(min(covariance22),
                     max(covariance22)), main="Second covariance",
                     xlab="Distance", ylab="Covariance",...)          
         }
     if(bivariate&&dyn){
            #par(mfrow=c(2,2))
       plot(lags_m, covariance11, type='l', ylim=c(min(covariance11),
                     max(covariance11)), main="First covariance",
                     xlab="Distance", ylab="Covariance",...)
       plot(lags_m, covariance22, type='l', ylim=c(min(covariance22),
                     max(covariance22)), main="Second covariance",
                     xlab="Distance", ylab="Covariance",...)          
         }

        
#######################################
#### space time case  covariance ######
#######################################
if(ispatim){
            # build the covariance matrix:
            plagt <- !is.null(fix.lags)
            plags <- !is.null(fix.lagt)
            numplot <- 1+plags+plagt
          
            par(mai=c(.2,.2,.2,.2))
            persp(lags_m, lagt_m, covariance, xlab="Distance", ylab="Time",
                  zlab="Covariance", ltheta=90,
                  shade=0.75, ticktype="detailed", phi=30,
                  theta=30,main="Fitted space-time covariance",
                  cex.axis=.8, cex.lab=.8)
            if(plagt){
                fixed_space <- lags[fix.lags]
                ii_fit <- .GeoCovariogram_nearest_index(lags_m, fixed_space)
                par(mai=c(.5,.5,.5,.5),mgp=c(1.6,.6,0))
                plot(lagt_m, covariance[ii_fit,], xlab="Time",
                     ylab="Covariance", type="l",cex.axis=.8,cex.lab=.8,
                     main=paste0("Space-time cov: temporal profile (h=",
                                 signif(fixed_space, 4), ")"),...)}
            if(plags){
                fixed_time <- lagt[fix.lagt]
                jj_fit <- .GeoCovariogram_nearest_index(lagt_m, fixed_time)
                par(mai=c(.5,.5,.5,.5),mgp=c(1.6,.6,0))
                plot(lags_m, covariance[,jj_fit], xlab="Distance",
                     ylab="Covariance", type="l",cex.axis=.8,cex.lab=.8,
                     main=paste0("Space-time cov: spatial profile (u=",
                                 signif(fixed_time, 4), ")"),...)
            }
        }
#######################################
#### spatial case  covariance #########
#######################################
       if(!ispatim && !bivariate){# spatial case:
            if(add.cov & dev.cur()!=1){
                lines(lags_m, covariance,...)
                if(show.range && is.finite(Range)) abline(v=Range)}
            else{
                plot(lags_m, covariance, type='l', #ylim=c(0,max(covariance)), 
                     main="Spatial covariance",
                     xlab="Distance", ylab="Covariance",...)
                if(show.range && is.finite(Range)) abline(v=Range)}}
         }
########################### 
########################### 
########################### 
    # display the semivariogram function
    if(show.vario){

#####################################
#### bivariate case semivariogram ###
##################################### 
 if(bivariate){
      if(bivariate&&!dyn){
       plot(vario$centers,vario$variograms[1,], main="First semivariogram",ylim=c(0,max(vario$variograms[1,])),
           xlim=c(0,max(vario$centers)),
                     xlab="Distance", ylab="SemiVariogram",...)
       lines(lags_m, variogram11, type='l',...)
       if(min(vario$variogramst)>0) {ll1=0;ll=max(vario$variogramst)}
       if(min(vario$variogramst)<0) {ll1=min(vario$variogramst);ll=-min(vario$variogramst)}
       plot(vario$centers,vario$variogramst, main="Cross semivariogram",ylim=c(ll1,ll),
         xlim=c(0,max(vario$centers)),
                     xlab="Distance", ylab="SemiVariogram",...)
       lines(lags_m, variogram12, type='l',...)
       plot(vario$centers,vario$variogramst, main="Cross semivariogram",ylim=c(ll1,ll),
         xlim=c(0,max(vario$centers)),
                     xlab="Distance", ylab="SemiVariogram",...)
       lines(lags_m, variogram12, type='l',...)
       plot(vario$centers,vario$variograms[2,], main="Second semivariogram",ylim=c(0,max(vario$variograms[2,])),
         xlim=c(0,max(vario$centers)),
                     xlab="Distance", ylab="SemiVariogram",...)
       lines(lags_m, variogram22, type='l',...)  

   }
 
   if(bivariate&&dyn){
       plot(vario$centers,vario$variograms[1,], main="First semivariogram",ylim=c(0,max(vario$variograms[1,])),
                     xlab="Distance", ylab="SemiVariogram",...)
       lines(lags_m, variogram11, type='l',...)
       plot(vario$centers,vario$variograms[2,], main="Second semivariogram",ylim=c(0,max(vario$variograms[2,])),
                     xlab="Distance", ylab="SemiVariogram",...)
       lines(lags_m, variogram22, type='l',...)  


        }
}
#######################################
#### space time case  semivariogram ###
#######################################

 if(ispatim){

 # spatio-temporal case:
            plagt <- !is.null(fix.lags)
            plags <- !is.null(fix.lagt)
            nrowp <- 1
            ncolp <- 1
            if(isvario){
                ncolp <- ncolp+1
                if(plags || plagt) nrowp <- nrowp+1}
            else ncolp <- ncolp+plags+plagt
            sup <- 0
            tup <- 0

            ## Profile bookkeeping.  For dynamic supports the empirical
            ## temporal margin at h=0 may be unavailable because no sites are
            ## collocated across distinct times.  In that special case, when
            ## fix.lags selects the zero-lag entry, GeoCovariogram falls back
            ## to the first positive-distance ST bin containing finite pairs.
            ## The fitted curve is evaluated at that bin center, while the
            ## title reports the actual empirical distance interval.
            profile_space_idx <- fix.lags
            profile_space_value <- if(plagt) lags[fix.lags] else NA_real_
            profile_space_fallback <- FALSE
            profile_space_lower <- NA_real_
            profile_space_upper <- NA_real_
            profile_space_first_bin <- FALSE

            if(isvario){
                nbins <- length(vario$centers)
                if(is.null(vario$centert)) nbint <- length(vario$bint)
                else nbint <- length(vario$centert)

                ## GeoVariogram stores the positive-distance/positive-time
                ## surface row-wise, with the spatial and temporal margins
                ## held separately.  Rebuild the rectangular empirical grid
                ## explicitly.  For dynamic supports gamma(0,u) is available
                ## only when locations are actually collocated across times;
                ## unavailable cells therefore remain NA rather than being
                ## silently represented as zero.
                est_cross <- matrix(vario$variogramst,
                                    nrow=nbins,ncol=nbint,byrow=TRUE)
                evario <- rbind(c(0,vario$variogramt),
                                cbind(vario$variograms,est_cross))
                evario.grid <- as.matrix(expand.grid(lags,lagt))
                if(nrow(evario) != length(lags) || ncol(evario) != length(lagt) ||
                   nrow(evario.grid) != length(evario))
                  stop("GeoCovariogram: inconsistent empirical space-time variogram grid.",
                       call.=FALSE)

                empirical_z <- as.numeric(evario)
                empirical_ok <- is.finite(empirical_z)
                if(any(empirical_ok) && requireNamespace("scatterplot3d", quietly=TRUE)) {
                  ## Vertical stems are useful for a fixed spatial support but
                  ## become a dense wall for dynamic supports.  Plot finite
                  ## binned values as points in the dynamic case instead.
                  scatterplot3d::scatterplot3d(
                    evario.grid[empirical_ok,1], evario.grid[empirical_ok,2],
                    empirical_z[empirical_ok],
                    type=if(dyn) "p" else "h",
                    pch=if(dyn) 20 else 1,
                    highlight.3d=!dyn,cex.axis=.7,cex.lab=.7,
                    main=if(dyn) "Empirical dynamic ST semivariogram" else paste("Empirical",vario.main),
                    xlab="Distance",ylab="Time",zlab=vario.zlab,
                    mar=c(2,2,2,2),mgp=c(0,0,0))
                } else {
                  plot.new()
                  title(if(dyn) "Empirical dynamic ST semivariogram" else paste("Empirical",vario.main))
                  if(any(empirical_ok)) {
                    graphics::text(0.5,0.5,"Install 'scatterplot3d' for the 3D empirical panel")
                  }
                }
                if(plagt) {
                  profile_info <- .GeoCovariogramDynamicTemporalProfile(
                    vario=vario,evario=evario,lags=lags,fix.lags=fix.lags,
                    dynamic=dyn)
                  profile_space_idx <- profile_info$index
                  profile_space_value <- profile_info$value
                  profile_space_fallback <- profile_info$fallback
                  profile_space_lower <- profile_info$lower
                  profile_space_upper <- profile_info$upper
                  profile_space_first_bin <- profile_info$first.bin
                }
                if(plagt) tup <- .GeoCovariogram_finite_max(evario[profile_space_idx,],0)
                if(plags) sup <- .GeoCovariogram_finite_max(evario[,fix.lagt],0)
            }

            par(mai=c(.2,.2,.2,.2))
            if(dyn) {
              persp(lags_m, lagt_m, variogram, xlab="Distance",
                    ylab="Time", zlab=vario.zlab, ltheta=90,
                    shade=0.35, col="grey85", border="grey45",
                    ticktype="detailed", phi=30, theta=30,
                    main="Fitted ST semivariogram", cex.axis=.8,cex.lab=.8)
            } else {
              persp(lags_m, lagt_m, variogram, xlab="Distance",
                    ylab="Time", zlab=vario.zlab, ltheta=90,
                    shade=0.75, ticktype="detailed", phi=30,
                    theta=30,main=vario.main, cex.axis=.8,cex.lab=.8)
            }

            vvv=nuisance["sill"]
            if(copula_cont) vvv <- copula_variance
            ########
            if(!copula_cont && gamma)         vvv=2*exp(mm["mean"])^2/as.numeric(nuisance["shape"])
            if(!copula_cont && weibull)       vvv=exp(mm["mean"])^2*(gamma(1+2/as.numeric(nuisance["shape"]))/gamma(1+1/as.numeric(nuisance["shape"]))^2-1)
            if(!copula_cont && loglogistic)   vvv=exp(mm["mean"])^2*
                               (2*as.numeric(nuisance['shape'])*sin(pi/nuisance['shape'])^2/(pi*sin(2*pi/nuisance['shape']))-1)
            if(!copula_cont && loggauss)      vvv=variance_disp(model,c(mm,nuisance))
            if(!copula_cont && binomial)      vvv=fitted$n*pnorm(mm['mean'])*(1-pnorm(mm['mean']))
            if(!copula_cont && poisson)       vvv=exp(mm['mean'])
            #if(poissongamma)  vvv=exp(mm['mean']*(1+1/nuisance["shape"]))
            if(!copula_cont && poissongamma)  {MM=exp(mm['mean']); vvv=MM*(1+MM/as.numeric(nuisance["shape"]))
        }


            if(!copula_cont && poissonZIP){        p=pnorm(nuisance['pmu']); MM=exp(mm['mean']); 
                              vvv=(1-p)*MM*(1+p*MM)}

            if(!copula_cont && poissongammaZIP){     p=pnorm(nuisance['pmu']); MM=exp(mm['mean']); kk=as.numeric(nuisance["shape"]);
                              vvv=(1-p)*MM*(1+p*MM+(MM/kk))}

            if(!copula_cont && geom)          vvv= (1-pnorm(mm['mean']))/pnorm(mm['mean'])^2
            if(!copula_cont && binomialneg)   vvv= fitted$n*(1-pnorm(mm['mean']))/pnorm(mm['mean'])^2
            if(!copula_cont && binomialnegZINB) {
                                  MM=pnorm(mm['mean']);pg=pnorm(nuisance['pmu'])
                                }
            if(!copula_cont && skewgausssian) vvv=variance_disp(model,c(mm,nuisance))
            if(!copula_cont && skewlaplace)   vvv=variance_disp(model,c(mm,nuisance))
            if(!copula_cont && studentT)      vvv=variance_disp(model,c(mm,nuisance))
            if(!copula_cont && tukeyh)        vvv=variance_disp(model,c(mm,nuisance))
            if(!copula_cont && sas)           { 
                              d <- as.numeric(nuisance['tail']); e <- as.numeric(nuisance['skew'])
                              MM <- sinh(e/d) * exp(0.25) * ( besselK(.25,(d+1)/(2*d)) + besselK(.25,(1-d)/(2*d)) ) / sqrt(8*pi)
                              vv_std <- cosh(2*e/d) * exp(0.25) * ( besselK(.25,(d+2)/(2*d)) + besselK(.25,(2-d)/(2*d)) ) / sqrt(32*pi) - 0.5 - MM^2
                              vvv <- as.numeric(nuisance['sill']) * vv_std   # <-- sill*vv
               }
            if(!copula_cont && tukeyh2)  {        hr=as.numeric(nuisance["tail1"]);hl=as.numeric(nuisance["tail2"]);
                                  mom_hh <- .GeoTukeyh2StandardMoments(hr, hl)
                                  vvv=as.numeric(nuisance["sill"]) * mom_hh$variance
                         }
            ########
            if(plagt){
                fixed_space <- profile_space_value
                ii_fit <- .GeoCovariogram_nearest_index(lags_m, fixed_space)
                par(mai=c(.5,.5,.3,.3),mgp=c(1.6,.6,0))

                temporal_title <- if(dyn && profile_space_fallback) {
                  if(is.finite(profile_space_upper)) {
                    ## Exact h=0 pairs are not included in the positive-
                    ## distance ST surface; they belong to the true temporal
                    ## margin when available.
                    if(profile_space_first_bin)
                      paste0("Temporal profile (0<h<",
                             signif(profile_space_upper,4), ")")
                    else
                      paste0("Temporal profile (",signif(profile_space_lower,4),
                             "<=h<",signif(profile_space_upper,4), ")")
                  } else {
                    paste0("Temporal profile (h~",signif(fixed_space,4),")")
                  }
                } else if(dyn) {
                  paste0("Temporal profile (h=",signif(fixed_space,4),")")
                } else {
                  paste0(vario.ylab, ": temporal profile (h=",
                         signif(fixed_space,4), ")")
                }

                plot(lagt_m, variogram[ii_fit,], xlab="Time",cex.axis=.8,cex.lab=.8,
                     ylab=vario.ylab, type="l", ylim=c(0,max(vvv,tup,na.rm=TRUE)),
                     main=temporal_title,...)
                if(isvario) {
                  emp_x <- lagt[-1]
                  emp_y <- evario[profile_space_idx,-1]
                  emp_ok <- is.finite(emp_x) & is.finite(emp_y)
                  if(any(emp_ok)) {
                    points(emp_x[emp_ok],emp_y[emp_ok],...)
                  } else if(dyn) {
                    msg <- if(abs(fixed_space) <= sqrt(.Machine$double.eps))
                      "No collocated dynamic pairs for empirical h=0 profile"
                    else "No empirical pairs in this dynamic profile"
                    mtext(msg,side=3,line=-1.15,cex=.62)
                  }
                }
               }

            if(plags){
                fixed_time <- lagt[fix.lagt]
                jj_fit <- .GeoCovariogram_nearest_index(lagt_m, fixed_time)
                par(mai=c(.5,.5,.3,.3),mgp=c(1.6,.6,0))
                plot(lags_m, variogram[,jj_fit], xlab="Distance",cex.axis=.8,cex.lab=.8,
                     ylab=vario.ylab, type="l", ylim=c(0,max(vvv,sup,na.rm=TRUE)),
                     main=if(dyn) paste0("Spatial profile (u=",signif(fixed_time,4),")")
                          else paste0(vario.ylab, ": spatial profile (u=",
                                      signif(fixed_time, 4), ")"),...)
                if(isvario) {
                  emp_x <- lags[-1]
                  emp_y <- evario[-1,fix.lagt]
                  emp_ok <- is.finite(emp_x) & is.finite(emp_y)
                  if(any(emp_ok)) points(emp_x[emp_ok],emp_y[emp_ok],...)
                  else if(dyn) mtext("No empirical pairs in this dynamic profile",
                                     side=3,line=-1.15,cex=.62)
                }
                }
  #################################

 } ####### end space time
          
################################
#### spaatial semivariogram ###
#################################
    if(!ispatim && !bivariate){

   #################################
  #################################
            if(add.vario & dev.cur()!=1){
                points(vario$centers, vario$variograms,...)
                lines(lags_m, variogram,...)
                if(show.range && is.finite(Range)) abline(v=Range)}
            else{
                bnds <- range(variogram)
                bnds[1] <- min(bnds[1], min(vario$variograms))
                bnds[2] <- max(bnds[2], max(vario$variograms))
        
                plot(lags_m, variogram, type='l',  
                     main=vario.main,xlab="Distance",
                     ylab=vario.ylab,...)
                points(vario$centers, vario$variograms,...)
                if(show.range && is.finite(Range)) abline(v=Range)}
        }
 #################################       
}
    if(do_plot && ispatim) par(mai=c(1.02 ,0.85 ,0.85 ,0.45),mgp=c(3,1,0))
    # return the estimated covariance function
    if(answer.cov) {
      if(!bivariate) {
        result <- list(lags=lags_m,lagt=lagt_m, covariance=covariance)
      } else if(gaussian) {
        result <- list(lags=lags_m,lagt=lagt_m,
                       covariance11=covariance11, covariance12=covariance12,
                       covariance22=covariance22)
      }
    }
    # return the estimated variogram function
    if(answer.vario) {
        if(!is.list(result)) {
            if(!bivariate) result <- list(lags=lags_m,lagt=lagt_m, variogram=variogram)
            if(bivariate && gaussian)
              result <- list(lags=lags_m,lagt=lagt_m, variogram11=variogram11,
                             variogram12=variogram12,variogram22=variogram22)
        } else {
            if(!bivariate) result$variogram <- variogram
            if(bivariate && gaussian) {
              result$variogram11 <- variogram11
              result$variogram12 <- variogram12
              result$variogram22 <- variogram22
            }
        }
    }
    if(answer.range) {
      if(!is.list(result)) result <- list(range = Range)
      else result$range <- Range
    }
   ########################################################## 
if(!is.null(result)) {
        invisible(result)
    } else if(!is.null(OLS)) {
        invisible(OLS)
    } else {
        invisible()
    }
  }
