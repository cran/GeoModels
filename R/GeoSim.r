####################################################
### File name: GeoSim.r
####################################################


# Simulate spatial and spatio-temporal random felds:
GeoSim <- function(coordx=NULL, coordy=NULL,coordz=NULL, coordt=NULL, coordx_dyn=NULL,corrmodel, distance="Eucl", grid=FALSE,
     method="cholesky",model='Gaussian', n=1, param, anisopars=NULL, radius=1,
      sparse=FALSE,X=NULL,spobj=NULL,nrep=1,progress=TRUE,check.duplicates=FALSE)
{
    call <- match.call()
####################################################################
############ internal function #####################################
####################################################################
##############################################################################
     RFfct1<- function(ccov,dime,nuisance,simd,X,ns)
    {
        numcoord=ccov$numcoord; numtime=ccov$numtime;grid=ccov$grid;
        spacetime=ccov$spacetime;bivariate=ccov$bivariate

        if(!bivariate) {
          external_mean <- .GeoMean_external(nuisance, n = dime)
          if(!is.null(external_mean)) {
            sim <- external_mean + as.numeric(simd)
          } else {
            X <- .GeoMean_design(X, n = dime, name = "X")
            beta <- .GeoMean_beta(nuisance, ncol(X))
            sim <- as.numeric(X %*% beta) + as.numeric(simd)
          }
        }
        if(bivariate) {
          X <- .GeoBivariateDesign(X, ns, context = "GeoSim")
          mu_biv <- .GeoBivariateObservedMean(
            nuisance, X, ns, context = "GeoSim"
          )
          sim <- mu_biv + as.numeric(simd)
        }
            if(!spacetime&&!bivariate) sim <- c(sim)
            else sim <- matrix(sim, nrow=numtime, ncol=numcoord,byrow=TRUE)
        return(sim)
    }
####################################################################
############# END internal functions ###############################
####################################################################


    if(!is.character(corrmodel) || length(corrmodel) != 1L) stop("the name of the correlation model is wrong")
    corrmodel <- .GeoSimulationCanonicalCorrModel(corrmodel)
    if(is.null(CkCorrModel(corrmodel))) stop("the name of the correlation model is wrong")
    model_requested <- gsub("[[:blank:]]", "",model)
    model <- .GeoSimulationCanonicalModel(model_requested)
    n <- .GeoSimulationEffectiveN(model_requested, model, n)
    nrep <- .GeoValidateNrep(nrep)
    .GeoValidateDirectSimulationModel(model, simulator = "GeoSim")

    distance=gsub("[[:blank:]]", "",distance)
    method <- .GeoValidateExactSimulationMethod(method)

if(is.null(coordz))
{
    if(grid) { xgrid=coordx;ygrid=coordy;
               numxgrid=length(xgrid);numygrid=length(ygrid) }
}
else {if(grid) stop("grid can not be coupled with coordz \n ")}


spacetime_dyn=FALSE


param <- lapply(param, unname) ## carlo
param_input <- param

# Enforce the finite-moment Tukey-hh parameter domain at the public
# simulation entry point.  Keep this explicit even though the common
# marginal validator performs the same check: GeoSim/GeoSimapprox must
# never proceed with h_l or h_r >= 1/2.
if (identical(model, "Tukeyh2")) {
    for (nm in c("tail1", "tail2")) {
        value <- suppressWarnings(as.numeric(param[[nm]]))
        if (length(value) != 1L || !is.finite(value) || value < 0 || value >= 0.5)
            stop(sprintf("param$%s must lie in [0, 0.5) for Tukeyh2 simulation.", nm),
                 call. = FALSE)
    }
}

param <- .GeoValidateLatentIntegerParameters(model, param)
param <- .GeoValidateSimulationMarginal(model, param)
##############################################################################
###### extracting sp object informations if necessary              ###########
##############################################################################
bivariate<-CheckBiv(CkCorrModel(corrmodel))
spacetime<-CheckST(CkCorrModel(corrmodel))
space=!spacetime&&!bivariate
if(bivariate && !(model %in% c("Gaussian", "Gauss"))) {
  stop("GeoSim currently supports bivariate simulation only for model = 'Gaussian'.",
       call. = FALSE)
}
if(!is.null(spobj)) {
   if(space||bivariate){
        a=sp2Geo(spobj,NULL); coordx=a$coords 
       if(!a$pj) {if(distance!="Chor") distance="Geod"}
    }
   if(spacetime){
        a=sp2Geo(spobj,NULL); coordx=a$coords ; coordt=a$coordt 
        if(!a$pj) {if(distance!="Chor") distance="Geod"}
     }
}
###############################################################
###############################################################


    .GeoValidate3D(
      coordx = coordx, coordy = coordy, coordz = coordz,
      coordx_dyn = coordx_dyn, distance = distance, grid = grid,
      context = "GeoSim"
    )

    if(spacetime && !is.null(coordx_dyn)) {
      coordx_dyn <- .GeoNormalizeDynamicCoordinates(
        coordx_dyn, coordt = coordt, context = "GeoSim"
      )
      .GeoValidateDynamicRows(X, coordx_dyn, name = "X", context = "GeoSim")
    }
    if(bivariate) {
      biv_info <- .GeoBivariateCoordinateInfo(
        coordx = coordx, coordy = coordy, coordz = coordz,
        coordx_dyn = coordx_dyn, grid = grid, context = "GeoSim"
      )
      if(biv_info$dynamic) coordx_dyn <- biv_info$blocks
      if(!is.null(X)) {
        X <- .GeoBivariateDesign(X, biv_info$ns, context = "GeoSim")
      }
    }

    .GeoCheckDuplicateLocations(
      coordx = coordx, coordy = coordy, coordz = coordz, coordt = coordt,
      coordx_dyn = coordx_dyn, spacetime = spacetime, bivariate = bivariate,
      grid = grid, context = "GeoSim", check.duplicates = check.duplicates
    )
    .dup_old <- .GeoSuppressDuplicateChecksForNestedCalls()
    on.exit(.GeoRestoreDuplicateCheckOption(.dup_old), add = TRUE)

    if(!is.null(coordx_dyn))  spacetime_dyn=TRUE
    unname(coordt);
    if(is.null(coordx_dyn)){
    unname(coordx);unname(coordy)}

if(!bivariate) { if(is.null(param$sill))  param$sill=1 }


################################################################################
################ starting spatial and spatiotemporal case#######################
################################################################################
 if(!bivariate)
    {
 ############### check and construct the mean specification
 ntime_mean <- if(is.null(coordt)) 1L else length(coordt)
 if(!is.null(coordx_dyn)) {
   nobs_mean <- sum(vapply(coordx_dyn, function(z) nrow(as.matrix(z)), integer(1L)))
 } else if(isTRUE(grid)) {
   nz_mean <- if(is.null(coordz)) 1L else length(coordz)
   nobs_mean <- length(coordx) * length(coordy) * nz_mean * ntime_mean
 } else {
   nspace_mean <- if(is.matrix(coordx) || is.data.frame(coordx)) nrow(as.matrix(coordx)) else length(coordx)
   nobs_mean <- nspace_mean * ntime_mean
 }
 if(model %in% c("Binomial", "BinomialLogistic")) {
   n <- .GeoValidateIntegerVector(n, "n", nobs = nobs_mean, minimum = 1L)
 }
 if(model %in% c("Geom", "Geometric")) {
   n <- .GeoValidateIntegerVector(n, "n", nobs = NULL, minimum = 1L)
   if(length(n) != 1L || as.integer(n[1L]) != 1L)
     stop("For direct Geometric simulation, n is fixed to 1 (the r = 1 Negative-Binomial case).", call.=FALSE)
   n <- 1L
 }
 if(model %in% c("BinomialNeg", "BinomialNegZINB")) {
   n <- .GeoValidateIntegerVector(n, "n", nobs = NULL, minimum = 1L)
   if(length(n) != 1L)
     stop("For direct Negative-Binomial simulation, n is the common number r of successes and must be a positive integer scalar.", call.=FALSE)
 }
 external_mean <- .GeoMean_external(param, n = nobs_mean)
 if(!is.null(external_mean) && !is.null(X))
   stop("A site-specific mean vector and X cannot be used together.\n")
 if(!is.null(X)) {
   X <- .GeoMean_design(X, n = nobs_mean, name = "X")
   num_betas <- ncol(X)
   invisible(.GeoMean_beta(param, num_betas))
 } else {
   supplied_mean <- .GeoMean_supplied_names(param)
   if(is.null(external_mean) && length(setdiff(supplied_mean, "mean")))
     stop("Multiple mean coefficients require an explicit design matrix X.\n")
   num_betas <- 1L
   if(is.null(external_mean)) invisible(.GeoMean_beta(param, 1L))
 }
 param_check <- param
 if(!is.null(external_mean)) param_check$mean <- 0
 if( !all(names(param_check) %in% c(CorrParam(corrmodel), NuisParam2(model,bivariate,num_betas=num_betas))) )
   stop("Nuisance and correlation parameters must be included in param\n")
#################################
    if(model %in% c("SkewGaussian","Beta",'Kumaraswamy','Kumaraswamy2','LogGaussian',
        "Binomial","BinomialNeg",#"BinomialNegZINB",
                    "StudentT","SkewStudentT","Poisson","TwoPieceTukeyh","PoissonZIP","PoissonGamma","PoissonGammaZIP","PoissonWeibull",
                     "TwoPieceBimodal", "TwoPieceStudentT","TwoPieceGaussian","TwoPieceGauss","Tukeyh","Tukeyh2","Tukeygh","SinhAsinh",
                    "Gamma","Weibull","LogLogistic","Logistic","SkewLaplace","BinomialLogistic"))
       {
           if(spacetime_dyn && is.list(X)) X <- do.call(rbind, lapply(X, as.matrix))
           mm <- .GeoMean_vector(param, X = X, n = nobs_mean, external = external_mean)
           param <- .GeoMean_zero(param, num_betas)

        if((model %in% c("SkewGaussian","TwoPieceGaussian","Logistic","SkewLaplace",
          "TwoPieceGauss","Gamma","Weibull","LogLogistic","Poisson","PoissonZIP","Tukeyh","Tukeyh2","PoissonGamma","PoissonGammaZIP","PoissonWeibull",
          'LogGaussian',"TwoPieceTukeyh","TwoPieceBimodal", "Tukeygh","SinhAsinh",
          "Binomial","BinomialNeg",
                    "StudentT","SkewStudentT","TwoPieceStudentT","Gaussian")))  
        {
          vv<-param$sill;
          param$sill=1
        }


        if(model%in% c("SkewGaussian","SkewStudentT","TwoPieceTukeyh","TwoPieceBimodal",
               "TwoPieceStudentT","TwoPieceGaussian","TwoPieceGauss"))
               { sk<-param$skew
               if(model%in% c("TwoPieceTukeyh")) tl<-param$tail
               if(model%in% c("TwoPieceBimodal")) bimo<-param$shape
               }
       }
#################################
  if(model %in% c("Tukeygh","SinhAsinh"))  {
          param$mean=0
          sk<-param$skew; tl<-param$tail}
    if(model %in% c("Tukeyh"))  {
          param$mean=0
          tl<-param$tail}
         if(model %in% c("Tukeyh2"))  {
          param$mean=0
          t1l<-param$tail1
          t2l<-param$tail2
          }
#################################
   if(model %in% c("Wrapped"))  {
            mm <- 2 * atan(.GeoMean_vector(param, X = X, n = nobs_mean,
                                           external = external_mean)) + pi
            param <- .GeoMean_zero(param, num_betas)
        }

if(model%in% c("SkewGaussian","StudentT","SkewStudentT","TwoPieceTukeyh",
               "TwoPieceStudentT","TwoPieceGaussian"))
     { nugget=param$nugget;param$nugget=0}  ### ojo!!


zz=param[grep("mean", names(param))]   
ll=list(sill=param$sill,nugget=param$nugget)
mc=append(zz,ll)

pc=param[CorrParam(corrmodel)]



ccov = GeoCovmatrix(coordx=coordx, coordy=coordy,coordz=coordz, coordt=coordt, coordx_dyn=coordx_dyn, corrmodel=corrmodel,
                   distance=distance,grid=grid,model="Gaussian", n=n,
                param=append(mc,pc), anisopars=anisopars, radius=radius, sparse=sparse,copula=NULL,
                X=X)

######################
# matrix decomposition and square root
if(!sparse) {
   decompvarcov <- MatDecomp(ccov$covmatrix,method)
   if(is.logical(decompvarcov)){print(" Covariance matrix is not positive definite");stop()}
   sqrtvarcov <- MatSqrt(decompvarcov,method)
}
else {
 cholS <- spam::chol(ccov$covmatrix)
 # cholS is the upper triangular part of the permutated matrix Sigma
 iord <- spam::ordering(cholS, inv=TRUE)
 R <- spam::as.spam(cholS)
}

######################
##############################################
## putting 2 nugget
 if(model%in% c("PoissonZIP","BinomialNegZINB","PoissonGammaZIP"))
{
  II=diag(nrow(ccov$covmatrix));
  II[lower.tri(II)] <- (1-param$nugget2)/(1-param$nugget1)
  II <- t(II);
  II[lower.tri(II)] <- (1-param$nugget2)/(1-param$nugget1)
  ccov_with_nug=(ccov$covmatrix)*(II)

if(!sparse) {
  decompvarcov1 <- MatDecomp(ccov_with_nug,method)
  if(is.logical(decompvarcov1)){print(" Covariance matrix is not positive definite");stop()}
  sqrtvarcov1 <- MatSqrt(decompvarcov1,method)
  }
  else {
   cholS1 <- spam::chol(ccov_with_nug)
   iord1 <- spam::ordering(cholS1, inv=TRUE)
   R1 <- spam::as.spam(cholS1)
  }
}
##############################################
## putting 1 nugget
if(model%in% c("SkewGaussian","StudentT","SkewStudentT","TwoPieceTukeyh",
               "TwoPieceStudentT","TwoPieceGaussian"))
{
  II=diag(nrow(ccov$covmatrix));
  II[lower.tri(II)] <- (1-nugget);
  II <- t(II);
  II[lower.tri(II)] <- (1-nugget)
  ccov_with_nug=ccov$covmatrix*II

if(!sparse) {
  decompvarcov1 <- MatDecomp(ccov_with_nug,method)
if(is.logical(decompvarcov1)){print(" Covariance matrix is not positive definite");stop()}
  sqrtvarcov1 <- MatSqrt(decompvarcov1,method)
  }
  else {   
    cholS1 <- spam::chol(ccov_with_nug)
    iord1 <- spam::ordering(cholS1, inv=TRUE)
    R1 <- spam::as.spam(cholS1)
   }
}
##########################################################################################################################################
 numtime=1
 numcoord=ccov$numcoord;
 if(spacetime) numtime=ccov$numtime;
 dime<-numcoord*numtime

########################

   if(spacetime_dyn) {
        coords=NULL
       coords=do.call(rbind,args=c(coordx_dyn))
      ns <- vapply(coordx_dyn, nrow, integer(1L))
      if(ncol(coords)==2) {coordx <- coords[,1]; coordy <- coords[,2];coordz=NULL}
      if(ncol(coords)==3) {coordx <- coords[,1]; coordy <- coords[,2];coordz=coords[,3]}
       dime=sum(ns)
      ccov$numtime=1
   }
   else { dime=nrow(ccov$covmatrix) }
########################
    
SIM=list()
###################################################
### starting number of replicates
###################################################
if(nrep>1){

cat("Performing", nrep, "simulations using Cholesky method...\n")}
if(progress){
restore_progress_handlers <- .GeoProgressHandlersPush(
  handler = "txtprogressbar", global = TRUE
)
on.exit(restore_progress_handlers(), add = TRUE)
pb <- progressr::progressor(
  along = seq_len(nrep), on_exit = FALSE, auto_finish = FALSE
)
## Finish the progressor before restoring the progressr handlers.  This is
## essential because progressr may temporarily capture stdout while a global
## handler is active; leaving the progressor unfinished can leave that sink
## open after GeoSim() returns.
on.exit(pb(type = "finish"), add = TRUE, after = FALSE)

}

for( L in seq_len(nrep)){
    spec <- .GeoSimulationLatentSpec(model, param, n)
    k <- spec$k
    npoi <- spec$npoi
    param <- spec$param
    

  ################################################################################
  ################################################################################
   if(model %in% c("BinomialNeg","BinomialNegZINB")) {
      dd <- NULL
   } else {
      dd <- array(0, dim = c(dime, 1, k))
   }

 #########################################

#########################################################
KK=1;ssp=double(dime)

if(model%in% c("SkewGaussian","StudentT","SkewStudentT","TwoPieceTukeyh",
               "TwoPieceStudentT","TwoPieceGaussian"))
{
  ss=matrix(rnorm(dime) , nrow=dime, ncol = 1)
  if(sparse) simD=as.numeric((array(c(ss),c(1,dime)) %*% R1) [,iord1] )
    else     #simD=crossprod(sqrtvarcov1,ss) 
             simD=sqrtvarcov1%*%ss

      ccov1=ccov
      ccov1$covmatrix=ccov_with_nug
      if(space) simDD <- c(simD)
      else simDD <- matrix(simD, nrow=ccov1$numtime, ncol=ccov1$numcoord,byrow=TRUE)
}

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # INIT efficient counters for BinomialNeg / BinomialNegZINB (univariate)
  if(model %in% c("BinomialNeg","BinomialNegZINB")){
    rvec <- rep.int(as.integer(n[1L]), dime)
    counts   <- integer(dime)     # successi cumulati
    finished <- logical(dime)     # chi ha raggiunto n successi
    nb_vals  <- integer(dime)     # output (tempo al successo r) = trial - r
    trial    <- 0L                # contatore di tentativi (iterazioni Bernoulli)
  }
  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  if(model %in% c("Poisson", "PoissonZIP", "PoissonGamma", "PoissonGammaZIP", "PoissonWeibull")) {
    poi_counts <- integer(dime)
  }

  while(KK<=npoi) {
  for(i in 1:k) {

    ss=matrix(rnorm(dime) , nrow=dime, ncol = 1)

    #### simulating N(0,1) with matrix decomposition using sparse or dense matrices without nugget!
      if(sparse)   simd=as.numeric((array(c(ss),c(1,dime)) %*% R) [,iord] )
      else         #simd=crossprod(sqrtvarcov,ss) 
                   simd=sqrtvarcov%*%ss



    #######################################################################
    nuisance<-param[ccov$namesnuis]
    sim<-RFfct1(ccov,dime,nuisance,simd,ccov$X,ns)

    ####################################
    ####### starting cases #############
    ####################################
     #if(model %in% c("Binomial", "BinomialNeg","BinomialNegZINB")) {
      #   simdim <- dim(sim)
       #  sim <- as.numeric(sim>0)
        # dim(sim) <- simdim
         # }
 if(model %in% c("Binomial", "BinomialNeg","BinomialNegZINB")) {
    simdim <- dim(sim)
    sim_vec <- if(is.null(simdim)) c(sim) else c(t(sim))
    eta <- if(length(mm) == 1) rep(mm, length(sim_vec)) else c(mm)
    sim_vec <- as.numeric(sim_vec + eta > 0)
    if(is.null(simdim)) sim <- sim_vec
    else sim <- matrix(sim_vec, nrow=simdim[1], ncol=simdim[2], byrow=TRUE)
    }    
    ####################################
    if(model %in% c("Weibull","SkewGaussian","SkewGauss","Binomial","BinomialLogistic","Poisson","PoissonGamma","PoissonGammaZIP","PoissonWeibull","PoissonZIP","Beta","Kumaraswamy","Kumaraswamy2",
              "LogGaussian","TwoPieceTukeyh",
                "Gamma","LogLogistic","Logistic","SkewLaplace","StudentT",
                "SkewStudentT","TwoPieceStudentT","TwoPieceGaussian","TwoPieceGauss","TwoPieceBimodal")) {
        dd[,,i]=t(sim)

     }
     ####################################
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if(model %in% c("BinomialNeg","BinomialNegZINB")){
                 # vec Bernoulli corrente (ordine identico a c(t(sim)) usato prima)
                 trial <- trial + 1L
                 vec <- c(t(sim))
                 idx <- !finished
                 if(any(idx)){
                   counts[idx] <- counts[idx] + vec[idx]
                   new_finished <- (counts >= rvec) & !finished
                   if(any(new_finished)){
                     nb_vals[new_finished] <- trial - rvec[new_finished]
                     finished <- finished | new_finished
                   }
                 }
                 if(all(finished)) {break;}### stopping rule (tutti chiusi)
               }
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    }
 ####################################
  if(model %in% c("Poisson","PoissonZIP"))   {
   pois1=0.5*(dd[,,1]^2+dd[,,2]^2)
   ssp=ssp+c(pois1)
   active <- ssp <= c(exp(mm))
   poi_counts <- poi_counts + as.integer(active)
   if(!any(active)) break  ## stopping rule
 }
  ####################################
 if(model %in% c("PoissonGamma","PoissonGammaZIP"))   {
   if(KK==1) sim3 <- .GeoSimulationSquareSums(dd, 3:k)
   #################################
   pois1=0.5*(dd[,,1]^2+dd[,,2]^2)
   ssp=ssp+c(pois1)
   threshold <- c(exp(mm)*sim3/(k-2))
   active <- ssp <= threshold
   poi_counts <- poi_counts + as.integer(active)
   if(!any(active)) break  ## stopping rule
 }
if(model %in% c("PoissonWeibull"))   {

  if(KK==1) sim3 <- .GeoSimulationSquareSums(dd, 3:k)
   #################################
   pois1=0.5*(dd[,,1]^2+dd[,,2]^2)
   ssp=ssp+c(pois1)
   threshold <- c(exp(mm)*(sim3/2)^(1/param$shape)/(gamma(1+1/param$shape)))
   active <- ssp <= threshold
   poi_counts <- poi_counts + as.integer(active)
   if(!any(active)) break  ## stopping rule
 }

 KK=KK+1
}

 ####### end for #########################
 ###############################################################################################
 #### simulation for discrete random field based on indipendent copies  of GRF ######
 ###############################################################################################
 if(model %in% c("Binomial","BinomialLogistic","Poisson","PoissonGamma","PoissonGammaZIP","PoissonWeibull","PoissonZIP","BinomialNeg","BinomialNegZINB"))   {

   if(model %in% c("Poisson","PoissonGamma","PoissonWeibull"))   {sim=poi_counts;byrow=TRUE}
    
    if(model %in% c("PoissonZIP","PoissonGammaZIP"))   {
      ss=matrix(rnorm(dime) , nrow=dime, ncol = 1)

      if(sparse)   a=as.numeric((array(c(ss),c(1,dime)) %*% R1) [,iord1] )
      else         #a=crossprod(sqrtvarcov1,ss) 
                       a=sqrtvarcov1%*%ss

     ###
      a[a<as.numeric(param$pmu)]=0;a[a!=0]=1
      sim=a*poi_counts;
      byrow=TRUE
      }
########################################
   if(model == "Binomial") {
                  sim <- .GeoSimulationBinomialCounts(dd, n)
                  byrow <- TRUE }
########################################
if(model == "BinomialLogistic") {
                  sim <- .GeoSimulationBinomialLogisticCounts(dd, n, mm)
                  byrow <- TRUE }
#######################################
   if(model %in% c("BinomialNeg"))   {
          if(any(!finished)) {
              stop("Some locations did not reach n successes. Increase the maximum number of trials.")
          }
          sim <- nb_vals
          byrow=TRUE
          }
  if(model %in% c("BinomialNegZINB"))   {
          if(any(!finished)) {
              stop("Some locations did not reach n successes. Increase the maximum number of trials.")
          }
          sim <- nb_vals
      ss=matrix(rnorm(dime) , nrow=dime, ncol = 1)

      if(sparse)   a=as.numeric((array(c(ss),c(1,dime)) %*% R1) [,iord1] )
      else         #a=crossprod(sqrtvarcov1,ss) 
                   a=sqrtvarcov1%*%ss
     ###
          a[a<as.numeric(param$pmu)]=0;a[a!=0]=1
          sim=a*sim
          byrow=TRUE
          }
   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#############################################
############### formatting data #############
#############################################

    sim <- .GeoSimulationFormat(
      sim, grid = grid, space = space, numxgrid = numxgrid, numygrid = numygrid,
      numtime = numtime, numcoord = numcoord, byrow = byrow
    )
}
#########################################################################################################
#### Common transformation of latent Gaussian fields for univariate models.
#########################################################################################################
    transformed <- .GeoSimulationTransformUnivariate(
      model = model, sim = sim, dd = dd, simDD = get0("simDD", ifnotfound = NULL, inherits = FALSE),
      mm = get0("mm", ifnotfound = 0, inherits = FALSE),
      vv = get0("vv", ifnotfound = 1, inherits = FALSE),
      sk = get0("sk", ifnotfound = 0, inherits = FALSE),
      tl = get0("tl", ifnotfound = 0, inherits = FALSE),
      t1l = get0("t1l", ifnotfound = 0, inherits = FALSE),
      t2l = get0("t2l", ifnotfound = 0, inherits = FALSE),
      bimo = get0("bimo", ifnotfound = NULL, inherits = FALSE),
      param = param, k = k, spacetime = spacetime
    )
    if(transformed$handled) {
      sim <- transformed$sim
      if(transformed$needs_format) {
        sim <- .GeoSimulationFormat(
          sim, grid = grid, space = space, numxgrid = numxgrid, numygrid = numygrid,
          numtime = numtime, numcoord = numcoord, byrow = transformed$byrow
        )
      }
    }
##################################################################
###########. formatting data for space time dynamic case. #########
    if(spacetime_dyn) sim <- .GeoSimulationSplitDynamic(sim, ns)
###################################################################

SIM[[L]]=sim
if(progress) pb()
} 
######################################################################
##########  end of replicates ########################################
######################################################################
}
###################################################################
################ end spatial and spatiotemporal case###############
###################################################################


################################################################
################ starting bivariate  case#######################
################################################################
if(bivariate){
  sel1 <- substr(names(param), 1, 6) == "mean_1"
  sel2 <- substr(names(param), 1, 6) == "mean_2"
  num_betas <- c(sum(sel1), sum(sel2))

  zz <- param[grep("mean", names(param))]
  pc <- param[CorrParam(corrmodel)]
  ccov <- GeoCovmatrix(
    coordx = coordx, coordy = coordy, coordz = coordz, coordt = coordt,
    coordx_dyn = coordx_dyn, corrmodel = corrmodel, distance = distance,
    grid = grid, model = "Gaussian", n = n, param = append(pc, zz),
    anisopars = anisopars, radius = radius, sparse = sparse, copula = NULL, X = X
  )

  if(!sparse) {
    decompvarcov <- MatDecomp(ccov$covmatrix, method)
    if(is.logical(decompvarcov)) {
      print(" Covariance matrix is not positive definite")
      stop()
    }
    sqrtvarcov <- MatSqrt(decompvarcov, method)
  } else {
    cholS <- spam::chol(ccov$covmatrix)
    iord <- spam::ordering(cholS, inv = TRUE)
    R <- spam::as.spam(cholS)
  }

  ns <- as.integer(ccov$ns)
  if(length(ns) != 2L || any(ns < 1L))
    stop("GeoSim: invalid bivariate block sizes returned by GeoCovmatrix.")
  dime <- nrow(ccov$covmatrix)
  if(dime != sum(ns))
    stop("GeoSim: covariance dimension does not match the bivariate block sizes.")

  SIM <- vector("list", nrep)
  for(L in seq_len(nrep)) {
    if(spacetime_dyn) {
      ccov$numtime <- 1L
      ccov$numcoord <- dime
    }
    numcoord <- ccov$numcoord
    numtime <- ccov$numtime
    ss <- matrix(rnorm(dime), nrow = dime, ncol = 1L)
    if(sparse) simd <- as.numeric((array(c(ss), c(1L, dime)) %*% R)[, iord])
    else simd <- sqrtvarcov %*% ss

    nuisance <- ccov$param[ccov$namesnuis]
    sim <- RFfct1(ccov, dime, nuisance, simd, ccov$X, ns)
    sim <- .GeoSimulationFormat(
      c(sim), grid = grid, space = FALSE, numxgrid = numxgrid, numygrid = numygrid,
      numtime = numtime, numcoord = numcoord, byrow = FALSE
    )
    if(spacetime_dyn) sim <- .GeoSimulationSplitDynamic(sim, ns)
    SIM[[L]] <- sim
  }
}

######################################################################
##########  end of bivariate  ########################################
######################################################################

if(nrep==1) SIM=SIM[[1]] 

##################################################################
#######################################
    if(ccov$bivariate)   ccov$numtime=1

    # Delete the global variables:
    # Return the objects list:
    GeoSim <- list(bivariate = bivariate,
    coordx = ccov$coordx,
    coordy = ccov$coordy,
    coordz = ccov$coordz,
    coordt = ccov$coordt,
    coordx_dyn =coordx_dyn,
    corrmodel = corrmodel,
    data = SIM,
    distance = distance,
    grid = grid,
    model = model_requested,
    method=method,
    n=n,
    numcoord = ccov$numcoord,
    numtime = ccov$numtime,
    param = param_input,
    radius = radius,
    #randseed=.Random.seed,
    spacetime = spacetime,
    sparse=ccov$sparse,
    nrep=nrep,
    X=X)
#}
##############################################
    structure(c(GeoSim, call = call), class = c("GeoSim"))
}
