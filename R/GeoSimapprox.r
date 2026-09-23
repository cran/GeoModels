####################################################
### File name: GeoSimapprox.r
####################################################

# Simulate approximate spatial and spatio-temporal random felds:
GeoSimapprox <- function(coordx=NULL, coordy=NULL,coordz=NULL,coordt=NULL, coordx_dyn=NULL,corrmodel, distance="Eucl", grid=FALSE,
 max.ext=1,method="TB", L=10000,model='Gaussian', parallel=FALSE,ncores=6,
     n=1, param, anisopars=NULL,radius=6371,X=NULL,spobj=NULL,nrep=1,progress=TRUE,check.duplicates=FALSE)
{
    call <- match.call()
####################################################################
############  starting internal function ###########################
####################################################################
ddim<-function(coordx,coordy,coordz,coordt,grid=FALSE)
{
  dimt <- if(is.null(coordt)) 1L else length(coordt)
  if(isTRUE(grid)) {
    dims <- length(coordx) * max(1L, length(coordy)) * max(1L, length(coordz))
  } else if(!is.null(dim(coordx))) {
    dims <- nrow(as.matrix(coordx))
  } else {
    dims <- length(coordx)
  }
  dims * dimt
}
###############################
RFfct1<- function(numcoord,numtime,spacetime,bivariate,dime,nuisance,simd,X,ns)
    {
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
          X <- .GeoBivariateDesign(X, ns, context = "GeoSimapprox")
          mean_vec <- .GeoBivariateObservedMean(
            nuisance, X, ns, context = "GeoSimapprox"
          )
          sim <- mean_vec + as.numeric(simd)
        }

        if(!spacetime&&!bivariate) sim <- c(sim)
        else if(bivariate) {
          ## The approximate bivariate kernel uses a common set of sites.
          sim <- matrix(sim, nrow=2L, ncol=ns[1L], byrow=TRUE)
        } else sim <- matrix(sim, nrow=numtime, ncol=numcoord,byrow=TRUE)
        return(sim)
    }
#########################################################

######################################################################
############## main  internal function  ##############################
######################################################################
simu_approx=function(numxgrid,numygrid,coordx,coordy,coordz,coords,coordt,method,corrmodel,param,M,L,bivariate,spacetime,parallel,ncores)
{
##### spatial case 


if(!spacetime){    
   ## Turning Bands
 if(method=="TB") { 
        if(bivariate)
        {
          tau1 <- as.numeric(param$nugget_1)
          tau2 <- as.numeric(param$nugget_2)
          sill1 <- as.numeric(param$sill_1)
          sill2 <- as.numeric(param$sill_2)
          if(any(lengths(list(tau1,tau2,sill1,sill2)) != 1L) ||
             any(!is.finite(c(tau1,tau2,sill1,sill2))) ||
             tau1 < 0 || tau1 >= 1 || tau2 < 0 || tau2 >= 1 ||
             sill1 <= 0 || sill2 <= 0)
            stop("Bivariate TB simulation requires positive sills and nugget parameters in [0,1).", call.=FALSE)
          raw <- tbm2d(coords,coordt, param, corrmodel,L,bivariate)
          if(is.null(raw)) return(NULL)
          raw <- as.matrix(raw)
          raw[,1] <- sqrt(sill1) * (sqrt(1-tau1) * raw[,1] + sqrt(tau1) * rnorm(nrow(raw)))
          raw[,2] <- sqrt(sill2) * (sqrt(1-tau2) * raw[,2] + sqrt(tau2) * rnorm(nrow(raw)))
          simu <- raw
        }
        else
        {   tau2=param$nugget
            vvv=param$sill
            ##----------###
            param$nugget=0
            param$sill=1
            ##----------###

         simu=tbm2d_uni(coords, coordt, param, corrmodel, L, bivariate,parallel,ncores)$simu


         if(is.null(simu)) {print("TB simulation method fails\n");return(simu)}
         simu=sqrt(vvv)*(sqrt(1-tau2)*simu+sqrt(tau2)*rnorm(length(simu)))
        }
           simu=c(simu)
    }     


 ## Circulant embeeding           for regular grid          
  if(method=="CE") simu= c(SimCE(numxgrid,numygrid,coordx,coordy,coordz,corrmodel,param,mean.val=0, max.ext)$X)        
}
##### spacetime case ######
if(spacetime)  {
 # if(method=="Vecchia") print("ciao") 
  if(method=="CE") simu=CE_Space_Time(
      coords=coords, time.seq=coordt, param=param,
      corrmodel=corrmodel, distance=distance, max.ext=max.ext
    )
    }
return(simu)
}
####################################################################
####################################################################
############# END internal functions ###############################
####################################################################
####################################################################
    M=NULL
    corrmodel <- .GeoSimulationCanonicalCorrModel(corrmodel)
    if(is.null(CkCorrModel(corrmodel))) stop("The name of the correlation model is not correct\n")
    if(!(method=="TB"||method=="CE")) stop("The method of simulation is not correct\n")
    model_requested <- gsub("[[:blank:]]", "",model)
    model <- .GeoSimulationCanonicalModel(model_requested)
    n <- .GeoSimulationEffectiveN(model_requested, model, n)
    nrep <- .GeoValidateNrep(nrep)
    .GeoValidateDirectSimulationModel(model, simulator = "GeoSimapprox")
    distance=gsub("[[:blank:]]", "",distance)
    method=gsub("[[:blank:]]", "",method)
    if(method == "TB") L <- .GeoValidateTBControls(L)
    numxgrid=numygrid=NULL
    spacetime_dyn <- !is.null(coordx_dyn)
    coordx_dyn_output <- coordx_dyn
    grid_output <- grid

param <- lapply(param, unname) 
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

### some check for parallel #######
n.cores <- .GeoResolveWorkers(parallel, ncores)
parallel <- isTRUE(parallel) && n.cores > 1L
##################################
## Flag consumed by simu_approx()/tbm2d_uni() calls below: whether the
## TB simulation should parallelize internally over spatial chunks.
## Its actual value is decided later (right before running the
## replicates), based on whether outer (nrep) or inner (points)
## parallelism is the better fit -- only one of the two is ever active.
parallel_inner <- FALSE
##################################

##############################################################################
###### extracting sp object informations if necessary              ###########
##############################################################################
bivariate<-CheckBiv(CkCorrModel(corrmodel))
spacetime<-CheckST(CkCorrModel(corrmodel))
space=!spacetime&&!bivariate
if(bivariate && model != "Gaussian") {
  stop("GeoSimapprox currently supports bivariate approximate simulation only for model = 'Gaussian'; use GeoSim() for supported bivariate simulation.",
       call. = FALSE)
}
if(method == "TB") .GeoValidateTBCorrelationModel(corrmodel, bivariate = bivariate)
biv_info <- NULL
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
if(method == "TB") .GeoValidateTBDistance(distance, context = "GeoSimapprox")
if(bivariate) {
  biv_info <- .GeoBivariateCoordinateInfo(
    coordx=coordx, coordy=coordy, coordz=coordz,
    coordx_dyn=coordx_dyn, grid=grid, context="GeoSimapprox"
  )
  if(biv_info$dynamic) {
    same_support <- identical(dim(biv_info$blocks[[1L]]), dim(biv_info$blocks[[2L]])) &&
      isTRUE(all.equal(biv_info$blocks[[1L]], biv_info$blocks[[2L]],
                       tolerance=sqrt(.Machine$double.eps), check.attributes=FALSE))
    if(!same_support)
      stop("GeoSimapprox: approximate bivariate simulation requires common locations for the two variables; use GeoSim for heterotopic coordinates.", call.=FALSE)
  }
  ## tbm2d simulates both components on one common coordinate matrix.
  coordx <- biv_info$blocks[[1L]]
  coordy <- coordz <- NULL
  coordx_dyn <- NULL
  spacetime_dyn <- FALSE
  grid <- FALSE
}
if(method == "CE" && bivariate) {
  stop("GeoSimapprox does not support bivariate simulation with method = 'CE'; use method = 'TB' or GeoSim().",
       call. = FALSE)
}
if(spacetime && method == "TB") {
  stop("GeoSimapprox does not implement turning-bands simulation for spatio-temporal models; use method = 'CE' for supported fixed-location separable models or use GeoSim().",
       call. = FALSE)
}
if(method == "CE" && spacetime_dyn) {
  stop("CE spatio-temporal simulation requires fixed spatial locations; use GeoSim() for dynamic locations.",
       call. = FALSE)
}
.GeoCheckDuplicateLocations(
  coordx = coordx, coordy = coordy, coordz = coordz, coordt = coordt,
  coordx_dyn = coordx_dyn, spacetime = spacetime, bivariate = bivariate,
  grid = grid, context = "GeoSimapprox", check.duplicates = check.duplicates
)
###############################################################
###############################################################

    if(spacetime_dyn) {
      coords <- do.call(rbind, lapply(coordx_dyn, as.matrix))
    } else if(grid) {
      xgrid=coordx;ygrid=coordy
      numxgrid=length(xgrid);numygrid=length(ygrid)
      if(is.null(coordz)) {
        coords=as.matrix(expand.grid(xgrid,ygrid))
      } else {
        zgrid=coordz;numzgrid=length(zgrid)
        coords=as.matrix(expand.grid(xgrid,ygrid,zgrid))
      }
    } else {
      coords=coordx
      if(!is.null(coordy)) coords=cbind(coordx,coordy,coordz)
      coords <- as.matrix(coords)
    }
    coords_orig=coords


    .GeoValidate3D(
      coordx = if (is.null(coordx_dyn)) coords else NULL,
      coordx_dyn = coordx_dyn, distance = distance, grid = grid,
      approximate = TRUE, method = method, context = "GeoSimapprox"
    )

####
    spatial_dimension <- ncol(coords)
    if(method == "TB" && identical(spatial_dimension, 3L)) {
      if(!isTRUE(space) || isTRUE(bivariate) || isTRUE(spacetime)) {
        stop("GeoSimapprox supports three-dimensional TB simulation only for purely spatial univariate fields.",
             call. = FALSE)
      }
      if(!is.null(anisopars)) {
        stop("GeoSimapprox does not support anisopars for three-dimensional TB simulation; transform the coordinates explicitly.",
             call. = FALSE)
      }
    }

    numtime <- if(is.null(coordt)) 1L else length(coordt)
    if(spacetime_dyn) numtime=1L
    numcoord=nrow(coords)
    dime <- if(spacetime_dyn) nrow(coords) else numcoord*numtime
    if(bivariate) numtime=2
    if(!bivariate && model %in% c("Binomial", "BinomialLogistic")) {
      n <- .GeoValidateIntegerVector(n, "n", nobs = dime, minimum = 1L)
    }
    if(!bivariate && model %in% c("BinomialNeg", "BinomialNegZINB")) {
      n <- .GeoValidateIntegerVector(n, "n", nobs = NULL, minimum = 1L)
      if(length(n) != 1L)
        stop("For direct Negative-Binomial simulation, n is the common number r of successes and must be a positive integer scalar.", call.=FALSE)
    }
  ####

    if(method == "CE" && space) {
      if(!isTRUE(grid))
        stop("The spatial CE method works only for a regular two-dimensional grid.", call.=FALSE)
      if(tolower(distance) != "eucl")
        stop("The spatial CE method requires distance = 'Eucl'.", call.=FALSE)
      if(!is.null(anisopars))
        stop("The spatial CE method does not support anisopars; transform the grid explicitly or use method = 'TB'.", call.=FALSE)
      if(is.null(coordy) || !is.null(coordz))
        stop("The spatial CE method requires separate two-dimensional grid axes coordx and coordy.", call.=FALSE)
    }

    if(!(is.null(anisopars))) coords=GeoAniso(coords,c(anisopars$angle,anisopars$ratio))
    ##############################################################################
    ##############################################################################

    if(spacetime)
       {if(!(corrmodel %in% c("Matern_Matern","GenWend_GenWend","GenWend_Matern_GenWend_Matern")))
          stop("GeoSimapprox(method='CE') does not implement this spatio-temporal correlation model.", call.=FALSE)
        }
    if(bivariate) {
         if(!(corrmodel %in% c("Bi_matern","Bi_Matern"))) stop("Not implemented for this correlation model  \n")}
   ################################################################################
    unname(coordt);
    if(is.null(coordx_dyn)){
    unname(coordx);unname(coordy);;unname(coordz)}
  
  ################################################################################
  ################ setting parameters for each model #############################
  ################################################################################
     if(!bivariate) {
       nobs_mean <- dime
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
       if(!all(names(param_check) %in% c(CorrParam(corrmodel),
                                         NuisParam(model,bivariate,num_betas=num_betas))))
         stop("only nuisance and correlation parameters must be included in param\n")
     } else {
       ns_biv <- biv_info$ns
       X <- .GeoBivariateDesign(X, ns_biv, context="GeoSimapprox")
       p_mean_biv <- ncol(X)
       invisible(.GeoBivariateValidateMeanParameters(
         param, p_mean_biv, context="GeoSimapprox"
       ))
       mean_biv <- .GeoBivariateObservedMean(
         param, X, ns_biv, context="GeoSimapprox"
       )
       num_betas <- c(p_mean_biv, p_mean_biv)
       allowed <- c(CorrParam(corrmodel),
                    NuisParam(model, bivariate=TRUE, num_betas=num_betas))
       if(!all(names(param) %in% allowed))
         stop("only nuisance and correlation parameters must be included in param\n")
     }
    
    if(!bivariate) { if(is.null(param$sill))  param$sill=1 }
if(!bivariate) { if(is.null(param$sill)) param$sill=1}
else   { if(is.null(param$sill_1)) param$sill_1=1
         if(is.null(param$sill_2)) param$sill_2=1}


#################################
    if(model %in% c("SkewGaussian","SkewGauss","Beta",'Kumaraswamy','Kumaraswamy2','LogGaussian',#"Binomial","BinomialNeg","BinomialNegZINB",
                    "StudentT","SkewStudentT","Poisson","TwoPieceTukeyh","PoissonZIP","PoissonGamma","PoissonGammaZIP","PoissonWeibull",
                     "TwoPieceBimodal", "TwoPieceStudentT","TwoPieceGaussian","TwoPieceGauss","Tukeyh","Tukeyh2","Tukeygh","SinhAsinh",
                    "Gamma","Weibull","LogLogistic","Logistic","SkewLaplace","BinomialLogistic"))
       {
        if(spacetime_dyn && is.list(X)) X <- do.call(rbind, lapply(X, as.matrix))

  if(!bivariate){
           mm <- .GeoMean_vector(param, X = X, n = nobs_mean, external = external_mean)
           param <- .GeoMean_zero(param, num_betas)
        

        if((model %in% c("SkewGaussian","SkewGauss","TwoPieceGaussian","Logistic","SkewLaplace",
          "TwoPieceGauss","Gamma","Weibull","LogLogistic","Poisson","PoissonZIP","Tukeyh","Tukeyh2","PoissonGamma",
          "PoissonGammaZIP","PoissonWeibull",
          'LogGaussian',"TwoPieceTukeyh","TwoPieceBimodal", "Tukeygh","SinhAsinh",
                    "StudentT","SkewStudentT","TwoPieceStudentT")))  {vv<-param$sill; param$sill=1}

        if(model%in% c("SkewGaussian","SkewGauss","SkewStudentT","TwoPieceTukeyh","TwoPieceBimodal",
               "TwoPieceStudentT","TwoPieceGaussian","TwoPieceGauss"))
               { sk<-param$skew
                 if(model%in% c("TwoPieceTukeyh")) tl<-param$tail
                 if(model%in% c("TwoPieceBimodal")) bimo<-param$shape
               }
        }
        else {
            mm <- mean_biv
            param <- .GeoMean_zero(param, p_mean_biv, "mean_1")
            param <- .GeoMean_zero(param, p_mean_biv, "mean_2")
            vv1<-param$sill_1;param$sill_1=1-param$nugget_1;
            vv2<-param$sill_2;param$sill_2=1-param$nugget_2;vv=c(vv1,vv2)
            sk1<-param$skew_1;sk2<-param$skew_2;sk=c(sk1,sk2)
        }
      }
#################################
  if(model %in% c("Tukeygh","SinhAsinh"))  {
         if(!bivariate){
          param$mean=0
          sk<-param$skew; tl<-param$tail}
         else {
            mm <- mean_biv
            param <- .GeoMean_zero(param, p_mean_biv, "mean_1")
            param <- .GeoMean_zero(param, p_mean_biv, "mean_2")
            vv1<-param$sill_1;param$sill_1=1;vv2<-param$sill_2;param$sill_2=1;vv=c(vv1,vv2)
            sk1<-param$skew_1;sk2<-param$skew_2;sk=c(sk1,sk2)
            tl1<-param$tail_1;tl2<-param$tail_2;tl=c(tl1,tl2)
        }}

    if(model %in% c("Tukeyh"))  {
         if(!bivariate){
          param$mean=0
          tl<-param$tail}
         else {
            mm <- mean_biv
            param <- .GeoMean_zero(param, p_mean_biv, "mean_1")
            param <- .GeoMean_zero(param, p_mean_biv, "mean_2")
            vv1<-param$sill_1;param$sill_1=1;vv2<-param$sill_2;param$sill_2=1;vv=c(vv1,vv2)
            tl1<-param$tail_1;tl2<-param$tail_2;sk=c(tl1,tl2)
        }}

    if(model %in% c("Tukeyh2"))  {
         if(!bivariate){
          param$mean=0
          t1l<-param$tail1
          t2l<-param$tail2
           }
     }
#################################
   if(model %in% c("Wrapped"))  {
        if(!bivariate){
            mm <- 2 * atan(.GeoMean_vector(param, X = X, n = nobs_mean,
                                           external = external_mean)) + pi
            param <- .GeoMean_zero(param, num_betas)}
        else {
            mm <- 2*atan(mean_biv)+pi
            param <- .GeoMean_zero(param, p_mean_biv, "mean_1")
            param <- .GeoMean_zero(param, p_mean_biv, "mean_2")
        }}

###################################################
############# start  do_one_sim ###################
###################################################

do_one_sim <- function(LL, pb = NULL) {
  if (!is.null(pb)) pb(sprintf("LL=%g", LL))
    spec <- .GeoSimulationLatentSpec(model, param, n)
    k <- spec$k
    npoi <- spec$npoi
    param <- spec$param
     #if(model %in% c("Beta")) {  k=round(param$shape1)+round(param$shape2)
  ################################################################################
  ################################################################################
   ns <- if(bivariate) biv_info$ns else NULL
   if(spacetime_dyn) {
        coords=NULL
       if(bivariate) coordt=c(1,2)
       coords=do.call(rbind,args=c(coordx_dyn))
       ns=lengths(coordx_dyn)/2
       if(ncol(coords)==2) {coordx <- coords[,1]; coordy <- coords[,2];coordz <- NULL}
         if(ncol(coords)==3) {coordx <- coords[,1]; coordy <- coords[,2];coordz <- coords[,3]}
       dime=sum(ns)
   }
   else { dime=ddim(coordx,coordy,coordz,coordt,grid=grid)
          if(bivariate) ns <- biv_info$ns
        }
   if(bivariate || model %in% c("BinomialNeg", "BinomialNegZINB")) {
     dd <- NULL
   } else {
     dd <- array(0, dim = c(dime, 1L, k))
   }
 #########################################

#########################################################
KK=1;ssp=double(dime)
if(model%in% c("SkewGaussian","StudentT","SkewStudentT","TwoPieceTukeyh",
               "TwoPieceStudentT","TwoPieceGaussian"))
{
      simD=simu_approx(numxgrid,numygrid,coordx,coordy,coordz,coords,coordt,method,corrmodel,param,M,L,bivariate,spacetime,parallel_inner,n.cores)
      if(!spacetime&&!bivariate) simDD <- c(simD)
      else simDD <- matrix(simD, nrow=numtime, ncol=numcoord,byrow=TRUE)
      param$nugget=0 #ojo
}

if(model %in% c("BinomialNeg","BinomialNegZINB")) {
  rvec <- rep.int(as.integer(n[1L]), dime)
  counts <- integer(dime)
  finished <- logical(dime)
  nb_vals <- integer(dime)
  trial <- 0L
}
if(model %in% c("Poisson", "PoissonZIP", "PoissonGamma", "PoissonGammaZIP", "PoissonWeibull")) {
  poi_counts <- integer(dime)
}

while(KK<=npoi) {
for(i in 1:k) {

################# here the approximated simulation  ##################################################
simd=simu_approx(numxgrid,numygrid,coordx,coordy,coordz,coords,coordt,method,corrmodel,param,M,L,bivariate,spacetime,parallel_inner,n.cores)
######################################################################################################

namesnuis<-NuisParam(model, bivariate,num_betas=num_betas,copula=NULL)
nuisance<-param[namesnuis]
if(i==1&&(model=="SkewGaussian"||model=="SkewGauss")&&bivariate) param["pcol"]=0
####################################
sim<- RFfct1(numcoord,numtime,spacetime,bivariate,
                  dime,nuisance,simd,X,ns)
     
####################################
####### starting cases #############
####################################
if(model %in% c("Binomial", "BinomialNeg","BinomialNegZINB")) {

        simdim <- dim(sim)
        sim <- as.numeric(sim>0)
        dim(sim) <- simdim
         }
####################################
if(model %in% c("Weibull","SkewGaussian","SkewGauss","Binomial","BinomialLogistic","Poisson","PoissonGamma","PoissonWeibull","PoissonZIP","PoissonGammaZIP","Beta","Kumaraswamy","Kumaraswamy2",
              "LogGaussian","TwoPieceTukeyh",
                "Gamma","LogLogistic","Logistic","SkewLaplace","StudentT",
                "SkewStudentT","TwoPieceStudentT","TwoPieceGaussian","TwoPieceGauss","TwoPieceBimodal")) {
       if(!bivariate) dd[,,i]=t(sim)
       if(bivariate)  dd[,,i]=t(sim)
     }
 ####################################
if(model %in% c("BinomialNeg","BinomialNegZINB")){
                 trial <- trial + 1L
                 vec <- c(t(sim))
                 idx <- !finished
                 if(any(idx)) {
                   counts[idx] <- counts[idx] + vec[idx]
                   new_finished <- (counts >= rvec) & !finished
                   if(any(new_finished)) {
                     nb_vals[new_finished] <- trial - rvec[new_finished]
                     finished <- finished | new_finished
                   }
                 }
                 if(all(finished)) break
               }
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
 if(model %in% c("Binomial","BinomialLogistic","Poisson","PoissonGamma","PoissonWeibull","PoissonGammaZIP","PoissonZIP","BinomialNeg","BinomialNegZINB"))   {
   if(model %in% c("Poisson","PoissonGamma","PoissonWeibull"))   {sim=poi_counts;byrow=TRUE}
    if(model %in% c("PoissonZIP"))   {
      a=simu_approx(numxgrid,numygrid,coordx,coordy,coordz,coords,coordt,method,corrmodel,param,M,L,bivariate,spacetime,parallel_inner,n.cores)
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
          if(any(!finished)) stop("Some locations did not reach n successes during BinomialNeg simulation.", call.=FALSE)
          sim <- nb_vals
          byrow=TRUE
          }
  if(model %in% c("BinomialNegZINB"))   {
          if(any(!finished)) stop("Some locations did not reach n successes during BinomialNegZINB simulation.", call.=FALSE)
          sim <- nb_vals
      a=simu_approx(numxgrid,numygrid,coordx,coordy,coordz,coords,coordt,method,corrmodel,param,M,L,bivariate,spacetime,parallel_inner,n.cores)
     ###
          a[a<as.numeric(param$pmu)]=0;a[a!=0]=1
          sim=a*sim
          byrow=TRUE
          }
#############################################
############### formatting data #############
#############################################

    sim <- .GeoSimulationFormat(
      sim, grid = grid, space = !spacetime && !bivariate,
      numxgrid = numxgrid, numygrid = numygrid, numtime = numtime,
      numcoord = numcoord, byrow = byrow
    )
}
#########################################################################################################
#### Common transformation of latent Gaussian fields for univariate models.
#########################################################################################################
if(!bivariate) {
  transformed <- .GeoSimulationTransformUnivariate(
    model = model, sim = sim, dd = dd, simDD = get0("simDD", ifnotfound = NULL, inherits = TRUE),
    mm = get0("mm", ifnotfound = 0, inherits = TRUE),
    vv = get0("vv", ifnotfound = 1, inherits = TRUE),
    sk = get0("sk", ifnotfound = 0, inherits = TRUE),
    tl = get0("tl", ifnotfound = 0, inherits = TRUE),
    t1l = get0("t1l", ifnotfound = 0, inherits = TRUE),
    t2l = get0("t2l", ifnotfound = 0, inherits = TRUE),
    bimo = get0("bimo", ifnotfound = NULL, inherits = TRUE),
    param = param, k = k, spacetime = spacetime
  )
  if(transformed$handled) {
    sim <- transformed$sim
    if(transformed$needs_format) {
      sim <- .GeoSimulationFormat(
        sim, grid = grid, space = !spacetime, numxgrid = numxgrid, numygrid = numygrid,
        numtime = numtime, numcoord = numcoord, byrow = transformed$byrow
      )
    }
  }
} else {
  ## GeoSimapprox rejects non-Gaussian bivariate marginals above.
  sim <- .GeoSimulationFormat(
    c(sim), grid = grid, space = FALSE, numxgrid = numxgrid, numygrid = numygrid,
    numtime = numtime, numcoord = numcoord, byrow = FALSE
  )
}
##################################################################
###########. formatting data for space time dynamic case. #########
    if(spacetime_dyn) sim <- .GeoSimulationSplitDynamic(sim, ns)
##################################################################
##################################################################
##################################################################

    sim 
  }

##################################################################################
##################### end do_one_sim #############################################
##################################################################################

## Decide ONCE which level of parallelism to use -- outer (across the
## nrep replicates, via future_lapply) or inner (within a single TB
## simulation, across spatial chunks, inside tbm2d_uni). The two are
## never active at the same time to avoid nested-parallelism overhead
## and oversubscription.
##
## For TB, the scheduler uses both N * L (work in one field) and
## N * L * nrep (aggregate work), plus the number of independent
## replicates.  This replaces the old rules based on N * nrep and the
## overly conservative N * L >= 2.5e9 inner threshold.
force_inner <- isTRUE(getOption("GeoModels.force_tb_inner_parallel", FALSE))

if (identical(method, "TB")) {
  tb_parallel_mode <- .GeoTBParallelMode(
    dime = dime,
    L = L,
    nrep = nrep,
    parallel = parallel,
    ncores = n.cores,
    force_inner = force_inner
  )
  use_outer <- identical(tb_parallel_mode, "outer")
  use_inner <- identical(tb_parallel_mode, "inner")
} else {
  ## Preserve the pre-existing scheduling rule for non-TB approximate
  ## simulators (e.g. CE).  The new calibration is specific to TB.
  cost_outer <- as.numeric(dime) * as.numeric(nrep)
  use_outer <- isTRUE(parallel) && nrep > 1 && cost_outer > 1000000
  use_inner <- FALSE
}

parallel_inner <- use_inner   # read by the simu_approx() calls inside do_one_sim

if (use_outer) {

         cat("Performing", nrep, "simulation using",method ," method"," ...\n")
  # Preserve the caller's warning policy; future relays worker warnings.

  old_plan_outer <- future::plan()
  on.exit(try(future::plan(old_plan_outer), silent = TRUE), add = TRUE)
  future::plan(future::multisession, workers = n.cores)

  run_outer <- function(pb = NULL) {
    future.apply::future_lapply(
      seq_len(nrep),
      function(LL) {
        if (!is.null(pb)) pb(sprintf("Simulation %d of %d", LL, nrep))
        do_one_sim(LL)
      },
      future.seed = TRUE
    )
  }

  ## Do not open an empty with_progress() context when progress=FALSE.
  ## This is important when GeoSimapprox() is called inside another
  ## with_progress() block (for example by GeoVarest): an unnecessary
  ## nested context invalidates the outer progressor.
  if (isTRUE(progress)) {
    SIM <- progressr::with_progress({
      pb <- progressr::progressor(steps = nrep)
      run_outer(pb)
    })
  } else {
    SIM <- run_outer(NULL)
  }

} else {

  if (use_inner) {
    ## The plan is opened ONCE here, before looping over replicates, and
    ## reused by every call to tbm2d_uni() (which no longer opens/closes
    ## its own plan). This avoids paying worker start/stop overhead on
    ## every single replicate.
    cat("Performing", nrep, "simulation using",method ,
        "method (parallel over spatial chunks) ...\n")

    ## GeoVarest opens one multisession plan for the complete simulation
    ## phase and reuses it across all K one-replicate calls.  Stand-alone
    ## GeoSimapprox calls still manage and restore their own plan.
    external_plan <- isTRUE(getOption("GeoModels.tb_inner_plan_active", FALSE))
    if (!external_plan) {
      old_plan_inner <- future::plan()
      future::plan(future::multisession, workers = n.cores)
      on.exit(try(future::plan(old_plan_inner), silent = TRUE), add = TRUE)
    }
  }

  ## GeoVarest can provide one deterministic seed per requested replicate.
  ## This internal hook allows several TB replicates to share the expensive
  ## GeoSimapprox setup while preserving exactly the same RNG stream that
  ## would be obtained from separate nrep=1 calls.  It is deliberately an
  ## option rather than a public argument: stand-alone GeoSimapprox semantics
  ## are unchanged.
  rep_seeds_internal <- getOption("GeoModels.GeoSimapprox_rep_seeds", NULL)
  if (!is.null(rep_seeds_internal)) {
    rep_seeds_internal <- suppressWarnings(as.numeric(rep_seeds_internal))
    if (length(rep_seeds_internal) != nrep ||
        any(!is.finite(rep_seeds_internal)) ||
        any(rep_seeds_internal < 0) ||
        any(rep_seeds_internal > .Machine$integer.max)) {
      stop("Internal GeoSimapprox replicate seeds are invalid.", call. = FALSE)
    }
    rep_seeds_internal <- as.integer(rep_seeds_internal)
  }

  with_rep_seed <- function(seed_value, expr) {
    if (is.null(seed_value)) return(force(expr))

    ## Preserve and restore the caller RNG state using base R.
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else NULL

    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)

    set.seed(seed_value)
    force(expr)
  }

  run_sequential_replicates <- function(pb = NULL) {
    results <- vector("list", nrep)
    for (LL in seq_len(nrep)) {
      if (!is.null(pb)) pb(sprintf("Simulation %d of %d", LL, nrep))
      if (is.null(rep_seeds_internal)) {
        results[[LL]] <- do_one_sim(LL)
      } else {
        results[[LL]] <- with_rep_seed(
          rep_seeds_internal[[LL]],
          do_one_sim(LL)
        )
      }
    }
    results
  }

  ## For nrep=1, and whenever progress=FALSE, evaluate directly.  Avoiding
  ## an empty nested with_progress() context keeps an outer GeoVarest
  ## progressor alive throughout its complete simulation phase.
  if (isTRUE(progress) && nrep > 1L) {
    SIM <- progressr::with_progress({
      pb <- progressr::progressor(steps = nrep)
      run_sequential_replicates(pb)
    })
  } else {
    SIM <- run_sequential_replicates(NULL)
  }
}



if (nrep == 1) SIM <- SIM[[1]]

################################################
################################################
################################################



    #######################################
    if(bivariate)   numtime=1
    # Delete the global variables:
    # Return the objects list:
    GeoSim <- list(bivariate = bivariate,
    coordx = coords_orig[,1],
    coordy = coords_orig[,2],
    coordz = if(ncol(coords_orig) >= 3L) coords_orig[,3] else NULL,
    coordt = coordt,
    coordx_dyn =coordx_dyn_output,
    corrmodel = corrmodel,
    data = SIM,
    distance = distance,
    grid = grid_output,
    model = model_requested,
    method=method,
    n=n,
    numcoord = numcoord,
    numtime = numtime,
    param = param_input,
    radius = radius,
    nrep=nrep,
    spacetime = spacetime,
    X=X)
#}
##############################################
    structure(c(GeoSim, call = call), class = c("GeoSim"))
}