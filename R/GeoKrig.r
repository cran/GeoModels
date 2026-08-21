
####################################################
### File name: GeoKrig.r
####################################################

GeoKrig= function(estobj=NULL,data, coordx=NULL, coordy=NULL,coordz=NULL, coordt=NULL, coordx_dyn=NULL, corrmodel, distance="Eucl", grid=FALSE, loc, 
                method="cholesky", model="Gaussian", n=1,nloc=NULL, mse=FALSE,  param, anisopars=NULL,
               radius=1, sparse=FALSE,time=NULL,  type_krig="Simple",weigthed=TRUE,
               which=1, copula=NULL,X=NULL,Xloc=NULL,Mloc=NULL,spobj=NULL,spdata=NULL,varcov=NULL,progress=FALSE,
               check.duplicates=FALSE)

{
##################################################################################################################
##################################################################################################################

############## internal functions are defined in KrigUtils.R ############################


######################################
########## START #####################
######################################
call <- match.call()

###############################
## checking if there is a  GeoFit object
###############################
if(!is.null(estobj)){
   if(!inherits(estobj,"GeoFit"))
               stop("need  a 'GeoFit' object as input\n")

data=estobj$data
################################################################
if(!estobj$grid){  #not regular grid 
 if(!estobj$bivariate){  if(is.null(estobj$coordx_dyn)) coordx=cbind(estobj$coordx,estobj$coordy,estobj$coordz) ####spatial and spacetime fixed loc
                          else                          coordx=NULL
                      } ## spatial (temporal) non regular case
 else  {    if(is.null(estobj$coordx_dyn))  { coordx=estobj$coordx[1:estobj$ns[1]]    # bivariate not dynamic    
                                              coordy=estobj$coordy[1:estobj$ns[2]] 
                                            }  
            else {coordx=NULL}                                      # bivariate  dynamic  
       }
 }                  # regular grid
else  { coordx=estobj$coordx; coordy=estobj$coordy; coordz=estobj$coordz} # grid(nunca)

   if(isTRUE(estobj$spacetime)) coordt=estobj$coordt
   else coordt=NULL

   coordx_dyn=estobj$coordx_dyn
   corrmodel=estobj$corrmodel
   model=estobj$model
   distance=estobj$distance
   grid=estobj$grid
   n=estobj$n
   param=append(estobj$param,estobj$fixed)
   radius=estobj$radius
   copula=estobj$copula
   anisopars=estobj$anisopars
   varcov=estobj$varcov
   if(!is.null(estobj$sparse)) sparse=estobj$sparse
   if(!is.null(anisopars))  {param$angle=NULL; param$ratio=NULL}

   X=estobj$X
}
##################################
###### end check geofit object####
##################################
model_requested <- gsub("[[:blank:]]", "", as.character(model)[1L])
model <- .GeoPredictionCanonicalModel(model_requested)
param_input <- param
n_input <- n
nloc_input <- nloc
if(is.null(CkModel(model))) stop("The name of the  model  is not correct\n")
if(!is.character(corrmodel) || is.null(CkCorrModel(corrmodel))) stop("the name of the correlation model is wrong")
corrmodel=gsub("[[:blank:]]", "",corrmodel)
if (!is.null(copula)) {
  copula <- as.character(copula)[1L]
  if (is.na(copula) || !nzchar(gsub("[[:blank:]]", "", copula))) {
    copula <- NULL
  } else {
    copula <- gsub("[[:blank:]]", "", copula)
  }
}
if(!is.null(copula) && !(copula %in% c("Gaussian", "Clayton", "SkewGaussian")))
  stop("GeoKrig currently implements copula-based linear prediction for 'Gaussian', 'Clayton', and 'SkewGaussian'.", call.=FALSE)
distance=gsub("[[:blank:]]", "",distance)
method=.GeoKrigMatchMethod(method, sparse=sparse, context="GeoKrig")
type_key=tolower(gsub("[[:blank:]]", "",type_krig))
if(!(type_key %in% c("simple","universal"))) stop("type_krig must be 'Simple' or 'Universal'")
type_krig=if(type_key=="simple") "Simple" else "Universal"
if(!is.logical(weigthed) || length(weigthed)!=1L || is.na(weigthed))
  stop("weigthed must be TRUE or FALSE")
if(!is.logical(progress) || length(progress)!=1L || is.na(progress))
  stop("progress must be TRUE or FALSE")
if(!isTRUE(weigthed))
  warning("weigthed is retained for backward compatibility and has no effect in standard kriging.", call.=FALSE)


########################################################
####### extracting sp objects if necessary 
########################################################
bivariate<-CheckBiv(CkCorrModel(corrmodel))
spacetime<-CheckST(CkCorrModel(corrmodel))
space=!spacetime&&!bivariate
.GeoKrigValidatePredictionModel(model, copula = copula, bivariate = bivariate, context = "GeoKrig")
model_code_copula <- as.integer(CkModel(model))
copula_cont <- !bivariate && !is.null(copula) &&
  copula %in% c("Gaussian", "Clayton", "SkewGaussian") &&
  (model_code_copula %in% .GeoGaussianCopulaContinuousModels ||
   (copula %in% c("Clayton", "SkewGaussian") &&
    model_code_copula %in% .GeoGaussianCopulaDiscreteModels))
if(!is.null(spobj)) {
   if(space||bivariate){
        a=sp2Geo(spobj,spdata); coordx=a$coords 
       if(!a$pj) {if(distance!="Chor") distance="Geod"}
    }
   if(spacetime){
        a=sp2Geo(spobj,spdata); coordx=a$coords ; coordt=a$coordt 
        if(!a$pj) {if(distance!="Chor") distance="Geod"}
     }
   if (!is.null(a$Y)) data <- a$Y
   if (!is.null(a$X)) X <- a$X
}

if(spacetime && !is.null(coordx_dyn)) {
  coordx_dyn <- .GeoNormalizeDynamicCoordinates(
    coordx_dyn, coordt = coordt, context = "GeoKrig"
  )
  .GeoValidateDynamicData(data, coordx_dyn, context = "GeoKrig")
  .GeoValidateDynamicRows(X, coordx_dyn, name = "X", context = "GeoKrig")
}
if(bivariate) {
  biv_info <- .GeoBivariateCoordinateInfo(
    coordx = coordx, coordy = coordy, coordz = coordz,
    coordx_dyn = coordx_dyn, grid = grid, context = "GeoKrig"
  )
  if(biv_info$dynamic) coordx_dyn <- biv_info$blocks
  data <- .GeoBivariateDataVector(data, biv_info$ns, context = "GeoKrig")
  if(!is.null(X)) X <- .GeoBivariateDesign(X, biv_info$ns, context = "GeoKrig")
}
if(!bivariate && is.null(coordx_dyn)) {
  coords_input <- .GeoKrigCoordinateMatrix(
    coordx = coordx, coordy = coordy, coordz = coordz, grid = grid,
    context = "GeoKrig"
  )
  data <- .GeoKrigNormalizeData(
    data, nspace = nrow(coords_input),
    coordt = if(spacetime) coordt else NULL,
    grid = grid, context = "GeoKrig"
  )
}

.GeoCheckDuplicateLocations(
  coordx = coordx, coordy = coordy, coordz = coordz, coordt = coordt,
  coordx_dyn = coordx_dyn, spacetime = spacetime, bivariate = bivariate,
  grid = grid, context = "GeoKrig", check.duplicates = check.duplicates
)
.dup_old <- .GeoSuppressDuplicateChecksForNestedCalls()
on.exit(.GeoRestoreDuplicateCheckOption(.dup_old), add = TRUE)

######################################
############ some checks##############
######################################
if (is.vector(loc)) {
  loc <- matrix(loc, nrow = 1)
} else if (!is.matrix(loc)) {
  loc <- as.matrix(loc)
}

if(!(ncol(loc)==2||ncol(loc)==3))   stop("loc parameter must be a matrix N x 2 or N x 3")
if(!is.numeric(loc) || !nrow(loc) || any(!is.finite(loc)))
  stop("loc must be a finite numeric N x 2 or N x 3 matrix")
if(!is.null(anisopars) && !identical(distance, "Eucl"))
  stop("GeoKrig: anisopars requires distance = 'Eucl'.")


.GeoValidate3D(
  coordx = coordx, coordy = coordy, coordz = coordz,
  coordx_dyn = coordx_dyn, distance = distance, grid = grid,
  context = "GeoKrig"
)
.GeoValidate3D(
  coordx = loc, distance = distance, grid = FALSE,
  context = "GeoKrig prediction locations"
)
.GeoValidateCoordinateCompatibility(
  coordx = coordx, coordy = coordy, coordz = coordz,
  coordx_dyn = coordx_dyn, loc = loc, context = "GeoKrig"
)

if (!is.null(X) && !is.matrix(X) && !is.list(X)) X <- as.matrix(X)

## Mean specification.  A vector param$mean is an external mean at the
## observations; mean, mean1, ... are regression coefficients otherwise.
Mobs <- .GeoMean_external(param)
if (!is.null(Mobs) && !is.null(X)) {
  stop("A fixed observation mean and X cannot be used together.\n")
}
if (!is.null(Mloc) && !is.null(Xloc)) {
  stop("Use either Mloc or Xloc, not both.\n")
}
if (!is.null(Mobs) && is.null(Mloc)) {
  stop("Mloc is required when the observation mean is supplied as a vector.\n")
}
if (!is.null(X) && is.null(Xloc) && is.null(Mloc)) {
  stop("Xloc or Mloc is required when X is supplied.\n")
}
if (is.null(X) && !is.null(Xloc)) {
  stop("X must be supplied when Xloc is used.\n")
}
if (is.null(X) && is.null(Mobs) &&
    length(setdiff(.GeoMean_supplied_names(param), "mean"))) {
  stop("Parameters mean1, mean2, ... require an explicit design matrix X; without X the mean is intercept-only.\n")
}
if(identical(type_krig, "Universal") && (!is.null(Mobs) || !is.null(Mloc))) {
  warning("External observation/prediction means are treated as known; ",
          "type_krig is therefore set to 'Simple'.", call. = FALSE)
  type_krig <- "Simple"
}
if(identical(type_krig, "Universal") && isTRUE(mse) && is.null(varcov)) {
  stop("Universal kriging with mse = TRUE requires varcov (the covariance matrix of the estimated parameters).")
}
if(spacetime) {
  if(is.null(time)) stop("At least one temporal instant is needed for space-time kriging")
  if(!is.numeric(time) || !length(time) || any(!is.finite(time)))
    stop("time must be a non-empty finite numeric vector")
  time=as.numeric(time)
}
loc_orig=loc
if(!is.null(anisopars)) {  loc=GeoAniso(loc,c(anisopars$angle,anisopars$ratio))}

if(!is.null(Mloc)) {
    if(!is.numeric(Mloc)) stop("Mloc must be numeric")
    Mloc <- as.numeric(Mloc)
    if(any(!is.finite(Mloc))) stop("Mloc contains non-finite values")
}
######################################
############ end some checks##########
######################################



#### locations points to predict ###
if(is.null(time)) time=0; numloc = nrow(loc); tloc = length(time);
if(!tloc)  tloc = 1;

## Normalize observation/prediction size parameters before covariance construction.
## A location-specific n can never be averaged to invent nloc.
count_sizes <- NULL
count_model_names <- c("Binary", "Bernoulli", "Binomial",
                       "Geom", "Geometric", "BinomialNeg", "BinomialNegZINB",
                       "Poisson", "PoissonZIP", "PoissonGamma",
                       "PoissonGammaZIP", "PoissonGammaZIP1")
n_cov <- n
if (!bivariate && model %in% count_model_names) {
  count_sizes <- .GeoKrigCountSizes(
    model, n = n, nloc = nloc, nobs = length(data),
    npred = numloc * tloc, context = "GeoKrig"
  )
  ## GeoCovmatrix distinguishes a site-specific Binomial trial vector from
  ## common field-size parameters.  Keep the expanded observation/prediction
  ## vectors in count_sizes for Corr_c_bin, but pass only the representation
  ## accepted by the covariance constructor.
  if (identical(model, "Binomial")) {
    n_cov <- count_sizes$obs
  } else if (model %in% c("Binary", "Bernoulli", "Geom", "Geometric")) {
    n_cov <- 1L
  } else if (model %in% c("BinomialNeg", "BinomialNegZINB")) {
    ## .GeoKrigCountSizes() has already validated that n is one common
    ## positive integer r (and that nloc, if supplied, agrees with it).
    n_cov <- as.integer(round(n[1L]))
  }
}

if(ncol(loc) == 2L) {
  locx = loc[,1]; locy = loc[,2]; locz = double(nrow(loc))
} else {
  locx = loc[,1]; locy = loc[,2]; locz = loc[,3]
}



type="Standard"

######################################
############ standard kriging ########
######################################
if(type %in% c("Standard","standard")) {

################################################
######################################################## 
############ optimal prediction cases ##################   
######################################################## 
logGausstemp=SinhAsinhtemp=Tukeyh2temp=Tukeyhtemp=FALSE 

######################################################## 
######################################################## 
## Remove an external observation mean before covariance construction.
## The covariance code must never interpret its entries as mean1, mean2, ... .
if (!is.null(Mobs) && !copula_cont) param$mean <- 0
Mtemp <- Mobs
   if(model %in% c("Weibull","Gamma","LogLogistic","LogGaussian") && !copula_cont)          # non-copula multiplicative kriging uses a standardized covariance matrix
{
     paramtemp=param; sel=substr(names(param),1,4)=="mean";  ## selecting mean values
     meantemp=names(param[sel])             ## saving mean parameters
     param=param[!sel];                     ## not mean parameters
     param$mean=0;        ## mean must be =0 when calling covariance matrix
     Xtemp=X;X=NULL                         ## saving X and setting X=NULL
}
#####################################################
##### computing covariance matrix  ##################
#####################################################

covmatrix = GeoCovmatrix(coordx=coordx, coordy=coordy,coordz=coordz, coordt=coordt, coordx_dyn=coordx_dyn,
            corrmodel=corrmodel, distance= distance,grid=grid,model=model,n=n_cov,
            param=param, anisopars=anisopars, radius=radius,sparse=sparse,copula=copula,X=X)

#####################################################
## Normalize legacy/no-copula GeoCovmatrix objects.  Some fitted objects may
## store an absent copula as NULL, NA, or an empty string; downstream kriging
## logic must treat all of these representations identically.
if (!is.null(covmatrix$copula)) {
  cov_copula <- as.character(covmatrix$copula)[1L]
  if (is.na(cov_copula) || !nzchar(gsub("[[:blank:]]", "", cov_copula)))
    covmatrix$copula <- NULL
}
covmatrix$param=unlist(covmatrix$param)
if(bivariate) tloc=1
spacetime_dyn=FALSE; if(!is.null(covmatrix$coordx_dyn)) spacetime_dyn=TRUE
########dimat is dimension  ######
if(!spacetime_dyn) dimat=covmatrix$numcoord*covmatrix$numtime
if(spacetime_dyn)  dimat =sum(covmatrix$ns)
dimat2=numloc*tloc

## GeoCovmatrix needs the external observation locations for Gaussian-copula
## covariance construction, but downstream parameter extraction expects a
## scalar regression coefficient.  MM retains the external mean vector.
if (copula_cont && !is.null(Mobs)) param$mean <- 0

MM <- Mtemp
if(!is.null(MM)) param$mean <- 0
###############
if(model %in% c("Weibull","Gamma","LogLogistic","LogGaussian") && !copula_cont) {
          if(is.null(Xtemp)) X=matrix(1, nrow=dimat, ncol=1L)
          else               X=Xtemp
          param=paramtemp
          if(!is.null(Mtemp)) param$mean=0
          covmatrix$namesnuis=unique(c(meantemp,covmatrix$namesnuis))   
}
else { if(!is.null(X)) X=covmatrix$X 
           else        X=matrix(1,nrow=dimat,ncol=1)

     }

###############
if (!bivariate) {
  X <- .GeoMean_design(X, n = dimat, name = "X")
  num_betas <- ncol(X)
  if (!is.null(MM) && length(MM) != dimat)
    stop("The fixed observation mean must have one value per observation.\n")
  if (!is.null(Mloc)) {
    if (length(Mloc) == 1L) Mloc <- rep(Mloc, dimat2)
    if (length(Mloc) != dimat2)
      stop("Mloc must have one value per prediction point/time combination.\n")
  }
  if (is.null(Mloc)) {
    Xloc <- .GeoMean_design(Xloc, n = dimat2, p = num_betas, name = "Xloc")
  } else {
    Xloc <- NULL
  }
} else {
  X <- .GeoBivariateDesign(
    if(is.null(X)) covmatrix$X else X,
    covmatrix$ns, context = "GeoKrig"
  )
  num_betas <- ncol(X)
  invisible(.GeoBivariateValidateMeanParameters(
    param, num_betas, context = "GeoKrig"
  ))
  if (!is.null(Mloc)) {
    Mloc <- as.numeric(Mloc)
    if (length(Mloc) == 1L) Mloc <- rep(Mloc, dimat2)
    if (length(Mloc) != dimat2 || any(!is.finite(Mloc)))
      stop("GeoKrig: Mloc must contain one finite value per prediction location.")
    Xloc <- NULL
  } else {
    Xloc <- .GeoMean_design(Xloc, n = dimat2, p = num_betas, name = "Xloc")
  }
}
XX <- X
NS=0
if(spacetime||bivariate) { NS=cumsum(covmatrix$ns); NS=c(0,NS)[-(length(covmatrix$ns)+1)] }
###########
nuisance = if (bivariate) covmatrix$param[covmatrix$namesnuis] else param[covmatrix$namesnuis]; nuisance=Filter(Negate(is.null),nuisance)
## n is an argument rather than a named nuisance parameter for count models,
## but the generic copula moment/quantile engines need the common size value.
if (!isTRUE(bivariate) && !is.null(copula) &&
    copula %in% c("Clayton", "SkewGaussian") &&
    as.integer(CkModel(model)) %in% c(11L, 16L)) {
  nn <- if (as.integer(CkModel(model)) == 11L && !is.null(count_sizes))
    unique(as.numeric(c(count_sizes$obs, count_sizes$pred))) else as.numeric(n_cov)[1L]
  if (length(nn) != 1L || !is.finite(nn))
    stop("GeoKrig: copula-based Binomial prediction currently requires one common n across observations and prediction locations.", call. = FALSE)
  nuisance["n"] <- nn[1L]
}
if(!bivariate) {
  sel <- names(nuisance) %in% .GeoMean_names(num_betas)
  betas <- .GeoMean_beta(nuisance, num_betas)
}
if(bivariate) {
  sel <- substr(names(nuisance),1,4)=="mean"
  beta_biv <- .GeoBivariateValidateMeanParameters(
    nuisance, num_betas, context = "GeoKrig"
  )
  betas1 <- beta_biv$beta1
  betas2 <- beta_biv$beta2
}
## Count cross-covariance kernels use the same alphabetically ordered
## nuisance vector as GeoCovmatrix().  Keeping one ordering is essential for
## models such as PoissonGamma/PoissonGammaZIP where `shape` is not adjacent
## to `nugget` in NuisParam().
other_nuis_named <- nuisance[!sel]
other_nuis <- as.numeric(other_nuis_named[order(names(other_nuis_named))])

## Marginal location predictors used by copula covariance and
## prediction.  These are eta = X beta (or externally supplied values), not
## necessarily the marginal means on the observed scale.
eta_obs_cop <- eta_loc_cop <- NULL
if (copula_cont) {
  eta_obs_cop <- if (is.null(MM)) as.numeric(X %*% betas) else as.numeric(MM)
  eta_loc_cop <- if (!is.null(Mloc)) as.numeric(Mloc) else as.numeric(Xloc %*% betas)
  if (length(eta_obs_cop) != dimat || length(eta_loc_cop) != dimat2 ||
      any(!is.finite(c(eta_obs_cop, eta_loc_cop))))
    stop("GeoKrig: inconsistent copula marginal locations.", call. = FALSE)
}

############
cmodel=corrmodel
cdistance=distance
corrmodel = CkCorrModel(covmatrix$corrmodel)
distance = CheckDistance(covmatrix$distance)
corrparam = covmatrix$param[covmatrix$namescorr]# selecting the correlation parameters
if(bivariate) {if(!(which==1 || which==2)) {stop("which  parameter must be 1 or 2")}}
pred = NULL
varpred=vv=vv2=NULL
k = 0

###### observed coordinates in covariance ordering
ccc <- .GeoKrigObservedCoordinates(covmatrix, dimat, context = "GeoKrig")

###### spatial anisotropy case
## Apply anisotropy only after the final prediction-coordinate matrix
## has been constructed for spatial, grid, space-time, and dynamic inputs.
if(!is.null(anisopars)) {
    ccc=GeoAniso(ccc,c(anisopars$angle,anisopars$ratio))
}
###############################################################

if(bivariate) dataT <- .GeoBivariateDataVector(data, covmatrix$ns, context = "GeoKrig")
else if(spacetime && spacetime_dyn) dataT <- as.numeric(unlist(data, use.names = FALSE))
else dataT <- as.numeric(data)

####################################################################
############### blocked global kriging ###############################
####################################################################
## The high-volume observed-scale copula path and the ordinary Gaussian
## baseline are evaluated in prediction-location blocks.  The observed
## covariance matrix is factorized once and the factorization is reused for
## every block.  The full weights matrix is still retained because it is part
## of the historical GeoKrig return object; therefore blocking reduces peak
## cross-covariance memory without changing the returned object.
blocked_continuous <- !bivariate && (copula_cont ||
  (is.null(covmatrix$copula) && as.integer(covmatrix$model) == 1L))

if (blocked_continuous) {
  if (copula_cont) {
    moments_obs_cop <- .GeoCopulaMarginalMoments(
      covmatrix$model, nuisance, eta_obs_cop
    )
    moments_loc_cop <- .GeoCopulaMarginalMoments(
      covmatrix$model, nuisance, eta_loc_cop
    )
    mu <- as.numeric(moments_obs_cop$mean)
    muloc <- as.numeric(moments_loc_cop$mean)
    vvar_all <- as.numeric(moments_loc_cop$variance)
    eta_range_cop <- range(c(eta_obs_cop, eta_loc_cop), finite = TRUE)
  } else {
    mu <- if (is.null(MM)) as.numeric(X %*% betas) else as.numeric(MM)
    muloc <- if (!is.null(Mloc)) as.numeric(Mloc) else as.numeric(Xloc %*% betas)
    vvar_all <- rep(as.numeric(covmatrix$param["sill"]), dimat2)
    eta_range_cop <- NULL
  }

  if (length(mu) != dimat || length(muloc) != dimat2 ||
      length(vvar_all) != dimat2 ||
      any(!is.finite(c(mu, muloc, vvar_all)))) {
    stop("GeoKrig: inconsistent blocked marginal moments.", call. = FALSE)
  }

  if (identical(type_krig, "Universal")) {
    if (copula_cont) {
      d_obs <- .GeoCopulaMarginalMeanDerivative(
        covmatrix$model, nuisance, eta_obs_cop
      )
      d_loc <- .GeoCopulaMarginalMeanDerivative(
        covmatrix$model, nuisance, eta_loc_cop
      )
      X_system <- X * d_obs
      Xloc_system <- Xloc * d_loc
    } else {
      X_system <- X
      Xloc_system <- Xloc
    }
    mean_names_system <- .GeoMean_names(ncol(X_system))
  } else {
    X_system <- Xloc_system <- mean_names_system <- NULL
  }

  krig_factor <- .GeoKrigFactor(
    covmatrix, method = method, context = "GeoKrig"
  )
  krig_blocks <- .GeoKrigLocationBlocks(
    dimat, numloc, tloc, context = "GeoKrig"
  )
  n_krig_blocks <- length(krig_blocks)
  if (isTRUE(progress) && n_krig_blocks > 1L) {
    restore_progress_handlers <- .GeoProgressHandlersPush(
      handler = "txtprogressbar", global = TRUE
    )
    on.exit(restore_progress_handlers(), add = TRUE)
    krig_pb <- progressr::progressor(
      steps = n_krig_blocks, on_exit = FALSE, auto_finish = FALSE
    )
    ## Finish before restoring the temporary global handler, including on
    ## errors or user interrupts.
    on.exit(krig_pb(type = "finish"), add = TRUE, after = FALSE)
  } else {
    krig_pb <- function(...) {}
  }

  ## Historical API: GeoKrig returns all kriging weights.  Keep that result
  ## exactly, but never materialize the complete cross-covariance matrix.
  krig_weights <- matrix(NA_real_, nrow = dimat, ncol = dimat2)
  pp <- numeric(dimat2)
  vv <- if (isTRUE(mse)) numeric(dimat2) else NULL
  datas <- as.matrix(as.numeric(dataT) - as.numeric(mu))
  nugget_scale <- 1 - as.numeric(covmatrix$param["nugget"])
  if (!is.finite(nugget_scale))
    stop("GeoKrig: invalid nugget parameter.", call. = FALSE)

  loc_matrix <- cbind(locx, locy, locz)
  block_id <- 0L
  for (loc_idx in krig_blocks) {
    block_id <- block_id + 1L
    first_task <- (loc_idx[1L] - 1L) * tloc + 1L
    last_task <- loc_idx[length(loc_idx)] * tloc
    task_idx <- seq.int(first_task, last_task)
    ntarget_block <- length(task_idx)

    rho_block <- .GeoKrigCorrelationBlock(
      covmatrix = covmatrix, observed_coords = ccc,
      corrmodel = corrmodel, corrparam = corrparam, distance = distance,
      loc_block = loc_matrix[loc_idx, , drop = FALSE],
      time = time, tloc = tloc, NS = NS, which = which,
      context = "GeoKrig"
    )
    latent_block <- nugget_scale * rho_block

    if (copula_cont) {
      eta_obs_pair <- rep(eta_obs_cop, times = ntarget_block)
      eta_loc_pair <- rep(eta_loc_cop[task_idx], each = dimat)
      CC_block <- .GeoCopulaCrossCovFast(
        covmatrix$copula, latent_block, covmatrix$model, nuisance,
        eta1 = eta_obs_pair, eta2 = eta_loc_pair,
        eta_range = eta_range_cop
      )
      dim(CC_block) <- c(dimat, ntarget_block)
    } else {
      CC_block <- matrix(
        latent_block * as.numeric(covmatrix$param["sill"]),
        nrow = dimat, ncol = ntarget_block
      )
    }

    Xloc_block <- if (identical(type_krig, "Universal")) {
      Xloc_system[task_idx, , drop = FALSE]
    } else NULL

    krig_system_block <- .GeoKrigSystem(
      covmatrix = covmatrix, CC = CC_block, method = method,
      type_krig = type_krig, X = X_system, Xloc = Xloc_block,
      target_variance = vvar_all[task_idx], mse = mse, varcov = varcov,
      mean_names = mean_names_system, factor = krig_factor,
      context = "GeoKrig"
    )
    wb <- krig_system_block$weights
    krig_weights[, task_idx] <- wb
    pp[task_idx] <- muloc[task_idx] + as.numeric(crossprod(datas, wb))
    if (isTRUE(mse)) vv[task_idx] <- krig_system_block$mse

    rm(rho_block, latent_block, CC_block, wb, krig_system_block)
    if (copula_cont) rm(eta_obs_pair, eta_loc_pair)
    krig_pb(sprintf("block %d/%d", block_id, n_krig_blocks))
  }

  cvv <- c(vv)
  if (mse) {
    if (anyNA(cvv)) cvv[is.na(cvv)] <- 1e-30
    selp <- cvv <= 0
    if (any(selp)) cvv[selp] <- 1e-30
  }
  if (spacetime) {
    pred <- matrix(t(pp), nrow = tloc, ncol = numloc)
    varpred <- if (isTRUE(mse)) matrix(cvv, nrow = tloc, ncol = numloc) else NULL
  } else {
    pred <- c(pp)
    varpred <- if (isTRUE(mse)) cvv else NULL
  }
}

####################################################################
############### computing  correlation vector ######################
####################################################################
if(!blocked_continuous && covmatrix$model %in% c(1,10,59,18,21,12,26,24,27,38,29,39,28,9,34,40,20,22,25,33,42,50))    
{
## gaussian=1
## skew gaussian=10
## skew Laplace=59
## student =12
## skewstudent =18
## gamma=21
## weibull=26
## loglogistic=24
## loggaussian=22
## twopieceStudentT=27
## beta=28
## twopieceGaussian=29
## twopieceTukeyh=38
## twopiecebimodal=39
## tukeygh=9
## tukey=34
## tukeyhh=40
## sihasin=20

if((type=="Standard"||type=="standard")) {

cc<-dotCall64::.C64(
"Corr_c",
SIGNATURE=c(
"double","double","double","double","double","integer",
"integer","double","double","double","integer","integer",
"integer","integer","integer","integer","double","integer",
"integer","double","integer","integer","double"
),
 corri = dotCall64::vector_dc("double",dimat * dimat2),
ccc[,1],ccc[,2],ccc[,3],covmatrix$coordt,corrmodel, 0,
locx,locy,locz,covmatrix$numcoord,numloc,tloc,
covmatrix$ns,NS,covmatrix$numtime,corrparam,
covmatrix$spacetime,covmatrix$bivariate,time,distance,which-1,covmatrix$radius,
INTENT=c("w",rep("r",22)),NAOK=TRUE,PACKAGE="GeoModels"
)

####  transforming gaussian correlations depending on the (non)Gaussian  model
if(bivariate){  corri=cc$corri  }

else {
        rho=cc$corri; cc=(1-as.numeric(covmatrix$param["nugget"]))*rho

if (copula_cont) {
  eta_obs_pair <- rep(eta_obs_cop, times = dimat2)
  eta_loc_pair <- rep(eta_loc_cop, each = dimat)
  corri <- .GeoCopulaCrossCovFast(
    covmatrix$copula, cc, covmatrix$model, nuisance,
    eta1 = eta_obs_pair, eta2 = eta_loc_pair
  )
  moments_obs_cop <- .GeoCopulaMarginalMoments(
    covmatrix$model, nuisance, eta_obs_cop
  )
  moments_loc_cop <- .GeoCopulaMarginalMoments(
    covmatrix$model, nuisance, eta_loc_cop
  )
} else {

################ Gaussian ##################################
if(covmatrix$model==1)   {  vv=as.numeric(covmatrix$param['sill']); 
          if(is.null(covmatrix$copula)) corri=cc 
          else {  if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov(cc,covmatrix$model, nuisance)}
} 
#################### skew gaussian #########################################
if(covmatrix$model==10) {    
vv=as.numeric(covmatrix$param['sill'])
if(is.null(covmatrix$copula)){
    corr2=rho^2; sk=as.numeric(covmatrix$param['skew']);sk2=sk^2; 
    corri=(2*sk2)*(sqrt(1-corr2) + rho*asin(rho)-1)/(pi*vv+sk2*(pi-2)) + (cc*vv)/(vv+sk2*(1-2/pi))
                        }
    else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov(cc,covmatrix$model, nuisance)}
}
if(covmatrix$model==59) {    
  vv = as.numeric(covmatrix$param['sill'])
  if(is.null(covmatrix$copula)){
    sk = as.numeric(covmatrix$param['skew'])
    corri = cc^2
  }
  else {
    if(covmatrix$copula=="Gaussian")
      corri = gaussian_copula_cov(cc,covmatrix$model, nuisance)
  }
}
##################### student T ########################################
if(covmatrix$model==12) 
{
vv=as.numeric(covmatrix$param['sill']); nu=1/as.numeric(covmatrix$param['df'])
if(is.null(covmatrix$copula)){
   cc2 <- cc^2;gamma_num <- gamma((nu - 1)/2)^2;gamma_den <- gamma(nu/2)^2;Fvals <- hypergeo::hypergeo(0.5, 0.5, nu/2, cc2)
   corri <- ((nu - 2) * gamma_num * Re(Fvals) * cc) / (2 * gamma_den)
}
else {
      if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov(cc,covmatrix$model, nuisance)
     }
}
########################## tukeyh ###################################
if(covmatrix$model==34) 
{   vv=as.numeric(covmatrix$param['sill']); h=as.numeric(covmatrix$param['tail'])
    if(h>0){ 
            if(is.null(covmatrix$copula)){ corri=(cc*(1-2*h)^(1.5))/((1-h)^2-(h*cc)^2)^(1.5) }
            else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov(cc,covmatrix$model, nuisance)}
            }
    else{ if(is.null(covmatrix$copula))  corri=cc
          else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov(cc,covmatrix$model, nuisance)}
        }
}
######################### tukeyh2###################################
if(covmatrix$model==40)
{ vv=as.numeric(covmatrix$param['sill']);
  tail1=as.numeric(covmatrix$param['tail1']);tail2=as.numeric(covmatrix$param['tail2'])
  if(is.null(covmatrix$copula)){
    mom_hh <- .GeoTukeyh2StandardMoments(tail1, tail2)
    mm <- mom_hh$mean
    vv1 <- mom_hh$variance
    corri <- .GeoTukeyh2Correlation(cc, tail1, tail2)
  }
  else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov(cc,covmatrix$model, nuisance)}
}
####################### sas ######################################
if (covmatrix$model == 20) {
  vv <- as.numeric(covmatrix$param['sill'])
  tail <- as.numeric(covmatrix$param['tail'])
  skew <- as.numeric(covmatrix$param['skew'])
  if(is.null(covmatrix$copula)){
    mom_sas <- .GeoSASStandardMoments(skew, tail)
    mm <- mom_sas$mean
    vv1 <- mom_sas$variance
    corri <- corrsas(cc, skew, tail)
  }
  else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov(cc,covmatrix$model, nuisance)}
}
######################## skew student T#####################################
if(covmatrix$model==18) 
{
                  vv=as.numeric(covmatrix$param['sill']); nu=1/as.numeric(covmatrix$param['df']);sk=as.numeric(covmatrix$param['skew'])
            if(is.null(covmatrix$copula)){ 
                  sk2=sk*sk;l=nu/2; f=(nu-1)/2; w=sqrt(1-sk2);y=rho;
                  CorSkew=(2*sk2/(pi*w*w+sk2*(pi-2)))*(sqrt(1-y*y)+y*asin(y)-1)+w*w*cc/(w*w+sk2*(1-2/pi)) ;
                  corri=(pi*(nu-2)*gamma(f)^2/(2*(pi*gamma(l)^2-sk2*(nu-2)*gamma(f)^2)))*(Re(hypergeo::hypergeo(0.5,0.5,l,y*y))*((1-2*sk2/pi)*CorSkew+2*sk2/pi)-2*sk2/pi);
              }
                else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov(cc,covmatrix$model, nuisance)}
}
######################## # two piece Student T#####################################
         if(covmatrix$model==27) {  
             nu=as.numeric(1/covmatrix$param['df']);
             sk=as.numeric(covmatrix$param['skew']);
             vv=as.numeric(covmatrix$param['sill'])
            if(is.null(covmatrix$copula)){ 
             corr2=cc^2;sk2=sk^2
             a1=Re(hypergeo::hypergeo(0.5,0.5,nu/2,corr2))
             a2=cc*asin(cc) + (1-corr2)^(0.5)
             ll=qnorm((1-sk)/2)
             p11=pbivnorm::pbivnorm(ll,ll, rho = rho, recycle = TRUE)
             a3=3*sk2 + 2*sk + 4*p11 - 1
             KK=( nu*(nu-2)*gamma((nu-1)/2)^2) / (nu*pi*gamma(nu/2)^2*(3*sk2+1)-4*sk2*nu*(nu-2)*gamma((nu-1)/2)^2 )
             corri= KK*(a1*a2*a3-4*sk2);
         }
          else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov(cc,covmatrix$model, nuisance)}
                      }
########################### two piece Tukey h ##################################
         if(covmatrix$model==38) {  
                        tail=as.numeric(covmatrix$param['tail']);sk=as.numeric(covmatrix$param['skew']);vv=as.numeric(covmatrix$param['sill'])
 if(is.null(covmatrix$copula)){ 
                        corr2=cc^2;sk2=sk^2;
                        gg2=(1-(1-corr2)*tail)^2
                        xx=corr2/gg2
                        A=(asin(sqrt(xx))*sqrt(xx)+sqrt(1-xx))/(1-xx)^(1.5)
                        ll=qnorm((1-sk)/2)
                        p11=pbivnorm::pbivnorm(ll,ll, rho = rho, recycle = TRUE)
                        a3=3*sk2 + 2*sk + 4*p11 - 1
                        mm=8*sk2/(pi*(1-tail)^2);
                        ff=(1+3*sk2)/(1-2*tail)^(1.5)
                        M=(2*(1-corr2)^(3/2))/(pi*gg2)
                        corri=  (M*A*a3-mm)/( ff- mm)
                    }
                     else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov(cc,covmatrix$model, nuisance)}
                      }
############################################################
         if(covmatrix$model==39) {  # bimodal
                                  nu=as.numeric(nuisance['df']); sk=as.numeric(nuisance['skew'])
                                  vv=as.numeric(covmatrix$param['sill']); delta=as.numeric(nuisance['shape'])
                                  alpha=2*(delta+1)/nu; nn=2^(1-alpha/2)
                                  ll=qnorm((1-sk)/2)
                                  p11=pbivnorm::pbivnorm(ll,ll, rho = rho, recycle = TRUE)
                                  corr2=cc^2;sk2=sk^2
                                  a1=Re(hypergeo::hypergeo(-1/alpha ,-1/alpha,nu/2,corr2))
                                  a3=3*sk2 + 2*sk + 4*p11 - 1
                                  MM=(2^(2/alpha)*(gamma(nu/2 + 1/alpha))^2)
                                  vari=2^(2/alpha)*(gamma(nu/2 + 2/alpha))*gamma(nu/2)* (1+3*sk2) - sk2*2^(2/alpha+2)*gamma(nu/2+1/alpha)^2
                                  corri= MM*(a1*a3-4*sk2)/vari
                      }
############################################################
        if(covmatrix$model==29) 
{  # two piece Gaussian
             if(is.null(covmatrix$copula)){ 
                          corr2=sqrt(1-cc^2)
                          vv=as.numeric(covmatrix$param['sill'])
                          sk=as.numeric(nuisance['skew']); sk2=sk^2
                          ll=qnorm((1-sk)/2)
                          p11=pbivnorm::pbivnorm(ll,ll, rho = rho, recycle = TRUE)
                          KK=3*sk2+2*sk+ 4*p11 - 1
                          corri=(2*((corr2 + cc*asin(cc))*KK)- 8*sk2)/(3*pi*sk2  -  8*sk2   +pi   )
                      }
                      else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov(cc,covmatrix$model, nuisance)}
}
############################################################
       if(covmatrix$model==28)
{

         shape1=as.numeric(covmatrix$param['shape1']);shape2=as.numeric(covmatrix$param['shape2']);

             if(is.null(covmatrix$copula)){ 
          corr2=cc^2       
         cc1=0.5*(shape1+shape2); vv=shape1*shape2/((cc1+1)*(shape1+shape2)^2)
         idx=which(abs(corr2)>1e-10);corr22=corr2[idx]; nu2=shape1/2;alpha2=shape2/2
         A <- 0; res <- 0; ss <- 0; kk <- 0
         while(kk<=100){
             p1=2*(lgamma(cc1+kk)-lgamma(cc1)+lgamma(nu2+1+kk)-lgamma(nu2+1));p2=lgamma(kk+1)+(lgamma(nu2+kk)-lgamma(nu2))+2*(lgamma(cc1+1+kk)-lgamma(cc1+1))
             b1=p1-p2; b2=log(hypergeo::genhypergeo(U=c(cc1+kk,cc1+kk,alpha2), L=c(cc1+kk+1,cc1+kk+1), polynomial=TRUE,maxiter=1000, z=corr22)); b3=kk*log(corr22)
             sum=exp(b1+b2+b3); res=res+sum
             if (all(sum<1e-6)){ break} else{ A=res}
         kk=kk+1
        }
         cc[idx]=A; corri=shape1*(cc1 + 1 ) * ((1-corr2)^(cc1) *cc -1)/shape2
         corri[-idx]=0
    }
       else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov(cc,covmatrix$model, nuisance)}
}
############################################################    

if(covmatrix$model==21) 
{ # gamma
                        sh=as.numeric(covmatrix$param["shape"])
                         if(is.null(covmatrix$copula)){  corri=cc^2 }
                         else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov(cc,covmatrix$model, nuisance)}
}
if(covmatrix$model==26) 
{  # weibull
                        sh=as.numeric(covmatrix$param['shape']);
                        if(is.null(covmatrix$copula)){ 
                        corri <- .GeoWeibullCorr(cc, sh)
                         }
                        else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov(cc,covmatrix$model, nuisance)}           
}
if(covmatrix$model==24) 
{  # loglogistic
                        sh=as.numeric(covmatrix$param['shape'])
                             if(is.null(covmatrix$copula)){ 
                            cc2 <- cc^2 ; num <- pi * sin(2 * pi / sh);den <- 2 * sh * (sin(pi / sh))^2 - pi * sin(2 * pi / sh);factor <- num / den
                            F1 <- hypergeo::hypergeo(-1/sh, -1/sh, 1, cc2)
                            F2 <- hypergeo::hypergeo( 1/sh,  1/sh, 1, cc2)
                            corri <- factor * (Re(F1) * Re(F2) - 1)
                        }
                        else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov(cc,covmatrix$model, nuisance)}
}
if(covmatrix$model==22)
 {  # loggaussian
                        ss=as.numeric(covmatrix$param['sill']);
                               if(is.null(covmatrix$copula)){  corri= (exp(ss*cc)-1)/(exp(ss)-1)}
                                     else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov(cc,covmatrix$model, nuisance)}
}
############################################################
} ## end legacy/non-copula marginal-correlation transformations
 }        
############################################################
if(bivariate)  {          if(which==1)    vvar=covmatrix$param["sill_1"]+covmatrix$param["nugget_1"]
                          if(which==2)    vvar=covmatrix$param["sill_2"]+covmatrix$param["nugget_2"]
               }
else    {
    ##### computing target marginal variance
    if (copula_cont) {
      vvar <- as.numeric(moments_loc_cop$variance)
      M <- 0
    } else {
    ##### legacy/non-copula mean and variances for each model
     if(covmatrix$model==1)   {vvar= vv;                M=0 }
     if(covmatrix$model==10)  {vvar= vv+sk^2*(1-2/pi); M=sk*sqrt(2/pi)}
     if(covmatrix$model==59)  {vvar = vv * (1 - 2*sk + 2*sk^2) / (sk^2 * (1-sk)^2)
                               M= sqrt(vv) * (1/sk - 1/(1-sk))
                              }
     if(covmatrix$model==12)  {vvar= vv*nu/(nu-2) ;     M=0 }
     if(covmatrix$model==34) {vvar= vv*(1-2*h)^(-1.5); M=0 }
     if(covmatrix$model==40) { vvar= vv* vv1 ; M=sqrt(vv)*mm}
     if(covmatrix$model==20) { vvar= vv* vv1; M=sqrt(vv)*mm}
   ##############################################################
     if(covmatrix$model==18)  {
                               D1=(nu-1)*0.5; D2=nu*0.5;
                               mm=sqrt(nu)*gamma(D1)*sk/(sqrt(pi)*gamma(D2));
                               vvar=vv*(nu/(nu-2)-mm^2); M=sqrt(vv)*mm
                              }
     if(covmatrix$model==27)  {
                              ttemp=gamma(0.5*(nu-1))/gamma(0.5*nu)
                              vvar= vv*((nu/(nu-2))*(1+3*sk2) - 4*sk2*(nu/pi)*ttemp^2); M= - sqrt(vv)*2*sk*sqrt(nu/pi)*ttemp
                              }
     if(covmatrix$model==28)  {
                              ssup=as.numeric((covmatrix$param["max"]-covmatrix$param["min"]))
                              vvar=(ssup^2)*(shape1*shape2/((cc1+1)*(shape1+shape2)^2))
                              M=ssup*shape1/(shape1+shape2)+as.numeric(covmatrix$param["min"])
                              muloc=0;mu=0
                              }
     if(covmatrix$model==29)  {vvar= vv*((1+3*sk2) - 8*sk2/pi ) ;M=-sqrt(vv)*2*sk*sqrt(2/pi)}
     if(covmatrix$model==38)  {vvar= vv*(ff - mm);M= -sqrt(vv)*2*sk*sqrt(2/pi)/(1-tail)}
     if(covmatrix$model==39)  {vvar= vv*vari/(nn^(2/alpha)*gamma(nu/2)^2);M=-sqrt(vv)*sk*2^(1/alpha+1)*gamma(nu/2+1/alpha)/(nn^(1/alpha)*gamma(nu*0.5)) }
     if(covmatrix$model==21)  {vvar= 2/sh }
     if(covmatrix$model==24)  {vvar= 2*sh*sin(pi/sh)^2/(pi*sin(2*pi/sh))-1 }
     if(covmatrix$model==26)  {vvar= gamma(1+2/sh)/gamma(1+1/sh)^2-1   }
     if(covmatrix$model==22)  {vvar= exp(ss)-1   }

     if(SinhAsinhtemp||Tukeyh2temp||Tukeyhtemp||logGausstemp) vvar=1
    } ## end legacy target variance

 }
###############################################################  

################################################################
############setting mean for the external mean case ################
################################################################
   if(!bivariate){
    if (copula_cont) {
      mu <- as.numeric(moments_obs_cop$mean)
      muloc <- as.numeric(moments_loc_cop$mean)
    } else {
      mu <- if (is.null(MM)) as.numeric(X %*% betas) else as.numeric(MM)
      muloc <- if (!is.null(Mloc)) as.numeric(Mloc) else as.numeric(Xloc %*% betas)
    }
                        }          
    else {
      i1 <- seq_len(covmatrix$ns[1L])
      i2 <- covmatrix$ns[1L] + seq_len(covmatrix$ns[2L])
      X11 <- X[i1, , drop = FALSE]
      X22 <- X[i2, , drop = FALSE]
      mu <- as.numeric(c(X11 %*% betas1, X22 %*% betas2))
      muloc <- .GeoBivariatePredictionMean(
        nuisance, Xloc, dimat2, which, num_betas, Mloc,
        context = "GeoKrig"
      )
    }
#### updating mean
        if(!bivariate && !copula_cont){
          #### additive model on the real line
         if(covmatrix$model %in% c(10,18,29,27,38,39,28, 34,40,20,59))
                                { muloc=muloc + M;   mu=mu + M }     
          ### multiplicative model on the positive real line
          if( (covmatrix$model %in% c(21,26,24,22)))  {emuloc=exp(muloc);emu=exp(mu) }
         }
         else{}

##################################################################
##########computing kriging weights##################################
##################################################################
if(bivariate) {
  ## Corr_c calls the bivariate covariance model directly; its output already
  ## contains sill, nugget and cross-covariance scaling.
  CC <- matrix(corri, nrow = dimat, ncol = dimat2)
} else if(is.null(covmatrix$copula)) {
  CC <- matrix(corri*vvar, nrow=dimat, ncol=dimat2)
} else if(covmatrix$copula %in% c("Gaussian", "Clayton", "SkewGaussian")) {
  CC <- matrix(corri, nrow=dimat, ncol=dimat2)
}
if(identical(type_krig, "Universal")) {
  if(bivariate) {
    X_system <- .GeoBivariateBlockDesign(X, covmatrix$ns)
    Xloc_system <- .GeoBivariateTargetBlockDesign(
      Xloc, nloc = dimat2, p = num_betas, which = as.integer(which)
    )
    mean_names_system <- colnames(X_system)
  } else {
    if (copula_cont) {
      d_obs <- .GeoCopulaMarginalMeanDerivative(covmatrix$model, nuisance, eta_obs_cop)
      d_loc <- .GeoCopulaMarginalMeanDerivative(covmatrix$model, nuisance, eta_loc_cop)
      X_system <- X * d_obs
      Xloc_system <- Xloc * d_loc
    } else {
      X_system <- X
      Xloc_system <- Xloc
    }
    mean_names_system <- .GeoMean_names(ncol(X_system))
  }
} else {
  X_system <- Xloc_system <- mean_names_system <- NULL
}
krig_system <- .GeoKrigSystem(
  covmatrix = covmatrix, CC = CC, method = method,
  type_krig = type_krig, X = X_system, Xloc = Xloc_system,
  target_variance = vvar, mse = mse, varcov = varcov,
  mean_names = mean_names_system, context = "GeoKrig"
)
krig_weights = krig_system$weights
BB = NULL

##################################################################
################# simple or universal kriging #################################
##################################################################
if(type_krig %in% c("Simple","Universal"))  {

if(!bivariate)
{  

if (copula_cont) {
  datas <- as.matrix(as.numeric(dataT) - as.numeric(mu))
  pp <- as.numeric(muloc) + as.numeric(crossprod(datas, krig_weights))
} else {

if(type_krig %in% c("Simple","Universal")) {
               ############################ optimal linear predictors #######################
               if(covmatrix$model %in% c(1,12,27,38,29,10,18,39,37,28,9, 34,40,20,59))
                {                        
                      datas=as.matrix(c(dataT)-c(mu))
                      pp = c(muloc)      +  crossprod(datas,krig_weights)
                }
             }
               ###################################################
               #### gamma weibull loglogistic loggaussian 
if(covmatrix$model %in% c(21,24,26,22))
                      {       
                 ones=rep(1,length(c(dataT))); one=rep(1,length(c(muloc)))    
                 datas=as.matrix(c(dataT)/emu-ones)
                 pp = c(emuloc) * ( one + crossprod(datas,krig_weights))  
                      }
} ## end legacy/non-copula predictor
}     
else {
  dat <- as.numeric(dataT) - as.numeric(mu)
  pp <- as.numeric(muloc) + as.numeric(crossprod(dat, krig_weights))
}
#### MSE COMPUTATION ##################

  if(mse) {
    if (copula_cont) {
      vv <- krig_system$mse
    } else {
      # Gaussian, StudentT, skew-Gaussian and two-piece linear kriging
      if(covmatrix$model %in% c(1,12,27,38,29,10,18,39,28,40,34,9,20,59))
        vv <- krig_system$mse

      # gamma, Weibull, log-logistic and log-Gaussian
      if(covmatrix$model %in% c(21,26,24,22))  {
        vv <- krig_system$mse
        vv <- emuloc^2 * vv
      }
    }
  }     


  }  

######################################################
####formatting data ##################################
######################################################
cvv=c(vv)

##### setting almost zero when zero mse
if(mse) { if(anyNA(cvv)) 
          {
              selna=is.na(cvv)
              cvv[selna]=1e-30
          }
          selp=(cvv<=0); 
          if(any(selp)) {cvv[selp]=1e-30}
    }
####   formatting data  
if(spacetime||bivariate) {
            pred=matrix(t(pp),nrow=tloc,ncol=numloc);
            varpred=if(isTRUE(mse)) matrix(cvv,nrow=tloc,ncol=numloc) else NULL;
            }
else {pred=c(pp);varpred=if(isTRUE(mse)) cvv else NULL}

}

} 

####################################################################################################################################
###################### binomial  negative binomial and poisson (inflated) #####################################
####################################################################################################################################
if(!copula_cont && covmatrix$model %in% c(2,11,14,19,30,36,16,43,44,45,46,47,57,58))
{

     if(type=="Standard"||type=="standard") {

      if(!bivariate){
               mu <- if (is.null(MM)) as.numeric(X %*% betas) else as.numeric(MM)
               mu0 <- if (!is.null(Mloc)) as.numeric(Mloc) else as.numeric(Xloc %*% betas)
      }
     if(bivariate)  {mu  = c(X11%*%betas1,X22%*%betas2)}
     ## Count-size vectors in observation and prediction ordering.
     n_obs_count <- if (is.null(count_sizes)) rep.int(1L, dimat) else count_sizes$obs
     n_pred_count <- if (is.null(count_sizes)) rep.int(1L, dimat2) else count_sizes$pred
     KK <- n_pred_count

## covariance between observed and prediction sites
if(covmatrix$model %in% c(2,11,14,16,19,30,36,43,44,45,46,47,57,58))
{ 
cop=0  
if(!is.null(covmatrix$copula)){ if(covmatrix$copula=="Gaussian") cop=1  }

ccorr = dotCall64::.C64(
  "Corr_c_bin",
  SIGNATURE = c(
    rep("double", 5),
    "integer", "integer", "double", "double", "double",
    rep("integer", 9),
    "double", "double",           # mean_obs, mean_loc
    "double", "double",           # other_nuis, corrparam
    "integer", "integer",         # spacetime, bivariate
    "double",                     # time
    "integer", "integer",         # distance, which
    "double",                     # radius
    "integer"                     # cop
  ),
  corri = dotCall64::vector_dc("double", dimat * dimat2),
  ccc[,1], ccc[,2], ccc[,3], covmatrix$coordt,
  corrmodel, as.integer(0),
  locx, locy, locz,
  covmatrix$numcoord, numloc,
  covmatrix$model, tloc,
  as.integer(KK), as.integer(n_obs_count), covmatrix$ns, NS, covmatrix$numtime,
  as.double(c(mu)),              # mean_obs
  as.double(c(mu0)),             # mean_loc
  as.double(other_nuis), as.double(corrparam),
  as.integer(covmatrix$spacetime), as.integer(bivariate),
  as.double(time), as.integer(distance), as.integer(which - 1),
  as.double(covmatrix$radius), as.integer(cop),
  INTENT = c("w", rep("r", 29)),
  PACKAGE = "GeoModels",
  VERBOSE = 0,
  NAOK = TRUE
)
}

if(cop==1) { 
nuisance$n=n_obs_count[1L]
corri= gaussian_copula_cov(ccorr$corri,covmatrix$model, nuisance)
     }
else {corri=ccorr$corri}

##################inverse of cov matrix##############################################
        if(!bivariate) cc = matrix(corri,nrow=dimat,ncol=dimat2)
        else           {}
        krig_weights = .GeoKrigSolve(covmatrix, cc, method=method,
                                         context="GeoKrig")
        BB = NULL
######################################################################################

       if(type_krig %in% c("Simple","Universal"))   {

      ##########################################################
       if(covmatrix$model==30||covmatrix$model==36){  ### poisson
        p0=exp(mu0); pmu=exp(mu)
            if(!bivariate) { 
            datas=c(dataT)-c(pmu)

            pp = c(p0) + crossprod(datas, krig_weights)
             }
            else{}
           if(mse)  vvar=p0
          }
      ##########################################################
      if(covmatrix$model==46||covmatrix$model==47){  ### poisson gamma
         p0=exp(mu0); pmu=exp(mu)
         if(!bivariate) { 
               datas=c(dataT)-c(pmu)
               pp = c(p0) + crossprod(datas, krig_weights)
          }
        else{}
        if(mse)  vvar=p0*(1+p0/covmatrix$param['shape'])
          }
      ##########################################################
       if(covmatrix$model==43||covmatrix$model==44){  ### poisson inflated
        p=pnorm(covmatrix$param['pmu'])
        p0=exp(mu0); pmu=exp(mu)
            if(!bivariate)  { 
             datas=c(dataT)-(1-p)*c(pmu)
             pp = (1-p)*c(p0) + crossprod(datas,krig_weights)
             }
            else{}
           if(mse)  vvar=(1-p)*p0*(1+p*p0)
        }
       ##########################################################
       if(covmatrix$model==2||covmatrix$model==11){  ### binomial
        p0=pnorm(mu0); pmu=pnorm(mu)
            if(!bivariate) { 
               datas=c(dataT)-n_obs_count*c(pmu)
               pp = n_pred_count*c(p0) + crossprod(datas,  krig_weights) 
            }
            else{}
           if(mse) vvar=KK*p0*(1-p0)
          }
      ###########################################################
           if(covmatrix$model==19){  ### binomial2
            p0=pnorm(mu0); pmu=pnorm(mu)
            if(!bivariate)   
              { datas=c(dataT)-c(n_obs_count*pmu)
                pp = n_pred_count*c(p0) + crossprod(datas,krig_weights)
              }
            else{}
          if(mse)    vvar=n_pred_count*p0*(1-p0)
          }
        ##########################################################
       if(covmatrix$model==14||covmatrix$model==16){    ###geometric or negative binomial
         p0=pnorm(mu0); pmu=pnorm(mu)
            if(!bivariate)
            {   k1=c(p0);k2=c(pmu);
                aa=n_pred_count*(1-k1)/k1
                datas=c(dataT)-n_obs_count*(1-k2)/k2
                pp = aa + crossprod(datas ,krig_weights)
          }
            else{}
            if(mse) vvar=aa/k1
          }
            ##########################################################
         if(covmatrix$model==57||covmatrix$model==58){  ### poisson gamma inflated
         p=pnorm(covmatrix$param['pmu'])
         p0=exp(mu0); pmu=exp(mu)
         if(!bivariate) { 
               datas=c(dataT)-(1-p)*c(pmu)
               pp = (1-p)*c(p0) + crossprod(datas, krig_weights)
          }
        else{}
        if(mse)  vvar=p0*(1+ (p0*(1+covmatrix$param['shape']*p)/((1-p)*covmatrix$param['shape'])) )
        }
        ##########################################################
        if(covmatrix$model==45){    ###inflated negative binomial
            p0=pnorm(mu0); pmu=pnorm(mu)
            p=as.numeric(pnorm(covmatrix$param['pmu']))
            if(!bivariate) { k1=c(p0);k2=c(pmu);
                             datas=c(dataT)-(1-p)*n_obs_count*(1-k2)/k2
                             pp = (1-p)*n_pred_count*(1-k1)/k1 + crossprod(datas,krig_weights)
                           }
            else{}
            if(mse) vvar=n_pred_count*(1-k1)*(1-p)*(1+n_pred_count*p*(1-k1))/k1^2
          }
        if(mse){
          ## Compute only the diagonal kriging variance; never allocate an
          ## Npred x Npred matrix.  Universal kriging additionally propagates
          ## uncertainty in the estimated regression coefficients.
          vv <- as.numeric(vvar) - colSums(cc * krig_weights)
          if (identical(type_krig, "Universal")) {
            d_obs <- .GeoKrigCountMeanDerivative(
              covmatrix$model, covmatrix$param, mu, n_obs_count
            )
            d_loc <- .GeoKrigCountMeanDerivative(
              covmatrix$model, covmatrix$param, mu0, n_pred_count
            )
            X_system_count <- X * d_obs
            Xloc_system_count <- Xloc * d_loc
            mean_names_count <- .GeoMean_names(ncol(X_system_count))
            var_mean_count <- .GeoKrigMeanVarcov(
              varcov, mean_names = mean_names_count, context = "GeoKrig"
            )
            trend_gap_count <- Xloc_system_count - t(krig_weights) %*% X_system_count
            vv <- vv + rowSums((trend_gap_count %*% var_mean_count) * trend_gap_count)
          }
        }
        }
##############################################################
if(!mse) vv <- NULL
cvv=c(vv)
##### setting almost zero when zero mse
if(mse) {
   if(anyNA(cvv)) 
     {selna=is.na(cvv)
     cvv[selna]=1e-30
    }
    selp=(cvv<=0); 
    if(any(selp)) {cvv[selp]=1e-30}
    }
####   formatting data  
if(spacetime||bivariate) {
            pred=matrix(t(pp),nrow=tloc,ncol=numloc);
            varpred=if(isTRUE(mse)) matrix(cvv,nrow=tloc,ncol=numloc) else NULL;
            }
          else {pred=c(pp);varpred=if(isTRUE(mse)) cvv else NULL}
    }}

if(tloc==1)  {pred=c(pred);varpred=c(varpred)}
if(!isTRUE(mse)) varpred <- NULL

    # Return the objects list:
   GeoKrig= list(     bivariate=bivariate,
                   coordx = covmatrix$coordx,coordy = covmatrix$coordy,coordz = covmatrix$coordz,
                   coordt = covmatrix$coordt,coordx_dyn=covmatrix$coordx_dyn,
                   covmatrix=covmatrix$covmatrix,corrmodel = corrmodel,copula=copula,data=data,distance = distance,
                   grid=covmatrix$grid,loc=loc_orig,nozero=covmatrix$nozero,ns=covmatrix$ns,X=XX,
                   numcoord = covmatrix$numcoord,numloc= numloc,numtime = covmatrix$numtime,numt = tloc,
                   model=model,fit_model=model_requested,param = param_input,pred=pred,n=n_input,nloc=nloc_input,
                   radius=radius,spacetime = covmatrix$spacetime,time=time,type=type,type_krig=type_krig,method=method,weigthed=weigthed,varcov=varcov,
                   weights=krig_weights,mse=varpred)
    structure(c(GeoKrig , call = call), class = c("GeoKrig"))
}

}