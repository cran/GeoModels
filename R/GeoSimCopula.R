####################################################
### File name: GeoSimCopula.r
####################################################


# Simulate spatial and spatio-temporal random felds:
GeoSimCopula <- function(coordx=NULL, coordy=NULL, coordz=NULL,coordt=NULL, coordx_dyn=NULL,corrmodel, distance="Eucl", grid=FALSE,
     method="cholesky",model='Gaussian', n=1, param,anisopars=NULL, radius=1, sparse=FALSE,
     copula="Gaussian",X=NULL,spobj=NULL,nrep=1,progress=FALSE,check.duplicates=FALSE,
     L=10000,parallel=FALSE,ncores=6)
{
call <- match.call()

    corrmodel <- .GeoSimulationCanonicalCorrModel(corrmodel)
    if(is.null(CkCorrModel(corrmodel))) stop("The name of the correlation model is not correct\n")
    model_requested <- gsub("[[:blank:]]", "",model)
    model <- .GeoSimulationCanonicalModel(model_requested)
    if(is.null(CkModel(model))) stop("The name of the model is not correct\n")
    .GeoValidateCopulaSimulationModel(model)
    n <- .GeoSimulationEffectiveN(model_requested, model, n)
    nrep <- .GeoValidateNrep(nrep)
    distance=gsub("[[:blank:]]", "",distance)
    method <- .GeoValidateCopulaSimulationMethod(method)
    copula=gsub("[[:blank:]]", "",as.character(copula)[1L])
    param <- lapply(param, unname)
    param_input <- param
    .GeoValidateCopulaMarginal(model, param)

##############################################################################
###### extracting sp object informations if necessary              ###########
##############################################################################
bivariate<-CheckBiv(CkCorrModel(corrmodel))
spacetime<-CheckST(CkCorrModel(corrmodel))
space=!spacetime&&!bivariate
if(bivariate) {
  stop("GeoSimCopula currently supports only univariate correlation models.",
       call. = FALSE)
}

.GeoCopula_flatten <- function(z, spacetime=FALSE, grid=FALSE) {
  if(is.list(z)) return(unlist(z, use.names=FALSE))
  if(is.matrix(z) && isTRUE(spacetime) && !isTRUE(grid))
    return(as.numeric(t(z)))
  as.numeric(z)
}
.GeoCopula_restore <- function(z, template, spacetime=FALSE, grid=FALSE) {
  z <- as.numeric(z)
  if(is.list(template)) {
    lens <- lengths(template)
    if(sum(lens) != length(z))
      stop("GeoSimCopula: simulated data length does not match dynamic coordinate blocks.",
           call. = FALSE)
    ends <- cumsum(lens)
    starts <- c(1L, utils::head(ends, -1L) + 1L)
    return(Map(function(a, b) z[a:b], starts, ends))
  }
  if(is.matrix(template) && isTRUE(spacetime) && !isTRUE(grid))
    return(matrix(z, nrow=nrow(template), ncol=ncol(template), byrow=TRUE))
  if(!is.null(dim(template))) return(array(z, dim=dim(template)))
  z
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
if(identical(method, "TB")) {
  L <- .GeoValidateTBControls(L)
  .GeoValidateTBDistance(distance, context = "GeoSimCopula")
  if(isTRUE(spacetime)) {
    stop("GeoSimCopula(method='TB') does not support spatio-temporal correlation models; use method = 'cholesky' or 'svd'.",
         call. = FALSE)
  }
  .GeoValidateTBCorrelationModel(corrmodel, bivariate = FALSE)
}
###############################################################
###############################################################

if((copula!="Clayton")&&(copula!="Gaussian")&&(copula!="AMH")&&(copula!="SkewGaussian")) stop("the type of copula is wrong")

## Use exact parameter-name matching. In R, param$nu may partially match
## param$nugget when nu is absent, which can silently corrupt copula checks.
has_param <- function(nm) !is.null(names(param)) && nm %in% names(param)
get_param <- function(nm) {
  if (!has_param(nm)) return(NULL)
  param[[nm, exact = TRUE]]
}

if(copula=="Clayton" && !has_param("nu")) stop("Clayton copula needs a nu parameter")
if(copula=="AMH" && !has_param("nu")) stop("AMH copula needs a nu parameter")
if(copula=="AMH") {
  nu_amh <- as.numeric(get_param("nu"))
  if(length(nu_amh) != 1L || !is.finite(nu_amh))
    stop("AMH nu must be a finite scalar.", call.=FALSE)
}
if(copula=="SkewGaussian" && !has_param("nu")) stop("Skew Gaussian copula needs a nu parameter")
if(!has_param("nugget")) stop("GeoSimCopula requires an explicit param$nugget.", call.=FALSE)

nugget_cop <- as.numeric(get_param("nugget"))
if(length(nugget_cop) != 1L || !is.finite(nugget_cop) ||
   nugget_cop < 0 || nugget_cop >= 1)
  stop("GeoSimCopula requires param$nugget to be a finite scalar in [0,1).", call.=FALSE)

#### corr parameters
paramcorr=param[CorrParam(corrmodel)]

.GeoCheckDuplicateLocations(
  coordx = coordx, coordy = coordy, coordz = coordz, coordt = coordt,
  coordx_dyn = coordx_dyn, spacetime = spacetime, bivariate = FALSE,
  grid = grid, context = "GeoSimCopula", check.duplicates = check.duplicates
)
.dup_old <- .GeoSuppressDuplicateChecksForNestedCalls()
on.exit(.GeoRestoreDuplicateCheckOption(.dup_old), add = TRUE)

## Keep copula construction separate from the numerical latent-field engine.
## Exact methods use GeoSim(), while turning bands uses GeoSimapprox().
.GeoCopula_latent_sim <- function(latent_model, latent_param, n_latent = 1L) {
  common_args <- list(
    coordx = coordx, coordy = coordy, coordz = coordz, coordt = coordt,
    coordx_dyn = coordx_dyn, corrmodel = corrmodel, distance = distance,
    grid = grid, model = latent_model, n = n_latent, param = latent_param,
    anisopars = anisopars, radius = radius, nrep = 1L,
    progress = FALSE, check.duplicates = FALSE
  )

  if(identical(method, "TB")) {
    return(do.call(
      GeoSimapprox,
      c(common_args, list(method = "TB", L = L,
                          parallel = parallel, ncores = ncores))
    ))
  }

  do.call(
    GeoSim,
    c(common_args, list(method = method, sparse = sparse))
  )
}

SIM=list()
###################################################
### starting number of replicates
###################################################
for(rep_idx in seq_len(nrep)){



####Gaussian copula #############################################
if(copula=="Gaussian")
{
param1=c(list(mean=0,sill=1,nugget=get_param("nugget")),paramcorr)

sim <- .GeoCopula_latent_sim("Gaussian", param1, n_latent = 1L)
latent_template <- sim$data
latent <- .GeoCopula_flatten(latent_template, sim$spacetime, sim$grid)
unif <- pnorm(latent, mean=0, sd=1)

}

####skewGaussian copula #############################################
if(copula=="SkewGaussian")
{

nu <- as.numeric(get_param("nu"))
if(length(nu) != 1L || !is.finite(nu) || abs(nu) >= 1)
  stop("nu parameter must be finite and between -1 and 1")
alpha <- nu / sqrt(1 - nu^2)

## The skew-Gaussian copula is generated from
##   Z(s) = alpha * abs(G1(s)) + G2(s),
## where G1 and G2 are independent Gaussian RFs sharing exactly the
## same correlation function.  Simulating the two latent fields
## explicitly also applies the nugget to both of them, consistently
## with the copula likelihood and its latent-field definition.
param_latent <- c(list(mean=0, sill=1, nugget=get_param("nugget")), paramcorr)
sim_g1 <- .GeoCopula_latent_sim("Gaussian", param_latent, n_latent = 1L)
sim_g2 <- .GeoCopula_latent_sim("Gaussian", param_latent, n_latent = 1L)
latent_template <- sim_g2$data
g1 <- .GeoCopula_flatten(sim_g1$data, sim_g1$spacetime, sim_g1$grid)
g2 <- .GeoCopula_flatten(sim_g2$data, sim_g2$spacetime, sim_g2$grid)
latent <- alpha * abs(g1) + g2
omega <- sqrt(alpha^2 + 1)
unif <- sn::psn(latent, xi=0, omega=omega, alpha=alpha)
sim <- sim_g2
}


#if (copula == "SkewGaussian") {
#  delta <- param$nu  
#  print("h")
#  print(delta)
#  eps   <- 1e-6;
#  alpha <- delta / sqrt(1 - delta^2)  
#  print(alpha) 
#  skew_eta <- alpha
#  param1 <- c(list(mean   = 0,sill   = 1,nugget = param$nugget,skew   = skew_eta),paramcorr)
#  sim <- GeoSim(
#    coordx = coordx, coordy = coordy,coordz = coordz, coordt = coordt,
#    coordx_dyn = coordx_dyn,corrmodel = corrmodel,distance  = distance, grid = grid,
#    method = method, model = "SkewGaussian",
#    n = 1, param = param1,anisopars = anisopars, radius = radius,sparse = sparse, nrep = 1
#  )
#  omega <- sqrt(1 + alpha^2)
#  unif  <- sn::psn(sim$data, xi = 0, omega = omega, alpha = alpha)
#}

####beta copula #############################################
if(copula=="Clayton")
{
nu_clayton <- as.numeric(get_param("nu"))
if(length(nu_clayton) != 1L || !is.finite(nu_clayton) || nu_clayton <= 0 ||
   abs(nu_clayton - round(nu_clayton)) > sqrt(.Machine$double.eps)) {
  stop("GeoSimCopula: the constructive Clayton simulator requires nu to be a positive integer; nu is not rounded automatically.",
       call. = FALSE)
}
pp <- as.integer(round(nu_clayton))
param1=c(list(shape1=pp,shape2=2,mean=0,min=0,max=1,nugget=get_param("nugget")),paramcorr)
sim <- .GeoCopula_latent_sim("Beta", param1, n_latent = 1L)
latent_template <- sim$data
latent <- .GeoCopula_flatten(latent_template, sim$spacetime, sim$grid)
unif <- latent^(pp/2)
}
####AMH copula #############################################
if(copula=="AMH")
{

param1=c(list(mean=0,nugget=as.numeric(get_param("nugget")),shape=2),paramcorr)
sim <- .GeoCopula_latent_sim("Gamma", param1, n_latent = 1L)
param2=c(list(mean=as.numeric(get_param("nu")),nugget=as.numeric(get_param("nugget"))),paramcorr)
sim2 <- .GeoCopula_latent_sim("BinomialNeg", param2, n_latent = 1L)
latent_template <- sim$data
latent1 <- .GeoCopula_flatten(sim$data, sim$spacetime, sim$grid)
latent2 <- .GeoCopula_flatten(sim2$data, sim2$spacetime, sim2$grid)
if(length(latent1) != length(latent2))
  stop("GeoSimCopula: incompatible latent simulations in the AMH construction.", call.=FALSE)
kk <- latent1/(latent2+1)
pp <- pnorm(as.numeric(get_param("nu")))
unif <- pp/(exp(kk)+pp-1)
}
####################################################################

## Keep inverse-CDF evaluations away from exact 0 and 1.
eps_u <- .Machine$double.eps
unif <- pmin(pmax(as.numeric(unif), eps_u), 1 - eps_u)

if(!sim$bivariate){
           ## length(data) is valid for fixed and dynamic space-time layouts.
           nobs_mean <- length(unif)
           external_mean <- .GeoMean_external(param, n = nobs_mean)
           if(!is.null(external_mean) && !is.null(X))
             stop("A site-specific mean vector and X cannot be used together.\n")
           if(!is.null(X)) {
             if(is.list(X)) X <- do.call(rbind, lapply(X, as.matrix))
             X <- .GeoMean_design(X, n = nobs_mean, name = "X")
             invisible(.GeoMean_beta(param, ncol(X)))
           } else if(is.null(external_mean) &&
                     length(setdiff(.GeoMean_supplied_names(param), "mean"))) {
             stop("Multiple mean coefficients require an explicit design matrix X.\n")
           }
           mm <- .GeoMean_vector(param, X = X, n = nobs_mean,
                                 external = external_mean)

           # Keep param unchanged across replications.  In particular, do not
           # reset mean, mean1, ... inside the nrep loop.
    }       
##############################
##############################
#######     models     #######
##############################
simcop <- NULL

####################################
############ discrete  RF ##########
####################################
.GeoCopula_size <- function(size, nobs, model_name) {
  size <- as.numeric(size)
  if(length(size) == 1L) size <- rep(size, nobs)
  if(length(size) != nobs || any(!is.finite(size)) ||
     any(size <= 0) || any(size != round(size))) {
    stop(sprintf("For %s, n must contain positive integers and have length 1 or the number of observations.", model_name),
         call. = FALSE)
  }
  size
}
if(model=="Binomial") {
 nn <- .GeoCopula_size(n, length(unif), model_requested)
 simcop=qbinom(unif, size=nn, prob=pnorm(mm))
}
if(model=="BinomialNeg") {
 nn <- .GeoCopula_size(n, length(unif), model_requested)
 simcop=qnbinom(unif, size=nn, prob=pnorm(mm))
}
if(model=="BinomialNegZINB") {
 nn <- .GeoCopula_size(n, length(unif), model_requested)
 simcop=qzinegbin(unif, size=nn, #prob = pnorm(mm),                # require package VGAM
        munb=   nn*(1-pnorm(mm))/pnorm(mm),#,n/pnorm(mm)-n,
        pstr0 = as.numeric(param$pmu))
}
if(model=="Poisson") {
 simcop=qpois(unif, lambda=exp(mm))
}
if(model=="PoissonZIP") {
 simcop=qzipois(unif, lambda=exp(mm), pstr0 = as.numeric(param$pmu) ) # require package VGAM
}
####################################
############ positive real  RF #####
####################################
if(model=="Gamma") 
         {
p2=as.numeric(param$shape)
if(length(p2) != 1L || !is.finite(p2) || p2 <= 0)
  stop("Gamma shape must be a finite positive scalar.", call.=FALSE)
simcop=exp(mm)*qgamma(unif, shape=p2/2, scale=2/p2)
         }
############
if(model=="Weibull") 
         {
p2=as.numeric(param$shape)
if(length(p2) != 1L || !is.finite(p2) || p2 <= 0)
  stop("Weibull shape must be a finite positive scalar.", call.=FALSE)
simcop= exp(mm)*qweibull(unif,shape=p2,scale=1/(gamma(1+1/p2 )))

         }
########################
if(model=="LogLogistic") 
         {

p2=as.numeric(param$shape)
if(length(p2) != 1L || !is.finite(p2) || p2 <= 1)
  stop("LogLogistic shape must be finite and greater than 1 under the mean parametrization used by GeoSimCopula.", call.=FALSE)
cc=gamma(1+1/p2)*gamma(1-1/p2)
simcop=.GeoMarginalQLogLogistic(unif,shape = p2,scale=exp(mm)/cc)            
         }
##############################

if(model=="LogGaussian")
{ 
    vv=as.numeric(param$sill)
   if(length(vv) != 1L || !is.finite(vv) || vv <= 0)
     stop("LogGaussian sill must be a finite positive scalar.", call.=FALSE)
   simcop = qlnorm(unif, mm-vv/2, sqrt(vv))
     #simcop = qlnorm(unif, 1+mm, sqrt(vv))
}
##############################
##############################
############ Real RF #########
##############################


if(model=="Gaussian") {
         vv <- as.numeric(param$sill)
         if(length(vv) != 1L || !is.finite(vv) || vv <= 0)
           stop("Gaussian sill must be a finite positive scalar.", call.=FALSE)
         simcop=qnorm(unif,mean=mm,sd=sqrt(vv))
         }
if(model=="Logistic") 
         {
         vv <- as.numeric(param$sill)
         if(length(vv) != 1L || !is.finite(vv) || vv <= 0)
           stop("Logistic sill must be a finite positive scalar.", call.=FALSE)
         simcop=qlogis(unif,location=mm,scale=sqrt(vv))
         }
if(model=="StudentT") {
vv=as.numeric(param$sill)
dd=as.numeric(param$df)
if(length(vv) != 1L || !is.finite(vv) || vv <= 0)
  stop("StudentT sill must be a finite positive scalar.", call.=FALSE)
if(length(dd) != 1L || !is.finite(dd) || dd <= 0)
  stop("StudentT df is the reciprocal degrees-of-freedom parameter and must be finite and positive.", call.=FALSE)
simcop=mm+sqrt(vv)*qt(unif,df=1/dd)

         }
############
if(model=="SkewGaussian") 
         {
vv=as.numeric(param$sill)  
sk=as.numeric(param$skew)  
omega=as.numeric(sqrt((sk^2 + vv)/vv)) 
alpha=as.numeric(sk/(vv^0.5))         
simcop=mm+sqrt(vv)*sn::qsn(unif,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))
         }
#######################################  
if(model=="SkewStudentT")
{
  vv=as.numeric(param$sill) 
  delta <- as.numeric(param$skew)
  if(length(delta) != 1L || !is.finite(delta) || abs(delta) >= 1)
    stop("For SkewStudentT, skew must lie strictly between -1 and 1.", call.=FALSE)
  alpha <- delta/sqrt(1-delta^2)
  eps_df <- as.numeric(param$df)
  if(length(eps_df) != 1L || !is.finite(eps_df) || eps_df <= 0)
    stop("SkewStudentT df is the reciprocal degrees-of-freedom parameter and must be finite and positive.", call.=FALSE)
  dd <- 1 / eps_df
  simcop=mm+sqrt(vv)*sn::qst(unif, xi=0, omega=1, alpha=alpha, nu=dd)
}
#######################################
if(model=="SinhAsinh")
{
vv=as.numeric(param$sill) 
tail = as.numeric(param$tail) 
skew = as.numeric(param$skew)
simcop=mm+sqrt(vv)*sinh(1/tail * asinh(qnorm(unif))+skew/tail)
}
#######################################   OK
if(model=="Tukeyh")
{
vv=as.numeric(param$sill) 
tail = as.numeric(param$tail) 
uu=qnorm(unif)
simcop=mm+sqrt(vv)*uu*exp(0.5*tail*uu^2);
}
#######################################  OK
if(model=="Tukeyh2")
{
tail1 =as.numeric(param$tail1);tail2 = as.numeric(param$tail2)
vv=as.numeric(param$sill) 
simcop =mm+sqrt(vv)*.GeoMarginalQTukeyH2(unif, tail1, tail2)
}

#######################################  OK
if(model=="SkewLaplace")
{  

# Quantile della Skew/Asymmetric Laplace
sk=as.numeric(param$skew) 
vv=as.numeric(param$sill)
if(length(vv) != 1L || !is.finite(vv) || vv <= 0)
  stop("SkewLaplace sill must be a finite positive scalar.", call.=FALSE)
if(length(sk) != 1L || !is.finite(sk) || sk <= 0 || sk >= 1)
  stop("SkewLaplace skew must lie strictly between 0 and 1.", call.=FALSE)
simcop =mm+sqrt(vv)*.GeoMarginalQSkewLaplace(unif, sk)
}
#######################################  OK
if(model=="Tukeygh")
{
tail = as.numeric(param$tail);skew =  as.numeric(param$skew) 
vv=as.numeric(param$sill) 
uu=qnorm(unif)
if(abs(skew) <= sqrt(.Machine$double.eps)) {
  core <- uu*exp(0.5*tail*uu^2)
} else {
  core <- expm1(skew*uu)*exp(0.5*tail*uu^2)/skew
}
simcop=mm+sqrt(vv)*core
}
#######################################   OK
if(model=="TwoPieceGaussian")
{
skew = as.numeric(param$skew) 
vv=as.numeric(param$sill) 
simcop =mm+sqrt(vv)*.GeoMarginalQTwoPieceGaussian(unif, skew)
}

#######################################  
if(model=="TwoPieceStudentT")
{
vv=as.numeric(param$sill) 
skew = as.numeric(param$skew)
eps_df <- as.numeric(param$df)
if(length(eps_df) != 1L || !is.finite(eps_df) || eps_df <= 0)
  stop("TwoPieceStudentT df is the reciprocal degrees-of-freedom parameter and must be finite and positive.", call.=FALSE)
df   = 1/eps_df
simcop=mm+sqrt(vv)*.GeoMarginalQTwoPieceStudentT(unif, skew, df)
}
#######################################   OK
if(model=="TwoPieceTukeyh")
{
vv=as.numeric(param$sill) 
skew = as.numeric(param$skew)
tail= as.numeric(param$tail)
simcop =mm+sqrt(vv)*.GeoMarginalQTwoPieceTukeyH(unif, skew, tail)
}




####################################
############ bounded  real RF    ###
####################################
if(model=="Kumaraswamy") 
{
p1=as.numeric(param$shape1);p2=as.numeric(param$shape2)
pmin=as.numeric(param$min);pmax=as.numeric(param$max);
if(any(!is.finite(c(p1,p2))) || any(c(p1,p2) <= 0)) stop("Kumaraswamy shape parameters must be finite and positive.", call.=FALSE)
if(!is.finite(pmin) || !is.finite(pmax) || pmin >= pmax) stop("Kumaraswamy requires finite min < max.", call.=FALSE)
simcop=pmin + (pmax-pmin)*(1-(1-unif)^(1/p2))^(1/p1)
}
############
if(model=="Kumaraswamy2") 
{ # parametrization using beta median  regression
mm=1/(1+exp(-mm))
p2=as.numeric(param$shape)
pmin=as.numeric(param$min);pmax=as.numeric(param$max);
if(!is.finite(p2) || p2 <= 0) stop("Kumaraswamy2 shape must be finite and positive.", call.=FALSE)
if(!is.finite(pmin) || !is.finite(pmax) || pmin >= pmax) stop("Kumaraswamy2 requires finite min < max.", call.=FALSE)
shape2 <- log(0.5)/log1p(-mm^p2)
simcop=pmin + (pmax-pmin)*(1-(1-unif)^(1/shape2))^(1/p2)
}
############
if(model=="Beta") 
{ # parametrization using beta regression
pmin=as.numeric(param$min);pmax=as.numeric(param$max);
sh1=as.numeric(param$shape1);sh2=as.numeric(param$shape2)
if(any(!is.finite(c(sh1,sh2))) || any(c(sh1,sh2) <= 0)) stop("Beta shape parameters must be finite and positive.", call.=FALSE)
if(!is.finite(pmin) || !is.finite(pmax) || pmin >= pmax) stop("Beta requires finite min < max.", call.=FALSE)
simcop=pmin + (pmax-pmin)*qbeta(unif,shape1=sh1,shape2=sh2)
}
############
if(model=="Beta2") 
{ # parametrization using beta mean  regression

mm=1/(1+exp(-mm))
p2=as.numeric(param$shape)
pmin=as.numeric(param$min);pmax=as.numeric(param$max);
if(!is.finite(p2) || p2 <= 0) stop("Beta2 shape must be finite and positive.", call.=FALSE)
if(!is.finite(pmin) || !is.finite(pmax) || pmin >= pmax) stop("Beta2 requires finite min < max.", call.=FALSE)
simcop=pmin + (pmax-pmin)*qbeta(unif,shape1=mm*p2,shape2=(1-mm)*p2)
}





############
if(is.null(simcop))
  stop(sprintf("GeoSimCopula does not implement marginal model '%s'.", model_requested),
       call. = FALSE)
simcop <- .GeoCopula_restore(simcop, latent_template, sim$spacetime, sim$grid)
##############################




SIM[[rep_idx]]=simcop
}
#### end numrep ##############
##############################
##############################

 if(nrep==1) SIM=SIM[[1]] 

    GeoSim_Copula <- list(bivariate = sim$bivariate,
    coordx = sim$coordx,
    coordy = sim$coordy,
    coordt = sim$coordt,
    coordz=sim$coordz,
    coordx_dyn =sim$coordx_dyn,
    corrmodel = corrmodel,
    data = SIM,
    distance = sim$distance,
    grid = sim$grid,
    model = model_requested,
    method=method,
    n=n,
    numcoord = sim$numcoord,
    numtime = sim$numtime,
    param = param_input,
    radius = radius,
    spacetime = sim$spacetime,
    sparse=if(is.null(sim$sparse)) FALSE else sim$sparse,
    copula=copula,
    nrep=nrep,
    X=X)
#}
##############################################
  structure(c(GeoSim_Copula, call = call), class = c("GeoSim_Copula"))
}
