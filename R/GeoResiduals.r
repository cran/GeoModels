####################################################
### File name: GeoResiduals.r
####################################################

GeoResiduals<-function(fit)

{
call <- match.call()
if(!inherits(fit,"GeoFit"))  stop("A GeoFit object is needed as input\n")
######
extmean=FALSE
if(!fit$bivariate) {
    mmext <- .GeoMean_external(fit$fixed)
    if(!is.null(mmext)) {
        extmean=TRUE; fit$fixed$mean=0
    }
}  ## external fixed mean
###
nm=names(fit$param); nf=names(fit$fixed)
fit$param=unlist(fit$param); fit$fixed=unlist(fit$fixed)
model=fit$model        #type of model
if(fit$bivariate && !(model %in% c("Gaussian", "Gauss"))) {
    stop("GeoResiduals currently supports bivariate residuals only for Gaussian models.",
         call. = FALSE)
}
num_betas=fit$numbetas  #number of mean parameters


temp_X=fit$X

if(!fit$bivariate)
{

res1 <- NULL

## extracting mean parameters
namescorr <- CorrParam(fit$corrmodel) 
namesnuis <- NuisParam(fit$model,fit$bivariate,num_betas, copula = fit$copula)
param <- c(fit$param, fit$fixed)
namesparam<- names(param)
paramcorr <- param[namescorr]
nuisance <- param[namesnuis]
sel <- names(nuisance) %in% .GeoMean_names(num_betas)


copula=fit$copula


#################################
#### extracting data ########
#################################
if(is.list(fit$coordx_dyn)) dd=unlist(fit$data)
else dd=c(t(fit$data))

##### if X is NULL, use an intercept-only design matrix with the correct length
##### This is important for purely spatial objects, where coordt can be NULL.
if(is.null(fit$X)) {
    fit$X <- matrix(1, nrow = length(dd), ncol = 1L)
} else {
    fit$X <- .GeoMean_design(fit$X, n = length(dd), name = "X")
}

#################################
#### computing mean ########
#################################
if(extmean) {
    if(length(mmext) != length(dd))
        stop("the fixed mean vector does not match the number of observations\n")
    mu <- as.numeric(mmext)
} else {
    beta2 <- .GeoMean_beta(param, ncol(fit$X))
    mu <- as.numeric(fit$X %*% beta2)
}


############################################################
#if(!is.null(copula)){    #### copula models
#if(copula=="Clayton"||copula=="Gaussian"){
if(model=="Beta2")
             {
              mm=c(1/(1+exp(-mu)))
              sh=as.numeric(param['shape'])
              pmin=as.numeric(param['min']); pmax=as.numeric(param['max'])
              if(!is.finite(pmin) || !is.finite(pmax) || pmax <= pmin)
                  stop("Beta2 residuals require finite min < max", call. = FALSE)
              zz=(dd-pmin)/(pmax-pmin)
              res1=pbeta(zz,mm*sh,(1-mm)*sh)
             }
if(model=="Kumaraswamy2") {
             mm=c(1/(1+exp(-mu)))
             sh=as.numeric(param['shape'])
             pmin=as.numeric(param['min']); pmax=as.numeric(param['max'])
             if(!is.finite(pmin) || !is.finite(pmax) || pmax <= pmin)
                 stop("Kumaraswamy2 residuals require finite min < max", call. = FALSE)
             zz=(dd-pmin)/(pmax-pmin)
             shape2=log(0.5)/log1p(-(mm^sh))
             if(any(!is.finite(shape2)) || any(shape2 <= 0))
                 stop("Kumaraswamy2 produced invalid fitted shape parameters", call. = FALSE)
             res1=numeric(length(zz))
             res1[zz >= 1] = 1
             inside=which(zz > 0 & zz < 1)
             if(length(inside))
                 res1[inside]=1-(1-zz[inside]^sh)^shape2[inside]
             }
          #}
#}
#else {           #### non-copula models
###  positive multiplicative models
if(model %in% c("Gamma","Weibull","LogLogistic","LogGaussian"))
{  
res1=dd/exp(c(mu))
}
### additive  models  on the real line
if(model %in% c("Gaussian","Gauss","SkewGaussian","Logistic", "SkewLaplace",
               "Tukeyh","Tukeyh2","SinhAsinh","Tukeygh","Gaussian_misp_Tukeygh",
               "StudentT",  "Gaussian_misp_StudentT","Gaussian_misp_SkewStudentT","SkewStudentT",
               "TwoPieceGaussian","TwoPieceTukeyh","TwoPieceGauss","TwoPieceStudentT","TwoPieceBimodal"))
{res1=(dd-c(mu))/sqrt(as.numeric(param['sill']))}


if(model=="Gaussian_misp_Binomial")
{
    aa=pnorm(c(mu))
    hh=fit$n*aa
    res1=(dd-hh)/sqrt(hh*(1-aa))
}
if(model=="Gaussian_misp_BinomialNeg")
{
    aa=pnorm(c(mu))
    hh=fit$n*(1-aa)/aa
    res1=(dd-hh)/sqrt(hh/aa)
}
if(model=="Gaussian_misp_Poisson")
{
    aa=exp(c(mu))
    res1=(dd-aa)/sqrt(aa)
}

## Pearson residuals for the discrete marginal models described in the paper.
## The formulas use their fitted marginal means and variances; they are
## deterministic and are intended mainly for residual variograms and
## standardized diagnostic plots.
.GeoResidual_n <- function(n, nobs, model_name) {
    if(length(n) == 1L) n <- rep(as.numeric(n), nobs)
    if(length(n) != nobs || any(!is.finite(n)) || any(n <= 0) || any(n != round(n)))
        stop(sprintf("For %s residuals, n must contain positive integers and have length 1 or the number of observations.", model_name),
             call. = FALSE)
    as.numeric(n)
}
.GeoResidual_pearson <- function(obs, fitted_mean, fitted_var) {
    if(any(!is.finite(fitted_mean)) || any(!is.finite(fitted_var)) ||
       any(fitted_var < 0))
        stop("Non-finite fitted marginal moments encountered while computing residuals.",
             call. = FALSE)
    den <- sqrt(pmax(fitted_var, .Machine$double.eps))
    (obs - fitted_mean) / den
}

if(model %in% c("Binary", "Bernoulli")) {
    prob <- pnorm(c(mu))
    res1 <- .GeoResidual_pearson(dd, prob, prob * (1 - prob))
}
if(model == "Binomial") {
    nn <- .GeoResidual_n(fit$n, length(dd), model)
    prob <- pnorm(c(mu))
    fitted_mean <- nn * prob
    res1 <- .GeoResidual_pearson(dd, fitted_mean, nn * prob * (1 - prob))
}
if(model == "BinomialLogistic") {
    nn <- .GeoResidual_n(fit$n, length(dd), model)
    prob <- plogis(c(mu))
    fitted_mean <- nn * prob
    res1 <- .GeoResidual_pearson(dd, fitted_mean, nn * prob * (1 - prob))
}
if(model %in% c("Geom", "Geometric", "BinomialNeg")) {
    nn <- .GeoResidual_n(if(model %in% c("Geom", "Geometric")) 1 else fit$n,
                         length(dd), model)
    prob <- pnorm(c(mu))
    fitted_mean <- nn * (1 - prob) / prob
    fitted_var <- nn * (1 - prob) / (prob^2)
    res1 <- .GeoResidual_pearson(dd, fitted_mean, fitted_var)
}
if(model == "Poisson") {
    lambda <- exp(c(mu))
    res1 <- .GeoResidual_pearson(dd, lambda, lambda)
}
if(model == "PoissonGamma") {
    fitted_mean <- exp(c(mu))
    shape <- as.numeric(param['shape'])
    if(length(shape) != 1L || !is.finite(shape) || shape <= 0)
        stop("PoissonGamma residuals require a positive finite shape parameter.",
             call. = FALSE)
    fitted_var <- fitted_mean * (1 + fitted_mean / shape)
    res1 <- .GeoResidual_pearson(dd, fitted_mean, fitted_var)
}
#########
#}

if(is.null(res1)) stop("residuals are not implemented for this model in GeoResiduals()\n")

fit$X <- matrix(1, nrow = length(dd), ncol = 1)

#### updating  object
mm=0
names(mm)="mean"
nuis_update=c(mm,nuisance[!sel])
fit$param=c(nuis_update,paramcorr)
fit$numbetas=1
fit$X <- matrix(1, nrow = length(dd), ncol = 1)

#if(!is.null(copula)){
#############################################
#if(copula=="Clayton"||copula=="Gaussian"){
if(model %in% c("Beta2")) {fit$param['shape']=2;fit$param['mean']=0; fit$param['max']=1; fit$param['min']=0}
if(model %in% c("Kumaraswamy2")) {fit$param['shape']=1;fit$param['mean']=0; fit$param['max']=1; fit$param['min']=0}
#}
#}else{
#####################################################
if(model %in% c("Gaussian","Logistic","Tukeyh","Tukeyh2","Tukeygh","SinhAsinh", "Gaussian_misp_StudentT","Gaussian_misp_Tukeygh",
         "StudentT","TwoPieceGauss","TwoPieceStudentT","TwoPieceGaussian","TwoPieceTukeyh","TwoPieceBimodal"))
{fit$param['sill']=1;fit$param['mean']=0}

if(model %in% c("SkewGaussian")) 
{
param['mean']=0;
 fit$param['skew']=as.numeric(param['skew'])/sqrt(as.numeric(param['sill']))
 fit$param['sill']=1
}


if(model %in% c("SkewLaplace"))

{
param['mean']=0;
fit$param['sill']=1
}




if(model %in% c("Gaussian_misp_SkewStudentT","SkewStudentT")) 
{param['mean']=0;
 fit$param['skew']=as.numeric(param['skew'])
 fit$param['df']= fit$param['df']
 fit$param['sill']=1
}
##
#}

## names of estimated and fixed parameters

obj=as.list(fit$param)
fit$param=obj[nm]
fit$fixed=obj[nf]

.drop_empty_na <- function(z) {
    if(length(z)==0) return(z)
    keep <- !vapply(z, function(v) length(v)==0L || all(is.na(v)), logical(1))
    z[keep]
}

fit$param=.drop_empty_na(fit$param)
fit$fixed=.drop_empty_na(fit$fixed)


###adding mean and variance if missing for some reason
if(model %in% c("Gaussian","SkewGaussian","SkewLaplace","Logistic","Tukeyh","Tukeyh2","Tukeygh","SinhAsinh",
 "Gaussian_misp_StudentT","Gaussian_misp_Tukeygh","StudentT","TwoPieceGauss",
 "TwoPieceStudentT","TwoPieceGaussian","TwoPieceTukeyh",
 "TwoPieceBimodal","Gaussian_misp_SkewStudentT","SkewStudentT")){
if(!sum(names(fit$param)=="mean")) fit$param["mean"]=0
if(!sum(names(fit$param)=="sill")) fit$param["sill"]=1

   sel=substr(names(unlist(fit$fixed)),1,4)=="mean"
            if(sum(sel)>=1) fit$fixed=fit$fixed[!sel]
}

if (model %in% c("Weibull", "Poisson", "PoissonGamma", "Binomial", "Gamma",  "LogGaussian", 
        "LogLogistic", "BinomialNeg", "Binary", "Bernoulli", "Geom", "Geometric", 
        "Gaussian_misp_Poisson", "PoissonZIP", "Gaussian_misp_PoissonZIP", 
        "BinomialNegZINB", "PoissonZIP1", "Gaussian_misp_PoissonZIP1", 
        "BinomialNegZINB1", "Beta2", "Kumaraswamy2", "Beta", 
        "Kumaraswamy")) {  
            if(!sum(names(fit$param)=="mean")) fit$param["mean"]=0              
            sel=substr(names(unlist(fit$fixed)),1,4)=="mean"
            if(sum(sel)>=1) fit$fixed=fit$fixed[!sel]
}


### formatting data
if(fit$spacetime)     ## spacetime
{if(!is.list(fit$coordx_dyn)) 
          data_res=matrix(res1,nrow=nrow(fit$data),ncol=ncol(fit$data),byrow=TRUE)
 else{ 
     ns=unlist(lapply(fit$coordx_dyn,nrow))
     sim_temp=list()
                    for(k in 1:length(fit$coordt))
                       { if(k==1) {indx=1:(sum(ns[1:k]))}
                         if(k>1)    {indx=(sum(ns[1:(k-1)])+1):(sum(ns[1:k]))}
                         sim_temp[[k]]=c(res1)[indx] }
    data_res=sim_temp  
 }
}   
else   {data_res=as.vector(res1)}   ## space
}
#####################################################################
#####################################################################
if(fit$bivariate)
{
 X=fit$X
 ns=fit$ns

namescorr <- CorrParam(fit$corrmodel) 
namesnuis <- NuisParam(fit$model,fit$bivariate,num_betas)
param <- c(fit$param, fit$fixed)
namesparam<- names(param)
paramcov <- param[namescorr]
nuisance <- param[namesnuis]
###
if(is.list(fit$coordx_dyn)) dd=unlist(fit$data)
else dd=c(t(fit$data))

dd1=dd[1:ns[1]]
dd2=dd[(ns[1]+1):(ns[2]+ns[1])]
###

###
 sel1=substr(names(nuisance),1,6)=="mean_1"
 beta1=as.numeric(nuisance[sel1])
 sel2=substr(names(nuisance),1,6)=="mean_2"
 beta2=as.numeric(nuisance[sel2])

if(!is.null(X)){
 X1=as.matrix(X[1:ns[1],]);
 X2=as.matrix(X[(ns[1]+1):(ns[2]+ns[1]),]); 
}
else {
     X1=matrix(1,nrow=ns[1],ncol=1)
     X2=matrix(1,nrow=ns[2],ncol=1)
}
mu1=X1%*%beta1
mu2=X2%*%beta2

if(model %in% c("Gaussian","Gauss")){
       res1=(dd1-c(mu1))/sqrt(as.numeric(paramcov['sill_1']))
       res2=(dd2-c(mu2))/sqrt(as.numeric(paramcov['sill_2']))
        }



if(!is.list(fit$coordx_dyn)) 
           {data_res=rbind(res1,res2)}
else {data_res=list(); data_res[[1]]=res1;data_res[[2]]=res2}
}



#####################################################################
#####################################################################

fit$data=data_res
## Preserve the complete fitted object.  Diagnostic metadata such as
## neighbourhood indices, thinning design, anisotropy, pair counts and
## convergence information must survive the residual transformation.
fit$param=as.list(fit$param)
fit$fixed=as.list(fit$fixed)
fit$X=NULL
fit$residual_type=if(model %in% c("Binary", "Bernoulli", "Binomial", "BinomialLogistic",
                                  "Geom", "Geometric", "BinomialNeg", "Poisson",
                                  "PoissonGamma"))
                    "Pearson" else "model-standardized"
fit$call=call
class(fit)="GeoFit"
fit
#########################
}
