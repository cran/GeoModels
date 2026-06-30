####################################################
### File name: GeoSimCopula.r
####################################################


# Simulate spatial and spatio-temporal random felds:
GeoSimCopula <- function(coordx, coordy=NULL, coordz=NULL,coordt=NULL, coordx_dyn=NULL,corrmodel, distance="Eucl", grid=FALSE,
     method="cholesky",model='Gaussian', n=1, param,anisopars=NULL, radius=1, sparse=FALSE,
     copula="Gaussian",X=NULL,spobj=NULL,nrep=1,progress=FALSE)
{

if(is.null(CkCorrModel (corrmodel))) stop("The name of the correlation model  is not correct\n")

if(is.null(CkModel(model))) stop("The name of the  model  is not correct\n")
    corrmodel=gsub("[[:blank:]]", "",corrmodel)
    model=gsub("[[:blank:]]", "",model)
    distance=gsub("[[:blank:]]", "",distance)
    method=gsub("[[:blank:]]", "",method)

##############################################################################
###### extracting sp object informations if necessary              ###########
##############################################################################
bivariate<-CheckBiv(CkCorrModel(corrmodel))
spacetime<-CheckST(CkCorrModel(corrmodel))
space=!spacetime&&!bivariate
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

if((copula!="Clayton")&&(copula!="Gaussian")&&(copula!="AMH")&&(copula!="SkewGaussian")) stop("the type of copula is wrong")

if(copula=="Clayton") {if(is.null(param$nu)) stop("Clayton copula need a nu parameter")}
if(copula=="AMH"    ) {if(is.null(param$nu)) stop("AMH copula need a nu parameter")}
if(copula=="SkewGaussian") {if(is.null(param$nu)) stop("Skew Gaussian copula need a nu parameter")}

#### corr parameters
paramcorr=param[CorrParam(corrmodel)]

    
SIM=list()
###################################################
### starting number of replicates
###################################################
for( L in 1:nrep){



####Gaussian copula #############################################
if(copula=="Gaussian")
{
param1=c(list(mean=0,sill=1,nugget=param$nugget),paramcorr)

sim=GeoSim(coordx=coordx, coordy=coordy,coordz=coordz,coordt=coordt, coordx_dyn=coordx_dyn,corrmodel=corrmodel, 
    distance=distance, grid=grid,
 method=method,model='Gaussian', n=1, param=param1,anisopars=anisopars, radius=radius, sparse=sparse,nrep=1)
unif=pnorm(sim$data,mean=0,sd=1);

}

####skewGaussian copula #############################################
if(copula=="SkewGaussian")
{

nu <- as.numeric(param$nu)
if(abs(nu) >= 1) stop("nu parameter must be between -1 and 1")
alpha <- nu / sqrt(1 - nu^2)
   
param1=c(list(mean=0,sill=1,nugget=param$nugget,skew=alpha),paramcorr)
sim=GeoSim(coordx=coordx, coordy=coordy,coordz=coordz,coordt=coordt, coordx_dyn=coordx_dyn,corrmodel=corrmodel, 
    distance=distance, grid=grid,
     method=method,model='SkewGaussian', n=1, param=param1,anisopars=anisopars, radius=radius, sparse=sparse,nrep=1)
omega=as.numeric(sqrt((alpha^2 + 1)/1))
unif=sn::psn(sim$data,xi=0,omega= omega,alpha= alpha)
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
pp=round(as.numeric(param['nu']))
param1=c(list(shape1=pp,shape2=2,sill=1,mean=0,min=0,max=1,nugget=param$nugget),paramcorr)
sim=GeoSim(coordx=coordx, coordy=coordy,coordz=coordz,coordt=coordt, coordx_dyn=coordx_dyn,corrmodel=corrmodel, 
    distance=distance, grid=grid,
     method=method,model='Beta', n=1, param=param1,anisopars=anisopars, radius=radius, sparse=sparse,nrep=1)
unif=(sim$data)^(pp/2)
}
####AMH copula #############################################
if(copula=="AMH")
{

param1=c(list(mean=0,nugget=as.numeric(param$nugget),shape=2),paramcorr)
sim=GeoSim(coordx=coordx, coordy=coordy,coordz=coordz,coordt=coordt, coordx_dyn=coordx_dyn,corrmodel=corrmodel, 
    distance=distance, grid=grid,
   method=method,model='Gamma', param=param1,anisopars=anisopars, radius=radius, sparse=sparse,nrep=1)
param2=c(list(mean=as.numeric(param['nu']),nugget=as.numeric(param$nugget)),paramcorr)
sim2=GeoSim(coordx=coordx, coordy=coordy,coordz=coordz,coordt=coordt, coordx_dyn=coordx_dyn,corrmodel=corrmodel, 
    distance=distance, grid=grid,
     method=method,model='BinomialNeg', n=1, param=param2,anisopars=anisopars, radius=radius, sparse=sparse,nrep=1)
kk=sim$data/(sim2$data+1)
pp=pnorm(as.numeric(param['nu']))
unif=pp/(exp(kk)+pp-1)
}
####################################################################

if(!sim$bivariate){
           if(is.null(dim(X))) {X=as.matrix(rep(1,sim$numcoord*sim$numtime))}
           mm <- rep(0, sim$numcoord*sim$numtime)
           sel=substr(names(param),1,4)=="mean";
           num_betas=sum(sel) 
           if(num_betas==1)  mm<-rep(as.numeric(param$mean), sim$numcoord*sim$numtime)
           if(num_betas>1)   mm<-as.numeric(X%*%as.numeric((param[sel])))

           # Do not reset param$mean, param$mean1, ... here.
           # This function may simulate several replications in the same call.
           # If param is modified inside the nrep loop, then from the second
           # replication onward the mean parameters are zero and the BinomialNeg
           # marginal becomes qnbinom(..., prob = pnorm(0)) = qnbinom(..., prob = 0.5).
    }       
##############################
##############################
#######     models     #######
##############################
if(!sim$bivariate) {}

####################################
############ discrete  RF ##########
####################################
if(model=="Binomial") {
 simcop=qbinom(unif, size=n, prob=pnorm(mm))
}
if(model=="BinomialNeg") {
 simcop=qnbinom(unif, size=n, prob=pnorm(mm))
}
if(model=="BinomialNegZINB") {
 simcop=qzinegbin(unif, size=n, #prob = pnorm(mm),                # require package VGAM
        munb=   n*(1-pnorm(mm))/pnorm(mm),#,n/pnorm(mm)-n,
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
simcop=exp(mm)*qgamma(unif,shape=p2/2,scale=p2/2)
         }
############
if(model=="Weibull") 
         {
p2=as.numeric(param$shape)
simcop= exp(mm)*qweibull(unif,shape=p2,scale=1/(gamma(1+1/p2 )))

         }
########################
if(model=="LogLogistic") 
         {
qllogis1 <- function(p, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE) {
  if (missing(shape)) stop("argument 'shape' is missing, with no default")
  if (shape <= 0) stop("'shape' must be positive")
  if (scale <= 0) stop("'scale' must be positive")
  if (log.p) { p <- exp(p)}
  if (any(p < 0 | p > 1)) stop("probabilities must be between 0 and 1")
  if (!lower.tail) {p <- 1 - p }
  result <- numeric(length(p)); result[p == 0] <- 0;result[p == 1] <- Inf
  valid <- p > 0 & p < 1
  if (any(valid)) {p_valid <- p[valid];result[valid] <- scale * (p_valid / (1 - p_valid))^(1/shape)}
  return(result)
}

p2=as.numeric(param$shape)
cc=gamma(1+1/p2)*gamma(1-1/p2)
simcop=qllogis1(unif,shape = p2,scale=exp(mm)/cc)            
         }
##############################

if(model=="LogGaussian")
{ 
    vv=as.numeric(param$sill)
   simcop = qlnorm(unif, mm-vv/2, sqrt(vv))
     #simcop = qlnorm(unif, 1+mm, sqrt(vv))
}
##############################
##############################
############ Real RF #########
##############################


if(model=="Gaussian") {
         simcop=qnorm(unif,mean=mm,sd=sqrt(as.numeric(param$sill)))
         }
if(model=="Logistic") 
         {
         simcop=qlogis(unif,location=mm,scale=sqrt(as.numeric(param$sill)))
         }
if(model=="StudentT") {
vv=as.numeric(param$sill)
dd=as.numeric(param$df)
simcop=mm+sqrt(vv)*qt(unif,df=round(1/dd))

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
  sk=as.numeric(param$skew)
  dd <- round(1 / as.numeric(param$df))
 simcop=mm+sqrt(vv)*sn::qst(unif, xi=0, omega=1, alpha=sk, nu=dd)
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
qtpTukeyh22= function(x,tail1,tail2){
  uu <- qnorm(x,0,1)
  res <- numeric(length(x))
  sel1 <- uu >= 0
  sel2 <- uu < 0
  res[sel1] <- uu[sel1]*exp(0.5*tail1*uu[sel1]^2)
  res[sel2] <- uu[sel2]*exp(0.5*tail2*uu[sel2]^2)
  return(res)
}
tail1 =as.numeric(param$tail1);tail2 = as.numeric(param$tail2)
vv=as.numeric(param$sill) 
simcop =mm+sqrt(vv)*qtpTukeyh22(unif,tail1,tail2)
}

#######################################  OK
if(model=="SkewLaplace")
{  

# Quantile della Skew/Asymmetric Laplace
qtpSkewLaplace22 <- function(u, sk) {
  mm=0
  vv=1
  s  <- sqrt(vv)
  res <- rep(NA_real_, length(u))
  hi <- u >= sk
  lo <- u <  sk
  res[hi] <- mm - (s / sk) * (log1p(-u[hi]) - log1p(-sk))
  res[lo] <- mm + (s / (1 - sk)) * (log(u[lo]) - log(sk))
  res
}
sk=as.numeric(param$skew) 
vv=as.numeric(param$sill)
simcop =mm+sqrt(vv)*qtpSkewLaplace22(unif,sk)
}
#######################################  OK
if(model=="Tukeygh")
{
tail = as.numeric(param$tail);skew =  as.numeric(param$skew) 
vv=as.numeric(param$sill) 
uu=qnorm(unif)
simcop=mm+sqrt(vv)*(exp(skew*uu)-1)*exp(0.5*tail*uu^2)/skew
}
#######################################   OK
if(model=="TwoPieceGaussian")
{
qtpGaussian = function(x,skew){
  res <- numeric(length(x))
  sel1 <- x > 0 & x < 0.5*(1+skew)
  sel2 <- x <= 1 & x >= 0.5*(1+skew)
  res[sel1] <- (1+skew)*qnorm(x[sel1]/(1+skew))
  res[sel2] <- (1-skew)*qnorm((x[sel2]-skew)/(1-skew))
  return(res)
}
skew = as.numeric(param$skew) 
vv=as.numeric(param$sill) 
simcop =mm+sqrt(vv)*qtpGaussian(unif,skew)
}

#######################################  
if(model=="TwoPieceStudentT")
{
qtpt = function(x,skew,df){
  res <- numeric(length(x))
  sel1 <- x > 0 & x < 0.5*(1+skew)
  sel2 <- x <= 1 & x >= 0.5*(1+skew)
  res[sel1] <- (1+skew)*qt(x[sel1]/(1+skew),df=df)
  res[sel2] <- (1-skew)*qt((x[sel2]-skew)/(1-skew),df=df)
  return(res)
}
vv=as.numeric(param$sill) 
skew = as.numeric(param$skew)
df   = 1/as.numeric(param$df)
simcop=mm+sqrt(vv)*qtpt(unif,skew,df)
}
#######################################   OK
if(model=="TwoPieceTukeyh")
{
qtukh=function(xx,tail)
{
  uu=qnorm(xx)
  q_t=uu*exp(0.5*tail*uu^2)
return(q_t)
}
qtptukey = function(x,skew,tail){
  res <- numeric(length(x))
  sel1 <- x > 0 & x < 0.5*(1+skew)
  sel2 <- x <= 1 & x >= 0.5*(1+skew)
  res[sel1] <- (1+skew)*qtukh(xx=x[sel1]/(1+skew),tail=tail)
  res[sel2] <- (1-skew)*qtukh(xx=(x[sel2]-skew)/(1-skew),tail=tail)
  return(res)
}
vv=as.numeric(param$sill) 
skew = as.numeric(param$skew)
tail= as.numeric(param$tail)
simcop =mm+sqrt(vv)*qtptukey(unif,skew,tail)
}




####################################
############ bounded  real RF    ###
####################################
if(model=="Kumaraswamy") 
{
p1=param$shape1;p2=param$shape2
pmin=as.numeric(param$min);pmax=as.numeric(param$max);
simcop=pmin + (pmax-pmin)*((1-unif^(1/p1))^(1/p2))
}
############
if(model=="Kumaraswamy2") 
{ # parametrization using beta median  regression
mm=1/(1+exp(-mm))
p2=as.numeric(param$shape)
pmin=as.numeric(param$min);pmax=as.numeric(param$max);
aa=log(1-mm^(p2))/log(0.5)
simcop=pmin + (pmax-pmin)*((1-unif^(aa))^(1/p2))
}
############
if(model=="Beta") 
{ # parametrization using beta regression
pmin=as.numeric(param$min);pmax=as.numeric(param$max);
simcop=pmin + (pmax-pmin)*qbeta(unif,shape1=param$shape1,shape2=param$shape2)
}
############
if(model=="Beta2") 
{ # parametrization using beta mean  regression

mm=1/(1+exp(-mm))
p2=as.numeric(param$shape)
pmin=as.numeric(param$min);pmax=as.numeric(param$max);
simcop=pmin + (pmax-pmin)*qbeta(unif,shape1=mm*p2,shape2=(1-mm)*p2)
}





############
if (!grid) simcop=c(simcop)
##############################




SIM[[L]]=simcop
}
#### end numrep ##############
##############################
##############################

 if(nrep==1) SIM=SIM[[1]] 

    GeoSim_Copula <- list(bivariate = sim$bivariate,
    coordx = sim$coordx,
    coordy = sim$coordy,
    coordt = sim$coordt,
    coordz=coordz,
    coordx_dyn =sim$coordx_dyn,
    corrmodel = corrmodel,
    data = SIM,
    distance = sim$distance,
    grid = sim$grid,
    model = model,
    method=method,
    n=n,
    numcoord = sim$numcoord,
    numtime = sim$numtime,
    param = param,
    radius = radius,
    spacetime = sim$spacetime,
    sparse=sim$sparse,
    copula=copula,
    nrep=nrep,
    X=X)
#}
##############################################
  structure(c(GeoSim_Copula, call = call), class = c("GeoSim_Copula"))
}
