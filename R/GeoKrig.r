####################################################
### File name: GeoKrig.r
####################################################


GeoKrig= function(estobj=NULL,data, coordx, coordy=NULL,coordz=NULL, coordt=NULL, coordx_dyn=NULL, corrmodel, distance="Eucl", grid=FALSE, loc, maxdist=NULL,
               maxtime=NULL, method="cholesky", model="Gaussian", n=1,nloc=NULL, mse=FALSE, lin_opt=TRUE, param, anisopars=NULL,
               radius=1, sparse=FALSE,taper=NULL, tapsep=NULL, time=NULL, type="Standard",type_mse=NULL, type_krig="Simple",weigthed=TRUE,
               which=1, copula=NULL,X=NULL,Xloc=NULL,Mloc=NULL,spobj=NULL,spdata=NULL)

{
##################################################################################################################
##################################################################################################################
##################################################################################################################
getInv <- function(covmatrix, b, mse = TRUE) {
  if (!covmatrix$sparse) {
    U <- MatDecomp(covmatrix$covmatrix, "cholesky")
    if (is.logical(U)) stop("Covariance matrix is not positive definite")
    vec <- forwardsolve(U, b)
    Invc <- forwardsolve(U, vec, transpose = TRUE)
    mse_val <- if (mse) crossprod(vec) else NULL
  } else {
    cc <- if (spam::is.spam(covmatrix$covmatrix)) covmatrix$covmatrix else spam::as.spam(covmatrix$covmatrix)
    U <- try(spam::chol.spam(cc), silent = TRUE)
    if (inherits(U, "try-error")) stop("Covariance matrix is not positive definite")
    vec <- spam::forwardsolve(U, b)
    Invc <- spam::backsolve(U, vec)
    mse_val <- if (mse) spam::crossprod.spam(vec) else NULL
  }
  list(a = Invc, b = mse_val)
}

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

   if(length(estobj$coordt)==1) coordt=NULL
   else coordt=estobj$coordt
   
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
   if(!is.null(anisopars))  {param$angle=NULL; param$ratio=NULL}
    
   X=estobj$X
}
##################################
###### end check geofit object####
##################################


if(is.null(CkModel(model))) stop("The name of the  model  is not correct\n")
    if( !is.character(corrmodel)|| is.null(CkCorrModel(corrmodel)))       stop("the name of the correlation model is wrong")
corrmodel=gsub("[[:blank:]]", "",corrmodel); model=gsub("[[:blank:]]", "",model)
distance=gsub("[[:blank:]]", "",distance); method=gsub("[[:blank:]]", "",method)
type_krig=gsub("[[:blank:]]", "",type_krig); type=gsub("[[:blank:]]", "",type)

########################################################
####### extracting sp objects if necessary 
########################################################
bivariate<-CheckBiv(CkCorrModel(corrmodel))
spacetime<-CheckST(CkCorrModel(corrmodel))
space=!spacetime&&!bivariate
if(!is.null(spobj)) {
   if(space||bivariate){
        a=sp2Geo(spobj,spdata); coordx=a$coords 
       if(!a$pj) {if(distance!="Chor") distance="Geod"}
    }
   if(spacetime){
        a=sp2Geo(spobj,spdata); coordx=a$coords ; coordt=a$coordt 
        if(!a$pj) {if(distance!="Chor") distance="Geod"}
     }
   if(!is.null(a$Y)&&!is.null(a$X)) {data=a$Y ; X=a$X }
}
######################################
############ some checks##############
######################################
#if(is.vector(loc))    loc=t(as.matrix(loc)) ## case of 1 location sites given as vector
#if(!is.matrix(loc))   loc=as.matrix(loc)
if (is.vector(loc)) {
  loc <- matrix(loc, nrow = 1)
} else if (!is.matrix(loc)) {
  loc <- as.matrix(loc)
}

if(!(ncol(loc)==2||ncol(loc)==3))   stop("loc parameter must be a matrix  N X 2 or N X 3 ")

#if(!is.null(Xloc))
#         { if(is.vector(Xloc)) Xloc=matrix(Xloc,nrow=1)
#           else                Xloc=as.matrix(Xloc)
#         }
#if(!is.null(X)) X=as.matrix(X)
if (!is.null(Xloc)) {
  if (is.vector(Xloc)) {
    Xloc <- matrix(Xloc, nrow = 1)
  } else if (!is.matrix(Xloc)) {
    Xloc <- as.matrix(Xloc)
  }
}
if (!is.null(X) && !is.matrix(X)) X <- as.matrix(X)

if(is.matrix(X) &&is.null(Xloc))  stop("Covariates for locations to predict are missing \n")
if(is.null(X) &&is.matrix(Xloc))  stop("Covariates  are missing \n")
if(CheckST(CkCorrModel(corrmodel))) {if(is.null(time)) stop("At least one temporal instants is needed for space-time kriging\n ")  } 
if(!is.null(Mloc)&&  !is.null(Xloc)) stop("Mloc or Xloc must be fixed\n")    
if((length(param$mean)>1)&&is.null(Mloc)) stop("Mloc must be fixed \n")  
loc_orig=loc
if(!is.null(anisopars)) {  loc=GeoAniso(loc,c(anisopars$angle,anisopars$ratio))}

if(!is.null(Mloc))
{ 
    if(!is.vector(Mloc)) stop("Mloc must be a vector")
    if(space){ 
       if(length(Mloc)==1) Mloc=rep(Mloc,nrow(loc));if(nrow(loc)!=length(Mloc)) stop("Lenght of the  mean vector fixed does not match the number of locations to predict \n")
         }
    if(spacetime){
       if(length(Mloc)==1) Mloc=rep(Mloc,nrow(loc)*length(time));if(nrow(loc)*length(time)!=length(Mloc)) stop("Lenght of the  mean vector fixed does not match the number of locations to predict \n")
          }
}
if(!(type_krig=="Simple"||type_krig=="Optim")) stop("type_krig must be equal to Simple or Optim")
######################################
############ end some checks##########
######################################



#### checking 2-d or 3-d compatibility ##########
if(!is.null(estobj)){
if(is.null(estobj$coordz)&&ncol(loc)==3) stop("locations to predict are three dimensional but location used to predict are two dimensional \n")
if(!is.null(estobj$coordz)&&ncol(loc)==2) stop("locations to predict are two dimensional but location used to predict are three dimensional \n")
}
else{ if(is.null(coordx_dyn))
       { 
    if(is.null(coordz)){
       if(ncol(coordx)==2&&ncol(loc)==3) stop("locations to predict are three dimensional but location used to predict are two dimensional \n")
       if(ncol(coordx)==3&&ncol(loc)==2) stop("locations to predict are two dimensional but location used to predict are three dimensional \n")
                        }
       else {
       if(ncol(loc)==2) stop("locations to predict are two dimensional but location used to predict are three dimensional \n")
           }
       }
      else {} 
}
######################################
######################################
#### locations points to predict ###
if(is.null(time)) time=0; numloc = nrow(loc); tloc = length(time);
if(!tloc)  tloc = 1; 
if(!is.null(estobj)){
if(is.null(estobj$coordz)) {locx = loc[,1];locy = loc[,2]; locz=double(length(loc[,1]))}
else                {locx = loc[,1];locy = loc[,2]; locz=loc[,3]}
}
else{
if(ncol(loc) ==2) {locx = loc[,1];locy = loc[,2]; locz=double(length(loc[,1]))}
if(ncol(loc) ==3) {locx = loc[,1];locy = loc[,2]; locz=loc[,3]}

}

bb=0

######################################
############ standard kriging ########
######################################
if(type %in% c("Standard","standard")) {

################################################
#if(!bivariate) {
######################################################## 
############ optimal prediction cases ##################   
######################################################## 
logGausstemp=SinhAsinhtemp=Tukeyh2temp=Tukeyhtemp=FALSE 
if((model %in% c("LogGaussian","SinhAsinh","Tukeyh2","Tukeyh"))&&type_krig=="Optim")
{
 if(!bivariate) 
  { ## saving means and variance ###
    uparam=unlist(param); sel=substr(names(uparam),1,4)=="mean"; gbetas=  as.numeric(uparam[sel])

    if(!is.null(X)) {me=X%*%gbetas; meloc=Xloc%*%gbetas}
    else { if(length(gbetas)>1) {me=param$mean;meloc=Mloc} ## external mean
           else{me=meloc=param$mean}}
    vvm=as.numeric(param['sill']); 
    ###############
    if(model=="LogGaussian") {logGausstemp=TRUE}
    if(model=="Tukeyh")      {Tukeyhtemp=TRUE;th=as.numeric(param['tail']); param['tail']=NULL;}
    if(model=="Tukeyh2")     {Tukeyh2temp=TRUE;t1=as.numeric(param['tail1']); t2=as.numeric(param['tail2']);param['tail1']=NULL; param['tail2']=NULL;}
    if(model=="SinhAsinh")  {SinhAsinhtemp=TRUE; sk=as.numeric(param['skew']); tail=as.numeric(param['tail']); param['skew']=NULL; param['tail']=NULL;}
   ## setting standard gaussian
    model="Gaussian" ; param['sill']=1; param['mean']=0   # we need a standard "Gaussian" covariance matrix
  }

    if(bivariate){}
}
######################################################## 
######################################################## 
Mtemp=NULL
   if(model %in% c("Weibull","Gamma","LogLogistic","LogGaussian")&&type_krig=="Simple")          # we need a x covariane matrix with mean=0   x=gamma,weibull,loglogistic
{
     paramtemp=param; sel=substr(names(param),1,4)=="mean";  ## selecting mean values
     meantemp=names(param[sel])             ## saving   mean parameters
     param=param[!sel];                     ## not mean paramters
     if(length(paramtemp$mean)>1)    Mtemp=paramtemp$mean
     param$mean=0;        ## mean must be =0 when calling covariance matrix
     Xtemp=X;X=NULL                         ## saving X and setting X=NULL
}
#####################################################
##### computing covariance matrix  ##################
#####################################################




covmatrix = GeoCovmatrix(coordx=coordx, coordy=coordy,coordz=coordz, coordt=coordt, coordx_dyn=coordx_dyn,
            corrmodel=corrmodel, distance= distance,grid=grid,maxdist= maxdist,maxtime=maxtime,model=model,n=n,
            param=param, anisopars=anisopars, radius=radius,sparse=sparse,taper=taper,tapsep=tapsep,type=type,copula=copula,X=X)
    
#####################################################
covmatrix$param=unlist(covmatrix$param)
if(bivariate) tloc=1
spacetime_dyn=FALSE; if(!is.null(covmatrix$coordx_dyn)) spacetime_dyn=TRUE
########dimat is dimension  ######
if(!spacetime_dyn) dimat=covmatrix$numcoord*covmatrix$numtime
if(spacetime_dyn)  dimat =sum(covmatrix$ns)
dimat2=numloc*tloc
if(is.null(X)) XX=rep(1,dimat)
else XX=X

MM=NULL
if(!is.null(Mtemp))  {MM=Mtemp;param$mean=0}
else  { if(length(param$mean)>1)  { MM=param$mean; param$mean=0 } } ## in the case of non constant  external mean        
###############
if(model %in% c("Weibull","Gamma","LogLogistic","LogGaussian")&&type_krig=="Simple") {
          if(is.null(Xtemp)) X=matrix(rep(1,dimat))
          else               X=Xtemp
          param=paramtemp
          if(!is.null(Mtemp)) param$mean=0
          covmatrix$namesnuis=unique(c(meantemp,covmatrix$namesnuis))   
}
else { if(!is.null(X)) X=covmatrix$X 
           else        X=matrix(1,nrow=dimat,ncol=1)
          
     }

###############
num_betas=ncol(X); NS=0
if(spacetime||bivariate) { NS=cumsum(covmatrix$ns); NS=c(0,NS)[-(length(covmatrix$ns)+1)] }
if(is.null(Xloc)) Xloc=as.matrix(rep(1,dimat2))
else {if(spacetime_dyn) Xloc=as.matrix(Xloc)}
###########
nuisance = param[covmatrix$namesnuis]; nuisance=Filter(Negate(is.null),nuisance)
sel=substr(names(nuisance),1,4)=="mean"
betas=as.numeric(nuisance[sel])   ## mean paramteres
if(length(betas)>1 && is.null(X)) stop("Covariates matrix X is missing\n")
if(bivariate) { sel1=substr(names(nuisance),1,6)=="mean_1"; betas1=as.numeric(nuisance[sel1])   ## mean1 paramteres
                    sel2=substr(names(nuisance),1,6)=="mean_2"; betas2=as.numeric(nuisance[sel2])   ## mean1 paramteres
              }
other_nuis=as.numeric(nuisance[!sel])
 
############
cmodel=corrmodel
cdistance=distance
corrmodel = CkCorrModel(covmatrix$corrmodel)
distance = CheckDistance(covmatrix$distance)
corrparam = covmatrix$param[covmatrix$namescorr]# selecting the correlation parametrs
if(bivariate) {if(!(which==1 || which==2)) {stop("which  parameter must be 1 or 2")}}
pred = NULL
varpred=varpred2=vv=vv2=NULL
k = 0
  
###### spatial case not grid
if(!grid){
  if(!spacetime&&!bivariate){
    if(is.null(covmatrix$coordz))  ccc=cbind(covmatrix$coordx,covmatrix$coordy,0) 
    else                           ccc=cbind(covmatrix$coordx,covmatrix$coordy,covmatrix$coordz)
    }
}
###### spatial anisotropy case 
if(!is.null(anisopars)) {  ccc=GeoAniso(ccc,c(anisopars$angle,anisopars$ratio))}

if(grid) {if(is.null(coordz)) { ccc=as.matrix(expand.grid(covmatrix$coordx,covmatrix$coordy));ccc=cbind(ccc,0) }
          else                { ccc=as.matrix(expand.grid(covmatrix$coordx,covmatrix$coordy,covmatrix$coordz)) }
          grid=FALSE
         }
else  {
      if((spacetime||bivariate)&&(!spacetime_dyn)) 
      {  ## space time not dynamic
          if(!is.null(covmatrix$coordz)) ccc=cbind(rep(covmatrix$coordx,covmatrix$numtime),rep(covmatrix$coordy,covmatrix$numtime),rep(covmatrix$coordz,covmatrix$numtime))
          else                          ccc=cbind(rep(covmatrix$coordx,covmatrix$numtime),rep(covmatrix$coordy,covmatrix$numtime),0)
      }
      if((spacetime||bivariate)&&(spacetime_dyn)) {
            ccc=do.call(rbind,args=c(covmatrix$coordx_dyn));
            if(ncol(ccc)==2) ccc=cbind(ccc,0)
                                                  }
}
###############################################################

if((spacetime||bivariate)&&spacetime_dyn)  dataT=t(unlist(data))  
else dataT=t(data)
  
####################################################################
############### computing  correlation vector ######################
####################################################################
if(covmatrix$model %in% c(1,10,18,21,12,26,24,27,38,29,39,28,9, 34,40,20,22))    
{
## gaussian=1
## skew gaussian=10
## student =12
## skewstudent =18
## gamma=21
## weibull=26
## loglogistic=24
## loggaussian=22 ojo
## twopieceStudentT=27
## beta=28
## twopieceGaussian=29
## twopieceTukeyh=38
## twopiecebimodal=39
## tukeygh=9
## tukey=34
## tukeyyh2=40
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
if(bivariate){  corri=cc$corri  }  #  adding models.....

else {  ## for each model..
        rho=cc$corri; cc=(1-as.numeric(covmatrix$param["nugget"]))*rho
############## fin qua ok

################ Gaussian ##################################
if(covmatrix$model==1)   {  vv=as.numeric(covmatrix$param['sill']); 
          if(is.null(covmatrix$copula)) corri=cc 
          else {  if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov_fast(cc,covmatrix$model, nuisance)}
} 
#################### skew gaussian #########################################
if(covmatrix$model==10) {    
vv=as.numeric(covmatrix$param['sill'])
if(is.null(covmatrix$copula)){
    corr2=rho^2; sk=as.numeric(covmatrix$param['skew']);sk2=sk^2; 
    corri=(2*sk2)*(sqrt(1-corr2) + rho*asin(rho)-1)/(pi*vv+sk2*(pi-2)) + (cc*vv)/(vv+sk2*(1-2/pi))
                        }
    else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov_fast(cc,covmatrix$model, nuisance)}
}
##################### student T ########################################
if(covmatrix$model==12) 
{
vv=as.numeric(covmatrix$param['sill']); nu=1/as.numeric(covmatrix$param['df'])
if(is.null(covmatrix$copula)){corri=((nu-2)*gamma((nu-1)/2)^2*Re(hypergeo::hypergeo(0.5,0.5 ,nu/2 ,rho^2))*cc)/(2*gamma(nu/2)^2)}
else {
      if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov_fast(cc,covmatrix$model, nuisance)
     }
}
########################## tukeyh ###################################
if(covmatrix$model==34&&type_krig=="Simple") 
{   vv=as.numeric(covmatrix$param['sill']); h=as.numeric(covmatrix$param['tail'])
    if(h>0){ 
            if(is.null(covmatrix$copula)){ corri=(cc*(1-2*h)^(1.5))/((1-h)^2-(h*cc)^2)^(1.5) }
            else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov_fast(cc,covmatrix$model, nuisance)}
            }
    else{ if(is.null(covmatrix$copula))  corri=cc
          else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov_fast(cc,covmatrix$model, nuisance)}
        }
}
#########################tukeyh2###################################
if(covmatrix$model==40&&type_krig=="Simple") # 
{ vv=as.numeric(covmatrix$param['sill']);
  tail1=as.numeric(covmatrix$param['tail1']);tail2=as.numeric(covmatrix$param['tail2'])
  if(is.null(covmatrix$copula)){ 
                           hr=tail1;hl=tail2
                             x1=1-(1-cc^2)*hr; y1=1-(1-cc^2)*hl; x2=(1-hr)^2-(cc*hr)^2;y2=(1-hl)^2-(cc*hl)^2
                             g=1-hl-hr+(1-cc^2)*hl*hr
                             h1=sqrt(1-cc^2/(x1^2))+(cc/x1)*asin(cc/x1);h2=sqrt(1-cc^2/(y1^2))+(cc/y1)*asin(cc/y1)
                             h3=sqrt(1-cc^2/(x1*y1))+sqrt(cc^2/(x1*y1))*asin(sqrt(cc^2/(x1*y1)))
                             p1=x1*h1/(2*pi*(x2)^(3/2))+cc/(4*(x2)^(3/2))
                             p2=y1*h2/(2*pi*(y2)^(3/2))+cc/(4*(y2)^(3/2))
                             p3=-(x1*y1)^(1/2)*h3/(2*pi*(g)^(3/2))+cc/(4*(g)^(3/2))
                             mm=(hr-hl)/(sqrt(2*pi)*(1-hl)*(1-hr))
                             vv1=0.5*(1-2*hl)^(-3/2)+0.5*(1-2*hr)^(-3/2)-(mm)^2
                             corri=(p1+p2+2*p3-mm^2)/vv1  # correlation
                           }
    else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov_fast(cc,covmatrix$model, nuisance)}
}
####################### sas ######################################
if (covmatrix$model == 20 && type_krig == "Simple") { 
  # Estrazione parametri una volta sola
  vv <- as.numeric(covmatrix$param['sill'])
  tail <- as.numeric(covmatrix$param['tail'])
  skew <- as.numeric(covmatrix$param['skew'])
    if(is.null(covmatrix$copula)){ 
  d <- tail
  e <- skew
  mm=sinh(e/d)*exp(0.25)*(besselK(.25,(d+1)/(2*d))+besselK(.25,(1-d)/(2*d)))/(sqrt(8*pi))
  vv1=cosh(2*e/d)*exp(0.25)*(besselK(.25,(d+2)/(2*d))+besselK(0.25,(2-d)/(2*d)))/(sqrt(32*pi))-0.5-mm^2
  corri=corrsas(cc,e,d)
}
      else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov_fast(cc,covmatrix$model, nuisance)}
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
                else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov_fast(cc,covmatrix$model, nuisance)}
}
######################## # two piece StudenT#####################################
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
          else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov_fast(cc,covmatrix$model, nuisance)}
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
                     else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov_fast(cc,covmatrix$model, nuisance)}
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
                      else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov_fast(cc,covmatrix$model, nuisance)}
}
############################################################
       if(covmatrix$model==28)   ##  beta case
{
        
         shape1=as.numeric(covmatrix$param['shape1']);shape2=as.numeric(covmatrix$param['shape2']);

             if(is.null(covmatrix$copula)){ 
          corr2=cc^2       
         cc1=0.5*(shape1+shape2); vv=shape1*shape2/((cc1+1)*(shape1+shape2)^2)
         idx=which(abs(corr2)>1e-10);corr22=corr2[idx]; nu2=shape1/2;alpha2=shape2/2
         res=0;ss=0;k=0
         while(k<=100){
             p1=2*(lgamma(cc1+k)-lgamma(cc1)+lgamma(nu2+1+k)-lgamma(nu2+1));p2=lgamma(k+1)+(lgamma(nu2+k)-lgamma(nu2))+2*(lgamma(cc1+1+k)-lgamma(cc1+1))
             b1=p1-p2; b2=log(hypergeo::genhypergeo(U=c(cc1+k,cc1+k,alpha2), L=c(cc1+k+1,cc1+k+1), polynomial=TRUE,maxiter=1000, z=corr22)); b3=k*log(corr22)
             sum=exp(b1+b2+b3); res=res+sum
             if (all(sum<1e-6)){ break} else{ A=res}
         k=k+1
        }
         cc[idx]=A; corri=shape1*(cc1 + 1 ) * ((1-corr2)^(cc1) *cc -1)/shape2 ## correlation
         corri[-idx]=0
    }
       else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov_fast(cc,covmatrix$model, nuisance)}
}
############################################################    


if(covmatrix$model==21) 
{ # gamma
                        sh=as.numeric(covmatrix$param["shape"])
                         if(is.null(covmatrix$copula)){  corri=cc^2 }
                         else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov_fast(cc,covmatrix$model, nuisance)}
}
if(covmatrix$model==26) 
{  # weibull
                        sh=as.numeric(covmatrix$param['shape']);
                        if(is.null(covmatrix$copula)){ 
                         bcorr=    (gamma(1+1/sh))^2/((gamma(1+2/sh))-(gamma(1+1/sh))^2)
                        corri=bcorr*((1-cc^2)^(1+2/sh)*Re(hypergeo::hypergeo(1+1/sh,1+1/sh ,1 ,cc^2)) -1)
                         }
                        else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov_fast(cc,covmatrix$model, nuisance)}
                        
}
if(covmatrix$model==24) 
{  # loglogistic
                        sh=as.numeric(covmatrix$param['shape'])
                             if(is.null(covmatrix$copula)){ 
                        corri=((pi*sin(2*pi/sh))/(2*sh*(sin(pi/sh))^2-pi*sin(2*pi/sh)))*(Re(hypergeo::hypergeo(-1/sh, -1/sh, 1,cc^2))*Re(hypergeo::hypergeo(1/sh, 1/sh, 1,cc^2)) -1)
                        }
                        else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov_fast(cc,covmatrix$model, nuisance)}
}
if(covmatrix$model==22&&type_krig=="Simple")
 {  # loggauusian
                        ss=as.numeric(covmatrix$param['sill']);
                               if(is.null(covmatrix$copula)){  corri= (exp(ss*cc)-1)/(exp(ss)-1)}
                                     else { if(covmatrix$copula=="Gaussian") corri= gaussian_copula_cov_fast(cc,covmatrix$model, nuisance)}
}
############################################################
 }        
############################################################
if(bivariate)  {          if(which==1)    vvar=covmatrix$param["sill_1"]+covmatrix$param["nugget_1"]
                          if(which==2)    vvar=covmatrix$param["sill_2"]+covmatrix$param["nugget_2"]
               }
else    {
    #####  computing mean and variances for each model
     if(covmatrix$model==1)   {vvar= vv;                M=0 }#gaussian
     if(covmatrix$model==10)  {vvar= vv+sk^2*(1-2/pi) ; M=sk*sqrt(2/pi)}## skewgaus
     if(covmatrix$model==12)  {vvar= vv*nu/(nu-2) ;     M=0 }     ## studentT
    #  if(covmatrix$model==9)  {                                      #tukeygh
    #                           tail2=tail*tail;skew2=skew*skew; u=1-tail;
    #                           mm=(exp(skew2/(2*u))-1)/(skew*sqrt(u));
    #                           vvar=vv* ((exp(2*skew2/(1-2*tail))-2*exp(skew2/(2*(1-2*tail)))+1)/(skew2*sqrt(1-2*tail))-mm^2)
    #                           M=sqrt(vv)*mm
    #                           }
     if(covmatrix$model==34&&type_krig=="Simple") {vvar= vv*(1-2*h)^(-1.5); M=0 }          ## tukey h
     if(covmatrix$model==40&&type_krig=="Simple") { vvar= vv* vv1 ; M=sqrt(vv)*mm}        ## tukey h2        
     if(covmatrix$model==20&&type_krig=="Simple") { vvar= vv* vv1; M=sqrt(vv)*mm}                        ## sas    
   ##############################################################
     if(covmatrix$model==18)  { #skew student T
                               D1=(nu-1)*0.5; D2=nu*0.5;
                               mm=sqrt(nu)*gamma(D1)*sk/(sqrt(pi)*gamma(D2));
                               vvar=vv*(nu/(nu-2)-mm^2); M=sqrt(vv)*mm
                              }
     if(covmatrix$model==27)  { # two piece studentT
                              ttemp=gamma(0.5*(nu-1))/gamma(0.5*nu)
                              vvar= vv*((nu/(nu-2))*(1+3*sk2) - 4*sk2*(nu/pi)*ttemp^2); M= - sqrt(vv)*2*sk*sqrt(nu/pi)*ttemp
                              }
     if(covmatrix$model==28)  { #beta
                              ssup=as.numeric((covmatrix$param["max"]-covmatrix$param["min"]))
                              vvar=(ssup^2)*(shape1*shape2/((cc1+1)*(shape1+shape2)^2))
                              M=ssup*shape1/(shape1+shape2)+as.numeric(covmatrix$param["min"])
                              muloc=0;mu=0
                              }
     if(covmatrix$model==29)  {vvar= vv*((1+3*sk2) - 8*sk2/pi ) ;M=-sqrt(vv)*2*sk*sqrt(2/pi)} # two piece Gaussian
     if(covmatrix$model==38)  {vvar= vv*(ff - mm);M= -sqrt(vv)*2*sk*sqrt(2/pi)/(1-tail)} # two piecetukeyh
     if(covmatrix$model==39)  {vvar= vv*vari/(nn^(2/alpha)*gamma(nu/2)^2);M=-sqrt(vv)*sk*2^(1/alpha+1)*gamma(nu/2+1/alpha)/(nn^(1/alpha)*gamma(nu*0.5)) }#twopiecebimodal
     if(covmatrix$model==21)  {vvar= 2/sh }          #gamma
     if(covmatrix$model==24)  {vvar= 2*sh*sin(pi/sh)^2/(pi*sin(2*pi/sh))-1 }  ##loglogistic
     if(covmatrix$model==26)  {vvar= gamma(1+2/sh)/gamma(1+1/sh)^2-1   }    ## weibull
     if(covmatrix$model==22&&type_krig=="Simple")  {vvar= exp(ss)-1   }    ## loggaussian

     if(SinhAsinhtemp||Tukeyh2temp||Tukeyhtemp||logGausstemp) vvar=1
     
 }
###############################################################  

################################################################
############setting mean for the exte mean case ################
################################################################
   if(!bivariate){
    if(is.null(MM)) {mu=X%*%betas; muloc=Xloc%*%betas}
                    else            {mu=MM;         muloc=Mloc}  
                        }          # for non constant external mean                
    else{ X11=X[1:covmatrix$ns[1],]; X22=X[(covmatrix$ns[1]+1):(covmatrix$ns[1]+covmatrix$ns[2]),]
               if(!is.null(Xloc)) 
                    { X11_loc=Xloc[(1:(nrow(Xloc)/2)),]; X22_loc=Xloc[(nrow(Xloc)/2+1):nrow(Xloc),]}
               mu=c(X11%*%matrix(betas1),X22%*%matrix(betas2))
                   if(!is.null(Xloc)) muloc=c(X11_loc%*%matrix(betas1),X22_loc%*%matrix(betas2))
        }  
#### updating mean
        if(!bivariate){
          #### additive model on the real line
         if(covmatrix$model %in% c(10,18,29,27,38,39,28, 34,40,20))
                                { muloc=muloc + M;   mu=mu + M }     
          ### multiplicative model on the positive real line
          if( (covmatrix$model %in% c(21,26,24,22)))  {emuloc=exp(muloc);emu=exp(mu) }
         }
         else{}
   
##################################################################
##########computing kriging weights##################################
##################################################################
if(is.null(covmatrix$copula))                 { CC = matrix(corri*vvar,nrow=dimat,ncol=dimat2)}
else  {if(covmatrix$copula=="Gaussian")       { CC = matrix(corri,nrow=dimat,ncol=dimat2)    }}
MM=getInv(covmatrix,CC,mse)  
krig_weights = MM$a       #compute (\Sigma^-1) %*% cc
BB=MM$b                   #compute t(cc)%*%(\Sigma^-1) %*% cc

##################################################################
################# simple kriging #################################
##################################################################
if(type_krig=='Simple'||type_krig=='Optim')  {
      
if(!bivariate) ## space and spacetime simple kringing
{  
   #### optimal median predictors ###
   if(type_krig=='Optim'){
               if(SinhAsinhtemp) # Sinh
               {
                     ss= (c(dataT)-c(me))/sqrt(vvm)
                     zz=  sinh(tail*asinh(ss)-sk)
                     #kk=Rfast::Crossprod(as.matrix(zz),krig_weights)
                     kk=crossprod(zz,krig_weights)
                     pp = c(meloc)      +  sqrt(vvm)*   sinh( (1/tail)*(asinh(kk)+sk))
               }
            if(Tukeyh2temp) # Tukeyh2
              {
                     ss= (c(dataT)-c(me))/sqrt(vvm)
                     zz=ifelse(ss>=0,(VGAM::lambertW(t1*ss^2)/t1)^(1/2),-(VGAM::lambertW(t2*ss^2)/t2)^(1/2))
                     kk=crossprod(zz,krig_weights) 
                     pp = c(meloc)      +  sqrt(vvm)* ifelse(kk>=0,kk*exp(0.5*t1*kk^2), kk*exp(0.5*t2*kk^2))
               }
            if(Tukeyhtemp) # Tukeyh
             {
                     ss= (c(dataT)-c(me))/sqrt(vvm)
                     zz=ifelse(ss>=0,(VGAM::lambertW(th*ss^2)/th)^(1/2),-(VGAM::lambertW(th*ss^2)/th)^(1/2))
                     kk=crossprod(zz,krig_weights)
                     pp = c(meloc)      +  sqrt(vvm)* ifelse(kk>=0,kk*exp(0.5*th*kk^2), kk*exp(0.5*th*kk^2))
              }

               if(logGausstemp) # Tukeyh
             {
              ## de oliveira 2006 (equation(2))
                rp=vvm
                datas=as.matrix(c(log(dataT))-c(me-0.5*rp))
                pp = c(meloc-0.5*rp) +  crossprod(datas, krig_weights)
                QQ=diag(as.matrix(diag(covmatrix$param['sill'],dimat2) - BB))
                pp=exp(pp+QQ/2) #/exp(covmatrix$param['sill']/2)
            }
        }
 if(type_krig=='Simple'){
               ############################ optimal linear predictors #######################
               if(covmatrix$model %in% c(1,12,27,38,29,10,18,39,37,28,9, 34,40,20))   ####gaussian, StudenT, two piece  skew gaussian bimodal tukeyh tukey hh
                {                        
                      datas=as.matrix(c(dataT)-c(mu))
                      pp = c(muloc)      +  crossprod(datas,krig_weights)
                }
             }
               ###################################################
               #### gamma weibull loglogistic loggaussian 
if(covmatrix$model %in% c(21,24,26,22)&&type_krig=="Simple")
                      {       ones=rep(1,length(c(dataT)))
                              one=rep(1,length(c(muloc)))
                 datas=as.matrix(c(dataT)/emu-ones)
                pp = c(emuloc) * ( one + crossprod(datas,krig_weights))
      
                      }
               ####log gaussian   simple kriging
#if(covmatrix$model==1&&logGausstemp)   {  
 #                pp = (c(emuloc)+covmatrix$param['sill']/2) +
  #                                        krig_weights %*% (c(dataT)-exp(c(mu)+covmatrix$param['sill']/2))
      #        }     
}     ####simple kriging
      
else  {   ## bivariate  case   cokriging
  dat = c(dataT) - as.numeric(c(rep(covmatrix$param['mean_1'],covmatrix$ns[1]),
                            rep(covmatrix$param['mean_2'],covmatrix$ns[2])))

  if(which==1) pp = c(param$mean_1) + crossprod(dat,krig_weights)
  if(which==2) pp = c(param$mean_2) + crossprod(dat,krig_weights)
      }
#### MSE COMPUTATION ##################

  if(mse) {
                bb=0
            #BB=crossprod(CC,krig_weights)
           # Gaussian,StudentT,skew-Gaussian,two piece linear kriging
        if(covmatrix$model %in% c(1,12,27,38,29,10,18,39,28,40,34,9,20))
                {vv=diag(as.matrix(diag(vvar,dimat2) - BB  + bb)) } ## simple variance  kriging predictor variance

             #gamma  #weibull     #loglogistic #loggaussian
           if(covmatrix$model %in% c(21,26,24,22))  
                    vv=emuloc^2*diag(as.matrix(diag(vvar,dimat2)- BB + bb))

           if(covmatrix$model %in% c(22)&&type_krig=="Optimal")
                    vv =    exp(muloc + covmatrix$param['sill']/2)^2 *diag(as.matrix(diag(exp(vvar),dimat2) - exp(BB+ bb))) 
               }     # end if(mse)
  }  #### end kriging

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
            varpred=matrix(cvv,nrow=tloc,ncol=numloc);
            }
else {pred=c(pp);varpred=cvv}

}#end gaussian standard kriging

} ####

####################################################################################################################################
###################### binomial  binomial negative and poisson poissongamma (inflated) #####################################
####################################################################################################################################
if(covmatrix$model %in% c(2,11,14,19,30,36,16,43,44,45,46,47,57,58))
{

     if(type=="Standard"||type=="standard") {
     
      if(!bivariate){
               if(is.null(MM)) {mu=X%*%betas; mu0=Xloc%*%betas}
               else {mu=MM;mu0=Mloc}        
      }
     if(bivariate)  {mu  = c(X11%*%betas1,X22%*%betas2)}
     KK=0
     if(covmatrix$model %in% c(2,11))    ### binomial
     { if(is.null(nloc)) {KK=rep(round(mean(n)),dimat2)}
          else {if(is.numeric(nloc)) {
                     if(length(nloc)==1) KK=rep(nloc,dimat2) 
                     else                KK=nloc
                                    }
               } 
         if(length(n)==1) {n=rep(n,dimat) }
         if(length(KK)!=dimat2) stop("dimension of nloc is wrong\n")
    }
     if(covmatrix$model==19) KK=min(nloc)
     if(covmatrix$model==16||covmatrix$model==45) KK=n
## ojo que es la covarianza
if(covmatrix$model %in% c(2,11,14,16,19,30,36,43,44,45,46,47,57,58))
{ 
cop=0  
if(!is.null(covmatrix$copula)){ if(covmatrix$copula=="Gaussian") cop=1  }

ccorr=dotCall64::.C64('Corr_c_bin',SIGNATURE = c(rep("double",5),
                                    "integer","integer", "double","double","double",
                                    rep("integer",9), "double", 
                                    "double", "double","integer","integer","double", "integer","integer","double","integer"),   
                                    corri=dotCall64::vector_dc("double",dimat*dimat2),  ccc[,1], ccc[,2], ccc[,3] ,covmatrix$coordt,corrmodel,as.integer(0),
                                    locx, locy,locz,covmatrix$numcoord,numloc,covmatrix$model,  tloc,
                                    KK,n,covmatrix$ns,NS,covmatrix$numtime, 
                                    rep(c(mu),dimat2),  other_nuis, corrparam,covmatrix$spacetime,
                                    bivariate, time,distance,as.integer(which-1), covmatrix$radius,cop,
                                    INTENT = c("w",rep("r",28)),PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
 }

if(cop==1) { 
nuisance$n=n[1]
    corri= gaussian_copula_cov_fast(ccorr$corri,covmatrix$model, nuisance)  ##it works in the stationary case..
     }
else {corri=ccorr$corri}


##################inverse of cov matrix##############################################
        if(!bivariate) cc = matrix(corri,nrow=dimat,ncol=dimat2)
        else           {}
        MM=getInv(covmatrix,cc,mse)  #compute (\Sigma^-1) %*% cc
        krig_weights = MM$a
        BB=MM$b
######################################################################################

       if(type_krig=='Simple'||type_krig=='simple')  {
          ##########################################################
       if(covmatrix$model==30||covmatrix$model==36){  ### poisson
        p0=exp(mu0); pmu=exp(mu)
            if(!bivariate) { 
            datas=c(dataT)-c(pmu)
            pp = c(p0) + crossprod(datas, krig_weights)
             }  ## simple kriging
            else{} #todo
           if(mse)  vvar=p0  ### variance (possibly no stationary)
          }
      ##########################################################
      if(covmatrix$model==46||covmatrix$model==47){  ### poisson gamma
         p0=exp(mu0); pmu=exp(mu)
         if(!bivariate) { 
               datas=c(dataT)-c(pmu)
               pp = c(p0) + crossprod(datas, krig_weights)
          }  ## simple kriging
        else{} #todo
        if(mse)  vvar=p0*(1+p0/covmatrix$param['shape'])  ### variance (possibly no stationary)     
          }
      ##########################################################
       if(covmatrix$model==43||covmatrix$model==44){  ### poisson  inflated
        p=pnorm(covmatrix$param['pmu'])
        p0=exp(mu0); pmu=exp(mu)
            if(!bivariate)  { 
             datas=c(dataT)-(1-p)*c(pmu)
             pp = (1-p)*c(p0) + crossprod(datas,krig_weights)
             }  ## simple kriging
            else{} #todo
           if(mse)  vvar=(1-p)*p0*(1+p*p0)  ### variance (possibly no stationary)
        }
       ##########################################################
       if(covmatrix$model==2||covmatrix$model==11){  ### binomial
        p0=pnorm(mu0); pmu=pnorm(mu)
            if(!bivariate) { 
           #    pp = kk*c(p0) + Rfast::Crossprod(datas,  krig_weights) 
               datas=c(dataT)-n*c(pmu)

               pp = KK*c(p0) + crossprod(datas,  krig_weights) 
            }  ## simple kriging
            else{} #todo
           if(mse) vvar=KK*p0*(1-p0)  ### variance (possibly no stationary
          }
      ###########################################################
           if(covmatrix$model==19){  ### binomial2
            p0=pnorm(mu0); pmu=pnorm(mu)
            if(!bivariate)   
              { datas=c(dataT)-c(n*pmu)
                pp = nloc*c(p0) + crossprod(datas,krig_weights)
              } 
            else{} #todo
          if(mse)    vvar=nloc*p0*t(1-p0) ### variance (possibly no stationary)
          }
        ##########################################################
       if(covmatrix$model==14||covmatrix$model==16){    ###geometric or negative binomial
         p0=pnorm(mu0); pmu=pnorm(mu)
            if(!bivariate) ## space and spacetime
            {   k1=c(p0);k2=c(pmu);
                aa=n*(1-k1)/k1  
                datas=c(dataT)-n*(1-k2)/k2
                pp = aa + crossprod(datas ,krig_weights)
          }
            else{}   #tood
            if(mse) vvar=aa/k1   ### variance (possibly no stationary)
          }
            ##########################################################
            if(covmatrix$model==57||covmatrix$model==58){  ### poisson gamma inflated
         p=pnorm(covmatrix$param['pmu'])
         p0=exp(mu0); pmu=exp(mu)
         if(!bivariate) { 
               datas=c(dataT)-(1-p)*c(pmu)
               pp = (1-p)*c(p0) + crossprod(datas, krig_weights)
          }  ## simple kriging
        else{} #todo
        if(mse)  vvar=p0*(1+ (p0*(1+covmatrix$param['shape']*p)/((1-p)*covmatrix$param['shape'])) )  ### variance (possibly no stationary)     
          }
        ##########################################################
        if(covmatrix$model==45){    ###inflated negative binomial
            p0=pnorm(mu0); pmu=pnorm(mu)
            p=as.numeric(pnorm(covmatrix$param['pmu']))
            if(!bivariate) { k1=c(p0);k2=c(pmu);
                             datas=c(dataT)-(1-p)*n*(1-k2)/k2
                             pp = (1-p)*n*(1-k1)/k1 + crossprod(datas,krig_weights)
                           }
            else{}   #tood
            if(mse) vvar=n*(1-k1)*(1-p)*(1+n*p*(1-k1))/k1^2
          }
        if(mse){
                  bb=0
                  vv = diag(sqrt(tcrossprod(vvar))  - BB  + bb)
                }
        }
##############################################################
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
            varpred=matrix(cvv,nrow=tloc,ncol=numloc);
            }
          else {pred=c(pp);varpred=cvv}
    }}  ##end binary or binomial or geometric kriging  poisson


if(tloc==1)  {pred=c(pred);varpred=c(varpred);varpred2=c(varpred2)}

    # Return the objects list:
   GeoKrig= list(     bivariate=bivariate,
                   coordx = covmatrix$coordx,coordy = covmatrix$coordy,coordz = covmatrix$coordz,
                   coordt = covmatrix$coordt,coordx_dyn=covmatrix$coordx_dyn,
                   covmatrix=covmatrix$covmatrix,corrmodel = corrmodel,copula=copula,data=data,distance = distance,
                   grid=covmatrix$grid,loc=loc_orig,nozero=covmatrix$nozero,ns=covmatrix$ns,X=XX,
                   numcoord = covmatrix$numcoord,numloc= numloc,numtime = covmatrix$numtime,numt = tloc,
                   maxdist=maxdist,maxtime=maxtime,model=model,param = covmatrix$param,pred=pred,n=n,
                   radius=radius,spacetime = covmatrix$spacetime,time=time,type=type,type_krig=type_krig,
                   mse=varpred,mse2=varpred2)
    structure(c(GeoKrig , call = call), class = c("GeoKrig"))
}

}

