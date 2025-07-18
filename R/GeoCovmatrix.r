

GeoCovmatrix <- function(estobj=NULL,coordx, coordy=NULL, coordz=NULL, coordt=NULL,coordx_dyn=NULL,corrmodel, distance="Eucl", grid=FALSE,
                       maxdist=NULL, maxtime=NULL, model="Gaussian", n=1, param, anisopars=NULL,radius=1,
                       sparse=FALSE,taper=NULL, tapsep=NULL, type="Standard",copula=NULL,X=NULL,spobj=NULL)

{
  ########################################################################################################
  ##########  Internal function: computing covariance matrix #############################################
  ########################################################################################################




    Cmatrix <- function(bivariate, coordx, coordy,coordz, coordt,corrmodel, dime, n, ns, NS, nuisance, numpairs,
                           numpairstot, model, paramcorr, setup, radius, spacetime, spacetime_dyn,type,copula,ML,other_nuis)
    {
###################################################################################
############### computing correlation #############################################
###################################################################################

if(model %in% c(1,9,34,12,20,18,39,27,38,29,21,26,24,10,22,40,28,33,42))
{

  if(type=="Standard") {
    if(spacetime) 
       cr=dotCall64::.C64('CorrelationMat_st_dyn2',SIGNATURE = c(rep("double",5),"integer","double","double","double","integer","integer"),
            corr=dotCall64::vector_dc("double",numpairstot),coordx,coordy,coordz,coordt,
            corrmodel,nuisance,paramcorr,radius,ns,NS,
            INTENT = c("w",rep("r",10)),
            PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
     if(bivariate) 
       cr=dotCall64::.C64('CorrelationMat_biv_dyn2',SIGNATURE = c(rep("double",5),"integer","double","double","double","integer","integer"),
            corr=dotCall64::vector_dc("double",numpairstot),coordx,coordy,coordz,coordt,
            corrmodel,nuisance,paramcorr,radius,ns,NS,
            INTENT = c("w",rep("r",10)),
            PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
   if(!bivariate&&!spacetime)  
        cr=dotCall64::.C64('CorrelationMat2',SIGNATURE = c(rep("double",5),"integer","double","double","double","integer","integer"),
            corr=dotCall64::vector_dc("double",numpairstot),coordx,coordy,coordz,coordt,
            corrmodel,nuisance,paramcorr,radius,ns,NS,
            INTENT = c("w",rep("r",10)),
            PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
    
}
###############################################################
if(type=="Tapering")  {
        fname <- 'CorrelationMat_tap'

        if(spacetime) fname <- 'CorrelationMat_st_tap'
       if(bivariate) fname <- 'CorrelationMat_biv_tap'
#############
cr=dotCall64::.C64(fname,SIGNATURE = c("double","double","double","double","double","integer","double","double","double","integer","integer"),
     corr=dotCall64::vector_dc("double",numpairs), coordx,coordy,coordz,coordt,corrmodel,nuisance, paramcorr,radius,ns,NS,
 INTENT = c("w","r","r","r","r","r","r","r", "r", "r", "r"),
            PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
#############
     ## deleting correlation equual  to 1 because there are problems  with hipergeometric function
        sel=(abs(cr$corr-1)<.Machine$double.eps);cr$corr[sel]=0
      }
    }  
    

###################################################################################
###################################################################################

if(model==1)   ## gaussian case
{
  if(!bivariate)
     {
      corr=cr$corr*(1-as.numeric(nuisance['nugget']))
      vv=nuisance['sill']
     }
if(bivariate){}
}
###############################################################
if(model==9)  {  ## TukeyGH

if(!bivariate)   {
          corr=cr$corr*(1-as.numeric(nuisance['nugget']))
          h=as.numeric(nuisance['tail'])
          g=as.numeric(nuisance['skew'])
          ss=as.numeric(nuisance['sill'])

  if(abs(g)<=1e-05 && h<=1e-05) { vv=ss }
  if(abs(g)>1e-05&&h<=1e-05) {  #
              aa=( -exp(g^2)+exp(g^2*2))*g^(-2)
              corr <- (( -exp(g^2)+exp(g^2*(1+corr)))*g^(-2))/aa
              vv=  aa*ss
              }
   if(abs(g)<=1e-05&& h>1e-05){ ##
              aa=(1-2*h)^(1.5) # variance
              corr <- aa*corr/( (1-h)^2-h^2*corr^2 )^(1.5)
              vv=aa*ss
               }
  if(h>1e-05&&abs(g)>1e-05){ # ok
                  rho=corr; rho2=rho*rho;
                  tail=h; eta=g
                  eta2=eta*eta; tail2=tail*tail;
                  u=1-tail;a=1+rho;
                  A1=exp(a*eta2/(1-tail*a));
                  A2=2*exp(0.5*eta2*  (1-tail*(1-rho2))  / (u*u- tail2*rho2)  );
                  A3=eta2*sqrt(u*u- rho2*tail2)
                  mu=(exp(eta2/(2*u))-1)/(eta*sqrt(u));
                  cova=ss*((A1-A2+1)/A3-mu*mu)
                  vari=ss*((exp(2*eta2/(1-2*tail))-2*exp(eta2/(2*(1-2*tail)))+1)/(eta2*
                           sqrt(1-2*tail))-mu*mu);
                  corr <-cova/vari
                  vv=vari
            }
}
if(bivariate){}
}
###############################################################
if(model==40)  {  ## TukeyH2

if(!bivariate)    {
    corr=cr$corr*(1-as.numeric(nuisance['nugget']))
    hr=as.numeric(nuisance['tail1'])
    hl=as.numeric(nuisance['tail2'])

x1=1-(1-corr^2)*hr
y1=1-(1-corr^2)*hl
x2=(1-hr)^2-(corr*hr)^2
y2=(1-hl)^2-(corr*hl)^2
g=1-hl-hr+(1-corr^2)*hl*hr
h1=sqrt(1-corr^2/(x1^2))+(corr/x1)*asin(corr/x1)
h2=sqrt(1-corr^2/(y1^2))+(corr/y1)*asin(corr/y1)
h3=sqrt(1-corr^2/(x1*y1))+sqrt(corr^2/(x1*y1))*asin(sqrt(corr^2/(x1*y1)))
p1=x1*h1/(2*pi*(x2)^(3/2))+corr/(4*(x2)^(3/2))
p2=y1*h2/(2*pi*(y2)^(3/2))+corr/(4*(y2)^(3/2))
p3=-(x1*y1)^(1/2)*h3/(2*pi*(g)^(3/2))+corr/(4*(g)^(3/2))

   mm=(hr-hl)/(sqrt(2*pi)*(1-hl)*(1-hr))
   vv1=0.5*(1-2*hl)^(-3/2)+0.5*(1-2*hr)^(-3/2)-(mm)^2
   corr=(p1+p2+2*p3-mm^2)/vv1
  vv=as.numeric(nuisance['sill'])*vv1
       }
if(bivariate){}
}

###############################################################
if(model==34)  {  ## TukeyH

if(!bivariate)    {
       corr=cr$corr*(1-as.numeric(nuisance['nugget']))
    h=nuisance['tail']
   if(!h)     vv=as.numeric(nuisance['sill'])
   if(h){     aa=(1-2*h)^(-1.5) # variance
              #corr <- aa*corr/( (1-h)^2-h^2*corr^2 )^(1.5)
              #corr <- (-corr/((1+h*(corr-1))*(-1+h+h*corr)*(1+h*(-2+h-h*corr^2))^0.5))/aa
              corr =(corr*(1-2*h)^(1.5))/((1-h)^2-(h*corr)^2)^(1.5)
              vv=aa*as.numeric(nuisance['sill'])
              }
       }
if(bivariate){}
}
######################################################################
if(model==12)   ##  student case
    {
if(!bivariate){
   corr=cr$corr*(1-as.numeric(nuisance['nugget']))
nu=as.numeric(1/nuisance['df'])
#corr=exp(log(nu-2)+2*lgamma(0.5*(nu-1))-log(2)-2*lgamma(nu/2)+log(Re(hypergeo::hypergeo(0.5,0.5, nu/2,corr^2)))+log(corr))
if(nu<170) corr=((nu-2)*gamma((nu-1)/2)^2*Re(hypergeo::hypergeo(0.5,0.5 ,nu/2 ,corr^2))*cr$corr)/(2*gamma(nu/2)^2)
vv=as.numeric(nuisance['sill'])*(nu)/(nu-2)
}
if(bivariate){}
}

###############################################################
  if(model==18)   ##  skew student case
    {
     if(!bivariate)
{
    corr=cr$corr*(1-as.numeric(nuisance['nugget']))
    nu=as.numeric(1/nuisance['df']); sk=as.numeric(nuisance['skew'])
    skew2=sk*sk;l=nu/2; f=(nu-1)/2; w=sqrt(1-skew2);y=corr;
    CorSkew=(2*skew2/(pi*w*w+skew2*(pi-2)))*(sqrt(1-y*y)+y*asin(y)-1)+w*w*cr$corr/(w*w+skew2*(1-2/pi)) ;
    mm=sqrt(nu)*gamma(f)*sk/(sqrt(pi)*gamma(l));
    corr=(pi*(nu-2)*gamma(f)^2/(2*(pi*gamma(l)^2-skew2*(nu-2)*gamma(f)^2)))*(Re(hypergeo::hypergeo(0.5,0.5,l,y*y))
                                *((1-2*skew2/pi)*CorSkew+2*skew2/pi)-2*skew2/pi);
    vv=as.numeric(nuisance['sill'])*(nu/(nu-2)-mm*mm);
}
if(bivariate){}
}
###############################################################

if (model == 20) {  ## SAS
  # Estrazione parametri
  nugget <- as.numeric(nuisance['nugget'])
  d <- as.numeric(nuisance['tail'])     # deve essere > 0
  e <- as.numeric(nuisance['skew'])     # può essere qualsiasi valore
  sill <- as.numeric(nuisance['sill'])
  mm=sinh(e/d)*exp(0.25)*(besselK(.25,(d+1)/(2*d))+besselK(.25,(1-d)/(2*d)))/(sqrt(8*pi))
  vs=cosh(2*e/d)*exp(0.25)*(besselK(.25,(d+2)/(2*d))+besselK(0.25,(2-d)/(2*d)))/(sqrt(32*pi))-0.5-mm^2

  corr_base <- cr$corr * (1 - nugget)
  corr=corrsas(corr_base, e, d)        #external function for sas
  vv <- sill * vs
}



###############################################################
 if(model==39)   ##  two piece bimodal case
    {
if(!bivariate)
{
 corr=cr$corr*(1-as.numeric(nuisance['nugget']))
 nu=as.numeric(nuisance['df']); sk=as.numeric(nuisance['skew'])
 delta=as.numeric(nuisance['shape'])
 alpha=2*(delta+1)/nu
 nn=2^(1-alpha/2)
 ll=qnorm((1-sk)/2)
 p11=pbivnorm::pbivnorm(ll,ll, rho = cr$corr, recycle = TRUE)
 corr2=corr^2;sk2=sk^2
 a1=Re(hypergeo::hypergeo(-1/alpha ,-1/alpha,nu/2,corr2))
 a3=3*sk2 + 2*sk + 4*p11 - 1

 MM=(2^(2/alpha)*(gamma(nu/2 + 1/alpha))^2)
 vari=2^(2/alpha)*(gamma(nu/2 + 2/alpha))*gamma(nu/2)* (1+3*sk2) - sk2*2^(2/alpha+2)*gamma(nu/2+1/alpha)^2
 corr= MM*(a1*a3-4*sk2)/vari
 vv=as.numeric(nuisance['sill'])*vari/(nn^(2/alpha)*gamma(nu/2)^2)
}

if(bivariate){}
}
###############################################################
if(model==27)   ##  two piece student case case
    {
  if(!bivariate)
{
          corr=cr$corr*(1-as.numeric(nuisance['nugget']))
          nu=as.numeric(1/nuisance['df']); sk=as.numeric(nuisance['skew'])
          corr2=corr^2;sk2=sk^2
          a1=Re(hypergeo::hypergeo(0.5,0.5,nu/2,corr2))
          a2=corr*asin(corr) + (1-corr2)^(0.5)
          ll=qnorm((1-sk)/2)
          p11=pbivnorm::pbivnorm(ll,ll, rho = cr$corr, recycle = TRUE)
          a3=3*sk2 + 2*sk + 4*p11 - 1
          KK=( nu*(nu-2)*gamma((nu-1)/2)^2) / (nu*pi*gamma(nu/2)^2*(3*sk2+1)-4*sk2*nu*(nu-2)*gamma((nu-1)/2)^2 )
          corr= KK*(a1*a2*a3-4*sk2);
          ttemp=gamma(0.5*(nu-1))/gamma(0.5*nu)
          vv=as.numeric(nuisance['sill'])*((nu/(nu-2))*(1+3*sk2) - 4*sk2*(nu/pi)*ttemp^2)
      }

if(bivariate){}
}
###############################################################
if(model==38)   ##  two piece tukey h  case
    {
      if(!bivariate)
{
          corr=cr$corr*(1-as.numeric(nuisance['nugget']))
          tail=as.numeric(nuisance['tail']); sk=as.numeric(nuisance['skew'])
          corr2=corr^2;sk2=sk^2;
          gg2=(1-(1-corr2)*tail)^2
          xx=corr2/gg2
          A=(asin(sqrt(xx))*sqrt(xx)+sqrt(1-xx))/(1-xx)^(1.5)
          ll=qnorm((1-sk)/2)
          p11=pbivnorm::pbivnorm(ll,ll, rho = cr$corr, recycle = TRUE)
          a3=3*sk2 + 2*sk + 4*p11 - 1
          mm=8*sk2/(pi*(1-tail)^2);
          ff=(1+3*sk2)/(1-2*tail)^(1.5)
          M=(2*(1-corr2)^(3/2))/(pi*gg2)
          corr=  (M*A*a3-mm)/( ff- mm)
          vv= as.numeric(nuisance['sill'])*((1-2*tail)^(-1.5)* (1+3*(sk2)) - 4*(sk2)*2/(pi*(1-tail)^2))
  }
if(bivariate){}
}
###############################################################
if(model==29)   ##  two piece gaussian case
    {
      if(!bivariate)
{
          corr1=cr$corr*(1-as.numeric(nuisance['nugget']))
          sk=as.numeric(nuisance['skew']);
          corr2=sqrt(1-corr1^2); sk2=sk^2
          ll=qnorm((1-sk)/2)
          p11=pbivnorm::pbivnorm(ll,ll, rho = cr$corr, recycle = TRUE)
          KK=3*sk2+2*sk+ 4*p11 - 1
          corr=(2*((corr2 + corr1*asin(corr1))*KK)- 8*sk2)/(3*pi*sk2  -  8*sk2   +pi   )
          vv=as.numeric(nuisance['sill'])*(1+3*sk2-8*sk2/pi)
  }
if(bivariate){}
}
###############################################################
if(model==10)  {  ##  skew Gaussian case

          if(!bivariate){
              corr=cr$corr*(1-as.numeric(nuisance['nugget']))
              sk=as.numeric(nuisance['skew'])
              corr2=cr$corr^2; ; sk2=sk^2; vv=as.numeric(nuisance['sill'])
              corr=(2*sk2)*(sqrt(1-corr2) + cr$corr*asin(cr$corr)-1)/(pi*vv+sk2*(pi-2)) + (corr*vv)/(vv+sk2*(1-2/pi))
              vv=vv+as.numeric(nuisance['skew'])^2*(1-2/pi)
     }
      if(bivariate){}
}
#################################################################################
################ covariance matrix for models defined on the real line ##########
#################################################################################
if(model %in% c(1,9,34,12,20,18,39,27,38,29,10,40)){
if(!bivariate)
{
 if(type=="Standard"){
     # Builds the covariance matrix:
        varcov <-  diag(dime)
        varcov[lower.tri(varcov)] <- corr
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr
        varcov=varcov*vv
      }
    if(type=="Tapering")  {
          vcov <- vv*corr;
          varcov <- new("spam",entries=vcov,colindices=setup$ja,
                             rowpointers=setup$ia,dimension=as.integer(rep(dime,2)))
          diag(varcov)=vv

        }
}
#####
    if(bivariate)      {
       if(type=="Standard"){
          corr <- cr$corr
          varcov<-diag(dime)
          varcov[lower.tri(varcov,diag=TRUE)] <- corr
          varcov <- t(varcov)
          varcov[lower.tri(varcov,diag=TRUE)] <- corr
        }
        if(type=="Tapering")  {

          varcov <-new("spam",entries=cr$corr,colindices=setup$ja,
                         rowpointers=setup$ia,dimension=as.integer(rep(dime,2)))
        }
      }
}
#################################################################################
############# covariance models defined on a bounded support of the real line ###
#################################################################################
if(model==28)   ##  beta case
    {
      if(!bivariate) {
         corr=cr$corr*(1-as.numeric(nuisance['nugget']))
         corr2=corr^2
         shape1=as.numeric(nuisance['shape1']);
         shape2=as.numeric(nuisance['shape2']);
         ssup=as.numeric(nuisance['max'])-as.numeric(nuisance['min'])
         cc=0.5*(shape1+shape2)
         vv=ssup^2*shape1*shape2/((cc+1)*(shape1+shape2)^2)        ## variance remember min and max!
         idx=which(abs(corr2)>1e-10);corr22=corr2[idx]
######################
  #nu=shape1;alpha=shape2
  nu2=shape1/2;alpha2=shape2/2
  res=0;ss=0;k=0
  while(k<=100){
    p1=2*(lgamma(cc+k)-lgamma(cc)+lgamma(nu2+1+k)-lgamma(nu2+1))
    p2=lgamma(k+1)+(lgamma(nu2+k)-lgamma(nu2))+2*(lgamma(cc+1+k)-lgamma(cc+1))
    b1=p1-p2
    b2=log(hypergeo::genhypergeo(U=c(cc+k,cc+k,alpha2), L=c(cc+k+1,cc+k+1), polynomial=TRUE,maxiter=1000, z=corr22))
    b3=k*log(corr22)
    sum=exp(b1+b2+b3)
    res=res+sum
    if (all(sum<1e-6)){
      break
    } else{
      A=res
    }
    k=k+1
  }
######################
         corr[idx]=A
         corr=shape1*(cc + 1 ) * ((1-corr2)^(cc) *corr -1)/shape2 ## correlation
         corr[-idx]=0
        }
     if(bivariate){}
}

if(model==33||model==42)   ##  Kumaraswamy case
    {
      if(!bivariate) {
#######################

      corr=cr$corr*(1-as.numeric(nuisance['nugget']))
      if(model==33){
                ga=as.numeric(nuisance['shape2'])
                eta=as.numeric(nuisance['shape1'])
                   }
      if(model==42){
                ga=as.numeric(nuisance['shape2'])
                eta=log1p(-(1+exp(mu))^(-ga))/log(0.5)
                   }

      mm=eta*beta(1+1/ga,eta)
      NN=length(corr);
      res=double(NN)
      bb=.C("corr_kuma_vec",as.double(corr),as.double(eta),as.double(ga),
      res=as.double(res),as.integer(NN),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
      corr=bb$res
      vv=eta*beta(1+2/ga,eta)-mm^2
      }
      if(bivariate){}
}

#################################################################################
# covariance matrix for models defined on a bounded support of the  real line  ##
#################################################################################
if(model %in% c(28,33,42)){
if(!bivariate)
{
 if(type=="Standard"){
     # Builds the covariance matrix:
        varcov <-  diag(dime)
        varcov[lower.tri(varcov)] <- corr
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr
        varcov=varcov*vv
      }
    if(type=="Tapering")  {
          vcov <- vv*corr;
          varcov <- new("spam",entries=vcov,colindices=setup$ja,
                             rowpointers=setup$ia,dimension=as.integer(rep(dime,2)))
          diag(varcov)=vv
        }
}
#####
    if(bivariate)      {
       if(type=="Standard"){
          corr <- cr$corr
          varcov<-diag(dime)
          varcov[lower.tri(varcov,diag=TRUE)] <- corr
          varcov <- t(varcov)
          varcov[lower.tri(varcov,diag=TRUE)] <- corr
        }
        if(type=="Tapering")  {
          varcov <-new("spam",entries=cr$corr,colindices=setup$ja,
                         rowpointers=setup$ia,dimension=as.integer(rep(dime,2)))
        }
      }
}

#################################################################################
################ wrapped Gaussian RF                     #######################
#################################################################################


if(model==13)   ##  wrapped
    {
      if(!bivariate) {
         corr=cr$corr*(1-as.numeric(nuisance['nugget']))
         vv=sinh(as.numeric(nuisance['sill']))
         corr=sinh(as.numeric(nuisance['sill'])*corr)/vv
        }
     if(bivariate){}
    }

#################################################################################
################ models defined on the positive real line #######################
#################################################################################
if(model %in% c(21,26,24,22)) 
    {

         mu = ML   ## mean  function

if(model==21)   ##  gamma case
    {
      if(!bivariate) {
         corr=cr$corr*(1-as.numeric(nuisance['nugget']))
         corr=corr^2   ### gamma correlation
         vv=exp(mu)^2 * 2/as.numeric(nuisance['shape'])
        }
     if(bivariate){}
    }
###############################################################
if(model==26)   ##  weibull case
    {
        if(!bivariate) {
         corr=cr$corr*(1-as.numeric(nuisance['nugget']))
         sh=as.numeric(nuisance['shape'])
         vv=exp(mu)^2 * (gamma(1+2/sh)/gamma(1+1/sh)^2-1)
         bcorr=    gamma(1+1/sh)^2/(gamma(1+2/sh)-gamma(1+1/sh)^2)
         corr=bcorr*((1-corr^2)^(1+2/sh)*Re(hypergeo::hypergeo(1+1/sh,1+1/sh ,1 ,corr^2))-1)
        }
  if(bivariate) {}
    }
###############################################################
    if(model==24)   ## log-logistic case
    {
       if(!bivariate) {
       corr=cr$corr*(1-as.numeric(nuisance['nugget']))
       sh=as.numeric(nuisance['shape'])
       vv=(exp(mu))^2* (2*sh*sin(pi/sh)^2/(pi*sin(2*pi/sh))-1)
       corr= ((pi*sin(2*pi/sh))/(2*sh*(sin(pi/sh))^2-pi*sin(2*pi/sh)))*
             (Re(hypergeo::hypergeo(-1/sh,-1/sh ,1 ,corr^2))* Re(hypergeo::hypergeo(1/sh,1/sh ,1 ,corr^2)) -1)
        }
   if(bivariate){}
    }
###############################################################
if(model==22)  {  ## Log Gaussian
      if(!bivariate) {
      corr=cr$corr*(1-as.numeric(nuisance['nugget']))
      vvar=as.numeric(nuisance['sill'])
      corr=(exp(vvar*corr)-1)/(exp(vvar)-1)
            vv=(exp(mu))^2*(exp(vvar)-1) 
             }
    if(bivariate){}
  }
#################################################################################
################ covariance matrix for models defined on the positive real line #
#################################################################################

       if(type=="Standard"){
     # Builds the covariance matrix:
        varcov <-  diag(dime)
        varcov[lower.tri(varcov)] <- corr
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr
      }
    if(type=="Tapering")  {
          varcov <- new("spam",entries=corr,colindices=setup$ja,
                             rowpointers=setup$ia,dimension=as.integer(rep(dime,2)))
          diag(varcov)=1
        }
    V=vv%*%t(vv)
   
    varcov=varcov*sqrt(V)
}

###############################################################
################################ end continous models #########
###############################################################

###############################################################
################################ start discrete #models ########
###############################################################
if(model %in% c(2,11,30,16,14,43,45,46,57,58)){ #  binomial (negative)Gaussian type , Poisson (inflated) poissongamma poisson gamma zip

if(model==2||model==11)
{if(length(n)==1) n=rep(n,dime)}
if(!bivariate){
mu=ML

if(type=="Standard")  {

    fname <-'CorrelationMat_dis2'
    if(spacetime) fname <- 'CorrelationMat_st_dyn_dis2'

  cr=dotCall64::.C64(fname,SIGNATURE =
      c("double","double","double","double","double","integer","double","integer","double","double","double","integer","integer","integer"),
        corr=dotCall64::vector_dc("double",numpairstot), coordx,coordy,coordz, coordt,corrmodel,c(mu), n,other_nuis,paramcorr,radius,ns,NS,model,
  INTENT = c("w","r","r","r","r","r","r","r", "r", "r","r", "r", "r", "r"),
             PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)

     corr=cr$corr # ojo que corr en este caso es una covarianza y ya va con el nugget
     varcov <-  diag(dime)
     varcov[lower.tri(varcov)] <- corr
     varcov <- t(varcov)
     varcov[lower.tri(varcov)] <- corr

}
############################
############################
 if(type=="Tapering")  {
        tap <-new("spam",entries=setup$taps,colindices=setup$ja,
                         rowpointers=setup$ia,dimension=as.integer(rep(dime,2)))
        idx=spam::triplet(tap)$indices
        
        fname <- "CorrelationMat_dis_tap"
        if(spacetime) fname <- "CorrelationMat_st_dis_tap"
     
if(spacetime){


cr<-dotCall64::.C64("CorrelationMat_st_dis_tap",corr=dotCall64::vector_dc("double",numpairs),coordx,coordy,coordz,coordt,
    corrmodel,other_nuis,paramcorr,radius,ns,NS,n[idx[,1]],n[idx[,2]],mu[idx[,1]],mu[idx[,2]],model,
    SIGNATURE=c("double","double","double","double","double","integer","double",
        "double","double","integer","integer","integer","integer","double","double","integer"),
    INTENT=c("w",rep("r",15)),
          PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
}
else
{
cr<-dotCall64::.C64("CorrelationMat_dis_tap",corr=dotCall64::vector_dc("double",numpairs),coordx,coordy,coordz,coordt,corrmodel,
    other_nuis,paramcorr,radius,ns,NS,n[idx[,1]],n[idx[,2]],mu[idx[,1]],mu[idx[,2]],model,
   SIGNATURE=c("double","double","double","double","double","integer","double",
        "double","double","integer","integer","integer","integer","double","double","integer"),
    INTENT=c("w",rep("r",15)),
           PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
}
    varcov <-new("spam",entries=cr$corr,colindices=setup$ja,
            rowpointers=setup$ia,dimension=as.integer(rep(dime,2)))

        }
            ## updating the diagonal with variance
  if(model %in% c(2,11)) { pg=pnorm(mu);  vv=pg*(1-pg)*n; diag(varcov)=vv }
  if(model %in% c(14))   { pg=pnorm(mu); vv=  (1-pg)/pg^2; diag(varcov)=vv }
  if(model %in% c(16))   { pg=pnorm(mu); vv=n*(1-pg)/pg^2; diag(varcov)=vv }
  if(model %in% c(30))   { vv=exp(mu); diag(varcov)=vv } ## poisson


  if(model %in% c(46))  
   { mm=exp(mu); vv=mm*(1+mm/as.numeric(param$shape)); diag(varcov)=vv } ## poissongamma

  if(model %in% c(43))   { mm=exp(mu); pg=pnorm(param$pmu)
                           vv=(1-pg)*mm*(1+pg*mm)
                           diag(varcov)=vv }

  if(model %in% c(45))   {
                            p=pnorm(param$pmu)
                            MM=pnorm(mu)
                            vv=n*(1-MM)*(1-p)*(1+n*p*(1-MM)) /MM^2
                            diag(varcov)=vv 
                        }

  if(model %in% c(57))   { mm=exp(mu); vv=mm*(1+mm/as.numeric(param$shape))
                            pg=pnorm(param$pmu)
                           vv=(1-pg)*vv*(1+pg*vv)
                           diag(varcov)=vv }
}

}
###############################################################
################################ end discrete #models #########
###############################################################
return(varcov)
}

#############################################################################################
#################### end internal function ##################################################
#############################################################################################


##############################################################################
###### extracting GeoFit object informations if necessary              #######
##############################################################################
if(!is.null(estobj)){
   if(!inherits(estobj,"GeoFit"))
               stop("need  a 'GeoFit' object as input\n")
if(!estobj$grid){  #not regular grid 

 if(!estobj$bivariate){  if(is.null(estobj$coordx_dyn)) coordx=cbind(estobj$coordx,estobj$coordy,estobj$coordz)
                         else cord=estobj$coordx_dyn
                      } ## spatial (temporal) non regular case
 else  {    if(is.null(estobj$coordx_dyn))  { coordx=estobj$coordx[1:estobj$ns[1]]    # bivariate not dynamic    
                                              coordy=estobj$coordy[1:estobj$ns[2]] 
                                            }  
            else {coordx_dyn=estobj$coordx_dyn}                                      # bivariate  dynamic  
       }
 }
else  { coordx=estobj$coordx; 
        coordy=estobj$coordy
        coordz=estobj$coordz
      }
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
################# end extracting information###################################################
if( !is.character(corrmodel)|| is.null(CkCorrModel(corrmodel)))       stop("the name of the correlation model is wrong")
 if(is.null(CkModel(model))) stop("The name of the  model  is not correct\n")

bivariate<-CheckBiv(CkCorrModel(corrmodel))
spacetime<-CheckST(CkCorrModel(corrmodel))
space=!spacetime&&!bivariate

##############################################################################
###### extracting sp object informations if necessary              ###########
##############################################################################

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



#######################################################
     ## setting zero mean and nugget if no mean or nugget is fixed
    if(!bivariate){
    if(is.null(param$mean)) param$mean<-0
    if(is.null(param$nugget)) param$nugget<-0  
    #if(length(param$mean)>1) MM=fixed$mean
    if(is.null(param$sill)) param$sill<-1  
}
    else{
    if(is.null(param$mean_1)) param$mean_1<-0
    if(is.null(param$mean_2)) param$mean_2<-0
    if(is.null(param$nugget_1)) param$nugget_1<-0
    if(is.null(param$nugget_2)) param$nugget_2<-0 
    if(is.null(param$sill_1)) param$sill_1<-0
    if(is.null(param$sill_2)) param$sill_2<-0 

}

    if(is.null(coordx_dyn)){
    unname(coordx);unname(coordy);unname(coordz)}
    else{coordx=NULL;coordy=NULL;coordz=NULL}


    #if the covariance is compact supported  and option sparse is used
    #then set the code as a tapering and an object spam is returned
if(sparse) {
    covmod=CkCorrModel(corrmodel)
    if(covmod %in% c(10,11,13,15,19,6,7,21,22,23,30,
                     63,64,65,66,67,68,
                     69,70,71,72,73,74,75,76,77,78,
                     111,112,129,113,114,131,132,130,134,
                     115,116,120))
    {
      type="Tapering"
      if(bivariate){
                  #taper="unit_matrix_biv"
                  if(covmod %in% c(111,113,115,130)) maxdist=c(param$scale,param$scale,param$scale)
                  if(covmod %in% c(112,114,116,134)) maxdist=c(param$scale_1,param$scale_12,param$scale_2)
                  if(covmod %in% c(120,129,131,132)) maxdist=c(param$scale_1,0.5*(param$scale_1+param$scale_2),param$scale_2)
      }
      if(spacetime)
      {

      maxdist=param$scale_s;maxtime=param$scale_t

       if(covmod==63||covmod==65||covmod==67) {  tapsep=c(param$power2_s,param$power_t,param$scale_s,param$scale_t,param$sep) }
       if(covmod==64||covmod==66||covmod==68) {  tapsep=c(param$power_s,param$power2_t,param$scale_s,param$scale_t,param$sep) }
    }
      if(space){  ### spatial Gen Wend (reparametrized)
        maxdist=param$scale
        if(covmod==6)   maxdist=as.numeric(param$scale*exp((lgamma(2*param$smooth+1/param$power2+1)-lgamma(1/param$power2))/ (1+2*param$smooth) ))
        if(covmod==7)   maxdist=as.numeric(param$scale*exp((lgamma(2*param$smooth+param$power2+1)-lgamma(param$power2))/ (1+2*param$smooth) ))
        if(covmod==23)
              {
              log_num <- log(2) * (2 * param$smooth + 1) + lgamma((1/param$power2 + 1) / 2 + param$smooth) + lgamma((1/param$power2 + 2 + 1) / 2 + 2 * param$smooth)
              log_den <- lgamma(1/param$power2 / 2) +lgamma((1/param$power2 + 2) / 2 + param$smooth)
              maxdist <- param$scale*exp(log_num - log_den)^(1/(1+2*param$smooth))
              }
        
        if(covmod==30)    {
              log_num <- log(2) * (2 * param$smooth + 1) + lgamma((param$power2 + 1) / 2 + param$smooth) + lgamma((param$power2 + 2 + 1) / 2 + 2 * param$smooth)
              log_den <- lgamma(param$power2 / 2) +lgamma((param$power2 + 2) / 2 + param$smooth)
              maxdist <- param$scale*exp(log_num - log_den)^(1/(1+2*param$smooth))
              }
      }
  }
  taper=corrmodel
}
##### sill parameter fixed for some models

if(!bivariate){

if(model %in% c("Weibull","Poisson","Binomial","Gamma","LogLogistic",
        "BinomialNeg","Bernoulli","Geometric","Gaussian_misp_Poisson",
        'PoissonZIP','Gaussian_misp_PoissonZIP','BinomialNegZINB', 'PoissonGammaZIP','PoissonGammaZIP1','PoissonGamma',
        'PoissonZIP1','Gaussian_misp_PoissonZIP1','BinomialNegZINB1',
        'Beta2','Kumaraswamy2','Beta','Kumaraswamy')){
if(is.null(param$sill)) param$sill=1
else                    param$sill=1
}
}


coords=NULL
coordz=NULL

if(is.null(coordx_dyn)){
 #######################################
    if(grid) { coords=as.matrix(expand.grid(coordx,coordy)); grid=FALSE }
    else     {  if(!is.null(coordy)) { coords=as.matrix(cbind(coordx,coordy,coordz))}
                else                 { coords=as.matrix(coordx);}
                if(ncol(coords)==3) coordz=coords[,3]
             }         
}
    coords_orig=coords
 #######################################  
 
    if(!is.null(anisopars)) {  coords=GeoAniso(coords,c(anisopars$angle,anisopars$ratio))}
#######################################


    checkinput <- CkInput(coords[,1], coords[,2],coordz, coordt, coordx_dyn, corrmodel, NULL, distance, "Simulation",
                             NULL, grid, NULL, maxdist, maxtime,  model=model, n,  NULL,
                              param, radius, NULL, taper, tapsep,  "Standard", NULL, NULL,copula,X)




    
    if(!is.null(checkinput$error)) stop(checkinput$error)
    spacetime_dyn=FALSE
    if(!is.null(coordx_dyn))  spacetime_dyn=TRUE
    # Initialising the parameters:
    
    initparam <- StartParam(coords[,1], coords[,2],coordz, coordt,coordx_dyn, corrmodel, NULL, distance, "Simulation",
                           NULL, grid, NULL, maxdist, NULL,maxtime, model, n,
                           param, NULL, NULL, radius, NULL, taper, tapsep,  type, 
                           type, FALSE,copula,X,FALSE,FALSE)

    if(is.null(coordz)) cc=cbind(initparam$coordx,initparam$coordy,0)
    else cc=cbind(initparam$coordx,initparam$coordy,initparam$coordz)

    if(!spacetime_dyn) dime=initparam$numcoord*initparam$numtime
    else               dime=sum(initparam$ns)

    if(!initparam$bivariate) numpairstot=dime*(dime-1)*0.5
    if(initparam$bivariate)  numpairstot=dime*(dime-1)*0.5+dime
    if(!is.null(initparam$error)) stop(initparam$error)
    setup<-initparam$setup
    if(initparam$type=="Tapering")
    {
       corr <- double(initparam$numpairs)
       if(ncol(cc)==2)  ccz=corr
       if(ncol(cc)==3)  ccz=cc[,3]

       #tapmod <- setup$tapmodel
       ### unit taperssss ####
       if(sparse){
           if(spacetime) tapmod=230
           if(bivariate) tapmod=147
           if(space) tapmod=36
           }
       else(tapmod=CkCorrModel(taper))
   
      if(initparam$spacetime) fname= "CorrelationMat_st_tap"
      else {
         if(!initparam$spacetime) fname= "CorrelationMat_tap"
         if(initparam$bivariate) fname= "CorrelationMat_biv_tap"
     }

if(is.null(tapsep)) tapsep=1
   tp=dotCall64::.C64(fname,SIGNATURE =
         c("double","double","double","double","double","integer","double","double","double","integer","integer"),
        tapcorr=dotCall64::vector_dc("double",initparam$numpairs),
        cc[,1],cc[,2],ccz,initparam$coordt,tapmod, 1,tapsep,1,1,1,
  INTENT = c("w","r","r","r","r","r","r","r", "r", "r","r"),
             PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)

        setup$taps<-tp$tapcorr

    }
 ### end tapering   
    if(is.null(X))  initparam$X=as.matrix(rep(1,dime))

    if(bivariate) {if(is.null(X))  initparam$X=as.matrix(rep(1,initparam$ns[1]+initparam$ns[2])) }
    if(!space){
          initparam$NS=cumsum(initparam$ns);
            if(spacetime_dyn){  initparam$NS=c(0,initparam$NS)[-(length(initparam$ns)+1)]}
            else{               initparam$NS=rep(0,initparam$numtime)}
    }
    if(is.null(initparam$NS)) initparam$NS=0

    if(initparam$model %in% c(43,45)) initparam$namesnuis=initparam$namesnuis[!initparam$namesnuis %in% "nugget"]



other_nuis=NULL
ML=0
##### sill parameter fixed for some models
if(!bivariate){

if(model %in% c("Weibull","Poisson","Binomial","Gamma","LogLogistic",
        "BinomialNeg","Bernoulli","Geometric","Gaussian_misp_Poisson",
        'PoissonZIP','Gaussian_misp_PoissonZIP','BinomialNegZINB', 'PoissonGammaZIP','PoissonGammaZIP1','PoissonGamma',
        'PoissonZIP1','Gaussian_misp_PoissonZIP1','BinomialNegZINB1',"LogGaussian",
        'Beta2','Kumaraswamy2','Beta','Kumaraswamy')) 
  {
  ### setting mean because mean affect variance in this case
  parc=initparam$param[initparam$namescorr]
  sel1=!(names(initparam$param) %in% names(parc))
  inip=initparam$param[sel1] ## deleting corr parameters

  sel=substr(names(inip),1,4)=="mean"
  beta=as.numeric(inip[sel]) ## mean parameters
  other_nuis2=inip[!sel]
  other_nuis=as.numeric(other_nuis2[order(names(other_nuis2))]) # other nuis parameters ordered by names
  if((dim(initparam$X)[2])>1){ ML=initparam$X%*%c(beta) }
  else                       { if(sum(sel)>1)   ML=beta 
                               else  ML=initparam$X*c(beta) }
                             } 
}


######  calling main correlation functions
if(is.null(copula)) {
    covmatrix<- Cmatrix(initparam$bivariate,cc[,1],cc[,2],cc[,3],initparam$coordt,initparam$corrmodel,dime,n,initparam$ns,
                        initparam$NS, initparam$param[initparam$namesnuis],
                        initparam$numpairs,numpairstot,initparam$model,
                        initparam$param[initparam$namescorr],setup,initparam$radius,initparam$spacetime,spacetime_dyn,initparam$type,copula,ML,other_nuis)
                    }
else {

 if(copula=="Gaussian") 
    covmatrix<- Cmatrix_copula(initparam$bivariate,cc[,1],cc[,2],cc[,3],initparam$coordt,initparam$corrmodel,dime,n,initparam$ns,
                        initparam$NS, initparam$param[initparam$namesnuis],
                        initparam$numpairs,numpairstot,initparam$model,
                        initparam$param[initparam$namescorr],setup,initparam$radius,initparam$spacetime,spacetime_dyn,initparam$type,copula,ML,other_nuis)
                    }


    if(type=="Tapering") sparse=TRUE

    # Delete the global variables:
    if(!space)
                      .C('DeleteGlobalVar', PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
    else                     
                      .C('DeleteGlobalVar2', PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
 ##### setting names nuis parameters
    if(initparam$bivariate)   initparam$numtime=2
    if(!initparam$bivariate){namesnuis = NuisParam(model, bivariate,ncol(initparam$X),copula)}
    else                    {namesnuis =initparam$namesnuis}

if(is.null(X)) initparam$X=NULL
#######################################
    # Return the objects list:
    CovMat <- list(bivariate =  initparam$bivariate,
                   coordx = coords_orig[,1],
                   coordy = coords_orig[,2],
                   coordz = coordz,
                   coordt = initparam$coordt,
                   coordx_dyn = coordx_dyn,
                   covmatrix=covmatrix,
                   copula=copula,
                   corrmodel = corrmodel,
                   distance = distance,
                   grid=   grid,
                   nozero=initparam$setup$nozero,
                   maxdist = maxdist,
                   maxtime = maxtime,
                   n=n,
                   ns=initparam$ns,
                   NS=initparam$NS,
                   model=initparam$model,
                   namescorr = initparam$namescorr,
                   namesnuis=namesnuis,
                   namessim = initparam$namessim,
                   numblock = initparam$numblock,
                   numcoord = initparam$numcoord,
                   numcoordx = initparam$numcoordx,
                   numcoordy = initparam$numcoordy,
                   numtime = initparam$numtime,
                   param = param,
                   radius = radius,
                   setup=setup,
                   spacetime = initparam$spacetime,
                   sparse=sparse,
                   tapmod=taper,
                   tapsep=tapsep,
                   X=initparam$X)
    structure(c(CovMat, call = call), class = c("GeoCovmatrix"))
}

print.GeoCovmatrix <- function(x, digits = max(3, getOption("digits") - 3), ...)
  {
  if(x$biv) {biv="Bivariate";x$numtime=1} else {biv="Univariate"}
  if(x$numtime==1) {type="Spatial"
                   if(x$biv) type="Spatial Biviariate"
             } else {type="Spatio-temporal"}
 
  if(!x$biv) {dime=x$numcoord*x$numtime}
  else       {dime=2*x$numcoord}
  if(x$model==1) model <- 'Gaussian'
  if(x$model==21) model <- 'Gamma'
  if(x$model==39) model <- 'TwoPieceBimodal'
  if(x$model==24) model <- 'LogLogistic'
  if(x$model==30) model <- 'Poisson'
  if(x$model==46) model <- 'PoissonGamma'
  if(x$model==43) model <- 'PoissonZIP'
  if(x$model==50) model <- 'Beta2'
  if(x$model==12) model <- 'StudentT'
  if(x$model==9) model <- 'Tukeygh '
  if(x$model==18) model <- 'SkewStudentT'
  if(x$model==25) model <- 'Logistic'
  if(x$model==34) model <- 'Tukeyh'
  if(x$model==40) model <- 'Tukeyh2'
  if(x$model==23) model <- 'Gamma2'
  if(x$model==22) model <- 'LogGaussian'
  if(x$model==10) model <- 'SkewGaussian'
  if(x$model==27) model <- 'TwoPieceStudentT'
  if(x$model==38) model <- 'TwoPieceTukeyh'
  if(x$model==29) model <- 'TwoPieceGaussian'
  if(x$model==20) model <- 'SinhAsinh'    
  if(x$model==13) model <- 'Wrapped'
  if(x$model==26) model <- 'Weibull'
  if(x$model==11) model <- 'Binomial'
  if(x$model==49) model <- 'BinomialLogistic'
  if(x$model==33) model <- 'Kumaraswamy'
  if(x$model==42) model <- 'Kumaraswamy2'
  if(x$model==28) model <- 'Beta'
  if(x$model==31) model <- 'Binomial_TwoPieceGauss'
  if(x$model==32) model <- 'BinomialNeg_TwoPieceGauss'
  if(x$model==19) model <- 'Binomial2'     
  if(x$model==16) model <- 'BinomialNeg'
  if(x$model==45) model <- 'BinomialNegZINB'
  if(x$model==14) model <- 'Geometric'
  if(x$model==15) model <- 'PoisBin'
  if(x$model==17) model <- 'PoisBinNeg'  
                      
    cat('\n##################################################################')
    cat('\n Covariance matrix type:', type,'\n')
    cat('\n Dimension:',dime,'X',dime,'\n')
    cat('\n Model:', model, '\n')
    cat('\n Correlation model:', x$corrmodel, '\n')
    cat('\n Number of spatial coordinates:', x$numcoord, '\n')
    cat('\n Number of dependent temporal realisations:', x$numtime, '\n')
    cat('\n Type of the random field:', biv, '\n')
    cat('\n Number of  parameters:', length(x$param), '\n')
    cat('\nParameters:\n')
    print.default(unlist(x$param), digits = digits, print.gap = 2,
                  quote = FALSE)
    cat('\n##################################################################\n')
    invisible(x)
  }

