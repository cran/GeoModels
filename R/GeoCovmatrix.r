

GeoCovmatrix <- function(estobj=NULL,coordx=NULL, coordy=NULL, coordz=NULL, coordt=NULL,coordx_dyn=NULL,corrmodel, distance="Eucl", grid=FALSE,
                     model="Gaussian", n=1, param, anisopars=NULL,radius=1,
                       sparse=FALSE,copula=NULL,X=NULL,spobj=NULL,check.duplicates=FALSE)

{
  call <- match.call()

  ########################################################################################################
  ##########  Internal function: computing covariance matrix #############################################
  ########################################################################################################




    Cmatrix <- function(bivariate, coordx, coordy,coordz, coordt,corrmodel, dime, n, ns, NS, nuisance, numpairs,
                           numpairstot, model, paramcorr, setup, radius, spacetime, spacetime_dyn,type,copula,ML,other_nuis)
    {
###################################################################################
############### computing correlation #############################################
###################################################################################

if(model %in% c(1,9,34,12,20,18,39,27,38,29,21,26,24,10,22,40,28,33,42,59,13))
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
     ## Sparse univariate kernels include the diagonal.  Its exact marginal
     ## variance is installed below, so skip singular transformation endpoints
     ## there only.  Bivariate kernels already return covariances, and a value
     ## of one can be a genuine variance or cross-covariance.
     if (!bivariate) {
       sparse_rows <- rep.int(seq_len(dime), diff(setup$ia))
       cr$corr[setup$ja == sparse_rows] <- 0
     }
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
              aa=(1-2*h)^(1.5) # reciprocal standardized variance
              corr <- aa*corr/( (1-h)^2-h^2*corr^2 )^(1.5)
              vv=ss/aa
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
    corr_lat=cr$corr*(1-as.numeric(nuisance['nugget']))
    hr=as.numeric(nuisance['tail1'])
    hl=as.numeric(nuisance['tail2'])
    mom_hh <- .GeoTukeyh2StandardMoments(hr, hl)
    corr <- .GeoTukeyh2Correlation(corr_lat, hr, hl)
    vv=as.numeric(nuisance['sill'])*mom_hh$variance
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
   rho_base=cr$corr
   rho_num=rho_base*(1-as.numeric(nuisance['nugget']))
   nu=as.numeric(1/nuisance['df'])
   # In the constructive Student-t RF the common scale-mixture component
   # depends on the nugget-free latent correlation, while the Gaussian
   # numerator carries the nugget attenuation.
   corr=.GeoStudentCorrelation(rho_base, rho_num, nu)
   vv=as.numeric(nuisance['sill'])/(1-2/nu)
}
if(bivariate){}
}

###############################################################
  if(model==18)   ##  skew student case
    {
     if(!bivariate)
{
    rho_num <- cr$corr * (1 - as.numeric(nuisance['nugget']))
    nu <- as.numeric(1 / nuisance['df'])
    sk <- as.numeric(nuisance['skew'])
    corr <- .GeoSkewStudentCorrelation(cr$corr, rho_num, nu, sk)
    vv <- as.numeric(nuisance['sill']) * .GeoSkewStudentMoments(nu, sk)$variance
}
if(bivariate){}
}
###############################################################

if (model == 20) {  ## SAS
  nugget <- as.numeric(nuisance['nugget'])
  d <- as.numeric(nuisance['tail'])
  e <- as.numeric(nuisance['skew'])
  sill <- as.numeric(nuisance['sill'])
  mom_sas <- .GeoSASStandardMoments(e, d)
  corr_base <- cr$corr * (1 - nugget)
  corr <- corrsas(corr_base, e, d)
  vv <- sill * mom_sas$variance
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
if(model==59)  {  ##  skewLaplace case

          if(!bivariate){
              corr=cr$corr*(1-as.numeric(nuisance['nugget']))
              sk=as.numeric(nuisance['skew'])
              corr=corr^2
              vv=as.numeric(nuisance['sill'])*(1/sk^2 + 1/(1-sk)^2)
     }
      if(bivariate){}
}
#################################################################################
################ covariance matrix for models defined on the real line ##########
#################################################################################
if(model %in% c(1,9,34,12,20,18,39,27,38,29,10,40,59)){
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
  res=0;ss=0;k=0;A=rep(0,length(corr22))
  while(k<=100){
    p1=2*(lgamma(cc+k)-lgamma(cc)+lgamma(nu2+1+k)-lgamma(nu2+1))
    p2=lgamma(k+1)+(lgamma(nu2+k)-lgamma(nu2))+2*(lgamma(cc+1+k)-lgamma(cc+1))
    b1=p1-p2
    b2=log(hypergeo::genhypergeo(U=c(cc+k,cc+k,alpha2), L=c(cc+k+1,cc+k+1), polynomial=TRUE,maxiter=1000, z=corr22))
    b3=k*log(corr22)
    sum=exp(b1+b2+b3)
    res=res+sum
    A=res
    if (all(sum<1e-6)){
      break
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
                ga=as.numeric(nuisance['shape'])
                mu_kuma2=as.numeric(ML)
                if(length(mu_kuma2) != 1L && (max(mu_kuma2)-min(mu_kuma2)) > sqrt(.Machine$double.eps))
                  stop("GeoCovmatrix: plain Kumaraswamy2 covariance is currently implemented only for a constant mean; use the copula covariance engine for location-dependent margins.", call.=FALSE)
                pmed=plogis(mu_kuma2[1L])
                eta=log(0.5)/log1p(-(pmed^ga))
                   }

      mm=eta*beta(1+1/ga,eta)
      NN=length(corr);
      res=double(NN)
      bb=.C("corr_kuma_vec",as.double(corr),as.double(eta),as.double(ga),
      res=as.double(res),as.integer(NN),PACKAGE='GeoModels', NAOK=TRUE)
      corr=bb$res
      ssup=as.numeric(nuisance['max'])-as.numeric(nuisance['min'])
      vv=ssup^2*(eta*beta(1+2/ga,eta)-mm^2)
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


if(model==13)   ## wrapped Gaussian
    {
      if(!bivariate) {
         nugget_wrap <- as.numeric(nuisance['nugget'])
         sill_wrap <- as.numeric(nuisance['sill'])
         if(!is.finite(sill_wrap) || sill_wrap <= 0 ||
            !is.finite(nugget_wrap) || nugget_wrap < 0 || nugget_wrap >= 1)
           stop("GeoCovmatrix: Wrapped requires sill > 0 and 0 <= nugget < 1.",
                call. = FALSE)

         ## Off-diagonal correlation of the latent Gaussian residuals.
         rho_wrap <- cr$corr * (1 - nugget_wrap)

         ## Jammalamadaka--Sarma circular correlation.  This is exactly
         ## sinh(sill * rho) / sinh(sill), evaluated by a stable formula.
         corr_wrap <- .GeoWrappedCircularCorrelation(rho_wrap, sill_wrap)

         ## GeoCovmatrix cannot use ordinary covariance of numerical angles,
         ## which depends on the arbitrary branch cut.  For Wrapped it returns
         ## the covariance of sin(Theta-mu), the component underlying the
         ## circular correlation above.  The cosine-component covariance is
         ## attached and exposed by the outer GeoCovmatrix object.
         cov_sin <- .GeoWrappedSinCovariance(rho_wrap, sill_wrap)
         cov_cos <- .GeoWrappedCosCovariance(rho_wrap, sill_wrap)
         vv_sin <- .GeoWrappedSinVariance(sill_wrap)
         vv_cos <- .GeoWrappedCosVariance(sill_wrap)

         if(type=="Standard") {
           varcov <- diag(vv_sin, dime)
           varcov[lower.tri(varcov)] <- cov_sin
           varcov <- t(varcov)
           varcov[lower.tri(varcov)] <- cov_sin

           coscov <- diag(vv_cos, dime)
           coscov[lower.tri(coscov)] <- cov_cos
           coscov <- t(coscov)
           coscov[lower.tri(coscov)] <- cov_cos

           circorr <- diag(1, dime)
           circorr[lower.tri(circorr)] <- corr_wrap
           circorr <- t(circorr)
           circorr[lower.tri(circorr)] <- corr_wrap
         }
         if(type=="Tapering") {
           varcov <- new("spam", entries=cov_sin, colindices=setup$ja,
                         rowpointers=setup$ia, dimension=as.integer(rep(dime,2)))
           diag(varcov) <- vv_sin
           coscov <- new("spam", entries=cov_cos, colindices=setup$ja,
                         rowpointers=setup$ia, dimension=as.integer(rep(dime,2)))
           diag(coscov) <- vv_cos
           circorr <- new("spam", entries=corr_wrap, colindices=setup$ja,
                          rowpointers=setup$ia, dimension=as.integer(rep(dime,2)))
           diag(circorr) <- 1
         }
         attr(varcov, "wrapped_cos_covmatrix") <- coscov
         attr(varcov, "wrapped_circular_correlation") <- circorr
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
         ## Use Eq. (6) directly.  The previous Euler-transformed form
         ## involved 2F1(1+1/sh,1+1/sh;1;rho^2), which is numerically
         ## singular for the hypergeo backend near rho^2=1.
         corr=.GeoWeibullCorr(corr, sh)
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
if(model==14) n=rep.int(1L,dime)
if(model %in% c(16,45)) n=rep.int(as.integer(n[1L]),dime)
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

  if(model %in% c(57,58))   { mm=exp(mu); kk=as.numeric(param$shape);
                            pg=pnorm(param$pmu)
                           vv=(1-pg)*mm*(1+pg*mm+(mm/kk))
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
                         else coordx_dyn=estobj$coordx_dyn
                      } ## spatial (temporal) non regular case
 else  {    if(is.null(estobj$coordx_dyn))  {
                                              n1 <- as.integer(estobj$ns[1L])
                                              idx1 <- seq_len(n1)
                                              coordx=estobj$coordx[idx1]    # bivariate, common static locations
                                              coordy=estobj$coordy[idx1]
                                              coordz=if(is.null(estobj$coordz)) NULL else estobj$coordz[idx1]
                                            }  
            else {coordx_dyn=estobj$coordx_dyn}                                      # bivariate  dynamic  
       }
 }
else  { coordx=estobj$coordx; 
        coordy=estobj$coordy
        coordz=estobj$coordz
      }
   coordt=estobj$coordt
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
model <- .GeoCovCanonicalModel(model)
if( !is.character(corrmodel)|| is.null(CkCorrModel(corrmodel)))       stop("the name of the correlation model is wrong")
 if(is.null(CkModel(model))) stop("The name of the  model  is not correct\n")


 compact_models <- c(
  "GenWend_Matern","GenWend_Matern2","Wend0", "wend0","Wend1", "wend1","Wend2", "wend2","Genwend", "GenWend",
  "Hypergeometric2","Hypergeometric","Hypergeometric_Matern","Hypergeometric_Matern2","GenWend_Hole", "GenWend_hole",
  "GenWend_Matern_Hole", "GenWend_Matern_hole","Wen0_space", "wen0_space","Wen0_time",  "wen0_time",
  "Wen1_space", "wen1_space","Wen1_time",  "wen1_time","Wen2_space", "wen2_space","Wen2_time",  "wen2_time",
  "Gneiting_wen_S", "gneiting_wen_S", "Gneiting_Wen_S","Gneiting_wen_T", "gneiting_wen_T","Genwend_Genwend_nosep", "GenWend_GenWend_nosep",
  "Wend0_Wend0","Wend0_Wend1","Wend0_Wend2","Wend1_Wend0","Wend1_Wend1","Wend1_Wend2","Wend2_Wend0","Wend2_Wend1","Wend2_Wend2",
  "GenWend_GenWend","Bi_wend0_sep", "Bi_Wend0_contr","Bi_wend0","Bi_wend1_sep","Bi_Wend1_contr","Bi_wend1",
  "Bi_wend2_sep","Bi_wend2","Bi_Wend2_contr","Bi_GenWend_sep", "Bi_genWend_sep","Bi_GenWend_contr", "Bi_genWend_contr", "Bi_GenWend", "Bi_genWend"
)

if (sparse && !(corrmodel %in% compact_models)) {
  stop("Sparse option should be used only for compactly supported covariance models\n")
}

bivariate<-CheckBiv(CkCorrModel(corrmodel))
spacetime<-CheckST(CkCorrModel(corrmodel))
space=!spacetime&&!bivariate

## Normalize the optional copula field before input checking.
if(!is.null(copula)) {
  copula <- as.character(copula)[1L]
  if(is.na(copula) || !nzchar(trimws(copula))) copula <- NULL
}
if(!is.null(copula) && !(copula %in% c("Gaussian", "Clayton", "SkewGaussian")))
  stop("GeoCovmatrix supports copula='Gaussian', 'Clayton', or 'SkewGaussian'.",
       call. = FALSE)
if(!is.null(copula) && copula %in% c("Clayton", "SkewGaussian") &&
   !(as.integer(CkModel(model)) %in% .GeoCopulaPredictionModels))
  stop("GeoCovmatrix supports copula='", copula,
       "' only for the validated copula prediction margins.", call. = FALSE)

if(is.null(copula)) {
  plain_supported <- c(
    "Gaussian","Tukeygh","Tukeyh","Tukeyh2","StudentT","SinhAsinh",
    "SkewStudentT","TwoPieceBimodal","TwoPieceStudentT","TwoPieceTukeyh",
    "TwoPieceGaussian","SkewGaussian","SkewLaplace","Beta","Kumaraswamy",
    "Kumaraswamy2","Wrapped","Gamma","Weibull","LogLogistic","LogGaussian",
    "Binary","Bernoulli","Binomial","Geometric","Geom","BinomialNeg",
    "Poisson","PoissonZIP","BinomialNegZINB","PoissonGamma",
    "PoissonGammaZIP","PoissonGammaZIP1"
  )
  if(!(model %in% plain_supported))
    stop("GeoCovmatrix: plain observed-scale covariance is not implemented for model = '",
         model, "'.", call.=FALSE)
  if(bivariate && model != "Gaussian")
    stop("GeoCovmatrix: plain bivariate observed-scale covariance is currently implemented only for model = 'Gaussian'.",
         call.=FALSE)
  if(model %in% c("PoissonGammaZIP","PoissonGammaZIP1") && (spacetime || sparse))
    stop("GeoCovmatrix: PoissonGammaZIP/PoissonGammaZIP1 covariance is currently implemented only for dense purely spatial designs.",
         call.=FALSE)
}

Mobs_cov <- NULL

tapsep=NULL ;   maxdist=NULL; maxtime=NULL; type="Standard"    #setted for older compatibility
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

## Preserve the user's grid declaration in the returned object.  Internally a
## regular grid is expanded to explicit locations and then processed as an
## irregular coordinate matrix.
grid_input <- isTRUE(grid)

.GeoValidate3D(
  coordx = coordx, coordy = coordy, coordz = coordz,
  coordx_dyn = coordx_dyn, distance = distance, grid = grid,
  context = "GeoCovmatrix"
)

  .st_old <- .GeoSTEnter(corrmodel,
    coordx = coordx, coordy = coordy,
    coordz = coordz, coordx_dyn = coordx_dyn,
    distance = distance, context = "GeoCovmatrix")
  on.exit(.GeoSTRestore(.st_old), add = TRUE)
  .GeoSTValidate(corrmodel, param, context = "GeoCovmatrix")

## Hole-effect generalized Wendland validity depends on the Euclidean
## dimension through the turning-bands construction.  Infer d from the
## coordinate representation before flattening it.
hole_d <- 2L
if(!is.null(coordz)) hole_d <- 3L
if(is.matrix(coordx) && ncol(coordx) >= 3L) hole_d <- 3L
if(is.list(coordx_dyn) && length(coordx_dyn) &&
   any(vapply(coordx_dyn, function(z) is.matrix(z) && ncol(z) >= 3L, logical(1L))))
  hole_d <- 3L
.GeoValidateHoleWendlandParameters(corrmodel, param, d=hole_d, context="GeoCovmatrix")

if(!is.null(coordx_dyn)) {
  coordx_dyn <- .GeoNormalizeDynamicCoordinates(
    coordx_dyn,
    coordt = if(spacetime) coordt else NULL,
    context = "GeoCovmatrix"
  )
  if(spacetime) {
    .GeoValidateDynamicRows(X, coordx_dyn, name = "X",
                            context = "GeoCovmatrix")
  }
}
if(bivariate) {
  biv_info <- .GeoBivariateCoordinateInfo(
    coordx = coordx, coordy = coordy, coordz = coordz,
    coordx_dyn = coordx_dyn, grid = grid, context = "GeoCovmatrix"
  )
  if(biv_info$dynamic) coordx_dyn <- biv_info$blocks
  if(!is.null(X)) {
    X <- .GeoBivariateDesign(X, biv_info$ns, context = "GeoCovmatrix")
  }
}

.GeoCheckDuplicateLocations(
  coordx = coordx, coordy = coordy, coordz = coordz, coordt = coordt,
  coordx_dyn = coordx_dyn, spacetime = spacetime, bivariate = bivariate,
  grid = grid, context = "GeoCovmatrix", check.duplicates = check.duplicates
)

#######################################################
     ## setting zero mean and nugget if no mean or nugget is fixed
    if(!bivariate){
    if(is.null(param$mean)) param$mean<-0
    ## CkInput(..., fcall = "Simulation") validates the union of the
    ## Gaussian nuisance parameters and the declared model nuisance
    ## parameters.  Hence the auxiliary generic nugget must exist here even
    ## for double-nugget ZIP models.  It is removed *after* StartParam(),
    ## before building the positional nuisance vector passed to
    ## CorrelationMat_dis2 (see cleanup below).
    if(is.null(param$nugget)) param$nugget<-0
    #if(length(param$mean)>1) MM=fixed$mean
    if(is.null(param$sill)) param$sill<-1  
}
    else{
    if(is.null(param$mean_1)) param$mean_1<-0
    if(is.null(param$mean_2)) param$mean_2<-0
    if(is.null(param$nugget_1)) param$nugget_1<-0
    if(is.null(param$nugget_2)) param$nugget_2<-0
    sill_names <- c("sill_1", "sill_2")
    missing_sills <- sill_names[vapply(sill_names, function(nm)
      is.null(param[[nm]]), logical(1L))]
    if(length(missing_sills))
      stop("GeoCovmatrix: bivariate models require parameter(s) ",
           paste(missing_sills, collapse=", "), ".", call.=FALSE)

}

    if(is.null(coordx_dyn)){
    coordx=unname(coordx);coordy=unname(coordy);coordz=unname(coordz)}
    else{coordx_dyn=lapply(coordx_dyn, unname)}


    #if the covariance is compact supported  and option sparse is used
    #then set the code as a tapering and an object spam is returned
if(sparse) {
    covmod=CkCorrModel(corrmodel)
    if(covmod %in% c(10,11,13,15,19,6,7,21,22,23,26,29,30,
                     63,64,65,66,67,68,
                     69,70,71,72,73,74,75,76,77,78,85,
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
        ## Reparameterized hole-effect GW is GW_{nu*a,xi-1/2,nu,d,k};
        ## hence its exact compact support is nu*a, not a.
        if(covmod==29) maxdist=as.numeric(param$power2*param$scale)
      }
  }
  taper=corrmodel
}
##### sill parameter fixed for some models

if(!bivariate){

if(model %in% c("Weibull","Poisson","Binomial","Gamma","LogLogistic",
        "BinomialNeg","Binary","Bernoulli","Geom","Geometric","Gaussian_misp_Poisson",
        'PoissonZIP','Gaussian_misp_PoissonZIP','BinomialNegZINB', 'PoissonGammaZIP','PoissonGammaZIP1','PoissonGamma',
        'PoissonZIP1','Gaussian_misp_PoissonZIP1','BinomialNegZINB1',
        'Beta2','Kumaraswamy2','Beta','Kumaraswamy')){
if(is.null(param$sill)) param$sill=1
else                    param$sill=1
}
}


coords=NULL
coordx_dyn_orig=coordx_dyn

if(is.null(coordx_dyn)){
 #######################################
    if(grid) { coords=as.matrix(expand.grid(coordx,coordy)); grid=FALSE }
    else     {  if(!is.null(coordy)) { coords=as.matrix(cbind(coordx,coordy,coordz))}
                else                 { coords=as.matrix(coordx);}
                if(ncol(coords)==3) coordz=coords[,3]
                else                coordz=NULL
             }
    coords_orig=coords
    if(!is.null(anisopars)) {  coords=GeoAniso(coords,c(anisopars$angle,anisopars$ratio))}
}
else {
    coords_orig=do.call(rbind, coordx_dyn_orig)
    coordx_dyn_metric=coordx_dyn_orig
    if(!is.null(anisopars)) {
      if(ncol(coords_orig)!=2L)
        stop("GeoCovmatrix: anisotropy for dynamic three-dimensional coordinates is not implemented.", call.=FALSE)
      coordx_dyn_metric=lapply(coordx_dyn_orig, function(z)
        GeoAniso(z,c(anisopars$angle,anisopars$ratio)))
    }
    coordx_dyn=coordx_dyn_metric
    coords=do.call(rbind, coordx_dyn)
    if(ncol(coords)==3) coordz=coords[,3]
    else                coordz=NULL
}
coordz_orig=if(ncol(coords_orig)==3L) coords_orig[,3] else NULL
 #######################################
#######################################

    if(!bivariate) {
      nobs_mean_cov <- if(!is.null(coordx_dyn)) {
        sum(vapply(coordx_dyn, function(z) nrow(as.matrix(z)), integer(1L)))
      } else {
        nrow(coords) * if(is.null(coordt)) 1L else length(coordt)
      }
      if(!is.numeric(n) || !length(n) || any(!is.finite(n)))
        stop("GeoCovmatrix: n must be finite numeric.", call.=FALSE)
      if(any(n < 1) || any(abs(n-round(n)) > sqrt(.Machine$double.eps)))
        stop("GeoCovmatrix: n must contain positive integer values.", call.=FALSE)
      if(is.null(copula) && model %in% c("Binary","Bernoulli")) {
        if(length(n) != 1L || as.integer(round(n[1L])) != 1L)
          stop("GeoCovmatrix: Binary/Bernoulli is the one-trial field; n must equal 1.", call.=FALSE)
        n <- 1L
      }
      if(is.null(copula) && model %in% c("Geom","Geometric")) {
        if(length(n) != 1L || as.integer(round(n[1L])) != 1L)
          stop("GeoCovmatrix: Geometric is the r = 1 Negative-Binomial field; n must equal 1.", call.=FALSE)
        n <- 1L
      }
      if(is.null(copula) && model %in% c("BinomialNeg","BinomialNegZINB")) {
        if(length(n) != 1L)
          stop("GeoCovmatrix: for direct Negative-Binomial fields n is the common number r of successes and must be a positive integer scalar.", call.=FALSE)
        n <- as.integer(round(n[1L]))
      } else if(model == "Binomial" &&
                length(n) != 1L && length(n) != nobs_mean_cov) {
        stop("GeoCovmatrix: Binomial n must be scalar or contain exactly one value per observation (",
             nobs_mean_cov, ").", call.=FALSE)
      }
      Mobs_cov <- .GeoMean_external(param, n = nobs_mean_cov)
      if(!is.null(Mobs_cov) && !is.null(X))
        stop("A site-specific mean vector and X cannot be used together.\n")
      if(!is.null(X)) {
        X <- .GeoMean_design(X, n = nobs_mean_cov, name = "X")
        invisible(.GeoMean_beta(param, ncol(X)))
      } else {
        supplied_mean <- .GeoMean_supplied_names(param)
        if(is.null(Mobs_cov) && length(setdiff(supplied_mean, "mean")))
          stop("Multiple mean coefficients require an explicit design matrix X.\n")
        if(is.null(Mobs_cov)) invisible(.GeoMean_beta(param, 1L))
      }
      if(!is.null(Mobs_cov)) param$mean <- 0
    }


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
                           type, FALSE,copula,X,FALSE)

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

    if(initparam$model %in% c(43,45,57)) initparam$namesnuis=initparam$namesnuis[!initparam$namesnuis %in% "nugget"]



other_nuis=NULL
ML=0
##### sill parameter fixed for some models
if(!bivariate){

if(model %in% c("Weibull","Poisson","Binomial","Gamma","LogLogistic",
        "BinomialNeg","Binary","Bernoulli","Geom","Geometric","Gaussian_misp_Poisson",
        'PoissonZIP','Gaussian_misp_PoissonZIP','BinomialNegZINB', 'PoissonGammaZIP','PoissonGammaZIP1','PoissonGamma',
        'PoissonZIP1','Gaussian_misp_PoissonZIP1','BinomialNegZINB1',"LogGaussian",
        'Beta2','Kumaraswamy2','Beta','Kumaraswamy')) 
  {
  ### setting mean because mean affects the covariance in this case
  parc <- initparam$param[initparam$namescorr]
  inip <- initparam$param[!(names(initparam$param) %in% names(parc))]
  mean_names <- .GeoMean_names(ncol(initparam$X))
  other_nuis2 <- inip[!(names(inip) %in% mean_names)]
  ## Build the ZIP nuisance vectors in the exact positional order expected by
  ## CorrelationMat_dis2.  StartParam() deliberately carries an auxiliary
  ## Gaussian `nugget` for Simulation-input validation; for double-nugget ZIP
  ## models that auxiliary entry must never reach the C kernel.
  zip_order <- switch(as.character(initparam$model),
    "43" = c("nugget1", "nugget2", "pmu"),
    "44" = c("nugget1", "nugget2", "pmu"),
    "45" = c("nugget1", "nugget2", "pmu"),
    "57" = c("nugget1", "nugget2", "pmu", "shape"),
    "58" = c("nugget", "pmu", "shape"),
    NULL)
  if(!is.null(zip_order)) {
    missing_zip <- setdiff(zip_order, names(other_nuis2))
    if(length(missing_zip))
      stop("GeoCovmatrix: internal ZIP nuisance construction is missing ",
           paste(missing_zip, collapse=", "), ".", call.=FALSE)
    other_nuis <- as.numeric(other_nuis2[zip_order])
  } else {
    other_nuis <- as.numeric(other_nuis2[order(names(other_nuis2))])
  }
  if(!is.null(Mobs_cov)) {
    ML <- Mobs_cov
  } else {
    beta <- .GeoMean_beta(initparam$param, ncol(initparam$X))
    ML <- as.numeric(initparam$X %*% beta)
  }
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
    covmatrix<- Cmatrix_copula(initparam$bivariate,cc[,1],cc[,2],cc[,3],initparam$coordt,initparam$corrmodel,dime,n,initparam$ns,
                        initparam$NS, initparam$param[initparam$namesnuis],
                        initparam$numpairs,numpairstot,initparam$model,
                        initparam$param[initparam$namescorr],setup,initparam$radius,initparam$spacetime,spacetime_dyn,initparam$type,copula,ML,other_nuis)
                    }

    wrapped_cos_covmatrix <- attr(covmatrix, "wrapped_cos_covmatrix", exact = TRUE)
    wrapped_circular_correlation <- attr(covmatrix, "wrapped_circular_correlation", exact = TRUE)

    if(type=="Tapering") sparse=TRUE

    # Delete the global variables:
    if(!space)
                      .C('DeleteGlobalVar', PACKAGE='GeoModels', NAOK=TRUE)
    else                     
                      .C('DeleteGlobalVar2', PACKAGE='GeoModels', NAOK=TRUE)
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
                   coordz = coordz_orig,
                   coordt = initparam$coordt,
                   coordx_dyn = coordx_dyn_orig,
                   covmatrix=covmatrix,
                   copula=copula,
                   corrmodel = corrmodel,
                   distance = distance,
                   grid=   grid_input,
                   nozero=initparam$setup$nozero,
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
                   ## Kept separate because param must remain a scalar parameter
                   ## vector for covariance and kriging internals.
                   mean_external = if(!initparam$bivariate) Mobs_cov else NULL,
                   radius = radius,
                   setup=setup,
                   spacetime = initparam$spacetime,
                   sparse=sparse,
    
                   X=initparam$X)
    if(initparam$model == 13) {
      CovMat$wrapped_cos_covmatrix <- wrapped_cos_covmatrix
      CovMat$wrapped_circular_correlation <- wrapped_circular_correlation
    }
    CovMat$call <- call
    structure(CovMat, class = c("GeoCovmatrix"))
}

print.GeoCovmatrix <- function(x, digits = max(3, getOption("digits") - 3), ...)
  {
  is_bivariate <- isTRUE(x$bivariate)
  if(is_bivariate) biv <- "Bivariate" else biv <- "Univariate"

  if(is_bivariate) type <- "Spatial Bivariate"
  else if(isTRUE(x$spacetime)) type <- "Spatio-temporal"
  else type <- "Spatial"

  dd <- dim(x$covmatrix)
  dime <- if(length(dd)) dd[1L] else length(x$covmatrix)

  model <- paste0("model code ", x$model)
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
  if(x$model==59) model <- 'SkewLaplace'
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
    if(!is.null(x$copula)) cat('\n Copula:', x$copula, '\n')
    cat('\n Correlation model:', x$corrmodel, '\n')
    if(!is.null(x$coordx_dyn)) {
      cat('\n Spatial block sizes:', paste(x$ns, collapse = ', '), '\n')
    } else {
      cat('\n Number of spatial coordinates:', x$numcoord, '\n')
    }
    if(isTRUE(x$spacetime))
      cat('\n Number of dependent temporal realisations:', length(x$coordt), '\n')
    cat('\n Type of the random field:', biv, '\n')
    cat('\n Number of parameters:', length(x$param), '\n')
    cat('\nParameters:\n')
    print.default(unlist(x$param), digits = digits, print.gap = 2,
                  quote = FALSE)
    cat('\n##################################################################\n')
    invisible(x)
  }
