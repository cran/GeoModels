####################################################
### File name: GeoCorrFct.r
####################################################

GeoCorrFct<- function(x,t=NULL,corrmodel, model="Gaussian",distance="Eucl",  
                                  param, radius=6371,n=1,covariance=FALSE,variogram=FALSE)

{
  
  #############################################################################################
CorrelationFct <- function(bivariate,corrmodel, lags, lagt, numlags, numlagt, mu,model, nuisance,param,N)
    {


       if(!bivariate) { 
         nn=numlags*numlagt
         #print(nn)
         p=dotCall64::.C64('VectCorrelation',SIGNATURE = c("double","integer","double","integer","integer","double",
                  "integer","double","double","double","integer"),  
    corr=dotCall64::vector_dc("double",nn),corrmodel,lags,numlags, numlagt,mu,model,nuisance,param,lagt,N,
         INTENT =c("rw",rep("r",10)),PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)

             cc=p$corr
                    }
        else    {
                nn=numlags*4
                            p=dotCall64::.C64('VectCorrelation_biv',SIGNATURE = c("double","double","integer","double","integer","integer","double",
                  "integer","double","double","double","integer"),  
    corr=dotCall64::vector_dc("double",nn),vario=double(nn),corrmodel,lags,numlags, numlagt,mu,model,nuisance,param,lagt,N,
         INTENT =c("rw",rep("r",11)),PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)

                 cc=c(p$corr,p$vario)   

                    }
        return(cc)
    }
  #############################################################################################
  #################### end internal function ##################################################
  #############################################################################################
    call <- match.call()

   # opar=par(no.readonly = TRUE)
#on.exit(par(opar))

   # if(covariance&&variogram||!covariance&&!variogram) covariance=TRUE;variogram=FALSE
 
    model <- .GeoCovCanonicalModel(model)
    if(is.null(CkCorrModel (corrmodel))) stop("The name of the correlation model  is not correct\n")
    if(is.null(CkModel(model)))   stop("The name of the  model  is not correct\n")
    if (!is.numeric(x) || is.complex(x) || !length(x) ||
        any(!is.finite(x)) || any(x < 0))
      stop("GeoCorrFct: distances must be finite, nonnegative numeric values.", call. = FALSE)
    if (!is.null(t) && (!is.numeric(t) || is.complex(t) || !length(t) ||
                       any(!is.finite(t)) || any(t < 0)))
      stop("GeoCorrFct: temporal lags must be finite and nonnegative.", call. = FALSE)
    if(!is.list(param)) stop("param  must be a list\n")
    .GeoValidateParameterNames(param, "GeoCorrFct")

  .st_old <- .GeoSTEnter(corrmodel, distance = distance,
                         context = "GeoCorrFct")
  on.exit(.GeoSTRestore(.st_old), add = TRUE)
  .GeoSTValidate(corrmodel, param, context = "GeoCorrFct")

    spacetime<-CheckST(CkCorrModel(corrmodel))
    bivariate<-CheckBiv(CkCorrModel(corrmodel))

    supported_models <- c(
      "Gaussian","SkewGaussian","StudentT","SkewStudentT","Tukeyh","Tukeyh2",
      "Tukeygh","Kumaraswamy","Kumaraswamy2","SinhAsinh","TwoPieceStudentT",
      "TwoPieceTukeyh","TwoPieceGaussian","Gamma","Weibull","LogLogistic",
      "LogGaussian","Binary","Bernoulli","Binomial","BinomialNeg","Geom","Geometric","Poisson",
      "PoissonGamma","SkewLaplace"
    )
    if(!(model %in% supported_models))
      stop("GeoCorrFct: observed-scale correlation/covariance is not implemented for model = '",
           model, "'.", call.=FALSE)
    if(bivariate && model != "Gaussian")
      stop("GeoCorrFct: bivariate observed-scale curves are currently implemented only for model = 'Gaussian'.",
           call.=FALSE)

    if(model %in% c("Binary","Bernoulli")) {
      if(!is.numeric(n) || length(n) != 1L || !is.finite(n) || n != 1)
        stop("GeoCorrFct: Binary/Bernoulli is the one-trial field; n must equal 1.", call.=FALSE)
      n <- 1L
    }
    if(model == "Binomial") {
      if(!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 || n > .Machine$integer.max ||
         abs(n-round(n)) > sqrt(.Machine$double.eps))
        stop("GeoCorrFct: a lag-only Binomial curve requires one common positive integer n; use GeoCovmatrix for site-specific n.", call.=FALSE)
      n <- as.integer(round(n))
    }
    if(model %in% c("BinomialNeg","Geom","Geometric")) {
      if(model %in% c("Geom","Geometric")) n <- 1L
      if(!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 || n > .Machine$integer.max ||
         abs(n-round(n)) > sqrt(.Machine$double.eps))
        stop("GeoCorrFct: Negative-Binomial n is the common number r of successes and must be one positive integer scalar.", call.=FALSE)
      n <- as.integer(round(n))
    }

    if(!bivariate) {
      if(is.null(param$sill)) param$sill <- 1
      if(is.null(param$nugget)) param$nugget <- 0
    }

if(is.null(t)) t=0
   
mu=0;nuisance=0
mm=0
num_beta=c(1,1)

nx=length(x)
if(spacetime) nt=length(t)
else nt=1

num_betas=c(1,1)
if(!bivariate) {
  supplied_mean <- .GeoMean_supplied_names(param)
  if(length(setdiff(supplied_mean, "mean")))
    stop("GeoCorrFct() has no design matrix X; only the intercept parameter 'mean' is allowed.\n")
  if(is.null(param$mean)) param$mean <- 0
  if(length(param$mean) != 1L || !is.numeric(param$mean) || !is.finite(param$mean))
    stop("GeoCorrFct() requires a finite scalar intercept 'mean'.\n")
}

mu=if(!bivariate) as.numeric(param$mean) else 0
if (bivariate && !covariance) {
  if (is.null(param[["sill_1"]])) param$sill_1 <- 1
  if (is.null(param[["sill_2"]])) param$sill_2 <- 1
}
.GeoValidateCorrelationParameters(corrmodel, param, context = "GeoCorrFct")
.GeoValidateCurveNuisance(param, model, bivariate)
  ## selecting nuisance mean and correlation parameters
      if(!bivariate){
        corr_names <- CorrelationPar(CkCorrModel(corrmodel))
        parcorr <- unlist(param[corr_names], use.names=FALSE)
        storage.mode(parcorr) <- "double"

        nuis_names <- NuisParam(model,FALSE,num_betas)
        nuisance <- unlist(param[nuis_names], use.names=TRUE)
        storage.mode(nuisance) <- "double"
        sel <- names(nuisance) %in% .GeoMean_names(1L)
        mm <- as.numeric(nuisance[sel])
        nuisance <- nuisance[!sel]
        }
      if(bivariate){
        if(!covariance) param$sill_1=param$sill_2=1
        corr_names <- CorrelationPar(CkCorrModel(corrmodel))
        parcorr <- unlist(param[corr_names], use.names=FALSE)
        storage.mode(parcorr) <- "double"

        nuis_names <- NuisParam(model,bivariate,num_betas)
        nuisance <- unlist(param[nuis_names], use.names=TRUE)
        storage.mode(nuisance) <- "double"
    }

correlation <- CorrelationFct(bivariate,CkCorrModel(corrmodel), x, t, nx, nt,mu,
                                     CkModel(model), nuisance,parcorr,n)

if(is.null(t)) t=0
if(length(t)>1)  A=as.matrix(expand.grid(x,t)) 


####### Gaussian
   if(model=="Gaussian"){
        if(bivariate){  

                        cova11 <- correlation[1:nx]
                        cova12 <- correlation[(nx+1):(2*nx)]
                        cova22 <- correlation[(3*nx+1):(4*nx)]
                        variogram11  <- correlation[(4*nx+1):(5*nx)]
                        variogram12  <- correlation[(5*nx+1):(6*nx)]
                        variogram22  <- correlation[(7*nx+1):(8*nx)]
                      
                               }
        else { 
             vs=as.numeric(nuisance["sill"])

             if(length(t)>1) cova <- correlation*(1-as.numeric(nuisance["nugget"])) + as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)
             else   cova <- correlation*(1-as.numeric(nuisance["nugget"]))+as.numeric(nuisance["nugget"])*I(x==0 & t[1L]==0)
        }
}



if(!bivariate){
##########################
###### non Gaussian cases
##########################

   if(model=="SkewGaussian") {    
            if(bivariate) {}
              else {

               if(length(t)>1) correlation1=(1-as.numeric(nuisance['nugget']) )*correlation  + as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)
               else correlation1=(1-as.numeric(nuisance['nugget']) )*correlation  +as.numeric(nuisance["nugget"])*I(x==0)

              vv=as.numeric(nuisance['sill']);sk=as.numeric(nuisance['skew']);sk2=sk^2;corr2=correlation^2;  
              cc=(2*sk2)*(sqrt(1-corr2) + correlation*asin(correlation)-1)/(pi*vv+sk2*(pi-2)) + (correlation1*vv)/(vv+sk2*(1-2/pi))
              vs=(vv+sk2*(1-2/pi))
              cova=cc; 
               }
                   }
  ##########################################
   if(model=="StudentT")        { if(bivariate) {}
                        else {
                              
                 if(length(t)>1) correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)
                 else            correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(x==0)

                              nu=1/as.numeric(nuisance['df']);sill=as.numeric(nuisance['sill'])
                              vs=sill/(1-2/nu)
                              cc=.GeoStudentCorrelation(correlation, correlation1, nu)
                              
                              cova=cc;
                               }
                  }

      ##########################################
   if(model=="SkewStudentT")        { if(bivariate) {}
                        else {
                              
                 if(length(t)>1) correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)
                 else            correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(x==0)

                              nu=1/as.numeric(nuisance['df']);sill=as.numeric(nuisance['sill']);sk=as.numeric(nuisance['skew'])
                              cc <- .GeoSkewStudentCorrelation(correlation, correlation1, nu, sk)
                              vs <- sill * .GeoSkewStudentMoments(nu, sk)$variance
                              cova=cc;
                               }
                  }
##########################################

 if(model=="Tukeyh")        { if(bivariate) {}
                        else {
                             

    if(length(t)>1) correlation=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)
             else   correlation=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(x==0)

                              h=as.numeric(nuisance['tail'])
                              sill=as.numeric(nuisance['sill'])
                              vs=  sill*(1-2*h)^(-1.5)    
                              cc=(correlation*(1-2*h)^(1.5))/((1-h)^2-(h*correlation)^2)^(1.5)
                             
                              cova=cc;
                             }
                  }
   if(model=="Tukeyh2")        { if(bivariate) {}
                        else {
                if(length(t)>1) correlation=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)
                else       correlation=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(x==0)
                              hr=as.numeric(nuisance['tail1']); hl=as.numeric(nuisance['tail2'])
                              mom_hh <- .GeoTukeyh2StandardMoments(hr, hl)
                              cova=.GeoTukeyh2Correlation(correlation, hr, hl)
                              vs=as.numeric(nuisance['sill'])*mom_hh$variance
                             }
                  }   
##########################################
  if(model=="Tukeygh")       { if(bivariate) {}
                        else {
                             
                          if(length(t)>1)
               correlation=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)  
                          else     correlation=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(x==0)

                              rho=correlation
                              tail=as.numeric(nuisance['tail'])
                              sill=as.numeric(nuisance['sill'])
                              eta=as.numeric(nuisance['skew'])
                              rho2=rho*rho; eta2=eta*eta; tail2=tail*tail;

                        if(tail>1e-05&&abs(eta)>1e-05){
                               rho2=rho*rho; eta2=eta*eta; tail2=tail*tail;
                               u=1-tail; a=1+rho;
                               mu=(exp(eta2/(2*u))-1)/(eta*sqrt(u));
                               var0=(exp(2*eta2/(1-2*tail))-2*exp(eta2/(2*(1-2*tail)))+1)/(eta2*sqrt(1-2*tail))-mu*mu
                               vs=sill*var0
                               A1=exp(a*eta2/(1-tail*a));
                               A2=2*exp(0.5*eta2*  (1-tail*(1-rho2))  / (u*u- tail2*rho2)  );
                               A3=eta2*sqrt(u*u- rho2*tail*tail);
                               cc=((A1-A2+1)/A3-mu*mu)/var0
                               
                               cova=cc;
                            } 
                        if(tail<=1e-05&&abs(eta)>1e-05){
                              var0=( -exp(eta^2)+exp(eta^2*2))*eta^(-2)
                              vs=sill*var0
                              cc=(( -exp(eta^2)+exp(eta^2*(1+rho)))*eta^(-2))/var0
                              
                              cova=cc;
                            } 
                        if(tail>1e-05&&abs(eta)<=1e-05){
                              ## When skew is approximately zero, Tukey-gh reduces to Tukey-h.
                              ## This avoids the undefined object `h` previously used here and
                              ## keeps the formula consistent with the Tukeyh block above.
                              vs=sill*(1-2*tail)^(-1.5)
                              cc=(rho*(1-2*tail)^(1.5))/((1-tail)^2-(tail*rho)^2)^(1.5)
                              
                              cova=cc;
                            } 
                        if(tail<=1e-05&&abs(eta)<=1e-05){
                            cova=rho;
                            vs=sill
                            }  
                }}


if(model=="Kumaraswamy"||model=="Kumaraswamy2")  { if(bivariate) {}
else {


if(length(t)>1)
               correlation=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0) 
else correlation=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(x==0)

if(model=="Kumaraswamy"){
ga=as.numeric(nuisance['shape2'])
eta=as.numeric(nuisance['shape1'])
}
if(model=="Kumaraswamy2"){
ga=as.numeric(nuisance['shape'])
pmed=plogis(mu)
eta=log(0.5)/log1p(-(pmed^ga))
}
mm=eta*beta(1+1/ga,eta)
ssup=as.numeric(nuisance['max'])-as.numeric(nuisance['min'])
sill=ssup^2*(eta*beta(1+2/ga,eta)-mm^2)
###
NN=length(correlation);
res=double(NN)
bb=.C("corr_kuma_vec",as.double(correlation),as.double(eta),as.double(ga), 
     res=as.double(res),
    as.integer(NN),PACKAGE='GeoModels', NAOK=TRUE)
rho=bb$res
cova=rho;
vs=sill
 }
}

##########################


if(model == "SinhAsinh") {
  if(bivariate) {} else {
    nugget <- as.numeric(nuisance['nugget'])
    d <- as.numeric(nuisance['tail'])
    e <- as.numeric(nuisance['skew'])
    sill <- as.numeric(nuisance['sill'])
    if(d <= 0 || !is.finite(d) || !is.finite(e))
      stop("GeoCorrFct: SinhAsinh requires finite skew and a positive finite tail.", call.=FALSE)
    if(length(t)>1)
      rho_lat <- correlation*(1-nugget) + nugget*I(A[,1]==0&A[,2]==0)
    else
      rho_lat <- correlation*(1-nugget) + nugget*I(x==0)
    mom_sas <- .GeoSASStandardMoments(e, d)
    cova <- corrsas(rho_lat,e,d)
    vs <- sill*mom_sas$variance
  }
}

##########################################
 if(model=="TwoPieceStudentT")        { if(bivariate) {}
                        else {

                            if(length(t)>1)
                                   correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)
                              else correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(x==0)

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
                              cova=cc;
                               }
                  } 
  ##########################################
  if(model=="TwoPieceTukeyh")        { if(bivariate) {}
                        else {
                              

                              if(length(t)>1) correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)
                              else correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(x==0)


                              tail=as.numeric(nuisance['tail']); sk=nuisance['skew'];sill=nuisance['sill']
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
                              cova=cc;
                               }
                  }  
##########################################
 if(model=="TwoPieceGaussian")        { 
                        if(bivariate) {}
                        else {        

                         if(length(t)>1) correlation1=correlation*(1-as.numeric(nuisance['nugget']))+as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)
                        else correlation1=correlation*(1-as.numeric(nuisance['nugget']))+as.numeric(nuisance["nugget"])*I(x==0)


                        corr2=sqrt(1-correlation1^(2))
                        sk=as.numeric(nuisance['skew']); sk2=sk^2
                        ll=qnorm((1-sk)/2)
                        p11=pbivnorm::pbivnorm(ll,ll, rho = correlation, recycle = TRUE)
                        KK=3*sk2+2*sk+ 4*p11 - 1
                        cc=(2*((corr2 + correlation1*asin(correlation1))*KK)- 8*sk2)/(3*pi*sk2  -  8*sk2   +pi   )
                        vs= as.numeric(nuisance['sill'])*(1+3*sk2-8*sk2/pi)
                        cova=cc;
                        }
                  } 
    ##########################################
 if(model=="SkewLaplace") { if(bivariate) {}
                        else {
                          if(length(t)>1) rho_lat=correlation*(1-as.numeric(nuisance['nugget']))+as.numeric(nuisance['nugget'])*I(A[,1]==0&A[,2]==0)
                          else rho_lat=correlation*(1-as.numeric(nuisance['nugget']))+as.numeric(nuisance['nugget'])*I(x==0)
                          sk=as.numeric(nuisance['skew'])
                          cova=rho_lat^2
                          vs=as.numeric(nuisance['sill'])*(1/sk^2+1/(1-sk)^2)
                        }
                  }
    ##########################################
 if(model=="Gamma")        { if(bivariate) {}
                        else {
                              
                          if(length(t)>1) correlation=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)
                          else   correlation=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(x==0)



                              vs=2*exp(mm)^2/as.numeric(nuisance['shape'])
                              cc=correlation^2
                              cova=cc;
                            }
                  }
##########################################
if(model=="Weibull")        { if(bivariate) {} 
                        else {
                         

                        if(length(t)>1)   correlation=correlation*(1-as.numeric(nuisance['nugget'] )  )+as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)
                        else correlation=correlation*(1-as.numeric(nuisance['nugget'] )  )+as.numeric(nuisance["nugget"])*I(x==0)


                        sh=as.numeric(nuisance["shape"])
                        vs=exp(mm)^2*(gamma(1+2/sh)/gamma(1+1/sh)^2-1)
                        cc=.GeoWeibullCorr(correlation, sh)
                        cova=cc;
                     }
                    }
##########################################
if(model=="LogLogistic")    { if(bivariate) {}  
                      else { 


                     if(length(t)>1)  correlation=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)
                     else correlation=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(x==0)


                     sh=as.numeric(nuisance["shape"])
                     vs=exp(mm)^2*(2*sh*sin(pi/sh)^2/(pi*sin(2*pi/sh))-1)
                     cc=((pi*sin(2*pi/sh))/(2*sh*(sin(pi/sh))^2-pi*sin(2*pi/sh)))*
                                    (Re(hypergeo::hypergeo(-1/sh, -1/sh, 1,correlation^2))*
                                     Re(hypergeo::hypergeo( 1/sh,  1/sh, 1,correlation^2)) -1)
                      cova=cc;
                    }
                    }
if(model=="LogGaussian")    { if(bivariate) {}  
                      else {    
                    
                    if(length(t)>1)  correlation=correlation*(1-as.numeric(nuisance["nugget"]))  +as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)
                    else  correlation=correlation*(1-as.numeric(nuisance["nugget"]))   +as.numeric(nuisance["nugget"])*I(x==0)


                    vvar=as.numeric(nuisance["sill"])
                    yy=vvar*correlation
                    cc<-(exp(yy)-1)/(exp(vvar)-1)
                    vs<-exp(2*mm)*(exp(vvar)-1)
                    cova=cc; 
                     }
                  }
#############################################
if(model %in% c("Binary","Bernoulli","Binomial","BinomialNeg","Geom","Geometric"))
   {
                    if(!bivariate) {
                           pp=pnorm(mu)
                           if(model %in% c("Binary","Bernoulli")) vs=pp*(1-pp)
                           if(model=="Binomial") vs=as.numeric(n)*pp*(1-pp)
                           if(model=="BinomialNeg") vs=as.numeric(n)*(1-pp)/pp^2
                           if(model %in% c("Geom","Geometric")) vs=(1-pp)/pp^2
                           ## VectCorrelation returns the observed-scale
                           ## correlation.  Do not rely on numerical evaluation
                           ## of the bivariate Gaussian cdf at latent rho = 1:
                           ## that is a singular boundary point.  Enforce the
                           ## exact same-site endpoint directly, so that
                           ## correlation(0)=1, covariance(0)=Var(Y), and
                           ## variogram(0)=0 by construction.
                           cova=correlation
                           if(spacetime) {
                             if(length(t)>1L) zero_lag <- (A[,1] == 0 & A[,2] == 0)
                             else zero_lag <- (x == 0 & t[1L] == 0)
                           } else zero_lag <- (x == 0)
                           cova[zero_lag] <- 1
                           }
   }
if(model=="Poisson") {
                    if(bivariate) {}
                    if(!bivariate) {   
                           if(length(t)>1) correlation=(1-as.numeric(nuisance['nugget']))*correlation  +as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)
                           else correlation=(1-as.numeric(nuisance['nugget']))*correlation  +as.numeric(nuisance["nugget"])*I(x==0)
                           corr2=correlation^2    
                           vs=exp(mu);
                           z=2*vs/(1-corr2)
                           cc=corr2*(1-(besselI(z,0,expon.scaled = TRUE)+besselI(z,1,expon.scaled = TRUE)))
                           cova=cc
                         }
                      }
if(model=="PoissonGamma") {
        if(bivariate) {}
        if(!bivariate) {
          ## PEM / Poisson-Erlang mixture correlation.  Corollary 5.1 is an
          ## off-diagonal expression and contains terms involving
          ## 1-rho^2.  At the same space-time point rho=1 is a singular
          ## numerical boundary for the hypergeometric representation, while
          ## the theoretical endpoint is exactly Corr{M(s),M(s)}=1.  Evaluate
          ## the closed form only away from zero lag and impose the endpoint
          ## directly.
          if(length(t)>1L) {
            zero_lag <- (A[,1] == 0 & A[,2] == 0)
          } else {
            zero_lag <- (x == 0 & t[1L] == 0)
          }

          nug <- as.numeric(nuisance['nugget'])
          latent_corr <- (1-nug)*correlation + nug*zero_lag
          a <- as.numeric(nuisance['shape'])
          mean_pg <- exp(mu)

          ## Use the compiled stationary PEM kernel shared with covariance
          ## matrices and kriging.  Integer kappa values use terminating 2F1
          ## polynomials; half-integers/general values use the compiled
          ## evaluator rather than R-level integration/continued fractions.
          cc <- .Call("GeoCorrPGsVecCall",
                      as.double(latent_corr), as.double(mean_pg), as.double(a),
                      PACKAGE="GeoModels")
          cc[zero_lag] <- 1

          ## Var{M(s)} = mu * (1 + mu/kappa), with kappa = shape.
          vv <- mean_pg*(1+mean_pg/a)
          vs <- vv
          cova <- cc
        }
    }


  } ## end not bivariate                 
#################################################################

if(!bivariate){
   if(!covariance) vs=1
   res=cova*vs; 
   if(variogram) res=vs*(1-cova)
 }
if(bivariate)
{
 if(!variogram) res=rbind(cova11,cova12,cova22)
 else  res=rbind(variogram11,variogram12,variogram22)
}

GeoCorrFct <- list(corr=res,
                   distances=x,
                    times=t,
                    model=model,
                    distance=distance,  
                    param=param,
                    radius=radius,
                    n=n,
                    covariance=covariance,
                    variogram=variogram,
                    spacetime=spacetime,
                    bivariate=bivariate)

structure(c(GeoCorrFct, call = call), class = c("GeoCorrFct"))




}




plot.GeoCorrFct<- function(x,type="p",...)
  {

if(!inherits(x,"GeoCorrFct"))       stop("Enter an object obtained from the function GeoCorrFct\n")

opar=par(no.readonly = TRUE)
on.exit(par(opar))

#old.par <- par(no.readonly = TRUE)

space=!x$bivariate&&!x$spacetime
################################# spatial ###################################################################
if(space) {
    if(type=="p") { plot(x$distances,x$corr,xlab="Distances",type="p",axes=FALSE,...);axis(1);axis(2)}
    if(type=="l") { 
                   if(as.numeric(x$param$nugget)>0){maxcor=max(x$corr)
                                                    plot.new(); plot.window(xlim=c(0,max(x$distances)),ylim= c(0,maxcor))
                                                    if(!x$variogram) points(0,maxcor,pch=20) else points(0,0,pch=20)
                                                    sel=x$distances>0
                                                    lines(x$distances[sel],x$corr[sel],...)
                                                    axis(1);axis(2)}
                   else  {   plot(x$distances,x$corr,type="l",axes=FALSE,...)        ;axis(1);axis(2)}                    
                  } 
   }

if(x$bivariate){
##################################  bivariate ##################################################################
  if(type=="p") {
par(mfrow=c(1,3))
 plot(x$distances,x$corr[1,],xlab="Distances",main="1st component",axes=FALSE,type="p",...);axis(1);axis(2)
 plot(x$distances,x$corr[2,],xlab="Distances",main="Cross component",axes=FALSE,type="p",...);axis(1);axis(2)
 plot(x$distances,x$corr[3,],xlab="Distances",main="2nd component",axes=FALSE,type="p",...);axis(1);axis(2)
}
  if(type=="l") {
  par(mfrow=c(1,3))  
###########
          if(as.numeric(x$param$nugget_1)>0){
                                                    maxcor=max(x$corr[1,])
                                                    plot.new(); plot.window(xlim=c(0,max(x$distances)),ylim= c(0,maxcor));title(main="1st component")
                                                    if(!x$variogram) points(0,maxcor,pch=20) else points(0,0,pch=20)
                                                    sel=x$distances>0
                                                    lines(x$distances[sel],x$corr[1,][sel],...)
                                                    axis(1);axis(2)
          }
          else   { plot(x$distances,x$corr[1,],xlab="Distances",main="1st component",axes=FALSE,type="l",...);axis(1);axis(2) }
##########
          if(as.numeric(x$param$nugget_1)>0||as.numeric(x$param$nugget_2)>0) {

                                                    maxcor=max(x$corr[2,])
                                                    mincor=min(x$corr[2,])
                                                    plot.new(); plot.window(xlim=c(0,max(x$distances)),ylim= c(mincor,maxcor));title(main="Cross component")
                                                    if(!x$variogram) {
                                                    bb=as.numeric(x$param$pcol*sqrt( (x$param$nugget_1+x$param$sill_1)*(x$param$nugget_2+x$param$sill_2)))
                                                    points(0,bb,pch=20) }
                                                    else points(0,0,pch=20)
                                                    sel=x$distances>0
                                                    lines(x$distances[sel],(x$corr[2,])[sel],...)
                                                    axis(1);axis(2)
          }
          else  plot(x$distances,x$corr[2,],xlab="Distances",main="Cross component",axes=FALSE,type="l",...);axis(1);axis(2)
##########
          if(as.numeric(x$param$nugget_2)>0){

                                                    maxcor=max(x$corr[3,])
                                                    plot.new(); plot.window(xlim=c(0,max(x$distances)),ylim= c(0,maxcor));title(main="2nd component")
                                                    if(!x$variogram) points(0,maxcor,pch=20) else points(0,0,pch=20)
                                                    sel=x$distances>0
                                                    lines(x$distances[sel],x$corr[3,][sel],...)
                                                    axis(1);axis(2)
          }
          else   { plot(x$distances,x$corr[3,],xlab="Distances",main="2nd component",axes=FALSE,type="l",...);axis(1);axis(2) }
  }
}
####################################################################################################
if(x$spacetime){
par(mfrow=c(1,3))
cc=matrix(x$corr,nrow=length(x$distances),ncol=length(x$times))

persp(cc,x= x$distances,y=x$times, theta = 20, phi = 30, 
     ticktype = "detailed",zlab="",xlab="Distance",ylab="Time")
# A=as.matrix(expand.grid(x$distances,x$times)) 
#scatterplot3d::scatterplot3d(A[,1],A[,2], c(cc),
#                              type="h",highlight.3d=TRUE,cex.axis=.7,cex.lab=.7,
#                              main="",
#                              ,xlab="Distance",ylab="Time")


 if(type=="p") {  plot( x$distances,cc[,1],type="p",xlab="Spatial Distances",main="Spatial Marginal",axes=FALSE,...)    ;axis(1);axis(2)
               
                  plot( x$times,cc[1,],    type="p",xlab="Times",main="Temporal marginal",axes=FALSE,...)    ;axis(1);axis(2)
              }

 if(type=="l") { 

    if(as.numeric(x$param$nugget)>0){

        maxcor=max(cc[1,])
        plot.new(); plot.window(xlim=c(0,max(x$distances)),ylim= c(0,maxcor));title(main="Spatial Marginal")
        if(!x$variogram) points(0,maxcor,pch=20) else points(0,0,pch=20)
        sel=x$distances>0
        lines(x$distances[sel],(cc[,1])[sel],...)
         axis(1);axis(2)

        maxcor=max(cc[,1])
        plot.new(); plot.window(xlim=c(0,max(x$times)),ylim= c(0,maxcor));title(main="Temporal marginal")
        if(!x$variogram) points(0,maxcor,pch=20) else points(0,0,pch=20)
        sel=x$times>0
        lines(x$times[sel],(cc[1,])[sel],...)
        axis(1);axis(2)
     }
    else { plot( x$distances,cc[,1],type="l",xlab="Spatial Distances",main="Spatial Marginal",axes=FALSE,...)    ;axis(1);axis(2)
            plot( x$times,cc[1,],    type="l",xlab="Times",main="Temporal marginal",axes=FALSE,...)    ;axis(1);axis(2)
         }   
   }
}
################################################################################################
invisible()
}  
