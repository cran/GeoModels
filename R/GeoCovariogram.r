####################################################
### File name: GeoCovVariogram.r
####################################################

### Compute and plot the (estimated) covariance function and the variogram
### from  GeoFit object
GeoCovariogram <- function(fitted, distance="Eucl", answer.cov=FALSE, answer.vario=FALSE,
                         answer.range=FALSE, fix.lags=NULL, fix.lagt=NULL,
                        show.cov=FALSE, show.vario=FALSE,  show.range=FALSE,
                        add.cov=FALSE, add.vario=FALSE, pract.range=95,
                        vario=NULL,invisible=FALSE,...)
  {
    result <- NULL

    CheckDistance<- function(distance)
    {
        CheckDistance <- NULL
        CheckDistance <- switch(distance,
                                eucl=0,
                                Eucl=0,
                                chor=1,  
                                Chor=1,
                                geod=2,
                                Geod=2,
                                proj=3,
                                Proj=3)
        return(CheckDistance)
    }
    # define the bivariate Gaussian distribution:
    vpbnorm <- function(corrmodel,lags,lagt,nuisance,numlags,numlagt,param,threshold)
    {
        p=.C("vpbnorm",  as.integer(corrmodel), as.double(lags), as.double(lagt),
           as.integer(numlags), as.integer(numlagt), as.double(nuisance),
           as.double(param), rho=double(numlags*numlagt), as.double(threshold), PACKAGE='GeoModels',
           NAOK=TRUE)
        return(p$rho)
    }            
    


CorrelationFct <- function(bivariate, corr_model, lags, lag_t, num_lags, num_lag_t, mu, model, nuisance, param, N) {
  if (!bivariate) {
    result <- dotCall64::.C64(
      "VectCorrelation",
      SIGNATURE = c("double", "integer", "double", "integer", "integer", "double", "integer", "double", "double", "double", "integer"),
      INTENT = c("w", rep("r", 10)),
      corr = dotCall64::vector_dc("double", num_lags * num_lag_t),
      corr_model = corr_model,
      lags = lags,num_lags = num_lags,num_lag_t = num_lag_t,mu = mu,
      model = model,nuisance = nuisance,param = param,lag_t = lag_t,N = N,PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
    
    cc <- result$corr
  } else {
    result <- dotCall64::.C64(
      "VectCorrelation_biv",
      SIGNATURE = c("double", "double", "integer", "double", "integer", "integer", "double", "integer", "double", "double", "double", "integer"),
      INTENT = c("w", "w", rep("r", 10)),
      corr = dotCall64::vector_dc("double", num_lags * 4),vario = dotCall64::vector_dc("double", num_lags * 4),corr_model = corr_model,
      lags = lags,num_lags = num_lags,num_lag_t = num_lag_t,mu = mu,
      model = model,nuisance = nuisance,param = param,lag_t = lag_t,N = N,PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
    
    cc <- c(result$corr, result$vario)
  }
  cc
}




if(show.range)  {
    # Pratical range in the Gaussian case:
    PracRangeNorm <- function(corrmodel, lags, lagt, nuisance, numlags, numlagt, param, pract.range)
    { 
        return(nuisance["sill"]*(1-nuisance["nugget"])*CorrelationFct(bivariate,corrmodel, lags, lagt, numlags, numlagt, mu,CkModel(fitted$model),nuisance, param,fitted$n)
            -(nuisance["sill"]*(1-pract.range)))
    
    }
}
###########################################################     
#############  starting ###################################
########################################################### 

    dyn=FALSE
    OLS=NULL
    isvario <- !is.null(vario) # is empirical variogram is passed?
    bivariate <- fitted$bivariate
    if(bivariate) fitted$numtime=1
    ispatim <- fitted$spacetime
    dyn<- is.list(fitted$coordx_dyn)

    
    if(!isvario) stop("an object vario is needed\n")
    if(!is.null(fitted$copula)) stop("copula models are not supported \n")

if (fitted$model %in% c("Weibull", "Poisson","PoissonGamma", "Binomial", "Gamma", 
        "LogLogistic", "BinomialNeg", "Bernoulli", "Geometric", 
        "Gaussian_misp_Poisson", "PoissonZIP", "Gaussian_misp_PoissonZIP", 
        "Gaussian_misp_PoissonGamma","PoissonGammaZIP","PoissonGammaZIP1",
        "BinomialNegZINB", "PoissonZIP1", "Gaussian_misp_PoissonZIP1", 
        "BinomialNegZINB1", "Beta2", "Kumaraswamy2", "Beta", 
        "Kumaraswamy")) { fitted$param$sill=1}


opar=par(no.readonly = TRUE)
on.exit(par(opar))

if(!ispatim && !bivariate){ if( (show.cov && show.vario) || (show.cov)) par(mfrow=c(1,2))}
if(show.vario && ispatim) par(mfrow=c(1,2))
if(show.vario && ispatim && !is.null(fix.lags) && !is.null(fix.lagt)) par(mfrow=c(2,2))
if(show.cov && ispatim && !is.null(fix.lags) && !is.null(fix.lagt)) par(mfrow=c(2,2))
if(show.cov && ispatim && is.null(fix.lags) && is.null(fix.lagt)) par(mfrow=c(1,1))
if(show.vario && bivariate) {par(mfrow=c(2,2))}
if(bivariate&&dyn) par(mfrow=c(1,2))

fitted$param=unlist(fitted$param)
fitted$fixed=unlist(fitted$fixed)
    
    # START ---- check input --- #

    if( !inherits(fitted,"GeoFit")  & !inherits(fitted,"GeoWLS") )
        stop("Enter an object obtained GeoFit or WLeastSquare")
  
    if(isvario & !inherits(vario,"GeoVariogram") )
        stop("Enter an object obtained from the function GeoVariogram")

    if(!is.numeric(pract.range) & answer.range)
        stop("Enter a number for the parameter % of sill")
    else{
        if(pract.range < 0 || pract.range > 100)
            stop("Entered an incorrect value for the % of sill")
        else
          pract.range <- pract.range / 100}
    # Set the type of process:
    model <- CkModel(fitted$model)
    gaussian <- model==1
    skewgausssian<- model==10
    gamma<- model==21
    skewstudentT<- model==18||model==37
    studentT<- model==12||model==35
    weibull<- model==26
    twopieceT<- model==27
    twopieceTukeyh<- model==38 
    twopieceGauss<- model==29
    twopiecebimodal<- model==39
    loggauss<- model==22
    binary <- model==2
    binomial <- model==11
    binomial2 <- model==19
    geom <- model==14
    binomialneg<-model==16
    binomialnegZINB<-model==45
    tukeyh<- model ==34
    tukeyh2<- model ==40
    sas<- model ==20
    poisson<- model==30
    poissongamma<- model==46
    poissonZIP<- model==43||model==44
    poissongammaZIP<- model==57
    loglogistic <- model==24
    tukeygh<- model==9||model==41
    Gaussian_misp_Binomial<-model==51
    Gaussian_misp_Poisson<-model==36
    Gaussian_misp_PoissonGamma<- model==47
    Gaussian_misp_BinomialNeg<-model==52
    zero <- 0;slow=1e-3;
    if(gaussian||skewgausssian||gamma||loggauss||binomial||Gaussian_misp_Binomial||Gaussian_misp_BinomialNeg||Gaussian_misp_Poisson||binomialneg||binomialnegZINB||geom||tukeyh||tukeyh2||sas||twopiecebimodal||skewstudentT
            ||twopieceGauss||twopieceTukeyh||twopieceT) slow=1e-6
    else slow=1e-5 
    # lags associated to empirical variogram estimation
    if(isvario){
    lags <- c(0,vario$centers);numlags <- length(lags)
    if(ispatim) {
     if(is.null(vario$centert)) lagt <-c(0,vario$bint) 
     else  lagt <-c(0,vario$centert)
    }
    else {lagt=0}
    numlagt <- length(lagt)
    }

    # check the fixed lags:
    if(!is.null(fix.lags)){
        if(!is.numeric(fix.lags) || fix.lags<0 || fix.lags>numlags){
            cat("invalid spatial fixed lag\n")
            return(result)}}
    if(!is.null(fix.lagt))
        if(!is.numeric(fix.lagt) || fix.lagt<0 || fix.lagt>numlagt){
            cat("invalid temporal fixed lag\n")
            return(result)}
    # END ---- check input --- #
    # set range (for practical range)
    if(fitted$corrmodel=="matern"){
        lower <- 1e-10
        upper <- 1e20}
    else{
        lower <- 0  
        upper <- 1e100}

    num_betas=fitted$numbetas
    mu=0;nuisance=0
    mm=0
    ## selecting nuisance mean annd corr parameters
      if(!bivariate){
       param <- c(fitted$fixed,fitted$param)[CorrelationPar(CkCorrModel(fitted$corrmodel))]
        nuisance <- c(fitted$fixed,fitted$param)[NuisParam2(fitted$model,FALSE,num_betas)]
        sel=substr(names(nuisance),1,4)=="mean"
        mm=nuisance[sel]
        nuisance=nuisance[!sel]
        }
      if(bivariate){
        param <- c(fitted$fixed,fitted$param)[CorrelationPar(CkCorrModel(fitted$corrmodel))]
        nuisance <- c(fitted$fixed,fitted$param)[NuisParam2(fitted$model,FALSE,num_betas)]
    }
    
     #############################################
     # computing the spatio-temporal distances where to compute the fitted model
     #############################################
    if(isvario) {

             if(!invisible) {   #### standard case
                              lags_m <- seq(slow,max(vario$centers),length.out =150)
                              if(ispatim)  lagt_m <-seq(slow,max(vario$bint),length.out =150)
                              else         lagt_m<-0
                            }
             else           {   ### invisible case
                                if (ispatim) {lags_m=c(slow,vario$centers)
                                              lagt_m <-c(slow,vario$bint)}
                                else         {lags_m=vario$centers
                                             lagt_m<-0}
                            }    
     }
    else{
        mmm <- double(2)
        type_dist <- CheckDistance(distance)
        p=.C("Maxima_Minima_dist",mmm=as.double(mmm),as.double(fitted$coordx),as.double(fitted$coordy),as.double(fitted$coordz),
        as.integer(fitted$numcoord),as.integer(type_dist),as.double(fitted$radius),
        PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
        lags_m <- seq(slow,p$mmm[2],length=150)
        if (ispatim) {
            tt <- double(2)
            lagt_m <- seq(slow,max(c(dist(fitted$coordt))) ,length=150)
        }
        else lagt_m<-0
    }
    numlags_m <- length(lags_m)
    numlagt_m <- length(lagt_m)
if(!bivariate) {
    if(is.null(fitted$X))  X=as.matrix(rep(1,numlags_m*numlagt_m))
    mu=X%*%as.numeric(mm)}
  else{}  
   ###############
    corrmodel <- CkCorrModel(fitted$corrmodel)

     nui=nuisance
     nui['sill']=1;
if(!(binomial||geom||binomialneg||binomialnegZINB||Gaussian_misp_Binomial||
       Gaussian_misp_BinomialNeg||Gaussian_misp_Poisson||Gaussian_misp_PoissonGamma)) nui['nugget']=0
else                                     nui['nugget']=nuisance['nugget']

    if(fitted$model=="Gaussian_misp_Binomial") fitted$model="Binomial"
    if(fitted$model=="Gaussian_misp_BinomialNeg") fitted$model="BinomialNeg"
    if(fitted$model=="Gaussian_misp_Poisson") fitted$model="Poisson"
    if(fitted$model=="Gaussian_misp_PoissonGamma") fitted$model="PoissonGamma"
  
    correlation <- CorrelationFct(bivariate,corrmodel, lags_m, lagt_m, numlags_m, numlagt_m,mu,
                                     CkModel(fitted$model), nui,param,fitted$n)
        #correlation <- CorrelationFct(bivariate,corrmodel, lags_m, lagt_m, numlags_m, numlagt_m,mu,
         #                            CkModel(fitted$model), nui,param,fitted$n)



##########################################
########### starting cases ###############
##########################################
    if(gaussian){
    llm=length(lags_m)
        if(bivariate){  covariance11 <- correlation[(0*llm+1):(1*llm)]
                        covariance12 <- correlation[(1*llm+1):(2*llm)]
                        covariance22 <- correlation[(3*llm+1):(4*llm)]
                        variogram11  <- correlation[(4*llm+1):(5*llm)]
                        variogram12  <- correlation[(5*llm+1):(6*llm)]
                        variogram22  <- correlation[(7*llm+1):(8*llm)]
                           }
        else { 

        #covariance <- nuisance["nugget"]+nuisance["sill"]*correlation
        covariance <- as.numeric(nuisance["sill"])*correlation*(1-as.numeric(nuisance["nugget"]))
        #variogram <- nuisance["nugget"]+nuisance["sill"]*(1-correlation)
        variogram <-as.numeric(nuisance["sill"])*(1-correlation*(1-as.numeric(nuisance["nugget"])))
       
        }
    }
##########################################
 if(skewgausssian) {    
            if(bivariate) {}
              else {
              correlation1=(1-as.numeric(nuisance['nugget']) )*correlation  
              vv=as.numeric(nuisance['sill']);sk=as.numeric(nuisance['skew']);sk2=sk^2;corr2=correlation^2;  
              cc=(2*sk2)*(sqrt(1-corr2) + correlation*asin(correlation)-1)/(pi*vv+sk2*(pi-2)) + (correlation1*vv)/(vv+sk2*(1-2/pi))
              vs=(vv+sk2*(1-2/pi))
          
              covariance=vs*cc;variogram=vs*(1-cc) }
                   }
##########################################
   if(twopieceT)        { if(bivariate) {}
                        else {
                              correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))
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
                              covariance=vs*cc;variogram=vs*(1-cc)  }
                  } 
  ##########################################
   if(twopieceTukeyh)        { if(bivariate) {}
                        else {
                              correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))
                              tail=as.numeric(nuisance['tail']); sk=as.numeric(nuisance['skew']);
                              sill=as.numeric(nuisance['sill'])
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
                              covariance=vs*cc;variogram=vs*(1-cc) 
                               }
                  }  
##########################################
   if(twopieceGauss)        { 
                        if(bivariate) {}
                        else {        
                        correlation1=correlation*(1-as.numeric(nuisance['nugget']))
                        corr2=sqrt(1-correlation1^(2))
                        sk=as.numeric(nuisance['skew']); sk2=sk^2
                        ll=qnorm((1-sk)/2)
                        p11=pbivnorm::pbivnorm(ll,ll, rho = correlation, recycle = TRUE)
                        KK=3*sk2+2*sk+ 4*p11 - 1
                        cc=(2*((corr2 + correlation1*asin(correlation1))*KK)- 8*sk2)/(3*pi*sk2  -  8*sk2   +pi   )
                        vs= as.numeric(nuisance['sill'])*(1+3*sk2-8*sk2/pi)
                        covariance=vs*cc;variogram=vs*(1-cc) 
                        }
                  } 
##########################################
 if(twopiecebimodal)        { if(bivariate) {}
                        else {                            
                                  correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))
                                  nu=as.numeric(nuisance['df']); sk=as.numeric(nuisance['skew'])
                                  delta=as.numeric(nuisance['shape'])
                                  alpha=2*(delta+1)/nu
                                  nn=2^(1-alpha/2)
                                  ll=qnorm((1-sk)/2)
                                  p11=pbivnorm::pbivnorm(ll,ll, rho = correlation, recycle = TRUE)
                                  corr2=correlation1^2;sk2=sk^2
                                  a1=Re(hypergeo::hypergeo(-1/alpha ,-1/alpha,nu/2,corr2))
                                  a3=3*sk2 + 2*sk + 4*p11 - 1
                                  MM=(2^(2/alpha)*(gamma(nu/2 + 1/alpha))^2) 
                                  vari=2^(2/alpha)*(gamma(nu/2 + 2/alpha))*gamma(nu/2)* (1+3*sk2) - sk2*2^(2/alpha+2)*gamma(nu/2+1/alpha)^2 
                                  cc= MM*(a1*a3-4*sk2)/vari
                                  vs=as.numeric(nuisance['sill'])*vari/(nn^(2/alpha)*gamma(nu/2)^2)
                                  covariance=vs*cc; variogram=vs*(1-cc)
                               }
                            }
##########################################
   if(studentT)        { if(bivariate) {}
                        else {
                              correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))
                              nu=1/as.numeric(nuisance['df']);sill=as.numeric(nuisance['sill'])
      
                              vs=sill*(nu)/(nu-2)
                              cc=((nu-2)*gamma((nu-1)/2)^2*Re(hypergeo::hypergeo(0.5,0.5 ,nu/2 ,correlation^2))*correlation1)/(2*gamma(nu/2)^2)
                              covariance=vs*cc;variogram=vs*(1-cc) 
                               }
                  }
##########################################  
   if(skewstudentT)        { if(bivariate) {}
                        else {
                               correlation1=correlation*(1-nuisance['nugget'] )
                               nu=as.numeric(1/nuisance['df']); sk=as.numeric(nuisance['skew'])
                               sill=as.numeric(nuisance['sill'])
                               skew2=sk*sk;l=nu/2; f=(nu-1)/2; w=sqrt(1-skew2);y=correlation;
                               CorSkew=(2*skew2/(pi*w*w+skew2*(pi-2)))*(sqrt(1-y*y)+y*asin(y)-1)+w*w*correlation1/(w*w+skew2*(1-2/pi)) ;
                               cc=(pi*(nu-2)*gamma(f)^2/(2*(pi*gamma(l)^2-skew2*(nu-2)*gamma(f)^2)))*(Re(hypergeo::hypergeo(0.5,0.5,l,y*y))
                                *((1-2*skew2/pi)*CorSkew+2*skew2/pi)-2*skew2/pi);
                               mm=sqrt(nu)*gamma(f)*sk/(sqrt(pi)*gamma(l));
                               vs=sill*(nu/(nu-2)-mm*mm);
                               covariance=vs*cc;variogram=vs*(1-cc)
                               }
                  }
##########################################
  if(tukeyh)        { if(bivariate) {}
                        else {
                              correlation=correlation*(1-nuisance['nugget'] )
                              h=as.numeric(nuisance['tail'])
                              sill=as.numeric(nuisance['sill'])
                              vs=  (1-2*h)^(-1.5)     
                              cc=correlation*(1-2*h)^(1.5)/ ((1-h)^2-(h*correlation)^2)^(1.5)
                              covariance=sill*vs*cc;variogram=sill*vs*(1-cc)  
                             } 
                  } 
    if(tukeyh2)        { if(bivariate) {}
                        else {
                              correlation=correlation*(1-nuisance['nugget'] )
                              hr=as.numeric(nuisance['tail1']); hl=as.numeric(nuisance['tail2'])
                              sill=as.numeric(nuisance['sill'])
                              corr=correlation
                              corr[corr>=0.99999999]=0.99999999
                              x1=1-(1-corr^2)*hr; x2=(1-hr)^2-(corr*hr)^2
                              y1=1-(1-corr^2)*hl; y2=(1-hl)^2-(corr*hl)^2
                              g=1-hl-hr+(1-corr^2)*hl*hr
                              h1=sqrt(1-corr^2/(x1^2))+(corr/x1)*asin(corr/x1);
                              h2=sqrt(1-corr^2/(y1^2))+(corr/y1)*asin(corr/y1)
                              
                              h3=sqrt(1-corr^2/(x1*y1))+sqrt(corr^2/(x1*y1))*asin(sqrt(corr^2/(x1*y1)))
                              p1=x1*h1/(2*pi*(x2)^(3/2))+corr/(4*(x2)^(3/2)); p2=y1*h2/(2*pi*(y2)^(3/2))+corr/(4*(y2)^(3/2))
                              p3=-(x1*y1)^(1/2)*h3/(2*pi*(g)^(3/2))+corr/(4*(g)^(3/2))
                              mm=(hr-hl)/(sqrt(2*pi)*(1-hl)*(1-hr))
                              vs=0.5*((1-2*hl)^(-3/2)+(1-2*hr)^(-3/2)) -(mm)^2
                              cc=(p1+p2+2*p3-mm^2)/vs
                          covariance=sill*vs*cc;variogram=sill*vs*(1-cc)  
                             } 
                  }   
##########################################
  if(tukeygh)        { if(bivariate) {}
                        else {
                              correlation=correlation*(1-nuisance['nugget'] )
                              rho=correlation
                              tail=as.numeric(nuisance['tail'])
                              sill=as.numeric(nuisance['sill'])
                              eta=as.numeric(nuisance['skew'])
                              rho2=rho*rho; eta2=eta*eta; tail2=tail*tail;

                        if(tail>1e-05&&abs(eta)>1e-05){
                               rho2=rho*rho; eta2=eta*eta; tail2=tail*tail;
                               u=1-tail; a=1+rho;
                               mu=(exp(eta2/(2*u))-1)/(eta*sqrt(u));
                               vs=(exp(2*eta2/(1-2*tail))-2*exp(eta2/(2*(1-2*tail)))+1)/(eta2*sqrt(1-2*tail))-mu*mu;
                               A1=exp(a*eta2/(1-tail*a));
                               A2=2*exp(0.5*eta2*  (1-tail*(1-rho2))  / (u*u- tail2*rho2)  );
                               A3=eta2*sqrt(u*u- rho2*tail*tail);
                               cc=((A1-A2+1)/A3-mu*mu)/vs;
                               covariance=sill*vs*cc;variogram=sill*vs*(1-cc)
                            } 
                        if(tail<=1e-05&&abs(eta)>1e-05){
                              vs=( -exp(eta^2)+exp(eta^2*2))*eta^(-2)
                              cc= (( -exp(eta^2)+exp(eta^2*(1+rho)))*eta^(-2))/vs
                              covariance=sill*vs*cc;variogram=sill*vs*(1-cc) 
                            } 
                        if(tail>1e-05&&abs(eta)<=1e-05){
                              vs=  (1-2*tail)^(-1.5)     
                              cc=(-rho/((1+h*(tail-1))*(-1+tail+tail*rho)*(1+tail*(-2+tail-tail*rho^2))^0.5))/vs
                              covariance=sill*vs*cc;variogram=sill*vs*(1-cc) 
                            } 
                        if(tail<=1e-05&&abs(eta)<=1e-05){
                            covariance=sill*rho;variogram=sill*(1-rho)
                            }  
                                
                      } 
                  }     
  if(sas) { if(bivariate) {}
                        else {
         
  nugget <- as.numeric(nuisance['nugget'])
  d <- as.numeric(nuisance['tail'])     # deve essere > 0
  e <- as.numeric(nuisance['skew'])     # può essere qualsiasi valore
  sill <- as.numeric(nuisance['sill'])
  correlation=correlation*(1-as.numeric(nuisance['nugget'] ))

mm=sinh(e/d)*exp(0.25)*(besselK(.25,(d+1)/(2*d))+besselK(.25,(1-d)/(2*d)))/(sqrt(8*pi))
  vv=cosh(2*e/d)*exp(0.25)*(besselK(.25,(d+2)/(2*d))+besselK(0.25,(2-d)/(2*d)))/(sqrt(32*pi))-0.5-mm^2
  corr=corrsas(correlation, e, d)
covariance=sill*vv*corr;variogram=sill*vv*(1-corr)  
                        }
          }
##########################################
 if(gamma)        { if(bivariate) {}
                        else {
                              correlation=correlation*(1-as.numeric(nuisance['nugget'] ))
                              vs=2*exp(mm)^2/as.numeric(nuisance['shape'])
                              cc=correlation^2
                              covariance=vs*cc;variogram=vs*(1-cc)  }
                  }
##########################################
 if(weibull)        { if(bivariate) {} 
                        else {
                            ssh=as.numeric(nuisance['shape'])
                        correlation=correlation*(1-as.numeric(nuisance['nugget'] ))  
                        vs=exp(mm)^2*(gamma(1+2/ssh)/gamma(1+1/ssh)^2-1)
                        auxcorr= (gamma(1+1/ssh))^2/((gamma(1+2/ssh))-(gamma(1+1/ssh))^2)
                        cc=auxcorr*(Re(hypergeo::hypergeo(-1/ssh, -1/ssh, 1,correlation^2)) -1)
                        covariance=vs*cc;variogram=vs*(1-cc) 
                      #    print(cc) 
                    }
                    }
##########################################
  if(loglogistic)    { if(bivariate) {}  
                      else { 
                     correlation=correlation*(1-as.numeric(nuisance['nugget'] ))
                     sh=as.numeric(nuisance["shape"])
                     vs=exp(mm)^2*(2*sh*sin(pi/sh)^2/(pi*sin(2*pi/sh))-1)
                     cc=((pi*sin(2*pi/sh))/(2*sh*(sin(pi/sh))^2-pi*sin(2*pi/sh)))*
                                    (Re(hypergeo::hypergeo(-1/sh, -1/sh, 1,correlation^2))*
                                     Re(hypergeo::hypergeo( 1/sh,  1/sh, 1,correlation^2)) -1)
                      covariance=vs*cc;variogram=vs*(1-cc)   }
                    }
##########################################
 # if(loggauss)    { if(bivariate) {}  
  #                    else {    
   #                 correlation=(1-nuisance['nugget'] )*correlation
    #                vvar=as.numeric(nuisance["sill"])
     #               yy=nuisance["sill"]*correlation
      #              cc<-(exp(yy-1))
       #             vs<-(exp(vvar)-1)*exp(2*mm)   # ok
        #            covariance=vs*cc;variogram=vs*(1-cc/((exp(vvar)-1)*exp(2*vvar) ))   
         #            }
          #        }
  if(loggauss)    { if(bivariate) {}  
                      else {    
                    vvar=as.numeric(nuisance["sill"])
                    yy=vvar*correlation*(1-as.numeric(nuisance["nugget"]))
                    cc<-(exp(yy)-1)/(exp(vvar)-1)
                    vs<-(exp(vvar)-1)
                    covariance=vs*cc;variogram=vs*(1-cc)   
                     }
                  }
   #                   
   if(binary||binomial||binomial2||geom||binomialneg||binomialnegZINB||Gaussian_misp_Binomial||Gaussian_misp_BinomialNeg) {
                    if(bivariate) {}
                    if(!bivariate) {      
                            
                           pp=pnorm(mu)
                           if(binary||binomial||binomial2) vv=min(fitted$n)*pp*(1-pp)
                           if(geom)             vv=(1-pp)/pp^2;
                           if(binomialneg)      vv=(fitted$n)*(1-pp)/pp^2;
                           if(binomialnegZINB) { 
                                     pg=pnorm(as.numeric(nuisance['pmu']))
                                     vv=fitted$n*(1-pp)*(1-pg)*(1+fitted$n*pg*(1-pp)) /pp^2
                                               }
                           if(Gaussian_misp_Binomial||Gaussian_misp_BinomialNeg) vv=1                    
                           covariance=vv*correlation
                           variogram=vv*(1-correlation)
                           }
    }
     if(poisson||Gaussian_misp_Poisson) {
                    if(bivariate) {}
                    if(!bivariate) {   
                           correlation=(1-as.numeric(nuisance['nugget']))*correlation   
                           corr2=correlation^2    
                           vv=exp(mu);
                           z=2*vv/(1-corr2)
                           cc=corr2*(1-(besselI(z,0,expon.scaled = TRUE)+besselI(z,1,expon.scaled = TRUE)))
                           if(Gaussian_misp_Poisson) vv=1   
                           covariance=vv*cc
                           variogram=vv*(1-cc)}
                   }
   if(poissongamma||Gaussian_misp_PoissonGamma) {

                    if(bivariate) {}
                    if(!bivariate) {
                       correlation=(1-as.numeric(nuisance['nugget']))*correlation
                       corr2=correlation^2
                       a=as.numeric(nuisance['shape'])
                       b=a/exp(mu);
                       KK=b*(1-corr2)
                       KK1=(a+1)/(2+KK)
                       dd=  exp(log(b)+0.5*log(KK)+a*log(2+KK)-log(1+b)-(a+0.5)*log(4+KK))
                       aa=hypergeo::hypergeo((1 - a)/2, -a/2, 1, 4/(2+KK)^2)
                       bb=KK1*hypergeo::hypergeo((2-a)/2, (1-a)/2, 2, 4/(2+KK)^2)
                       cc=Re(corr2*(1-dd*(aa+bb)))
                       vv=exp(mu)*(1+1/b)
                       covariance=vv*cc
                       variogram=vv*(1-cc)
                           }
                   }
     if(poissonZIP) {
                    if(bivariate) {}
                    if(!bivariate) {   
                           p=pnorm(as.numeric(nuisance['pmu']));MM=exp(mu)
                           vv=(1-p)*MM*(1+p*MM) 
                           p1=1-2*p+pbivnorm::pbivnorm(as.numeric(nuisance['pmu']),as.numeric(nuisance['pmu']), rho =
                           (1-as.numeric(nuisance['nugget2']))*correlation, recycle = TRUE)
      
                           corr2=((1-as.numeric(nuisance['nugget1']))*correlation)^2    
                           z=2*vv/(1-corr2)
                           cc1=corr2*(1-(besselI(z,0,expon.scaled = TRUE)+besselI(z,1,expon.scaled = TRUE)))
                           cc=(p1*cc1*MM+MM^2*(p1-(1-p)^2))/vv
                           variogram=vv*(1-cc)
                      }  
            }
       if(poissongammaZIP) {
                    if(bivariate) {}
                    if(!bivariate) {   
                           p=pnorm(as.numeric(nuisance['pmu']));

                           ##....
                           cc=1
                           variogram=vv*(1-cc)
                      }  
            }




##########################################
############ starting graphics############
##########################################
      vario.main <- "Spatial semi-variogram"
      vario.ylab <- "Semi-Variogram"
        if(ispatim){
            dim(covariance) <- c(numlags_m, numlagt_m)
            dim(variogram) <- c(numlags_m, numlagt_m)
            vario.main <- "Space-time semi-variogram"
            vario.zlab <- "Semi-Variogram" }
        # compute the practical range:
        if(show.range || answer.range){
            Range <- uniroot(PracRangeNorm, c(lower, upper), corrmodel=corrmodel,
                             nuisance=nuisance, numlags=1, numlagt=1, lagt=lagt,
                             param=param, pract.range=pract.range)$root}



if(show.cov){

#####################################
#### bivariate case covariance ######
##################################### 
        if(bivariate&&!dyn){
            #par(mfrow=c(2,2))
       plot(lags_m, covariance11, type='l', ylim=c(min(covariance11),
                     max(covariance11)), main="First covariance",
                     xlab="Distance", ylab="Covariance",...)
       plot(lags_m, covariance12, type='l', ylim=c(min(covariance12),
                     max(covariance12)), main="Cross covariance",
                     xlab="Distance", ylab="Covariance",...)      
       plot(lags_m, covariance12, type='l', ylim=c(min(covariance12),
                     max(covariance12)), main="Cross covariance",
                     xlab="Distance", ylab="Covariance",...)   
       plot(lags_m, covariance22, type='l', ylim=c(min(covariance22),
                     max(covariance22)), main="Second covariance",
                     xlab="Distance", ylab="Covariance",...)          
         }
     if(bivariate&&dyn){
            #par(mfrow=c(2,2))
       plot(lags_m, covariance11, type='l', ylim=c(min(covariance11),
                     max(covariance11)), main="First covariance",
                     xlab="Distance", ylab="Covariance",...)
       plot(lags_m, covariance22, type='l', ylim=c(min(covariance22),
                     max(covariance22)), main="Second covariance",
                     xlab="Distance", ylab="Covariance",...)          
         }

        
#######################################
#### space time case  covariance ######
#######################################
if(ispatim){
            # build the covariance matrix:
            plagt <- !is.null(fix.lags)
            plags <- !is.null(fix.lagt)
            numplot <- 1+plags+plagt
          
            par(mai=c(.2,.2,.2,.2))
            persp(lags_m, lagt_m, covariance, xlab="Distance", ylab="Time",
                  zlab="Covariance", ltheta=90,
                  shade=0.75, ticktype="detailed", phi=30,
                  theta=30,main="Fitted space-time covariance",
                  , cex.axis=.8, cex.lab=.8,zlim=c(0,max(covariance))) #
            if(plagt){
                par(mai=c(.5,.5,.5,.5),mgp=c(1.6,.6,0))
                plot(lagt_m, covariance[fix.lags,], xlab="Time",
                     ylab="Covariance", type="l",cex.axis=.8,cex.lab=.8,
                     main="Space-time cov: temporal profile",...)}
                points()
            if(plags){
                par(mai=c(.5,.5,.5,.5),mgp=c(1.6,.6,0))
                plot(lags_m, covariance[,fix.lagt], xlab="Distance",
                     ylab="Covariance", type="l",cex.axis=.8,cex.lab=.8,
                     main="Space-time cov: spatial profile",...)
            }
        }
#######################################
#### spatial case  covariance #########
#######################################
       if(!ispatim && !bivariate){# spatial case:
            if(add.cov & dev.cur()!=1){
                lines(lags_m, covariance,...)
                if(show.range) abline(v=Range)}
            else{
                plot(lags_m, covariance, type='l', #ylim=c(0,max(covariance)), 
                     main="Spatial covariance",
                     xlab="Distance", ylab="Covariance",...)
                if(show.range) abline(v=Range)}}
         }
OLS=NULL
########################### 
########################### 
########################### 
    # display the semivariogram function
    if(show.vario){

#####################################
#### bivariate case semivariogram ###
##################################### 
 if(bivariate){
      if(bivariate&&!dyn){
       plot(vario$centers,vario$variograms[1,], main="First semi-variogram",ylim=c(0,max(vario$variograms[1,])),
           xlim=c(0,max(vario$centers)),
                     xlab="Distance", ylab="Semi-Variogram",...)
       lines(lags_m, variogram11, type='l',...)
       if(min(vario$variogramst)>0) {ll1=0;ll=max(vario$variogramst)}
       if(min(vario$variogramst)<0) {ll1=min(vario$variogramst);ll=-min(vario$variogramst)}
       plot(vario$centers,vario$variogramst, main="Cross semi-variogram",ylim=c(ll1,ll),
         xlim=c(0,max(vario$centers)),
                     xlab="Distance", ylab="Semi-Variogram",...)
       lines(lags_m, variogram12, type='l',...)
       plot(vario$centers,vario$variogramst, main="Cross semivariogram",ylim=c(ll1,ll),
         xlim=c(0,max(vario$centers)),
                     xlab="Distance", ylab="Semi-Variogram",...)
       lines(lags_m, variogram12, type='l',...)
       plot(vario$centers,vario$variograms[2,], main="Second semi-variogram",ylim=c(0,max(vario$variograms[2,])),
         xlim=c(0,max(vario$centers)),
                     xlab="Distance", ylab="Semi-Variogram",...)
       lines(lags_m, variogram22, type='l',...)  
  if(invisible)  OLS= sum( (vario$variograms[2,]-variogram22)^2) + sum((vario$variograms[1,]-variogram11))^2 + sum((vario$variogramst-variogram12))^2

   }
 
   if(bivariate&&dyn){
       plot(vario$centers,vario$variograms[1,], main="First semi-variogram",ylim=c(0,max(vario$variograms[1,])),
                     xlab="Distance", ylab="Semi-Variogram",...)
       lines(lags_m, variogram11, type='l',...)
       plot(vario$centers,vario$variograms[2,], main="Second semi-variogram",ylim=c(0,max(vario$variograms[2,])),
                     xlab="Distance", ylab="Semi-Variogram",...)
       lines(lags_m, variogram22, type='l',...)  

    if(invisible)  OLS= sum( (vario$variograms[2,]-variogram22)^2) + sum((vario$variograms[1,]-variogram11))^2

        }
}
#######################################
#### space time case  semivariogram ###
#######################################

 if(ispatim){

 # spatio-temporal case:
            plagt <- !is.null(fix.lags)
            plags <- !is.null(fix.lagt)
            nrowp <- 1
            ncolp <- 1
            if(isvario){
                ncolp <- ncolp+1
                if(plags || plagt) nrowp <- nrowp+1}
            else ncolp <- ncolp+plags+plagt
            sup <- 0
            tup <- 0
            if(isvario){
                nbins <- length(vario$centers)
                
            if(is.null(vario$centert)) nbint <- length(vario$bint)
            else  nbint <- length(vario$centert)
         # if(!dyn) {

            evario=matrix(vario$variogramst,nrow=length(vario$centers),ncol=nbint,byrow=TRUE)
             evario=rbind(c(0,vario$variogramt),cbind(vario$variograms,evario))

                  #evario <- matrix(vario$variogramst,nrow=nbins,ncol=nbint,byrow=TRUE)
                  #evario <- rbind(c(zero,vario$variogramt),cbind(vario$variograms,evario))
         
                  if(is.null(vario$centert)) 
                          evario.grid <- as.matrix(expand.grid(c(0,vario$centers),c(0,vario$bint)))
                  else  evario.grid <- as.matrix(expand.grid(c(0,vario$centers),c(0,vario$centert)))
            #      }
         ## else   {
           ##       evario <- matrix(vario$variogramst,nrow=nbins-1,ncol=nbint,byrow=TRUE)
             #     evario <- rbind(c(zero,vario$variogramt),cbind(vario$variograms,evario))
              #    evario.grid <- expand.grid(c(vario$centers),c(0,vario$bint))
               #   }
             
                scatterplot3d::scatterplot3d(evario.grid[,1],evario.grid[,2], c(evario),
                              type="h",highlight.3d=TRUE,cex.axis=.7,cex.lab=.7,
                              main=paste("Empirical",vario.main),xlab="Distance",
                              ylab="Time",zlab=vario.zlab,mar=c(2,2,2,2),mgp=c(0,0,0))
                if(plagt) tup <- max(evario[fix.lags,],na.rm=TRUE)
                if(plags) sup <- max(evario[,fix.lagt],na.rm=TRUE)
            } 

            par(mai=c(.2,.2,.2,.2))
            #print(variogram);print(lags_m);print(lagt_m)
            persp(lags_m, lagt_m, variogram, xlab="Distance",
                  ylab="Time", zlab=vario.zlab, ltheta=90,
                  shade=0.75, ticktype="detailed", phi=30,
                  theta=30,main=vario.main, cex.axis=.8,
                   cex.lab=.8)  #zlim=c(0,max(variogram))

            vvv=nuisance["sill"]
            ########
            if(gamma)         vvv=2*exp(mm["mean"])^2/as.numeric(nuisance["shape"])
            if(weibull)       vvv=exp(mm["mean"])^2*(gamma(1+2/as.numeric(nuisance["shape"]))/gamma(1+1/as.numeric(nuisance["shape"]))^2-1)
            if(loglogistic)   vvv=exp(mm["mean"])^2*
                               (2*as.numeric(nuisance['shape'])*sin(pi/nuisance['shape'])^2/(pi*sin(2*pi/nuisance['shape']))-1)
            if(loggauss)      vvv=(exp(nuisance["sill"])-1)#*(exp(mm['mean']))^2
            if(binomial)      vvv=fitted$n*pnorm(mm['mean'])*(1-pnorm(mm['mean']))
            if(poisson)       vvv=exp(mm['mean'])
            #if(poissongamma)  vvv=exp(mm['mean']*(1+1/nuisance["shape"]))
            if(poissongamma)  {MM=exp(mm['mean']); vvv=MM*(1+MM/as.numeric(nuisance["shape"]))
        }


            if(poissonZIP){        p=pnorm(nuisance['pmu']); MM=exp(mm['mean']); 
                              vvv=(1-p)*MM*(1+p*MM)}

            if(poissongammaZIP){     p=pnorm(nuisance['pmu']); MM=exp(mm['mean']);      ### to fix 
                              vvv=(1-p)*MM*(1+p*MM)}

            if(geom)          vvv= (1-pnorm(mm['mean']))/pnorm(mm['mean'])^2
            if(binomialneg)   vvv= fitted$n*(1-pnorm(mm['mean']))/pnorm(mm['mean'])^2
            if(binomialnegZINB) {
                                  MM=pnorm(mm['mean']);pg=pnorm(nuisance['pmu'])
                                  vvv=fitted$n*(1-MM)*(1-pg)*(1+fitted$n*pg*(1-MM)) /MM^2
                                }
            if(skewgausssian) vvv=(nuisance["sill"]+as.numeric(nuisance["skew"]))^2*(1-2/pi)
            if(studentT)      vvv=as.numeric(nuisance["df"])/(as.numeric(nuisance["df"])-2)
            if(tukeyh)        vvv=(1-2*as.numeric(nuisance["tail"]))^(-1.5)
            if(sas)           { d=as.numeric(nuisance['tail']); e=as.numeric(nuisance['skew'])
                                MM=sinh(e/d)*exp(0.25)*(besselK(.25,(d+1)/(2*d))+besselK(.25,(1-d)/(2*d)))/(sqrt(8*pi))
                                vvv=cosh(2*e/d)*exp(0.25)*(besselK(.25,(d+2)/(2*d))+besselK(0.25,(2-d)/(2*d)))/(sqrt(32*pi))-0.5-(MM)^2
                            }
            if(tukeyh2)  {        hr=as.numeric(nuisance["tail1"]);hl=as.numeric(nuisance["tail2"]);
                                  mm=(hr-hl)/(sqrt(2*pi)*(1-hl)*(1-hr))
                                  vvv=0.5*(1-2*hl)^(-3/2)+0.5*(1-2*hr)^(-3/2)-(mm)^2
                         }
            ########
            if(plagt){
                par(mai=c(.5,.5,.3,.3),mgp=c(1.6,.6,0))
                plot(lagt_m, variogram[fix.lags,], xlab="Time",cex.axis=.8,cex.lab=.8,
                     ylab=vario.ylab, type="l", ylim=c(0,max(vvv,tup)), main=paste(vario.ylab,": temporal profile",
                     sep=""),...)
                if(isvario) points(lagt[-1], evario[fix.lags,][-1],...)
               }

            if(plags){
                par(mai=c(.5,.5,.3,.3),mgp=c(1.6,.6,0))
                plot(lags_m, variogram[,fix.lagt], xlab="Distance",cex.axis=.8,cex.lab=.8,
                     ylab=vario.ylab, type="l", ylim=c(0,max(vvv,sup)), main=paste(vario.ylab,": spatial profile",
                     sep=""),...)
                if(isvario) points(lags[-1], evario[,fix.lagt][-1],...)
                }
  #################################

                      if(invisible) OLS=sum(evario-variogram)^2
 } ####### end space time
          
################################
#### spaatial semivariogram ###
#################################
    if(!ispatim && !bivariate){

   #################################
        if(invisible) OLS=sum(vario$variograms-variogram)^2
  #################################
            if(add.vario & dev.cur()!=1){
                points(vario$centers, vario$variograms,...)
                lines(lags_m, variogram,...)
                if(show.range) abline(v=Range)}
            else{
                bnds <- range(variogram)
                bnds[1] <- min(bnds[1], min(vario$variograms))
                bnds[2] <- max(bnds[2], max(vario$variograms))
        
                plot(lags_m, variogram, type='l',  
                     main=vario.main,xlab="Distance",
                     ylab=vario.ylab,...)
                points(vario$centers, vario$variograms,...)
                if(show.range) abline(v=Range)}
        }
 #################################       
}
    if(ispatim) par(mai=c(1.02 ,0.85 ,0.85 ,0.45),mgp=c(3,1,0))
    # return the estimated covariance function
    if(answer.cov) {result <- list(lags=lags_m,lagt=lagt_m, covariance=covariance)}
    # return the estimated variogram function
    if(answer.vario) {
        if(!is.list(result)) {if(!bivariate) if(gaussian) result <- list(lags=lags_m,lagt=lagt_m, variogram=variogram)
                              if(bivariate)  if(gaussian) result <- list(lags=lags_m,lagt=lagt_m, variogram11=variogram11,
                                                                        variogram12=variogram12,variogram22=variogram22)          
                              }
        else {
            if(!bivariate){if(gaussian) result$variogram <- variogram}
            if(bivariate){
                if(gaussian) {result$variogram11 <- variogram11;result$variogram12 <- variogram12;result$variogram22 <- variogram22}
                }}}
   ########################################################## 
if(!is.null(result)) {
        invisible(result)
    } else if(!is.null(OLS)) {
        invisible(OLS)
    } else {
        invisible()
    }
  }
