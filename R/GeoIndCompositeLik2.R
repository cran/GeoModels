####################################################
### File name: CompIndLik2.r
####################################################

### Optim call for Indipendence Composite log-likelihood maximization

CompIndLik2 <- function(bivariate, coordx, coordy ,coordz,coordt,coordx_dyn, data, flagcorr, flagnuis, fixed,grid,
                           lower, model, n, namescorr, namesnuis, namesparam,
                           numparam,  optimizer, onlyvar, param, spacetime, type,
                           upper,namesupper, varest, ns, X,sensitivity,copula,MM)
  {

############# utility functions ############
lambertW0_custom <- function(x, tol = 1e-12, maxiter = 100) {
  if (length(x) > 1) {
    return(sapply(x, function(xi) lambertW0_custom(xi, tol, maxiter)))
  }
  if (is.na(x) || is.infinite(x)) return(x)
  if (x < -1/exp(1)) return(NaN)
  if (x == 0) return(0)
  if (x == -1/exp(1)) return(-1)
  if (x > 3) {
    # Per valori grandi: W(x) ≈ log(x) - log(log(x))
    w <- log(x) - log(log(x))
  } else if (x > -0.32358) {
    if (x < 1.5) {
      # Serie di potenze modificata
      w <- x * (1 - x + (3/2) * x^2 - (8/3) * x^3 + (125/24) * x^4)
    } else {
      # Approssimazione logaritmica migliorata
      lx <- log(x)
      w <- lx - log(lx) + lx/(1 + lx)
    }
  } else {
    sigma <- -1 - exp(1) * x
    if (sigma > 0) {
      p <- sqrt(2 * sigma)
      w <- -1 + p * (1 - p/3 + 11*p^2/72 - 43*p^3/540 + 769*p^4/17280)
    } else {
      w <- log(-x) - log(-log(-x))
    }
  }
  for (iter in 1:maxiter) {
    ew <- exp(w)
    wew <- w * ew
    f <- wew - x
    if (abs(f) <= tol * (1 + abs(x))) break
    fprime <- ew * (1 + w)
    if (abs(fprime) < 1e-20) break  # Evita divisione per zero
    delta <- f / fprime
    w_new <- w - delta
    if (abs(delta) <= tol * (1 + abs(w))) break
    w <- w_new
  }
  return(w)
}
one_log_tukeyh <- function(data, mm, sill, tail) {
  q <- (data - mm) / sqrt(sill)
  aa <- tail * q * q
  lambert_vals <- lambertW0_custom(aa)
  x <- sign(q) * sqrt(lambert_vals / tail)
  extra <- 1 / (1 + lambert_vals)
  return(log(dnorm(x) * x * extra / (q * sqrt(sill))))
}
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

indloglik<- function(fan,data,mm,nuis){

## gaussian and misspecified gaussian
    if(fan=="Ind_Pair_Gauss")                  { sill=nuis[1];
                                                 res=sum(dnorm(data,mean =mm,sd=sqrt(sill), log = TRUE))
                                               }

    if(fan=="Ind_Pair_Gauss_misp_T")           {    df=1/nuis[1];sill=nuis[2]
                                                    vv=sill*df/(df-2)
                                                  res=sum(dnorm(data,mean =mm,sd=sqrt(vv), log = TRUE))
                                               }

    if(fan=="Ind_Pair_Gauss_misp_SkewT")        {    df=1/nuis[1];sill=nuis[2];skew=nuis[3]; 
                                                     kk=gamma(0.5*(df-1))/gamma(0.5*(df)) 
                                                     me=  sqrt(sill)*sqrt(df/pi)*kk*skew  ;    
                                                     vv=  sill*(df/(df-2) - me^2)
                                                    res=sum(dnorm(data,mean=(mm+me),sd=sqrt(vv), log = TRUE)) 
                                                }

    if(fan=="Ind_Pair_Gauss_misp_Tukeygh")      {   sill =nuis[1];eta  = nuis[2];tail = nuis[3]; 
                                                    eta2=eta*eta; u=1-tail;
                                                    me=sqrt(sill)*(exp(eta2/(2*u))-1)/(eta*sqrt(u));
                                                    vv=sill*((exp(2*eta2/(1-2*tail))-2*exp(eta2/(2*(1-2*tail)))+1)/(eta2*sqrt(1-2*tail)))-me*me;
                                                  res=sum(dnorm(data, mean=mm+me,sd=sqrt(vv),  log = TRUE))
                                                }

## non gaussian over R
    if(fan=="Ind_Pair_T")                     {  sill=nuis[2]; df=1/nuis[1]
                                                res=sum(dt((data-mm)/sqrt(sill), df, log = TRUE)-0.5*log(sill))
                                              }  
    if(fan=="Ind_Pair_Logistic")              {  sill=nuis[1]; 
                                                res=sum(dlogis(data,mm, sqrt(sill) ,log = TRUE))
                                              }  
    if(fan=="Ind_Pair_SkewGauss")   {
                                                sk=nuis[2];sill=nuis[1]
                                                omega=sk*sk + sill
                                                #alpha=sk/(sill^0.5)
                                               #res=sum( sn::dsn((data-mm)/sqrt(omega),xi=0,omega= 1,alpha= alpha,log=TRUE)-0.5*log(omega)) 
                                                q=data-mm
                                               res=sum(log(2)-0.5*log(omega)+dnorm(q/(sqrt(omega)),log=TRUE)+pnorm( (sk*q)/(sqrt(sill)*sqrt(omega)),log.p=TRUE))
                                   }


   if(fan=="Ind_Pair_SinhGauss")   {           
                                                skew=nuis[2];sill=nuis[1];tail=nuis[3]
                                                q=(data-mm)/(sqrt(sill));
                                                b1=tail*asinh(q)-skew;Z1=sinh(b1);
                                                res=sum(-0.5*log(q^2+1)-0.5*log(2*pi*sill)+log(cosh(b1))+log(tail)-Z1*Z1/2);
                                       }

  if(fan=="Ind_Pair_Tukeyh")        { 
                                           sill=nuis[1];tail=nuis[2]
                                           res=sum(one_log_tukeyh(data,mm,sill,tail)) 
                                    }
   if(fan=="Ind_Pair_Tukeyhh")        { sill=nuis[1];tail1=nuis[3];tail2=nuis[2];
    res=sum(  log(  exp(one_log_tukeyh(data,mm,sill,tail2))*I(data>=mm) +
                    exp(one_log_tukeyh(data,mm,sill,tail1))*I(data<mm) )) 
                                      }
   if(fan=="Ind_Pair_TWOPIECEGauss") { sill=nuis[1];eta=nuis[2]
                                        y=(data-mm)/sqrt(sill)
                    res=sum( log( dnorm(y/(1-eta))*I(y>=0) + dnorm(y/(1+eta))*I(y<0)) -0.5*log(sill))
                                      }  
   if(fan=="Ind_Pair_TWOPIECETukeyh") { sill=nuis[1];eta=nuis[2];tail=nuis[3]
                                        y=(data-mm)/sqrt(sill);
                    res= sum(  log( exp(one_log_tukeyh(y/(1-eta),0,1,tail))*I(y>=0) + 
                                    exp(one_log_tukeyh(y/(1+eta),0,1,tail))*I(y<0) )     -0.5*log(sill) )
                                       }
   if(fan=="Ind_Pair_TWOPIECET") { sill=nuis[2];eta=nuis[3];tail=1/nuis[1]
                                  y=(data-mm)/sqrt(sill);
                    res=sum( log(  dt(y/(1-eta),df=tail)*I(y>=0) + dt(y/(1+eta),df=tail)*I(y<0)) -0.5*log(sill))
                                       }
## non gaussian over R^+

if(fan== "Ind_Pair_Gamma")                 {
                                              shape=nuis[2]
                                         
                                              res=sum(dgamma(data, shape=shape/2, scale = 1/(shape/(2*exp(mm))), log = TRUE))
                                            }
if(fan== "Ind_Pair_Weibull")              {

                                              shape=nuis[2] 
                                              res=sum(dweibull(data, shape=shape,scale=exp(mm)/(gamma(1+1/shape)) , log = TRUE))
                                            }
 if(fan== "Ind_Pair_LogGauss")              {
                                              sill=nuis[1]
                                              c1=mm-sill/2
                                              res=sum(dlnorm(data, meanlog =c1, sdlog = sqrt(sill), log = TRUE))
                                            }
if(fan== "Ind_Pair_LogLogistic")            {    
                                             dllogis1 <- function(x, shape, rate = 1, scale = 1/rate, log = FALSE) {
                                               z <- x / scale
                                               num <- (shape / scale) * z^(shape - 1)
                                               denom <- (1 + z^shape)^2
                                               dens <- num / denom
                                               if (log) dens <- log(dens)
                                               return(dens)
                                             }
                                              shape=nuis[2]
                                              ci=gamma(1+1/shape)*gamma(1-1/shape)
                                              res=sum(dllogis1(data, shape, scale = exp(mm)/ci, log = TRUE))
                                            }
## non Gassian bounded support
  if(fan== "Ind_Pair_Beta2")                {    mmax=nuis[4];mmin=nuis[3]
                                                 shape=nuis[2]
                                    
                                                 me=1/(1+exp(-mm))
                                                 res=sum(dbeta((data-mmin)/(mmax-mmin), me*shape, (1-me)*shape,log=TRUE)-log(mmax-mmin))
                                            }
 if(fan== "Ind_Pair_Kumaraswamy2")                {  
                                             mmax=nuis[4];mmin=nuis[3]
                                             shape=nuis[2]
                                             q=(data-mmin)/(mmax-mmin);k=1-q^shape
                                             m1=1/(1+exp(-mm));
                                             shapei=log(0.5)/log1p(-m1^shape);
                                             res=sum(log(shapei)+log(shape)+(shape-1)*log(q)+(shapei-1)*log(k)-log(mmax-mmin))            
                                             }

#### discrete
    if(fan=="Ind_Pair_Pois")              {mm=exp(mm);res=sum(dpois(data, mm, log = TRUE)) }
    if(fan=="Ind_Pair_Gauss_misp_Pois")   {mm=exp(mm);res=sum(dnorm(data, mean = mm, sd =sqrt(mm), log = TRUE))}
    if(fan=="Ind_Pair_BinomGauss")        res=sum(dbinom(data, n, pnorm(mm), log = TRUE))
    if(fan=="Ind_Pair_BinomGauss_misp")   { pp=pnorm(mm); mm=n*pp;vv=mm*(1-pp)
                                           res=sum(dnorm(data, mean =mm , sd =sqrt(vv), log = TRUE))}
    if(fan=="Ind_Pair_BinomnegGauss")     res=sum(dnbinom(data, n, pnorm(mm), log = TRUE))
    if(fan=="Ind_Pair_PoisGamma")         {mm=exp(mm);res=sum(dnbinom(data, nuis[2], mu=mm, log = TRUE))}

    if(fan=="Ind_Pair_PoisGammaZIP")      {mm=exp(mm);pp=pnorm(nuis[3])
                                           res1=sum(log(pp+(1-pp)*dnbinom(data[data==0], nuis[2],mu=mm[data==0],0)));
                                           res2=sum(log(1-pp)+dnbinom(data[data!=0],nuis[2],mu=mm[data!=0],1)); res=res1+res2}


    if(fan=="Ind_Pair_Gauss_misp_PoisGamma") {mm=exp(mm);res=sum(dnorm(data, mean = mm, sd =sqrt(mm*(1+mm/nuis[2])), log = TRUE))}
    if(fan=="Ind_Pair_PoisZIP")           {mm=exp(mm);pp=pnorm(nuis[3])
                                            res1=sum(log(pp+(1-pp)*dpois(data[data==0],mm[data==0],0)));
                                            res2=sum(log(1-pp)+dpois(data[data!=0],mm[data!=0],1)); res=res1+res2}
    if(fan=="Ind_Pair_Gauss_misp_PoisZIP"){mm=exp(mm);pp=pnorm(nuis[3])
                                            res1=sum(log(pp+(1-pp)*dnorm(data[data==0],mean=mm[data==0],sd =sqrt(mm[data==0]), log = FALSE)));
                                            res2=sum(log(1-pp)+dnorm(data[data!=0],mean=mm[data!=0], sd =sqrt(mm[data!=0]), log = TRUE)); res=res1+res2}
    if(fan=="Ind_Pair_BinomnegGaussZINB") {pm=pnorm(mm);pp=pnorm(nuis[3])
                                            res1=sum(log(pp+(1-pp)*dnbinom(data[data==0],n,pm[data==0],log = FALSE)));
                                            res2=sum(log(1-pp)+dnbinom(data[data!=0],n, pm[data!=0], log = TRUE)); res=res1+res2}
### .......
return(-res)
}

 compindloglik2 <- function(param,  data,fixed, fan, n, 
                              namesnuis,namesparam,X,MM)
      {
        names(param) <- namesparam
        param <- c(param, fixed)
        nuisance <- param[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
        mm=as.numeric(nuisance[sel])   ## mean paramteres
        other_nuis=as.numeric(nuisance[!sel])   ## or nuis parameters (nugget sill skew df)
        if((is.null(MM))) Mean=c(X%*%mm)
        else Mean=c(MM)
        result=indloglik(fan,data,Mean,other_nuis)
        return(result)
      }

 compindloglik_biv2 <- function(param, data1,data2,fixed, fan, n,   ## to do...
                           namesnuis,namesparam,X,MM)
      {

        names(param) <- namesparam
        param <- c(param, fixed)
        nuisance <- param[namesnuis]
        sel1=substr(names(nuisance),1,6)=="mean_1"
        mm1=as.numeric(nuisance[sel1])
        sel2=substr(names(nuisance),1,6)=="mean_2"
        mm2=as.numeric(nuisance[sel2])
        sel=substr(names(nuisance),1,4)=="mean"
        X1=as.matrix(X[1:ns[1],]);X2=as.matrix(X[(ns[1]+1):(ns[2]+ns[1]),]); 
        other_nuis=as.numeric(nuisance[!sel]) 
        Mean=c(X1%*%mm1,X2%*%mm2)
        res=double(1)
     
        result=4
        return(-result)
      }

   ##################################################################################################
   ############### starting function ###############################################################
   ##################################################################################################

    numcoord=length(coordx);numtime=1;spacetime_dyn=FALSE
    if(spacetime) numtime=length(coordt)
    if(bivariate) numtime=2
    if(!is.null(coordx_dyn)) spacetime_dyn=TRUE

    dimat=numcoord*numtime#
    if((spacetime||bivariate)) dimat <- sum(ns)
    NS=cumsum(ns)
    if(is.null(dim(X))){X=as.matrix(rep(1,dimat))}
    #else(if(bivariate) X=rbind(X,X))


    fname <- NULL; hessian <- FALSE
     
###################### pairwise ###############################################

    # --- lookup vector ---------------------------------------------------------
#  (keep only the models that are actually used)
lookup <- c(
  Gauss                    =  1,
  BinomGauss               = 11,          # duplicate key: last one wins
  BinomnegGauss            = c(14, 16),   # both map to the same name
  WrapGauss                = 13,
  SkewGauss                = 10,
  Gamma                    = 21,
  Kumaraswamy2             = 42,
  Beta2                    = 50,
  Weibull                  = 26,
  LogLogistic              = 24,
  Logistic                 = 25,
  LogGauss                 = 22,
  TWOPIECET                = 27,
  TWOPIECEGauss            = 29,
  T                        = 12,
  Tukeyh                   = 34,
  Tukeyhh                  = 40,
  Gauss_misp_Tukeygh       = 41,
  Gauss_misp_Pois          = 36,
  Gauss_misp_T             = 35,
  Gauss_misp_SkewT         = 37,
  SinhGauss                = 20,
  TWOPIECETukeyh           = 38,
  Pois                     = 30,
  PoisGamma                = 46,
  PoisGammaZIP             = 57,
  PoisZIP                  = 43,
  Gauss_misp_PoisZIP       = 44,
  BinomnegGaussZINB        = 45,
  Gauss_misp_PoisGamma     = 47,
  BinomGauss_misp          = 51
)

# reverse the map: model numbers → names
model2name <- setNames(rep(names(lookup), lengths(lookup)),
                       unlist(lookup, use.names = FALSE))
fname <- paste0("Ind_Pair_", model2name[ as.character(model) ])
if (is.na(fname)) stop("Unknown model: ", model)
   
########################################################################                                            
    if(sensitivity) hessian=TRUE
    
     if((spacetime||bivariate)&&(!spacetime_dyn))    data=c(t(data))
     if((spacetime||bivariate)&&(spacetime_dyn))     data=unlist(data)          
     if(spacetime||bivariate)   NS=c(0,NS)[-(length(ns)+1)]

####     
tot=c(param,fixed) ## all the parameters
## deleting corr para
tot=tot[is.na(pmatch(names(tot),namescorr))]
## deleting nugget
tot=tot[names(tot)!='nugget'];namesnuis=namesnuis[namesnuis!="nugget"]
param=tot[pmatch(namesparam,names(tot))];param=param[!is.na(names(param))];param=param[!is.na(param)];

if(model  %in%  c(2,14,16,21,42,50,26,24,25,30,46,43,11))  ## model where sill must be fixed = to 1
 {param=param[names(param)!='sill'];a=1; names(a)="sill";
  if(is.null(unlist(fixed['sill']))) fixed=c(fixed,a)}

if(!is.null(copula))
    {if(copula=="Clayton")
           {param=param[names(param)!='nu'];a=2; names(a)="nu";
           if(is.null(unlist(fixed['nu']))) fixed=c(fixed,a)}}

namesparam=names(param)
###updating upper and lower bound if necessary
sel=pmatch(namesparam,namesupper)
lower=lower[sel]
upper=upper[sel]
#### not exactly zero for the mean parameters starting values
sel=substr(names(param),1,4)=="mean"&param==0
param[sel]=0.1
param=as.numeric(param)

   if(!onlyvar){
  ##############################.  spatial or space time ############################################
   if(!bivariate)           {
    if(length(param)==1) {
      
         optimizer="optimize"  
         if(is.na(lower)||is.na(upper))  {
            if(model %in% c(2,14,16,45,11,30,36)) {lower=-5;upper=5}
            else                            {lower=-1e+10;upper=1e+10}
           }
         else{
             if(model %in% c(2,14,16,45,11,30,36)) {lower=-5;upper=5}
         }  
     CompLikelihood <- optimize(f= compindloglik2,    
                              data=data, fixed=fixed, fan=fname,  lower=lower, n=n,
                               namesnuis=namesnuis,namesparam=namesparam, maximum = FALSE,
                              upper= upper,  X=X,MM=MM)}
   if(length(param)>1) {
    
    if(optimizer=='L-BFGS-B'){
      CompLikelihood <- optim(par=param,fn= compindloglik2, 
                              control=list(factr=1e-10,pgtol=1e-14, maxit=100000), 
                                data=data, fixed=fixed,
                              fan=fname, lower=lower, method='L-BFGS-B',n=n,
                               namesnuis=namesnuis,namesparam=namesparam, 
                              upper=upper,  X=X,MM=MM,   hessian=FALSE)
  }

    if(optimizer=='BFGS') 
        CompLikelihood <- optim(par=param, fn= compindloglik2,     
                           control=list(factr=1e-10,
                             reltol=1e-14, maxit=100000),data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='BFGS',n=n,
                              namesnuis=namesnuis,namesparam=namesparam,  X=X,MM=MM)
   if(optimizer=='SANN'){ 
      CompLikelihood <- optim(par=param, fn= compindloglik2,     
                           control=list(factr=1e-10,
                             reltol=1e-14, maxit=100000),data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='SANN',n=n,
                              namesnuis=namesnuis,namesparam=namesparam,  X=X,MM=MM)
  }
      if(optimizer=='Nelder-Mead')
        CompLikelihood <- optim(par=param, fn= compindloglik2,     
          control=list( reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='Nelder-Mead',n=n,
                                  namesnuis=namesnuis,namesparam=namesparam,  X=X,MM=MM)
    if(optimizer=='nlm')
    CompLikelihood <- nlm(f= compindloglik2,p=param,steptol = 1e-4,    data=data, fixed=fixed,
                               fan=fname,hessian=FALSE,n=n, namesnuis=namesnuis,namesparam=namesparam, 
                               iterlim=100000,   X=X,MM=MM)
  
    if(optimizer=='nlminb')
     CompLikelihood <-nlminb(objective= compindloglik2,start=param,   data=data, fixed=fixed,
                                control = list( iter.max=100000),
                              lower=lower,upper=upper,
                               fan=fname,n=n, namesnuis=namesnuis,namesparam=namesparam, 
                                 X=X,MM=MM)                   
    }}
######################################################################################
############################## bivariate  ############################################ 
######################################################################################                          
    if(bivariate)           { 
     if(length(param)==1)
        {
         optimizer="optimize" 

       CompLikelihood <- optimize(f= compindloglik_biv2,     
                              data=data, fixed=fixed,fan=fname,  lower=lower,n=n,
                               namesnuis=namesnuis, namesparam=namesparam,maximum = FALSE,
                              upper=upper,  X=X,MM=MM )}
      if(length(param)>1) {   
    if(optimizer=='L-BFGS-B'){
      CompLikelihood <- optim(param, compindloglik_biv2, control=list(pgtol=1e-14, maxit=100000),
                              method='L-BFGS-B',hessian=FALSE,lower=lower, upper=upper,
                                  
                              data=data, fixed=fixed,fan=fname,n=n,
                               namesnuis=namesnuis, namesparam=namesparam,
                               X=X ,MM=MM )}

    
      if(optimizer=='BFGS')
      CompLikelihood <- optim(param, compindloglik_biv2,     control=list(
                              reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='BFGS',n=n,
                               namesnuis=namesnuis,namesparam=namesparam ,  X=X,MM=MM )
   if(optimizer=='Nelder-Mead')
      CompLikelihood <- optim(param, compindloglik_biv2,     control=list(
                              reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='Nelder-Mead',n=n,
                               namesnuis=namesnuis,namesparam=namesparam ,  X=X,MM=MM )
     if(optimizer=='nlm') 
        CompLikelihood <- nlm( f= compindloglik_biv2,p=param,     data=data, fixed=fixed,
                               fan=fname,hessian=FALSE,n=n, namesnuis=namesnuis,namesparam=namesparam, 
                                 X=X,MM=MM )
    if(optimizer=='nlminb') 
        CompLikelihood <- nlminb( objective= compindloglik_biv2,start=param, 
                                     control = list( iter.max=100000),
                              lower=lower,upper=upper,
                                   data=data, fixed=fixed,
                               fan=fname,n=n, namesnuis=namesnuis,namesparam=namesparam, 
                                 X=X ,MM=MM)
   }
 }  

                   
      ########################################################################################   
      ########################################################################################
    # check the optimisation outcome
      if(optimizer=='Nelder-Mead'||optimizer=='SANN'){
        CompLikelihood$value = -CompLikelihood$value
        names(CompLikelihood$par)<- namesparam
        if(CompLikelihood$convergence == 0)
        CompLikelihood$convergence <- 'Successful'
        else
        if(CompLikelihood$convergence == 1)
        CompLikelihood$convergence <- 'Iteration limit reached'
        else
        CompLikelihood$convergence <- "Optimization may have failed"
        if(CompLikelihood$value==-1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
    }
        if(optimizer=='nmk'||optimizer=='nmkb'){
        CompLikelihood$value = -CompLikelihood$value
        names(CompLikelihood$par)<- namesparam
        if(CompLikelihood$convergence == 0)
        CompLikelihood$convergence <- 'Successful'
        else CompLikelihood$convergence <- "Optimization may have failed"
        if(CompLikelihood$value==-1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
    }
      
    if(optimizer=='L-BFGS-B'||optimizer=='BFGS'){
        CompLikelihood$value = -CompLikelihood$value
        names(CompLikelihood$par)<- namesparam
        if(CompLikelihood$convergence == 0)
        CompLikelihood$convergence <- 'Successful'
        else
        if(CompLikelihood$convergence == 1)
        CompLikelihood$convergence <- 'Iteration limit reached'
        else
        CompLikelihood$convergence <- "Optimization may have failed"
        if(CompLikelihood$value==-1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
    }

     if(optimizer=='nlm'){
        CompLikelihood$par <- CompLikelihood$estimate
        names(CompLikelihood$par)<- namesparam
        CompLikelihood$value <- -CompLikelihood$minimum
        if(CompLikelihood$code == 1|| CompLikelihood$code == 2)
        CompLikelihood$convergence <- 'Successful'
        else
        if(CompLikelihood$code == 4)
        CompLikelihood$convergence <- 'Iteration limit reached'
        else
        CompLikelihood$convergence <- "Optimization may have failed"
        if(CompLikelihood$value==-1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
    }

    if(optimizer=='nlminb'||optimizer=='multinlminb'){
        CompLikelihood$par <- CompLikelihood$par
        names(CompLikelihood$par)<- namesparam
        CompLikelihood$value <- -CompLikelihood$objective
        if(CompLikelihood$convergence == 0) { CompLikelihood$convergence <- 'Successful' }
        else {CompLikelihood$convergence <- "Optimization may have failed" }
        if(CompLikelihood$objective==-1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
    }
    if(optimizer=='optimize'){
    param<-CompLikelihood$minimum
    CompLikelihood$par<-param  
    names(CompLikelihood$par)<- namesparam
    maxfun <- -CompLikelihood$objective
    CompLikelihood$value <- maxfun
    CompLikelihood$convergence <- 'Successful'
    }
  } ##### end if !onlyvar
    else {
          CompLikelihood=as.list(0)
          names(CompLikelihood)="value"
          CompLikelihood$par <- param
          CompLikelihood$claic <- NULL;CompLikelihood$clbic <- NULL;
          CompLikelihood$convergence <- 'Successful'
          if(!bivariate) CompLikelihood$value = -  compindloglik2(param=CompLikelihood$par ,     
                              data=data, fixed=fixed, fan=fname,
                             n=n,namesnuis=namesnuis,namesparam=namesparam,  X=X,MM=MM)
          else CompLikelihood$value = - compindloglik_biv2(param=CompLikelihood$par ,     
                data=data, fixed=fixed, fan=fname,
                             n=n,namesnuis=namesnuis,namesparam=namesparam,  X=X,MM=MM)

          if(hessian) 
          {
               if(!bivariate)  
                CompLikelihood$hessian=numDeriv::hessian(func= compindloglik2,x=param,method="Richardson",     
                              data=data,fixed=fixed,fan=fname,n=n,
                               namesnuis=namesnuis, namesparam=namesparam,
                                X=X,MM=MM )
               if(bivariate)  
               CompLikelihood$hessian=numDeriv::hessian(func= compindloglik_biv2,x=param,method="Richardson",   
                             data=data, fixed=fixed,fan=fname,n=n, namesnuis=namesnuis,namesparam=namesparam, 
                               X=X ,MM=MM)
               rownames(CompLikelihood$hessian)=namesparam
               colnames(CompLikelihood$hessian)=namesparam
          }
  }

#####################################
if((sensitivity||varest))
  {
if(!bivariate)  
CompLikelihood$hessian=numDeriv::hessian(func= compindloglik2,x=CompLikelihood$par,method="Richardson",      
                              data=data, fixed=fixed,fan=fname,n=n,
                               namesnuis=namesnuis, namesparam=namesparam,
                                X=X ,MM=MM)
if(bivariate)  
CompLikelihood$hessian=numDeriv::hessian(func= compindloglik_biv2,x=CompLikelihood$par,method="Richardson",    
                          data=data, fixed=fixed,fan=fname,n=n, namesnuis=namesnuis,namesparam=namesparam, 
                                 X=X ,MM=MM)
rownames(CompLikelihood$hessian)=namesparam
colnames(CompLikelihood$hessian)=namesparam
  }

if(hessian) CompLikelihood$sensmat=CompLikelihood$hessian

    return(CompLikelihood)
  }

