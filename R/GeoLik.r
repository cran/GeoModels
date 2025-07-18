####################################################
### File name: GeoLik.r
####################################################


### Optim call for log-likelihood maximization 
Lik <- function(copula,bivariate,coordx,coordy,coordz,coordt,coordx_dyn,corrmodel,data,fixed,flagcor,flagnuis,grid,lower,
                       mdecomp,model,namescorr,namesnuis,namesparam,numcoord,numpairs,numparamcor,numtime,
                       optimizer,onlyvar,param,radius,setup,spacetime,sparse,varest,taper,type,upper,ns,X,neighb,MM,aniso)
{
 ######### computing upper trinagular of covariance matrix   
    matr<- function(corrmat,corr,coordx,coordy,coordz,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius)
    {


  #      cc <- .C(as.character(corrmat),
  #  cr=double(length(corr)),as.double(coordx), as.double(coordy),as.double(coordz),as.double(coordt), as.integer(corrmodel),
  #              as.double(nuisance),as.double(paramcorr),as.double(radius), as.integer(ns),as.integer(NS),
  #          PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)$cr

           cc=dotCall64::.C64(as.character(corrmat),
        SIGNATURE =c(rep("double",5),"integer","double","double","double","integer","integer"),
                          cr=dotCall64::vector_dc("double",length(corr)), coordx, coordy, coordz,coordt, corrmodel, nuisance,paramcorr,radius, ns,NS,
                INTENT = c("w",rep("r",10)),
              PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$cr
        return(cc)
    }
   ######### computing upper trinagular of covariance matrix   it is for poisson
    matr2 <- function(corrmat,corr,coordx,coordy,coordz,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius,model,mu)
    {

  #     cc <- .C(as.character(corrmat),
  #     cr=double(length(corr),as.double(coordx), as.double(coordy), as.double(coordz),as.double(coordt), as.integer(corrmodel), as.double(c(mu)),
  #                        as.integer(ns), as.double(nuisance),as.double(paramcorr),as.double(radius), as.integer(ns),as.integer(NS),as.integer(model),
  #      PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)$cr
        
  cc=dotCall64::.C64(as.character(corrmat),
         SIGNATURE = c(rep("double",5),"integer","double", "integer","double","double","double","integer","integer","integer"),  
                          cr=dotCall64::vector_dc("double",length(corr)), coordx, coordy, coordz,coordt, corrmodel, c(mu),
                          ns, nuisance,
                          paramcorr,
                          radius, ns,NS,model,
         INTENT =    c("rw",rep("r",13)),
             PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$cr
        return(cc)
    }
    ### START Defining the objective functions

######### 
build_correlation_matrix <- function(corr_vec, ident) {
    corrmat_full <- ident  
    # Riempimento efficiente
    lt <- lower.tri(ident)
    corrmat_full[lt] <- corr_vec
    corrmat_full <- t(corrmat_full)  # Trasponi
    corrmat_full[lt] <- corr_vec  # Ri-riempi il triangolo inferiore
    
    return(corrmat_full)
}
######### Restricted log-likelihood for multivariate normal density:
    LogNormDenRestr <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        llik <- 1.0e8
        # Computes the covariance matrix:
        varcov <- (nuisance['sill'])*ident
        varcov[lower.tri(varcov)] <- cova
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- cova
        #  decomposition of the covariance matrix:
        decompvarcov <- MatDecomp(varcov,mdecomp)
        if(is.logical(decompvarcov)) return(llik)       
        logdetvarcov <- MatLogDet(decompvarcov,mdecomp) 
        ivarcov <- MatInv(varcov)
        sumvarcov <- sum(ivarcov)
        p <- ivarcov-array(rowSums(ivarcov),c(dimat,1))%*%colSums(ivarcov)/sumvarcov
        llik <- 0.5*(const+logdetvarcov+log(sumvarcov)+crossprod(t(crossprod(stdata,p)),stdata))
        return(llik)
    }
######### Restricted log-likelihood for bivariate multivariate normal density:
         LogNormDenRestr_biv <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        llik <- 1.0e8
        # Computes the covariance matrix:
        ident[lower.tri(ident,diag=T)] <- cova
        ident <- t(ident)
        ident[lower.tri(ident,diag=T)] <- cova
        #  decomposition of the covariance matrix:
        decompvarcov <- MatDecomp(ident,mdecomp)
        if(is.logical(decompvarcov)) return(llik)       
        logdetvarcov <- MatLogDet(decompvarcov,mdecomp) 
        ivarcov <- MatInv(ident)
        sumvarcov <- sum(ivarcov)
        p <- ivarcov-array(rowSums(ivarcov),c(dimat,1))%*%colSums(ivarcov)/sumvarcov
        llik <- 0.5*(const+ logdetvarcov +log(sumvarcov)+crossprod(t(crossprod(stdata,p)),stdata))
        return(llik)
    }
    
    
######### Tapering 2 log-likelihood for multivariate normal density:
    LogNormDenTap <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        lliktap <- 1.0e8
        # Computes the vector of the correlations:
        varcovtap <- new("spam",entries=cova*setup$taps,colindices=setup$ja,
        rowpointers=setup$ia,dimension=as.integer(rep(dimat,2)))
        cholvctap <- spam::update.spam.chol.NgPeyton(setup$struct, varcovtap)
        #if(class(cholvctap)=="try-error") return(lliktap)
        logdet <- c(spam::determinant(cholvctap)$modulus)
        inv <- spam::solve.spam(cholvctap)
        slot(varcovtap,"entries") <- inv[setup$idx]*setup$taps
        lliktap= 0.5*(const+2*logdet+drop(t(stdata)%*%varcovtap%*%stdata))
        return(lliktap)
    }

######### Tapering 1 sas log-likelihood 
LogshDenTap1<- function(const,cova,ident,dimat,mdecomp,nuisance,setup,sill,stdata)
{

    varcovtap <- new("spam",entries=setup$taps*cova,colindices=setup$ja,
    rowpointers=setup$ia,dimension=as.integer(rep(dimat,2)))
   
          #cholvctap <-spam::chol.spam(varcovtap)
           cholvctap <- spam::update.spam.chol.NgPeyton(setup$struct, varcovtap)
          logdet <- 2*c(spam::determinant(cholvctap)$modulus)

    ################################################
        skew=as.numeric(nuisance["skew"])
        delta=as.numeric(nuisance["tail"])
        Z=sinh(delta * asinh(stdata)-skew)
        C=delta*sqrt((1+Z^2)/(stdata^2+1))
    llik <- 0.5*( const + const*log(sill)/log(2*pi)+logdet - 2*sum(log(C)) +sum(Z* spam::solve.spam(cholvctap, Z)))
                    
    return(llik)
}


######### Tapering 1 normal log-likelihood 
LogNormDenTap1 <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
{

    varcovtap <- new("spam",entries=cova*setup$taps,colindices=setup$ja,
    rowpointers=setup$ia,dimension=as.integer(rep(dimat,2)))
    cholvctap <- spam::update.spam.chol.NgPeyton(setup$struct, varcovtap)
    #cholvctap <-spam::chol.spam(varcovtap)
    logdet <- c(spam::determinant(cholvctap)$modulus)
    lliktap= 0.5*(const+2*logdet+sum(stdata* spam::solve.spam(cholvctap, stdata)))
    return(lliktap)
}

    
######### Tapering 2 log-likelihood for bivariate multivariate normal density:
    LogNormDenTap_biv <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        lliktap <- 1.0e8
        # Computes the vector of the correlations:
        varcovtap <- try(new("spam",entries=cova*setup$taps,colindices=setup$ja,
        rowpointers=setup$ia,dimension=as.integer(rep(dimat,2))),silent=TRUE)
        cholvctap <- try(spam::chol.spam(varcovtap),silent=TRUE)
        if(inherits(cholvctap,"try-error")) {return(lliktap)}
        logdet <- c(spam::determinant(cholvctap)$modulus)
        inv <- spam::solve.spam(cholvctap)
        slot(varcovtap,"entries") <- inv[setup$idx]*setup$taps
        lliktap= 0.5*(const+2*logdet+drop(t(stdata)%*%varcovtap%*%stdata))
        return(lliktap)
    }
######### Tapering 1 log-likelihood for bivariate multivariate normal density:
    LogNormDenTap1_biv <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        lliktap <- 1.0e8
        # Computes the vector of the correlations:
        cova[cova==(nuisance['sill'])] <- nuisance['sill']+nuisance['nugget']
        varcovtap <- new("spam",entries=cova*setup$taps,colindices=setup$ja,
        rowpointers=setup$ia,dimension=as.integer(rep(dimat,2)))
        cholvctap <- spam::update.spam.chol.NgPeyton(setup$struct, varcovtap)
        if(inherits(cholvctap, "try-error"))  {return(lliktap)}
        logdet <- c(spam::determinant.spam.chol.NgPeyton(cholvctap)$modulus)
        lliktap <- 0.5*(const+2*logdet+sum(stdata* spam::solve.spam(cholvctap, stdata)))
        return(lliktap)
    }

######### Standard log-likelihood function for multivariate normal density with sparse alg matrices
    LogNormDenStand_spam <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        llik <- 1.0e8
        # Computes the covariance matrix:
        varcov <- (nuisance['sill'])*ident
        varcov[lower.tri(varcov)] <- cova
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- cova
        mcov=spam::as.spam(varcov)#,silent=TRUE);if(class(mcov)=="try-error")   return(llik)
        cholS <- spam::chol.spam(mcov)# ,silent=TRUE);if(class(cholS)=="try-error")   return(llik)         
        llik=0.5*( const+2*c(spam::determinant.spam.chol.NgPeyton(cholS)$modulus) 
                                    + sum(stdata* spam::solve.spam(cholS,stdata)))
        return(llik)
    }    

    LogNormDenStand22 <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        llik <- 1.0e8
    
        # decomposition of the covariance matrix:
       decompvarcov <- MatDecomp(cova, mdecomp)
    if(is.logical(decompvarcov) && !decompvarcov) {
        return(llik)
    }
    logdetvarcov <- MatLogDet(decompvarcov, mdecomp)
    if(is.na(logdetvarcov)) {
        return(llik)
    }
         llik <- 0.5*(const+logdetvarcov+  sum((forwardsolve(decompvarcov, stdata, transpose = FALSE))^2))
        return(llik)
    }
    ######### CVV mdecomp:
    CVV <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        llik <- 1.0e8
        # Computes the covariance matrix:
        varcov <- (nuisance['sill'])*ident
        varcov[lower.tri(varcov)] <- cova
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- cova      
        # decomposition of the covariance matrix:
        decompvarcov <- MatDecomp(varcov,mdecomp)
        if(is.logical(decompvarcov)) return(llik)  
        invarcov <- MatInv(varcov)
        D=diag(1/diag(invarcov))
        M=crossprod(invarcov,D);C=tcrossprod(M,M)
        llik <- mean(crossprod(t(crossprod(stdata,C)),stdata))
        return(llik)
    }

    ######### CVV mdecomp in the bivariate case
CVV_biv <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        llik <- 1.0e8
        # Computes the covariance matrix:
        ident[lower.tri(ident,diag=T)] <- cova
        ident <- t(ident)
        ident[lower.tri(ident,diag=T)] <- cova
        # decomposition of the covariance matrix:
        decompvarcov <- MatDecomp(ident,mdecomp)
        if(is.logical(decompvarcov)) return(llik)
        invarcov <- MatInv(ident)
        D=diag(1/diag(invarcov))
        M=crossprod(invarcov,D);C=tcrossprod(M,M)
        llik <- mean(crossprod(t(crossprod(stdata,C)),stdata))
        return(llik)
    }

 ######## Standard log-likelihood function for log gaussian random fields   


 #LogNormDenStand_LG <- function(const, cova, ident, dimat, mdecomp, nuisance, det, sill, setup, stdata) {
  #  llik <- 1.0e8  
   # varcov <- sill * ident
   # varcov[lower.tri(varcov)] <- cova
   # varcov <- t(varcov)
   # varcov[lower.tri(varcov)] <- cova  
   # decompvarcov <- MatDecomp(varcov, mdecomp)
   # if (is.logical(decompvarcov)) return(llik)  # Fallimento
   # logdetvarcov <- MatLogDet(decompvarcov, mdecomp)
  #  quad_form <-sum((forwardsolve(decompvarcov, stdata, transpose = FALSE))^2)
  #  llik <- 0.5 * (const + logdetvarcov + 2 * det + quad_form)
  #  return(llik)
#}







       
       
######### Standard log-likelihood function for multivariate bivariate normal density:
 LogNormDenStand_biv <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        llik <- 1.0e8
        # Computes the covariance matrix:
        ident[lower.tri(ident,diag=T)] <- cova
        ident <- t(ident)
        ident[lower.tri(ident,diag=T)] <- cova
        # decomposition of the covariance matrix:
       decompvarcov <- MatDecomp(ident, mdecomp)
    if(is.logical(decompvarcov) && !decompvarcov) {
        return(llik)
    }
    logdetvarcov <- MatLogDet(decompvarcov, mdecomp)
    if(is.na(logdetvarcov)) {
        return(llik)
    }
        #invarcov <- MatInv(decompvarcov,mdecomp)
        #llik <- 0.5*(const+logdetvarcov+crossprod(t(crossprod(stdata,invarcov)),stdata))
        llik <- 0.5*(const+logdetvarcov+  sum((forwardsolve(decompvarcov, stdata, transpose = FALSE))^2))
        return(llik)

        
    }
    ######### Standard log-likelihood function for multivariate bivariate normal density with sparse alg matrices
 LogNormDenStand_biv_spam <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        llik <- 1.0e8
        # Computes the covariance matrix:
        ident[lower.tri(ident,diag=T)] <- cova
        ident <- t(ident)
        ident[lower.tri(ident,diag=T)] <- cova
        # decomposition of the covariance matrix:
        mcov=try(spam::as.spam(ident),silent=TRUE)
        if(inherits(mcov,"try-error")){return(llik)}
        cholS <- try(chol(mcov) ,silent=T);if(inherits(cholS,"try-error")){ return(llik)}
        llik=0.5*( const+2*c(spam::determinant.spam.chol.NgPeyton(cholS)$modulus) 
                                    + sum(stdata* spam::solve.spam(cholS,stdata)))
         return(llik)
    }
    ### END Defining the objective functions
################################################################################################
################################################################################################
   # Call to the objective functions:
        loglik_miss_skewT<- function(param,const,coordx,coordy,coordz,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM,aniso,namesaniso)
    {
        llik <- 1.0e8
        names(param) <- namesparam
        # Set the parameter vector:
        pram <- c(param, fixed)
        paramcorr <- pram[namescorr]
        nuisance <- pram[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
                  mm=as.numeric(nuisance[sel])
          Mean=c(X%*%mm)
        if(!is.null(MM)) Mean=MM

      if(aniso){     ### anisotropy
            anisopar<-pram[namesaniso]
            coords1=GeoAniso (cbind(coordx,coordy,coordz), anisopars=as.numeric(anisopar))       
             if(ncol(coords1)==2) {coordx=coords1[,1];coordy=coords1[,2];coordz=NULL}
             if(ncol(coords1)==3) {coordx=coords1[,1];coordy=coords1[,2];coordz=coords1[,3]}
         }
   
        # Computes the vector of the correlations:
        corr=matr(corrmat,corr,coordx,coordy,coordz,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius)
        nu=1/nuisance['df']; skew2=nuisance['skew']^2
        #if(nu<2||abs(nuisance['skew'])>1)  return(llik)
        w=sqrt(1-skew2);
        KK=2*skew2/pi
        D1=(nu-1)/2; D2=nu/2
        CorSkew<-(2*skew2/(pi*w^2+skew2*(pi-2)))*(sqrt(1-corr^2)+corr*asin(corr)-1)+w^2*corr/(w^2+skew2*(1-2/pi))   
        corr3<-(pi*(nu-2)*gamma(D1)^2/(2*(pi*gamma(D2)^2-skew2*(nu-2)*gamma(D1)^2)))*(Re(hypergeo::hypergeo(0.5,0.5,D2,corr^2))*((1-KK)*CorSkew+KK)-KK)
        cova <- corr3*nuisance['sill'] *(1-nuisance['nugget'])
       #nuisance['nugget']=0
      loglik_u <- do.call(what="LogNormDenStand",args=list(stdata=(data-c(Mean)),const=const,cova=cova,dimat=dimat,ident=ident,
            mdecomp=mdecomp,nuisance=nuisance,setup=setup))
        return(loglik_u)
      }
################################################################################################
    # Call to the objective functions:
    loglik_miss_T <- function(param,const,coordx,coordy,coordz,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM,aniso,namesaniso)
    {
        llik <- 1.0e8
        names(param) <- namesparam
        # Set the parameter vector:
        pram <- c(param, fixed)
        paramcorr <- pram[namescorr]
        nuisance <- pram[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
           mm=as.numeric(nuisance[sel])
                 Mean=c(X%*%mm)
        if(!is.null(MM)) Mean=MM

     if(aniso){     ### anisotropy
            anisopar<-pram[namesaniso]
            coords1=GeoAniso (cbind(coordx,coordy,coordz), anisopars=as.numeric(anisopar))       
             if(ncol(coords1)==2) {coordx=coords1[,1];coordy=coords1[,2];coordz=NULL}
             if(ncol(coords1)==3) {coordx=coords1[,1];coordy=coords1[,2];coordz=coords1[,3]}
         }
      
        # Computes the vector of the correlations:
        corr=matr(corrmat,corr,coordx,coordy,coordz,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius)
        df=1/nuisance['df']
        #if(df<2)  return(llik)
         #if(df<170) corr=(df-2)*gamma((df-1)/2)^2/(2*gamma(df/2)^2)* corr *Re(hypergeo::hypergeo(0.5,0.5,df/2,corr^2)) 
         #else      
    corr=exp(log(df-2)+2*lgamma(0.5*(df-1))-(log(2)+2*lgamma(df/2))+log(Re(hypergeo::hypergeo(0.5,0.5, df/2,corr^2)))+log(corr))
        #if(is.nan(corr[1])||nuisance['sill']<0||nuisance['nugget']<0||nuisance['nugget']>1) return(llik)
    cova <- corr*(nuisance['sill'])*(1-nuisance['nugget'])
        #nuisance['nugget']=0
      loglik_u <- do.call(what="LogNormDenStand",args=list(stdata=(data-c(Mean)),const=const,cova=cova,dimat=dimat,ident=ident,
            mdecomp=mdecomp,nuisance=nuisance,setup=setup))
        return(loglik_u)
      }
       # Call to the objective functions:
    loglik_miss_Poisgamma <- function(param,const,coordx,coordy,coordz,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM,aniso,namesaniso)
    {
        llik <- 1.0e8
        names(param) <- namesparam
        # Set the parameter vector:
        pram <- c(param, fixed)
        paramcorr <- pram[namescorr]
        nuisance <- pram[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
        mm=as.numeric(nuisance[sel])
        Mean=c(X%*%mm)
        nuisance=nuisance[!sel]
        
        nuisance=nuisance[order(names(nuisance))]
        if(!is.null(MM)) Mean=MM

 if(aniso){     ### anisotropy
            anisopar<-pram[namesaniso]
            coords1=GeoAniso (cbind(coordx,coordy,coordz), anisopars=as.numeric(anisopar))       
             if(ncol(coords1)==2) {coordx=coords1[,1];coordy=coords1[,2];coordz=NULL}
             if(ncol(coords1)==3) {coordx=coords1[,1];coordy=coords1[,2];coordz=coords1[,3]}
         }
        # Computes the vector of the correlations:
        mu=Mean
        
        model=46
        corr=matr2(corrmat,corr,coordx,coordy,coordz,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius,model,mu)
        cova <-  ident
        cova[lower.tri(cova)] <- corr   
        cova <- t(cova)
        cova[lower.tri(cova)] <- corr 

        ff=exp(mu)
        diag(cova)=ff*(1+ff/as.numeric(nuisance["shape"]))
    if(nuisance['nugget']<0||nuisance['nugget']>1) return(llik)
        
      loglik_u <- do.call(what="LogNormDenStand22",args=list(stdata=data-c(ff),
           const=const,cova=cova,dimat=dimat,ident=ident,
            mdecomp=mdecomp,nuisance=nuisance,setup=setup))
        return(loglik_u)
    
      }
################################################################################################
    # Call to the objective functions:
    loglik_miss_Pois <- function(param,const,coordx,coordy,coordz,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM,aniso,namesaniso)
    {
        llik <- 1.0e8
        names(param) <- namesparam
        # Set the parameter vector:
        pram <- c(param, fixed)
        paramcorr <- pram[namescorr]
        nuisance <- pram[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
           mm=as.numeric(nuisance[sel])
           Mean=c(X%*%mm)
        if(!is.null(MM)) Mean=MM
          nuisance=nuisance[-sel]
        

    if(aniso){     ### anisotropy
            anisopar<-pram[namesaniso]
            coords1=GeoAniso (cbind(coordx,coordy,coordz), anisopars=as.numeric(anisopar))       
             if(ncol(coords1)==2) {coordx=coords1[,1];coordy=coords1[,2];coordz=NULL}
             if(ncol(coords1)==3) {coordx=coords1[,1];coordy=coords1[,2];coordz=coords1[,3]}
         }
      
        # Computes the vector of the correlations:
        mu=Mean
        model=30
     
        corr=matr2(corrmat,corr,coordx,coordy,coordz,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius,model,mu)
        cova <-  ident
        cova[lower.tri(cova)] <- corr   
        cova <- t(cova)
        cova[lower.tri(cova)] <- corr 
        diag(cova)=exp(mu)

    if(nuisance['nugget']<0||nuisance['nugget']>1) return(llik)
        #nuisance['nugget']=0
      loglik_u <- do.call(what="LogNormDenStand22",args=list(stdata=data-c(exp(mu)),
           const=const,cova=cova,dimat=dimat,ident=ident,
            mdecomp=mdecomp,nuisance=nuisance,setup=setup))
        return(loglik_u)
    
      }



##################################################################################
################## Tukey hh  ####################################################
##################################################################################
LogNormDenStand_Tukey2H <- function(const, cova, dimat, mdecomp, nuisance, sill, setup, stdata) {
    llik <- 1.0e8 # Valore di default in caso di errore
    delta1 <- as.numeric(nuisance["tail1"])
    delta2 <- as.numeric(nuisance["tail2"])
    
    if(delta1 < 0 || delta1 >= 0.5 || delta2 < 0 || delta2 >= 0.5 || sill <= 0) return(llik) 

    decompvarcov <- MatDecomp(cova, mdecomp)
    if(is.logical(decompvarcov) && !decompvarcov) {
        return(llik)
    }
    logdetvarcov <- MatLogDet(decompvarcov, mdecomp)
    if(is.na(logdetvarcov)) {
        return(llik)
    }
    
    # 3. Partizionamento dati positivi/negativi
    pos_idx <- which(stdata >= 0)
    neg_idx <- which(stdata < 0)
    

    if(length(pos_idx) > 0) {
        st_pos <- stdata[pos_idx]
        W_pos <- VGAM::lambertW(delta1 * st_pos^2)
        g_pos <- sqrt(W_pos / delta1)
        jac_pos <- g_pos / (st_pos * (1 + W_pos))
    }
    
    if(length(neg_idx) > 0) {
        st_neg <- stdata[neg_idx]
        W_neg <- VGAM::lambertW(delta2 * st_neg^2)
        g_neg <- sqrt(W_neg / delta2)
        jac_neg <- -g_neg / (st_neg * (1 + W_neg))  
    }
    
    tau_inv <- numeric(dimat)
    if(length(pos_idx) > 0) tau_inv[pos_idx] <- g_pos
    if(length(neg_idx) > 0) tau_inv[neg_idx] <- -g_neg  
    
    # 7. Calcolo termini della log-verosimiglianza
    quad_form_calc <- try(forwardsolve(decompvarcov, tau_inv, transpose = FALSE), silent = TRUE)
    if(inherits(quad_form_calc, "try-error")) {
        return(llik)
    }
    
    quad_form <- sum(quad_form_calc^2)
    log_jacobian <- if(length(pos_idx) > 0 && length(neg_idx) > 0) {
        sum(log(jac_pos)) + sum(log(abs(jac_neg)))
    } else if(length(pos_idx) > 0) {
        sum(log(jac_pos))
    } else {
        sum(log(abs(jac_neg)))
    }
    
    # Verifica finale dei calcoli
    if(is.na(quad_form) || is.infinite(quad_form) || 
       is.na(log_jacobian) || is.infinite(log_jacobian)) {
        return(llik)
    }
     print(logdetvarcov);
    llik <- 0.5 * (const + dimat * log(sill) + logdetvarcov + quad_form) - log_jacobian
    
    return(llik)
}
################################################################################################
loglik_tukey2h <- function(param, const, coordx, coordy, coordz, coordt, corr, corrmat, corrmodel, 
                          data, dimat, fixed, fname, grid, ident, mdecomp, model, namescorr, 
                          namesnuis, namesparam, radius, setup, X, ns, NS, MM, aniso, namesaniso) {

    llik <- 1.0e8
    names(param) <- namesparam
    pram <- c(param, fixed)
    paramcorr <- pram[namescorr]
    nuisance <- pram[namesnuis]
    sel <- substr(names(nuisance), 1, 4) == "mean"
    mm <- as.numeric(nuisance[sel])
    Mean <- c(X %*% mm)
    if(!is.null(MM)) Mean <- MM
    if(aniso) {
        anisopar <- pram[namesaniso]
        coords1 <- GeoAniso(cbind(coordx, coordy, coordz), anisopars = as.numeric(anisopar))       
        if(ncol(coords1) == 2) {
            coordx <- coords1[,1]
            coordy <- coords1[,2]
            coordz <- NULL
        }
        if(ncol(coords1) == 3) {
            coordx <- coords1[,1]
            coordy <- coords1[,2]
            coordz <- coords1[,3]
        }
    }
    sill <- as.numeric(nuisance['sill'])
    nugget <- as.numeric(nuisance['nugget'])
    tail1 <- as.numeric(nuisance['tail1'])
    tail2 <- as.numeric(nuisance['tail2'])
    if(tail1 < 0 || tail1 >= 0.5 || tail2 < 0 || tail2 >= 0.5 || nugget < 0 || nugget >= 1 || sill <= 0) {
        return(llik)
    }
     nuisance_temp=nuisance
        nuisance_temp['sill']=1
        nuisance_temp['nugget']=0
    corr <- matr(corrmat, corr, coordx, coordy, coordz, coordt, corrmodel, nuisance_temp, paramcorr, ns, NS, radius)
      if (is.nan(corr[1])) return(llik)
    corrmat_full =build_correlation_matrix(corr,ident)
    cova <- (1 - nugget) * corrmat_full + nugget * ident
    stdata <- (data - c(Mean)) / sqrt(sill)
    
 
    loglik_u <- LogNormDenStand_Tukey2H(const = const, cova = cova, dimat = dimat,
                                      mdecomp = mdecomp, nuisance = nuisance,
                                      sill = sill, setup = setup, stdata = stdata)
    return(loglik_u)
}
##################################################################################
################## Tukey h  ####################################################
##################################################################################

LogNormDenStand_TukeyH<- function(const, cova, ident, dimat, mdecomp, delta, sill, setup, stdata) {
llik <- 1.0e8
  # Decomposizione della matrice di covarianza
   decompvarcov <- MatDecomp(cova, mdecomp)
    if(is.logical(decompvarcov) && !decompvarcov) {
        return(llik)
    }
    logdetvarcov <- MatLogDet(decompvarcov, mdecomp)
    if(is.na(logdetvarcov)) {
        return(llik)
    }


  # Inversa della trasformazione Tukey-h (usando Lambert W)
  Wval <- VGAM::lambertW(delta * stdata^2)
  IL <- sign(stdata) * sqrt(Wval / delta)
  IW <- 1 / (stdata * (1 + Wval))  # Derivata della trasformazione inversa
  quad_form <- sum((forwardsolve(decompvarcov, IL, transpose = FALSE))^2)
  log_jacobian <- -2 * sum(log(abs(IL * IW)))
  # Calcolo della -logverosimiglianza
  llik <- 0.5 * (
    const +               # n * log(2π)
    dimat * log(sill) +       # n * log(σ²)
    logdetvarcov +        # log|Σ|
    quad_form +           # termine quadratico
    log_jacobian          # cambio di variabile
  )
  return(llik)
}
loglik_tukeyh <- function(param,const,coordx,coordy,coordz,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM,aniso,namesaniso)
    {
   llik <- 1.0e8
        names(param) <- namesparam
        # Set the parameter vector:
        pram <- c(param, fixed)
        paramcorr <- pram[namescorr]
        nuisance <- pram[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
          mm=as.numeric(nuisance[sel])
               Mean=c(X%*%mm)
        if(!is.null(MM)) Mean=MM

    if(aniso){     ### anisotropy
            anisopar<-pram[namesaniso]
            coords1=GeoAniso (cbind(coordx,coordy,coordz), anisopars=as.numeric(anisopar))       
             if(ncol(coords1)==2) {coordx=coords1[,1];coordy=coords1[,2];coordz=NULL}
             if(ncol(coords1)==3) {coordx=coords1[,1];coordy=coords1[,2];coordz=coords1[,3]}
         }
      
        # Computes the vector of the correlations:
        sill=as.numeric(nuisance['sill']); 
        tail=as.numeric(nuisance['tail'])
        nugget=as.numeric(nuisance['nugget'])
        if(tail<0||tail>0.5||nugget<0||nugget>=1||sill<0) return(llik)
        nuisance_temp=nuisance
        nuisance_temp['sill']=1
        nuisance_temp['nugget']=0
         corr=matr(corrmat,corr,coordx,coordy,coordz,coordt,corrmodel,nuisance_temp,paramcorr,ns,NS,radius)
           if (is.nan(corr[1])) return(llik)
    corrmat_full =build_correlation_matrix(corr,ident)
      cova<-(1 - nugget) * corrmat_full + nugget * ident

        loglik_u <- do.call(what="LogNormDenStand_TukeyH",
            args=list(stdata=((data-c(Mean))/sqrt(sill)),const=const,cova=cova,dimat=dimat,ident=ident,
            mdecomp=mdecomp,delta=tail,sill=sill,setup=setup))
        return(loglik_u)
      }
##################################################################################
################## Log Gauss  ####################################################
##################################################################################
LogNormDenStand_LG <- function(const, cova, ident, dimat, mdecomp, nuisance, 
                               sill, setup, stdata, V_data, mu_s) {
llik <- 1.0e8

 decompvarcov <- MatDecomp(cova, mdecomp)
    if(is.logical(decompvarcov) && !decompvarcov) {
        return(llik)
    }
    logdetvarcov <- MatLogDet(decompvarcov, mdecomp)
    if(is.na(logdetvarcov)) {
        return(llik)
    }
  y <- forwardsolve(decompvarcov, stdata, transpose = FALSE, upper.tri = FALSE)
    quadterm <- sum(y^2)
  sigma <- sqrt(sill)
  jacobian_term <- -dimat * log(sigma) - sum(log(V_data)) - sum(log(mu_s))
  llik <- -0.5 * (const + logdetvarcov + quadterm) + jacobian_term
  return(-llik)
}

loglik_loggauss <- function(param, const, coordx, coordy, coordz, coordt, corr, corrmat, corrmodel,
                            data, dimat, fixed, fname, grid, ident, mdecomp, model,
                            namescorr, namesnuis, namesparam, radius, setup, X,
                            ns, NS, MM, aniso, namesaniso) {
  
 llik <- 1.0e8
  names(param) <- namesparam
  pram <- c(param, fixed)
  # Parametri
  paramcorr <- pram[namescorr]
  nuisance   <- pram[namesnuis]
  # Termini di media
  sel <- substr(names(nuisance), 1, 4) == "mean"
  mm <- as.numeric(nuisance[sel])
  Mean <- as.vector(X %*% mm)
  if (!is.null(MM)) Mean <- MM
  # Anisotropia
  if (aniso) {
    anisopar <- as.numeric(pram[namesaniso])
    coords1 <- GeoAniso(cbind(coordx, coordy, coordz), anisopars = anisopar)
    coordx <- coords1[, 1]
    coordy <- coords1[, 2]
    if (ncol(coords1) == 3) coordz <- coords1[, 3] else coordz <- NULL
  }
  # Parametri del modello log-Gaussiano
  sill <- as.numeric(nuisance["sill"])    
  nugget <- as.numeric(nuisance["nugget"]) 
  # Controlli di validità
  if (sill <= 0 || is.na(nugget) || nugget < 0 || nugget >= 1) return(llik)
  # Media del modello: μ(s) = exp(X(s)^T β)
  mu_s <- exp(Mean)
  V_data <- data / mu_s
  sigma <- sqrt(sill)
  Z_data <- (log(V_data) + sill/2) / sigma
  nuisance_temp <- nuisance
  nuisance_temp["sill"] <- 1  # Standardizzazione per correlazione
  corr <- matr(corrmat, corr, coordx, coordy, coordz, coordt, corrmodel,
               nuisance_temp, paramcorr, ns, NS, radius)
  if (is.nan(corr[1])) return(llik)
       corrmat_full =build_correlation_matrix(corr,ident)
      cova<-(1 - nugget) * corrmat_full + nugget * ident
  loglik_gaussian <- LogNormDenStand_LG(
    stdata = Z_data,const = const,cova = cova,dimat = dimat,ident = ident,
    mdecomp = mdecomp,nuisance = nuisance,sill = sill,setup = setup,V_data = V_data,
    mu_s = mu_s)
  
  return(loglik_gaussian)
}
##################################################################################
################## Sin H  ########################################################
##################################################################################
    LogNormDenStand_SH <- function(const,cova,mdecomp,nuisance,sill,setup,stdata)
    {
llik <- 1.0e8
         decompvarcov <- MatDecomp(cova, mdecomp)
    if(is.logical(decompvarcov) && !decompvarcov) {
        return(llik)
    }
    logdetvarcov <- MatLogDet(decompvarcov, mdecomp)
    if(is.na(logdetvarcov)) {
        return(llik)
    }
        skew=as.numeric(nuisance["skew"])
        delta=as.numeric(nuisance["tail"])
        Z=sinh(delta * asinh(stdata)-skew)
        C=delta*sqrt((1+Z^2)/(stdata^2+1))
        llik <- 0.5*( const + const*log(sill)/log(2*pi)
                      +logdetvarcov - 2*sum(log(C))
                      +sum((forwardsolve(decompvarcov, Z, transpose = FALSE))^2))
        return(llik)
    }
loglik_sh <- function(param,const,coordx,coordy,coordz,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM,aniso,namesaniso)
    {

      llik <- 1.0e8
  names(param) <- namesparam
  pram <- c(param, fixed)
  paramcorr <- pram[namescorr]
  nuisance   <- pram[namesnuis]
  # Termini di media
  sel <- substr(names(nuisance), 1, 4) == "mean"
  mm <- as.numeric(nuisance[sel])
  Mean <- as.vector(X %*% mm)
  if (!is.null(MM)) Mean <- MM
  # Anisotropia
  if (aniso) {
    anisopar <- as.numeric(pram[namesaniso])
    coords1 <- GeoAniso(cbind(coordx, coordy, coordz), anisopars = anisopar)
    coordx <- coords1[, 1]
    coordy <- coords1[, 2]
    if (ncol(coords1) == 3) coordz <- coords1[, 3] else coordz <- NULL
  }
         # Computes the vector of the correlations:
        sill=as.numeric(nuisance['sill']);
        nugget=as.numeric(nuisance['nugget'])
        tail=as.numeric(nuisance['tail'])
        nuisance_temp=nuisance
        nuisance_temp['sill']=1
        if(tail<=0||sill<=0||nugget<0||nugget>=1) return(llik)
        corr=matr(corrmat,corr,coordx,coordy,coordz,coordt,corrmodel,nuisance_temp,paramcorr,ns,NS,radius)
     if (is.nan(corr[1])) return(llik)
    corrmat_full =build_correlation_matrix(corr,ident)
      cova<-(1 - nugget) * corrmat_full + nugget * ident

        loglik_u <- do.call(what="LogNormDenStand_SH",
            args=list(stdata=((data-c(Mean))/(sqrt(sill))),const=const,cova=cova,
            mdecomp=mdecomp,nuisance=nuisance,sill=sill,setup=setup))
        return(loglik_u)
      }
####################################################################################
################## Gaussian ########################################################
####################################################################################

LogNormDenStand <- function(const, cova, dimat, mdecomp, stdata) {
   llik <- 1.0e8
     decompvarcov <- MatDecomp(cova, mdecomp)
    if(is.logical(decompvarcov) && !decompvarcov) {
        return(llik)
    }
    logdetvarcov <- MatLogDet(decompvarcov, mdecomp)
    if(is.na(logdetvarcov)) {
        return(llik)
    }
    quad <- sum((forwardsolve(decompvarcov, stdata, transpose = FALSE))^2)
    llik <- 0.5 * (const + logdetvarcov + quad)
    return(llik)
}

loglik <- function(param, const, coordx, coordy, coordz, coordt, corr, corrmat, corrmodel, data, dimat, fixed, fname,
                   grid, ident, mdecomp, model, namescorr, namesnuis, namesparam, radius, setup, X, ns, NS, MM, aniso, namesaniso) 
    {
    llik <- 1.0e8
    names(param) <- namesparam
    pram <- c(param, fixed)
    paramcorr <- pram[namescorr]
    nuisance <- pram[namesnuis]
    
    # Estrazione della media
    sel <- substr(names(nuisance), 1, 4) == "mean"
    mm <- as.numeric(nuisance[sel])
    Mean <- as.vector(X %*% mm)
    if (!is.null(MM)) Mean <- MM
    
    # Gestione anisotropia se richiesta
    if (aniso) {
        anisopar <- pram[namesaniso]
        coords1 <- GeoAniso(cbind(coordx, coordy, coordz), anisopars = as.numeric(anisopar))
        if (ncol(coords1) == 2) { coordx <- coords1[,1]; coordy <- coords1[,2]; coordz <- NULL }
        if (ncol(coords1) == 3) { coordx <- coords1[,1]; coordy <- coords1[,2]; coordz <- coords1[,3] }
    }
    
    sill <- as.numeric(nuisance['sill'])
    nugget <- as.numeric(nuisance['nugget'])
    
    # Controllo parametri
    if (is.nan(sill) || is.nan(nugget) || sill <= 0 || nugget < 0 || nugget >= 1) return(llik)
    #nuisance['sill']=1; nuisance['nugget']=0;
    corr <- matr(corrmat, corr, coordx, coordy, coordz, coordt, corrmodel, nuisance, paramcorr, ns, NS, radius)
    if (is.nan(corr[1])) return(llik)
    corrmat_full <- build_correlation_matrix(corr, ident)
    cova <- sill*((1 - nugget) * corrmat_full + nugget * ident)  # matrice di correlazione
    
    stdata <- (data - c(Mean))
    # Calcolo della log-likelihood
    loglik_u <- do.call(
        what = fname,
        args = list(
            stdata = stdata,
            const = const,
            cova = cova,
            dimat = dimat,
            mdecomp = mdecomp
        )
    )
    return(loglik_u)
}

  #####################################################    

# loglikvecchia <- function(param,vecchia.approx,data,fixed,dimat,
 #                   model,namescorr,namesnuis,namesparam,X,MM,aniso,namesaniso)
 #   {
 #       llik <- 1.0e8
 #       names(param) <- namesparam
 #       # Set the parameter vector:
 #       pram <- c(param, fixed)
 #       paramcorr <- pram[namescorr]
 #       nuisance <- pram[namesnuis]
 #       sel=substr(names(nuisance),1,4)=="mean"
 #       mm=as.numeric(nuisance[sel])
 #       Mean=c(X%*%mm)
 #       if(!is.null(MM)) Mean=MM
 #  
 #         nuggets=as.numeric(nuisance['nugget'])+ 0.00001
 #       data=c(data-X%*%mm)
 #       ppar=as.numeric(c(nuisance['sill'], paramcorr[1], paramcorr[2]))
 #   if(ppar[2]<0|| ppar[3]<0||nuisance['sill']<0||nuisance['nugget']<0||nuisance['nugget']>1){return(llik)}
 #       loglik_u=GPvecchia::vecchia_likelihood(data,vecchia.approx,
 #                        covparms=ppar,nuggets=nuggets,covmodel ="matern")
 #       return(-loglik_u)
 #     }

      
################################################################################################                     
     loglik_biv <- function(param,const,coordx,coordy,coordz,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM,aniso,namesaniso)
      {

        # Set the parameter vector:
        names(param) <- namesparam
        pram <- c(param, fixed)
        paramcorr <- pram[namescorr]
        nuisance <- pram[namesnuis]

        sel1=substr(names(nuisance),1,6)=="mean_1"
        mm1=as.numeric(nuisance[sel1])
        sel2=substr(names(nuisance),1,6)=="mean_2"
        mm2=as.numeric(nuisance[sel2])
       
        X1=as.matrix(X[1:ns[1],]);X2=as.matrix(X[(ns[1]+1):(ns[2]+ns[1]),]); 
        mm=as.double(c(X1%*%mm1,X2%*%mm2))
        # Standardizes the data:
         stdata <- data-mm 
      # Computes the vector of the correlations
         corr=matr(corrmat,corr,coordx,coordy,coordz,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius)
      # Computes the log-likelihood
       loglik_b <- do.call(what=fname,args=list(stdata=stdata,const=const,cova=corr,ident=ident,dimat=dimat,
            mdecomp=mdecomp,nuisance=nuisance,setup=setup))

        return(loglik_b)
      }


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
    ### START the main code of the function:

    spacetime_dyn=FALSE; NS=0;fname=NULL
    if(!is.null(coordx_dyn)) spacetime_dyn=TRUE
    if(grid)     {a=expand.grid(coordx,coordy);coordx=a[,1];coordy=a[,2]; }

    ####################################
    if(!spacetime_dyn) dimat <- numcoord*numtime# length of data
    if(spacetime_dyn)  dimat =sum(ns)

    if(is.null(dim(X))) {
    if(!bivariate) X=as.matrix(rep(1,dimat))  # matrix of covariates
    if( bivariate) {X=as.matrix(rep(1,ns[1]+ns[2]));# X=rbind(X,X)
        }}
    else{
    if(!bivariate) num_betas=ncol(X)  
    if( bivariate) num_betas=c(ncol(X),ncol(X)) }
     
    corrmat<-"CorrelationMat"# set the type of correlation matrix   

    if(model==36||model==47) corrmat<-"CorrelationMat_dis"# set the type of correlation matrix 
    if(spacetime)  { corrmat<-"CorrelationMat_st_dyn"
                     if(model==36||model==47)  corrmat="CorrelationMat_st_dyn_dis"
                    } 
    if(bivariate)  corrmat<-"CorrelationMat_biv_dyn"  
     if(spacetime||bivariate){
          NS=cumsum(ns);
          if(spacetime_dyn){ data=t(unlist(data));NS=c(0,NS)[-(length(ns)+1)]}
          else {data=matrix(t(data),nrow=1);NS=rep(0,numtime)}  
    }       
    ####################################
    dd=0
     
    if(bivariate)  dd=dimat
    numpairstot <- dimat*(dimat-1)/2+dd
    const<-dimat*log(2*pi)# set the likelihood constant

     if(is.null(neighb))
     {
     corr<-double(numpairstot)# initialize the correlation
     ident <- diag(dimat)# set the identity matrix}
     }

 #########################################################################  
 ######################################################################### 
if(model==1||model==20){  ## gaussian case
    lname <- 'loglik'
    if(!is.null(neighb))  lname <-'loglikvecchia'
    if(bivariate)  {lname <- 'loglik_biv'}
  
    # detects the type of likelihood:
    if(type==3){
        if(bivariate) fname<-"LogNormDenRestr_biv"
        else          fname<-"LogNormDenRestr"
        const <- const-log(2*pi)}
    if(type==4) { if(bivariate)   fname <- 'LogNormDenStand_biv'
                  else            fname <- 'LogNormDenStand'
                }
    if(type==5||type==6){  #tapering
           corrmat<-"CorrelationMat_tap"
           if(spacetime) corrmat<-"CorrelationMat_st_tap"
           if(bivariate) {if(type==5) fname <- 'LogNormDenTap_biv'
                       if(type==6) fname <- 'LogNormDenTap1_biv'
                      corrmat<-"CorrelationMat_biv_tap" }
            else          {if(type==5) fname <- 'LogNormDenTap'
                       if(type==6)  fname <- 'LogNormDenTap1'
                                   
                      }             
        corr <- double(numpairs)
        tapcorr <- double(numpairs)
        tcor=matr(corrmat,tapcorr,coordx,coordy,coordz,coordt,setup$tapmodel,c(0,0,1),1,ns,NS,radius)
        tape <- new("spam",entries=tcor,colindices=setup$ja,rowpointers=setup$ia,dimension=as.integer(rep(dimat,2)))
        setup$struct <- try(spam::chol.spam(tape,silent=TRUE))
        setup$taps<-tcor
        }
        if(type==8)
        { if(bivariate) fname<-"CVV_biv"
          else          fname<-"CVV"
       }
     if(sparse) fname <- paste(fname,"_spam",sep="")
 }
 ############################################################################
############################################################################
hessian=FALSE
 if(model==20){   ## SAS case
     lname <- 'loglik_sh'
    if(bivariate)  {lname <- 'loglik_biv_sh'}
    fname='LogNormDenStand_SH'
    if(type==6) fname <- 'LogshDenTap1'

}
 if(model==34){   ## Tukeyh case
     lname <- 'loglik_tukeyh'
    if(bivariate)  {lname <- 'loglik_biv_tukeyh'}

}

 if(model==40){   ## Tukey2h case
     lname <- 'loglik_tukey2h'
    if(bivariate)  {lname <- 'loglik_biv_tukey2h'}
   
}
 if(model==47){   ## gaussian misspecified poissongamma
     lname <- 'loglik_miss_Poisgamma'
    if(bivariate)  {lname <- 'loglik_biv_miss_Poisgamma'}

}
 if(model==35){   ## gaussian misspecified t
     lname <- 'loglik_miss_T'
    if(bivariate)  {lname <- 'loglik_biv_miss_T'}

}
 if(model==36){   ## Poisson misspecified t
     lname <- 'loglik_miss_Pois'
    if(bivariate)  {lname <- 'loglik_biv_miss_Pois'}

}
 if(model==37){   ## gaussian misspecified skewt
     lname <- 'loglik_miss_skewT'
    if(bivariate)  {lname <- 'loglik_biv_miss_skewT'}

}
 if(model==22){   ## loggaussian  case
     lname <- 'loglik_loggauss'
    if(bivariate)  {lname <- 'loglik_biv_loggauss'}

}

#############
 if(type!=5&&type!=6){ corrmat <- paste(corrmat,"2",sep="") }
#############
##################
namesaniso=c("angle","ratio")
##################



if(is.null(coordz)) coordz=double(length(coordx))

if(!onlyvar){   # performing optimization

 

    maxit=10000
    # Optimize the log-likelihood:
   if(length(param)==1)
        {
         optimizer="optimize"         
  Likelihood <- optimize(f=eval(as.name(lname)),const=const,coordx=coordx,coordy=coordy,coordz=coordz, coordt=coordt,corr=corr,corrmat=corrmat,
                          corrmodel=corrmodel,data=t(data),dimat=dimat,fixed=fixed,
                          fname=fname,grid=grid,ident=ident,lower=lower,maximum = FALSE,mdecomp=mdecomp,
                          model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,
                          upper=upper,radius=radius,setup=setup,X=X,ns=ns,NS=NS,MM=MM,aniso=aniso,namesaniso=namesaniso)
        }
  if(length(param)>1)
        {
  
  #### vecchia  approxxx
    #if(!is.null(neighb)) { 
     #       locs=cbind(coordx,coordy)
    #
     #       vecchia.approx=GPvecchia::vecchia_specify(locs,m=neighb,ordering="maxmin")#,conditioning="NN",cond.yz="SGV")
     #     if(optimizer=="nlminb")
     #        Likelihood <- nlminb(objective=eval(as.name(lname)),start=param,vecchia.approx=vecchia.approx,
     #                        control = list( iter.max=100000),dimat=dimat,
     #                    lower=lower,upper=upper, hessian=hessian,
     #                      data=t(data),fixed=fixed,
     #                     model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,X=X,MM=MM,aniso=aniso,namesaniso=namesaniso) 
     #       if(optimizer=="Nelder-Mead")
     #             Likelihood <- optim(param,eval(as.name(lname)),vecchia.approx=vecchia.approx,
     #                        control=list(reltol=1e-14, maxit=maxit),dimat=dimat,hessian=hessian, data=t(data),fixed=fixed,
     #                        model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,X=X,MM=MM,aniso=aniso,namesaniso=namesaniso)
     #                  }
     #   else{  ### no vecchia
if(optimizer=='L-BFGS-B')
                        Likelihood <- optim(param,fn=eval(as.name(lname)),const=const,coordx=coordx,coordy=coordy,coordz=coordz,coordt=coordt,corr=corr,corrmat=corrmat,
                          corrmodel=corrmodel,control=list(
                          pgtol=1e-14,maxit=maxit),data=t(data),dimat=dimat,fixed=fixed,
                          fname=fname,grid=grid,ident=ident,lower=lower,mdecomp=mdecomp,method=optimizer,
                          model=model,namescorr=namescorr,hessian=hessian,
                          namesnuis=namesnuis,upper=upper,namesparam=namesparam,radius=radius,setup=setup,X=X,ns=ns,NS=NS,MM=MM,aniso=aniso,namesaniso=namesaniso) 


  #if(optimizer=='BFGS')
   #                   Likelihood <- optim(param,eval(as.name(lname)),const=const,coordx=coordx,coordy=coordy,coordt=coordt,corr=corr,corrmat=corrmat,
    #                      corrmodel=corrmodel,control=list(
     #                     pgtol=1e-14,maxit=maxit),data=t(data),dimat=dimat,
      #                   fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,method=optimizer,
       #                   model=model,namescorr=namescorr,hessian=hessian,
       #                   namesnuis=namesnuis,namesparam=namesparam,radius=radius,setup=setup,X=X,ns=ns,NS=NS,MM=MM,aniso=aniso,namesaniso=namesaniso)
  if(optimizer=='Nelder-Mead'||optimizer=='SANN'||optimizer=="BFGS")
                      Likelihood <- optim(param,eval(as.name(lname)),const=const,coordx=coordx,coordy=coordy,coordz=coordz,coordt=coordt,corr=corr,corrmat=corrmat,
                          corrmodel=corrmodel,control=list(
                          reltol=1e-14, maxit=maxit),data=t(data),dimat=dimat,
                          fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,method=optimizer,
                          model=model,namescorr=namescorr,hessian=hessian,
                          namesnuis=namesnuis,namesparam=namesparam,radius=radius,setup=setup,X=X,ns=ns,NS=NS,MM=MM,aniso=aniso,namesaniso=namesaniso)
  if(optimizer=='nlm')
                      Likelihood <- nlm(eval(as.name(lname)),param,const=const,coordx=coordx,coordy=coordy,coordz=coordz,coordt=coordt,corr=corr,corrmat=corrmat,
                          corrmodel=corrmodel,data=t(data),dimat=dimat,fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,hessian=hessian,
                          model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,radius=radius,setup=setup,iterlim = maxit,X=X,ns=ns,NS=NS,MM=MM,aniso=aniso,namesaniso=namesaniso)
    if(optimizer=='bobyqa')
                      Likelihood <-minqa::bobyqa(fn = get(lname),par=param, 
                                     control = list(maxfun=maxit),
                              lower=lower,upper=upper, #hessian=hessian,
                          const=const,coordx=coordx,coordy=coordy,coordz=coordz,coordt=coordt,corr=corr,corrmat=corrmat,
                          corrmodel=corrmodel,data=t(data),dimat=dimat,fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,
                          model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,radius=radius,
                          setup=setup,X=X,ns=ns,NS=NS,MM=MM,aniso=aniso,namesaniso=namesaniso)
                           

  if(optimizer=='nlminb')
                      Likelihood <-nlminb(objective=eval(as.name(lname)),start=param,
                             control = list( iter.max=100000),
                          lower=lower,upper=upper, hessian=hessian,
                          const=const,coordx=coordx,coordy=coordy,coordz=coordz,coordt=coordt,corr=corr,corrmat=corrmat,
                          corrmodel=corrmodel,data=t(data),dimat=dimat,fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,
                          model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,radius=radius,
                          setup=setup,X=X,ns=ns,NS=NS,MM=MM,aniso=aniso,namesaniso=namesaniso)
                               
    #}
}





 if(optimizer %in% c('Nelder-Mead','L-BFGS-B','BFGS','nmk','nmkb','multiNelder-Mead','bobyqa','SANN'))
                   {names(Likelihood$par)=namesparam
                    param <- Likelihood$par
                    if(optimizer=='bobyqa')  maxfun <- -Likelihood$fval
                    else                     maxfun <- -Likelihood$value
                   Likelihood$value <- maxfun
    if(optimizer %in% c('Nelder-Mead','L-BFGS-B','BFGS','SANN'))  Likelihood$counts=as.numeric(Likelihood$counts[1])
     if(optimizer %in% c('nmk','nmkb','bobyqa'))                    Likelihood$counts=as.numeric(Likelihood$feval)
               }

    #if(optimizer=='ucminf'){
    #               names(Likelihood$par)=namesparam
    #               param <- Likelihood$par
    #               maxfun <- -Likelihood$value
     #              Likelihood$value <- maxfun
     # }
       if(optimizer %in% c('nlminb','multinlminb')){
                   names(Likelihood$par)=namesparam
                   param <- Likelihood$par
                   maxfun <- -Likelihood$objective
                   Likelihood$value <- maxfun
                   Likelihood$counts=as.numeric(Likelihood$iterations)
      }
    if(optimizer=='nlm')
           { names(Likelihood$estimate)=namesparam
             param <- Likelihood$estimate
             #names(param)<-namesparam
             maxfun <- -as.numeric(Likelihood$minimum)
             Likelihood$value <- maxfun
             Likelihood$param <- param
             Likelihood$counts=as.numeric(Likelihood$iterations)
          }
     if(optimizer=='optimize')  
          {param<-Likelihood$minimum
           names(param)<-namesparam
           maxfun <- -Likelihood$objective
           Likelihood$value <- maxfun
           Likelihood$param<-param}
    numparam<-length(param)# parameter size
    Likelihood$claic <- NULL; Likelihood$clbic <- NULL;
    if(type==3||type==4) 
             {Likelihood$claic <- -2*(maxfun-numparam) 
              Likelihood$clbic <- -2*maxfun+numparam*log(dimat)    }
    ### Some checks of the output from the optimization procedure:
     if(optimizer=='bobyqa'){
    if(Likelihood$ierr == 0)
      Likelihood$convergence <- 'Successful'
    else
        Likelihood$convergence <- 'Optimization may have failed'
    }

    if(optimizer=='Nelder-Mead' ||  optimizer=='L-BFGS-B'||  optimizer=='BFGS'||  optimizer=='SANN'){
    if(Likelihood$convergence == 0)
      Likelihood$convergence <- 'Successful'
    else
      if(Likelihood$convergence == 1)
        Likelihood$convergence <- 'Iteration limit reached'
      else
        Likelihood$convergence <- 'Optimization may have failed'}
    if(optimizer=='optimize'){  Likelihood$convergence <- 'Successful'}
    if(optimizer=='nmk' || optimizer=='nmkb'){
                   if(Likelihood$convergence == 0) Likelihood$convergence <- 'Successful'
                   else Likelihood$convergence <- 'Optimization may have failed'}
    if(optimizer=='nlm'){
               if(Likelihood$code == 1||Likelihood$code == 2)
               Likelihood$convergence <- 'Successful'
               else
               if(Likelihood$code == 4)
               Likelihood$convergence <- 'Iteration limit reached'
               else
               Likelihood$convergence <- 'Optimization may have failed'}
        if(optimizer=='nlminb'||optimizer=='multinlminb'){
               if(Likelihood$convergence == 0)
               Likelihood$convergence <- 'Successful'
               else
               Likelihood$convergence <- 'Optimization may have failed'}
      # if(optimizer=='ucminf'){
       #        if(Likelihood$convergence== 1||Likelihood$convergence== 2||Likelihood$convergence == 4)
        #       Likelihood$convergence <- 'Successful'
         #      else
          #     Likelihood$convergence <- 'Optimization may have failed'}
    if(maxfun==-1.0e8) Likelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
     }
     else {Likelihood=as.list(0)   # in the case of computing just the asym variance
           names(Likelihood)="value"
           maxfun=Likelihood$value
           Likelihood$param <- param
           Likelihood$claic <- NULL;Likelihood$clbic <- NULL;
           Likelihood$convergence <- 'None'
           numparam<-length(param)
            }
   


if(varest) 
  {
 #    aa=try(abs(det(Likelihood$hessian)),silent=T)
 #  if(aa<1e-08||is.null(Likelihood$hessian)||min(eigen(Likelihood$hessian)$values)<0)
 # {  

Likelihood$hessian=numDeriv::hessian(func=eval(as.name(lname)),x=Likelihood$par,method="Richardson",  const=const,coordx=coordx,coordy=coordy,coordz=coordz,
            coordt=coordt,corr=corr,corrmat=corrmat,
            corrmodel=corrmodel,data=t(data),dimat=dimat,fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,
            model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,radius=radius,setup=setup,X=X,ns=ns,NS=NS,MM=MM,aniso=aniso,namesaniso=namesaniso)

Likelihood$score=numDeriv::grad(func=eval(as.name(lname)),x=Likelihood$par,method="Richardson",  const=const,coordx=coordx,coordy=coordy,coordz=coordz,
            coordt=coordt,corr=corr,corrmat=corrmat,
            corrmodel=corrmodel,data=t(data),dimat=dimat,fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,
            model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,radius=radius,setup=setup,X=X,ns=ns,NS=NS,MM=MM,aniso=aniso,namesaniso=namesaniso)
rownames(Likelihood$hessian)=namesparam
colnames(Likelihood$hessian)=namesparam
names(Likelihood$score)=namesparam
}
#}


   if(Likelihood$convergence == 'Successful' || Likelihood$convergence =='None'||Likelihood$convergence =='Optimization may have failed'){  

    # if optimization has failed it does not compute stderr

    ### START Computing the asymptotic variance-covariance matrices:  


    if(varest){
    if((model==20||model==22||model==1||model==34||model==35||model==37)&& !(type==5||type==6))      {
       # if(is.null(Likelihood$hessian)) {Likelihood$hessian=numDeriv::hessian(func=eval(as.name(lname)),x=param,method="Richardson",
       #     const=const,coordx=coordx,coordy=coordy,coordt=coordt,corr=corr,corrmat=corrmat,
       #                   corrmodel=corrmodel,data=t(data),dimat=dimat,
       #                   fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,
       #                   model=model,namescorr=namescorr,
       #                   namesnuis=namesnuis,namesparam=namesparam,radius=radius,setup=setup,X=X,ns=ns,NS=NS)
       # rownames(Likelihood$hessian)=namesparam; colnames(Likelihood$hessian)=namesparam
        #             }

         aa=try(abs(det(Likelihood$hessian)),silent=T)
     if(inherits(aa,"try-error")) {
                        #options(warn = -1)
                        warning("Asymptotic information matrix is singular",immediate.=TRUE)
                        Likelihood$varcov <- NULL  } 
        else 
        {
        
         ii=solve(Likelihood$hessian)
         mm=min(eigen(ii)$values)
            if(mm<=0)   { 
                        warning("Asymptotic information matrix is not positive-definite")
                        Likelihood$varcov <- NULL  } 
            else       {Likelihood$varcov <-  ii
                        Likelihood$stderr <- sqrt(diag(( Likelihood$varcov))) } 
         }
     }
     else
     {

        #Checks if the resulting variance and covariance matrix:
        mm=min(eigen(Likelihood$varcov)$values)
        if(is.null(Likelihood$varcov)||mm<0){
            if(mm<0)                       warning("Asymptotic information matrix is not positive-definite")
            if(is.null(Likelihood$varcov)) warning("Asymptotic information matrix is singular")
            Likelihood$varcov <- 'none'
            Likelihood$stderr <- 'none'}
        else
        {
            dimnames(Likelihood$varcov)<-list(namesparam,namesparam)
            Likelihood$stderr<-diag(Likelihood$varcov)
        if(any(Likelihood$stderr < 0)) Likelihood$stderr <- 'none'
        else{
            Likelihood$stderr<-sqrt(Likelihood$stderr)
            names(Likelihood$stderr)<-namesparam}}
    }}}
    ### END the main code of the function:
if(varest) if(is.null(Likelihood$varcov)){Likelihood$varcov <- 'none';Likelihood$stderr <- 'none'}

    return(Likelihood)
  }
