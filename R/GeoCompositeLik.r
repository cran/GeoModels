####################################################
### File name: CompLik.r
####################################################



CompLik <- function(copula,bivariate, coordx, coordy ,coordz,coordt,coordx_dyn,corrmodel, data, distance, flagcorr, flagnuis, fixed,grid,
                           likelihood,lower, model, n, namescorr, namesnuis, namesparam,
                           numparam, numparamcorr, optimizer, onlyvar, param, spacetime, type,
                           upper, varest, weigthed, ns, X,sensitivity,MM,aniso)
  {
    ### Define the object function:
    comploglik <- function(param,coords,coordt, corrmodel, data, fixed, fan, n, namescorr, 
                              namesnuis,namesparam,namesaniso,weigthed,X,ns,NS,MM)
      {
        names(param) <- namesparam
        param <- c(param, fixed)
        paramcorr <- param[namescorr]
        nuisance <- param[namesnuis]

        #aniso<-param[namesaniso]

        sel=substr(names(nuisance),1,4)=="mean"
        mm=as.numeric(nuisance[sel])
        Mean=c(X%*%mm)
        if(!is.null(MM)) Mean=MM

       other_nuis=as.numeric(nuisance[!sel])   ## or nuis parameters (nugget sill skew df)

        result=dotCall64::.C64(as.character(fan),
        SIGNATURE = c("integer","double","double","double","double",
                         "integer","double","integer","double","double","double","double",
                          "integer","integer","integer","integer"),  
                         corrmodel ,coords[,1],coords[,2] , coordt ,  data , n , paramcorr ,  weigthed , 
                   res=dotCall64::vector_dc("double",1), Mean,0,other_nuis ,
                     ns , NS , 
         INTENT =    c("r","r","r","r","r","r","r","r","w","r", "r","r","r","r","r","r"),
         PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$res
 
         return(-result)
      }
     comploglik_biv <- function(param,coords,coordt, corrmodel, data, fixed, fan, n, namescorr, namesnuis,namesparam,namesaniso,weigthed,X,ns,NS,MM)
      {

        names(param) <- namesparam
        param <- c(param, fixed)
        paramcorr <- param[namescorr]
        nuisance <- param[namesnuis]
        
        #aniso<-param[namesaniso]

        sel1=substr(names(nuisance),1,6)=="mean_1"
        mm1=as.numeric(nuisance[sel1])
        sel2=substr(names(nuisance),1,6)=="mean_2"
        mm2=as.numeric(nuisance[sel2])
        sel=substr(names(nuisance),1,4)=="mean"
        X1=as.matrix(X[1:ns[1],]);X2=as.matrix(X[(ns[1]+1):(ns[2]+ns[1]),]); 
        other_nuis=as.numeric(nuisance[!sel]) 



        result=dotCall64::.C64(as.character(fan),
        SIGNATURE = c("integer","double","double","double","double","double",
                         "integer","double","integer","double","double","double","double",
                          "integer","integer","integer","integer"),  
                         corrmodel ,coordx,coordy , coordz,coordt ,  data , n , paramcorr ,  weigthed , 
                   res=dotCall64::vector_dc("double",1), c(X1%*%mm1,X2%*%mm2),0,other_nuis ,
                     ns , NS ,
         INTENT =    c("r","r","r","r","r","r","r","r","w","r", "r","r","r","r","r","r","r"),
         PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$res
        return(-result)
      }
   ############################################################################################################################################
   ############################################################################################################################################
   ############################################################################################################################################
 
    numcoord=length(coordx);numtime=1;spacetime_dyn=FALSE
    if(spacetime) numtime=length(coordt)
    if(bivariate) numtime=2
    if(!is.null(coordx_dyn)) spacetime_dyn=TRUE
    if(!spacetime_dyn)  dimat <- numcoord*numtime#
    if(spacetime_dyn)  dimat <- sum(ns)

    NS=cumsum(ns)
    if(is.null(dim(X)))
    {
       X=as.matrix(rep(1,dimat))
       if((spacetime||bivariate)&& spacetime_dyn)  X=as.matrix(rep(1,NS[numtime]))
    }
    else(if(bivariate) X=rbind(X,X))

  
    fname <- NULL; hessian <- FALSE 
    namesaniso=c("angle","ratio")

    if(all(model==1,likelihood==1,type==2)) fname <- 'Comp_Cond_Gauss'
    if(all(model==1,likelihood==3,type==1)) fname <- 'Comp_Diff_Gauss'
    if(all(model==1,likelihood==3,type==2)) {fname <- 'Comp_Pair_Gauss'
                                              if(varest ) hessian <- TRUE}
    if(all(model==2,likelihood==3,type==2)){ fname <- 'Comp_Pair_BinomGauss'
                                             if(varest ) hessian <- TRUE}
    if(all(model==11,likelihood==3,type==2)){ fname <- 'Comp_Pair_BinomGauss'
                                              if(varest ) hessian <- TRUE}
    if(all(model==19,likelihood==3,type==2)){ namesnuis=c(namesnuis,"z")
                                              fixed<- c(fixed, list(z=min(n)))
                                              fname <- 'Comp_Pair_Binom2Gauss'
                                              if(varest ) hessian <- TRUE}
    if(all(model==14,likelihood==3,type==2)){ fname <- 'Comp_Pair_BinomnegGauss'
                                              if(varest ) hessian <- TRUE}
    if(all(model==16,likelihood==3,type==2)){ fname <- 'Comp_Pair_BinomnegGauss'
                                              if(varest ) hessian <- TRUE}
    if(all(model==15,likelihood==3,type==2)){ fname <- 'Comp_Pair_PoisbinGauss'
                                              if(varest ) hessian <- TRUE}
    if(all(model==17,likelihood==3,type==2)){ fname <- 'Comp_Pair_PoisbinnegGauss'
                                              if(varest ) hessian <- TRUE}    
    if(all(model==13,likelihood==3,type==2)){ fname <- 'Comp_Pair_WrapGauss'
                                              if(varest ) hessian <- TRUE}
    if(all(model==10,likelihood==3,type==2)){ fname <- 'Comp_Pair_SkewGauss'
                                              if(varest ) hessian <- TRUE}
     if(all(model==21,likelihood==3,type==2)){ fname <- 'Comp_Pair_Gamma'
                                              if(varest ) hessian <- TRUE}
      if(all(model==21,likelihood==1,type==2)){ fname <- 'Comp_Cond_Gamma'
                                              if(varest ) hessian <- TRUE}
    if(all(model==33,likelihood==3,type==2)){ fname <- 'Comp_Pair_Kumaraswamy'
                                              if(varest ) hessian <- TRUE}
      if(all(model==42,likelihood==3,type==2)){ fname <- 'Comp_Pair_Kumaraswamy2'
                                              if(varest ) hessian <- TRUE}
   if(all(model==26,likelihood==3,type==2)){ fname <- 'Comp_Pair_Weibull'
                                              if(varest ) hessian <- TRUE}    
   if(all(model==26,likelihood==1,type==2)){ fname <- 'Comp_Cond_Weibull'
                                              if(varest ) hessian <- TRUE}    
    if(all(model==28,likelihood==3,type==2)){ fname <- 'Comp_Pair_Beta'
                                              if(varest ) hessian <- TRUE}                                   
    if(all(model==24,likelihood==3,type==2)){ fname <- 'Comp_Pair_LogLogistic'
                                              if(varest ) hessian <- TRUE}     
    if(all(model==25,likelihood==3,type==2)){ fname <- 'Comp_Pair_Logistic'
                                              if(varest ) hessian <- TRUE}                                                                              
    if(all(model==23,likelihood==3,type==2)){ fname <- 'Comp_Pair_2Gamma'
                                              if(varest ) hessian <- TRUE}                                        
    if(all(model==22,likelihood==3,type==2)){ fname <- 'Comp_Pair_LogGauss';
                                              if(varest ) hessian <- TRUE}
    if(all(model==18,likelihood==3,type==2)){ fname <- 'Comp_Pair_SkewTGauss'
                                              if(varest ) hessian <- TRUE}  
    if(all(model==27,likelihood==3,type==2)){ fname <- 'Comp_Pair_TWOPIECET'
                                              if(varest ) hessian <- TRUE} 
    if(all(model==39,likelihood==3,type==2)){ fname <- 'Comp_Pair_TWOPIECEBIMODAL'
                                              if(varest ) hessian <- TRUE} 
    if(all(model==29,likelihood==3,type==2)){ fname <- 'Comp_Pair_TWOPIECEGauss'
                                              if(varest ) hessian <- TRUE} 
    if(all(model==31,likelihood==3,type==2)){ fname <- 'Comp_Pair_BinomTWOPIECEGauss'
                                              if(varest ) hessian <- TRUE} 
    if(all(model==32,likelihood==3,type==2)){ fname <- 'Comp_Pair_BinomnegTWOPIECEGauss'
                                              if(varest ) hessian <- TRUE}
    if(all(model==12,likelihood==3,type==2)){ fname <- 'Comp_Pair_T'
                                              if(varest ) hessian <- TRUE}  
    if(all(model==34,likelihood==3,type==2)){ fname <- 'Comp_Pair_Tukeyh' 
                                              if(varest ) hessian <- TRUE} 
    if(all(model==41,likelihood==3,type==2)){ fname <- 'Comp_Pair_Gauss_misp_Tukeygh' 
                                              if(varest ) hessian <- TRUE} 
    if(all(model==40,likelihood==3,type==2)){ fname <- 'Comp_Pair_Tukeyhh' 
                                              if(varest ) hessian <- TRUE} 
    if(all(model==36,likelihood==3,type==2)){ fname <- 'Comp_Pair_Gauss_misp_Pois'
                                              if(varest ) hessian <- TRUE}
    if(all(model==35,likelihood==3,type==2)){ fname <- 'Comp_Pair_Gauss_misp_T'
                                              if(varest ) hessian <- TRUE}
    if(all(model==37,likelihood==3,type==2)){ fname <- 'Comp_Pair_Gauss_misp_SkewT'
                                              if(varest ) hessian <- TRUE}
    if(all(model==20,likelihood==3,type==2)){ fname <- 'Comp_Pair_SinhGauss'
                                              if(varest ) hessian <- TRUE}
    if(all(model==38,likelihood==3,type==2)){ fname <- 'Comp_Pair_TWOPIECETukeyh'
                                              if(varest ) hessian <- TRUE} 
    if(all(model==30,likelihood==3,type==2)){ fname <- 'Comp_Pair_Pois'
                                              if(varest ) hessian <- TRUE}
    if(all(model==46,likelihood==3,type==2)){ fname <- 'Comp_Pair_PoisGamma'
                                              if(varest ) hessian <- TRUE}
    if(all(model==46,likelihood==3,type==2)){ fname <- 'Comp_Pair_PoisGamma'
                                              if(varest ) hessian <- TRUE}
    if(all(model==43,likelihood==3,type==2)){ fname <- 'Comp_Pair_PoisZIP'
                                              if(varest ) hessian <- TRUE}
    if(all(model==44,likelihood==3,type==2)){ fname <- 'Comp_Pair_Gauss_misp_PoisZIP'
                                              if(varest ) hessian <- TRUE}
    if(all(model==45,likelihood==3,type==2)){ fname <- 'Comp_Pair_BinomnegGaussZINB'
                                              if(varest ) hessian <- TRUE}
    if(sensitivity) hessian=TRUE
    if(spacetime) fname <- paste(fname,"_st",sep="")
    if(bivariate) fname <- paste(fname,"_biv",sep="")
    fname <- paste(fname,"2",sep="")

    if(!is.null(copula))
    {
        if(copula=="Gaussian") fname <- paste(fname,"GCop",sep="")
        if(copula=="Clayton")     fname <- paste(fname,"BCop",sep="")
    }

    if(aniso) fname <- paste(fname,"_aniso",sep="")

   # path.parent <- getwd()
  
    
   
    
    if((spacetime||bivariate)&&(!spacetime_dyn)){
                                  data=c(t(data))
                                  coordx=rep(coordx,numtime);coordy=rep(coordy,numtime);
                                  if(!is.null) coordz=rep(coordz,numtime);
                }
                if((spacetime||bivariate)&&(spacetime_dyn)) data=unlist(data)          
    
   if(spacetime||bivariate)   NS=c(0,NS)[-(length(ns)+1)]



coords=cbind(coordx,coordy,coordz)

   if(!onlyvar){
  ##############################.  spatial or space time ############################################
   if(!bivariate)           {
    if(length(param)==1) {
         optimizer="optimize"  
     CompLikelihood <- optimize(f=comploglik, coords=coords, coordt=coordt,corrmodel=corrmodel, 
                              data=data, fixed=fixed, fan=fname,  lower=lower, n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso, maximum = FALSE,
                              upper=upper,weigthed=weigthed,X=X,ns=ns,NS=NS,MM=MM)
          }
   if(length(param)>1) {
    if(optimizer=='L-BFGS-B')
      CompLikelihood <- optim(par=param,fn=comploglik, 
                              control=list(factr=1e-10,pgtol=1e-14, maxit=100000), 
                              coords=coords, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
                              fan=fname, lower=lower, method='L-BFGS-B',n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso, 
                              upper=upper,weigthed=weigthed,X=X,ns=ns,NS=NS, hessian=FALSE,MM=MM)
   
    if(optimizer=='BFGS') 
        CompLikelihood <- optim(par=param, fn=comploglik,  coords=coords, coordt=coordt,corrmodel=corrmodel, 
                           control=list(factr=1e-10,
                             reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='BFGS',n=n,namescorr=namescorr,
                                  namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso, weigthed=weigthed,X=X,ns=ns,NS=NS,MM=MM)

      if(optimizer=='Nelder-Mead')
        CompLikelihood <- optim(par=param, fn=comploglik,  coords=coords, coordt=coordt,corrmodel=corrmodel, 
          control=list( reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='Nelder-Mead',n=n,namescorr=namescorr,
                                  namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso, weigthed=weigthed,X=X,ns=ns,NS=NS,MM=MM)
    
 #    if(optimizer=='lbfgsb3')
 #      CompLikelihood <-lbfgsb3::lbfgsb3(prm=param,fn=comploglik,gr=NULL,lower=lower,upper=upper,
 #                             #control=list(factr=1e-10,pgtol=1e-14, maxit=100000), 
 #                             coords=coords, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
 #                             fan=fname,n=n,
 #                             namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, 
 #                             weigthed=weigthed,X=X,ns=ns,NS=NS)

    #if(optimizer=='nmk')
     # CompLikelihood <-dfoptim::nmk(par=param, fn=comploglik, control = list(maxfeval=100000,tol=1e-10),
      #                     coords=coords, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed, fan=fname,
       #                    n=n,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso, weigthed=weigthed,X=X,ns=ns,NS=NS,MM=MM)
    #if(optimizer=='nmkb')
    #{
     # CompLikelihood <-dfoptim::nmkb(par=param, fn=comploglik, control = list(maxfeval=100000,tol=1e-10),
      #                   lower=lower,upper=upper,
       #                  coords=coords, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed, fan=fname,
        #                 n=n,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso, weigthed=weigthed,X=X,ns=ns,NS=NS,MM=MM)
    #}
    if(optimizer=='nlm')
    CompLikelihood <- nlm(f=comploglik,p=param,steptol = 1e-4, coords=coords, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
                               fan=fname,hessian=FALSE,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso,  
                               iterlim=100000, weigthed=weigthed,X=X,ns=ns,NS=NS,MM=MM)
  
    if(optimizer=='nlminb')
     CompLikelihood <-nlminb(objective=comploglik,start=param,coords=coords, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
                                control = list( iter.max=100000),
                              lower=lower,upper=upper,
                               fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso, 
                               weigthed=weigthed,X=X,ns=ns,NS=NS,MM=MM)

   #      if(optimizer=='multinlminb'){
    #   CompLikelihood <- mcGlobaloptim::multiStartoptim(objectivefn=comploglik,
     #         coords=coords, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
      #                         fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso, 
       #                        weigthed=weigthed,X=X,ns=ns,NS=NS,MM=MM,
        #                            lower=lower,upper=upper,method = "nlminb", nbtrials = 500, 
         #                     control = list( iter.max=100000),
          #                 typerunif = "sobol")
           #                    }
     #if(optimizer=='multiNelder-Mead'){
      # CompLikelihood <- mcGlobaloptim::multiStartoptim(objectivefn=comploglik,
       #  coords=coords, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
        #                       fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso, MM=MM, 
         #                      weigthed=weigthed,X=X,ns=ns,NS=NS,lower=lower,upper=upper,
         # method = "Nelder-Mead", nbtrials = 500, 
          #                    control=list( reltol=1e-14, maxit=100000),
           #                typerunif = "sobol")
 # }

                                   
    }}
     ############################## bivariate  ############################################                           
    if(bivariate)           {
      
    
     if(length(param)==1)
        {
         optimizer="optimize" 
       CompLikelihood <- optimize(f=comploglik_biv, coords=coords, coordt=coordt,corrmodel=corrmodel, 
                              data=data, fixed=fixed,fan=fname,  lower=lower,n=n,MM=MM,
                              namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,namesaniso=namesaniso, maximum = FALSE,
                              upper=upper,weigthed=weigthed,X=X,ns=ns,NS=NS)}
      if(length(param)>1) {   
    if(optimizer=='L-BFGS-B'){
      CompLikelihood <- optim(param,comploglik_biv, control=list(pgtol=1e-14, maxit=100000),
                              method='L-BFGS-B',hessian=FALSE,lower=lower, upper=upper,MM=MM,
                              coords=coords, coordt=coordt,corrmodel=corrmodel, 
                              data=data, fixed=fixed,fan=fname,n=n,
                              namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,namesaniso=namesaniso, 
                             weigthed=weigthed,X=X,ns=ns,NS=NS )}
  
      if(optimizer=='BFGS')
      CompLikelihood <- optim(param,comploglik_biv, coords=coords, coordt=coordt,corrmodel=corrmodel, control=list(
                              reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='BFGS',n=n,MM=MM,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam ,namesaniso=namesaniso, weigthed=weigthed,X=X,ns=ns,NS=NS)
   if(optimizer=='Nelder-Mead')
      CompLikelihood <- optim(param,comploglik_biv, coords=coords, coordt=coordt,corrmodel=corrmodel, control=list(
                              reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='Nelder-Mead',n=n,MM=MM,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam ,namesaniso=namesaniso, weigthed=weigthed,X=X,ns=ns,NS=NS)
       if(optimizer=='nlm') 
        CompLikelihood <- nlm( f=comploglik_biv,p=param, coords=coords, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
                               fan=fname,hessian=FALSE,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso, MM=MM, 
                               weigthed=weigthed,X=X,ns=ns,NS=NS)
 
    if(optimizer=='nlminb') 
        CompLikelihood <- nlminb( objective=comploglik_biv,start=param, 
                                     control = list( iter.max=100000),
                              lower=lower,upper=upper,
                               coords=coords, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
                               fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso, MM=MM,
                               weigthed=weigthed,X=X,ns=ns,NS=NS)

    
   }}                    
      ########################################################################################   
      ########################################################################################
    # check the optimisation outcome
      if(optimizer=='Nelder-Mead'||optimizer=='multiNelder-Mead'){
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

    if(optimizer=='L-BFGS-B'||optimizer=='BFGS'||optimizer=='lbfgsb3c'){
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

    if(optimizer=='nlminb'||optimizer=='multinlminb' ){
     
        CompLikelihood$par <- CompLikelihood$par
        names(CompLikelihood$par)<- namesparam
        CompLikelihood$value <- -CompLikelihood$objective
        if(CompLikelihood$convergence == 0) { CompLikelihood$convergence <- 'Successful' }
        else {CompLikelihood$convergence <- "Optimization may have failed" }
        #if(CompLikelihood$objective==-1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
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
          if(!bivariate)CompLikelihood$value = - comploglik(param=CompLikelihood$par ,  coords=coords, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed, fan=fname,
                             n=n,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso, weigthed=weigthed,X=X,ns=ns,NS=NS,MM=MM)
          else CompLikelihood$value = -comploglik_biv(param=CompLikelihood$par ,  coords=coords, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed, fan=fname,
                             n=n,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso, weigthed=weigthed,X=X,ns=ns,NS=NS,MM=MM)
          

          if(hessian) 
          {
               if(!bivariate)  
                CompLikelihood$hessian=numDeriv::hessian(func=comploglik,x=param,method="Richardson",  coords=coords, coordt=coordt,corrmodel=corrmodel, 
                              data=data, fixed=fixed,fan=fname,n=n,
                              namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,namesaniso=namesaniso, 
                              weigthed=weigthed,X=X,ns=ns,NS=NS,MM=MM)
               if(bivariate)  
               CompLikelihood$hessian=numDeriv::hessian(func=comploglik_biv,x=param,method="Richardson",coords=coords, coordt=coordt,corrmodel=corrmodel, 
                             data=data, fixed=fixed,fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso, 
                             weigthed=weigthed,X=X,ns=ns,NS=NS,MM=MM)
               rownames(CompLikelihood$hessian)=namesparam
               colnames(CompLikelihood$hessian)=namesparam
          }
  }

#####################################
#if((sensitivity||varest)&&(is.null(CompLikelihood$hessian)||min(eigen(CompLikelihood$hessian)$values)<0))
if((sensitivity||varest))
  {
if(!bivariate)  

CompLikelihood$hessian=numDeriv::hessian(func=comploglik,x=CompLikelihood$par,method="Richardson",  coords=coords, coordt=coordt,corrmodel=corrmodel, 
                              data=data, fixed=fixed,fan=fname,n=n,
                              namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,namesaniso=namesaniso, 
                              weigthed=weigthed,X=X,ns=ns,NS=NS,MM=MM)
if(bivariate)  
CompLikelihood$hessian=numDeriv::hessian(func=comploglik_biv,x=CompLikelihood$par,method="Richardson",coords=coords, coordt=coordt,corrmodel=corrmodel, 
                             data=data, fixed=fixed,fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso, 
                               weigthed=weigthed,X=X,ns=ns,NS=NS,MM=MM)
rownames(CompLikelihood$hessian)=namesparam
colnames(CompLikelihood$hessian)=namesparam
  }


####################################
#if( (CompLikelihood$convergence!='Successful')||CompLikelihood$value==-1e+15)  print("Optimization failed: try with other starting values ")
if(hessian) CompLikelihood$sensmat=CompLikelihood$hessian

return(CompLikelihood)
  }

