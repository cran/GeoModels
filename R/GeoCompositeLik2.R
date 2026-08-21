####################################################
### File name: CompLik2.R
####################################################

### Optim call for Composite log-likelihood maximization


CompLik2 <- function(copula,bivariate, coordx, coordy ,coordz,coordt,coordx_dyn,corrmodel, data, distance, flagcorr, flagnuis, fixed,grid,
                           likelihood,lower, model, n, namescorr, namesnuis, namesparam,
                           numparam, numparamcorr, optimizer, onlyvar, param, spacetime, type,
                           upper, varest,  weigthed, ns, X,sensitivity,
                           colidx,rowidx,neighb,MM,aniso,score)
  {

comploglik2 <- function(param, colidx, rowidx, corrmodel, coords, data1, data2, fixed, fan, n, namescorr,
                              namesnuis, namesparam, namesaniso, weigthed, X, MM, aniso, type_cop, cond_pair) {

    names(param) <- namesparam
    param <- c(param, fixed)

    ## Guard copula and marginal parameter domains at every objective
    ## evaluation. This is essential for unconstrained optimizers such as
    ## Nelder-Mead, especially for the skew-Gaussian eta parameter.
    if (type_cop > 0L) {
        bad <- FALSE
        get1 <- function(nm) {
            if (!(nm %in% names(param))) return(NA_real_)
            suppressWarnings(as.numeric(param[[nm]]))
        }

        nug <- get1("nugget")
        if (!is.finite(nug) || nug < 0 || nug >= 1) bad <- TRUE

        if ("sill" %in% names(param)) {
            sv <- get1("sill")
            if (!is.finite(sv) || sv <= 0) bad <- TRUE
        }
        if (model == 12L) {
            dfv <- get1("df")
            if (!is.finite(dfv) || dfv <= 0) bad <- TRUE
        }
        if (model %in% c(21L, 26L, 50L, 42L)) {
            sh <- get1("shape")
            if (!is.finite(sh) || sh <= 0) bad <- TRUE
        }
        if (model %in% c(28L, 33L)) {
            sh1 <- get1("shape1"); sh2 <- get1("shape2")
            if (!is.finite(sh1) || !is.finite(sh2) || sh1 <= 0 || sh2 <= 0) bad <- TRUE
        }
        if (model %in% c(28L, 33L, 50L, 42L)) {
            lo <- get1("min"); hi <- get1("max")
            if (!is.finite(lo) || !is.finite(hi) || lo >= hi) bad <- TRUE
        }
        if (model == 59L) {
            sk <- get1("skew")
            if (!is.finite(sk) || sk <= 0 || sk >= 1) bad <- TRUE
        }

        if (type_cop == 2L) {
            nu_cop <- get1("nu")
            if (!is.finite(nu_cop) || nu_cop <= 0 ||
                abs(nu_cop - round(nu_cop)) > sqrt(.Machine$double.eps)) bad <- TRUE
        }
        if (type_cop == 3L) {
            eta_cop <- get1("nu")
            if (!is.finite(eta_cop) || abs(eta_cop) >= 1) bad <- TRUE
        }

        if (bad) return(1e20)
    }

    paramcorr <- param[namescorr]
    nuisance <- param[namesnuis]
    sel <- sel_mean

    if (is.null(MM)) {
        beta <- .GeoMean_beta(nuisance, ncol(X))
        Mean <- as.numeric(X %*% beta)
    } else {
        Mean <- as.numeric(MM)
    }
    other_nuis <- nuisance[!sel]

    if (aniso) {
        coords1 <- GeoAniso(coords, anisopars = param[namesaniso])
        c1 <- as.vector(t(coords1[colidx, , drop = FALSE]))
        c2 <- as.vector(t(coords1[rowidx, , drop = FALSE]))

        result <- dotCall64::.C64(
            as.character(fan),
            SIGNATURE = c("integer", "double", "double", "double", "double",
                         "integer", "integer", "double", "integer", "double",
                         "double", "double", "double", "integer", "integer"),
            corrmodel, c1, c2, data1, data2,
            n_col, n_row, paramcorr, weigthed,
            res = dotCall64::vector_dc("double", 1),
            Mean[colidx], Mean[rowidx], other_nuis, type_cop, cond_pair,
            INTENT = c(rep("r", 9), "w", rep("r", 5)),
            PACKAGE = 'GeoModels',
            VERBOSE = 0,
            NAOK = TRUE
        )$res

    } else {
        result <- dotCall64::.C64(
            as.character(fan),
            SIGNATURE = c("integer", "double", "double", "integer", "integer",
                         "double", "integer", "double", "double", "double",
                         "double", "integer", "integer"),
            corrmodel, data1, data2,
            n_col, n_row, paramcorr, weigthed,
            res = dotCall64::vector_dc("double", 1),
            Mean[colidx], Mean[rowidx], other_nuis, type_cop, cond_pair,
            INTENT = c(rep("r", 7), "w", rep("r", 5)),
            PACKAGE = 'GeoModels',
            VERBOSE = 0,
            NAOK = TRUE
        )$res
    }

    ## Optimizers work on the average negative composite log-likelihood.
    ## obj_scale is defined after the pair selection and is constant in theta.
    return(-result / obj_scale)
}
##################################################################
##################################################################
comploglik_biv2 <- function(param,colidx,rowidx, corrmodel, coords,data1,data2,fixed, fan, n,
                          namescorr, namesnuis,namesparam,namesaniso,weigthed,X,MM,aniso,type_cop,cond_pair)
      {

        names(param) <- namesparam
        param <- c(param, fixed)
        paramcorr <- param[namescorr]
        nuisance <- param[namesnuis]
        anisopar <- param[namesaniso]

        sel1=substr(names(nuisance),1,6)=="mean_1"
        mm1=as.numeric(nuisance[sel1])
        sel2=substr(names(nuisance),1,6)=="mean_2"
        mm2=as.numeric(nuisance[sel2])
        sel=substr(names(nuisance),1,4)=="mean"
        X1=as.matrix(X[1:ns[1],]);X2=as.matrix(X[(ns[1]+1):(ns[2]+ns[1]),]);
        other_nuis=as.numeric(nuisance[!sel])
        Mean=c(X1%*%mm1,X2%*%mm2)
        res=double(1)

        result=dotCall64::.C64(as.character(fan),
          SIGNATURE = c("integer","double","double", "integer","double","integer","double","double","double","double"),
                        corrmodel,data1, data2, n,paramcorr,weigthed, res=res,Mean[colidx],Mean[rowidx],other_nuis,
          INTENT =    c(rep("r",6),"w", rep("r", 3)),
             PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$res
        ## Optimizers work on the average negative composite log-likelihood.
        ## obj_scale is defined after the pair selection and is constant in theta.
        return(-result / obj_scale)
      }

   ##################################################################################################
   ############### starting function ################################################################
   ##################################################################################################

    numcoord=length(coordx);numtime=1;spacetime_dyn=FALSE
    if(spacetime) numtime=length(coordt)
    if(bivariate) numtime=2
    if(!is.null(coordx_dyn)) spacetime_dyn=TRUE

    dimat=numcoord*numtime
    if((spacetime||bivariate)) dimat <- sum(ns)
    NS=cumsum(ns)
    if(is.null(dim(X))){X=as.matrix(rep(1,dimat))}

    fname <- NULL; hessian <- FALSE
    namesaniso=c("angle","ratio")

    if(length(n)==1) n=rep(n,dimat)
####################### conditional ##############################################
    if(all(model==1, likelihood==1,type==2)) fname <- 'Comp_Cond_Gauss'
    if(all(model==21,likelihood==1,type==2)) fname <- 'Comp_Cond_Gamma'
    if(all(model==26,likelihood==1,type==2)) fname <- 'Comp_Cond_Weibull'
    if(all(model==22,likelihood==1,type==2)) fname <- 'Comp_Cond_LogGauss'
    if(all(model==28,likelihood==1,type==2)) fname <- 'Comp_Cond_Beta'
    if(all(model==50,likelihood==1,type==2)) fname <- 'Comp_Cond_Beta2'
    if(all(model==13,likelihood==1,type==2)) fname <- 'Comp_Cond_WrapGauss'
    if(all(model==33,likelihood==1,type==2)) fname <- 'Comp_Cond_Kumaraswamy'
    if(all(model==42,likelihood==1,type==2)) fname <- 'Comp_Cond_Kumaraswamy2'
    if(all(model==20,likelihood==1,type==2)) fname <- 'Comp_Cond_SinhGauss'
    if(all(model==38,likelihood==1,type==2)) fname <- 'Comp_Cond_TWOPIECETukeyh'
    if(all(model==27,likelihood==1,type==2)) fname <- 'Comp_Cond_TWOPIECET'
    if(all(model==29,likelihood==1,type==2)) fname <- 'Comp_Cond_TWOPIECEGauss'
    if(all(model==10,likelihood==1,type==2)) fname <- 'Comp_Cond_SkewGauss'
    if(all(model==12,likelihood==1,type==2)) fname <- 'Comp_Cond_T'
    if(all(model==34,likelihood==1,type==2)) fname <- 'Comp_Cond_Tukeyh'
    if(all(model==40,likelihood==1,type==2)) fname <- 'Comp_Cond_Tukeyhh'
    if(all(model==35,likelihood==1,type==2)) fname <- 'Comp_Cond_Gauss_misp_T'
    if(all(model==39,likelihood==1,type==2)) fname <- 'Comp_Cond_TWOPIECEBIMODAL'
    if(all(model==31,likelihood==1,type==2)) fname <- 'Comp_Cond_BinomTWOPIECEGauss'
    if(all(model==32,likelihood==1,type==2)) fname <- 'Comp_Cond_BinomnegTWOPIECEGauss'
    if(all(model==30,likelihood==1,type==2)) fname <- 'Comp_Cond_Pois'
    if(all(model==46,likelihood==1,type==2)) fname <- 'Comp_Cond_PoisGamma'
    if(all(model==36,likelihood==1,type==2)) fname <- 'Comp_Cond_Gauss_misp_Pois'
    if(all(model==2,likelihood==1,type==2))  fname <- 'Comp_Cond_BinomGauss'
    if(all(model==11,likelihood==1,type==2)&length(n)==1) fname <- 'Comp_Cond_BinomGauss'
    if(all(model==51,likelihood==1,type==2)&length(n)>=1) fname <- 'Comp_Cond_BinomNNGauss_misp'
    if(all(model==11,likelihood==1,type==2)&length(n)>1) fname <- 'Comp_Cond_BinomNNGauss'
    if(all(model==49,likelihood==1,type==2)&length(n)==1) fname <- 'Comp_Cond_BinomLogi'
    if(all(model==49,likelihood==1,type==2)&length(n)>1) fname <- 'Comp_Cond_BinomNNLogi'
    if(all(model==14,likelihood==1,type==2)) fname <- 'Comp_Cond_BinomnegGauss'
    if(all(model==16,likelihood==1,type==2)) fname <- 'Comp_Cond_BinomnegGauss'
    if(all(model==54,likelihood==1,type==2)) fname <- 'Comp_Cond_BinomnegBinary'
    if(all(model==43,likelihood==1,type==2)) fname <- 'Comp_Cond_PoisZIP'
    if(all(model==57,likelihood==1,type==2)) fname <- 'Comp_Cond_PoisGammaZIP'
    if(all(model==44,likelihood==1,type==2)) fname <- 'Comp_Cond_Gauss_misp_PoisZIP'
    if(all(model==45,likelihood==1,type==2)) fname <- 'Comp_Cond_BinomnegGaussZINB'
    if(all(model==56,likelihood==1,type==2)) fname <- 'Comp_Cond_BinomnegGaussZINB1'
    if(all(model==24,likelihood==1,type==2)) fname <- 'Comp_Cond_LogLogistic'
    if(all(model==25,likelihood==1,type==2)) fname <- 'Comp_Cond_Logistic'
    if(all(model==59,likelihood==1,type==2)) fname <- 'Comp_Cond_SkewLaplace'
    if(all(model==41,likelihood==1,type==2)) fname <- 'Comp_Cond_Gauss_misp_Tukeygh'
    if(all(model==37,likelihood==1,type==2)) fname <- 'Comp_Cond_Gauss_misp_SkewT'
    if(all(model==47,likelihood==1,type==2)) fname <- 'Comp_Cond_Gauss_misp_PoisGamma'
###################### pairwise ###############################################
    if(all(model==1,likelihood==3,type==2)) fname <- 'Comp_Pair_Gauss'
    if(all(model==2,likelihood==3,type==2)) fname <- 'Comp_Pair_BinomGauss'
    if(all(model==11,likelihood==3,type==2)&length(n)==1) fname <- 'Comp_Pair_BinomGauss'
    if(all(model==51,likelihood==3,type==2)&length(n)>=1) fname <- 'Comp_Pair_BinomNNGauss_misp'
    if(all(model==11,likelihood==3,type==2)&length(n)>1) fname <- 'Comp_Pair_BinomNNGauss'
    if(all(model==49,likelihood==3,type==2)&length(n)==1) fname <- 'Comp_Pair_BinomLogi'
    if(all(model==49,likelihood==3,type==2)&length(n)>1) fname <- 'Comp_Pair_BinomNNLogi'
    if(all(model==19,likelihood==3,type==2)) { namesnuis=c(namesnuis,"z")
                                              fixed<- c(fixed, list(z=min(n)))
                                              fname <- 'Comp_Pair_Binom2Gauss'}
    if(all(model==14,likelihood==3,type==2)) fname <- 'Comp_Pair_BinomnegGauss'
    if(all(model==16,likelihood==3,type==2)) fname <- 'Comp_Pair_BinomnegGauss'
    if(all(model==54,likelihood==3,type==2)) fname <- 'Comp_Pair_BinomnegBinary'
    if(all(model==15,likelihood==3,type==2)) fname <- 'Comp_Pair_PoisbinGauss'
    if(all(model==17,likelihood==3,type==2)) fname <- 'Comp_Pair_PoisbinnegGauss'
    if(all(model==13,likelihood==3,type==2)) fname <- 'Comp_Pair_WrapGauss'
    if(all(model==10,likelihood==3,type==2)) fname <- 'Comp_Pair_SkewGauss'
    if(all(model==21,likelihood==3,type==2)) fname <- 'Comp_Pair_Gamma'
    if(all(model==33,likelihood==3,type==2)) fname <- 'Comp_Pair_Kumaraswamy'
    if(all(model==42,likelihood==3,type==2)) fname <- 'Comp_Pair_Kumaraswamy2'
    if(all(model==28,likelihood==3,type==2)) fname <- 'Comp_Pair_Beta'
    if(all(model==50,likelihood==3,type==2)) fname <- 'Comp_Pair_Beta2'
    if(all(model==26,likelihood==3,type==2)) fname <- 'Comp_Pair_Weibull'
    if(all(model==24,likelihood==3,type==2)) fname <- 'Comp_Pair_LogLogistic'
    if(all(model==25,likelihood==3,type==2)) fname <- 'Comp_Pair_Logistic'
    if(all(model==59,likelihood==3,type==2)) fname <- 'Comp_Pair_SkewLaplace'
    if(all(model==23,likelihood==3,type==2)) fname <- 'Comp_Pair_2Gamma'
    if(all(model==22,likelihood==3,type==2)) fname <- 'Comp_Pair_LogGauss'
    if(all(model==18,likelihood==3,type==2)) fname <- 'Comp_Pair_SkewTGauss'
    if(all(model==27,likelihood==3,type==2)) fname <- 'Comp_Pair_TWOPIECET'
    if(all(model==39,likelihood==3,type==2)) fname <- 'Comp_Pair_TWOPIECEBIMODAL'
    if(all(model==29,likelihood==3,type==2)) fname <- 'Comp_Pair_TWOPIECEGauss'
    if(all(model==31,likelihood==3,type==2)) fname <- 'Comp_Pair_BinomTWOPIECEGauss'
    if(all(model==32,likelihood==3,type==2)) fname <- 'Comp_Pair_BinomnegTWOPIECEGauss'
    if(all(model==12,likelihood==3,type==2)) fname <- 'Comp_Pair_T'
    if(all(model==34,likelihood==3,type==2)) fname <- 'Comp_Pair_Tukeyh'
    if(all(model==40,likelihood==3,type==2)) fname <- 'Comp_Pair_Tukeyhh'
    if(all(model==41,likelihood==3,type==2)) fname <- 'Comp_Pair_Gauss_misp_Tukeygh'
    if(all(model==36,likelihood==3,type==2)) fname <- 'Comp_Pair_Gauss_misp_Pois'
    if(all(model==35,likelihood==3,type==2)) fname <- 'Comp_Pair_Gauss_misp_T'
    if(all(model==37,likelihood==3,type==2)) fname <- 'Comp_Pair_Gauss_misp_SkewT'
    if(all(model==20,likelihood==3,type==2)) fname <- 'Comp_Pair_SinhGauss'
    if(all(model==38,likelihood==3,type==2)) fname <- 'Comp_Pair_TWOPIECETukeyh'
    if(all(model==30,likelihood==3,type==2)) fname <- 'Comp_Pair_Pois'
    if(all(model==46,likelihood==3,type==2)) fname <- 'Comp_Pair_PoisGamma'
    if(all(model==43,likelihood==3,type==2)) fname <- 'Comp_Pair_PoisZIP'
    if(all(model==57,likelihood==3,type==2)) fname <- 'Comp_Pair_PoisGammaZIP'
    if(all(model==44,likelihood==3,type==2)) fname <- 'Comp_Pair_Gauss_misp_PoisZIP'
    if(all(model==45,likelihood==3,type==2)) fname <- 'Comp_Pair_BinomnegGaussZINB'
    if(all(model==56,likelihood==3,type==2)) fname <- 'Comp_Pair_BinomnegGaussZINB1'
    if(all(model==47,likelihood==3,type==2)) fname <- 'Comp_Pair_Gauss_misp_PoisGamma'
 #############################################################################
    if(sensitivity) hessian=TRUE
    if(spacetime) fname <- paste(fname,"_st",sep="")
    if(bivariate) fname <- paste(fname,"_biv",sep="")

    type_cop=0
    cond_pair=0

  if(!is.null(copula))
    {
        fname <- paste(fname,"Cop",sep="")
        if(copula=="Gaussian")      {type_cop=1; }
        if(copula=="Clayton")       {type_cop=2; }
        if(copula=="SkewGaussian")  {type_cop=3; }
        if(all(likelihood==1,type==2)) {cond_pair=1;  fname=gsub("Cond", "Pair", fname)}
        if(all(likelihood==3,type==2)) {cond_pair=0}
    }

    fname <- paste(fname,"2mem",sep="")
    if(aniso) fname <- paste(fname,"_aniso",sep="")

     if((spacetime||bivariate)&&(!spacetime_dyn))    data=c(t(data))
     if((spacetime||bivariate)&&(spacetime_dyn))     data=unlist(data)
     if(spacetime||bivariate)   NS=c(0,NS)[-(length(ns)+1)]

###### selecting data with indexes from composite likelihood

   if(is.null(neighb)) {colidx=colidx+1; rowidx=rowidx+1}
   data1=data[colidx]; data2=data[rowidx]

   ## Count observations do not change during optimization.  Validate and
   ## convert them once, then reuse the integer cache in every Poisson-Gamma
   ## pairwise likelihood evaluation (base, ZIP, and Gaussian-misspecified).
   if(!bivariate && model %in% c(46L,47L,57L)) {
     pg_cache <- .C("SetPGPairCounts",
       as.double(data1),as.double(data2),as.integer(length(data1)),
       ok=integer(1),PACKAGE="GeoModels",NAOK=TRUE)$ok
     if(pg_cache[1L] != 1L)
       stop("Poisson-Gamma data must be non-negative integer counts or NA")
   }

   ## Constant scaling of the objective.  The optimizer and numerical score
   ## use the average pairwise objective, while CompLikelihood$value is
   ## converted back to the total composite log-likelihood below.
   obj_scale <- max(1, length(data1))

   if((model==11||model==49||model==51)&&length(n)>1) {n1=n[colidx];n2=n[rowidx];n=c(n1,n2)}

   ## Quantities below are invariant across objective evaluations.
   ## Cache them once in the CompLik2 frame; comploglik2 accesses them
   ## through lexical scoping without changing its call signature.
   if(!bivariate) {
     sel_mean <- namesnuis %in% .GeoMean_names(ncol(X))
     n_col <- n[colidx]
     n_row <- n[rowidx]
   }

##################
lname="comploglik2"

coords=cbind(coordx,coordy,coordz)

## Original optimizer controls.
## Keep these unchanged to preserve the same numerical path as the previous code.
## The changes below only reclassify optimizer warnings after convergence.
opt_reltol <- 1e-12
opt_pgtol  <- 1e-10
opt_factr  <- 1e4
opt_maxit  <- 100000L

if(!onlyvar){
  ############################## spatial or space time ############################################
   if(!bivariate)           {
    if(length(param)==1) {
         optimizer="optimize"
     CompLikelihood <- .GeoModels_optimize_1d(
                              f=eval(as.name(lname)), param=param, lower=lower, upper=upper,
                              namesparam=namesparam,
                              colidx=colidx,rowidx=rowidx,corrmodel=corrmodel, coords=coords,
                              data1=data1,data2=data2, fixed=fixed, fan=fname, n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesaniso=namesaniso,
                              weigthed=weigthed,X=X, MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)}
   if(length(param)>1) {

    if(optimizer=='L-BFGS-B')
      CompLikelihood <- optim(par=param,fn=eval(as.name(lname)),
                              control=list(factr=opt_factr, pgtol=opt_pgtol, maxit=opt_maxit),
                              colidx=colidx,rowidx=rowidx,corrmodel=corrmodel, coords=coords,data1=data1,data2=data2, fixed=fixed,
                              fan=fname, lower=lower, method='L-BFGS-B',n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso,
                              upper=upper,weigthed=weigthed,X=X,  hessian=FALSE,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)

    if(optimizer=='BFGS')
        CompLikelihood <- optim(par=param, fn=eval(as.name(lname)),  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
                           control=list(reltol=opt_reltol, maxit=opt_maxit),data1=data1,data2=data2, fixed=fixed, fan=fname,
                              hessian=FALSE, method='BFGS',n=n,namescorr=namescorr,
                                  namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso,weigthed=weigthed,X=X,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
    if(optimizer=='SANN')
        CompLikelihood <- optim(par=param, fn=eval(as.name(lname)),  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
                           control=list(reltol=opt_reltol, maxit=opt_maxit),data1=data1,data2=data2, fixed=fixed, fan=fname,
                              hessian=FALSE, method='SANN',n=n,namescorr=namescorr,
                                  namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso,weigthed=weigthed,X=X, MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)

      if(optimizer=='Nelder-Mead')
        CompLikelihood <- optim(par=param, fn=eval(as.name(lname)),  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
          control=list(reltol=opt_reltol, maxit=opt_maxit), data1=data1,data2=data2, fixed=fixed, fan=fname,
                              hessian=FALSE, method='Nelder-Mead',n=n,namescorr=namescorr,
                                  namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso,weigthed=weigthed,X=X, MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)


    if(optimizer=='nlm')
    CompLikelihood <- nlm(f=eval(as.name(lname)),p=param,steptol = 1e-4, colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,data1=data1,data2=data2, fixed=fixed,
                               fan=fname,hessian=FALSE,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso,
                               iterlim=opt_maxit, weigthed=weigthed,X=X,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)


     if(optimizer=='bobyqa')
  {
     CompLikelihood <-minqa::bobyqa(par=param, fn = get(lname),lower=lower,upper=upper,
                        control = list(maxfun=opt_maxit),
                        colidx=colidx,rowidx=rowidx,corrmodel=corrmodel, coords=coords, data1=data1,data2=data2, fixed=fixed,fan=fname,n=n,namescorr=namescorr,
                       namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso, weigthed=weigthed,X=X,MM=MM,
                        aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
    }

    if(optimizer=='nlminb'){

     CompLikelihood <-nlminb(objective=eval(as.name(lname)),start=param,colidx=colidx,rowidx=rowidx,corrmodel=corrmodel, coords=coords, data1=data1,data2=data2, fixed=fixed,
                                control = list(iter.max=opt_maxit, eval.max=opt_maxit, rel.tol=opt_reltol, x.tol=1e-8),
                              lower=lower,upper=upper,hessian=FALSE,
                               fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso,
                               weigthed=weigthed,X=X, MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
    }
  }
}
######################################################################################
############################## bivariate  ############################################
######################################################################################
    if(bivariate)           {
     if(length(param)==1)
        {
         optimizer="optimize"
       CompLikelihood <- .GeoModels_optimize_1d(
                              f=comploglik_biv2, param=param, lower=lower, upper=upper,
                              namesparam=namesparam,
                              colidx=colidx,rowidx=rowidx,corrmodel=corrmodel, coords=coords,
                              data1=data1,data2=data2, fixed=fixed,fan=fname,n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesaniso=namesaniso,
                              weigthed=weigthed,X=X,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)}
      if(length(param)>1) {
    if(optimizer=='L-BFGS-B'){
      CompLikelihood <- optim(param,comploglik_biv2, control=list(factr=opt_factr, pgtol=opt_pgtol, maxit=opt_maxit),
                              method='L-BFGS-B',hessian=FALSE,lower=lower, upper=upper,
                               colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
                              data1=data1,data2=data2, fixed=fixed,fan=fname,n=n,
                              namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,namesaniso=namesaniso,
                             weigthed=weigthed,X=X,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)}

      if(optimizer=='BFGS')
      CompLikelihood <- optim(param,comploglik_biv2,  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,control=list(reltol=opt_reltol, maxit=opt_maxit), data1=data1,data2=data2, fixed=fixed, fan=fname,
                              hessian=FALSE, method='BFGS',n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam ,namesaniso=namesaniso,weigthed=weigthed,X=X,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
   if(optimizer=='Nelder-Mead')
      CompLikelihood <- optim(param,comploglik_biv2,  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,control=list(reltol=opt_reltol, maxit=opt_maxit), data1=data1,data2=data2, fixed=fixed, fan=fname,
                              hessian=FALSE, method='Nelder-Mead',n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam ,namesaniso=namesaniso,weigthed=weigthed,X=X,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
    if(optimizer=='nlm')
        CompLikelihood <- nlm( f=comploglik_biv2,p=param,  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel, coords=coords, data1=data1,data2=data2, fixed=fixed,
                               fan=fname,hessian=FALSE,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso,
                               weigthed=weigthed,X=X,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)

   if(optimizer=='bobyqa')   {
        CompLikelihood <- minqa::bobyqa(fn=comploglik_biv2,par=param,
                                     control = list(maxfun=opt_maxit),
                              lower=lower,upper=upper,
                                colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,data1=data1,data2=data2, fixed=fixed,
                               fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso,
                               weigthed=weigthed,X=X,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
                      }
    if(optimizer=='nlminb')
        CompLikelihood <- nlminb( objective=comploglik_biv2,start=param,
                                     control = list(iter.max=opt_maxit, eval.max=opt_maxit, rel.tol=opt_reltol, x.tol=1e-8),
                              lower=lower,upper=upper, hessian=FALSE,
                                colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,data1=data1,data2=data2, fixed=fixed,
                               fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso,
                               weigthed=weigthed,X=X,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)


   }
 }
########################################################################################
########################################################################################
    # Check the optimisation outcome.
    # Important: all optimizers minimize -loglik.  Therefore checks for penalized
    # values must be done on the minimized objective BEFORE changing sign.

    CompLikelihood$optimizer <- optimizer
    CompLikelihood$opt_code <- NA_integer_
    CompLikelihood$opt_message <- NULL
    CompLikelihood$opt_warning <- NULL
    CompLikelihood$warning <- NULL

    penalized_minimum <- function(obj) {
      is.null(obj) || length(obj) == 0L || any(!is.finite(obj)) || any(obj >= 1.0e8)
    }

    bad_parameter <- function(par) {
      is.null(par) || length(par) == 0L || any(!is.finite(par)) || any(abs(par) >= 1.0e14)
    }

    usable_fit <- function(min_obj, par) {
      !penalized_minimum(min_obj) && !bad_parameter(par)
    }

    make_opt_warning <- function(opt_name, code, msg = NULL) {
      out <- paste0(opt_name, ' warning: code=', code)
      if(!is.null(msg) && length(msg) > 0L && !is.na(msg) && nzchar(as.character(msg))) {
        out <- paste0(out, '; ', as.character(msg))
      }
      out
    }

    set_success_or_failure <- function(obj, min_obj, par, opt_name, code, msg = NULL,
                                       success_codes = 0L, iter_codes = integer(0L)) {
      obj$opt_code <- code
      obj$opt_message <- msg

      if(usable_fit(min_obj, par)) {
        obj$convergence <- 'Successful'
        if(!(length(code) == 1L && code %in% success_codes)) {
          if(length(code) == 1L && code %in% iter_codes) {
            obj$opt_warning <- paste0(opt_name, ': iteration limit reached')
          } else {
            obj$opt_warning <- make_opt_warning(opt_name, code, msg)
          }
        }
      } else {
        obj$convergence <- 'Optimization may have failed: Try with other starting parameters'
      }
      obj
    }

    if(optimizer=='Nelder-Mead'||optimizer=='multiNelder-Mead'||optimizer=='SANN'){
        opt_code <- CompLikelihood$convergence
        opt_msg  <- if(!is.null(CompLikelihood$message)) CompLikelihood$message else NULL
        min_obj  <- CompLikelihood$value

        CompLikelihood$value = -min_obj * obj_scale
        names(CompLikelihood$par)<- namesparam

        CompLikelihood <- set_success_or_failure(CompLikelihood, min_obj, CompLikelihood$par,
                                                 optimizer, opt_code, opt_msg,
                                                 success_codes = 0L, iter_codes = 1L)

        CompLikelihood$counts=as.numeric(CompLikelihood$counts[1])
    }

    if(optimizer=='nmk'||optimizer=='nmkb'){
        opt_code <- CompLikelihood$convergence
        min_obj  <- CompLikelihood$value

        CompLikelihood$value = -min_obj * obj_scale
        names(CompLikelihood$par)<- namesparam

        CompLikelihood <- set_success_or_failure(CompLikelihood, min_obj, CompLikelihood$par,
                                                 optimizer, opt_code, NULL,
                                                 success_codes = 0L)

        CompLikelihood$counts=as.numeric(CompLikelihood$feval)
    }

    if(optimizer=='L-BFGS-B'||optimizer=='BFGS'||optimizer=='lbfgsb3c'){
        opt_code <- CompLikelihood$convergence
        opt_msg  <- if(!is.null(CompLikelihood$message)) CompLikelihood$message else NULL
        min_obj  <- CompLikelihood$value

        CompLikelihood$value = -min_obj * obj_scale
        names(CompLikelihood$par)<- namesparam

        CompLikelihood <- set_success_or_failure(CompLikelihood, min_obj, CompLikelihood$par,
                                                 optimizer, opt_code, opt_msg,
                                                 success_codes = 0L, iter_codes = 1L)

        CompLikelihood$counts=as.numeric(CompLikelihood$counts[1])
    }

    if(optimizer=='nlm'){
        opt_code <- CompLikelihood$code
        min_obj  <- CompLikelihood$minimum

        CompLikelihood$par <- CompLikelihood$estimate
        names(CompLikelihood$par)<- namesparam
        CompLikelihood$value <- -min_obj * obj_scale

        CompLikelihood <- set_success_or_failure(CompLikelihood, min_obj, CompLikelihood$par,
                                                 optimizer, opt_code, NULL,
                                                 success_codes = c(1L, 2L),
                                                 iter_codes = 4L)

        CompLikelihood$counts=as.numeric(CompLikelihood$iterations)
    }

    if(optimizer=='bobyqa'){
        opt_code <- CompLikelihood$ierr
        opt_msg  <- if(!is.null(CompLikelihood$msg)) CompLikelihood$msg else NULL
        min_obj  <- CompLikelihood$fval

        names(CompLikelihood$par)<- namesparam
        CompLikelihood$value <- -min_obj * obj_scale

        CompLikelihood <- set_success_or_failure(CompLikelihood, min_obj, CompLikelihood$par,
                                                 optimizer, opt_code, opt_msg,
                                                 success_codes = 0L)

        CompLikelihood$counts=as.numeric(CompLikelihood$feval)
    }

    if(optimizer=='nlminb'||optimizer=='multinlminb'){
        opt_code <- CompLikelihood$convergence
        opt_msg  <- CompLikelihood$message
        min_obj  <- CompLikelihood$objective

        names(CompLikelihood$par)<- namesparam
        CompLikelihood$value <- -min_obj * obj_scale

        CompLikelihood <- set_success_or_failure(CompLikelihood, min_obj, CompLikelihood$par,
                                                 optimizer, opt_code, opt_msg,
                                                 success_codes = 0L)

        CompLikelihood$counts=as.numeric(CompLikelihood$iterations)
    }

    if(optimizer=='optimize'){
        param <- CompLikelihood$minimum
        min_obj <- CompLikelihood$objective

        CompLikelihood$par <- param
        names(CompLikelihood$par)<- namesparam
        CompLikelihood$value <- -min_obj * obj_scale
        CompLikelihood$counts <- NULL

        CompLikelihood <- set_success_or_failure(CompLikelihood, min_obj, CompLikelihood$par,
                                                 optimizer, 0L, NULL,
                                                 success_codes = 0L)
    }
  } ##### end if !onlyvar
    else {
          CompLikelihood=as.list(0)
          names(CompLikelihood)="value"
          CompLikelihood$par <- param
          CompLikelihood$claic <- NULL;CompLikelihood$clbic <- NULL;
          CompLikelihood$convergence <- 'Successful'
          CompLikelihood$opt_code <- 0L
          CompLikelihood$opt_message <- NULL
          CompLikelihood$opt_warning <- NULL
          CompLikelihood$warning <- NULL
          if(!bivariate) CompLikelihood$value = - obj_scale * comploglik2(param=CompLikelihood$par ,  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
                              data1=data1,data2=data2, fixed=fixed, fan=fname,type_cop=type_cop,
                             n=n,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso,weigthed=weigthed,X=X,MM=MM,aniso=aniso,cond_pair=cond_pair)
          else CompLikelihood$value = - obj_scale * comploglik_biv2(param=CompLikelihood$par ,  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
                data1=data1,data2=data2, fixed=fixed, fan=fname,type_cop=type_cop,
                             n=n,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso,weigthed=weigthed,X=X,MM=MM,aniso=aniso,cond_pair=cond_pair)
          ## The sensitivity matrix, when requested, is computed once below
          ## using the numerically scaled parameterization.
  }

#####################################
if (sensitivity || varest)
{
  ## The objective used by the optimizer is the average negative
  ## composite log-likelihood, -CL / obj_scale.  The sensitivity matrix
  ## must be returned on the scale of the total negative composite
  ## log-likelihood.
  ##
  ## Numerical differentiation on the original parameter scale can be
  ## unstable when parameters have very different magnitudes.  We first
  ## obtain a preliminary Hessian, use its diagonal to define a local
  ## rescaling, recompute the Hessian in the scaled parameterization, and
  ## finally transform it back to the original parameter scale.  In exact
  ## arithmetic this gives the same Hessian; the change is only numerical.
  ##
  ## The scaled Hessian is the primary operational estimate.  When its
  ## unit-free diagnostics indicate weak curvature or poor conditioning, a
  ## second Richardson evaluation with a different step is computed in the
  ## same scaled parameterization.  This second matrix is diagnostic only:
  ## it never replaces or regularizes the primary sensitivity matrix.

  theta_hat <- as.numeric(CompLikelihood$par)
  names(theta_hat) <- namesparam
  ntheta <- length(theta_hat)

  if (!bivariate) {
    objective_theta <- function(theta) {
      comploglik2(
        param = theta,
        colidx = colidx,
        rowidx = rowidx,
        corrmodel = corrmodel,
        coords = coords,
        data1 = data1,
        data2 = data2,
        fixed = fixed,
        fan = fname,
        n = n,
        namescorr = namescorr,
        namesnuis = namesnuis,
        namesparam = namesparam,
        namesaniso = namesaniso,
        weigthed = weigthed,
        X = X,
        MM = MM,
        aniso = aniso,
        type_cop = type_cop,
        cond_pair = cond_pair
      )
    }
  } else {
    objective_theta <- function(theta) {
      comploglik_biv2(
        param = theta,
        colidx = colidx,
        rowidx = rowidx,
        corrmodel = corrmodel,
        coords = coords,
        data1 = data1,
        data2 = data2,
        fixed = fixed,
        fan = fname,
        n = n,
        namescorr = namescorr,
        namesnuis = namesnuis,
        namesparam = namesparam,
        namesaniso = namesaniso,
        weigthed = weigthed,
        X = X,
        MM = MM,
        aniso = aniso,
        type_cop = type_cop,
        cond_pair = cond_pair
      )
    }
  }

  symmetrize_hessian <- function(HH) {
    (HH + t(HH)) / 2
  }

  relative_matrix_difference <- function(A, B) {
    denom <- max(1, max(abs(A)), max(abs(B)))
    max(abs(A - B)) / denom
  }

  ## Preliminary Hessian of the average negative composite
  ## log-likelihood.  It is used only to select numerically appropriate
  ## perturbation scales.  Multiplying by obj_scale here would make every
  ## perturbation unnecessarily smaller by sqrt(obj_scale).
  H0_average <- numDeriv::hessian(
    func = objective_theta,
    x = theta_hat,
    method = "Richardson"
  )
  H0_average <- symmetrize_hessian(H0_average)

  if (any(!is.finite(H0_average))) {
    stop(
      "The preliminary numerical sensitivity matrix contains non-finite values.",
      call. = FALSE
    )
  }

  ## Local diagonal rescaling.  Use the estimated curvature whenever it is
  ## finite and positive.  A parameter-based fallback is used only for
  ## invalid or zero preliminary diagonal entries; no global curvature floor
  ## is imposed, because legitimate curvatures can differ by many orders of
  ## magnitude.
  hdiag <- abs(diag(H0_average))
  valid_diag <- is.finite(hdiag) & hdiag > 0

  fallback_scale <- pmax(abs(theta_hat), 1)
  scale_vec <- numeric(ntheta)
  scale_vec[valid_diag] <- 1 / sqrt(hdiag[valid_diag])
  scale_vec[!valid_diag] <- fallback_scale[!valid_diag]

  if (any(!is.finite(scale_vec)) || any(scale_vec <= 0)) {
    stop(
      "Unable to construct a finite positive parameter scaling for the sensitivity matrix.",
      call. = FALSE
    )
  }

  A <- diag(scale_vec, nrow = ntheta, ncol = ntheta)

  objective_scaled <- function(eta) {
    theta <- theta_hat + as.vector(A %*% eta)
    objective_theta(theta)
  }

  compute_scaled_hessian_average <- function(method_args = NULL) {
    hessian_args <- list(
      func = objective_scaled,
      x = rep(0, ntheta),
      method = "Richardson"
    )
    if (!is.null(method_args)) hessian_args$method.args <- method_args

    HH <- do.call(numDeriv::hessian, hessian_args)
    HH <- symmetrize_hessian(HH)

    if (any(!is.finite(HH))) {
      stop(
        "The numerical sensitivity matrix contains non-finite values.",
        call. = FALSE
      )
    }
    HH
  }

  ## Primary Hessian with respect to eta, where
  ## theta = theta_hat + A eta.  Restore the total composite-likelihood
  ## scale only after numerical differentiation.
  H_scaled_average <- compute_scaled_hessian_average()
  H_scaled <- obj_scale * H_scaled_average
  H_scaled <- symmetrize_hessian(H_scaled)

  ## Since H_eta = A H_theta A, transform back as
  ## H_theta = A^{-1} H_eta A^{-1}.
  Ainv <- diag(1 / scale_vec, nrow = ntheta, ncol = ntheta)
  CompLikelihood$hessian <- Ainv %*% H_scaled %*% Ainv
  CompLikelihood$hessian <- symmetrize_hessian(CompLikelihood$hessian)

  rownames(CompLikelihood$hessian) <- namesparam
  colnames(CompLikelihood$hessian) <- namesparam

  names(scale_vec) <- namesparam
  rownames(H_scaled_average) <- namesparam
  colnames(H_scaled_average) <- namesparam
  rownames(H_scaled) <- namesparam
  colnames(H_scaled) <- namesparam

  eig_primary <- eigen(
    H_scaled,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  cond_primary <- tryCatch(
    kappa(H_scaled),
    error = function(e) Inf
  )

  max_abs_eig <- max(abs(eig_primary))
  eigen_ratio_primary <- if (is.finite(max_abs_eig) && max_abs_eig > 0) {
    min(eig_primary) / max_abs_eig
  } else {
    NA_real_
  }

  ## Adaptive diagnostic check.  Thresholds determine only whether a second
  ## numerical evaluation is useful; they never alter the returned Hessian.
  hessian_condition_threshold <- 1000
  hessian_eigenratio_threshold <- 1e-3
  hessian_relative_tolerance <- 1e-2
  hessian_check_d <- 0.05

  needs_stability_check <-
    !is.finite(cond_primary) ||
    cond_primary > hessian_condition_threshold ||
    !is.finite(eigen_ratio_primary) ||
    eigen_ratio_primary <= 0 ||
    abs(eigen_ratio_primary) < hessian_eigenratio_threshold

  hessian_check_performed <- FALSE
  hessian_stable <- TRUE
  hessian_check_error <- NULL
  H_check_scaled_average <- NULL
  H_check_scaled <- NULL
  eig_check <- rep(NA_real_, ntheta)
  cond_check <- NA_real_
  relative_difference_check <- NA_real_
  smallest_eigenvalue_sign_change <- NA

  if (isTRUE(needs_stability_check)) {
    hessian_check_performed <- TRUE

    H_check_result <- tryCatch(
      compute_scaled_hessian_average(
        method_args = list(d = hessian_check_d)
      ),
      error = function(e) e
    )

    if (inherits(H_check_result, "error")) {
      hessian_stable <- FALSE
      hessian_check_error <- conditionMessage(H_check_result)
      warning(
        paste0(
          "The adaptive sensitivity-matrix check failed: ",
          hessian_check_error,
          ". The primary scaled Richardson Hessian has been retained."
        ),
        call. = FALSE
      )
    } else {
      H_check_scaled_average <- H_check_result
      H_check_scaled <- obj_scale * H_check_scaled_average
      H_check_scaled <- symmetrize_hessian(H_check_scaled)

      eig_check <- eigen(
        H_check_scaled,
        symmetric = TRUE,
        only.values = TRUE
      )$values
      cond_check <- tryCatch(
        kappa(H_check_scaled),
        error = function(e) Inf
      )

      relative_difference_check <- relative_matrix_difference(
        H_scaled_average,
        H_check_scaled_average
      )
      smallest_eigenvalue_sign_change <-
        sign(min(eig_primary)) != sign(min(eig_check))

      hessian_stable <-
        is.finite(relative_difference_check) &&
        relative_difference_check <= hessian_relative_tolerance &&
        !isTRUE(smallest_eigenvalue_sign_change)

      rownames(H_check_scaled_average) <- namesparam
      colnames(H_check_scaled_average) <- namesparam
      rownames(H_check_scaled) <- namesparam
      colnames(H_check_scaled) <- namesparam

      if (!isTRUE(hessian_stable)) {
        warning(
          paste0(
            "The numerical sensitivity matrix is unstable across ",
            "Richardson step sizes. The primary scaled Hessian has been ",
            "retained; standard errors and composite-likelihood information ",
            "criteria should be interpreted with caution."
          ),
          call. = FALSE
        )
      } else if (min(eig_primary) <= 0 || min(eig_check) <= 0) {
        warning(
          paste0(
            "The sensitivity matrix is reproducibly non-positive-definite ",
            "under the adaptive Richardson check. The primary scaled ",
            "Hessian has been retained."
          ),
          call. = FALSE
        )
      }
    }
  }

  ## Backward-compatible attributes plus adaptive stability diagnostics.
  attr(CompLikelihood$hessian, "parameter_scaling") <- scale_vec
  attr(CompLikelihood$hessian, "scaled_hessian") <- H_scaled
  attr(CompLikelihood$hessian, "scaled_hessian_average") <- H_scaled_average
  attr(CompLikelihood$hessian, "scaled_condition") <- cond_primary
  attr(CompLikelihood$hessian, "scaled_eigenvalues") <- eig_primary
  attr(CompLikelihood$hessian, "hessian_method") <-
    "scaled_parameterization_Richardson_with_adaptive_check"
  attr(CompLikelihood$hessian, "objective_scale") <- obj_scale
  attr(CompLikelihood$hessian, "hessian_check_performed") <-
    hessian_check_performed
  attr(CompLikelihood$hessian, "hessian_stable") <- hessian_stable
  attr(CompLikelihood$hessian, "hessian_positive_definite") <-
    all(is.finite(eig_primary)) && min(eig_primary) > 0
  attr(CompLikelihood$hessian, "hessian_stability") <- list(
    condition_threshold = hessian_condition_threshold,
    eigenratio_threshold = hessian_eigenratio_threshold,
    relative_tolerance = hessian_relative_tolerance,
    check_d = hessian_check_d,
    primary_condition = cond_primary,
    primary_eigenvalues = eig_primary,
    primary_min_eigen_ratio = eigen_ratio_primary,
    check_performed = hessian_check_performed,
    check_error = hessian_check_error,
    check_condition = cond_check,
    check_eigenvalues = eig_check,
    relative_scaled_difference = relative_difference_check,
    smallest_eigenvalue_sign_change = smallest_eigenvalue_sign_change,
    stable = hessian_stable,
    check_scaled_hessian_average = H_check_scaled_average,
    check_scaled_hessian = H_check_scaled
  )
}
CompLikelihood$score=NULL
if(score)
{
  if(!bivariate)
    CompLikelihood$score=numDeriv::grad(func=comploglik2,
                                        x=as.numeric(CompLikelihood$par),
                                        method="Richardson",
                                        colidx=colidx, rowidx=rowidx,
                                        corrmodel=corrmodel, coords=coords,
                                        data1=data1, data2=data2, fixed=fixed,
                                        fan=fname, n=n, type_cop=type_cop,
                                        namescorr=namescorr, namesnuis=namesnuis,
                                        namesparam=namesparam, namesaniso=namesaniso,
                                        weigthed=weigthed, X=X, MM=MM, aniso=aniso,
                                        cond_pair=cond_pair)

  if(bivariate)
    CompLikelihood$score=numDeriv::grad(func=comploglik_biv2,
                                        x=as.numeric(CompLikelihood$par),
                                        method="Richardson",
                                        colidx=colidx, rowidx=rowidx,
                                        corrmodel=corrmodel, coords=coords,
                                        data1=data1, data2=data2, fixed=fixed,
                                        fan=fname, n=n, type_cop=type_cop,
                                        namescorr=namescorr, namesnuis=namesnuis,
                                        namesparam=namesparam, namesaniso=namesaniso,
                                        weigthed=weigthed, X=X, MM=MM, aniso=aniso,
                                        cond_pair=cond_pair)

  ## numDeriv differentiates the normalized objective -CL / obj_scale.
  ## Restore the gradient to the scale of the total negative composite
  ## log-likelihood. GeoVarest will only need to reverse its sign.
  CompLikelihood$score <- obj_scale * CompLikelihood$score
  names(CompLikelihood$score) <- namesparam
}
####################################
       serious_failure <- is.null(CompLikelihood$value) ||
                          length(CompLikelihood$value) == 0L ||
                          !is.finite(CompLikelihood$value) ||
                          CompLikelihood$value <= -1.0e14 ||
                          grepl('^Optimization may have failed', CompLikelihood$convergence)

       if(serious_failure) {
           CompLikelihood$warning <- 'Optimization may have failed: try with other starting values'
           warning(CompLikelihood$warning, call. = FALSE)
       } else {
           ## Soft optimizer warnings are kept in CompLikelihood$opt_warning,
           ## but no R warning is emitted. This avoids false alarms when the
           ## optimum is finite, non-penalized, and usable.
           if(varest){
              stop("subsambling is not working")
           }
       }
      if(hessian) CompLikelihood$sensmat=CompLikelihood$hessian
    return(CompLikelihood)
  }