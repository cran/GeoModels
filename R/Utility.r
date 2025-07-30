####################################################
### File name: Utility.r
####################################################

# Check if the correlation is bivariate
CheckBiv <- function(numbermodel)
{
    CheckBiv <- NULL
    if(!is.null(numbermodel)){
    if(numbermodel >= 101 & numbermodel <= 140) CheckBiv <- TRUE
    else CheckBiv <- FALSE
    }
    return(CheckBiv)
}
# Check if the correlation is spatial or spatial-temporal
CheckST <- function(numbermodel)
{
    CheckST <- NULL
     if(!is.null(numbermodel)){
    if(numbermodel >40 & numbermodel <= 100) CheckST <- TRUE
    else  CheckST <- FALSE
     }
    return(CheckST)
}


# Check the type of distances
CheckDistance <- function(distance) {
  out <- switch(distance,
                eucl = 0, Eucl = 0, 
                chor = 1, Chor = 1, 
                geod = 2, Geod = 2,
                stop("Invalid type of  di distance: use 'Eucl', 'Chor' o 'Geod'."))
  return(out)
}
### Procedures are in alphabetical order.
CkCorrModel <- function(corrmodel)
  {
    CkCorrModel <- NULL
    # Correlation function are in alphabetical order
    CkCorrModel <- switch(corrmodel,
            #spatial models
                             cauchy=1,Cauchy=1,
                             Matern1=2,
                             Matern2=3,
                             exponential=4,Exponential=4,Exp=4,exp=4,Matern0=4,
                             dagum = 5, Dagum = 5,
                             GenWend_Matern=6,Genwend_Matern=6,
                             GenWend_Matern2=7,Genwend_Matern2=7,
                             Gencauchy=8,GenCauchy=8,
                             Shkarofski=10,shkarofski=10,
                             Wend0=11,wend0=11,
                             stable=12,Stable=12,
                             Wend1=13,wend1=13,
                             matern=14,Matern=14,
                             Wend2=15,wend2=15,
                             wave=16,Wave=16,
                             Multiquadric=17,multiquadric=17,
                             Sinpower=18,sinpower=18,
                             Genwend=19,GenWend=19,
                             F_sphere=20,F_Sphere=20,
                             Hypergeometric2=21,HyperGeometric2=21, hypergeometric2=21,
                             Hypergeometric=22,HyperGeometric=22, hypergeometric=22,
                             Hypergeometric_Matern=23,HyperGeometric_Matern=23, hypergeometric_Matern=23,
                             Hypergeometric_Matern2=30,Hypergeometric_matern2=30, hypergeometric_Matern2=30,
                             Kummer=24,Kummer=24,
                             Kummer_Matern=25,Kummer_matern=25,
                             GenWend_Hole=26,GenWend_hole=26,
                             Matern_Hole=27,Matern_hole=27,
                             Schoenberg=28,schoenberg=28,
                             GenWend_Matern_Hole=29,GenWend_Matern_hole=29,
             # spatial-temporal non-separable models
                             gneiting=42,Gneiting=42,  #ok
                             iacocesare=44,Iacocesare=44, #ok
                             porcu=46,Porcu=46,
                             stein=48,Stein=48,          #ok
                             porcu1=50,Porcu1=50,
                             gneiting_GC=52,Gneiting_GC=52, #ok
                             gneiting_GC2=54,Gneiting_GC2=54,
                             sinpower_st=56,Sinpower_st=56,    #ok
                             multiquadric_st=58,Multiquadric_st=58,   #ok
                             gneiting_mat_T=61,Gneiting_mat_T=61, #ok
                             gneiting_mat_S=62,Gneiting_mat_S=62,Gneiting_Mat_S=62, #ok
                             Wen0_space=63,wen0_space=63,  #ok
                             Wen0_time=64,wen0_time=64,    #ok
                             Wen1_space=65,wen1_space=65,  #ok
                             Wen1_time=66,wen1_time=66,    #ok
                             Wen2_space=67,wen2_space=67,  #ok
                             Wen2_time=68,wen2_time=68,    #ok
                             Matern_Matern_sep=88,
                             Gneiting_wen_S=87,gneiting_wen_S=87,Gneiting_Wen_S=87,
                             Gneiting_wen_T=88,gneiting_wen_T=88, #ok
                             Matern_Matern_nosep=89,
                             Genwend_Genwend_nosep=85,
              # spatial-temporal separable models
                             Wend0_Wend0=69,wend0_wend0=69, #ok
                             Wend0_Wend1=70,wend0_wend1=70, #ok
                             Wend0_Wend2=71,wend0_wend2=71, #ok
                             Wend1_Wend0=72,wend1_wend0=72, #ok
                             Wend1_Wend1=73,wend1_wend1=73, #ok
                             Wend1_Wend2=74,wend1_wend2=74, #ok
                             Wend2_Wend0=75,wend2_wend0=75, #ok
                             Wend2_Wend1=76,wend2_wend1=76, #ok
                             Wend2_Wend2=77,wend2_wend2=77, #ok 
                             GenWend_GenWend=78,Genwend_Genwend=78, #ok 
                             exp_cauchy=82,Exp_Cauchy=82,
                             exp_exp=84,Exp_Exp=84,
                             Matern_Matern=86, matern_matern=86,  #ok
                             stable_stable=94,Stable_Stable=94,
                             prove=96,
              # Bivariate models
                             Bi_wend0_sep=111,Bi_Wend0_sep=111,
                             Bi_Wend0_contr=129,Bi_wend0_contr=129,
                             Bi_wend0=112,Bi_Wend0=112,
                             Bi_F_sphere=117,Bi_F_sphere=117,
                              Bi_F_sphere_sep=119,Bi_F_sphere_sep=119,
                              Bi_F_sphere_contr=121,Bi_F_sphere_contr=121,
                             Bi_wend1_sep=113,Bi_Wend1_sep=113,
                             Bi_Wend1_contr=131,Bi_wend1_contr=131,
                             Bi_wend1=114,Bi_Wend1=114,  
                             Bi_wend2_sep=115,Bi_Wend2_sep=115,
                             Bi_wend2=116,Bi_Wend2=116,
                             Bi_Wend2_contr=120,Bi_wend2_contr=120,
                             Bi_matern_contr=118,Bi_Matern_contr=118, 
                             Bi_matern_sep=122,  Bi_Matern_sep=122, 
                             Bi_matern=128,Bi_Matern=128,
                             Bi_LMC_contr=124,
                             Bi_LMC=126,
                             Bi_GenWend_sep=130,Bi_genWend_sep=130,
                             Bi_GenWend_contr=134,Bi_genWend_contr=134,
                             Bi_GenWend=132,Bi_genWend=132,
                             Bi_Matern_Cauchy=136,  
                             Bi_GenMatern_Cauchy=137,   Bi_genMatern_Cauchy=137,
                             Bi_genCauchy=138,Bi_GenCauchy=138,
                             Bi_Stable=139,Bi_stable=139,   
            ########Tapers
                   ##spatial
                             Bohman=28,
                             Wendland0=30,wendland0=30,
                             Wendland1=32,wendland1=32,
                             Wendland2=34, wendland2=34,
                             unit_matrix=36,
                             Spherical=38, spherical=38,
                   ##bivariate
                             Bi_Wendland1=140,Bi_wendland1=140,
                             Bi_Wendland2=142,Bi_wendland2=142,
                             Bi_Wendland3=144,Bi_wendland3=144,
                             Bi_wendland_asy=146,Bi_wendland_asy=146,
                             unit_matrix_biv=147,
                   ##spacetime separable
                             Wendland0_Wendland0=200,
                             Wendland0_Wendland1=202,
                             Wendland0_Wendland2=204,
                             Wendland1_Wendland0=206,
                             Wendland1_Wendland1=208,
                             Wendland1_Wendland2=210,
                             Wendland2_Wendland0=212,
                             Wendland2_Wendland1=214,
                             Wendland2_Wendland2=216,
                    ##spacetime no separable
                             Wendland0_time=218,
                             Wendland1_time=220,
                             Wendland2_time=222,
                             Wendland0_space=224,
                             Wendland1_space=226,
                             Wendland2_space=228,
                             unit_matrix_st=230)
    return(CkCorrModel)
  }



# checking models valid only on the sphere
CheckSph<- function(numbermodel)
  {
    Check <- FALSE
    if(numbermodel %in% c(17,18,56,58,20,117))  Check=TRUE    
    return(Check)
  }


CkInput <- function(coordx, coordy, coordz,coordt, coordx_dyn, corrmodel, data, distance, fcall, fixed, grid,
                      likelihood, maxdist, maxtime,  model, n,  optimizer, param,
                       radius, start, taper, tapsep, type, varest, weighted,
                       copula,X)
  {
    error <- NULL
    replicates=1


    # START Include internal functions:
    CheckParamRange <- function(param)
    { 
      #  if(!is.na(param['df'])) if(param['df'] > 1/2 || param['df'] < 0 ) return(FALSE)
        #if(!is.na(param['tail'])) if(param['tail'] >0.5) return(FALSE)
        if(!is.na(param['shape'])) if(param['shape'] <0) return(FALSE)
        if(!is.na(param['nugget'])) if(param['nugget'] < 0||param['nugget'] >= 1) return(FALSE)
        if(!is.na(param['nugget_1'])) if(param['nugget_1'] < 0) return(FALSE)
        if(!is.na(param['nugget_2'])) if(param['nugget_2'] < 0) return(FALSE)
        if(!is.na(param['power'])) if(param['power'] <=0 || param['power'] > 2) return(FALSE)
        if(!is.na(param['power_s'])) if(param['power_s'] <=0 || param['power_s'] > 2) return(FALSE)
        if(!is.na(param['power_t'])) if(param['power_t'] <=0 || param['power_t'] > 2) return(FALSE)
       # if(!is.na(param['power1'])) if(param['power1'] <=0 || param['power1'] > 2) return(FALSE)
        if(!is.na(param['power2'])) if(param['power2'] <= 0) return(FALSE)
        if(!is.na(param['power2_1'])) if(param['power2_1'] <= 0) return(FALSE)
        if(!is.na(param['power2_12'])) if(param['power2_12'] <= 0) return(FALSE)
        if(!is.na(param['power2_2'])) if(param['power2_2'] <= 0) return(FALSE)
        if(!is.na(param['sep'])) if(param['sep'] < 0 || param['sep'] > 1) return(FALSE)
        if(!is.na(param['scale'])) if(param['scale'] <= 0) return(FALSE)
        if(!is.na(param['scale_s'])) if(param['scale_s'] <= 0) return(FALSE)
        if(!is.na(param['scale_t'])) if(param['scale_t'] <= 0) return(FALSE)
        if(!is.na(param['scale_1']))   if(param['scale_1'] < 0) return(FALSE)
        if(!is.na(param['scale_2']))   if(param['scale_2'] < 0) return(FALSE)
        if(!is.na(param['scale_12']))   if(param['scale_12'] <= 0) return(FALSE)
        if(!is.na(param['sill'])) if(param['sill'] <= 0) return(FALSE)
        if(!is.na(param['sill_1'])) if(param['sill_1'] <= 0) return(FALSE)
        if(!is.na(param['sill_2'])) if(param['sill_2'] <= 0) return(FALSE)
        #if(!is.na(param['smooth'])) if(param['smooth'] <= 0) return(FALSE)
        if(!is.na(param['smooth_s'])) if(param['smooth_s'] < 0) return(FALSE)
        if(!is.na(param['smooth_t'])) if(param['smooth_t'] < 0) return(FALSE)
        #if(!is.na(param['smooth_1'])) if(param['smooth_1'] < 0) return(FALSE)
        if(!is.na(param['smooth_2'])) if(param['smooth_2'] < 0) return(FALSE)
        #if(!is.na(param['smooth_12'])) if(param['smooth_12'] < 0) return(FALSE)
        if(!is.na(param['pcol'])) if(param['pcol'] > 1 || param['pcol'] < -1  ) return(FALSE)
        return(TRUE)  
    }
   
    bivariate=CheckBiv(CkCorrModel(corrmodel))
    if(!bivariate)  if(length(coordt)>0&&is.list(X)) X=X[[1]]
    if(!bivariate) {if(is.null(X))  {X=1;num_betas=1} 
                    else num_betas=ncol(X)  }
    if( bivariate) {if(is.null(X))  {X=1;num_betas=c(1,1)} 
                    else { 
                               if(is.list(X))  num_betas=c(ncol(X[[1]]),ncol(X[[2]]))
                               else  num_betas=c(ncol(X),ncol(X))
                          } }
    #bivariate<-CheckBiv(CheckCorrModel(corrmodel))
    #spacetime<-CheckST(CheckCorrModel(corrmodel))
    #if(is.null(bivariate)&&is.null(bivariate))
    # END Include internal functions
    if(!is.null(copula)) {
       {if(copula!="Gaussian"&&copula!="Clayton"&&copula!="SkewGaussian") 
        error <- 'Copula model is wrong\n'
        return(list(error=error))}
       }

    ### check not kriging inputs   
    if(fcall!="Kriging"){

    ### START Checks inserted input
    # START common check fitting and simulation
    if(missing(coordx_dyn)){
        if(missing(coordx) || !is.numeric(coordx)){
        error <- 'insert a suitable set of numeric coordinates\n'
        return(list(error=error))}
        if(!is.null(coordy) & !is.numeric(coordy)){
        error <- 'insert a suitable set of numeric coordinates\n'
        return(list(error=error))}
          if(!is.null(coordz) & !is.numeric(coordz)){
        error <- 'insert a suitable set of numeric coordinates\n'
        return(list(error=error))}
    }
    if(missing(corrmodel) || !is.character(corrmodel)){
        error <- 'insert the correlation model\n'
        return(list(error=error))}
    if(!is.null(grid) & !is.logical(grid)){
        error <- 'the parameter grid need to be a logic value\n'
        return(list(error=error))}
    if(!is.null(model) & !is.character(model)){
        error <- 'insert the name of the random field\n'
        return(list(error=error))}
    if(is.null(CkModel(model))){
        error <- 'the model name of the random field is not correct\n'
        return(list(error=error))}
    if(is.null(replicates) || (abs(replicates-round(replicates))>0) || replicates<1){
        error <- 'the parameter replicates need to be a positive integer\n'
        return(list(error=error))}
    if(is.null(CheckDistance(distance))){
        error <- 'the name of the distance is not correct\n'
        return(list(error=error))}
    if(is.null(CkCorrModel(corrmodel))){
        error <- 'the name of the correlation model is not correct\n'
        return(list(error=error))}
    if(radius<0){
        error <- 'the radius of the sphere must be positive\n'
        return(list(error=error))}
    if(CkModel(model)==11&&(all(n<1)||!all(is.numeric(n))))
    {error <- 'the parameter n for the Binomial RF is wrong \n'
        return(list(error=error))}

    # START check fitting
    if(fcall=="Fitting"){
         if(!is.null(coordx_dyn))
            {
                if(!is.list(coordx_dyn)){
                error<-"Dynamical coordinates must be a list object"
              return(list(error=error))}
              if(!is.list(data)){
                error<-"Data must be a list"
              return(list(error=error))}
          }
        else {if(missing(data) || !is.numeric(data)){
              error <- 'insert a numeric vector or matrix of data\n'
              return(list(error=error))}}
        if(!is.null(fixed) & !is.list(fixed)){
            error <- 'insert fixed values as a list of parameters\n'
            return(list(error=error))}
        if(!is.null(fixed)){ 
            namfixed <- names(fixed)
     
        if(!all(namfixed %in% c(NuisParam2(model,CheckBiv(CkCorrModel(corrmodel)),num_betas,copula),CorrelationPar(CkCorrModel(corrmodel))))){
                error <- 'some names of the fixed parameters is/are not correct\n'
                return(list(error=error))}
        if(!CheckParamRange(unlist(fixed))){
            error <- 'some fixed values are out of the range\n'
            return(list(error=error))}}
           else {namfixed=NULL}   
        if(!is.null(likelihood) & !is.character(likelihood)){
            error <- 'insert the type of likelihood objects\n'
            return(list(error=error))}
        for(i in 1:length(maxdist)){
        if(!is.null(maxdist[i])){
            error <- "insert a positive numeric value for the maximum spatial distance\n"
            if(!is.numeric(maxdist)) return(list(error=error))
            else if(maxdist[i]<0) return(list(error=error))}}
        if(!is.null(maxtime)){
            error <- "insert a positive numeric value for the maximum time interval\n"
            if(!is.numeric(maxtime)) return(list(error=error))
            else if(maxtime<0) return(list(error=error))}
        if(!is.null(optimizer) & !is.character(optimizer)){
            error <- 'insert the type of maximising algorithm\n'
            return(list(error=error))}
        if(!is.null(varest) & !is.logical(varest)){
            error <- 'the parameter std.err need to be a logical value\n'
            return(list(error=error))}


        if(!is.null(type) & !is.character(type)){
            error <- 'insert the configuration of the likelihood objects\n'
            return(list(error=error))}
        if(is.null(CkType(type))){
            error <- 'the type name of the likelihood objects is not correct\n'
            return(list(error=error))}
        if(is.null(CkLikelihood(likelihood))){
           error <- 'the setting name of the likelihood objects is not correct\n'
           return(list(error=error))}
        if(likelihood == "Full"){
            if(!any(type == c("Restricted", "Standard", "Tapering","Tapering1","CV","Tapering2"))){
                error <- 'insert a type name of the likelihood objects compatible with the full likelihood\n'
                return(list(error=error))}}
        if(likelihood == "Marginal"){
            if(!any(type == c("Difference", "Pairwise","Independence"))){
                error <- 'insert a type name of the likelihood objects compatible with the composite-likelihood\n'
                return(list(error=error))}}
        if(varest & (likelihood == "Conditional" || likelihood == "Difference" ||
            likelihood == "Marginal"|| likelihood == "Marginal_2")){
            error <- 'insert the type of estimation method for the variances\n'
            return(list(error=error))}
        if(varest  & (likelihood == "Conditional" || likelihood == "Difference" ||
            likelihood == "Marginal"|| likelihood == "Marginal_2")){
            error <- 'the name of the estimation type for the variances is not correct\n'
            return(list(error=error))}

  ### check on starting parameters         
  if(!is.null(start)){
            if(!is.list(start)){
                error <- 'insert starting values as a list of parameters\n'
                return(list(error=error))}
 namstart <- names(start)

  if(!bivariate){
       if(num_betas>1){
       if(ncol(X)!=sum(substr(c(namstart,namfixed),1,4)=="mean"))
       {error <- 'number of covariates must be equal to to the regressian mean parameters\n'
        return(list(error=error)) }}
    }
  if(bivariate){
       if(num_betas[1]>1&&num_betas[2]>1){
       if(is.list(X)) 
            if(ncol(X[[1]])!=sum(substr(c(namstart,namfixed),1,6)=="mean_1")&&
          ncol(X[[2]])!=sum(substr(c(namstart,namfixed),1,6)=="mean_2"))
       { error <- 'number of covariates must be equal to to the regressian mean parameters\n'
                return(list(error=error)) }}
    }
  if(!all(namstart %in% c(NuisParam2(model,CheckBiv(CkCorrModel(corrmodel)),num_betas,copula), CorrelationPar(CkCorrModel(corrmodel)))))
    { error <- 'some names of the starting parameters is/are not correct\n'
      return(list(error=error))}

  if(any(namstart=='mean') & (type=='Difference' || type=='Restricted'))
  { error <- 'the mean parameter is not allow with the difference composite likelihood\n'
    return(list(error=error))}
  if(!CheckParamRange(unlist(start))){
                error <- 'some starting values are out of the range\n'
                return(list(error=error))
  }

   if(!is.null(fixed))
    if(any(namstart %in% namfixed)){
        error <- 'some fixed parameter name/s is/are matching with starting parameter name/s\n'
        return(list(error=error))}
  }
    ### end check on starting parameters  

  #### START - checks the format of the coordinates and dataset##
  if(!is.null(coordx_dyn))
     {AS=1}
  else
    {
    dimdata <- dim(data) # set the data dimension

    if(is.null(coordt)) # START 1) spatial random field
      {   
        if(CheckST(CkCorrModel(corrmodel)))
          {
            error <- 'temporal coordinates are missing\n'
            return(list(error=error))}
            if(grid) # START regular grid
              {  if(CheckBiv(CkCorrModel(corrmodel))) 
                    {
                      if(is.null(dimdata))
                      {error <- c('insert an array of d x d x 2  spatial observations\n')
                          return(list(error=error))}
                      if(length(dimdata)!=3)
                      {error <- c('the dimension of the data matrix is not correct\n')
                          return(list(error=error))}
                      if(length(coordx)!=dimdata[1] || length(coordy)!=dimdata[2]|| length(coordz)!=dimdata[3])
                      {error <- c('the number of coordinates does not match with the number of spatial observations\n')
                      return(list(error=error))}
                    }
                  else {                  
                     if(is.null(dimdata))
                     {error <- c('insert a matrix d x d of spatial observations\n')
                     return(list(error=error))}
                     if(length(dimdata)!=2)
                     {error <- c('the dimension of the data matrix is not correct\n')
                     return(list(error=error))}
                     if(length(coordx)!=dimdata[1] || length(coordy)!=dimdata[2]|| length(coordz)!=dimdata[3])
                     {error <- c('the number of coordinates does not match with the number of spatial observations\n')
                     return(list(error=error))}
                    }
              } # END regular grid
            else # START irregular grid
              {
                numsite <- length(data)
                if(CheckBiv(CkCorrModel(corrmodel))) numsite<- length(data)/2
                if(is.null(numsite))
                  {error <- c('insert a vector of spatial observations\n')
                    return(list(error=error))}
                if(is.null(coordy))
                  { dimcoord <- dim(coordx)
                    if(is.null(dimcoord))
                      {error <- c('insert a suitable set of coordinates\n')
                      return(list(error=error))}
                    else
                      { if(dimcoord[1]!=numsite || (dimcoord[2]!=2 && dimcoord[2]!=3 ))
                          {
                            error <- c('the number of coordinates does not match with the number of spatial observations\n')
                            return(list(error=error))
                          }}
                   }
                else
                  { if(length(coordx)!=length(coordy))
                      {error <- c('the number of the two coordinates does not match\n')
                        return(list(error=error))}
                    if(length(coordx)!=numsite)
                      {error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))}
                  }
              }
      } # END 1) spatial random field
    else # 2) case: spatial-temporal random field
      {
        if(!is.numeric(coordt))
          { error <- 'insert a numerical vector of temporal coordinates\n'
            return(list(error=error))}
        if(length(coordt)<=1)
          { error <- 'insert a numerical vector of temporal coordinates\n'
            return(list(error=error))}
    if(grid) # START regular grid
              {
                if(is.null(dimdata))
                  { error <- c('insert an array of d x d x t spatial-temporal observations\n')
                    return(list(error=error))}
                if(length(dimdata)!=3)
                  { error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))}
                if(length(coordx)!=dimdata[1] || length(coordy)!=dimdata[2]|| length(coordz)!=dimdata[3])
                  { error <- c('the number of coordinates does not match with the number of spatial observations\n')
                    return(list(error=error))}
                if(dimdata[3]!=length(coordt))
                  { error <- c('the time coordinate does not match with the third dimension of the data array\n')
                    return(list(error=error))}
              } # END regular grid
            else # START irregular grid
              {
                if(is.list(coordx_dyn)) {
                 if(length(coordx_dyn)!= length(coordt))
                      { error <- c('the time coordinate does not match with spatial dynamic coordinates number\n')
                    return(list(error=error))}
                }
                if(is.null(dimdata))
                  { error <- c('insert a matrix of t x d spatial-temporal observations\n')
                    return(list(error=error))}
                if(length(dimdata)!=2)
                  { error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))}
                if(is.null(coordy))
                  {  
                    if(dimdata[2]!=nrow(coordx) || dimdata[2]!=nrow(coordx))
                      { error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))  }
                  }
                else
                  {
                    if(length(coordx)!=length(coordy))
                      {error <- c('the number of the  spatial coordinates does not match\n')
                        return(list(error=error))}
                    if(length(coordx)!=dimdata[2])
                      { error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))}
                  }
                if(dimdata[1]!=length(coordt))
                  { error <- c('the time coordinate does not match with the number of the matrix rows\n')
                    return(list(error=error))}
              } # END irregular grid
      }
      # END - check the format of the inserted coordinates and dataset
      } 
    }
    # END check fitting
    # START check simulation
    if(fcall=="Simulation"){
        if(type=="Tapering")
          {
  
         if(!is.null(coordx_dyn))
            {if(!is.list(coordx_dyn)){
                error<-"Dynamical coordinates mst be a list object"
              return(list(error=error))}}
          if(is.list(coordx_dyn)) {
                 if(length(coordx_dyn)!= length(coordt))
                      { error <- c('the time coordinate does not match with spatial dynamic coordinates number\n')
                    return(list(error=error))}}
            if(!is.null(maxdist)){
            error <- "insert a positive numeric value for the maximum spatial distance\n"
            if(!is.numeric(maxdist)) return(list(error=error))
            else if(maxdist<0) return(list(error=error))}
            if(!is.null(maxtime)){
            error <- "insert a positive numeric value for the maximum temporal distance\n"
            if(!is.numeric(maxtime)) return(list(error=error))
            else if(maxtime<0) return(list(error=error))}
            if(!is.null(tapsep)){
            error <- "separability parameter of spacetime taper must be between 0  and 1\n"
            if(!is.numeric(tapsep)) return(list(error=error))
            else if(tapsep<0||tapsep>1) return(list(error=error))}
          }  ## end tapering

        if(is.null(param) || !is.list(param)){
            error <- 'insert the parameters as a list\n'
            return(list(error=error))}
       biv<-CheckBiv(CkCorrModel(corrmodel))
 
       a1=unique(c(NuisParam2("Gaussian",biv,num_betas,NULL),
                    NuisParam2(model,biv,num_betas,copula)))
       a2=CorrelationPar(CkCorrModel(corrmodel))
  
        if(length(param)!= length(c(a1,a2)))
             {
            error <- "some parameters are missing or does not match with the declared model\n"
            return(list(error=error))}

        if(!all( names(param) %in% c(unique(c(NuisParam2("Gaussian",biv,num_betas,copula),NuisParam2(model,biv,num_betas,copula))),
                                      CorrelationPar(CkCorrModel(corrmodel))))){
            error <- 'some names of the parameters are not correct\n'
            return(list(error=error))}
           
        if(is.list(coordx_dyn)) {
                if(biv) coordt=c(1,2)
                 if(length(coordx_dyn)!= length(coordt)){
                      { error <- c('the number of temporal instants does not match with  the dynamic  spatial coordinates number\n')
                    return(list(error=error))}}    
        }
        if(!CheckParamRange(unlist(param))){
            error <- 'some parameters are out of the range\n'
          return(list(error=error))}
     }# END check simulation
  
    }
 else{   

    if(missing(coordx)) {
    error <- "spatial locations must be a matrix of dimension 2\n"
    return(list(error=error))}
    else     {
    if(is.vector(coordx)&&!length(coordx)==2)    {
           error <- "spatial locations must be a vector of dimension 2\n"
           return(list(error=error))}
    if(is.matrix(coordx)&&!ncol(coordx)==2)       {
           error <- "spatial locations must be  a matrix of dimension 2\n"
           return(list(error=error))}
    }
   if((!is.null(coordt)&&!is.numeric(coordt))  ){
           error <- "time  must be a vector\n"
           return(list(error=error))}
   if(!type %in% c("simple","ordinary","Simple","Ordinary")){
           error <-"kriging type can be  Simple or Ordinary\n"
   return(list(error=error))}
   if(missing(data) || !is.numeric(data)){
           error <- "insert a numeric vector of data\n"
           return(list(error=error))}
  if(is.list(coordx_dyn)) {
                 if(length(coordx_dyn)!= length(coordt))
                      { error <- c('the time coordinate does not match with spatial dynamic coordinates number\n')
                    return(list(error=error))}}
        }

      # END check kriging
  }

CkLikelihood <- function(likelihood)
  {
    CkLikelihood <- switch(likelihood,
                              None=0,
                              Conditional=1,
                              Full=2,
                              Marginal=3,
                              Marginal_2=1,
                              Difference=4)
    return(CkLikelihood)
  }

CkModel <- function(model)
  {
    CkModel <- switch(model,
                         None=0,
                         Gaussian=1,Gauss=1,
                         Binary=2,
                         Tukeygh=9,
                         SkewGaussian=10,SkewGauss=10,
                         Binomial=11,
                         StudentT=12,
                         Wrapped=13,
                         Geom=14,Geometric=14,
                         PoisBin=15,
                         BinomialNeg=16,
                         PoisBinNeg=17,
                         SkewStudentT=18,
                         Binomial2=19,
                         SinhAsinh=20,
                         Gamma=21,
                         LogGaussian=22,LogGauss=22,
                         Gamma2=23,
                         LogLogistic=24,
                         Logistic=25,
                         Weibull=26,
                         TwoPieceStudentT=27,
                         Beta=28,
                         TwoPieceGaussian=29,TwoPieceGauss=29,
                         Poisson=30,poisson=30,
                         Binomial_TwoPieceGaussian=31,
                         Binomial_TwoPieceGauss=31,
                         BinomialNeg_TwoPieceGaussian=32,
                         BinomialNeg_TwoPieceGauss=32,
                         Kumaraswamy=33,
                         Tukeyh=34,tukeyh=34,
                         Gaussian_misp_StudentT=35,
                         Gaussian_misp_Poisson=36,
                         Gaussian_misp_SkewStudentT=37,
                         TwoPieceTukeyh=38,
                         TwoPieceBimodal=39,
                         Tukeyh2=40,tukeyh2=40,
                         Gaussian_misp_Tukeygh=41,
                         Kumaraswamy2=42,
                         PoissonZIP=43,
                         Gaussian_misp_PoissonZIP=44,
                         BinomialNegZINB=45,
                         PoissonGamma=46,poissongamma=46,
                         Gaussian_misp_PoissonGamma=47,
                         BinomialLogistic=49,Binomiallogistic=49,
                         Beta2=50,
                         Gaussian_misp_Binomial=51,
                         Gaussian_misp_BinomialNeg=52,
                         PoissonZIP1=53,
                         Binary_misp_BinomialNeg=54,
                         BinomialNegZINB1=56,
                         PoissonGammaZIP=57,
                         PoissonGammaZIP1=58,
                         )
    return(CkModel)
  }
  
CkType <- function(type)
  {
    CkType <- switch(type,
                        Difference=1,
                        Pairwise=2,
                        Restricted=3,
                        Standard=4,
                        Tapering=5,
                        Tapering2=5,
                        Tapering1=6,
                        GeoWLS=7,
                        CV=8,
                        Independence=9)
    return(CkType)
  }

#####  names of the correlation models ###############
CorrParam <- function(corrmodel)
{
  val <- CorrelationPar(CkCorrModel(corrmodel))
  if (is.null(val)) warning("Unrecognized correlation model ")
  return(val)
}

#####  names of the correlation models ###############
CorrelationPar <- function(corrmodel) {
  param_map <- list(
    "1" = c("power2", "scale"),
    "2" = c("scale"),
    "3" = c("scale"),
    "4" = c("scale"),
    "5" = c("power1", "power2", "scale"),
    "6" = c("power2", "scale", "smooth"),
    "7" = c("power2", "scale", "smooth"),
    "8" = c("power1", "power2", "scale"),
    "10" = c("scale_1", "scale_2", "smooth"),
    "11" = c("power2", "scale"),
    "12" = c("power", "scale"),
    "13" = c("power2", "scale"),
    "14" = c("scale", "smooth"),
    "15" = c("power2", "scale"),
    "16" = c("scale"),
    "17" = c("power", "scale"),
    "18" = c("power"),
    "19" = c("power2", "scale", "smooth"),
    "20" = c("scale", "smooth"),
    "21" = c("power1", "power2", "scale", "smooth"),
    "22" = c("power2", "scale", "smooth"),
    "23" = c("power2", "scale", "smooth"),
    "24" = c("power2", "scale", "smooth"),
    "25" = c("power2", "scale", "smooth"),
    "26" = c("power1", "power2", "scale", "smooth"),
    "27" = c("power1", "scale", "smooth"),
    "28" = c("scale"),
    "29" = c("power1", "power2", "scale", "smooth"),
    "30" = c("power2", "scale", "smooth"),
    "42" = c("power_s", "power_t", "scale_s", "scale_t", "sep"),
    "44" = c("power2", "power_s", "power_t", "scale_s", "scale_t"),
    "46" = c("power_s", "power_t", "scale_s", "scale_t", "sep"),
    "48" = c("power_t", "scale_s", "scale_t", "smooth_s"),
    "50" = c("power_s", "power_t", "scale_s", "scale_t", "sep"),
    "52" = c("power_s", "power_t", "scale_s", "scale_t", "sep"),
    "54" = c("power_s", "power_t", "scale_s", "scale_t", "sep"),
    "56" = c("scale_s", "scale_t", "smooth_t"),
    "58" = c("power_s", "power_t", "scale_s", "scale_t"),
    "60" = c("power_s", "power_t", "scale_s", "scale_t", "sep"),
    "61" = c("power_s", "power2_s", "scale_s", "scale_t", "sep", "smooth_t"),
    "62" = c("power_t", "power2_t", "scale_s", "scale_t", "sep", "smooth_s"),
    "63" = c("power_t", "power2_s", "power2_t", "scale_s", "scale_t", "sep"),
    "64" = c("power_s", "power2_s", "power2_t", "scale_s", "scale_t", "sep"),
    "65" = c("power_t", "power2_s", "power2_t", "scale_s", "scale_t", "sep"),
    "66" = c("power_s", "power2_s", "power2_t", "scale_s", "scale_t", "sep"),
    "67" = c("power_t", "power2_s", "power2_t", "scale_s", "scale_t", "sep"),
    "68" = c("power_s", "power2_s", "power2_t", "scale_s", "scale_t", "sep"),
    "69" = c("power2_s", "power2_t", "scale_s", "scale_t"),
    "70" = c("power2_s", "power2_t", "scale_s", "scale_t"),
    "71" = c("power2_s", "power2_t", "scale_s", "scale_t"),
    "72" = c("power2_s", "power2_t", "scale_s", "scale_t"),
    "73" = c("power2_s", "power2_t", "scale_s", "scale_t"),
    "74" = c("power2_s", "power2_t", "scale_s", "scale_t"),
    "75" = c("power2_s", "power2_t", "scale_s", "scale_t"),
    "77" = c("power2_s", "power2_t", "scale_s", "scale_t"),
    "78" = c("scale_s", "scale_t", "smooth_s", "smooth_t", "power2_s", "power2_t"),
    "84" = c("scale_s", "scale_t"),
    "85" = c("scale_s", "scale_t", "power2_s", "power2_t", "smooth_s", "smooth_t", "sep"),
    "86" = c("scale_s", "scale_t", "smooth_s", "smooth_t"),
    "87" = c("power_t", "power2_s", "power2_t", "scale_s", "scale_t", "sep", "smooth_s"),
    "88" = c("power_s", "power2_s", "power2_t", "scale_s", "scale_t", "sep", "smooth_t"),
    "89" = c("scale_s", "scale_t", "smooth_s", "smooth_t", "sep"),
    "94" = c("power_s", "power_t", "scale_s", "scale_t"),
    "96" = c("power2_s", "power2_t", "scale_s", "scale_t"),
    "111" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "power2", "scale"),
    "112" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "power2_1", "power2_12", "power2_2", "scale_1", "scale_12", "scale_2"),
    "113" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "power2", "scale"),
    "114" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "power2_1", "power2_12", "power2_2", "scale_1", "scale_12", "scale_2"),
    "115" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "power2", "scale"),
    "116" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "power2_1", "power2_12", "power2_2", "scale_1", "scale_12", "scale_2"),
    "117" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "scale_1", "scale_12", "scale_2", "smooth_1", "smooth_12", "smooth_2"),
    "118" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "scale_1", "scale_2", "smooth_1", "smooth_2"),
    "119" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "scale", "smooth"),
    "120" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "power2_1", "power2_2", "scale_1", "scale_2"),
    "121" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "scale_1", "scale_2", "smooth_1", "smooth_2"),
    "122" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "scale", "smooth"),
    "124" = c("a_1", "a_12", "a_2", "nugget_1", "nugget_2", "scale_1", "scale_2"),
    "126" = c("a_1", "a_12", "a_2", "a_21", "nugget_1", "nugget_2", "scale_1", "scale_2"),
    "128" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "scale_1", "scale_12", "scale_2", "smooth_1", "smooth_12", "smooth_2"),
    "129" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "power2_1", "power2_2", "scale_1", "scale_2"),
    "130" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "power2", "scale", "smooth"),
    "131" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "power2_1", "power2_2", "scale_1", "scale_2"),
    "132" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "power2_1", "power2_12", "power2_2", "scale_1", "scale_12", "scale_2", "smooth_1", "smooth_12", "smooth_2"),
    "134" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "power2_1", "power2_2", "scale_1", "scale_2", "smooth_1", "smooth_2"),
    "136" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "scale_1", "scale_12", "scale_2", "smooth_1", "smooth_12", "smooth_2", "power2_2"),
    "137" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "scale_1", "scale_12", "scale_2", "smooth_1", "smooth_12", "smooth_2", "power2_1", "power2_12", "power2_2"),
    "138" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "scale_1", "scale_12", "scale_2", "power1_1", "power1_12", "power1_2", "power2_1", "power2_12", "power2_2"),
    "139" = c("sill_1", "sill_2", "nugget_1", "nugget_2", "pcol", "scale_1", "scale_12", "scale_2", "power_1", "power_12", "power_2")
  )
  
  corrmodel_char <- as.character(corrmodel)
  if (corrmodel_char %in% names(param_map)) {
    return(param_map[[corrmodel_char]])
  } else {
    return(NULL)
  }
}

  #############################################################  
  #############################################################
NuisParam2 <- function(model, bivariate = FALSE, num_betas = c(1, 1), copula = NULL) {
  if (is.null(CkModel(model))) {
    stop("The name of the model is not correct.\n")
  }

  # Normalizza num_betas per univariati
  if (!bivariate && all(num_betas == c(1, 1))) {
    num_betas <- 1
  }

  #############################################################
  # Univariate case
  #############################################################
  if (!bivariate) {
    mm <- if (num_betas == 1) "mean" else {
      c("mean", paste0("mean", 1:(num_betas - 1)))
    }

    base_params <- c(mm, "nugget", "sill")

    # Aggiungi nu se copula è Clayton o SkewGaussian
    add_nu <- !is.null(copula) && copula %in% c("Clayton", "SkewGaussian")

    # Lista modelli con parametri specifici
    model_list <- list(
      simple = c("Gaussian", "Gauss", "Binomial", "Gaussian_misp_Binomial",
                 "BinomialLogistic", "Binomial2", "BinomialNeg", "Bernoulli",
                 "Poisson", "Gaussian_misp_Poisson", "Geom", "Geometric", "Wrapped",
                 "PoisBin", "PoisBinNeg", "LogGaussian", "LogGauss", "Logistic"),
      
      double_nugget = c("PoissonZIP", "Gaussian_misp_PoissonZIP", "BinomialNegZINB"),
      single_nugget_zip = c("PoissonZIP1", "Gaussian_misp_PoissonZIP1", "BinomialNegZINB1"),
      
      shape1 = c("PoissonGammaZIP1", "Gaussian_misp_PoissonGammaZIP1"),
      shape2 = c("PoissonGammaZIP", "Gaussian_misp_PoissonGammaZIP"),
      
      shape_only = c("Weibull", "weibull", "Gamma", "gamma", "LogLogistic",
                     "Loglogistic", "PoissonGamma", "PoissonWeibull",
                     "Gaussian_misp_PoissonGamma"),
      
      beta_like = c("Beta2", "Kumaraswamy2"),
      beta_full = c("Beta", "Kumaraswamy"),
      
      skew = c("SkewGaussian", "SkewGauss", "TwoPieceGaussian", "TwoPieceGauss",
               "Binomial_TwoPieceGaussian", "Binomial_TwoPieceGauss",
               "BinomialNeg_TwoPieceGaussian", "BinomialNeg_TwoPieceGauss"),
      
      skew_t = c("SkewStudentT", "TwoPieceStudentT", "Gaussian_misp_SkewStudentT"),
      bimodal = "TwoPieceBimodal",
      
      student_t = c("StudentT", "Gaussian_misp_StudentT"),
      
      tukey_h = c("Tukeyh", "tukeyh"),
      tukey_h2 = c("Tukeyh2", "tukeyh2"),
      tukey_gh = c("Tukeygh", "SinhAsinh", "TwoPieceTukeyh", "Gaussian_misp_Tukeygh")
    )

    # Assegnazione parametri
    if (model %in% model_list$simple) {
      param <- base_params
    } else if (model %in% model_list$double_nugget) {
      param <- c(mm, "nugget1", "nugget2", "pmu", "sill")
    } else if (model %in% model_list$single_nugget_zip) {
      param <- c(mm, "nugget", "pmu", "sill")
    } else if (model %in% model_list$shape1) {
      param <- c(mm, "nugget", "pmu", "sill", "shape")
    } else if (model %in% model_list$shape2) {
      param <- c(mm, "nugget1", "nugget2", "pmu", "sill", "shape")
    } else if (model %in% model_list$shape_only) {
      param <- c(mm, "nugget", "sill", "shape")
    } else if (model %in% model_list$beta_like) {
      param <- c(mm, "nugget", "sill", "shape", "min", "max")
    } else if (model %in% model_list$beta_full) {
      param <- c(mm, "nugget", "sill", "shape1", "shape2", "min", "max")
    } else if (model %in% model_list$skew) {
      param <- c(mm, "nugget", "sill", "skew")
    } else if (model %in% model_list$skew_t) {
      param <- c(mm, "df", "nugget", "sill", "skew")
    } else if (model %in% model_list$bimodal) {
      param <- c(mm, "df", "nugget", "sill", "shape", "skew")
    } else if (model %in% model_list$student_t) {
      param <- c(mm, "df", "nugget", "sill")
    } else if (model %in% model_list$tukey_h) {
      param <- c(mm, "nugget", "sill", "tail")
    } else if (model %in% model_list$tukey_h2) {
      param <- c(mm, "nugget", "sill", "tail1", "tail2")
    } else if (model %in% model_list$tukey_gh) {
      param <- c(mm, "nugget", "sill", "skew", "tail")
    } else {
      stop("Model not recognized in univariate case.")
    }

    if (add_nu) param <- c(param, "nu")
    return(param)
  }

  #############################################################
  # Bivariate case
  #############################################################
  if (bivariate) {
    mm1 <- if (num_betas[1] == 1) "mean_1" else c("mean_1", paste0("mean_1", 1:(num_betas[1] - 1)))
    mm2 <- if (num_betas[2] == 1) "mean_2" else c("mean_2", paste0("mean_2", 1:(num_betas[2] - 1)))
    mm <- c(mm1, mm2)

    model_list <- list(
      biv_simple = c("Gaussian", "Gauss", "Binomial", "BinomialLogistic", "Binomial2",
                     "BinomialNeg", "Geom", "Geometric", "Poisson", "Wrapped", "PoisBin",
                     "PoisBinNeg", "LogGaussian", "LogGauss", "Logistic"),
      
      biv_skew = c("SkewGaussian", "SkewGauss", "TwoPieceGaussian", "TwoPieceGauss"),
      biv_shape = c("Weibull", "Gamma", "LogGauss", "LogGaussian", "LogLogistic"),
      biv_t = c("StudentT", "Gaussian_misp_StudentT"),
      biv_tukeyh = c("Tukeyh", "tukeyh"),
      biv_tukeygh = c("Tukeygh", "SinhAsinh", "TwoPieceTukeyh")
    )

    if (model %in% model_list$biv_simple) {
      param <- mm
    } else if (model %in% model_list$biv_skew) {
      param <- c(mm, "skew_1", "skew_2")
    } else if (model %in% model_list$biv_shape) {
      param <- c(mm, "shape_1", "shape_2")
    } else if (model %in% model_list$biv_t) {
      param <- c(mm, "df_1", "df_2")
    } else if (model %in% model_list$biv_tukeyh) {
      param <- c(mm, "tail_1", "tail_2")
    } else if (model %in% model_list$biv_tukeygh) {
      param <- c(mm, "skew_1", "skew_2", "tail_1", "tail_2")
    } else {
      stop("Model not recognized in bivariate case.")
    }

    return(param)
  }

  stop("Unexpected error: check inputs.")
}

# Wrapper NuisParam
NuisParam <- function(model, bivariate = FALSE, num_betas = c(1, 1), copula = NULL) {
  a <- NuisParam2(model, bivariate, num_betas, copula)
  models_no_sill <- c(
    "Weibull", "Poisson", "Binomial", "Gamma", "LogLogistic",
    "BinomialNeg", "Bernoulli", "Geometric", "Gaussian_misp_Poisson",
    "PoissonGammaZIP", "PoissonGamma", "PoissonGammaZIP1", "Binary_misp_BinomialNeg",
    "PoissonZIP", "Gaussian_misp_PoissonZIP", "BinomialNegZINB",
    "PoissonZIP1", "Gaussian_misp_PoissonZIP1", "BinomialNegZINB1",
    "Beta2", "Kumaraswamy2", "Beta", "Kumaraswamy"
  )
  if (model %in% models_no_sill) {
    a <- a[a != "sill"]
  }
  return(a)
}
####################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

StartParam <- function(coordx, coordy,coordz ,coordt,coordx_dyn, corrmodel, data, distance, fcall,
                      fixed, grid,likelihood,  maxdist, neighb,maxtime, model, n, 
                      param, parscale,paramrange, radius, start, taper, tapsep, 
                      type,typereal,  weighted,copula, X,memdist,nosym)
{

####################################  
### START internal functions:
#################################### 
    replicates=1
    # Check if the correlation is bivariate
    CheckBiv <- function(corrmodel)
    {
        CheckBiv <- NULL
        if(corrmodel >= 101 & corrmodel <= 140) CheckBiv <- TRUE
        else CheckBiv <- FALSE
        return(CheckBiv)
    }
    # Check if the correlation is spatial or spatial-temporal
    CheckST <- function(corrmodel)
    {
        CheckST <- NULL
        if(corrmodel >40 & corrmodel <= 100) CheckST <- TRUE
        else  CheckST <- FALSE
        return(CheckST)
    }
    # Check the type of distances
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
##############################################################################################    
newtap<- function(coords,numcoord, coordt,numtime, distance,maxdist,maxtime,spacetime,bivariate,radius)
### using nearest.dist...
    {
      if(distance==0) method1="euclidean"
      if(distance==2||distance==1) method1="greatcircle"
if(method1=="greatcircle"){
      gb=spam::nearest.dist( x=coords,method = method1,
             delta = maxdist*360/(radius*2*pi), upper = NULL,miles=FALSE, R=radius)
      if(distance==2) gb@entries=radius*gb@entries             ##GC
      if(distance==1) gb@entries=2*radius*sin(0.5*gb@entries)  ##CH
      }

if(method1=="euclidean")
      gb=spam::nearest.dist(x=coords,method = method1,
                         delta = maxdist, upper =NULL,miles=FALSE, R=1)
      numpairs=length(gb@entries)
 ##loading only good distances..
   dotCall64::.C64("SetGlobalVar2",
        SIGNATURE = c("integer","integer",
             "double","integer","double",
            "double","double","integer",
            "integer","integer","integer","integer"),  
       numcoord,  numtime,  
       gb@entries,numpairs,srange[2],
       1,1,1, # to change for spacetime sparse
       spacetime,bivariate,1,1, INTENT =    c(rep("r",12)),
         PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
   
    colidx=gb@colindices
    rowidx=gb@rowpointers
    nozero=numpairs/(numcoord)^2
  return(list(colidx=colidx,rowidx=rowidx,numpairs=numpairs,nozero=nozero))
    }
####################################
### END Includes internal functions
####################################

    # Set the correlation and  if the correlation is space-time(spacetime=T and bivariate=F) or bivariate (F o T)  or univariate (case spacetime=F and bivariate=F)p
    corrmodel<-CkCorrModel(corrmodel)
    bivariate <- CheckBiv(corrmodel); if(bivariate) coordt=c(1,2)
    spacetime <- CheckST(corrmodel)
    isdyn=!is.null(coordx_dyn)
    space=!(spacetime||bivariate)
    ####################
    if(!bivariate)
    {
        if(is.null(X))  {X=1;num_betas=1}
        else 
          {if(is.list(X))  num_betas=ncol(X[[1]])
           else  num_betas=ncol(X) 
          }
    }
    if(bivariate)
    {
        if(is.null(X))  {X=1;num_betas=c(1,1)}
        else
                        { if(is.list(X))  num_betas=c(ncol(X[[1]]),ncol(X[[2]]))
                          else  num_betas=c(ncol(X),ncol(X)) }
    }
    ####################
    namesnuis <- NuisParam2(model,bivariate,num_betas,copula)
    ltimes=length(coordt)


    ### case gris handling coordz
    if(grid) {
              if(is.null(coordz)) { cc=as.matrix(expand.grid(coordx,coordy));coordx=cc[,1];coordy=cc[,2]; coordz=NULL}
              else {  
                    # if(is.null(coordz)) cc=as.matrix(expand.grid(coordx,coordy,0));
                     cc=as.matrix(expand.grid(coordx,coordy,coordz));
                     coordx=cc[,1];coordy=cc[,2]; coordz=cc[,3];
            

                  }
             }

    ### Set returning variables and initialize the model parameters:
    # Initialises the starting and fixed parameters' names
    error <- NULL
    ns<-NULL
    namesfixed <- namesstart <- namessim <- NULL
    numfixed <- numstart <- 0
    # Set the model, likelihood, correlation and the nuisance parameters:
    model <- CkModel(model)
    flagnuis <- NULL
    namescorr <- CorrelationPar(corrmodel)
    numparamcorr <- length(namescorr)
    paramcorr <- rep(1, numparamcorr)
    names(paramcorr) <- namescorr
    flagcorr <- NULL

    ### START settings the data structure:

    # set the coordinates sizes:
    if(is.null(coordx_dyn))  
    {
    
      if(is.null(coordy)){
                         if(ncol(coordx)==2) {
                                 coordy=coordx[,2];
                                 coordx=coordx[,1];
                                 coordz=NULL}
                        else { coordz=coordx[,3];coordy=coordx[,2];coordx=coordx[,1];}
                          }

      numcoord <- numcoordx <- numcoordy <-numcoordz <-length(coordx)

      if(bivariate && !is.null(nrow(coordx)) && !is.null(nrow(coordy))&& !is.null(nrow(coordz))  )   {  #heterotopic case
        numcoordx=nrow(coordx); 
        numcoordy=nrow(coordy);
        numcoordz=nrow(coordz);
        numcoord=numcoordy+numcoordx+numcoordz}
      
      ns<-rep(numcoord,ltimes)
    }
    else     ##  dynamic coords
    {
       env <- new.env()
       coords=do.call(rbind,args=c(coordx_dyn),envir = env) 

       if(is.list(X))  X=do.call(rbind,args=c(X),envir = env)
       if(ncol(coords)==2) ns=lengths(coordx_dyn)/2 
       if(ncol(coords)==3) ns=lengths(coordx_dyn)/3
       
       if(ncol(coords)==2) {coordx <- coords[,1]; coordy <- coords[,2];coordz=NULL}
       if(ncol(coords)==3) { coordx <- coords[,1]; coordy <- coords[,2];coordz <- coords[,3]}
       numcoord <- numcoordx <- numcoordy <-numcoordz <- length(coordx)
    }


   if(!space && is.null(coordx_dyn)) {coordx=rep(coordx,ltimes);coordy=rep(coordy,ltimes);coordz=rep(coordz,ltimes);}
    
    NS=cumsum(ns)
    if(!space)   NS=c(0,NS)[-(length(ns)+1)]

    # initialize tapering variables:
    tapering=ia=idx=ja=colidx=rowidx=integer(1)
    nozero<-NULL
    tapmodel=0
    cutoff <- FALSE
    distance<-CheckDistance(distance)
    ### END settings the data structure
    # START code for the simulation procedure


##################################################################################################################
    if(fcall=="Fitting"){
        ### Parameters' settings:
        nuisance=nuisance1=nuisance2=NULL
        likelihood <- CkLikelihood(likelihood)
        type <- CkType(type)
     if((!bivariate&&num_betas==1)||(bivariate&all(num_betas==c(1,1))))
     {
          if(model %in% c(1,10,12,18,9,20,13,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,39,40,41,42,46,47,48,50))   ### continous model 
          {
           if(!bivariate) {
                           mu <- mean(unlist(data))
                           if(any(type==c(1, 3, 7,8,4))){    # Checks the type of likelihood
                           if(is.list(fixed)) fixed$mean <- mu# Fixs the mean
                           else fixed <- list(mean=mu)}
                           nuisance <- c(mu, 0, var(c(unlist(data))))

                           if(likelihood==2 && (CkType(typereal)==5 || CkType(typereal)==7) ) tapering <- 1
                           if(model %in% c(10,29,31,32))         nuisance <- c(nuisance,0)
                           if(model %in% c(18,20,27,37,38,40,41))      nuisance <- c(0,nuisance,0)
                            if(model %in% c(39))      nuisance <- c(0,0,nuisance,0)
                           if(model %in% c(21,24,12,26,34,35,46,47,48))   nuisance <- c(0,nuisance)
                           #if(model %in% c(23,28,33))  nuisance <- c(0,0,0,nuisance)
                           if(model %in% c(23,28,33))  nuisance <- c(0,0,0,nuisance,0,0)
                           if(model %in% c(42,50))  nuisance <- c(0,nuisance,0,0)
                       }
             if(bivariate) {
                           if(is.null(coordx_dyn)) { mu1 <- mean(data[1,]); mu2 <- mean(data[2,])}
                           else                   { mu1 <- mean(data[[1]]); mu2 <- mean(data[[2]])}
                           if(any(type==c(1, 3, 7, 8,4))) {   # Checks the type of likelihood
                           if(is.list(fixed)) {fixed$mean_1 <- mu1;fixed$mean_2<- mu2}
                           else               {fixed <- list(mean_1=mu1,mean_2=mu2)  }
                                                         }
                           nuisance <- c(mu1,mu2)
                           if(model %in% c(10,29,31,32))  {nuisance <- c(nuisance,0.1,0.2)}
                           if(model %in% c(26,46,47,48,42,50))  {nuisance <- c(nuisance,0.1,0.2)}
                           if(model %in% c(21))  {nuisance <- c(nuisance,0.1)}
                          
                           if(likelihood==2 && (CkType(typereal)==5 || CkType(typereal)==7)) tapering <- 1
                 }
        }

if(model %in% c(11,13,14,15,16,19,17,30,45,49,51,52,53,54,56,58)){                                                        #discrete
            mu=0
            if(any(type==c(1, 3, 7,8,4))){    # Checks the type of likelihood
                           if(is.list(fixed)) 
                                              { fixed$mean <- mu}# Fixs the mean}
                                            else      {fixed <- list(mean=mu)}
                           
                           }                        

            nuisance <- c(mu, 0, 1)

            if(model==45) nuisance <- c(mu, 0, 0,0,1)
            if(model==53||model==56) nuisance <- c(mu, 0,0,1)
            if(model==58) nuisance <- c(mu, 0,0,0,1)
        }
 if(model %in% c(43,44)) nuisance <- c(0, 0, 0,0, 1)
 if(model %in% c(57)) nuisance <-    c(0, 0, 0,0, 0,1)

}


 if((!bivariate&&num_betas>1)||(bivariate&&num_betas[1]>1&&num_betas[2]>1) )
 {

   if(model %in% c(1,10,12,18,9,20,13,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,39,40,41,42,46,47,48,50))  ### continous models
{
    if(!bivariate)
    {
         if(any(type==c(1, 3, 7,8,4)))# Checks the type of likelihood
            if(is.list(fixed)) {
                               mu <- mean(unlist(data));fixed$mean <- mu# Fixs the mean
                               for(i in 1:(num_betas-1)) fixed[[paste("mean",i,sep="")]]=1 
                               }
            else               {mu <- mean(unlist(data));fixed <- list(mean=mu)}
            for(i in 1:num_betas) nuisance=c(nuisance,1);
            nuisance=c(nuisance,0,var(c(unlist(data))))
             if(model %in% c(10,29,31,32))        nuisance=c(nuisance,1)  
             if(model %in% c(21,24,12,26,34,35,46,47,48))  nuisance=c(nuisance,1) 
             if(model %in% c(18,20,27,37,38,40,41))     nuisance=c(1,nuisance,1) 
             if(model %in% c(39))     nuisance=c(1,1,nuisance,1) 
          #  if(model %in% c(23,28,33))         nuisance=c(nuisance,1,1,1)  
            if(model %in% c(23,28,33))         nuisance=c(nuisance,1,1,1,1,1)  
            if(model %in% c(42,50))  nuisance=c(nuisance,1,1,1) 
    }
    if(bivariate)
    {
            if(any(type==c(1, 3, 7,8,4)))# Checks the type of likelihood
            if(is.list(fixed)) 
        {
                if(!is.list(data))
                {
                                mu1 <- rowMeans(unlist(data))[1];fixed$mean_1 <- mu1# Fixs the mean
                                mu2 <- rowMeans(unlist(data))[2];fixed$mean_2 <- mu2# Fixs the mean
                                for(i in 1:(num_betas[1]-1)) fixed[[paste("mean_1",i,sep="")]]=1 
                                for(i in 1:(num_betas[2]-1)) fixed[[paste("mean_2",i,sep="")]]=1  # fixed$meani=1
                }
                  if(is.list(data))
                {
                                mu1 <- mean(data[[1]]);fixed$mean_1 <- mu1# Fixs the mean
                                mu2 <- mean(data[[2]]);fixed$mean_2 <- mu2# Fixs the mean
                                for(i in 1:(num_betas[1]-1)) fixed[[paste("mean_1",i,sep="")]]=1 
                                for(i in 1:(num_betas[2]-1)) fixed[[paste("mean_2",i,sep="")]]=1  # fixed$meani=1
                }
        }
            else  fixed <- list(mean_1=mu1,mean_2=mu2)
            for(i in 1:num_betas[1]) nuisance1=c(nuisance1,1);
            for(i in 1:num_betas[2]) nuisance2=c(nuisance2,1);
            nuisance=c(nuisance1,nuisance2 )
             if(model %in% c(10,29,31,32))        nuisance=c(nuisance,1,1)  
             if(model %in% c(21,24,12,26,34,35,46,47,48,42,50))  nuisance=c(nuisance,1,1) 
             if(model %in% c(18,20,27,37,38,40,41))     nuisance=c(1,nuisance,1) 
              if(model %in% c(39))     nuisance=c(1,1,nuisance,1) 
            #if(model %in% c(23,28,33))         nuisance=c(nuisance,1,1,1)
            if(model %in% c(23,28,33))         nuisance=c(nuisance,1,1,1,1,1)  

    }
}  ### end continous models
     if(model %in% c(2,11,13,14,15,16,19,17,30,49,51,52,54,43,45,53,56,57,58)) {  # discrete models

        if(any(type==c(1, 3, 7,8,4)))# Checks the type of likelihood
            if(is.list(fixed)) {
                               mu <- mean(unlist(data));fixed$mean <- mu# Fixs the mean
                               for(i in 1:(num_betas-1)) fixed[[paste("mean",i,sep="")]]=1  # fixed$meani=1
                               }
            else               {mu <- mean(unlist(data));fixed <- list(mean=mu)}

         if(model %in% c(2,11,13,14,15,16,19,17,30,49,51,52,54))   nuisance <- c(0,rep(0,num_betas-1) ,0, 1) 

       
  
     if(model %in% c(43,45)) nuisance <- c(0,rep(1,num_betas-1) ,0, 0,0,1)
     if(model %in% c(57)) nuisance <-    c(0,rep(1,num_betas-1) ,0, 0,0,0,1)

     if(model %in% c(53,56)) nuisance <- c(0,rep(1,num_betas-1) , 0,0,1)
     if(model %in% c(58)) nuisance <- c(0,rep(1,num_betas-1) , 0,0,0,1)

        }
     #
}

 if(!is.null(copula)) {if(copula=="Clayton"||copula=="SkewGaussian") nuisance=c(nuisance,2)}
# Update the parameter vector     
        names(nuisance) <- namesnuis
        namesparam <- sort(c(namescorr, namesnuis))
        param <- c(nuisance, paramcorr)
        param <- param[namesparam]
        numparam <- length(param)
        flag <- rep(1, numparam)
        namesflag <- namesparam
        names(flag) <- namesflag
        # Update the parameters with fixed values:
     
        if(!is.null(fixed)){
            fixed <- unlist(fixed)
            namesfixed <- names(fixed)
            numfixed <- length(namesfixed)

            flag[pmatch(namesfixed, namesflag)] <- 0
            param <- param[-pmatch(namesfixed, namesparam)]
                          numparamcorr <- numparamcorr-sum(namesfixed %in% namescorr)
            namesparam <- names(param)
            numparam <- length(param)   

        }
        else { }

        flagcorr <- flag[namescorr]
        flagnuis <- flag[namesnuis]
        # Update the parameters with starting values:
        if(!is.null(start)){
            start <- unlist(start)
            namesstart <- names(start)

            #if(any(type == c(1, 3, 7))){
                 if(any(type==c(1, 3, 7,8,4))){    # Checks the type of likelihood

                if(!bivariate) {   # univariate case
                       if(any(model==c(1,10,12,18,20,9,13,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,39,40,41,42,46,47,48,50,
                        11,14,15,16,19,17,30,45,49,51,52,54,53,56,58)))
                       if(any(namesstart == 'mean'))  start <- start[!namesstart == 'mean']
                       if(num_betas>1)
                       for(i in 1:(num_betas-1)) 
                       {  
                         if(any(namesstart == paste("mean",i,sep="")))  {
                              namesstart <- names(start) 
                              if(any(model==c(1,10,12,18,20,9,13,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,39,40,41,42,46,47,48,50,
                        11,14,15,16,19,17,30,45,49,51,52,54,53,56,58)))
                                                 start <- start[!namesstart == paste("mean",i,sep="")]
                            }
                       }
                }
                if(bivariate) {          
                                  if(any(namesstart == 'mean_1'))  start <- start[!namesstart == 'mean_1']        
                                  if(any(namesstart == 'mean_2'))  {namesstart <- names(start) ; start <- start[!namesstart == 'mean_2']}
                      
                       if(num_betas[1]>1)
                       for(i in 1:(num_betas[1]-1)) {  if(any(namesstart == paste("mean_1",i,sep="")))  {namesstart <- names(start) ; 
                       if(any(model==c(1,10,12,18,20,9,13,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,39,40,41,42,46,47,48,50)))
                                                 start <- start[!namesstart == paste("mean_1",i,sep="")]}}            
                       if(num_betas[2]>1)
                       for(i in 1:(num_betas[2]-1)) {  if(any(namesstart == paste("mean_2",i,sep="")))  {namesstart <- names(start) ; 
                       if(any(model==c(1,10,12,18,20,9,13,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,39,40,41,42,46,47,48,50)))
                                                 start <- start[!namesstart == paste("mean_2",i,sep="")]}}  
                                  }
                }

            namesstart <- names(start)
            numstart <- length(start)
            param[pmatch(namesstart,namesparam)] <- start
            }
        ### set the scale of the parameters:
        # Insert here!
        # set the range of the parameters if its the case
        paramrange=TRUE
        #if(paramrange) paramrange <- SetRangeParam(namesparam, numparam)
        #else 
        paramrange <- list(lower=NULL, upper=NULL)
        ### Set the data format:
        if(!space){ # setting spam indexes
            if(spacetime) numtime <- ltimes
            if(bivariate) numtime <- 2
                }
        else{              #    setting spam indexes
            numtime <- 1
            coordt <- 0
            data <- matrix(data, ncol=numcoord, nrow=replicates)
            }
      
       #if((typereal=="Tapering"&&type=="Tapering")||(typereal=="Tapering1"&&type=="Tapering1")||(typereal=="Tapering2"&&type=="Tapering2")){
     
        if(typereal=="Tapering"||typereal=="Tapering1"||typereal=="Tapering2"){
        tapering<-1

        if(!space){
        idx<-integer((numcoord*numtime)^2)
        ja<-integer((numcoord*numtime)^2)
        ia<-integer(numcoord*numtime+1)}
          else {idx=ja=ia=0}
        tapmodel<-CkCorrModel(taper)
                }
              ###### ojo aca!! if conditional then I used nn2 using "all" the indeces
        if(likelihood==1&&is.numeric(maxdist)&&is.null(neighb))   
        neighb= min(numcoord*numtime,500)

                                             
    }
# END code for the fitting procedure
    
##################################################################################################################
# START code for the simulation procedure
    if(fcall=="Simulation"){

        neighb=NULL;likelihood=2
        namesnuis <- sort(unique(c(namesnuis,NuisParam2("Gaussian",bivariate,num_betas,copula))))
        param <- unlist(param)
        numparam <- length(param)
        namesparam <- names(param)

        if(!bivariate) if(any(model!=c(43,45,53,56,57,58)))  namessim <- c("mean","sill","nugget","scale",namescorr[!namescorr=="scale"])
       if(bivariate)  namessim <- c("mean_1","mean_2","scale",namescorr[!namescorr=="scale"])  

        if(spacetime) numtime <- ltimes
        else {numtime <- 1; coordt <- 0}
        if(bivariate) numtime <- 2

        if((typereal=="Tapering"&&type=="Tapering")||(typereal=="Tapering1"&&type=="Tapering1")||(typereal=="Tapering2"&&type=="Tapering2")){
                 tapering<-1
                 nt=numcoord*numtime
               if(!space){  
                idx<-integer(nt^2)
                ja<-integer(nt^2)
                ia<-integer(nt+1)
                 }
                else {idx=ja=ia=0}
                tapmodel<-CkCorrModel(taper)
        }
        K=neighb
}  # END code for the simulation procedure
#####################################################################################
numpairs <- integer(1); srange <- double(1); trange <- double(1)

if(typereal=="Independence"){ maxdist=NULL;maxtime=NULL;K=neighb}
  
#################
distC=FALSE
if(!tapering) { if(is.null(neighb)&is.numeric(maxdist)) distC=TRUE  }### just for maxdist parameter
################
  
    if(is.null(maxdist)) srange<-c(srange,double(1)) else {srange<-c(srange,as.double(maxdist))}                # cutoff<-TRUE
    if(is.null(maxtime)) trange<-c(trange,double(1)) else {trange<-c(trange,as.double(maxtime))}                # cutoff<-TRUE
    isinit <- as.integer(1)
    if(is.null(tapsep))  tapsep=c(0.5,0.5)
    else  {if(length(tapsep)==1) tapsep=c(tapsep,0)}
    mem=FALSE
    if(tapering||memdist)  { mem=TRUE }   #### NB

    if(mem&&!tapering)  
      {        
         
                nn=numcoord*numtime
                if(spacetime&&isdyn)  nn=sum(ns)
                if(is.null(neighb)){
                    if(typereal=="Independence") colidx=rowidx=0
                    else         colidx=rowidx=integer(nn*(nn-1)/2)}
           
    }
    if(bivariate) {
    if(!srange[1]&&!srange[2])  srange=c(srange,0,0)
    if(is.na(srange[3])) srange[3]=srange[2];
    if(is.na(srange[4])) srange[4]=srange[2];}
    ###
    if(CheckSph(corrmodel))   radius=1
    ###
    aa=double(5);for(i in 1:length(tapsep)) aa[i]=tapsep[i];tapsep=aa

if(fcall=="Fitting"&likelihood==2&!is.null(neighb)) mem=FALSE # Vecchia gp case
if(fcall=="Fitting"&likelihood==2||fcall=="Simulation") mem=FALSE 
if(tapering) mem=TRUE

###################### using new "spam" with neighdist(just for space)#########################################################

if(tapering&space){
    cc=cbind(coordx, coordy,coordz);numcoord=nrow(cc);numtime=length(coordt)
    atap=newtap(cc,numcoord, coordt,numtime, distance,maxdist,maxtime,spacetime,bivariate,radius)
    ja=colidx=atap$colidx
    ia=rowidx=atap$rowidx
    numpairs=atap$numpairs
    nozero=atap$nozero
    K=NULL
}
else{          # all the rest
##############################################################
## loading distances in memory using brute force C routine ###
#############################################################
### aca paso solo para  simular o maximum likelihood o variogram 
### o si hay CL with  maxdist!!!
  
if(distC||fcall=="Simulation"||(fcall=="Fitting"&likelihood==2)||(fcall=="Fitting"&typereal=="GeoWLS")) {
if(fcall=="Fitting"&mem==TRUE&(!space)&!tapering)   {vv=length(NS); numcoord=NS[vv]+ns[vv]} # number of space time point in the case of coordxdyn

if(is.null(coordz)) coordz=double(numcoordx*numtime) ## is it necessary?

srange[which(srange==Inf)]=1e+50;trange[which(trange==Inf)]=1e+50
gb=.C('SetGlobalVar',as.integer(bivariate), as.double(coordx), as.double(coordy),as.double(coordz), as.double(coordt),as.integer(grid),ia=as.integer(ia),idx=as.integer(idx),  #7
           isinit=as.integer(isinit),ja=as.integer(ja), as.integer(mem), as.integer(numcoord),as.integer( numcoordx),  as.integer(numcoordy), as.integer(numcoordz), #6
           numpairs=as.integer(numpairs), as.double(radius),as.double(srange), as.double(tapsep),  as.integer(spacetime), #5
            as.integer(numtime),as.double(trange), as.integer(tapering), as.integer(tapmodel),as.integer(distance),as.integer(weighted), #6
           colidx= as.integer(colidx),rowidx= as.integer(rowidx), # 2
            as.integer(ns), as.integer(NS), as.integer(isdyn),
      PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)

  #if(!sum(coordz)) coordz=NULL

rm(colidx);rm(rowidx)
if(type=="Tapering") {rm(idx);rm(ja);rm(ia)}
##
## number  of selected pairs
numpairs <- gb$numpairs
## indexes for composite 
    colidx=gb$colidx
    rowidx=gb$rowidx
    colidx <- colidx[1:numpairs]
    rowidx  <- rowidx[1:numpairs]
## indexes for tapering 
    idx <- gb$idx
    ja <- gb$ja 
    ia <- gb$ia;
    isinit <- gb$isinit
    nozero <- numpairs/(numcoord*numtime)^2
    idx <- idx[1:numpairs]
    ja  <- ja[1:numpairs]
    K=neighb
  
}
#######################################################################
else   
###############################################################
################### loading distances in memory using nabor package 
#### it works when CL  using neighb  or maxdist AND neighb 
#############################################################
{ 
if(typereal!="Independence") {
  ########################## 
if(distance==0) distance1="Eucl";
if(distance==2) distance1="Geod";
if(distance==1) distance1="Chor";

if(all(neighb==0.5)) neighb=NULL ## ojo!!
if(sum(!is.finite(maxdist))) maxdist=NULL

if(space)   #  spatial case
{
##########################################
  
  K=neighb
  #x=cbind(coordx, coordy)
   x=cbind(coordx, coordy,coordz)

  sol=GeoNeighIndex(coordx=x,distance=distance1,maxdist=maxdist,neighb=K,radius=radius)

 ###    deleting symmetric indexes with associate distances
 if(nosym){
  aa=GeoNosymindices(cbind(sol$colidx,sol$rowidx),sol$lags)
  sol$rowidx=c(aa$xy[,1])
  sol$colidx=c(aa$xy[,2])
  sol$lags=c(aa$d) }

###
  nn = length(sol$lags)
  sol$lagt=0
  gb=list(); gb$colidx=sol$colidx;
             gb$rowidx=sol$rowidx ;
             gb$numpairs=nn

             

  ## loading space distances in memory 
  mmm=1;ttt=1
if(weighted)  mmm=max(sol$lags)
  
    ss <- dotCall64::.C64("SetGlobalVar2",
           SIGNATURE = c("integer", "integer",
            "double",
             "integer", "double",
                         "double", 
                         "integer", "double", "integer", "integer",
                         "integer", "integer"),
           as.integer(numcoord),
           as.integer(numtime),
           as.double(sol$lags),
           as.integer(nn),
           as.double(mmm),
           as.double(sol$lagt),
           as.integer(nn),
           as.double(ttt),
           as.integer(spacetime),
           as.integer(bivariate),
           as.integer(1),
           as.integer(1),
           INTENT = rep("r", 12), 
           PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)

} 



##############################################   
if(spacetime)   #  space time  case
{      
  K=neighb
  x=cbind(coordx, coordy,coordz)
 
  sol=GeoNeighIndex(coordx=x[1:numcoord,],
    coordx_dyn=coordx_dyn,
    coordt=coordt,distance=distance1,maxdist=maxdist,neighb=K,maxtime=maxtime,radius=radius)

    #sol=GeoNeighIndex(coordx=x,
    #coordx_dyn=coordx_dyn,
    #coordt=coordt,distance=distance1,maxdist=maxdist,neighb=K,maxtime=maxtime,radius=radius)

 # ###    deleting symmetric indexes with associate distances #unuseful
  if(nosym){ aa=GeoNosymindices(cbind(sol$colidx,sol$rowidx),sol$lags)
             sol$rowidx=c(aa$xy[,1])
             sol$colidx=c(aa$xy[,2])
             sol$lags=c(aa$d)
           }

  gb=list(); gb$colidx=sol$colidx;
             gb$rowidx=sol$rowidx ;
             nn=length(gb$colidx)
             gb$numpairs=nn
  ## loading space time distances in memory   
  mmm=1;ttt=1
if(weighted) { mmm=max(sol$lags) ;ttt=max(sol$lagt)}

  ss=dotCall64::.C64("SetGlobalVar2",
        SIGNATURE = c("integer","integer",
             "double",
             "integer","double",
              "double",
              "integer","double",
            "integer","integer","integer","integer"),  
       numcoord,  numtime, 
        sol$lags,
        nn,mmm, 
      sol$lagt,
      nn,ttt,
       spacetime,bivariate,1,1, INTENT =    c(rep("r",12)),
         PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
} 
##############################################  
if(bivariate)   # bivariate case 
{ 
  K=neighb
  x=cbind(coordx, coordy,coordz)

  sol=GeoNeighIndex(coordx=x, coordx_dyn=coordx_dyn, distance=distance1,maxdist=maxdist,neighb=K,maxtime=maxtime,radius=radius,bivariate=TRUE)
  ###    deleting symmetric indexes with associate distances
  if(nosym){ aa=GeoNosymindices(cbind(sol$colidx,sol$rowidx),sol$lags)
             sol$rowidx=c(aa$xy[,1])
             sol$colidx=c(aa$xy[,2])
             sol$lags=c(aa$d)}
  gb=list(); gb$colidx=sol$colidx;
             gb$rowidx=sol$rowidx ;
             #gb$first=sol$first
             #gb$second=sol$second
             nn=length(gb$colidx)
             gb$numpairs=nn
## loading space time distances in memory   
  mmm=1
if(weighted) { mmm=max(sol$lags)}
  
  
  ss=dotCall64::.C64("SetGlobalVar2",
        SIGNATURE = c("integer","integer",
             "double",
             "integer","double",
              "double","integer","double",
            "integer","integer","integer","integer"),  
       numcoord,  2,  
       sol$lags,
       nn,mmm, 
       1,nn,1,spacetime,bivariate,sol$first,sol$second,
        INTENT =    c(rep("r",12)),
         PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)

  
  
} #### end bivariate case

    numpairs <- gb$numpairs
    colidx=gb$rowidx 
    rowidx=gb$colidx
    idx <- 0;ja <- 0;ia <- 0
    isinit <- 1
    nozero <- numpairs/(numcoord*numtime)^2
    idx <- 0;ja  <- 0
}  
############################################## 
# end neighboord case
##############################################
if(is.null(coordt)) coordt=1
}}

########################################################################################
########################################################################################
########################################################################################
########################################################################################
    ### Returned list of objects:
    return(list(bivariate=bivariate,coordx=coordx,coordy=coordy,coordz=coordz,coordt=coordt,corrmodel=corrmodel,
                colidx = colidx ,rowidx=rowidx,
                data=data,distance=distance,
                error=error,flagcorr=flagcorr,flagnuis=flagnuis,fixed=fixed,likelihood=likelihood,
                lower=paramrange$lower,model=model,n=n,namescorr=namescorr,namesfixed=namesfixed,
                namesnuis=namesnuis,namesparam=namesparam,namessim=namessim,namesstart=namesstart,ns=ns,NS=NS,
                num_betas=num_betas,
                numcoord=numcoord,numcoordx=numcoordx,numcoordy=numcoordy,neighb=K,
                numfixed=numfixed,numpairs=numpairs,numparam=numparam,numparamcorr=numparamcorr,
                numstart=numstart,numtime=numtime,param=param,
                setup=list(                ## setup is a list
                ia=ia,idx=idx,ja=ja,nozero=nozero,tapmodel=tapmodel,tapsep=tapsep),  radius=radius,                            ## with tapered matrix informations
                spacetime=spacetime,srange=srange,start=start,upper=paramrange$upper,type=type,
                trange=trange,weighted=weighted,X=X))
}