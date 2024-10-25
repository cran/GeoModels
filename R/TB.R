#########################################################
#########################################################
############ functions for TB############################
#########################################################
#########################################################

spectral_density_1dR <- function(param, corrmodel, u_vec = seq(-1,1,l=100)){


    av <- param$scale
    nu1v <- param$smooth
    params_other <- param$power2
    if(is.null(params_other)) params_other <- 0
    norm_u <- u_vec
    N <- length(u_vec)
    model=CkCorrModel(corrmodel)

  
   result=dotCall64::.C64("spectral_density_1d",
        SIGNATURE = c("double","integer","double", "double","double","integer","double"),
                     norm_u, N, av, params_other, nu1v, model, simu1 = dotCall64::numeric_dc(N),
        INTENT =    c("r","r","r", "r","r", "r","w"),PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$simu1

  return(result)
}

##################################################### 
######### bivarate case #############################
#####################################################
tbm2d <- function(coord, coordt, param, corrmodel,L,bivariate){

  N=1; n <- dim(coord)[1]; d <- 2
    if(corrmodel == "Matern"){
       CC=1; N=1
       a <- as.numeric(param['scale']);  nu1 <- as.numeric(param['smooth'])
       a0 = a; nu0 =nu1
       parameters <- list("CC" = CC, "a" = a,"nu1" = nu1,"nu2" = 0 )
       P <- 1;  vtype = 0}

  if (corrmodel=="Bi_matern"||corrmodel=="Bi_Matern"){
    corrmodel <- "Matern"
    a <- matrix(0,2,2)
    a[1,1] <- as.numeric(param['scale_1'])
    a[2,2] <- as.numeric(param['scale_2'])
    a[1,2] <- a[2,1] <- as.numeric(param['scale_12'])
    
    nu1 <- matrix(0,2,2)
    nu1[1,1] <- as.numeric(param['smooth_1'])
    nu1[2,2] <- as.numeric(param['smooth_2'])
    nu1[1,2] <- nu1[2,1] <- as.numeric(param['smooth_12'])
    CC <- matrix(0,2,2)
    CC[1,1] <- CC[2,2] <- 1
    CC[1,2] <- CC[2,1] <- as.numeric(param['pcol'])
    a0 <- min(a)
    nu0 <- min(nu1)
    parameters <- list("C" = CC, "a" = a,"nu1" = nu1,"nu2" =matrix(0,2,2) )
    P<-2; vtype = 0  

  }
   parametersg <- list("a" = a0,"nu1" = nu0)

      A <- matrix(0, P, P*L*N); B <- matrix(0, P, P*L*N)
      #S <- ceiling(1e7*runif(3)); set.seed(S[1])
      G <- matrix(rgamma(P*L*N, nu0, scale = 1),P*L*N,d); #set.seed(S[2])
      u <- matrix(rnorm(P*L*N*d), P*L*N, d)/sqrt(G*2)/a0/(2*pi); #set.seed(S[3])
      phi <- 2*pi*runif(P*L*N)
      sequen <- c(seq(0,n-0.5, by = ceiling(1e6/P/N)),n)

    m = c()
   for (i in 1:(length(sequen)-1)){ m1 <- sequen[i+1]-sequen[i]; m = c(m1,m)} 
simu11 = as.numeric( rep(0,N*P*sum(m)*(length(sequen)-1)))

result=dotCall64::.C64("for_c",
         SIGNATURE = c("integer","double","double", "double",
                     "double","integer","integer","integer","integer",
                     "double","double","double","double","double",
                     "integer","integer","integer",
                     "double","double","integer","integer","double",
                     "double"), 
                         d,c(a),c(nu1), c(CC), c(parameters$nu2),P,N,L,CkCorrModel(corrmodel),
                         c(u),a0,nu0,c(A),c(B),c(sequen),length(sequen),n,coord,phi,vtype,m,
                         simu1=dotCall64::numeric_dc(length(simu11)),L,
         INTENT =    c("r","r","r","r","r","r","r","r","r", "r",
                       "r","r","r","r","r","r","r","r","r", "r",
                       "r", "rw","r"),
             PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$simu1
  simu =  matrix(result,n,P)
  return(simu)
}
#######################################################
######### univariate case #############################
#######################################################
tbm2d_uni <- function(coord, coordt, param, corrmodel, L, bivariate){
  
  N=1; n <- dim(coord)[1]; d <- 2


#########################
if(corrmodel=="GenWend")  
{
rep=(gamma(param$power2+2*param$smooth+1)/gamma(param$power2))^(-1/(1+2*param$smooth)) # inverse parametrization
param$scale=param$scale*rep
corrmodel="GenWend_Matern"
param$power2=1/param$power2
}
#########################
if(corrmodel=="Kummer")  
{
  rep= sqrt(2*(param$power2+1)) # inverse parametrization
param$scale=param$scale/rep
corrmodel="Kummer_Matern"
}
#########################


LIM1=8 # after this limit then GenWend_Matern is a matern approx..
LIM2=8 # after this limit then kummer_Matern is a matern approx..
  if(corrmodel == "GenWend_Matern"&&(1/param$power2)>=LIM1)
  {
    corrmodel = "Matern"
    param$smooth=param$smooth+0.5
    param$power2=NULL
  }
   if(corrmodel == "Kummer_Matern"&&param['power2']>=LIM2)
  {
    corrmodel = "Matern"
    param$power2=NULL
  }

 ################################### 
  if(corrmodel == "Matern"){
    CC=as.double(param['sill']); N=1
    a <- as.double(param['scale']);  
    nu1 <- as.double(param['smooth'])
    other <- as.double(0)
    a0 = a; nu0 = nu1
    parameters <- list("CC" = CC, "a" = a,"nu1" = nu1, "nu2" = 0, "other" = other)
    P <- 1;  vtype = 0
    parametersg <- list("a" = a0, "nu1" = nu0, other = other, CC = CC)
    u <- frequency_sampler(L, parametersg, corrmodel)
    model_num <- CkCorrModel(corrmodel)
  }
  ##################################################
  #if(corrmodel == "Kummer"){
  #  CC=as.double(param['sill']); N=1
  #  a <- as.double(param['scale']);  
  #  nu1 <- as.double(param['smooth']);
  #  other <- as.double(param['power2'])
  #  a0 = a; nu0 = nu1
  #  a0 <- a <- param$scale
  #  parameters <- list("CC" = CC, "a" = a,"nu1" = nu1,"nu2" = 0, "other" = other)
  #  P <- 1;  vtype = 0
  #  parametersg <- list("a" = a0, "nu1" = nu0, other = other, CC = CC)
  #  u <- frequency_sampler(L, parametersg, corrmodel)
  #  model_num <- CkCorrModel(corrmodel)
  #}
   if(corrmodel == "Kummer_Matern"){
    CC=as.double(param['sill']); N=1
    a <- as.double(param['scale']);  
    nu1 <- as.double(param['smooth']);
    other <- as.double(param['power2'])
    a0 = a; nu0 = nu1
    a0 <- a <- param$scale    
    parameters <- list("CC" = CC, "a" = a,"nu1" = nu1,"nu2" = 0, "other" = other)
    P <- 1;  vtype = 0
    parametersg <- list("a" = a0, "nu1" = nu0, other = other, CC = CC)
    u <- frequency_sampler(L, parametersg, corrmodel)
    model_num <- CkCorrModel(corrmodel)
  }
    ##################################################
 # if(corrmodel == "GenWend"){
 #   CC <- as.double(param['sill']); N=1;
 #   nu1 <- as.numeric( param['smooth'] )
 #   other <- as.numeric(param['power2'])
 #   a <- as.numeric(param['scale'])
 #   a0 = a; nu0 = nu1
 #   parameters <- list("CC" = CC, "a" = a,"nu1" = nu1, "nu2" = 0, "other" = other)
 #   P <- 1;  vtype = 0
 #   parametersg <- list("a" = a0, "nu1" = nu0, other = other, CC = CC)
 #   u <- frequency_sampler(L, parametersg, corrmodel)
 #   model_num <- CkCorrModel(corrmodel)
 # }
  if(corrmodel == "GenWend_Matern"){
    CC <- as.double(param['sill']); N=1;
    nu1 <- as.numeric( param['smooth'] )

    other <- 1/as.numeric(param['power2'])

    a <- as.numeric(param['scale'])
    a0 = a; nu0 = nu1
    parameters <- list("CC" = CC, "a" = a,"nu1" = nu1, "nu2" = 0, "other" = other)
    P <- 1;  vtype = 0
    parametersg <- list("a" = a0, "nu1" = nu0, other = other, CC = CC)
    u <- frequency_sampler(L, parametersg, corrmodel)
    model_num <- CkCorrModel(corrmodel)
  }

    ##################################################
  if(!is.null(u)){  
  A <- matrix(0, P, P*L*N); B <- matrix(0, P, P*L*N)
  phi <- 2*pi*runif(P*L*N)
  sequen <- c(seq(0,n-0.5, by = ceiling(1e6/P/N)),n)
  m = c()
  for (i in 1:(length(sequen)-1)){ m1 <- sequen[i+1]-sequen[i]; m = c(m1,m)}
  simu11 = as.numeric( rep(0,N*P*sum(m)*(length(sequen)-1)))

  result=dotCall64::.C64("for_c",
                         SIGNATURE = c("integer","double", "double","double",
                                       "double","integer","integer","integer","integer","double",
                                       "double","double","double","double",
                                       "integer","integer","integer",
                                       "double","double","integer","integer","double",
                                       "double","double"),
                         d_v = d, a_v = c(a), nu1_v = c(nu1), C_v = c(CC), nu2_v = c(parameters$nu2), 
                         P = P, N = N, L = L, model = model_num,
                         u = c(u), a0 = a0, nu0 = nu0, A = c(A), B = c(B), sequen = c(sequen),
                         largo_sequen = length(sequen), n = n, coord = coord, phi = phi, vtype = vtype,
                         m1 = m,
                         simu1=dotCall64::numeric_dc(length(simu11)), L1 = L, other,
                         INTENT = c("r", "r", "r", "r", "r","r","r","r", "r",
                                    "r", "r", "r", "r", "r", "r","r","r","r", "r",
                                    "r", "r", "rw","r","r"),
                         PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$simu1
 
  simu =  matrix(sqrt(CC)*result,n,P)
  }
  else
  {simu=u=parametersg=NULL }
  return( list( simu = simu, freqs = u, parametersg = parametersg) )
}
######################################################################
######################################################################
frequency_sampler <- function(L, parametersg, corrmodel){
  
  if(corrmodel == "Matern"){
    #Simulando frecuencias de una densidad espectral Matern
    G <- matrix(rgamma(L, parametersg$nu1, scale = 1),L,2)
    ini.sample <- matrix(rnorm(L*2), L, 2)/sqrt(G*2)/parametersg$a/(2*pi)
  }
  
  #if(corrmodel == "GenWend"){
  #  #Simulando frecuencias desde la densidad espectral
  #  param_geomodels.int <- list(scale = parametersg$a, smooth = parametersg$nu, sill = 1,
  #                              power2 = parametersg$other, nugget = 0)
  #  Q.particles <- 15
  #   dff=0.1#2*(1.5 + parametersg$nu)-1
    
  #  u_proposal <- cbind( rt(Q.particles*L, df = dff), rt(Q.particles*L, df = dff))
  #  u_norm <- sqrt( rowSums(u_proposal^2) )
  #  w0 <- dt(u_proposal[,1], df = dff)*dt(u_proposal[,2], df = dff)
  #  w1 <- pmax( (spectral_density_1dR(param_geomodels.int, corrmodel, u_vec = u_norm) ), 0)
  #  www <- w1/w0
  #   if(any(is.na(www))){ ini.sample=NULL}
  #  else{
  #  u_index <- sample(1:(Q.particles*L), size = L, replace = T, prob = www)
  #  ini.sample <- matrix(u_proposal[u_index,], ncol = 2 )
  #}
  #}

if(corrmodel == "GenWend_Matern"){
    #Simulando frecuencias desde la densidad espectral
    param_geomodels.int <- list(scale = parametersg$a, smooth = parametersg$nu, sill = 1,
                                power2 = parametersg$other, nugget = 0)
    Q.particles <- 15

    #dff <- 2*(1.5 + parametersg$nu)
     dff=0.05#2*(1.5 + parametersg$nu)-1
    
    u_proposal <- cbind( rt(Q.particles*L, df = dff), rt(Q.particles*L, df = dff))
    u_norm <- sqrt( rowSums(u_proposal^2) )
    w0 <- dt(u_proposal[,1], df = dff)*dt(u_proposal[,2], df = dff)
    w1 <- pmax( (spectral_density_1dR(param_geomodels.int, corrmodel, u_vec = u_norm) ), 0)
    www <- w1/w0
     if(any(is.na(www))){ ini.sample=NULL}
    else{
    u_index <- sample(1:(Q.particles*L), size = L, replace = T, prob = www)
    ini.sample <- matrix(u_proposal[u_index,], ncol = 2 )
  }
  }
  
  #if(corrmodel == "Kummer"){
  #  param_geomodels.int <- list(scale = parametersg$a, smooth = parametersg$nu, sill = 1, 
   #                             power2 = parametersg$other, nugget = 0)

   # Q.particles <- 15
   # dff=0.2#min(0.05,parametersg$other-1)
   # u_proposal <- cbind( rt(Q.particles*L,  df=dff), rt(Q.particles*L, df=dff))
   # u_norm <- sqrt( rowSums(u_proposal^2) )
   # w0 <- dt(u_proposal[,1], df=dff)*dt(u_proposal[,2], df=dff)
   # w1 <- pmax( (spectral_density_1dR(param_geomodels.int, corrmodel, u_vec = u_norm) ), 0)
   # www <- w1/w0
   # if(any(is.na(www))){ ini.sample=NULL}
   # else
   # {
   # u_index <- sample(1:(Q.particles*L), size = L, replace = T, prob = www)
   # ini.sample <- matrix(u_proposal[u_index,], ncol = 2 )
   # }
  #}

if(corrmodel == "Kummer_Matern"){
    #Simulando frecuencias desde la densidad espectral
    param_geomodels.int <- list(scale = parametersg$a, smooth = parametersg$nu, sill = 1, 
                                power2 = parametersg$other, nugget = 0)

    Q.particles <- 15
    dff=0.2#min(0.05,parametersg$other-1)
    u_proposal <- cbind( rt(Q.particles*L,  df=dff), rt(Q.particles*L, df=dff))
    u_norm <- sqrt( rowSums(u_proposal^2) )
    w0 <- dt(u_proposal[,1], df=dff)*dt(u_proposal[,2], df=dff)
    w1 <- pmax( (spectral_density_1dR(param_geomodels.int, corrmodel, u_vec = u_norm) ), 0)
    www <- w1/w0
    if(any(is.na(www))){ ini.sample=NULL}
    else{
    u_index <- sample(1:(Q.particles*L), size = L, replace = T, prob = www)
    ini.sample <- matrix(u_proposal[u_index,], ncol = 2 )
  }
  }
  
  u <- ini.sample
  return(u)
}