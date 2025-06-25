####################################################
### File name: GeoCorrFct_Cop_Fixed.r
####################################################

GeoCorrFct_Cop <- function(x, t=NULL, corrmodel, model="Gaussian", copula="Gaussian", 
                          distance="Eucl", param, radius=6371, n=1, 
                          covariance=FALSE, variogram=FALSE) {

############################################################################
## VERSIONE CORRETTA - OTTIMIZZAZIONI CONSERVATIVE
############################################################################

############################################################################
## Funzioni C (mantenute identiche all'originale)
############################################################################
biv_unif_CopulaClayton <- function(a, b, c, d) {
  sol <- .C("biv_unif_CopulaClayton_call", as.double(a), as.double(b),
           as.double(c), as.double(d), ress = double(1),
           PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
  return(exp(sol$ress))
}

biv_unif_CopulaGauss <- function(a, b, c) {
  sol <- .C("biv_unif_CopulaGauss_call", as.double(a), as.double(b),
           as.double(c), ress = double(1),
           PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
  return(sol$ress)
}

# Vectorizzazione (mantenuta identica)
biv_unif_CopulaClayton <- Vectorize(biv_unif_CopulaClayton, vectorize.args=c("a","b"))
biv_unif_CopulaGauss <- Vectorize(biv_unif_CopulaGauss, vectorize.args=c("a","b"))

# Funzioni argomento (identiche all'originale)
arg_clayton <- function(x, y, rho, nu, q1, q2) return(q1(x) * q2(y) * biv_unif_CopulaClayton(x, y, rho, nu))
arg_gaussian <- function(x, y, rho, q1, q2) return(q1(x) * q2(y) * biv_unif_CopulaGauss(x, y, rho))

############################################################################
## Funzione correlazione copula CORRETTA
############################################################################
corr_copula <- function(rho, copula, q1, q2, e1, e2, v1, v2, nu) {
  
  # OTTIMIZZAZIONE 1: Pre-calcolo di sqrt(v1*v2) e e1*e2 (fuori dal loop)
  sqrt_v1v2 <- sqrt(v1 * v2)
  e1_e2 <- e1 * e2
  
  if (copula == "Clayton") {
    # OTTIMIZZAZIONE 2: Pre-allocazione del vettore risultato
    res <- numeric(length(rho))
    
    for (i in seq_along(rho)) {
      if (rho[i] > 0.9999 & rho[i] <= 1) {
        res[i] <- 1
      } else {
        # Calcolo integrale (mantenuto identico)
        integral_result <- pracma::integral2(arg_clayton, 0, 1, 0, 1, 
                                           rho=rho[i], nu=nu, q1=q1, q2=q2)$Q
        res[i] <- (integral_result - e1_e2) / sqrt_v1v2
      }
    }
    return(res)
  }
  
  if (copula == "Gaussian") {
    res <- numeric(length(rho))
    
    for (i in seq_along(rho)) {
      if (rho[i] > 0.9999 & rho[i] <= 1) {
        res[i] <- 1
      } else {
        integral_result <- pracma::integral2(arg_gaussian, 0, 1, 0, 1, 
                                           rho=rho[i], q1=q1, q2=q2)$Q
        res[i] <- (integral_result - e1_e2) / sqrt_v1v2  
      }
    }
    return(res)
  }
}

############################################################################
## Funzione CorrelationFct (identica all'originale)
############################################################################
CorrelationFct <- function(bivariate, corrmodel, lags, lagt, numlags, numlagt, mu, model, nuisance, param, N) {
  if (!bivariate) {
    p <- dotCall64::.C64('VectCorrelation', SIGNATURE = c("double","integer","double","integer", "integer",
                        "double","integer","double","double","double","integer"),  
                        corr=dotCall64::vector_dc("double",numlags*numlagt), corrmodel, lags,numlags,numlagt, mu,model,nuisance,param,lagt,N,
                        INTENT = c("rw","r","r","r", "r", "r","r", "r","r", "r","r"),
                        PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
  } else {
    p <- dotCall64::.C64('VectCorrelation_biv',SIGNATURE = c("double","double","integer", "double",
                        "integer", "integer","double","integer","double","double",
                        "double","integer"),  
                        corr=dotCall64::vector_dc("double",numlags*4),vario=numlags*4, corrmodel, lags,numlags,numlagt,mu,model,nuisance, param,lagt, N,
                        INTENT = c("rw","r","r","r", "r", "r","r", "r","r", "r","r"),
                        PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE) 
  }
  return(p$corr)
}

############################################################################
## CORPO PRINCIPALE - CORREZIONI APPLICATE
############################################################################

# Validazione input (identica)
call <- match.call()
if(is.null(CkCorrModel(corrmodel))) stop("The name of the correlation model is not correct\n")
if(is.null(CkModel(model))) stop("The name of the model is not correct\n")
if(!is.numeric(x)) stop("Distances must be numeric\n")
if(sum(x<0)>=1) stop("Distances must be positive\n")

spacetime <- CheckST(CkCorrModel(corrmodel))
bivariate <- CheckBiv(CkCorrModel(corrmodel))

# Inizializzazione (identica)
mu <- 0
nuisance <- 0
mm <- 0
num_beta <- c(1,1)
nx <- length(x)

if(spacetime) {
  nt <- length(t) 
  A <- as.matrix(expand.grid(x,t)) 
} else {
  t <- 0
  nt <- 1
}

num_betas <- c(1,1)
if(sum((names(param)=='mean'))==0) param$mean <- 0

mu <- as.numeric(param$mean)

# Selezione parametri (identica)
if(!bivariate) {
  parcorr <- c(param)[CorrelationPar(CkCorrModel(corrmodel))]
  nuisance <- c(param)[NuisParam(model,FALSE,num_betas)]
  sel <- substr(names(nuisance),1,4)=="mean"
  mm <- as.numeric(nuisance[sel])
  nuisance <- nuisance[!sel]
}
if(bivariate) {
  parcorr <- c(param)[CorrelationPar(CkCorrModel(corrmodel))]
  nuisance <- c(param)[NuisParam(model,FALSE,num_betas)]
}

# Calcolo correlazione
correlation <- CorrelationFct(bivariate, CkCorrModel(corrmodel), x, t, nx, nt, mu,
                             CkModel(model), nuisance, parcorr, n)

# OTTIMIZZAZIONE 3: Calcolo nugget una sola volta
nugget_val <- as.numeric(nuisance['nugget'])
cc <- correlation * (1 - nugget_val)

# Aggiunta nugget effect (identica alla logica originale)
if(length(t) > 1) {
  correlation <- correlation * (1 - nugget_val) + nugget_val * I(A[,1]==0 & A[,2]==0)
} else {   
  correlation <- correlation * (1 - nugget_val) + nugget_val * I(x==0)
}

############################################################################
## Gestione modelli (corretta)
############################################################################

if(model == "Beta2") { 
  if(!bivariate) {
    delta <- as.numeric(nuisance["shape"])
    mu1 <- mu2 <- 1/(1+exp(-mm))  # CORREZIONE: assegnazione corretta
    
    # OTTIMIZZAZIONE 4: Pre-calcolo parametri
    a1 <- a2 <- mu1 * delta
    b1 <- b2 <- (1 - mu1) * delta
    
    # Funzioni quantile (identiche)
    q1 <- function(x) qbeta(x, a1, b1)
    q2 <- function(x) qbeta(x, a2, b2)
    
    # Momenti (identici)
    e1 <- a1/(a1+b1)
    e2 <- a2/(a2+b2)   
    v1 <- a1*b1/((a1+b1)^2*(a1+b1+1))
    v2 <- a2*b2/((a2+b2)^2*(a2+b2+1))
    
    min_val <- as.numeric(param['min'])
    max_val <- as.numeric(param['max'])
    dd <- max_val - min_val
    vs <- sqrt(v1*v2) * dd^2
  }
}        

if(model == "Gaussian") { 
  if(!bivariate) {
    sill <- as.numeric(nuisance["sill"])
    v1 <- v2 <- sill
    e1 <- e2 <- mm
    q1 <- function(x) qnorm(x, e1, sqrt(sill))
    q2 <- function(x) qnorm(x, e2, sqrt(sill))
    vs <- sill
  }
}       

############################################################################
## Calcolo copula (corretto)
############################################################################
if(copula == "Clayton") {
  cova <- corr_copula(cc, "Clayton", q1, q2, e1, e2, v1, v2, as.numeric(param['nu']))
}
if(copula == "Gaussian") {
  cova <- corr_copula(cc, "Gaussian", q1, q2, e1, e2, v1, v2, 0)
}

# Calcolo finale (identico)
if(!covariance) vs <- 1
res <- cova * vs
if(variogram) res <- vs * (1 - cova)

# Costruzione risultato (identica)
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