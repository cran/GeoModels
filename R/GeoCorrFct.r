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
           #p=.C('VectCorrelation', corr=double(numlags*numlagt), as.integer(corrmodel), as.double(lags),
            #                 as.integer(numlags), as.integer(numlagt), as.double(mu),as.integer(model)
             #                ,as.double(nuisance),as.double(param), as.double(lagt),as.integer(N), PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
                            
                    }
        else    {
                nn=numlags*4
                         #    p=.C('VectCorrelation_biv', corr=double(nn),vario=double(nn), as.integer(corrmodel), as.double(lags),
                         #    as.integer(numlags), as.integer(numlagt),  as.double(mu),as.integer(model),as.double(nuisance), as.double(param),
                         #    as.double(lagt), as.integer(N),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)

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
 
    if(is.null(CkCorrModel (corrmodel))) stop("The name of the correlation model  is not correct\n")
    if(is.null(CkModel(model)))   stop("The name of the  model  is not correct\n")
    if(!is.numeric(x)) stop("Distances must be numeric\n")
    if(sum(x<0)>=1) stop("Distances must be positive\n")
    if(!is.list(param)) stop("param  must be a list\n")

    spacetime<-CheckST(CkCorrModel(corrmodel))
    bivariate<-CheckBiv(CkCorrModel(corrmodel))

    if(!bivariate) {if(is.null(param$sill)) param$sill=1}

if(is.null(t)) t=0
   
mu=0;nuisance=0
mm=0
num_beta=c(1,1)

nx=length(x)
if(spacetime) nt=length(t)
else nt=1

num_betas=c(1,1)
if(!bivariate) if(sum((names(param)=='mean'))==0) param$mean=0 # adding mean if missing

mu=as.numeric(param$mean)
  ## selecting nuisance mean annd corr parameters
      if(!bivariate){
        parcorr <- c(param)[CorrelationPar(CkCorrModel(corrmodel))]
        nuisance <- c(param)[NuisParam(model,FALSE,num_betas)]
        sel=substr(names(nuisance),1,4)=="mean"
        mm=as.numeric(nuisance[sel])
        nuisance=nuisance[!sel]
        }
      if(bivariate){
        if(!covariance) param$sill_1=param$sill_2=1
        parcorr <- c(param)[CorrelationPar(CkCorrModel(corrmodel))]
        nuisance <- c(param)[NuisParam(model,bivariate,num_betas)]
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
             else   cova <- correlation*(1-as.numeric(nuisance["nugget"]))+as.numeric(nuisance["nugget"])*I(x==0)
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
                              vs=sill*(nu)/(nu-2)
                              cc=((nu-2)*gamma((nu-1)/2)^2*Re(hypergeo::hypergeo(0.5,0.5 ,nu/2 ,correlation^2))*correlation1)/(2*gamma(nu/2)^2)
                              
                              cova=cc;
                               }
                  }

      ##########################################
   if(model=="SkewStudentT")        { if(bivariate) {}
                        else {
                              
                 if(length(t)>1) correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)
                 else            correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))+as.numeric(nuisance["nugget"])*I(x==0)

                              nu=1/as.numeric(nuisance['df']);sill=as.numeric(nuisance['sill']);sk=as.numeric(nuisance['skew'])
                              skew2=sk*sk; f=(nu-1)/2; w=sqrt(1-skew2);y=correlation1;
                              CorSkew=(2*skew2/(pi*w*w+skew2*(pi-2)))*(sqrt(1-y*y)+y*asin(y)-1)+w*w*correlation/(w*w+skew2*(1-2/pi)) ;
                              mm=sqrt(nu)*gamma(f)*sk/(sqrt(pi)*gamma(nu/2));
                cc=(pi*(nu-2)*gamma(f)^2/(2*(pi*gamma(nu/2)^2-skew2*(nu-2)*gamma(f)^2)))*(Re(hypergeo::hypergeo(0.5,0.5,nu/2,y*y))*((1-2*skew2/pi)*CorSkew+2*skew2/pi)-2*skew2/pi);
                               vs=sill*(nu/(nu-2)-mm*mm);
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
                              sill=as.numeric(nuisance['sill'])
                              corr=correlation
                              corr[corr>=0.9999999999]=0.9999999999
                             x1=1-(1-corr^2)*hr; y1=1-(1-corr^2)*hl
                             x2=(1-hr)^2-(corr*hr)^2;y2=(1-hl)^2-(corr*hl)^2
                             g=1-hl-hr+(1-corr^2)*hl*hr
                             h1=sqrt(1-corr^2/(x1^2))+(corr/x1)*asin(corr/x1);
                             h2=sqrt(1-corr^2/(y1^2))+(corr/y1)*asin(corr/y1)
                             h3=sqrt(1-corr^2/(x1*y1))+sqrt(corr^2/(x1*y1))*asin(sqrt(corr^2/(x1*y1)))
                             p1=x1*h1/(2*pi*(x2)^(3/2))+corr/(4*(x2)^(3/2))
                             p2=y1*h2/(2*pi*(y2)^(3/2))+corr/(4*(y2)^(3/2))
                             p3=-(x1*y1)^(1/2)*h3/(2*pi*(g)^(3/2))+corr/(4*(g)^(3/2))

                             mm=(hr-hl)/(sqrt(2*pi)*(1-hl)*(1-hr))
                             vv1=0.5*(1-2*hl)^(-3/2)+0.5*(1-2*hr)^(-3/2)-(mm)^2
                             cc=(p1+p2+2*p3-mm^2)/vv1  # correlation
                             vs=as.numeric(nuisance['sill'])*vv1
                            
                             cova=cc; 
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
                               vs=sill*(exp(2*eta2/(1-2*tail))-2*exp(eta2/(2*(1-2*tail)))+1)/(eta2*sqrt(1-2*tail))-mu*mu;
                               A1=exp(a*eta2/(1-tail*a));
                               A2=2*exp(0.5*eta2*  (1-tail*(1-rho2))  / (u*u- tail2*rho2)  );
                               A3=eta2*sqrt(u*u- rho2*tail*tail);
                               cc=((A1-A2+1)/A3-mu*mu)/vs;
                               
                               cova=cc;
                            } 
                        if(tail<=1e-05&&abs(eta)>1e-05){
                              vs=sill*( -exp(eta^2)+exp(eta^2*2))*eta^(-2)
                              cc= (( -exp(eta^2)+exp(eta^2*(1+rho)))*eta^(-2))/vs
                              
                              cova=cc;
                            } 
                        if(tail>1e-05&&abs(eta)<=1e-05){
                              vs=  sill*(1-2*tail)^(-1.5)     
                              cc=(-rho/((1+h*(tail-1))*(-1+tail+tail*rho)*(1+tail*(-2+tail-tail*rho^2))^0.5))/vs
                              
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
ga=as.numeric(nuisance['shape2'])
eta=log1p(-(1+exp(mu))^(-ga))/log(0.5)
}
mm=eta*beta(1+1/ga,eta)
sill=eta*beta(1+2/ga,eta)-mm^2
###
NN=length(correlation);
res=double(NN)
bb=.C("corr_kuma_vec",as.double(correlation),as.double(eta),as.double(ga), 
     res=as.double(res),
    as.integer(NN),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
rho=bb$res
cova=rho;
vs=sill
 }
}

##########################


if(model == "SinhAsinh") {
  if(bivariate) {
    # Implementazione bivariata da definire
  } else {
    # Estrazione parametri
    nugget <- as.numeric(nuisance['nugget'])
    d <- as.numeric(nuisance['tail'])     # deve essere > 0
    e <- as.numeric(nuisance['skew'])     # può essere qualsiasi valore
    sill <- as.numeric(nuisance['sill'])
    
    # Applicazione nugget effect
    if(length(t) > 1) {
      correlation <- correlation * (1 - nugget) + nugget * I(A[,1] == 0 & A[,2] == 0)
    } else {
      correlation <- correlation * (1 - nugget) + nugget * I(x == 0)
    }
    corr_base <- correlation
    
    # Imposta dinamicamente max_coeff in base a |e| e d
    max_coeff <- if (abs(e) > 2 || d > 3) 10 else if (abs(e) > 1 || d > 2) 8 else 6
    j_vec <- 1:max_coeff
    
    # Calcolo robusto delle costanti
    sqrt_8pi <- sqrt(8 * pi)
    sqrt_32pi <- sqrt(32 * pi)
    exp_0.25 <- exp(0.25)
    
    # Gestione robusta delle funzioni iperboliche per evitare overflow
    sinh_e_d <- if (abs(e / d) > 20) {
      sign(e / d) * 0.5 * exp(abs(e / d))
    } else {
      sinh(e / d)
    }
    
    cosh_2e_d <- if (abs(2 * e / d) > 20) {
      0.5 * exp(abs(2 * e / d))
    } else {
      cosh(2 * e / d)
    }
    
    # Calcolo robusto delle funzioni di Bessel con gestione errori
    besselK_1 <- suppressWarnings(besselK(0.25, (d + 1) / (2 * d)))
    besselK_2 <- suppressWarnings(besselK(0.25, abs(1 - d) / (2 * d)))
    besselK_3 <- suppressWarnings(besselK(0.25, (d + 2) / (2 * d)))
    besselK_4 <- suppressWarnings(besselK(0.25, abs(2 - d) / (2 * d)))
    
    # Controllo validità dei risultati di Bessel
    if(any(is.na(c(besselK_1, besselK_2, besselK_3, besselK_4))) || 
       any(is.infinite(c(besselK_1, besselK_2, besselK_3, besselK_4)))) {
      warning("Problemi nel calcolo delle funzioni di Bessel, usando approssimazioni")
      # Fallback con approssimazioni per casi estremi
      besselK_1 <- ifelse(is.na(besselK_1) || is.infinite(besselK_1), 0, besselK_1)
      besselK_2 <- ifelse(is.na(besselK_2) || is.infinite(besselK_2), 0, besselK_2)
      besselK_3 <- ifelse(is.na(besselK_3) || is.infinite(besselK_3), 0, besselK_3)
      besselK_4 <- ifelse(is.na(besselK_4) || is.infinite(besselK_4), 0, besselK_4)
    }
    
    # Calcolo momento e varianza
    mm <- sinh_e_d * exp_0.25 * (besselK_1 + besselK_2) / sqrt_8pi
    vs <- cosh_2e_d * exp_0.25 * (besselK_3 + besselK_4) / sqrt_32pi - 0.5 - mm^2
    
    # Verifica che vs sia positiva
    if(vs <= 0) {
      warning("Varianza negativa o nulla, usando valore minimo")
      vs <- 1e-10
    }
    
    # Funzione integrand ottimizzata
    integrand <- function(z, alpha, kappa, j, r) {
      z_sq <- z^2
      z_pow <- if (j - 2 * r == 0) 1 else z^(j - 2 * r)
      aa <- z + sqrt(z_sq + 1)
      
      # Gestione robusta delle potenze per evitare overflow
      alpha_kappa <- alpha / kappa
      pow_kappa <- 1 / kappa
      
      aa_pow1 <- if (abs(log(aa) * pow_kappa) > 700) {
        sign(aa)^pow_kappa * exp(700 * sign(pow_kappa))
      } else {
        aa^pow_kappa
      }
      
      aa_pow2 <- if (abs(log(aa) * (-pow_kappa)) > 700) {
        sign(aa)^(-pow_kappa) * exp(700 * sign(-pow_kappa))
      } else {
        aa^(-pow_kappa)
      }
      
      exp_term1 <- if (abs(alpha_kappa) > 700) {
        sign(alpha_kappa) * exp(700)
      } else {
        exp(alpha_kappa)
      }
      
      exp_term2 <- if (abs(-alpha_kappa) > 700) {
        sign(-alpha_kappa) * exp(700)
      } else {
        exp(-alpha_kappa)
      }
      
      return(exp(-z_sq / 2) * z_pow * (exp_term1 * aa_pow1 - exp_term2 * aa_pow2))
    }
    
    # Cache per gli integrali
    II_cache <- new.env(hash = TRUE)
    
    # Funzione II ottimizzata con cache e gestione errori
    II <- function(alpha, kappa, j, r) {
      key <- paste(alpha, kappa, j, r, sep = "_")
      if (exists(key, envir = II_cache)) {
        return(get(key, envir = II_cache))
      }
      
      result <- tryCatch({
        integrate(integrand, lower = -Inf, upper = Inf,
                  alpha = alpha, kappa = kappa, j = j, r = r,
                  rel.tol = 1e-6, abs.tol = 1e-10)$value
      }, error = function(err) {
        warning(paste("Errore integrazione per j =", j, "r =", r, ". Usando range limitato."))
        tryCatch({
          integrate(integrand, lower = -10, upper = 10,
                    alpha = alpha, kappa = kappa, j = j, r = r,
                    rel.tol = 1e-4, abs.tol = 1e-8)$value
        }, error = function(err2) {
          warning(paste("Fallback fallito per j =", j, "r =", r, ". Usando 0."))
          return(0)
        })
      })
      
      assign(key, result, envir = II_cache)
      result
    }
    
    # Funzione coeff_j ottimizzata
    coeff_j <- function(alpha, kappa, j) {
      max_r <- floor(j / 2)
      if (max_r < 0) return(0)
      
      rr <- 0:max_r
      II_vals <- numeric(length(rr))
      
      for (i in seq_along(rr)) {
        II_vals[i] <- II(alpha, kappa, j, rr[i])
      }
      
      # Controllo validità dei risultati
      if (any(is.na(II_vals)) || any(is.infinite(II_vals))) {
        warning(paste("Valori non validi nell'integrale per j =", j))
        II_vals[is.na(II_vals) | is.infinite(II_vals)] <- 0
      }
      
      terms <- II_vals * (-1)^rr / (2^(rr + 1) * gamma(rr + 1) * gamma(j - 2 * rr + 1))
      res <- gamma(j + 1) * sum(terms) / sqrt(2 * pi)
      
      return(res)
    }
    
    # Calcolo coefficienti con controllo validità
    coeffs <- sapply(j_vec, function(j) coeff_j(e, d, j))
    
    # Verifica coefficienti validi
    if (any(is.na(coeffs)) || any(is.infinite(coeffs))) {
      warning("Coefficienti non validi rilevati, applicando correzioni")
      coeffs[is.na(coeffs) | is.infinite(coeffs)] <- 0
    }
    
    # Funzione di correlazione SAS ottimizzata
    corrsas <- function(rho) {
      if (abs(rho) < 1e-12) return(1)
      
      # Filtra coefficienti significativi
      non_zero_idx <- which(abs(coeffs) > 1e-12)
      if (length(non_zero_idx) == 0) return(1)
      
      # Controllo convergenza della serie
      if (abs(rho) >= 1 && any(j_vec[non_zero_idx] > 0)) {
        warning("Possibile divergenza della serie per |rho| >= 1")
      }
      
      numerator <- sum(coeffs[non_zero_idx]^2 * rho^j_vec[non_zero_idx] / 
                         gamma(j_vec[non_zero_idx] + 1))
      
      result <- numerator / vs
      
      # Controllo risultato valido
      if (is.na(result) || is.infinite(result)) {
        warning("Risultato non valido nella correlazione, usando fallback")
        return(1)
      }
      
      return(result)
    }
    
    # Applica la funzione di correlazione
    corr <- if (length(corr_base) > 1) {
      sapply(corr_base, corrsas)
    } else {
      corrsas(corr_base)
    }
    
    # Scaling finale
    vs <- vs * sill
    cova <- corr
  }
}
##########################################
 if(model=="TwopieceT")        { if(bivariate) {}
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
  if(model=="TwopieceTukeyh")        { if(bivariate) {}
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
 if(model=="TwopieceGauss")        { 
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
                        auxcorr= (gamma(1+1/sh))^2/((gamma(1+2/sh))-(gamma(1+1/sh))^2)
                        cc=auxcorr*(Re(hypergeo::hypergeo(-1/sh, -1/sh, 1,correlation^2)) -1)
                        cova=cc;
                     }
                    }
##########################################
if(model=="Loglogistic")    { if(bivariate) {}  
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
                    vs<-(exp(vvar)-1)
                    cova=cc; 
                     }
                  }
#############################################
if(model=="Binomial"||model=="BinomialNeg"||model=="Geometric"||model=="BinomialnegZINB")
   {
                    if(bivariate) {}
                    if(!bivariate) {          
                           pp=pnorm(mu)
                           if(model=="Binomial") vs=min(n)*pp*(1-pp)
                           if(model=="BinomialNeg")      vs=(n)*(1-pp)/pp^2;
                           if(model=="BinomialnegZINB") { 
                                     pg=pnorm(as.numeric(nuisance['pmu']))
                                     vs=n*(1-pp)*(1-pg)*(1+n*pg*(1-pp)) /pp^2
                                               }
                          ## en este caso el nugget ya se computo
                           if(length(t)>1)  cova=correlation     +as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)
                           else cova=correlation     +as.numeric(nuisance["nugget"])*I(x==0)
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
                       if(length(t)>1) correlation=(1-as.numeric(nuisance['nugget']))*correlation  +as.numeric(nuisance["nugget"])*I(A[,1]==0&A[,2]==0)
                       else correlation=(1-as.numeric(nuisance['nugget']))*correlation  +as.numeric(nuisance["nugget"])*I(x==0)
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
                       cova=vv*cc
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

