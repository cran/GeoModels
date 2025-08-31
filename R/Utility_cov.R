######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################


### decomposition of a square  matrix
MatDecomp <- function(mtx, method) {
  # Controlli preliminari veloci
  if (!is.matrix(mtx) || anyNA(mtx) || any(!is.finite(mtx))) return(FALSE)
  
  if (method == "cholesky") {
    if (!isSymmetric(mtx)) return(FALSE)
    mat.decomp <- tryCatch(
      FastGP::rcppeigen_get_chol(mtx),
      error = function(e) FALSE
    )
    if (isFALSE(mat.decomp)) return(FALSE)
    d <- diag(mat.decomp)
    if (any(d <= 0, na.rm = TRUE)) return(FALSE) 
  } else if (method == "svd") {
    mat.decomp <- tryCatch(
      svd(mtx),
      error = function(e) FALSE
    )
    if (isFALSE(mat.decomp)) return(FALSE)
    d <- mat.decomp$d
    if (any(d <= 0, na.rm = TRUE)) return(FALSE)
    if (!is.finite(sum(log(d)))) return(FALSE)
  }
  mat.decomp
}

MatSqrt <- function(mat.decomp, method) {
  if (method == "cholesky") return(mat.decomp)
  if (method == "svd") {
    d_sqrt <- sqrt(mat.decomp$d)
    return(t(mat.decomp$v * rep(d_sqrt, each = nrow(mat.decomp$v))))
  }
  stop("Invalid method")
}
MatInv <- function(mtx) {
  # Aggiunto controllo preliminare
  if (!is.matrix(mtx) || anyNA(mtx)) return(FALSE)
  tryCatch(
    FastGP::rcppeigen_invert_matrix(mtx),
    error = function(e) FALSE
  )
}
MatLogDet <- function(mat.decomp, method) {
  if (method == "cholesky") {
    d <- diag(mat.decomp)
    if (any(d <= 0, na.rm = TRUE)) return(NA)
    2 * sum(log(d))
    
  } else if (method == "svd") {
    d <- mat.decomp$d
    if (any(d <= 0, na.rm = TRUE)) return(NA)
    sum(log(d))
    
  } else {
    NA
  }
}

# utility function for geokrig
getInvC <- function(covmatrix, CC, mse = TRUE) {
    if (!covmatrix$sparse) {
      U <- tryCatch({
        FastGP::rcppeigen_get_chol(covmatrix$covmatrix)
      }, error = function(e) {
        stop("Covariance matrix is not positive definite")
      })
      vec <- forwardsolve(U, CC)
      Invc <- forwardsolve(U, vec, transpose = TRUE)
      mse_val <- if (mse) as.numeric(crossprod(vec)) else NULL
    } else {
      cc <- if (spam::is.spam(covmatrix$covmatrix)) covmatrix$covmatrix else spam::as.spam(covmatrix$covmatrix)
      U <- tryCatch({
        spam::chol.spam(cc)
      }, error = function(e) {
        stop("Covariance matrix is not positive definite")
      })
      vec <- spam::forwardsolve(U, CC)
      Invc <- spam::backsolve(U, vec)
      mse_val <- if (mse) as.numeric(spam::crossprod.spam(vec)) else NULL
    }
    list(a = Invc, b = mse_val)
  }

  
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

corrsas <- function(corr, skew, tail, max_coeff = NULL) {
     d <- tail
     e <- skew
    if (is.null(max_coeff)) {
        stability_indicator <- abs(e/d) + abs(d - 1) + max(abs(corr))
        is_stable <- stability_indicator < 2.0
        
        if (is_stable) {
            base_terms <- 6
            adjustment <- min(ceiling(stability_indicator), 6)
            max_coeff <- base_terms + adjustment
        } else {
            base_terms <- 8
            adjustment <- min(ceiling(stability_indicator - 2), 6)
            max_coeff <- base_terms + adjustment
        }
        max_coeff <- min(max_coeff, 12)
        max_coeff <- max(max_coeff, 6)
    }
    
    # Pre-calcolo di costanti
    sqrt_8pi <- sqrt(8 * pi)
    sqrt_32pi <- sqrt(32 * pi)
    sqrt_2pi <- sqrt(2 * pi)
    exp_0.25 <- exp(0.25)
    
    # Calcolo funzioni di Bessel con fallback
    tryCatch({
        besselK_1 <- besselK(0.25, (d + 1)/(2*d))
        besselK_2 <- besselK(0.25, (1 - d)/(2*d))
        besselK_3 <- besselK(0.25, (d + 2)/(2*d))
        besselK_4 <- besselK(0.25, (2 - d)/(2*d))
    }, error = function(e) {
        return(rep(NA, length(corr)))
    })
    
    # Calcolo mm e vv
    sinh_e_d <- sinh(e/d)
    cosh_2e_d <- cosh(2*e/d)
    mm <- sinh_e_d * exp_0.25 * (besselK_1 + besselK_2) / sqrt_8pi
    vv <- cosh_2e_d * exp_0.25 * (besselK_3 + besselK_4) / sqrt_32pi - 0.5 - mm^2
    
    if (abs(vv) < .Machine$double.eps) {
        return(rep(NA, length(corr)))
    }
    
    # Pre-calcolo di j_vec e gamma terms (evita ricalcoli)
    j_vec <- 1:max_coeff
    gamma_j_plus_1 <- gamma(j_vec + 1)  # Pre-calcolato
    
    # Funzione integranda semplificata
    integrand <- function(z, alpha, kappa, j, r) {
        z_sq <- z^2
        if (j - 2*r == 0) {
            z_pow <- 1
        } else {
            z_pow <- z^(j - 2*r)
        }
        aa <- z + sqrt(z_sq + 1)
        exp_term <- exp(-z_sq/2 + alpha/kappa)
        exp_term * z_pow * (aa^(1/kappa) - exp(-2*alpha/kappa) * aa^(-1/kappa))
    }
    
    # Cache globale per evitare ricreazione
   # if (!exists("II_global_cache", envir = .GlobalEnv)) {
   #     assign("II_global_cache", new.env(hash = TRUE), envir = .GlobalEnv)
   # }
   # II_cache <- get("II_global_cache", envir = .GlobalEnv)

    II_cache <- .GeoModels_env$II_global_cache
    
    # Funzione II ottimizzata 
    II <- function(alpha, kappa, j, r) {
        key <- paste(round(alpha, 8), round(kappa, 8), j, r, sep = "_")
        if (exists(key, envir = II_cache)) {
            return(get(key, envir = II_cache))
        }
        
        val <- tryCatch({
            integrate(integrand, lower = -Inf, upper = Inf, 
                     alpha = alpha, kappa = kappa, j = j, r = r,
                     rel.tol = 1e-6, # Leggermente meno rigoroso ma più veloce
                     subdivisions = 100)$value
        }, error = function(e) 0)
        
        assign(key, val, envir = II_cache)
        val
    }
    
    # Calcolo coefficienti ottimizzato
    coeffs <- numeric(max_coeff)
    
    for (j in j_vec) {
        max_r <- floor(j/2)
        if (max_r < 0) {
            coeffs[j] <- 0
            next
        }
        
        # Vettorizzazione del loop interno
        rr <- 0:max_r
        II_vals <- sapply(rr, function(r_val) II(e, d, j, r_val))
        
        # Calcolo gamma terms vettorizzato
        gamma_r_plus_1 <- gamma(rr + 1)
        gamma_j_minus_2r_plus_1 <- gamma(j - 2*rr + 1)
        
        # Controllo validità
        valid_gamma <- is.finite(gamma_r_plus_1) & is.finite(gamma_j_minus_2r_plus_1) & 
                      gamma_r_plus_1 > 0 & gamma_j_minus_2r_plus_1 > 0
        
        if (any(valid_gamma)) {
            terms <- II_vals[valid_gamma] * (-1)^rr[valid_gamma] / 
                    (2^(rr[valid_gamma] + 1) * gamma_r_plus_1[valid_gamma] * gamma_j_minus_2r_plus_1[valid_gamma])
            
            if (is.finite(gamma_j_plus_1[j]) && gamma_j_plus_1[j] > 0) {
                coeffs[j] <- gamma_j_plus_1[j] * sum(terms) / sqrt_2pi
            }
        }
    }
    
    # Filtra coefficienti validi
    valid_idx <- is.finite(coeffs) & !is.na(coeffs) & coeffs != 0
    if (!any(valid_idx)) {
        return(rep(NA, length(corr)))
    }
    
    coeffs <- coeffs[valid_idx]
    j_vec_valid <- j_vec[valid_idx]
    gamma_terms_valid <- gamma_j_plus_1[valid_idx]
    coeffs_sq <- coeffs^2  # Pre-calcolo
    
    # Funzione interna vettorizzata e ottimizzata
    corrsas_inner <- function(rho) {
        if (!is.finite(rho) || is.na(rho)) return(NA)
        if (rho == 0) return(0)
        
        # Calcolo diretto senza tryCatch per velocità
        rho_powers <- rho^j_vec_valid
        if (any(is.infinite(rho_powers))) return(NA)
        
        numerator <- sum(coeffs_sq * rho_powers / gamma_terms_valid)
        result <- numerator / vv
        
        if (is.finite(result)) result else NA
    }
    
    # Applicazione vettorizzata ottimizzata
    if (length(corr) == 1) {
        return(corrsas_inner(corr))
    } else {
        # Per vettori lunghi, usa vapply che è più veloce
        return(vapply(corr, corrsas_inner, numeric(1)))
    }
}




######################################################################################################
############## functions for gaussiam copula #########################################################
######################################################################################################
# marginal variances
variance_disp <- function(model_type, params) {
  switch(as.character(model_type),
    "1"  = as.numeric(params["sill"]),  # Gaussian
    "12" = as.numeric(params["sill"]) * (1/as.numeric(params["df"])) / (1/as.numeric(params["df"]) - 2),  # StudentT
    "22" = exp(2 * as.numeric(params["mean"])) * (exp(as.numeric(params["sill"])) - 1),  # LogGaussian
    "30" = exp(as.numeric(params["mean"])),  # Poisson
    "21" = 2 * exp(2 * as.numeric(params["mean"])) / as.numeric(params["shape"]),  # Gamma
    "26" = exp(2 * as.numeric(params["mean"])) * (gamma(1 + 2 / as.numeric(params["shape"])) /(gamma(1 + 1 / as.numeric(params["shape"]))^2) - 1),  # Weibull
    "11" = {n <- as.numeric(params["n"]);pp <- pnorm(as.numeric(params["mean"])); n * pp * (1 - pp)},  # Binomial
    "16" = {n <- as.numeric(params["n"]); pp<- pnorm(as.numeric(params["mean"]));n  * (1 - pp )/ pp^2},  # BinomialNeg
    "18" = {df <- round(1 / as.numeric(params["df"]));skew <- as.numeric(params["skew"]); mu <- as.numeric(params["mean"])
            mu * ((df / (df - 2)) * (1 + skew^2) - df * skew^2 * gamma(0.5 * (df - 1)) / (pi * gamma(0.5 * df)))
           },  # SkewStudentT
    "10" = as.numeric(params["sill"]) + as.numeric(params["skew"])^2 * (1 - 2 / pi),  # SkewGaussian
    "34" = { h <- as.numeric(params["h"]); as.numeric(params["sill"]) * (1 - 2 * h)^(-1.5)},  # Tukeyh
    stop(paste("Model not supported:", model_type))
  )
}

# quantiles
quantile_disp <- function(p, model_type, params) {
  switch(as.character(model_type),
    "1" = qnorm(p, mean = as.numeric(params["mean"]), sd = sqrt(as.numeric(params["sill"]))), # gaussian
    "12" = as.numeric(params["mean"]) + sqrt(as.numeric(params["sill"])) * qt(p, df = 1 / as.numeric(params["df"])), # StudentT
    "22" = qlnorm(p, meanlog = as.numeric(params["mean"]) - as.numeric(params["sill"]) / 2, sdlog = sqrt(as.numeric(params["sill"]))), # loggauss
    "30" = qpois(p, lambda = exp(as.numeric(params["mean"]))),  #poisson
    "21" = exp(as.numeric(params["mean"])) * qgamma(p, shape = as.numeric(params["shape"]) / 2, rate = as.numeric(params["shape"]) / 2), # gamma
    "26" = exp(as.numeric(params["mean"])) * qweibull(p, shape = as.numeric(params["shape"]), scale = 1 / gamma(1 + 1 / as.numeric(params["shape"]))), #weibull
    "11" = { n <- as.numeric(params["n"]); qbinom(p, size = n, prob = pnorm(as.numeric(params["mean"]))) },   #   binom
    "16" = { n <- as.numeric(params["n"]); qnbinom(p, size = n, prob = pnorm(as.numeric(params["mean"])))  }, # neg  binom
    "18" = as.numeric(params["mean"]) + sqrt(as.numeric(params["sill"])) * sn::qst(p, xi = 0, omega = 1, alpha = as.numeric(params["skew"]), nu = as.numeric(round(1 / params["df"]))), ## SkewStudentT
    "10" = as.numeric(params["mean"]) + sqrt(as.numeric(params["sill"])) * sn::qsn(p, xi = 0, omega = sqrt((as.numeric(params["skew"])^2 + as.numeric(params["sill"])) / as.numeric(params["sill"])), 
                alpha = as.numeric(params["skew"]) / sqrt(as.numeric(params["sill"]))), #skewgaussian
    "34" = as.numeric(params["mean"]) + sqrt(as.numeric(params["sill"])) * p * exp(0.5 * as.numeric(params["tail"]) * p^2),# Tukeyh
    stop(paste("Model not supported:", model_type))
  )
}

#####################################################################
# coefficients a_k in copula gaussian covariance 
####################################################################

  hermite_single <- function(k, t) {
    if (k == 0) return(rep(1, length(t)))
    if (k == 1) return(t)
    H_prev <- rep(1, length(t))
    H_curr <- t
    if (k > 1) {
      for (i in 2:k) {
        H_next <- t * H_curr - (i - 1) * H_prev
        H_prev <- H_curr
        H_curr <- H_next
      }
    }
    return(H_curr)
  }
compute_ak_vectorized <- function(M, model_type, params) {
  ak <- numeric(M)
  failed_integrations <- c()


  for (k in 1:M) {
    integrand <- function(t) {
      p <- pnorm(t)
      q <- quantile_disp(p, model_type, params)
      q[!is.finite(q)] <- 0
      result <- q * hermite_single(k, t) * dnorm(t)
      # Controllo aggiuntivo per valori non finiti
      result[!is.finite(result)] <- 0
      return(result)
    }
    # Strategia di integrazione progressiva
    result <- NA
    # Tentativo 1: Integrazione standard con più suddivisioni
    result <- tryCatch({
      integrate(integrand, lower = -12, upper = 12, 
                rel.tol = 1e-6, abs.tol = 1e-10, 
                subdivisions = 2000L)$value
    }, error = function(e) NA, warning = function(w) NA)

    # Tentativo 2: Limiti più ristretti se il primo fallisce
    if (is.na(result)) {
      result <- tryCatch({
        integrate(integrand, lower = -8, upper = 8, 
                  rel.tol = 1e-5, abs.tol = 1e-8, 
                  subdivisions = 1000L)$value
      }, error = function(e) NA, warning = function(w) NA)
    }
    # Tentativo 3: Integrazione adattiva con quadratura di Gauss-Hermite
  #  if (is.na(result)) {
  #    result <- tryCatch({
  #      integrate_gauss_hermite(integrand, k)
  #    }, error = function(e) NA)
  #  }
    # Tentativo 4: Approssimazione per k elevati
    if (is.na(result)) {
      if (k > 20) {
        # Per k elevati, i coefficienti tendono rapidamente a zero
        result <- 0
        #message(paste("Coefficient a_", k, " approximated to 0 (high k)", sep = ""))
      } else {
        # Integrazione molto conservativa
        result <- tryCatch({
          integrate(integrand, lower = -6, upper = 6, 
                    rel.tol = 1e-4, abs.tol = 1e-6,
                    subdivisions = 500L, stop.on.error = FALSE)$value
        }, error = function(e) {
          failed_integrations <<- c(failed_integrations, k)
          0
        })
      }
    }
    
    ak[k] <- ifelse(is.finite(result), result, 0)
  }
  
  # Report sui fallimenti
  if (length(failed_integrations) > 0) {
    message(paste("Integration failed for k =", paste(failed_integrations, collapse = ", "), 
                  "- coefficients set to 0"))
  }
  
  return(ak)
}
# main fnction


# Funzione principale migliorata
gaussian_copula_cov_fast <- function(rho, model_type, nuisance, M=30, cache_ak = TRUE, 
                                     auto_reduce_M = TRUE, verbose = FALSE) {
  

  # Cache management
  cache_key <- paste(model_type, paste(nuisance, collapse = "_"), M, sep = "_")
  
  ak_cache <- .GeoModels_env$ak_cache
  if (cache_ak && cache_key %in% names(ak_cache)) {
    ak_coeffs <- ak_cache[[cache_key]]
    if (verbose) message("Coefficients a_k loaded from cache")
  } else {
    if (verbose) message(paste("Computing coefficients a_k for M =", M))
    
    # Calcolo coefficienti con gestione automatica di M
    ak_coeffs <- compute_ak_vectorized(M, model_type, nuisance)
    
    # Riduzione automatica di M se molti coefficienti sono problematici
    if (auto_reduce_M) {
      # Trova l'ultimo coefficiente significativo
      significant_coeffs <- which(abs(ak_coeffs) > 1e-10)
      if (length(significant_coeffs) > 0) {
        effective_M <- max(significant_coeffs)
        if (effective_M < M * 0.7) {  
          M_new <- min(effective_M + 5, M)  
          ak_coeffs <- ak_coeffs[1:M_new]
          M <- M_new
          #message(paste("M automatically reduced to", M, "to improve stability"))
        }
      }
    }
    
    if (cache_ak) {
      ak_cache[[cache_key]] <- ak_coeffs
      .GeoModels_env$ak_cache <- ak_cache
    }
  }
  
  # Controllo convergenza della serie
  if (length(ak_coeffs) >= 10) {
    last_coeffs <- abs(ak_coeffs[(length(ak_coeffs)-4):length(ak_coeffs)])
    if (all(last_coeffs < 1e-8)) {
      if (verbose) message("Series converges well - last coefficients are small")
    } #else {
     # warning("Series may not converge well - consider reducing M")
    #}
  }
  
  # Calcolo covarianze
  M_eff <- length(ak_coeffs)
  k_vec <- 1:M_eff
  fact_k <- factorial(k_vec)
  outer_term <- (ak_coeffs^2) / fact_k
  
  covariances <- sapply(rho, function(r) {
    if (abs(r) < 1e-12) return(0)  # Correlazione quasi zero → Covarianza zero
    
    # Calcolo della serie con controllo overflow
    terms <- outer_term * r^k_vec
    
    # Controlla per overflow/underflow
    valid_terms <- is.finite(terms) & abs(terms) > .Machine$double.eps
    
    result <- sum(terms[valid_terms])
    
    if (!is.finite(result)) {
      #warning(paste("Non-finite result for rho =", r))
      return(0)
    }
    
    return(result)
  })
  
  if (verbose) {
    message(paste("Computed with", M_eff, "coefficients"))
    message(paste("Covariance range:", round(min(covariances), 6), "-", round(max(covariances), 6)))
  }
  
  return(covariances)
}



###############################################
Cmatrix_copula <- function(bivariate, coordx, coordy,coordz, coordt,corrmodel, dime, n, ns, NS, nuisance, numpairs,
                           numpairstot, model, paramcorr, setup, radius, spacetime, spacetime_dyn,type,copula,ML,other_nuis)
    {
###################################################################################
############### computing correlation #############################################
###################################################################################

  if(type=="Standard") {
    if(spacetime) 
       cr=dotCall64::.C64('CorrelationMat_st_dyn2',SIGNATURE = c(rep("double",5),"integer","double","double","double","integer","integer"),
            corr=dotCall64::vector_dc("double",numpairstot),coordx,coordy,coordz,coordt,
            corrmodel,nuisance,paramcorr,radius,ns,NS,
            INTENT = c("w",rep("r",10)),
            PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
     if(bivariate) 
       cr=dotCall64::.C64('CorrelationMat_biv_dyn2',SIGNATURE = c(rep("double",5),"integer","double","double","double","integer","integer"),
            corr=dotCall64::vector_dc("double",numpairstot),coordx,coordy,coordz,coordt,
            corrmodel,nuisance,paramcorr,radius,ns,NS,
            INTENT = c("w",rep("r",10)),
            PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
   if(!bivariate&&!spacetime)  
        cr=dotCall64::.C64('CorrelationMat2',SIGNATURE = c(rep("double",5),"integer","double","double","double","integer","integer"),
            corr=dotCall64::vector_dc("double",numpairstot),coordx,coordy,coordz,coordt,
            corrmodel,nuisance,paramcorr,radius,ns,NS,
            INTENT = c("w",rep("r",10)),
            PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
    
}
###############################################################
if(type=="Tapering")  {
        fname <- 'CorrelationMat_tap'
        if(spacetime) fname <- 'CorrelationMat_st_tap'
       if(bivariate) fname <- 'CorrelationMat_biv_tap'
#############
cr=dotCall64::.C64(fname,SIGNATURE = c("double","double","double","double","double","integer","double","double","double","integer","integer"),
     corr=dotCall64::vector_dc("double",numpairs), coordx,coordy,coordz,coordt,corrmodel,nuisance, paramcorr,radius,ns,NS,
 INTENT = c("w","r","r","r","r","r","r","r", "r", "r", "r"),
            PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
#############
     ## deleting correlation equual  to 1 because there are problems  with hipergeometric function
        sel=(abs(cr$corr-1)<.Machine$double.eps);cr$corr[sel]=0
      }
  
corr=cr$corr*(1-as.numeric(nuisance['nugget'])) 
if(model==11||model==16)  nuisance["n"]=n

cova=gaussian_copula_cov_fast(corr, model,nuisance  )
vv=variance_disp(model,nuisance)

  if(!bivariate)
{
 if(type=="Standard"){
     # Builds the covariance matrix:
        varcov <-  diag(dime)
        varcov[lower.tri(varcov)] <- cova
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- cova
        diag(varcov)=vv
    }
    if(type=="Tapering")  {
          vcov <- cova;
          varcov <- new("spam",entries=vcov,colindices=setup$ja,
                             rowpointers=setup$ia,dimension=as.integer(rep(dime,2)))
          diag(varcov)=vv
        }
}
else  { }

return(varcov)
} 
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
