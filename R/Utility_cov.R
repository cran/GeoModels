######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################


### Internal cache environment for covariance utilities
if (!exists(".GeoModels_env", inherits = FALSE)) {
  .GeoModels_env <- new.env(parent = emptyenv())
}

.GeoModels_get_cache <- function(name) {
  if (is.null(.GeoModels_env[[name]]) ||
      !is.environment(.GeoModels_env[[name]])) {
    .GeoModels_env[[name]] <- new.env(parent = emptyenv())
  }
  .GeoModels_env[[name]]
}


### Observed-scale correlation for the constructive Weibull random field.
###
### For shape = kappa, GeoModels uses the correlation from Eq. (6) of the
### Weibull-process paper,
###
###   C(kappa) * { 2F1(-1/kappa,-1/kappa;1;rho^2) - 1 }.
###
### The exact endpoint rho^2 = 1 is one.  Evaluating hypergeo::hypergeo()
### at that boundary is unnecessary and, for some shapes (notably kappa=2),
### can make the numerical backend attempt a divergent integral.  Handle the
### endpoints analytically and evaluate the hypergeometric function only for
### 0 < rho^2 < 1.
.GeoWeibullCorr <- function(rho, shape) {
  shape <- suppressWarnings(as.numeric(shape))
  if (length(shape) != 1L || !is.finite(shape) || shape <= 0)
    stop("Weibull shape must be a finite positive scalar.", call. = FALSE)

  rho <- as.numeric(rho)
  if (!length(rho)) return(numeric())
  if (any(!is.finite(rho)))
    stop("Weibull latent correlations must be finite.", call. = FALSE)

  ## Correlations supplied by the covariance core should lie in [-1,1].
  ## Clamp only round-off excursions before squaring.
  rho <- pmax(-1, pmin(1, rho))
  z <- rho * rho

  out <- numeric(length(z))
  at_zero <- (z == 0)
  at_one <- (z == 1)
  mid <- !(at_zero | at_one)

  out[at_zero] <- 0
  out[at_one] <- 1

  if (any(mid)) {
    ## C(kappa) = Gamma(1+1/kappa)^2 /
    ##            {Gamma(1+2/kappa)-Gamma(1+1/kappa)^2}.
    ## Compute it through lgamma/expm1 to avoid avoidable overflow.
    delta <- lgamma(1 + 2 / shape) - 2 * lgamma(1 + 1 / shape)
    denom_ratio <- expm1(delta)
    if (!is.finite(denom_ratio) || denom_ratio <= 0)
      stop("Invalid Weibull correlation normalizing constant.", call. = FALSE)

    factor <- 1 / denom_ratio
    hz <- Re(hypergeo::hypergeo(-1 / shape, -1 / shape, 1, z[mid]))
    out[mid] <- factor * (hz - 1)
  }

  ## Protect covariance matrices against tiny numerical overshoots while
  ## retaining the exact theoretical endpoints.
  out[abs(out) < 32 * .Machine$double.eps] <- 0
  out[abs(out - 1) < 64 * .Machine$double.eps] <- 1
  pmax(0, pmin(1, out))
}

### decomposition of a square  matrix
MatDecomp <- function(mtx, method) {
  # Controlli preliminari veloci
  if (!is.matrix(mtx) || anyNA(mtx) || any(!is.finite(mtx))) return(FALSE)
  
  if (method == "cholesky") {
    if (!isSymmetric(mtx, tol = 1e-8)) return(FALSE)
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
    return(t(t(mat.decomp$v) * d_sqrt))
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
    mse_val <- if (mse) crossprod(vec) else NULL
  } else {
    cc <- if (spam::is.spam(covmatrix$covmatrix)) covmatrix$covmatrix else spam::as.spam(covmatrix$covmatrix)
    U <- tryCatch({
      spam::chol.spam(cc)
    }, error = function(e) {
      stop("Covariance matrix is not positive definite")
    })
    vec <- spam::forwardsolve(U, CC)
    Invc <- spam::backsolve(U, vec)
    #mse_val <- if (mse) as.matrix(spam::crossprod.spam(vec)) else NULL
    mse_val <- if (mse) as.matrix(crossprod(vec)) else NULL
    
  }
  list(a = Invc, b = mse_val)
}
  
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
.GeoTukeyh2StandardMoments <- function(tail1, tail2) {
  hr <- suppressWarnings(as.numeric(tail1))[1L]
  hl <- suppressWarnings(as.numeric(tail2))[1L]
  if (!is.finite(hr) || !is.finite(hl) || hr < 0 || hl < 0 ||
      hr >= 0.5 || hl >= 0.5) {
    stop("Tukeyh2 tail1 and tail2 must lie in [0, 0.5).", call. = FALSE)
  }
  mm <- (hr - hl) / (sqrt(2 * pi) * (1 - hl) * (1 - hr))
  vv <- 0.5 * (1 - 2 * hl)^(-1.5) +
        0.5 * (1 - 2 * hr)^(-1.5) - mm^2
  list(mean = mm, variance = vv)
}

.GeoTukeyh2Correlation <- function(rho, tail1, tail2) {
  mom <- .GeoTukeyh2StandardMoments(tail1, tail2)
  hr <- as.numeric(tail1)[1L]
  hl <- as.numeric(tail2)[1L]
  rho <- as.numeric(rho)
  if (!length(rho)) return(numeric(0))
  out <- rep(NA_real_, length(rho))
  good <- is.finite(rho) & !is.na(rho)
  if (!any(good)) return(out)

  rr0 <- pmax(-1, pmin(1, rho[good]))
  ans <- numeric(length(rr0))
  zero <- abs(rr0) < 1e-15
  one <- abs(rr0 - 1) < 1e-15
  work <- !(zero | one)
  ans[zero] <- 0
  ans[one] <- 1

  if (any(work)) {
    eps <- 1e-12
    rr <- pmax(-1 + 1e-12, pmin(1 - 1e-12, rr0[work]))
    x1 <- 1 - (1 - rr^2) * hr
    y1 <- 1 - (1 - rr^2) * hl
    x2 <- pmax((1 - hr)^2 - (rr * hr)^2, eps)
    y2 <- pmax((1 - hl)^2 - (rr * hl)^2, eps)
    gg <- pmax(1 - hl - hr + (1 - rr^2) * hl * hr, eps)

    ax1 <- pmin(1, pmax(-1, rr / x1))
    ay1 <- pmin(1, pmax(-1, rr / y1))
    x1y1 <- pmax(x1 * y1, eps)
    axy <- pmin(1, pmax(-1, rr / sqrt(x1y1)))

    h1 <- sqrt(pmax(1 - rr^2 / x1^2, 0)) + (rr / x1) * asin(ax1)
    h2 <- sqrt(pmax(1 - rr^2 / y1^2, 0)) + (rr / y1) * asin(ay1)
    h3 <- sqrt(pmax(1 - rr^2 / x1y1, 0)) +
          (rr / sqrt(x1y1)) * asin(axy)

    p1 <- x1 * h1 / (2 * pi * x2^(3/2)) + rr / (4 * x2^(3/2))
    p2 <- y1 * h2 / (2 * pi * y2^(3/2)) + rr / (4 * y2^(3/2))
    p3 <- -sqrt(x1y1) * h3 / (2 * pi * gg^(3/2)) + rr / (4 * gg^(3/2))
    ans[work] <- (p1 + p2 + 2 * p3 - mom$mean^2) / mom$variance
  }
  ans <- pmax(-1, pmin(1, ans))
  out[good] <- ans
  out
}

.GeoSASStandardMoments <- function(skew, tail) {
  e <- suppressWarnings(as.numeric(skew))[1L]
  d <- suppressWarnings(as.numeric(tail))[1L]
  if (!is.finite(e) || !is.finite(d) || d <= 0)
    stop("SinhAsinh requires finite skew and a positive finite tail.", call. = FALSE)
  mm <- sinh(e / d) * exp(0.25) *
    (besselK(0.25, (d + 1) / (2 * d)) +
       besselK(0.25, abs(1 - d) / (2 * d))) / sqrt(8 * pi)
  vv <- cosh(2 * e / d) * exp(0.25) *
    (besselK(0.25, (d + 2) / (2 * d)) +
       besselK(0.25, abs(2 - d) / (2 * d))) / sqrt(32 * pi) - 0.5 - mm^2
  if (!is.finite(vv) || vv <= 0)
    stop("SinhAsinh marginal variance is not finite and positive.", call. = FALSE)
  list(mean = mm, variance = vv)
}

.GeoPlainTransformedMarginalMoments <- function(model_type, params, eta = 0) {
  model_type <- as.integer(model_type)
  eta <- as.numeric(eta)
  if (!length(eta) || any(!is.finite(eta)))
    stop("Marginal locations must be finite.", call. = FALSE)
  sill <- suppressWarnings(as.numeric(params["sill"]))[1L]
  if (!is.finite(sill) || sill <= 0)
    stop("sill must be finite and positive.", call. = FALSE)

  if (model_type == 34L) {
    h <- suppressWarnings(as.numeric(params["tail"]))[1L]
    if (!is.finite(h) || h < 0 || h >= 0.5)
      stop("Tukeyh tail must lie in [0, 0.5).", call. = FALSE)
    return(list(mean = eta,
                variance = rep(sill * (1 - 2 * h)^(-1.5), length(eta))))
  }
  if (model_type == 40L) {
    mom <- .GeoTukeyh2StandardMoments(params["tail1"], params["tail2"])
    return(list(mean = eta + sqrt(sill) * mom$mean,
                variance = rep(sill * mom$variance, length(eta))))
  }
  if (model_type == 20L) {
    mom <- .GeoSASStandardMoments(params["skew"], params["tail"])
    return(list(mean = eta + sqrt(sill) * mom$mean,
                variance = rep(sill * mom$variance, length(eta))))
  }
  stop("Unsupported plain transformed model id ", model_type, ".", call. = FALSE)
}

.GeoPlainTransformedCrossCov <- function(rho, model_type, params) {
  model_type <- as.integer(model_type)
  rho <- as.numeric(rho)
  sill <- suppressWarnings(as.numeric(params["sill"]))[1L]
  if (!is.finite(sill) || sill <= 0)
    stop("sill must be finite and positive.", call. = FALSE)

  if (model_type == 34L) {
    h <- suppressWarnings(as.numeric(params["tail"]))[1L]
    if (!is.finite(h) || h < 0 || h >= 0.5)
      stop("Tukeyh tail must lie in [0, 0.5).", call. = FALSE)
    if (h == 0) return(sill * rho)
    v0 <- (1 - 2 * h)^(-1.5)
    cc <- rho * (1 - 2 * h)^(1.5) /
      pmax((1 - h)^2 - (h * rho)^2, .Machine$double.eps)^(1.5)
    return(sill * v0 * cc)
  }
  if (model_type == 40L) {
    mom <- .GeoTukeyh2StandardMoments(params["tail1"], params["tail2"])
    return(sill * mom$variance *
             .GeoTukeyh2Correlation(rho, params["tail1"], params["tail2"]))
  }
  if (model_type == 20L) {
    ## `params` may be either a named numeric vector or a list.  Using
    ## params["skew"] / params["tail"] preserves list type in the latter
    ## case, while corrsas() deliberately requires numeric scalar inputs.
    ## Normalize once here so plain and Gaussian-copula SAS kriging share
    ## exactly the same covariance implementation.
    skew <- suppressWarnings(as.numeric(params[["skew"]]))[1L]
    tail <- suppressWarnings(as.numeric(params[["tail"]]))[1L]
    mom <- .GeoSASStandardMoments(skew, tail)
    return(sill * mom$variance * corrsas(rho, skew, tail))
  }
  stop("Unsupported plain transformed model id ", model_type, ".", call. = FALSE)
}

corrsas <- function(corr, skew, tail, max_coeff = NULL) {
    d <- tail
    e <- skew
    # ========================================================================
    if (!is.numeric(skew) || !is.numeric(tail)) {
        stop("skew and tail must be numeric")
    }
    if (length(skew) != 1 || length(tail) != 1) {
        stop("skew and tail must be scalars")
    }
    if (tail <= 0) {
        stop("tail must be positive")
    }
    if (!is.numeric(corr)) {
        stop("corr must be numeric")
    }
    if (length(corr) == 0) {
        return(numeric(0))
    }
    ## Exact Gaussian identity and independence do not require the
    ## Hermite/integration machinery below.
    if (abs(skew) < .Machine$double.eps && abs(tail - 1) < .Machine$double.eps)
        return(as.numeric(corr))
    if (all(corr == 0)) return(rep(0, length(corr)))
    
    # Determinazione automatica max_coeff
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
    
    # ========================================================================
    # ========================================================================
    compute_bessel_terms <- function(d) {
        tryCatch({
            k1 <- besselK(0.25, (d + 1)/(2*d))
            k2 <- besselK(0.25, abs(1 - d)/(2*d))
            k3 <- besselK(0.25, (d + 2)/(2*d))
            k4 <- besselK(0.25, abs(2 - d)/(2*d))
            if (any(!is.finite(c(k1, k2, k3, k4)))) {return(NULL)}
            list(k1 = k1, k2 = k2, k3 = k3, k4 = k4)
        }, error = function(e) {
            warning("Bessel function computation failed: ", e$message)
            NULL
        })
    }

    bessel_vals <- compute_bessel_terms(d)
    if (is.null(bessel_vals)) {
        return(rep(NA_real_, length(corr)))
    }
    
    #
    besselK_1 <- bessel_vals$k1
    besselK_2 <- bessel_vals$k2
    besselK_3 <- bessel_vals$k3
    besselK_4 <- bessel_vals$k4
    
    # 
    if (abs(e/d) < 0.01) {
        # Per valori molto piccoli, usa serie di Taylor
        x <- e/d
        sinh_e_d <- x * (1 + x^2/6)
        cosh_2e_d <- 1 + 2*x^2 * (1 + x^2/3)
    } else {
        sinh_e_d <- sinh(e/d)
        cosh_2e_d <- cosh(2*e/d)
    }
    
    mm <- sinh_e_d * exp_0.25 * (besselK_1 + besselK_2) / sqrt_8pi
    vv <- cosh_2e_d * exp_0.25 * (besselK_3 + besselK_4) / sqrt_32pi - 0.5 - mm^2
    
    if (!is.finite(vv) || abs(vv) < .Machine$double.eps) {
        return(rep(NA_real_, length(corr)))
    }
    
    # Pre-calcolo di j_vec e gamma terms
    j_vec <- 1:max_coeff
    gamma_j_plus_1 <- gamma(j_vec + 1)
    # ========================================================================
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
    # ========================================================================
    # ========================================================================
    II_cache <- .GeoModels_get_cache("II_global_cache")
    MAX_CACHE_SIZE <- 5000
    # 
    cache_size <- length(names(II_cache))
    if (cache_size >= MAX_CACHE_SIZE) {
        # Rimuovi metà degli elementi (strategia semplice)
        keys_to_remove <- head(sort(ls(II_cache, all.names = TRUE),
                                    method = "radix"), floor(MAX_CACHE_SIZE / 2))
        rm(list = keys_to_remove, envir = II_cache)
    }
    
    # =======================================================================
    # ========================================================================
    II <- function(alpha, kappa, j, r) {
        key <- paste(sprintf("%.17g", alpha), sprintf("%.17g", kappa),
                     j, r, sep = "_")
        if (exists(key, envir = II_cache, inherits = FALSE)) {
            return(get(key, envir = II_cache, inherits = FALSE))
        }
        val <- tryCatch({
            # Limiti finiti: l'integranda decade come exp(-z²/2)
            # Oltre ±6 il contributo è trascurabile (< 1e-9)
            integrate(integrand, 
                     lower = -6, 
                     upper = 6, 
                     alpha = alpha, 
                     kappa = kappa, 
                     j = j, 
                     r = r,
                     rel.tol = 1e-7,
                     abs.tol = 1e-9,
                     subdivisions = 150,
                     stop.on.error = FALSE)$value
        }, error = function(e) {
            warning("Integration failed for key ", key, ": ", e$message, call. = FALSE)
            0
        })
        # 
        if (is.finite(val)) {
            assign(key, val, envir = II_cache)
        } else {
            val <- 0
        }
        
        val
    }
    kahan_sum <- function(x) {
        if (length(x) == 0) return(0)
        s <- 0
        c <- 0
        for (xi in x) {
            y <- xi - c
            t <- s + y
            c <- (t - s) - y
            s <- t
        }
        s
    }
    # ========================================================================
    # ========================================================================
    coeffs <- numeric(max_coeff)
    for (j in j_vec) {
        max_r <- floor(j/2)
        if (max_r < 0) {
            coeffs[j] <- 0
            next
        }
        
        rr <- 0:max_r
        II_vals <- sapply(rr, function(r_val) II(e, d, j, r_val))
        gamma_r_plus_1 <- gamma(rr + 1)
        gamma_j_minus_2r_plus_1 <- gamma(j - 2*rr + 1)
        valid_gamma <- is.finite(gamma_r_plus_1) & is.finite(gamma_j_minus_2r_plus_1) & 
                      gamma_r_plus_1 > 0 & gamma_j_minus_2r_plus_1 > 0
        if (any(valid_gamma)) {
            terms <- II_vals[valid_gamma] * (-1)^rr[valid_gamma] / 
                    (2^(rr[valid_gamma] + 1) * gamma_r_plus_1[valid_gamma] * 
                     gamma_j_minus_2r_plus_1[valid_gamma])
            
            if (is.finite(gamma_j_plus_1[j]) && gamma_j_plus_1[j] > 0) {
                coeffs[j] <- gamma_j_plus_1[j] * kahan_sum(terms) / sqrt_2pi
            }
        }
    }
    valid_idx <- is.finite(coeffs) & !is.na(coeffs) & coeffs != 0
    if (!any(valid_idx)) {
        return(rep(NA_real_, length(corr)))
    }
    coeffs <- coeffs[valid_idx]
    j_vec_valid <- j_vec[valid_idx]
    gamma_terms_valid <- gamma_j_plus_1[valid_idx]
    coeffs_sq <- coeffs^2
    
    # ========================================================================
    # ========================================================================
    corrsas_inner <- function(rho) {
        if (!is.finite(rho) || is.na(rho)) return(NA_real_)
        if (rho == 0) return(0)
        if (rho == 1) return(1)
        
        # Gestione overflow con log-scale per correlazioni alte
        max_j <- max(j_vec_valid)
        if (abs(rho) > 0.99 && max_j > 50) {
            log_rho <- log(abs(rho))
            if (max_j * abs(log_rho) > 700) return(NA_real_)
            
            # Calcolo in log-scale
            log_terms <- 2*log(abs(coeffs)) + j_vec_valid * log_rho - log(gamma_terms_valid)
            max_log <- max(log_terms[is.finite(log_terms)])
            if (!is.finite(max_log)) return(NA_real_)
            
            terms_scaled <- exp(log_terms - max_log)
            numerator <- exp(max_log) * kahan_sum(terms_scaled)
        } else {
            rho_powers <- rho^j_vec_valid
            if (any(is.infinite(rho_powers))) return(NA_real_) 
            numerator <- kahan_sum(coeffs_sq * rho_powers / gamma_terms_valid)
        }
        result <- numerator / vv
        if (is.finite(result)) result else NA_real_
    }
    
    # ========================================================================
    # ========================================================================
    if (length(corr) == 1) {
        return(corrsas_inner(corr))
    } else {
        return(vapply(corr, corrsas_inner, numeric(1)))
    }
}



.GeoCovCanonicalModel <- function(model) {
  if(!is.character(model) || length(model) != 1L || is.na(model)) return(model)
  model <- gsub("[[:blank:]]", "", model)
  aliases <- c(
    Gauss = "Gaussian", SkewGauss = "SkewGaussian", LogGauss = "LogGaussian",
    TwoPieceGauss = "TwoPieceGaussian", TwopieceGauss = "TwoPieceGaussian",
    TwopieceT = "TwoPieceStudentT", TwopieceTukeyh = "TwoPieceTukeyh",
    Loglogistic = "LogLogistic",
    poisson = "Poisson", poissongamma = "PoissonGamma",
    Gaussian_misp_StudentT = "StudentT",
    Gaussian_misp_Poisson = "Poisson",
    Gaussian_misp_SkewStudentT = "SkewStudentT",
    Gaussian_misp_Tukeygh = "Tukeygh",
    Gaussian_misp_PoissonZIP = "PoissonZIP",
    Gaussian_misp_PoissonZIP1 = "PoissonZIP1",
    Gaussian_misp_PoissonGamma = "PoissonGamma",
    Gaussian_misp_Binomial = "Binomial",
    Gaussian_misp_BinomialNeg = "BinomialNeg",
    Binary_misp_BinomialNeg = "BinomialNeg"
  )
  if(model %in% names(aliases)) unname(aliases[[model]]) else model
}

######################################################################################################
############## functions for gaussiam copula #########################################################
######################################################################################################
# marginal variances
variance_disp <- function(model_type, params) {
  model_type <- as.integer(model_type)
  eta <- suppressWarnings(as.numeric(params["mean"]))
  if (!length(eta) || !is.finite(eta[1L])) eta <- 0
  eta <- eta[1L]

  switch(as.character(model_type),
    "1"  = as.numeric(params["sill"]),  # Gaussian
    "12" = {
      df_t <- 1 / as.numeric(params["df"])
      as.numeric(params["sill"]) * df_t / (df_t - 2)
    },  # StudentT
    "22" = exp(2 * eta) * (exp(as.numeric(params["sill"])) - 1),  # LogGaussian
    "30" = exp(eta),  # Poisson
    "21" = 2 * exp(2 * eta) / as.numeric(params["shape"]),  # Gamma
    "26" = exp(2 * eta) * (gamma(1 + 2 / as.numeric(params["shape"])) /
      (gamma(1 + 1 / as.numeric(params["shape"]))^2) - 1),  # Weibull
    "11" = {
      n <- as.numeric(params["n"])
      pp <- pnorm(eta)
      n * pp * (1 - pp)
    },  # Binomial
    "16" = {
      n <- as.numeric(params["n"])
      pp <- pnorm(eta)
      n * (1 - pp) / pp^2
    },  # BinomialNeg
    "28" = {
      a <- as.numeric(params["shape1"]); b <- as.numeric(params["shape2"])
      rr <- as.numeric(params["max"] - params["min"])
      rr^2 * a * b / ((a + b)^2 * (a + b + 1))
    },  # Beta
    "50" = {
      m <- plogis(eta); phi <- as.numeric(params["shape"])
      rr <- as.numeric(params["max"] - params["min"])
      rr^2 * m * (1 - m) / (phi + 1)
    },  # Beta2
    "33" = {
      a <- as.numeric(params["shape1"]); b <- as.numeric(params["shape2"])
      rr <- as.numeric(params["max"] - params["min"])
      m1 <- b * beta(1 + 1 / a, b)
      m2 <- b * beta(1 + 2 / a, b)
      rr^2 * (m2 - m1^2)
    },  # Kumaraswamy
    "42" = {
      a <- as.numeric(params["shape"]); med <- plogis(eta)
      b <- log(0.5) / log1p(-med^a)
      rr <- as.numeric(params["max"] - params["min"])
      m1 <- b * beta(1 + 1 / a, b)
      m2 <- b * beta(1 + 2 / a, b)
      rr^2 * (m2 - m1^2)
    },  # Kumaraswamy2
    "25" = pi^2 * as.numeric(params["sill"]) / 3,  # Logistic
    "59" = {
      sk <- as.numeric(params["skew"]); ss <- as.numeric(params["sill"])
      ss * (1 - 2 * sk + 2 * sk^2) / (sk^2 * (1 - sk)^2)
    },  # SkewLaplace
    "18" = {
      df <- 1 / as.numeric(params["df"])
      skew <- as.numeric(params["skew"])
      mu <- as.numeric(params["mean"])
      mu * ((df / (df - 2)) * (1 + skew^2) -
              df * skew^2 * gamma(0.5 * (df - 1)) / (pi * gamma(0.5 * df)))
    },  # legacy SkewStudentT expression retained for compatibility
    "10" = as.numeric(params["sill"]) + as.numeric(params["skew"])^2 * (1 - 2 / pi),
    "34" = {
      h <- suppressWarnings(as.numeric(params["tail"]))
      if (!length(h) || is.na(h)) h <- suppressWarnings(as.numeric(params["h"]))
      if (!is.finite(h) || h < 0 || h >= 0.5)
        stop("Tukeyh tail must lie in [0, 0.5).", call. = FALSE)
      as.numeric(params["sill"]) * (1 - 2 * h)^(-1.5)
    },
    "40" = {
      mom <- .GeoTukeyh2StandardMoments(params["tail1"], params["tail2"])
      as.numeric(params["sill"]) * mom$variance
    },
    "20" = {
      mom <- .GeoSASStandardMoments(params["skew"], params["tail"])
      as.numeric(params["sill"]) * mom$variance
    },
    stop(paste("Model not supported:", model_type))
  )
}

# Quantile of the observed margin evaluated at a standard-uniform probability.
# The location parameter is taken from params["mean"].  This function is used
# by the Hermite representation of Gaussian-copula covariances.
quantile_disp <- function(p, model_type, params) {
  model_type <- as.integer(model_type)
  eta <- suppressWarnings(as.numeric(params["mean"]))
  if (!length(eta) || !is.finite(eta[1L])) eta <- 0
  eta <- eta[1L]
  p <- pmin(pmax(as.numeric(p), .Machine$double.eps), 1 - .Machine$double.eps)

  switch(as.character(model_type),
    "1" = qnorm(p, mean = eta, sd = sqrt(as.numeric(params["sill"]))),
    "12" = eta + sqrt(as.numeric(params["sill"])) *
      qt(p, df = 1 / as.numeric(params["df"])),
    "22" = qlnorm(p, meanlog = eta - as.numeric(params["sill"]) / 2,
                   sdlog = sqrt(as.numeric(params["sill"]))),
    "30" = qpois(p, lambda = exp(eta)),
    "21" = exp(eta) * qgamma(p, shape = as.numeric(params["shape"]) / 2,
                              rate = as.numeric(params["shape"]) / 2),
    "26" = exp(eta) * qweibull(
      p, shape = as.numeric(params["shape"]),
      scale = 1 / gamma(1 + 1 / as.numeric(params["shape"]))
    ),
    "11" = {
      n <- as.numeric(params["n"])
      qbinom(p, size = n, prob = pnorm(eta))
    },
    "16" = {
      n <- as.numeric(params["n"])
      qnbinom(p, size = n, prob = pnorm(eta))
    },
    "28" = {
      lo <- as.numeric(params["min"]); hi <- as.numeric(params["max"])
      lo + (hi - lo) * qbeta(p, shape1 = as.numeric(params["shape1"]),
                             shape2 = as.numeric(params["shape2"]))
    },
    "50" = {
      lo <- as.numeric(params["min"]); hi <- as.numeric(params["max"])
      m <- plogis(eta); phi <- as.numeric(params["shape"])
      lo + (hi - lo) * qbeta(p, shape1 = m * phi, shape2 = (1 - m) * phi)
    },
    "33" = {
      lo <- as.numeric(params["min"]); hi <- as.numeric(params["max"])
      a <- as.numeric(params["shape1"]); b <- as.numeric(params["shape2"])
      lo + (hi - lo) * (1 - (1 - p)^(1 / b))^(1 / a)
    },
    "42" = {
      lo <- as.numeric(params["min"]); hi <- as.numeric(params["max"])
      a <- as.numeric(params["shape"]); med <- plogis(eta)
      b <- log(0.5) / log1p(-med^a)
      lo + (hi - lo) * (1 - (1 - p)^(1 / b))^(1 / a)
    },
    "25" = qlogis(p, location = eta, scale = sqrt(as.numeric(params["sill"]))),
    "59" = {
      sk <- as.numeric(params["skew"]); ss <- sqrt(as.numeric(params["sill"]))
      out <- numeric(length(p)); right <- p >= sk
      out[right] <- eta - (ss / sk) *
        (log1p(-p[right]) - log1p(-sk))
      out[!right] <- eta + (ss / (1 - sk)) *
        (log(p[!right]) - log(sk))
      out
    },
    "18" = eta + sqrt(as.numeric(params["sill"])) *
      sn::qst(p, xi = 0, omega = 1, alpha = as.numeric(params["skew"]),
              nu = 1 / as.numeric(params["df"])),
    "10" = eta + sqrt(as.numeric(params["sill"])) *
      sn::qsn(p, xi = 0,
              omega = sqrt((as.numeric(params["skew"])^2 + as.numeric(params["sill"])) /
                             as.numeric(params["sill"])),
              alpha = as.numeric(params["skew"]) / sqrt(as.numeric(params["sill"]))),
    "34" = {
      z <- qnorm(p)
      h <- suppressWarnings(as.numeric(params["tail"]))
      if (!length(h) || is.na(h)) h <- suppressWarnings(as.numeric(params["h"]))
      if (!is.finite(h) || h < 0 || h >= 0.5)
        stop("Tukeyh tail must lie in [0, 0.5).", call. = FALSE)
      eta + sqrt(as.numeric(params["sill"])) * z * exp(0.5 * h * z^2)
    },
    "40" = {
      z <- qnorm(p)
      hr <- as.numeric(params["tail1"]); hl <- as.numeric(params["tail2"])
      .GeoTukeyh2StandardMoments(hr, hl)
      h <- ifelse(z >= 0, hr, hl)
      eta + sqrt(as.numeric(params["sill"])) * z * exp(0.5 * h * z^2)
    },
    "20" = {
      z <- qnorm(p)
      d <- as.numeric(params["tail"]); e <- as.numeric(params["skew"])
      .GeoSASStandardMoments(e, d)
      eta + sqrt(as.numeric(params["sill"])) * sinh((asinh(z) + e) / d)
    },
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


# Cached standard-normal Gauss-Hermite rule and Hermite design matrix.
# This avoids rebuilding the same M x nquad basis in repeated kriging calls.
.GeoGaussianHermiteRule <- function(M, nquad) {
  M <- as.integer(M); nquad <- as.integer(nquad)
  key <- paste0("M", M, "|n", nquad)
  cache <- .GeoModels_get_cache("gaussian_hermite_rule")
  if (exists(key, envir = cache, inherits = FALSE))
    return(get(key, envir = cache, inherits = FALSE))

  rule <- .GeoNormalQuadrature(nquad)
  H <- .GeoHermiteBasis(rule$nodes, M)
  out <- list(
    nodes = rule$nodes,
    weights = rule$weights,
    prob = pnorm(rule$nodes),
    H = H,
    WH = H * rule$weights
  )
  assign(key, out, envir = cache)
  out
}

# Fast Gauss-Hermite coefficients for continuous Gaussian-copula margins.
# The rule returned by .GeoNormalQuadrature() integrates expectations under
# a standard normal law directly, so all coefficients are obtained in one
# matrix product instead of M adaptive one-dimensional integrations.
.GeoGaussianCopulaAKGH <- function(M, model_type, params, nquad = NULL) {
  M <- as.integer(M)
  model_type <- as.integer(model_type)
  if (length(M) != 1L || is.na(M) || M < 1L)
    stop("M must be a positive integer.", call. = FALSE)

  if (is.null(nquad)) {
    ## Skew-Laplace has a kink in its quantile function at p = skew and needs
    ## a denser rule.  Smooth continuous margins are already very accurate
    ## with a substantially smaller rule.
    nquad <- if (model_type == 59L) 256L else 80L
  }
  nquad <- as.integer(nquad)

  rule <- .GeoGaussianHermiteRule(M, nquad)
  qq <- quantile_disp(rule$prob, model_type, params)
  if (any(!is.finite(qq)))
    stop("Non-finite marginal quantiles while computing Gaussian-copula covariance.",
         call. = FALSE)

  as.numeric(crossprod(rule$WH, qq))
}

# Vectorized polynomial evaluator for
#   sum_{k=1}^M coefficient[k] * rho^k.
# Horner evaluation uses M vector operations, rather than one R callback for
# every covariance entry.
.GeoGaussianCopulaPolynomial <- function(rho, coefficient) {
  rho <- as.numeric(rho)
  if (!length(rho)) return(numeric(0))

  good <- is.finite(rho) & !is.na(rho)
  out <- rep(NA_real_, length(rho))
  if (!any(good)) return(out)

  rr <- rho[good]
  coefficient <- as.numeric(coefficient)
  M <- length(coefficient)
  if (!M) {
    out[good] <- 0
    return(out)
  }

  acc <- rep(coefficient[M], length(rr))
  if (M > 1L) {
    for (k in seq.int(M - 1L, 1L))
      acc <- coefficient[k] + rr * acc
  }
  val <- rr * acc
  val[abs(rr) < 1e-12] <- 0
  val[!is.finite(val)] <- 0
  out[good] <- val
  out
}

# Gaussian-copula covariance.  Continuous margins use one cached
# Gauss-Hermite coefficient calculation and vectorized Horner evaluation.
# Legacy/discrete margins retain the adaptive-integration coefficient path.
gaussian_copula_cov <- function(rho, model_type, nuisance, M = 30,
                                cache_ak = TRUE, auto_reduce_M = TRUE,
                                verbose = FALSE, nquad = NULL) {
  model_type <- as.integer(model_type)
  M <- as.integer(M)
  continuous_models <- c(1L, 12L, 20L, 21L, 22L, 25L, 26L,
                         28L, 33L, 34L, 40L, 42L, 50L, 59L)
  use_gh <- model_type %in% continuous_models

  ## For a Gaussian margin the Gaussian-copula covariance is exactly
  ## sill * rho, so no quadrature or series evaluation is needed.
  if (model_type == 1L) {
    sill <- suppressWarnings(as.numeric(unlist(nuisance["sill"], use.names = FALSE)))
    if (length(sill) == 1L && is.finite(sill))
      return(.GeoGaussianCopulaPolynomial(rho, sill))
  }

  if (use_gh && is.null(nquad))
    nquad <- if (model_type == 59L) 256L else 80L

  cache_key <- paste(
    if (use_gh) "gh" else "adaptive",
    model_type,
    paste(names(nuisance), nuisance, sep = "=", collapse = ";"),
    M,
    if (use_gh) as.integer(nquad) else "",
    sep = "|"
  )
  ak_cache <- .GeoModels_get_cache("ak_cache")

  if (cache_ak && exists(cache_key, envir = ak_cache, inherits = FALSE)) {
    ak_coeffs <- get(cache_key, envir = ak_cache, inherits = FALSE)
    if (verbose) message("Coefficients a_k loaded from cache")
  } else {
    if (verbose) {
      method <- if (use_gh) paste0("Gauss-Hermite (n = ", nquad, ")") else "adaptive integration"
      message("Computing coefficients a_k with ", method, " for M = ", M)
    }

    ak_coeffs <- if (use_gh) {
      .GeoGaussianCopulaAKGH(M, model_type, nuisance, nquad = nquad)
    } else {
      compute_ak_vectorized(M, model_type, nuisance)
    }

    if (auto_reduce_M) {
      significant_coeffs <- which(abs(ak_coeffs) > 1e-10)
      if (length(significant_coeffs) > 0L) {
        effective_M <- max(significant_coeffs)
        if (effective_M < M * 0.7) {
          M_new <- min(effective_M + 5L, M)
          ak_coeffs <- ak_coeffs[seq_len(M_new)]
        }
      }
    }

    if (cache_ak)
      assign(cache_key, ak_coeffs, envir = ak_cache)
  }

  if (length(ak_coeffs) >= 10L) {
    last_coeffs <- abs(tail(ak_coeffs, 5L))
    if (all(last_coeffs < 1e-8) && verbose)
      message("Series converges well - last coefficients are small")
  }

  k <- seq_along(ak_coeffs)
  coefficient <- (ak_coeffs^2) / factorial(k)
  covariances <- .GeoGaussianCopulaPolynomial(rho, coefficient)

  if (verbose && any(is.finite(covariances))) {
    rr <- range(covariances[is.finite(covariances)])
    message("Computed with ", length(ak_coeffs), " coefficients")
    message("Covariance range: ", round(rr[1L], 6), " - ", round(rr[2L], 6))
  }

  covariances
}


################################################################################
## Gaussian-copula covariance helpers used by GeoCovmatrix/GeoKrig.           ##
################################################################################
.GeoGaussianCopulaContinuousModels <- c(1L, 12L, 20L, 21L, 22L, 25L, 26L,
                                         28L, 33L, 34L, 40L, 42L, 50L, 59L)
## Discrete count margins supported by the observed-scale copula covariance
## engines.  They are kept separate because, unlike most continuous margins,
## all three have location-dependent quantile functions.
.GeoGaussianCopulaDiscreteModels <- c(11L, 16L, 30L) # Binomial, BinomialNeg, Poisson
.GeoCopulaPredictionModels <- c(.GeoGaussianCopulaContinuousModels,
                                .GeoGaussianCopulaDiscreteModels)

.GeoCopulaParamScalar <- function(params, name, required = TRUE, default = NULL) {
  nms <- names(params)
  if (is.null(nms) || !(name %in% nms)) {
    if (required) stop("Missing copula marginal parameter '", name, "'.", call. = FALSE)
    return(default)
  }
  x <- suppressWarnings(as.numeric(params[[name]]))
  if (length(x) != 1L || !is.finite(x))
    stop("Copula marginal parameter '", name, "' must be a finite scalar.", call. = FALSE)
  x
}

.GeoCopulaEtaVector <- function(eta, n, default = 0) {
  if (is.null(eta) || !length(eta)) eta <- default
  eta <- as.numeric(eta)
  if (length(eta) == 1L) eta <- rep(eta, n)
  if (length(eta) != n || any(!is.finite(eta)))
    stop("Internal copula location vector has incompatible length or non-finite values.",
         call. = FALSE)
  eta
}

.GeoCopulaMarginalMoments <- function(model_type, params, eta = 0) {
  model_type <- as.integer(model_type)
  eta <- as.numeric(eta)
  if (!length(eta) || any(!is.finite(eta)))
    stop("Copula marginal locations must be finite.", call. = FALSE)

  sill <- function() .GeoCopulaParamScalar(params, "sill")
  lohi <- function() {
    lo <- .GeoCopulaParamScalar(params, "min")
    hi <- .GeoCopulaParamScalar(params, "max")
    if (lo >= hi) stop("Copula bounded margins require min < max.", call. = FALSE)
    c(lo, hi)
  }

  if (model_type == 1L) {                         # Gaussian
    ss <- sill()
    return(list(mean = eta, variance = rep(ss, length(eta))))
  }
  if (model_type %in% c(20L, 34L, 40L)) {         # SAS / Tukeyh / Tukeyh2
    return(.GeoPlainTransformedMarginalMoments(model_type, params, eta))
  }
  if (model_type == 12L) {                        # StudentT
    ss <- sill(); eps <- .GeoCopulaParamScalar(params, "df"); df_t <- 1 / eps
    if (df_t <= 2) stop("StudentT linear prediction requires more than 2 degrees of freedom.", call. = FALSE)
    return(list(mean = eta, variance = rep(ss * df_t / (df_t - 2), length(eta))))
  }
  if (model_type == 22L) {                        # LogGaussian
    ss <- sill(); mm <- exp(eta)
    return(list(mean = mm, variance = mm^2 * (exp(ss) - 1)))
  }
  if (model_type == 30L) {                        # Poisson
    mm <- exp(eta)
    return(list(mean = mm, variance = mm))
  }
  if (model_type == 11L) {                        # Binomial
    nn <- .GeoCopulaParamScalar(params, "n")
    if (nn < 1 || abs(nn - round(nn)) > sqrt(.Machine$double.eps))
      stop("Binomial copula prediction requires a positive integer n.", call. = FALSE)
    pp <- pnorm(eta)
    return(list(mean = nn * pp, variance = nn * pp * (1 - pp)))
  }
  if (model_type == 16L) {                        # BinomialNeg
    nn <- .GeoCopulaParamScalar(params, "n")
    if (nn < 1 || abs(nn - round(nn)) > sqrt(.Machine$double.eps))
      stop("BinomialNeg copula prediction requires a positive integer n.", call. = FALSE)
    pp <- pnorm(eta)
    mm <- nn * (1 - pp) / pp
    return(list(mean = mm, variance = nn * (1 - pp) / pp^2))
  }
  if (model_type == 21L) {                        # Gamma
    sh <- .GeoCopulaParamScalar(params, "shape")
    if (sh <= 0) stop("Gamma shape must be positive.", call. = FALSE)
    mm <- exp(eta)
    return(list(mean = mm, variance = 2 * mm^2 / sh))
  }
  if (model_type == 26L) {                        # Weibull
    sh <- .GeoCopulaParamScalar(params, "shape")
    if (sh <= 0) stop("Weibull shape must be positive.", call. = FALSE)
    mm <- exp(eta)
    vv0 <- gamma(1 + 2 / sh) / gamma(1 + 1 / sh)^2 - 1
    return(list(mean = mm, variance = mm^2 * vv0))
  }
  if (model_type == 25L) {                        # Logistic
    ss <- sill()
    return(list(mean = eta, variance = rep(pi^2 * ss / 3, length(eta))))
  }
  if (model_type == 59L) {                        # SkewLaplace
    ss <- sill(); sk <- .GeoCopulaParamScalar(params, "skew")
    if (sk <= 0 || sk >= 1) stop("SkewLaplace skew must lie in (0,1).", call. = FALSE)
    shift <- sqrt(ss) * (1 / sk - 1 / (1 - sk))
    vv <- ss * (1 - 2 * sk + 2 * sk^2) / (sk^2 * (1 - sk)^2)
    return(list(mean = eta + shift, variance = rep(vv, length(eta))))
  }
  if (model_type == 28L) {                        # Beta
    lr <- lohi(); a <- .GeoCopulaParamScalar(params, "shape1")
    b <- .GeoCopulaParamScalar(params, "shape2")
    if (a <= 0 || b <= 0) stop("Beta shape parameters must be positive.", call. = FALSE)
    rr <- lr[2L] - lr[1L]
    mm <- lr[1L] + rr * a / (a + b)
    vv <- rr^2 * a * b / ((a + b)^2 * (a + b + 1))
    return(list(mean = rep(mm, length(eta)), variance = rep(vv, length(eta))))
  }
  if (model_type == 50L) {                        # Beta2
    lr <- lohi(); phi <- .GeoCopulaParamScalar(params, "shape")
    if (phi <= 0) stop("Beta2 shape must be positive.", call. = FALSE)
    m <- plogis(eta); rr <- lr[2L] - lr[1L]
    return(list(mean = lr[1L] + rr * m,
                variance = rr^2 * m * (1 - m) / (phi + 1)))
  }
  if (model_type == 33L) {                        # Kumaraswamy
    lr <- lohi(); a <- .GeoCopulaParamScalar(params, "shape1")
    b <- .GeoCopulaParamScalar(params, "shape2")
    if (a <= 0 || b <= 0) stop("Kumaraswamy shape parameters must be positive.", call. = FALSE)
    m1 <- b * beta(1 + 1 / a, b); m2 <- b * beta(1 + 2 / a, b)
    rr <- lr[2L] - lr[1L]
    return(list(mean = rep(lr[1L] + rr * m1, length(eta)),
                variance = rep(rr^2 * (m2 - m1^2), length(eta))))
  }
  if (model_type == 42L) {                        # Kumaraswamy2
    lr <- lohi(); a <- .GeoCopulaParamScalar(params, "shape")
    if (a <= 0) stop("Kumaraswamy2 shape must be positive.", call. = FALSE)
    med <- plogis(eta)
    b <- log(0.5) / log1p(-med^a)
    if (any(!is.finite(b)) || any(b <= 0))
      stop("Kumaraswamy2 produced invalid location-specific shape parameters.", call. = FALSE)
    m1 <- b * beta(1 + 1 / a, b); m2 <- b * beta(1 + 2 / a, b)
    rr <- lr[2L] - lr[1L]
    return(list(mean = lr[1L] + rr * m1,
                variance = rr^2 * (m2 - m1^2)))
  }

  stop("Copula marginal moments are not implemented for model id ", model_type, ".",
       call. = FALSE)
}

.GeoNormalQuadrature <- function(n = 60L) {
  n <- as.integer(n)
  key <- paste0("n", n)
  cache <- .GeoModels_get_cache("normal_quadrature")
  if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache))
  if (n < 10L) stop("Internal Gaussian quadrature requires at least 10 nodes.", call. = FALSE)
  J <- matrix(0, n, n)
  off <- sqrt(seq_len(n - 1L))
  J[cbind(seq_len(n - 1L), 2:n)] <- off
  J[cbind(2:n, seq_len(n - 1L))] <- off
  ee <- eigen(J, symmetric = TRUE)
  oo <- order(ee$values)
  out <- list(nodes = ee$values[oo], weights = ee$vectors[1L, oo]^2)
  assign(key, out, envir = cache)
  out
}

.GeoHermiteBasis <- function(x, M) {
  x <- as.numeric(x); M <- as.integer(M)
  H <- matrix(0, nrow = length(x), ncol = M)
  if (M < 1L) return(H)
  H[, 1L] <- x
  if (M >= 2L) {
    H[, 2L] <- x^2 - 1
    if (M >= 3L) {
      hm2 <- rep(1, length(x)); hm1 <- x
      for (k in 2:M) {
        hk <- x * hm1 - (k - 1) * hm2
        H[, k] <- hk
        hm2 <- hm1; hm1 <- hk
      }
    }
  }
  H
}

.GeoCopulaHermiteGHIndexed <- function(model_type, params, eta,
                                             M = 30L, nquad = 60L) {
  eta <- as.numeric(eta)
  uq <- unique(eta)
  rule <- .GeoGaussianHermiteRule(M, nquad)
  p <- rule$prob
  WH <- rule$WH
  out_u <- matrix(NA_real_, nrow = length(uq), ncol = M)

  cache <- .GeoModels_get_cache("copula_hermite_gh")
  base_names <- setdiff(names(params), .GeoMean_names(max(1L, sum(grepl("^mean", names(params))))))
  base_key <- paste(model_type,
                    paste(base_names, unlist(params[base_names]), sep = "=", collapse = ";"),
                    M, nquad, sep = "|")

  for (i in seq_along(uq)) {
    key <- paste0(base_key, "|eta=", format(uq[i], digits = 17, scientific = TRUE))
    if (exists(key, envir = cache, inherits = FALSE)) {
      out_u[i, ] <- get(key, envir = cache)
    } else {
      pp <- params
      pp["mean"] <- uq[i]
      q <- quantile_disp(p, model_type, pp)
      if (any(!is.finite(q)))
        stop("Non-finite marginal quantiles while computing Gaussian-copula covariance.", call. = FALSE)
      aa <- as.numeric(crossprod(WH, q))
      out_u[i, ] <- aa
      assign(key, aa, envir = cache)
    }
  }
  list(coef = out_u, index = match(eta, uq))
}

.GeoCopulaHermiteGH <- function(model_type, params, eta, M = 30L, nquad = 60L) {
  z <- .GeoCopulaHermiteGHIndexed(model_type, params, eta, M = M, nquad = nquad)
  z$coef[z$index, , drop = FALSE]
}

# Native location-dependent Gaussian-copula series evaluator.  The historical
# R implementation expanded one M-column Hermite matrix for every covariance
# pair.  For dense covariance and kriging blocks that can be hundreds of MB
# even though the same endpoint appears in many pairs.  Here only the unique
# endpoint coefficient rows are retained and native code follows integer row
# indices while accumulating the identical Hermite series.
.GeoGaussianCovSeriesIndexedC <- function(rho, z1, z2) {
  rho <- as.numeric(rho)
  n <- length(rho)
  if (!n) return(numeric(0))
  A1 <- as.matrix(z1$coef); A2 <- as.matrix(z2$coef)
  idx1 <- as.integer(z1$index); idx2 <- as.integer(z2$index)
  degree <- ncol(A1)
  if (degree < 1L || ncol(A2) != degree || length(idx1) != n ||
      length(idx2) != n || anyNA(idx1) || anyNA(idx2) ||
      any(idx1 < 1L | idx1 > nrow(A1)) || any(idx2 < 1L | idx2 > nrow(A2)))
    stop("Internal indexed Gaussian-copula series dimensions are inconsistent.",
         call. = FALSE)

  dotCall64::.C64(
    "GeoGaussianCovSeriesIndexed",
    SIGNATURE = c("double", "integer", "double", "integer", "integer",
                  "double", "integer", "integer", "integer", "double",
                  "double"),
    rho, as.integer(n), as.double(A1), as.integer(nrow(A1)), idx1,
    as.double(A2), as.integer(nrow(A2)), idx2, as.integer(degree),
    factorial(seq_len(degree)),
    out = dotCall64::vector_dc("double", n),
    INTENT = c(rep("r", 10), "w"), NAOK = TRUE, PACKAGE = "GeoModels"
  )$out
}

.GeoGaussianCopulaCrossCov <- function(rho, model_type, params,
                                       eta1 = 0, eta2 = eta1,
                                       M = 30L, nquad = 60L) {
  model_type <- as.integer(model_type)
  rho <- as.numeric(rho)
  n <- length(rho)
  if (!n) return(numeric(0))
  eta1 <- .GeoCopulaEtaVector(eta1, n)
  eta2 <- .GeoCopulaEtaVector(eta2, n)
  if (any(!is.finite(rho)) || any(abs(rho) > 1 + 1e-10))
    stop("Invalid latent correlations in Gaussian-copula covariance.", call. = FALSE)
  rho <- pmax(-1, pmin(1, rho))

  ## For additive location families the covariance does not depend on eta.
  additive <- c(1L, 12L, 20L, 25L, 34L, 40L, 59L)
  ## Positive multiplicative families have q_eta = exp(eta) q_0.
  multiplicative <- c(21L, 22L, 26L)
  ## Fixed-shape bounded families do not depend on the regression location.
  fixed_bounded <- c(28L, 33L)

  pp <- params
  pp["mean"] <- 0

  ## Tukey-h, Tukey-hh and sinh-arcsinh are monotone transforms of a
  ## single Gaussian RF.  Under a Gaussian copula their covariance is
  ## therefore exactly the same as in the corresponding plain model; do
  ## not approximate it with a truncated Hermite expansion.
  if (model_type %in% c(20L, 34L, 40L)) {
    return(.GeoPlainTransformedCrossCov(rho, model_type, pp))
  }
  if (model_type %in% c(additive, fixed_bounded)) {
    return(gaussian_copula_cov(rho, model_type, pp, M = M))
  }
  if (model_type %in% multiplicative) {
    base <- gaussian_copula_cov(rho, model_type, pp, M = M)
    return(exp(eta1 + eta2) * base)
  }
  if (model_type %in% c(50L, 42L, .GeoGaussianCopulaDiscreteModels)) {
    ## Count quantiles are step functions.  A denser normal quadrature reduces
    ## aliasing while preserving the fast indexed/cached covariance path.
    nquad_loc <- if (model_type %in% .GeoGaussianCopulaDiscreteModels)
      max(as.integer(nquad), 120L) else as.integer(nquad)
    if (isTRUE(getOption("GeoModels.gaussian.indexed_series", TRUE))) {
      z1 <- .GeoCopulaHermiteGHIndexed(
        model_type, params, eta1, M = M, nquad = nquad_loc
      )
      z2 <- .GeoCopulaHermiteGHIndexed(
        model_type, params, eta2, M = M, nquad = nquad_loc
      )
      return(.GeoGaussianCovSeriesIndexedC(rho, z1, z2))
    }

    ## Historical reference path retained behind an internal option for exact
    ## regression comparisons.  It intentionally materializes pair-expanded
    ## Hermite matrices and should not be used for production-sized problems.
    A1 <- .GeoCopulaHermiteGH(model_type, params, eta1, M = M, nquad = nquad_loc)
    A2 <- .GeoCopulaHermiteGH(model_type, params, eta2, M = M, nquad = nquad_loc)
    fac <- factorial(seq_len(M))
    out <- numeric(n)
    rp <- rho
    for (k in seq_len(M)) {
      term <- (A1[, k] * A2[, k] / fac[k]) * rp
      term[!is.finite(term)] <- 0
      out <- out + term
      rp <- rp * rho
    }
    return(out)
  }

  ## Backward-compatible fallback for legacy Gaussian-copula margins.
  gaussian_copula_cov(rho, model_type, params, M = M)
}

.GeoCopulaMarginalMeanDerivative <- function(model_type, params, eta) {
  model_type <- as.integer(model_type); eta <- as.numeric(eta)
  if (model_type %in% c(1L, 12L, 20L, 25L, 34L, 40L, 59L)) return(rep(1, length(eta)))
  if (model_type %in% c(21L, 22L, 26L, 30L)) return(exp(eta))
  if (model_type == 11L) {
    nn <- .GeoCopulaParamScalar(params, "n")
    return(nn * dnorm(eta))
  }
  if (model_type == 16L) {
    nn <- .GeoCopulaParamScalar(params, "n")
    pp <- pnorm(eta)
    return(-nn * dnorm(eta) / pp^2)
  }
  if (model_type == 50L) {
    lo <- .GeoCopulaParamScalar(params, "min"); hi <- .GeoCopulaParamScalar(params, "max")
    m <- plogis(eta)
    return((hi - lo) * m * (1 - m))
  }
  if (model_type %in% c(28L, 33L)) return(rep(0, length(eta)))
  if (model_type == 42L) {
    ## Stable centered finite-difference derivative of the Kumaraswamy2 mean.
    h <- 1e-5 * pmax(1, abs(eta))
    mp <- .GeoCopulaMarginalMoments(model_type, params, eta + h)$mean
    mm <- .GeoCopulaMarginalMoments(model_type, params, eta - h)$mean
    return((mp - mm) / (2 * h))
  }
  rep(1, length(eta))
}




################################################################################
## Fast Clayton-like copula covariance helpers used by GeoCovmatrix/GeoKrig.  ##
##                                                                            ##
## The bivariate Clayton-like density is expensive because it contains an      ##
## Appell F4 evaluation.  Kriging must therefore never call adaptive 2-D       ##
## integration separately for every covariance entry.  The implementation     ##
## below builds and caches a deterministic Gauss-Legendre representation of   ##
## the copula joint law on a modest grid of rho^2 values.  Marginal quantile   ##
## transforms are then applied with BLAS matrix products and the resulting     ##
## covariance table is interpolated.                                           ##
################################################################################

.GeoUnitQuadrature <- function(n = 36L) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 12L)
    stop("Internal unit-interval quadrature requires at least 12 nodes.", call. = FALSE)
  key <- paste0("n", n)
  cache <- .GeoModels_get_cache("unit_quadrature")
  if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache))

  k <- seq_len(n - 1L)
  off <- k / sqrt(4 * k^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1L)] <- off
  J[cbind(k + 1L, k)] <- off
  ee <- eigen(J, symmetric = TRUE)
  oo <- order(ee$values)
  x <- ee$values[oo]
  w <- 2 * ee$vectors[1L, oo]^2
  out <- list(nodes = (x + 1) / 2, weights = w / 2)
  assign(key, out, envir = cache)
  out
}

.GeoBalanceCopulaQuadrature <- function(P, w, tol = 1e-11, maxit = 80L) {
  P <- (P + t(P)) / 2
  P[!is.finite(P) | P < 0] <- 0
  tiny <- .Machine$double.xmin
  if (!any(P > 0)) stop("Copula quadrature produced zero mass.", call. = FALSE)

  for (it in seq_len(maxit)) {
    rs <- rowSums(P)
    s <- sqrt(w / pmax(rs, tiny))
    P <- sweep(P, 1L, s, `*`)
    P <- sweep(P, 2L, s, `*`)
    if (max(abs(rowSums(P) - w)) < tol) break
  }
  ## One final symmetric correction substantially reduces accumulated roundoff.
  P <- (P + t(P)) / 2
  rs <- rowSums(P)
  s <- sqrt(w / pmax(rs, tiny))
  P <- sweep(P, 1L, s, `*`)
  P <- sweep(P, 2L, s, `*`)
  P / sum(P)
}

.GeoClaytonKernelTable <- function(nu, nquad = 36L, nrho = 55L) {
  nu <- suppressWarnings(as.numeric(nu))
  if (length(nu) != 1L || !is.finite(nu) || nu <= 0 ||
      abs(nu - round(nu)) > sqrt(.Machine$double.eps))
    stop("Clayton-like linear prediction requires a positive integer nu.", call. = FALSE)
  nu <- as.integer(round(nu))
  nquad <- as.integer(nquad); nrho <- as.integer(nrho)
  if (nrho < 15L) stop("Internal Clayton covariance grid is too small.", call. = FALSE)

  ## The interior grid is concentrated near both independence and the
  ## comonotone limit.  rho^2=1 is handled analytically, avoiding evaluation
  ## of the singular limiting density.
  tt <- seq(0, 1, length.out = nrho - 1L)
  xgrid <- c(0.995 * 0.5 * (1 - cos(pi * tt)), 1)
  use_symmetry <- isTRUE(getOption("GeoModels.clayton.symmetric_kernel", TRUE))
  kernel_mode <- if (use_symmetry) "tri" else "full"
  key <- paste0("nu=", nu, "|q=", nquad, "|r=", nrho,
                "|v=3|mode=", kernel_mode)
  cache <- .GeoModels_get_cache("clayton_kernel_table")
  if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache))

  rule <- .GeoUnitQuadrature(nquad)
  ww <- outer(rule$weights, rule$weights)

  ## The Clayton-like copula density is exchangeable in (u, v).  Evaluate
  ## only one triangular half of the quadrature grid, then mirror it before
  ## balancing.  For the default 36-point rule this reduces each expensive
  ## Appell-F4 density slice from 1296 to 666 native evaluations.
  if (use_symmetry) {
    ij <- which(upper.tri(matrix(FALSE, nquad, nquad), diag = TRUE),
                arr.ind = TRUE)
    u <- rule$nodes[ij[, 1L]]
    v <- rule$nodes[ij[, 2L]]
    density_evaluations <- nrow(ij)
  } else {
    u <- rep(rule$nodes, times = nquad)
    v <- rep(rule$nodes, each = nquad)
    ij <- NULL
    density_evaluations <- nquad * nquad
  }

  joint <- array(0, dim = c(nquad, nquad, nrho))
  joint[, , 1L] <- ww
  joint[, , nrho] <- diag(rule$weights, nrow = nquad, ncol = nquad)

  if (nrho > 2L) {
    for (j in 2:(nrho - 1L)) {
      rho <- sqrt(xgrid[j])
      ans <- .C(
        "GeoCopulaDensity",
        as.integer(2L), as.double(nu),
        as.double(u), as.double(v), as.integer(length(u)),
        as.double(rho), res = double(length(u)),
        PACKAGE = "GeoModels", NAOK = TRUE
      )$res

      if (use_symmetry) {
        dens <- matrix(0, nrow = nquad, ncol = nquad)
        dens[ij] <- ans
        dens[cbind(ij[, 2L], ij[, 1L])] <- ans
      } else {
        dens <- matrix(ans, nrow = nquad, ncol = nquad)
      }

      P <- dens * ww
      joint[, , j] <- .GeoBalanceCopulaQuadrature(P, rule$weights)
    }
  }

  out <- list(x = xgrid, joint = joint, nodes = rule$nodes,
              weights = rule$weights, nu = nu,
              symmetric_kernel = use_symmetry,
              density_evaluations_per_slice = density_evaluations)
  assign(key, out, envir = cache)
  out
}

.GeoClaytonQuantileRows <- function(model_type, params, eta, nodes, weights) {
  eta <- as.numeric(eta)
  Q <- matrix(NA_real_, nrow = length(eta), ncol = length(nodes))
  V <- .GeoCopulaMarginalMoments(model_type, params, eta)$variance
  if (length(V) == 1L) V <- rep(V, length(eta))

  for (i in seq_along(eta)) {
    pp <- params
    pp["mean"] <- eta[i]
    q <- quantile_disp(nodes, model_type, pp)
    if (any(!is.finite(q)))
      stop("Non-finite marginal quantiles in Clayton-like covariance.", call. = FALSE)
    qc <- q - sum(weights * q)
    vq <- sum(weights * qc^2)
    if (!is.finite(vq) || vq <= 0 || !is.finite(V[i]) || V[i] <= 0)
      stop("Invalid marginal variance in Clayton-like covariance.", call. = FALSE)
    ## Calibrate the deterministic quadrature to the exact marginal variance.
    Q[i, ] <- qc * sqrt(V[i] / vq)
  }
  Q
}

.GeoInterp1Grid <- function(grid, values, x) {
  x <- pmax(grid[1L], pmin(grid[length(grid)], as.numeric(x)))
  ii <- findInterval(x, grid, all.inside = TRUE)
  ii <- pmin(ii, length(grid) - 1L)
  x0 <- grid[ii]; x1 <- grid[ii + 1L]
  a <- ifelse(x1 > x0, (x - x0) / (x1 - x0), 0)
  (1 - a) * values[ii] + a * values[ii + 1L]
}

.GeoInterp3Grid <- function(g1, g2, gx, tab, z1, z2, x) {
  z1 <- pmax(g1[1L], pmin(g1[length(g1)], as.numeric(z1)))
  z2 <- pmax(g2[1L], pmin(g2[length(g2)], as.numeric(z2)))
  x  <- pmax(gx[1L], pmin(gx[length(gx)], as.numeric(x)))

  i1 <- pmin(findInterval(z1, g1, all.inside = TRUE), length(g1) - 1L)
  i2 <- pmin(findInterval(z2, g2, all.inside = TRUE), length(g2) - 1L)
  ir <- pmin(findInterval(x,  gx, all.inside = TRUE), length(gx) - 1L)
  a1 <- (z1 - g1[i1]) / pmax(g1[i1 + 1L] - g1[i1], .Machine$double.eps)
  a2 <- (z2 - g2[i2]) / pmax(g2[i2 + 1L] - g2[i2], .Machine$double.eps)
  ar <- (x  - gx[ir]) / pmax(gx[ir + 1L] - gx[ir], .Machine$double.eps)

  at <- function(d1, d2, dr) tab[cbind(i1 + d1, i2 + d2, ir + dr)]
  c000 <- at(0L,0L,0L); c100 <- at(1L,0L,0L)
  c010 <- at(0L,1L,0L); c110 <- at(1L,1L,0L)
  c001 <- at(0L,0L,1L); c101 <- at(1L,0L,1L)
  c011 <- at(0L,1L,1L); c111 <- at(1L,1L,1L)
  c00 <- (1-a1)*c000 + a1*c100
  c10 <- (1-a1)*c010 + a1*c110
  c01 <- (1-a1)*c001 + a1*c101
  c11 <- (1-a1)*c011 + a1*c111
  c0 <- (1-a2)*c00 + a2*c10
  c1 <- (1-a2)*c01 + a2*c11
  (1-ar)*c0 + ar*c1
}



# Native Clayton interpolation helpers.  The historical R interpolators above
# are retained as regression references.  These wrappers move interval search,
# interpolation, and (for multiplicative margins) exp(eta1 + eta2) scaling into
# one native pass, avoiding large temporary vectors on warm-cache covariance
# calls.
.GeoClaytonInterp1GridC <- function(grid, values, x,
                                     eta1 = NULL, eta2 = NULL,
                                     multiplicative = FALSE) {
  grid <- as.numeric(grid)
  values <- as.numeric(values)
  x <- as.numeric(x)
  ng <- length(grid); n <- length(x)
  if (!n) return(numeric(0))
  if (ng < 2L || length(values) != ng)
    stop("Internal Clayton interpolation table has incompatible dimensions.",
         call. = FALSE)

  multiplicative <- isTRUE(multiplicative)
  if (multiplicative) {
    eta1 <- .GeoCopulaEtaVector(eta1, n)
    eta2 <- .GeoCopulaEtaVector(eta2, n)
  } else {
    ## C64 needs valid addresses, although the native routine does not
    ## dereference eta1/eta2 when multiplicative == 0.
    eta1 <- 0
    eta2 <- 0
  }

  dotCall64::.C64(
    "GeoClaytonInterp1",
    SIGNATURE = c("double", "integer", "double", "double", "double",
                  "double", "integer", "integer", "double"),
    grid, as.integer(ng), values, x, eta1, eta2, as.integer(n),
    as.integer(multiplicative),
    out = dotCall64::vector_dc("double", n),
    INTENT = c("r", "r", "r", "r", "r", "r", "r", "r", "w"),
    NAOK = TRUE, PACKAGE = "GeoModels"
  )$out
}

.GeoClaytonInterp3GridC <- function(eta_grid, x_grid, tab,
                                     eta1, eta2, x) {
  eta_grid <- as.numeric(eta_grid)
  x_grid <- as.numeric(x_grid)
  x <- as.numeric(x)
  n <- length(x); ne <- length(eta_grid); nx <- length(x_grid)
  if (!n) return(numeric(0))
  if (ne < 2L || nx < 2L || length(tab) != ne * ne * nx)
    stop("Internal Clayton location interpolation table has incompatible dimensions.",
         call. = FALSE)
  eta1 <- .GeoCopulaEtaVector(eta1, n)
  eta2 <- .GeoCopulaEtaVector(eta2, n)

  dotCall64::.C64(
    "GeoClaytonInterp3",
    SIGNATURE = c("double", "integer", "double", "integer", "double",
                  "double", "double", "double", "integer", "double"),
    eta_grid, as.integer(ne), x_grid, as.integer(nx), tab,
    eta1, eta2, x, as.integer(n),
    out = dotCall64::vector_dc("double", n),
    INTENT = c("r", "r", "r", "r", "r", "r", "r", "r", "r", "w"),
    NAOK = TRUE, PACKAGE = "GeoModels"
  )$out
}

.GeoClaytonInterp1Grid <- function(grid, values, x,
                                    eta1 = NULL, eta2 = NULL,
                                    multiplicative = FALSE) {
  if (isTRUE(getOption("GeoModels.clayton.native_lookup", TRUE))) {
    return(.GeoClaytonInterp1GridC(
      grid, values, x, eta1 = eta1, eta2 = eta2,
      multiplicative = multiplicative
    ))
  }
  out <- .GeoInterp1Grid(grid, values, x)
  if (isTRUE(multiplicative))
    out <- exp(.GeoCopulaEtaVector(eta1, length(out)) +
               .GeoCopulaEtaVector(eta2, length(out))) * out
  out
}

.GeoClaytonInterp3Grid <- function(eta_grid, x_grid, tab,
                                    eta1, eta2, x) {
  if (isTRUE(getOption("GeoModels.clayton.native_lookup", TRUE))) {
    return(.GeoClaytonInterp3GridC(
      eta_grid, x_grid, tab, eta1 = eta1, eta2 = eta2, x = x
    ))
  }
  .GeoInterp3Grid(eta_grid, eta_grid, x_grid, tab, eta1, eta2, x)
}

.GeoClaytonCovLookupTable <- function(model_type, params, nu,
                                           nquad = 36L, nrho = 55L,
                                           neta = 51L, eta_range = NULL) {
  model_type <- as.integer(model_type)
  kernel <- .GeoClaytonKernelTable(nu, nquad = nquad, nrho = nrho)
  additive <- c(1L, 12L, 20L, 25L, 34L, 40L, 59L)
  multiplicative <- c(21L, 22L, 26L)
  fixed_bounded <- c(28L, 33L)
  location_dependent <- c(50L, 42L, .GeoGaussianCopulaDiscreteModels)

  nms <- names(params)
  keep <- if (is.null(nms)) seq_along(params) else !grepl("^mean[0-9]*$", nms)
  ppkey <- params[keep]
  if (!is.null(names(ppkey))) ppkey <- ppkey[order(names(ppkey))]
  base_key <- paste(
    "m", model_type,
    paste(names(ppkey), format(as.numeric(unlist(ppkey, use.names = FALSE)), digits = 17, scientific = TRUE),
          sep = "=", collapse = ";"),
    "nu", as.integer(round(nu)), "q", nquad, "r", nrho, "e", neta,
    "kernel", if (isTRUE(kernel$symmetric_kernel)) "tri" else "full",
    "v5", sep = "|"
  )
  cache <- .GeoModels_get_cache("clayton_cov_lookup")

  if (model_type %in% c(additive, multiplicative, fixed_bounded)) {
    key <- paste0(base_key, "|fixed")
    if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache))
    Q <- .GeoClaytonQuantileRows(model_type, params, 0,
                                 kernel$nodes, kernel$weights)
    cv <- numeric(length(kernel$x))
    for (j in seq_along(kernel$x))
      cv[j] <- as.numeric(Q %*% kernel$joint[, , j] %*% t(Q))
    out <- list(kind = if (model_type %in% multiplicative) "multiplicative" else "fixed",
                x = kernel$x, cv = cv)
    assign(key, out, envir = cache)
    return(out)
  }

  if (model_type %in% location_dependent) {
    ## Beta2/Kumaraswamy2 and the count margins have location-dependent
    ## quantile transforms.  Do not
    ## build a huge fixed grid reaching probabilities arbitrarily close to 0/1:
    ## low-order uniform quadrature cannot resolve such nearly-degenerate beta
    ## margins and can return zero numerical variance.  Instead build a compact
    ## eta grid around the locations actually requested.  The expensive Clayton
    ## copula kernel (Appell F4 evaluations) is still shared session-wide; only
    ## the inexpensive marginal transform table depends on this eta band.
    er <- suppressWarnings(as.numeric(eta_range))
    er <- er[is.finite(er)]
    if (!length(er)) er <- 0
    elo <- floor(min(er)) - 1
    ehi <- ceiling(max(er)) + 1
    if (!is.finite(elo) || !is.finite(ehi)) stop("Invalid Clayton marginal-location range.", call. = FALSE)
    if (ehi <= elo) ehi <- elo + 2
    ## Keep a minimum two-unit band so nearby GeoKrigloc neighborhoods reuse
    ## the same cached table rather than generating a new table for each site.
    if ((ehi - elo) < 2) {
      mid <- 0.5 * (elo + ehi)
      elo <- floor(mid - 1)
      ehi <- ceiling(mid + 1)
    }
    key <- paste0(base_key, "|eta=", format(elo, scientific = FALSE), ":",
                  format(ehi, scientific = FALSE))
    if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache))

    ## Keep the interpolation spacing reasonably fine on wide regression
    ## ranges without changing the expensive copula kernel resolution.
    neta_eff <- max(as.integer(neta), as.integer(4L * ceiling(ehi - elo) + 1L))
    eg <- seq(elo, ehi, length.out = neta_eff)
    Q <- .GeoClaytonQuantileRows(model_type, params, eg,
                                 kernel$nodes, kernel$weights)
    ng <- length(eg); nr <- length(kernel$x)
    tab <- array(0, dim = c(ng, ng, nr))
    for (j in seq_len(nr)) {
      tmp <- Q %*% kernel$joint[, , j] %*% t(Q)
      tab[, , j] <- (tmp + t(tmp)) / 2
    }
    out <- list(kind = "location", x = kernel$x, eta = eg, tab = tab,
                eta_range = c(elo, ehi))
    assign(key, out, envir = cache)
    return(out)
  }

  stop("Unsupported margin in Clayton-like covariance.", call. = FALSE)
}

.GeoClaytonCopulaCrossCovFast <- function(rho, model_type, params,
                                           eta1 = 0, eta2 = eta1,
                                           nquad = 36L, nrho = 55L,
                                           neta = 51L, eta_range = NULL) {
  model_type <- as.integer(model_type)
  if (!(model_type %in% .GeoCopulaPredictionModels))
    stop("Clayton-like linear prediction is not implemented for model id ",
         model_type, ".", call. = FALSE)
  if (model_type %in% .GeoGaussianCopulaDiscreteModels)
    nquad <- max(as.integer(nquad), 48L)
  rho <- as.numeric(rho); n <- length(rho)
  if (!n) return(numeric(0))
  eta1 <- .GeoCopulaEtaVector(eta1, n)
  eta2 <- .GeoCopulaEtaVector(eta2, n)
  if (any(!is.finite(rho)) || any(abs(rho) > 1 + 1e-10))
    stop("Invalid latent correlations in Clayton-like covariance.", call. = FALSE)
  x <- pmin(1, pmax(0, rho^2))

  nu <- .GeoCopulaParamScalar(params, "nu")
  if (nu <= 0 || abs(nu - round(nu)) > sqrt(.Machine$double.eps))
    stop("Clayton-like linear prediction requires a positive integer nu.", call. = FALSE)
  lookup_eta_range <- if (is.null(eta_range)) {
    range(c(eta1, eta2), finite = TRUE)
  } else {
    eta_range
  }
  lookup <- .GeoClaytonCovLookupTable(
    model_type, params, nu,
    nquad = nquad, nrho = nrho, neta = neta,
    eta_range = lookup_eta_range
  )

  if (identical(lookup$kind, "fixed")) {
    out <- .GeoClaytonInterp1Grid(lookup$x, lookup$cv, x)
  } else if (identical(lookup$kind, "multiplicative")) {
    out <- .GeoClaytonInterp1Grid(
      lookup$x, lookup$cv, x, eta1 = eta1, eta2 = eta2,
      multiplicative = TRUE
    )
  } else {
    out <- .GeoClaytonInterp3Grid(
      lookup$eta, lookup$x, lookup$tab, eta1, eta2, x
    )
  }
  out[x < 1e-15] <- 0
  out
}


################################################################################
## Fast skew-Gaussian copula covariance helpers used by GeoCovmatrix/GeoKrig. ##
##                                                                            ##
## Write the observed-scale transform as                                      ##
##   g_s(G1,G2) = F_s^{-1[F_Zeta{gamma_eta |G1| + G2}].                      ##
## A bivariate Hermite expansion of g_s gives                                 ##
##   Cov{Y_i,Y_j} = sum_{a+b>=1} f_i(a,b) f_j(a,b) rho^(a+b).                ##
## Truncating this expansion is much cheaper than evaluating the bivariate     ##
## skew-Gaussian density for every covariance entry. It also preserves a       ##
## positive-semidefinite covariance construction: each total-degree term is   ##
## a diagonally scaled Hadamard power of the latent correlation matrix.        ##
##                                                                            ##
## The expensive univariate skew-normal transform is computed only once for   ##
## each copula eta. Location-dependent Beta2/Kumaraswamy2 coefficients are    ##
## tabulated on a compact eta grid; interpolation is done on the coefficient  ##
## Gram matrices, which is equivalent to linearly interpolating the Hermite   ##
## feature vectors and therefore retains the PSD construction.                 ##
################################################################################

.GeoSkewGaussianUniformGH <- function(eta, nquad = 60L) {
  eta <- suppressWarnings(as.numeric(eta))
  if (length(eta) != 1L || !is.finite(eta) || abs(eta) >= 1)
    stop("Skew-Gaussian linear prediction requires eta in (-1,1).", call. = FALSE)
  nquad <- as.integer(nquad)
  if (nquad < 20L)
    stop("Internal skew-Gaussian quadrature requires at least 20 nodes.", call. = FALSE)

  key <- paste0("eta=", format(eta, digits = 17, scientific = TRUE),
                "|q=", nquad, "|v=2")
  cache <- .GeoModels_get_cache("skewgaussian_uniform_gh")
  if (exists(key, envir = cache, inherits = FALSE))
    return(get(key, envir = cache, inherits = FALSE))

  ## g(G1,G2) depends on G1 only through |G1|.  Integrating that coordinate
  ## directly under its half-normal distribution is substantially more
  ## accurate for |eta| close to one than applying ordinary Gauss-Hermite
  ## quadrature across the kink at G1=0.
  half_rule <- .GeoUnitQuadrature(nquad)
  a_nodes <- qnorm((1 + half_rule$nodes) / 2)
  normal_rule <- .GeoNormalQuadrature(nquad)

  gamma_eta <- eta / sqrt(1 - eta^2)
  omega_eta <- sqrt(1 + gamma_eta^2)
  z <- outer(gamma_eta * a_nodes, normal_rule$nodes, `+`)
  u <- sn::psn(as.numeric(z), xi = 0, omega = omega_eta, alpha = gamma_eta)
  u <- matrix(pmin(1 - 1e-12, pmax(1e-12, u)), nrow = nquad, ncol = nquad)

  out <- list(
    half_nodes = a_nodes, half_weights = half_rule$weights,
    normal_nodes = normal_rule$nodes, normal_weights = normal_rule$weights,
    u = u, eta = eta, nquad = nquad
  )
  assign(key, out, envir = cache)
  out
}

.GeoHermiteBasis0 <- function(x, M) {
  M <- as.integer(M)
  cbind(`0` = rep(1, length(x)), .GeoHermiteBasis(x, M))
}

.GeoSkewGaussianFeatureGrams <- function(model_type, params, eta_copula,
                                          eta_grid = 0,
                                          degree = 36L, nquad = 60L) {
  model_type <- as.integer(model_type)
  degree <- as.integer(degree); nquad <- as.integer(nquad)
  if (degree < 4L || nquad <= degree)
    stop("Internal skew-Gaussian Hermite settings are inconsistent.", call. = FALSE)
  eta_grid <- as.numeric(eta_grid)
  if (!length(eta_grid) || any(!is.finite(eta_grid)))
    stop("Invalid skew-Gaussian marginal-location grid.", call. = FALSE)

  gh <- .GeoSkewGaussianUniformGH(eta_copula, nquad = nquad)
  B1 <- .GeoHermiteBasis0(gh$half_nodes, degree)
  ## Averaging over the unobserved sign of G1 eliminates all odd Hermite
  ## orders exactly because the response transform depends on |G1| only.
  if (degree >= 1L) B1[, seq.int(2L, degree + 1L, by = 2L)] <- 0
  B2 <- .GeoHermiteBasis0(gh$normal_nodes, degree)
  WB1 <- B1 * gh$half_weights
  WB2 <- B2 * gh$normal_weights
  sfac <- sqrt(factorial(0:degree))
  normmat <- outer(sfac, sfac)
  degmat <- outer(0:degree, 0:degree, `+`)
  use <- degmat > 0 & degmat <= degree

  nms <- names(params)
  keep <- if (is.null(nms)) seq_along(params) else !grepl("^mean[0-9]*$", nms)
  ppkey <- params[keep]
  if (!is.null(names(ppkey))) ppkey <- ppkey[order(names(ppkey))]
  key <- paste(
    "m", model_type,
    paste(names(ppkey),
          format(as.numeric(unlist(ppkey, use.names = FALSE)), digits = 17, scientific = TRUE),
          sep = "=", collapse = ";"),
    "eta", format(eta_copula, digits = 17, scientific = TRUE),
    "loc", paste(format(eta_grid, digits = 15, scientific = TRUE), collapse = ","),
    "d", degree, "q", nquad, "v1", sep = "|"
  )
  cache <- .GeoModels_get_cache("skewgaussian_feature_grams")
  if (exists(key, envir = cache, inherits = FALSE))
    return(get(key, envir = cache, inherits = FALSE))

  ng <- length(eta_grid)
  feature_by_degree <- vector("list", degree)
  for (k in seq_len(degree)) {
    idx <- which(degmat == k & use, arr.ind = TRUE)
    feature_by_degree[[k]] <- matrix(0, nrow = ng, ncol = nrow(idx))
    attr(feature_by_degree[[k]], "index") <- idx
  }

  for (ii in seq_len(ng)) {
    pp <- params
    pp["mean"] <- eta_grid[ii]
    qv <- quantile_disp(as.numeric(gh$u), model_type, pp)
    if (any(!is.finite(qv)))
      stop("Non-finite marginal quantiles in skew-Gaussian covariance.", call. = FALSE)
    Q <- matrix(qv, nrow = nquad, ncol = nquad)
    mom <- .GeoCopulaMarginalMoments(model_type, params, eta_grid[ii])
    Q <- Q - as.numeric(mom$mean)[1L]

    ## E[g H_a(G1) H_b(G2)] / sqrt(a! b!) are orthonormal Hermite
    ## feature coefficients.  Only total degree <= degree is retained.
    A <- crossprod(WB1, Q %*% WB2) / normmat
    A[!use] <- 0

    ## Calibrate the truncated feature norm to the exact marginal variance.
    ## This makes the rho=1 limit exact for equal margins while preserving the
    ## PSD feature representation (a common positive scale is applied at each
    ## marginal location).
    vtr <- sum(A^2)
    vex <- as.numeric(mom$variance)[1L]
    if (!is.finite(vtr) || vtr <= 0 || !is.finite(vex) || vex <= 0)
      stop("Invalid marginal variance in skew-Gaussian covariance.", call. = FALSE)
    A <- A * sqrt(vex / vtr)

    for (k in seq_len(degree)) {
      idx <- attr(feature_by_degree[[k]], "index")
      feature_by_degree[[k]][ii, ] <- A[idx]
    }
  }

  grams <- lapply(feature_by_degree, function(F) tcrossprod(F))
  out <- list(eta = eta_grid, grams = grams, degree = degree,
              eta_copula = eta_copula, nquad = nquad)
  assign(key, out, envir = cache)
  out
}

.GeoInterp2Weights <- function(g, z1, z2) {
  g <- as.numeric(g)
  z1 <- pmax(g[1L], pmin(g[length(g)], as.numeric(z1)))
  z2 <- pmax(g[1L], pmin(g[length(g)], as.numeric(z2)))
  if (length(g) < 2L)
    stop("Internal interpolation grid must contain at least two points.", call. = FALSE)
  i1 <- pmin(findInterval(z1, g, all.inside = TRUE), length(g) - 1L)
  i2 <- pmin(findInterval(z2, g, all.inside = TRUE), length(g) - 1L)
  a1 <- (z1 - g[i1]) / pmax(g[i1 + 1L] - g[i1], .Machine$double.eps)
  a2 <- (z2 - g[i2]) / pmax(g[i2 + 1L] - g[i2], .Machine$double.eps)
  list(i1 = as.integer(i1), i2 = as.integer(i2), a1 = a1, a2 = a2)
}

.GeoInterp2Grid <- function(g, tab, z1, z2) {
  z1 <- pmax(g[1L], pmin(g[length(g)], as.numeric(z1)))
  z2 <- pmax(g[1L], pmin(g[length(g)], as.numeric(z2)))
  if (length(g) == 1L) return(rep(tab[1L, 1L], max(length(z1), length(z2))))
  w <- .GeoInterp2Weights(g, z1, z2)
  i1 <- w$i1; i2 <- w$i2; a1 <- w$a1; a2 <- w$a2
  c00 <- tab[cbind(i1, i2)]
  c10 <- tab[cbind(i1 + 1L, i2)]
  c01 <- tab[cbind(i1, i2 + 1L)]
  c11 <- tab[cbind(i1 + 1L, i2 + 1L)]
  (1 - a1) * (1 - a2) * c00 +
    a1 * (1 - a2) * c10 +
    (1 - a1) * a2 * c01 +
    a1 * a2 * c11
}

.GeoSkewGaussianCovLookupTable <- function(model_type, params, eta_copula,
                                            degree = 36L, nquad = 60L,
                                            neta = 51L, eta_range = NULL) {
  model_type <- as.integer(model_type)
  additive <- c(1L, 12L, 20L, 25L, 34L, 40L, 59L)
  multiplicative <- c(21L, 22L, 26L)
  fixed_bounded <- c(28L, 33L)
  location_dependent <- c(50L, 42L, .GeoGaussianCopulaDiscreteModels)

  nms <- names(params)
  keep <- if (is.null(nms)) seq_along(params) else !grepl("^mean[0-9]*$", nms)
  ppkey <- params[keep]
  if (!is.null(names(ppkey))) ppkey <- ppkey[order(names(ppkey))]
  base_key <- paste(
    "m", model_type,
    paste(names(ppkey),
          format(as.numeric(unlist(ppkey, use.names = FALSE)), digits = 17, scientific = TRUE),
          sep = "=", collapse = ";"),
    "eta", format(eta_copula, digits = 17, scientific = TRUE),
    "d", degree, "q", nquad, "e", neta, "v2", sep = "|"
  )
  cache <- .GeoModels_get_cache("skewgaussian_cov_lookup")

  if (model_type %in% c(additive, multiplicative, fixed_bounded)) {
    key <- paste0(base_key, "|fixed")
    if (exists(key, envir = cache, inherits = FALSE))
      return(get(key, envir = cache, inherits = FALSE))
    ff <- .GeoSkewGaussianFeatureGrams(
      model_type, params, eta_copula, eta_grid = 0,
      degree = degree, nquad = nquad
    )
    coef <- vapply(ff$grams, function(G) G[1L, 1L], numeric(1L))
    out <- list(
      kind = if (model_type %in% multiplicative) "multiplicative" else "fixed",
      coef = coef, degree = seq_along(coef)
    )
    assign(key, out, envir = cache)
    return(out)
  }

  if (model_type %in% location_dependent) {
    er <- suppressWarnings(as.numeric(eta_range))
    er <- er[is.finite(er)]
    if (!length(er)) er <- 0
    ## Quantize the requested band to two-unit boundaries.  This makes the
    ## same table reusable by many GeoKrigloc neighborhoods while still
    ## avoiding the numerically degenerate extreme-location grid that a fixed
    ## global range would create for Beta2/Kumaraswamy2.
    elo <- 2 * floor((min(er) - 0.5) / 2)
    ehi <- 2 * ceiling((max(er) + 0.5) / 2)
    if (!is.finite(elo) || !is.finite(ehi))
      stop("Invalid skew-Gaussian marginal-location range.", call. = FALSE)
    if (ehi <= elo) ehi <- elo + 4
    if ((ehi - elo) < 4) {
      mid <- 0.5 * (elo + ehi)
      elo <- 2 * floor((mid - 1) / 2)
      ehi <- elo + 4
    }
    key <- paste0(base_key, "|eta=", format(elo, scientific = FALSE), ":",
                  format(ehi, scientific = FALSE))
    if (exists(key, envir = cache, inherits = FALSE))
      return(get(key, envir = cache, inherits = FALSE))

    neta_eff <- max(as.integer(neta), as.integer(4L * ceiling(ehi - elo) + 1L))
    eg <- seq(elo, ehi, length.out = neta_eff)
    ff <- .GeoSkewGaussianFeatureGrams(
      model_type, params, eta_copula, eta_grid = eg,
      degree = degree, nquad = nquad
    )
    out <- list(kind = "location", eta = eg, grams = ff$grams,
                gram_flat = unlist(ff$grams, use.names = FALSE),
                degree = seq_along(ff$grams), eta_range = c(elo, ehi))
    assign(key, out, envir = cache)
    return(out)
  }

  stop("Unsupported margin in skew-Gaussian covariance.", call. = FALSE)
}


# Reference evaluator retained for numerical regression tests.  It mirrors the
# historical skew-Gaussian covariance loop exactly at R level.
.GeoCovSeriesForwardR <- function(rho, coefficient, scale = NULL) {
  rho <- as.numeric(rho)
  n <- length(rho)
  if (!n) return(numeric(0))
  coefficient <- as.numeric(coefficient)
  out <- numeric(n)
  rp <- rho
  scaled <- !is.null(scale)
  if (scaled) {
    scale <- as.numeric(scale)
    if (length(scale) == 1L) scale <- rep(scale, n)
    if (length(scale) != n)
      stop("Internal covariance-series scale has incompatible length.", call. = FALSE)
  }
  for (k in seq_along(coefficient)) {
    ck <- rep(coefficient[k], n)
    if (scaled) ck <- ck * scale
    out <- out + ck * rp
    rp <- rp * rho
  }
  out[abs(rho) < 1e-15] <- 0
  out
}

# Native evaluator for fixed-coefficient covariance series.  One native pass
# replaces degree full-vector R passes and avoids the temporary ck/rp vectors
# created by the historical implementation.  The C loop uses the same
# low-to-high power accumulation order as .GeoCovSeriesForwardR().
.GeoCovSeriesForwardC <- function(rho, coefficient, scale = NULL) {
  rho <- as.numeric(rho)
  n <- length(rho)
  if (!n) return(numeric(0))
  coefficient <- as.numeric(coefficient)
  degree <- length(coefficient)
  if (!degree) return(numeric(n))

  scaled <- !is.null(scale)
  if (scaled) {
    scale <- as.numeric(scale)
    if (length(scale) == 1L) scale <- rep(scale, n)
    if (length(scale) != n)
      stop("Internal covariance-series scale has incompatible length.", call. = FALSE)
  } else {
    ## C64 still needs a valid address even though the native routine does not
    ## dereference scale when use_scale == 0.
    scale <- 1
  }

  dotCall64::.C64(
    "GeoCovSeriesForward",
    SIGNATURE = c("double", "integer", "double", "integer",
                  "double", "integer", "double"),
    rho, as.integer(n), coefficient, as.integer(degree),
    scale, as.integer(scaled),
    out = dotCall64::vector_dc("double", n),
    INTENT = c("r", "r", "r", "r", "r", "r", "w"),
    NAOK = TRUE, PACKAGE = "GeoModels"
  )$out
}

.GeoCovSeriesForward <- function(rho, coefficient, scale = NULL) {
  if (isTRUE(getOption("GeoModels.covseries.native", TRUE)))
    .GeoCovSeriesForwardC(rho, coefficient, scale = scale)
  else
    .GeoCovSeriesForwardR(rho, coefficient, scale = scale)
}

# Reference evaluator for location-dependent covariance coefficients.  This is
# the historical R path retained for regression tests and for disabling the
# native kernels through GeoModels.covseries.native.
.GeoCovSeriesInterp2ForwardR <- function(rho, eta_grid, grams, eta1, eta2) {
  rho <- as.numeric(rho)
  n <- length(rho)
  if (!n) return(numeric(0))
  eta1 <- .GeoCopulaEtaVector(eta1, n)
  eta2 <- .GeoCopulaEtaVector(eta2, n)
  out <- numeric(n)
  rp <- rho
  for (k in seq_along(grams)) {
    ck <- .GeoInterp2Grid(eta_grid, grams[[k]], eta1, eta2)
    out <- out + ck * rp
    rp <- rp * rho
  }
  out[abs(rho) < 1e-15] <- 0
  out
}

# Native evaluator for Beta2/Kumaraswamy2 skew-Gaussian covariance tables.
# Interpolation indices and weights are computed once in R, rather than once
# per Hermite degree, and the degree loop plus bilinear interpolation are then
# evaluated in one native pass over all covariance pairs.
.GeoCovSeriesInterp2ForwardC <- function(rho, eta_grid, grams, eta1, eta2,
                                             gram_flat = NULL) {
  rho <- as.numeric(rho)
  n <- length(rho)
  if (!n) return(numeric(0))
  eta_grid <- as.numeric(eta_grid)
  ng <- length(eta_grid)
  if (ng < 2L)
    stop("Internal covariance interpolation grid is too short.", call. = FALSE)
  degree <- length(grams)
  if (!degree) return(numeric(n))
  eta1 <- .GeoCopulaEtaVector(eta1, n)
  eta2 <- .GeoCopulaEtaVector(eta2, n)
  w <- .GeoInterp2Weights(eta_grid, eta1, eta2)

  if (is.null(gram_flat)) {
    gram_flat <- unlist(lapply(grams, function(z) {
      z <- as.matrix(z)
      if (!identical(dim(z), c(ng, ng)))
        stop("Internal covariance interpolation table has inconsistent dimensions.",
             call. = FALSE)
      as.numeric(z)
    }), use.names = FALSE)
  } else {
    gram_flat <- as.numeric(gram_flat)
    if (length(gram_flat) != degree * ng * ng)
      stop("Internal flattened covariance table has incompatible length.",
           call. = FALSE)
  }

  dotCall64::.C64(
    "GeoCovSeriesInterp2Forward",
    SIGNATURE = c("double", "integer", "double", "integer", "integer",
                  "integer", "integer", "double", "double", "double"),
    rho, as.integer(n), gram_flat, as.integer(ng), as.integer(degree),
    w$i1, w$i2, w$a1, w$a2,
    out = dotCall64::vector_dc("double", n),
    INTENT = c("r", "r", "r", "r", "r", "r", "r", "r", "r", "w"),
    NAOK = TRUE, PACKAGE = "GeoModels"
  )$out
}

.GeoCovSeriesInterp2Forward <- function(rho, eta_grid, grams, eta1, eta2,
                                            gram_flat = NULL) {
  if (isTRUE(getOption("GeoModels.covseries.native", TRUE)))
    .GeoCovSeriesInterp2ForwardC(
      rho, eta_grid, grams, eta1, eta2, gram_flat = gram_flat
    )
  else
    .GeoCovSeriesInterp2ForwardR(rho, eta_grid, grams, eta1, eta2)
}

.GeoSkewGaussianCopulaCrossCovFast <- function(rho, model_type, params,
                                                eta1 = 0, eta2 = eta1,
                                                degree = 36L, nquad = 60L,
                                                neta = 51L, eta_range = NULL) {
  model_type <- as.integer(model_type)
  if (!(model_type %in% .GeoCopulaPredictionModels))
    stop("Skew-Gaussian linear prediction is not implemented for model id ",
         model_type, ".", call. = FALSE)
  if (model_type %in% .GeoGaussianCopulaDiscreteModels) {
    degree <- max(as.integer(degree), 48L)
    nquad <- max(as.integer(nquad), 80L)
  }
  rho <- as.numeric(rho); n <- length(rho)
  if (!n) return(numeric(0))
  eta1 <- .GeoCopulaEtaVector(eta1, n)
  eta2 <- .GeoCopulaEtaVector(eta2, n)
  if (any(!is.finite(rho)) || any(abs(rho) > 1 + 1e-10))
    stop("Invalid latent correlations in skew-Gaussian covariance.", call. = FALSE)
  rho <- pmax(-1, pmin(1, rho))

  eta_copula <- .GeoCopulaParamScalar(params, "nu")
  if (abs(eta_copula) >= 1)
    stop("Skew-Gaussian linear prediction requires nu (eta) in (-1,1).", call. = FALSE)

  ## eta=0 is exactly the Gaussian copula.
  if (abs(eta_copula) < 1e-10) {
    return(.GeoGaussianCopulaCrossCov(
      rho, model_type, params, eta1 = eta1, eta2 = eta2
    ))
  }

  lookup_eta_range <- if (is.null(eta_range)) {
    range(c(eta1, eta2), finite = TRUE)
  } else {
    eta_range
  }
  lookup <- .GeoSkewGaussianCovLookupTable(
    model_type, params, eta_copula,
    degree = degree, nquad = nquad, neta = neta,
    eta_range = lookup_eta_range
  )

  ## Fixed and multiplicative margins have one coefficient vector shared by
  ## all pairs.  Location-dependent Beta2/Kumaraswamy2 margins use the same
  ## cached Hermite Gram tables, with bilinear interpolation and series
  ## accumulation fused in a native pass.
  if (!identical(lookup$kind, "location")) {
    scale <- if (identical(lookup$kind, "multiplicative"))
      exp(eta1 + eta2) else NULL
    return(.GeoCovSeriesForward(rho, lookup$coef, scale = scale))
  }

  .GeoCovSeriesInterp2Forward(
    rho, lookup$eta, lookup$grams, eta1 = eta1, eta2 = eta2,
    gram_flat = lookup$gram_flat
  )
}


## Common observed-scale covariance dispatcher for the three copula families
## supported by GeoModels, including the audited count margins used by linear
## prediction.  Keeping the dispatch in one place
## guarantees that GeoCovmatrix(), GeoCovariogram(), and GeoKrig() use the
## same numerical covariance engine.
.GeoCopulaCrossCovFast <- function(copula, rho, model_type, params,
                                    eta1 = 0, eta2 = eta1, eta_range = NULL) {
  if (is.null(copula) || length(copula) != 1L || is.na(copula))
    stop("A single copula name is required for copula covariance.", call. = FALSE)
  copula <- as.character(copula)
  if (!(copula %in% c("Gaussian", "Clayton", "SkewGaussian")))
    stop("Copula covariance is implemented for 'Gaussian', 'Clayton', and 'SkewGaussian'.",
         call. = FALSE)

  model_type <- as.integer(model_type)
  if (!(model_type %in% .GeoCopulaPredictionModels))
    stop("Observed-scale covariance for copula='", copula,
         "' is not implemented for model id ", model_type, ".", call. = FALSE)

  if (identical(copula, "Gaussian")) {
    return(.GeoGaussianCopulaCrossCov(
      rho, model_type, params, eta1 = eta1, eta2 = eta2
    ))
  }
  if (identical(copula, "Clayton")) {
    return(.GeoClaytonCopulaCrossCovFast(
      rho, model_type, params, eta1 = eta1, eta2 = eta2,
      eta_range = eta_range
    ))
  }
  .GeoSkewGaussianCopulaCrossCovFast(
    rho, model_type, params, eta1 = eta1, eta2 = eta2,
    eta_range = eta_range
  )
}


################################################################################
## Memory-bounded dense assembly for observed-scale copula covariance matrices.##
##                                                                            ##
## The native correlation routines return the strict lower triangle in         ##
## column-major order.  The former implementation materialized, at the same    ##
## time, a full pair-index matrix, two pairwise eta vectors, the complete       ##
## transformed covariance vector, and a lower.tri() logical matrix.  For large ##
## n those temporary objects can be comparable in size to the covariance       ##
## matrix itself.  Here the same ordered correlation vector is consumed in     ##
## bounded contiguous column blocks and written directly to the final matrix.  ##
################################################################################
.GeoDenseCopulaCovarianceBlocked <- function(raw_corr, dime, copula, model_type,
                                              params, eta_obs, variances,
                                              nugget_scale,
                                              block_pairs = 1000000L) {
  dime <- as.integer(dime)
  if (length(dime) != 1L || is.na(dime) || dime < 1L)
    stop("Invalid covariance dimension.", call. = FALSE)

  expected <- as.double(dime) * as.double(dime - 1L) / 2
  if (length(raw_corr) != expected)
    stop("Internal copula covariance ordering mismatch.", call. = FALSE)

  variances <- as.numeric(variances)
  if (length(variances) == 1L) variances <- rep(variances, dime)
  if (length(variances) != dime || any(!is.finite(variances)))
    stop("Invalid marginal variances in copula covariance assembly.", call. = FALSE)

  eta_obs <- as.numeric(eta_obs)
  if (length(eta_obs) == 1L) eta_obs <- rep(eta_obs, dime)
  if (length(eta_obs) != dime || any(!is.finite(eta_obs)))
    stop("Invalid marginal locations in copula covariance assembly.", call. = FALSE)

  nugget_scale <- as.numeric(nugget_scale)[1L]
  if (!is.finite(nugget_scale))
    stop("Invalid nugget scaling in copula covariance assembly.", call. = FALSE)

  block_pairs <- suppressWarnings(as.double(block_pairs)[1L])
  if (!is.finite(block_pairs) || block_pairs < 1)
    block_pairs <- 1000000
  block_pairs <- floor(block_pairs)

  varcov <- matrix(0, dime, dime)
  diag(varcov) <- variances
  if (dime <= 1L) return(varcov)

  ## A single global location range makes blockwise Clayton/SkewGaussian
  ## interpolation numerically identical to a one-shot covariance transform.
  eta_range <- range(eta_obs, finite = TRUE)
  location_sensitive <- as.integer(model_type) %in% c(21L, 22L, 26L, 42L, 50L,
                                                       .GeoGaussianCopulaDiscreteModels)

  pos <- 1
  col_start <- 1L
  while (col_start < dime) {
    first_len <- dime - col_start
    ncols <- max(1, floor(block_pairs / first_len))
    ncols <- min(ncols, dime - col_start)
    cols <- seq.int(col_start, length.out = ncols)
    lens <- dime - cols
    np <- sum(lens)

    idx <- seq.int(pos, length.out = np)
    rho <- as.numeric(raw_corr[idx]) * nugget_scale

    if (location_sensitive) {
      ## Generate only block-sized eta vectors; no global n(n-1)/2 pair index
      ## matrix (or list of per-column eta vectors) is ever materialized.
      eta1 <- numeric(np)
      eta2 <- numeric(np)
      eta_pos <- 1L
      for (ii in seq_along(cols)) {
        j <- cols[ii]
        nr <- lens[ii]
        take <- seq.int(eta_pos, length.out = nr)
        eta1[take] <- eta_obs[seq.int(j + 1L, dime)]
        eta2[take] <- eta_obs[j]
        eta_pos <- eta_pos + nr
      }
    } else {
      eta1 <- 0
      eta2 <- 0
    }

    cv <- .GeoCopulaCrossCovFast(
      copula, rho, model_type, params,
      eta1 = eta1, eta2 = eta2, eta_range = eta_range
    )
    if (length(cv) != np || any(!is.finite(cv)))
      stop("Invalid block in copula covariance assembly.", call. = FALSE)

    ## Write each strict-lower column and its transpose directly.  This avoids
    ## lower.tri() logical matrices and avoids a second full covariance vector.
    local_pos <- 1L
    for (ii in seq_along(cols)) {
      j <- cols[ii]
      nr <- lens[ii]
      rows <- seq.int(j + 1L, dime)
      take <- seq.int(local_pos, length.out = nr)
      vals <- cv[take]
      varcov[rows, j] <- vals
      varcov[j, rows] <- vals
      local_pos <- local_pos + nr
    }

    pos <- pos + np
    col_start <- col_start + ncols
  }

  if ((pos - 1) != length(raw_corr))
    stop("Internal copula covariance block accounting mismatch.", call. = FALSE)
  varcov
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
  
raw_corr <- cr$corr
cr <- NULL
nugget_scale <- 1 - as.numeric(nuisance['nugget'])
if (!length(nugget_scale) || !is.finite(nugget_scale[1L])) nugget_scale <- 1
nugget_scale <- nugget_scale[1L]
if(model == 11L) {
  ## GeoKrig expands a scalar Binomial trial count to one value per
  ## observation before covariance construction.  The Gaussian-copula
  ## Hermite covariance below currently uses one common marginal size, so
  ## collapse an equal vector back to its scalar value.  Do not silently
  ## use the first value when trial counts are genuinely location-specific.
  n_bin <- as.numeric(n)
  if (!length(n_bin) || any(!is.finite(n_bin)))
    stop("Cmatrix_copula: Binomial n must contain finite values.", call. = FALSE)
  if (length(n_bin) > 1L && any(n_bin != n_bin[1L]))
    stop("Cmatrix_copula: copula-based Binomial covariance currently requires a common n across locations; location-specific n is not yet implemented.", call. = FALSE)
  nuisance["n"] <- n_bin[1L]
}
if(model == 16L) {
  ## Negative-Binomial n is a common field parameter by construction.
  nuisance["n"] <- as.numeric(n)[1L]
}

if (!(copula %in% c("Gaussian", "Clayton", "SkewGaussian"))) {
  stop("Cmatrix_copula() currently supports Gaussian, Clayton-like, and skew-Gaussian copulas.")
}

common_copula_model <- as.integer(model) %in% .GeoGaussianCopulaContinuousModels ||
  (copula %in% c("Clayton", "SkewGaussian") &&
   as.integer(model) %in% .GeoGaussianCopulaDiscreteModels)
if (!bivariate && common_copula_model) {
  ## ML contains the location predictor X beta for models whose covariance
  ## changes with the marginal location (positive-scale margins, Beta2 and
  ## Kumaraswamy2).  For additive/fixed-shape margins it is harmless to use
  ## the intercept because their covariance does not depend on location.
  eta0 <- suppressWarnings(as.numeric(nuisance["mean"]))
  if (!length(eta0) || !is.finite(eta0[1L])) eta0 <- 0
  eta_obs <- if (is.numeric(ML) && length(ML) == dime && all(is.finite(ML))) {
    as.numeric(ML)
  } else {
    rep(eta0[1L], dime)
  }

  moments <- .GeoCopulaMarginalMoments(model, nuisance, eta_obs)
  vv <- as.numeric(moments$variance)

  if (type == "Standard") {
    varcov <- .GeoDenseCopulaCovarianceBlocked(
      raw_corr = raw_corr, dime = dime, copula = copula, model_type = model,
      params = nuisance, eta_obs = eta_obs, variances = vv,
      nugget_scale = nugget_scale
    )
  } else if (type == "Tapering") {
    corr <- as.numeric(raw_corr) * nugget_scale
    nr <- diff(setup$ia)
    row_id <- rep(seq_len(dime), nr)
    col_id <- as.integer(setup$ja)
    if (length(row_id) != length(corr) || length(col_id) != length(corr) ||
        any(col_id < 1L | col_id > dime))
      stop("Internal sparse copula covariance ordering mismatch.", call. = FALSE)
    cova <- .GeoCopulaCrossCovFast(
      copula, corr, model, nuisance,
      eta1 = eta_obs[row_id], eta2 = eta_obs[col_id],
      eta_range = range(eta_obs, finite = TRUE)
    )
    varcov <- new("spam", entries = cova, colindices = setup$ja,
                  rowpointers = setup$ia, dimension = as.integer(rep(dime, 2)))
    diag(varcov) <- vv
  }
} else {
  ## Legacy fallback for older Gaussian-copula models not handled by the
  ## common observed-scale covariance engine above.
  if (copula %in% c("Clayton", "SkewGaussian"))
    stop(sprintf("%s covariance is not implemented for model id %s.", copula, model), call. = FALSE)
  corr <- as.numeric(raw_corr) * nugget_scale
  cova <- gaussian_copula_cov(corr, model, nuisance)
  vv <- variance_disp(model, nuisance)
  if(!bivariate) {
    if(type=="Standard") {
      varcov <- diag(dime)
      varcov[lower.tri(varcov)] <- cova
      varcov <- t(varcov)
      varcov[lower.tri(varcov)] <- cova
      diag(varcov) <- vv
    }
    if(type=="Tapering") {
      varcov <- new("spam", entries=cova, colindices=setup$ja,
                    rowpointers=setup$ia, dimension=as.integer(rep(dime,2)))
      diag(varcov) <- vv
    }
  }
}

return(varcov)
}
######################################################################################################