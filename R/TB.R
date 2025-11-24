#########################################################
#########################################################
############ functions for TB############################
#########################################################
#########################################################



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
      G <- matrix(rgamma(P*L*N, nu0, scale = 1),P*L*N,d)
      u <- matrix(rnorm(P*L*N*d), P*L*N, d)/sqrt(G*2)/a0/(2*pi)
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
                         simu1=dotCall64::vector_dc("double",length(simu11)),L,
         INTENT =    c(rep("r",20),
                       "r", "rw","r"),
             PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$simu1
  simu =  matrix(result,n,P)
  return(simu)
}

#######################################################
######### univariate case #############################
#######################################################

safe_besselJ <- function(x, nu) {
  out <- numeric(length(x))
  x <- as.numeric(x)
  nu <- as.numeric(nu)

  tiny_thr <- 1e-6
  big_thr  <- 200

  small <- is.finite(x) & x <= tiny_thr & x >= 0
  if (any(small)) {
    z <- x[small]
    out[small] <- (pmax(z, .Machine$double.eps)/2)^nu / gamma(nu + 1)
  }

  mid <- is.finite(x) & x > tiny_thr & x <= big_thr
  if (any(mid)) {
    out[mid] <- base::besselJ(x[mid], nu)
  }

  big <- is.finite(x) & x > big_thr
  if (any(big)) {
    z <- x[big]
    out[big] <- sqrt(2/(pi*z)) * cos(z - nu*pi/2 - pi/4)
  }

  out[!is.finite(out)] <- 0
  out
}


sphere_area <- function(d) {
  2 * pi^(d / 2) / gamma(d / 2)
}

#########################################################
# beta  MIXTURE METHOD 
#########################################################
C_H <- function(kappa, d) {
  num <- gamma((d + 1) / 2 + kappa) * gamma(1 + kappa) * gamma(d / 2 + 1 + 2 * kappa)
  den <- 2^d * pi^(d / 2) * gamma(0.5 + kappa) * gamma(d / 2 + 1 + kappa) * gamma(d + 1 + 2 * kappa)
  num / den
}

f_T <- function(t, kappa, d, log_scale = FALSE) {
  t <- as.numeric(t)
  out <- numeric(length(t))
  out[t <= 0] <- if (log_scale) -Inf else 0
  if (all(t <= 0)) return(out)

  delta <- (d + 1) / 2 + kappa
  S_d <- sphere_area(d)
  CH  <- C_H(kappa, d)
  const <- S_d * CH * gamma(delta + 0.5)^2 * 4^(2 * delta - 1)

  pos <- t > 0
  x <- t[pos] / 2

  #bval <- besselJ(x, delta - 0.5)
  bval <-safe_besselJ(x, delta - 0.5)
  #big <- x > 200
  #if (any(big)) {
  #  nu <- delta - 0.5
  #  bval[big] <- sqrt(2 / (pi * x[big])) * cos(x[big] - nu * pi / 2 - pi / 4)
  #}

  val <- const * pmax(t[pos], 1e-12)^(d - 2 * delta) * (bval^2)

  if (log_scale) {
    out[pos] <- log(pmax(val, .Machine$double.eps))
  } else {
    out[pos] <- val
  }
  out
}

sample_fT_final <- function(n, kappa, d, T0 = NULL, safety = 1.15, batch = 100000, 
                             max_iter = 1000, optimize_T0 = TRUE) {
  delta <- (d + 1) / 2 + kappa
  S_d   <- sphere_area(d)
  CH    <- C_H(kappa, d)
  B_small <- S_d * CH
  C0 <- S_d * CH * gamma(delta + 0.5)^2 * 4^(2 * delta - 1)
  A_tail <- (4 / pi) * C0
  alpha  <- 2 * kappa + 2
  shape  <- alpha - 1
  
  M1_from_T0 <- function(T0) safety * (B_small * T0^d) / d
  M2_from_T0 <- function(T0) safety * (A_tail / (shape * T0^shape))
  
  if (is.null(T0) && optimize_T0) {
    objective <- function(T0) M1_from_T0(T0) + M2_from_T0(T0)
    T0 <- optimize(objective, interval = c(0.5, 50), maximum = FALSE)$minimum
  } else if (is.null(T0)) {
    t_mode <- tryCatch(
      optimize(function(t) f_T(t, kappa, d), c(1e-3, 40), maximum = TRUE)$maximum,
      error = function(e) 2
    )
    T0 <- max(6, 2 * t_mode)
  }
  
  M1 <- M1_from_T0(T0)
  M2 <- M2_from_T0(T0)
  w1 <- M1 / (M1 + M2)
  w2 <- 1 - w1
  
  g1_rnd <- function(n) T0 * runif(n)^(1 / d)
  g1_pdf <- function(t) ifelse(t >= 0 & t <= T0, d * t^(d - 1) / (T0^d), 0)
  g2_rnd <- function(n) T0 * (1 - runif(n))^(-1 / shape)
  g2_pdf <- function(t) ifelse(t > T0, shape * (T0^shape) / (t^(shape + 1)), 0)
  
  samples <- numeric(0)
  n_total_proposed <- 0
  n_iter <- 0
  
  while (length(samples) < n && n_iter < max_iter) {
    n_iter <- n_iter + 1
    n_batch <- min(batch, (n - length(samples)) * 5)
    
    use_body <- runif(n_batch) < w1
    n_body <- sum(use_body)
    n_tail <- n_batch - n_body
    
    proposals <- numeric(n_batch)
    M_vals <- numeric(n_batch)
    g_vals <- numeric(n_batch)
    
    if (n_body > 0) {
      proposals[use_body] <- g1_rnd(n_body)
      M_vals[use_body] <- M1
      g_vals[use_body] <- g1_pdf(proposals[use_body])
    }
    
    if (n_tail > 0) {
      proposals[!use_body] <- g2_rnd(n_tail)
      M_vals[!use_body] <- M2
      g_vals[!use_body] <- g2_pdf(proposals[!use_body])
    }
    
    f_vals <- f_T(proposals, kappa, d)
    u <- runif(n_batch)
    accept <- u < f_vals / (M_vals * g_vals)
    
    samples <- c(samples, proposals[accept])
    n_total_proposed <- n_total_proposed + n_batch
  }
  
  samples[1:min(n, length(samples))]
}

#####################################################################
simulate_GH_beta <- function(L, d, kappa, mu, l, a) {
  # Caso limite mu = 1: U -> 1 quasi sicuramente
  if (abs(mu - 1) < 1e-8) {
    shapeV1 <- d / 2 + 2 * kappa + 1
    shapeV2 <- mu / 2 - d / 2 - 0.5 - kappa + l
    V <- if (shapeV2 <= 1e-10) rep(1, L) else rbeta(L, shapeV1, shapeV2)
    U <- rep(1, L)   # limite U → 1
  } else if (mu > 1) {
    shapeU1 <- kappa + 1
    shapeU2 <- mu / 2 - 0.5
    shapeV1 <- d / 2 + 2 * kappa + 1
    shapeV2 <- mu / 2 - d / 2 - 0.5 - kappa + l
    U <- rbeta(L, shapeU1, shapeU2)
    V <- if (shapeV2 <= 1e-10) rep(1, L) else rbeta(L, shapeV1, shapeV2)
  } else {
    stop("mu must be >= 1 for the Beta mixture representation.")
  }

  b <- a * sqrt(U * V)
  Tvals <- sample_fT_final(L, kappa, d)
  R <- Tvals / b

  Theta <- matrix(rnorm(L * d), L, d)
  Theta <- Theta / sqrt(rowSums(Theta^2))
  R * Theta
}
###################################################################

#########################################################
# GASPER MIXTURE METHOD 
#########################################################

log_C_H_n <- function(eta, n, d) {
  kappa <- eta + n - d/2
  lgamma(d/2) + lgamma((d + 1)/2 + kappa) +
    lgamma(1 + kappa - n) + lgamma(d/2 + 1 + 2*kappa - n) -
    ((d + 2*n) * log(2) + (d/2) * log(pi) +
       lgamma(d/2 + n) + lgamma(0.5 + kappa - n) +
       lgamma(d/2 + 1 + kappa) + lgamma(d + 1 + 2*kappa))
}

C_4F3_terminating <- function(n, eta, delta, beta, gamma) {
  if (n == 0) return(1.0)
  tk <- 1.0
  s <- 1.0
  for (k in 0:(n - 1)) {
    num <- (-n + k) * (n + 2*eta + k) * (eta + 1 + k) * (delta + k)
    den <- (eta + 0.5 + k) * (beta + k) * (gamma + k) * (k + 1)
    tk  <- tk * (num / den)
    s   <- s + tk
  }
  s
}

w_n_stable <- function(nmax, delta, beta, gamma, d) {
  eta <- 0.5 * (beta + gamma - delta - 1.5)
  logL <- lgamma(delta) + lgamma(beta - d/2) + lgamma(gamma - d/2) -
          (d*log(2) + (d/2)*log(pi) + lgamma(delta - d/2) + lgamma(beta) + lgamma(gamma))

  w_log <- numeric(nmax + 1)
  for (n in 0:nmax) {
    Cn_val <- C_4F3_terminating(n, eta, delta, beta, gamma)
    if (!is.finite(Cn_val) || Cn_val < 0) Cn_val <- 0
    logCH <- log_C_H_n(eta, n, d)
    term1 <- 2 * (lgamma(eta + 1) - lgamma(eta + n + 1))
    term2 <- log(2 * n + 2 * eta) + (lgamma(2 * eta + 1 + n) - lgamma(2 * eta + 1)) -
             log(n + 2 * eta) - lgamma(n + 1)
    w_log[n + 1] <- logL + term1 - 4 * n * log(2) + term2 +
                    log(max(Cn_val, .Machine$double.eps)) - logCH
  }
  m <- max(w_log)
  w <- exp(w_log - m)
  w <- pmax(w, 0)
  w / sum(w)
}

fT_n <- function(t, eta, n, d) {
  t <- as.numeric(t)
  Sd <- sphere_area(d)
  CH <- exp(log_C_H_n(eta, n, d))
  const <- Sd * CH * 4^(2*eta + 2*n) * gamma(1 + eta + n)^2
  #const * t^(d - 1 - 2*eta) * besselJ(t / 2, eta + n)^2  
   const * t^(d - 1 - 2*eta) * safe_besselJ(t / 2, eta + n)^2            
}

sample_fT_n_adaptive <- function(N, eta, n, d, safety = 1.2) {
  Sd <- sphere_area(d)
  CH <- exp(log_C_H_n(eta, n, d))
  const <- Sd * CH * 4^(2*eta + 2*n) * gamma(1 + eta + n)^2
  nu <- eta + n

  A_body <- Sd * CH
  A_tail <- const * (4 / pi)
  alpha  <- 2*eta + 1 - d

  M1_from_T0 <- function(T0) safety * (A_body / d) * T0^(d + 2*n)
  M2_from_T0 <- function(T0) safety * (A_tail / (alpha * T0^alpha))
  T0 <- optimize(function(T0) M1_from_T0(T0) + M2_from_T0(T0),
                 interval = c(0.5, 60), maximum = FALSE)$minimum

  M1 <- M1_from_T0(T0)
  M2 <- M2_from_T0(T0)
  w1 <- M1 / (M1 + M2)

  g1_rnd <- function(m) T0 * runif(m)^(1 / d)
  g1_pdf <- function(t) ifelse(t >= 0 & t <= T0, d * t^(d - 1) / (T0^d), 0)
  g2_rnd <- function(m) T0 * (1 - runif(m))^(-1 / alpha)
  g2_pdf <- function(t) ifelse(t > T0, alpha * T0^alpha / (t^(alpha + 1)), 0)

  f_target <- function(t) const * t^(d - 1 - 2*eta) * besselJ(t / 2, nu)^2

  samples <- numeric(0)
  while (length(samples) < N) {
    m <- min(10000, (N - length(samples)) * 5)
    use_body <- runif(m) < w1
    t_prop <- numeric(m)
    gval <- numeric(m)
    Mval <- numeric(m)

    if (any(use_body)) {
      nb <- sum(use_body)
      t_prop[use_body] <- g1_rnd(nb)
      gval[use_body]   <- g1_pdf(t_prop[use_body])
      Mval[use_body]   <- M1
    }
    if (any(!use_body)) {
      nt <- sum(!use_body)
      t_prop[!use_body] <- g2_rnd(nt)
      gval[!use_body]   <- g2_pdf(t_prop[!use_body])
      Mval[!use_body]   <- M2
    }

    fval <- f_target(t_prop)
    acc  <- runif(m) < fval / (Mval * gval)
    samples <- c(samples, t_prop[acc])
  }
  samples[1:N]
}

simulate_GH_gasper <- function(L, d, kappa, mu, l, a, nmax = 15) {
  delta <- (d + 1) / 2 + kappa
  beta  <- delta + mu / 2
  gamma <- delta + mu / 2 + l
  eta   <- 0.5 * (beta + gamma - delta - 1.5)
  
  weights <- w_n_stable(nmax, delta, beta, gamma, d)
  n_indices <- sample(0:nmax, size = L, replace = TRUE, prob = weights)
  
  T_samples <- numeric(L)
  
  for (n in 0:nmax) {
    indices_n <- which(n_indices == n)
    n_samples <- length(indices_n)
    
    if (n_samples == 0) next
    
    T_samples[indices_n] <- sample_fT_n_adaptive(
      N = n_samples, 
      eta = eta, 
      n = n, 
      d = d,
      safety = 1.2
    )
  }
  
  R_samples <- T_samples / a
  
  Theta <- matrix(rnorm(L * d), nrow = L, ncol = d)
  Theta <- Theta / sqrt(rowSums(Theta^2))
  
  Omega <- R_samples * Theta
  
  return(Omega)
}

#########################################################
# MAIN SIMULATION FUNCTION
#########################################################

tbm2d_uni <- function(coord, coordt, param, corrmodel, L, bivariate, parallel, ncores){
  
  N <- 1
  n <- dim(coord)[1]
  d <- dim(coord)[2]
  condition <- FALSE

#########################################################################
  if (corrmodel %in% c("GenWend","Genwend","Hypergeometric","GenWend_Matern","Genwend_Matern","Hypergeometric_Matern")) {
    if (corrmodel %in% c("GenWend","Genwend","GenWend_Matern","Genwend_Matern"))   l <- 0.5
    if (corrmodel %in% c("Hypergeometric","Hypergeometric_Matern"))                l <- d/2 + param$smooth
    
    if (corrmodel %in% c("GenWend_Matern","Genwend_Matern","Hypergeometric_Matern"))  mu <- 1/as.numeric(param['power2'])
    else                                                                              mu <- as.numeric(param['power2'])
    
    C2 = mu/2 - d/2 - 0.5 - param$smooth + l >= 0
    C1 = mu >= 1
    condition <- isTRUE(C1 && C2)
    
    CC <- as.double(param['sill'])
    N <- 1
    a  <- as.double(param['scale'])
    nu1 <- as.double(param['smooth'])
    
    if (corrmodel %in% c("GenWend_Matern","Hypergeometric_Matern"))
      other <- 1/as.numeric(param['power2'])
    else
      other <- as.numeric(param['power2'])
    
    a0 <- a
    nu0 <- nu1
    parameters <- list("CC" = CC, "a" = a,"nu1" = nu1, "nu2" = 0, "other" = other)
    P <- 1
    vtype <- 0
    parametersg <- list("a" = a0, "nu1" = nu0, other = other, CC = CC,
                        mu = mu, kappa = param$smooth, l = l)
    model_num <- CkCorrModel(corrmodel)
  }
#########################################################################
  if (corrmodel == "Matern") {
    CC <- as.double(param['sill'])
    N <- 1
    a  <- as.double(param['scale'])
    nu1 <- as.double(param['smooth'])
    other <- as.double(0)
    a0 <- a
    nu0 <- nu1
    parameters <- list("CC" = CC, "a" = a,"nu1" = nu1, "nu2" = 0, "other" = other)
    P <- 1
    vtype <- 0
    parametersg <- list("a" = a0, "nu1" = nu0, other = other, CC = CC)
    model_num <- CkCorrModel(corrmodel)
  }
#########################################################################
  if (corrmodel %in% c("Kummer","Kummer_Matern")) {
    CC <- as.double(param['sill'])
    N <- 1
    a  <- as.double(param['scale'])
    nu1 <- as.double(param['smooth'])
    other <- as.double(param['power2'])
    a0 <- a
    nu0 <- nu1
    parameters <- list("CC" = CC, "a" = a,"nu1" = nu1,"nu2" = 0, "other" = other)
    P <- 1
    vtype <- 0
    parametersg <- list("a" = a0, "nu1" = nu0, other = other, CC = CC)
    model_num <- CkCorrModel(corrmodel)
  }
######################################################################### 
  u <- frequency_sampler(L, parametersg, corrmodel, d, condition)
  #######
  phi <- runif(L, 0, (2*pi))
  
  Nlocs <- nrow(coord)
  simu <- NULL
  
  if (!is.null(u)) {
    if (Nlocs > 100000) {
      Chunks <- ceiling(Nlocs/100000)
      limits <- unique(c(seq(0, Nlocs, by = 100000), Nlocs))
      low_limits <- limits[1:(length(limits)-1)]
      upp_limits <- limits[2:length(limits)]
      simu <-  matrix(NA, Nlocs, 1)
      
      if (!parallel) {
        for (i in 1:length(low_limits)) {
          indexes <- (low_limits[i]+1):upp_limits[i]
          Ninner <- length(indexes)
          result <- dotCall64::.C64("TBD1d",
                     SIGNATURE = c("double","double","double","double",
                                   "double","integer","integer","double"),
                     ux = c(u[,1]), uy = c(u[,2]), 
                     sx = c(coord[indexes,1]), 
                     sy = c(coord[indexes,2]), 
                     phi = phi, L = as.integer(L), n = as.integer(Ninner), 
                     result = dotCall64::vector_dc("double",Ninner), 
                     INTENT = c("r","r","r","r","r","r","r","rw"),
                     PACKAGE='GeoModels', VERBOSE=0, NAOK=TRUE)$result
          
          simu[indexes,1] <- sqrt(CC)*result/sqrt(2*L)
        }
      } else {
        future::plan(future::multisession, workers = ncores)
        chunk_list <- lapply(1:Chunks, function(k) {
          list(indexes = (low_limits[k]+1):upp_limits[k], chunk_id = k)
        })
        chunk_results <- future.apply::future_lapply(chunk_list, function(chunk_info) {
          indexes <- chunk_info$indexes
          Ninner <- length(indexes)
          result <- dotCall64::.C64("TBD1d",
                     SIGNATURE = c("double","double","double","double",
                                   "double","integer","integer","double"),
                     ux = c(u[,1]), uy = c(u[,2]), 
                     sx = c(coord[indexes,1]), 
                     sy = c(coord[indexes,2]), 
                     phi = phi, L = as.integer(L), n = as.integer(Ninner), 
                     result = dotCall64::vector_dc("double",Ninner), 
                     INTENT = c("r","r","r","r","r","r","r","rw"),
                     PACKAGE='GeoModels', VERBOSE=0, NAOK=TRUE)$result
          result <- sqrt(CC)*result/sqrt(2*L)
          return(list(indexes=indexes, result=result))
        }, future.seed=TRUE, future.packages="GeoModels")
        
        for (chunk_result in chunk_results) {
          simu[chunk_result$indexes,1] <- chunk_result$result
        }
        future::plan(sequential)
      }
    } else {
      result <- dotCall64::.C64("TBD1d",
                   SIGNATURE = c("double","double","double","double",
                                 "double","integer","integer","double"),
                   ux = c(u[,1]), uy = c(u[,2]), 
                   sx = c(coord[,1]), sy = c(coord[,2]), 
                   phi = phi, L = as.integer(L), n = as.integer(Nlocs), 
                   result = dotCall64::vector_dc("double",Nlocs), 
                   INTENT = c("r","r","r","r","r","r","r","rw"),
                   PACKAGE='GeoModels', VERBOSE=0, NAOK=TRUE)$result
      simu <- matrix(sqrt(CC)*result/sqrt(2*L), Nlocs, 1)
    }
  }
  
  return(list(simu=simu, freqs=u, parametersg=parametersg))
}
#########################################################
#########################################################
#########################################################
frequency_sampler <- function(L, parametersg, corrmodel, d, condition) {
  
  ini.sample <- NULL
  #########################################################
  if (corrmodel == "Matern") {
    spl <- rgamma(L, parametersg$nu1, scale = 1)
    G <- matrix(spl, nrow = L, ncol = d)
    N <- matrix(rnorm(L * d), nrow = L, ncol = d)
    ini.sample <- (N / sqrt(G*2)) / parametersg$a
    ini.sample <- ini.sample / (2*pi)
  }
  ######################################################### 
  if (corrmodel %in% c("Kummer", "Kummer_Matern")) {
    nu    <- as.numeric(parametersg$nu1)
    alpha <- as.numeric(parametersg$other)

    if (corrmodel == "Kummer") {
      betaK <- as.numeric(parametersg$a)
    } else {
      betaK <- as.numeric(parametersg$a) * sqrt(2 * (alpha + 1))
    }
    spl <- rf(L, df1 = 2*nu, df2 = 2*alpha)
    G <- matrix( (nu / alpha) * spl, nrow = L, ncol = d)
    G <- pmax(G, .Machine$double.eps)
    #G <- matrix( (nu / alpha) * rf(L * d, df1 = 2*nu, df2 = 2*alpha), nrow = L, ncol = d)
    N <- matrix(rnorm(L * d), nrow = L, ncol = d)
    ini.sample <- (N / sqrt(G)) / betaK
    ini.sample <- ini.sample / (2*pi)
  }
#########################################################
  if (corrmodel %in% c("GenWend","GenWend_Matern","Genwend","Genwend_Matern",
                       "Hypergeometric","Hypergeometric_Matern")) {
    
    kappa <- as.numeric(parametersg$kappa)
    mu    <- as.numeric(parametersg$mu)
    l     <- as.numeric(parametersg$l)
    
    if (corrmodel %in% c("GenWend","Genwend","Hypergeometric")) {
      a <- as.numeric(parametersg$a)
    } else if (corrmodel %in% c("GenWend_Matern","Genwend_Matern")) {
      a <- as.numeric(parametersg$a) * 
           (gamma(mu + 2*kappa + 1) / gamma(mu))^(1 / (1 + 2*kappa))
    } else if (corrmodel == "Hypergeometric_Matern") {
      num <- (2*kappa + 1) * gamma((mu + 1)/2 + kappa) * 
             gamma((mu + d + 1)/2 + 2*kappa)
      den <- gamma(mu/2) * gamma((mu + d)/2 + kappa)
      a <- as.numeric(parametersg$a) * (num / den)^(1 / (1 + 2*kappa))
    }
    
    if (condition) {
      #message("==> Using Beta-mixture method")
      ini.sample <- simulate_GH_beta(L, d, kappa, mu, l, a)
      ini.sample <- ini.sample / (2 * pi)
    } else {
      #message("==> Using Gasper mixture method")
      ini.sample <- simulate_GH_gasper(L, d, kappa, mu, l, a, nmax = 15)
      ini.sample <- ini.sample / (2 * pi)
    }
  }
  
  return(ini.sample)
}