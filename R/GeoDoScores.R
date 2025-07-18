####################################################
### File name: GeoDoScores.r
####################################################

GeoDoScores <- function(data, method = "cholesky", matrix) {
  
  if (!inherits(matrix, "GeoCovmatrix")) {
    stop("A GeoCovmatrix object is needed as input\n")
  }
  nsites <- length(matrix$coordx)
  ntime <- if (matrix$spacetime) length(matrix$coordt) else if (matrix$bivariate) 2L else 1L
  dime <- nsites * ntime
  varcov <- matrix$covmatrix
  if (nrow(varcov) != length(data)) {
    stop("The dimension of the covariance matrix and/or the vector data are not correct\n")
  }
  # Rimozione nomi per efficienza
  dimnames(varcov) <- NULL
  data <- as.numeric(data)
  
########## Funzione interna ottimizzata ############################
  getInv2 <- function(covmatrix, b) {
    if (!covmatrix$sparse) {
      # Case non sparse
      U <- MatDecomp(covmatrix$covmatrix, method)
      if (is.logical(U)) {
        stop("Covariance matrix is not positive definite")
      }
      
      vec <- forwardsolve(U, b)
      Invc <- forwardsolve(U, vec, transpose = TRUE)
      Inv <- FastGP::rcppeigen_invert_matrix(covmatrix$covmatrix)
      diagInv <- diag(Inv)
      
    } else {
      # Case sparse
      cc <- covmatrix$covmatrix
      
      # Conversione efficiente a spam se necessario
      if (!spam::is.spam(cc)) {
        cc <- spam::as.spam(cc)
      }
      
      U <- tryCatch(
        spam::chol.spam(cc),
        error = function(e) stop("Covariance matrix is not positive definite")
      )
      
      vec <- spam::forwardsolve(U, b)
      Invc <- spam::backsolve(U, vec)
      Inv <- spam::solve.spam(U)
      diagInv <- diag(Inv)
    }
    
    list(a = Invc, b = Inv, c = diagInv)
  }
################################################
  
  param <- matrix$param
  if (is.null(matrix$X)) {
    matrix$X <- matrix(1, ncol = 1, nrow = dime)
  }
  MM <- 0
  namesnuis <- matrix$namesnuis
  nuisance <- param[namesnuis]
  
  if (length(param$mean) == 1L) {
    sel <- startsWith(names(nuisance), "mean")
    mm <- as.numeric(nuisance[sel])
    MM <- matrix$X %*% mm
  } else {
    MM <- param$mean
  }
  model <- matrix$model
  data <- switch(
    as.character(model),
    # Gaussian, StudentT, Tukey, Logistic
    "1" = , "35" = , "12" = , "34" = , "25" = data - MM,
    # SkewGaussian
    "10" = {
      kk <- as.numeric(param['skew'])
      data - (MM + kk * sqrt(2/pi))
    },
    # Binomial
    "11" = data - matrix$n * pnorm(MM),
    # Poisson
    "30" = , "36" = data - exp(MM),
    # Poisson inflated
    "43" = , "44" = {
      p <- pnorm(as.numeric(param['pmu']))
      data - (1 - p) * exp(MM)
    },
    # Two piece t models
    "27" = {
      kk <- as.numeric(param['skew'])
      dd <- as.numeric(param['df'])
      ss <- as.numeric(param['sill'])
      adjustment <- (2 * kk * sqrt(ss * dd) * gamma((dd - 1) / 2)) / 
                   (gamma(dd / 2) * sqrt(pi))
      data - (MM - adjustment)
    },
    # Two piece gaussian
    "29" = {
      kk <- as.numeric(param['skew'])
      ss <- as.numeric(param['sill'])
      data - (MM - 2 * kk * sqrt(2 * ss / pi))
    },
    # Two piece tukeyh
    "38" = {
      kk <- as.numeric(param['skew'])
      ss <- as.numeric(param['sill'])
      tt <- as.numeric(param['tail'])
      data - (MM - 2 * kk * sqrt(2 * ss / pi) / (1 - tt))
    },
    # Tukeyh2
    "40" = {
      ss <- as.numeric(param['sill'])
      t1 <- as.numeric(param['tail1'])
      t2 <- as.numeric(param['tail2'])
      data - (MM + sqrt(ss) * (t1 - t2) / (sqrt(2 * pi) * (1 - t1) * (1 - t2)))
    },
    # Binomial negative
    "16" = data - matrix$n * (1 - pnorm(MM)) / pnorm(MM),
    
    # Binomial negative inflated
    "45" = {
      p <- pnorm(as.numeric(param['pmu']))
      data - (1 - p) * matrix$n * (1 - pnorm(MM)) / pnorm(MM)
    },
    
    # SAS
    "20" = {
      ss <- as.numeric(param['sill'])
      kk <- as.numeric(param['skew'])
      tt <- as.numeric(param['tail'])
      bessel_term <- besselK(0.25, (tt + 1) / (2 * tt)) + 
                     besselK(0.25, (1 - tt) / (2 * tt))
      adjustment <- sqrt(ss) * sinh(kk / tt) * exp(0.25) * bessel_term / sqrt(8 * pi)
      data - (MM + adjustment)
    },
    data
  )
  # computing inverse
  cc <- getInv2(matrix, data)
  temp <- cc$a  # inv %*% data
  inv <- cc$b   # inv
  vv <- cc$c    # diag inv
  inv_sqrt_vv <- 1 / sqrt(vv)
  inv_vv <- 1 / vv
  z <- inv_vv * temp
  zz <- inv_sqrt_vv * temp
  dime_inv <- 1 / dime

  MAD <- median(z)
  RMSE <- sqrt(dime_inv * sum(z^2))
  MAE <- dime_inv * sum(abs(z))
  LSCORE <- 0.5 * dime_inv * (sum(log(2 * pi * vv)) + sum(zz^2))
  pnorm_zz <- pnorm(zz)
  CRPS <- dime_inv * (
    sum(inv_sqrt_vv * zz * (2 * pnorm_zz - 1)) +
    2 * sum(inv_sqrt_vv * pnorm_zz) +
    sum(inv_sqrt_vv) / sqrt(pi)
  )
  
  # Results
  list(
    RMSE   = RMSE,
    LSCORE = LSCORE,
    MAD    = MAD,
    CRPS   = CRPS,
    MAE    = MAE
  )
}