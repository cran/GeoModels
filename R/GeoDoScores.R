####################################################
### File name: GeoDoScores.r
####################################################

GeoDoScores <- function(data, method = "cholesky", matrix) {

  if (!inherits(matrix, "GeoCovmatrix")) {
    stop("A GeoCovmatrix object is needed as input\n")
  }

  method <- tolower(gsub("[[:blank:]]", "", as.character(method)[1L]))
  if (!(method %in% c("cholesky", "svd"))) {
    stop("method must be 'cholesky' or 'svd'", call. = FALSE)
  }

  ## The data length is the authoritative number of observations.  This also
  ## covers dynamic space-time layouts, where length(coordx) is not nsite.
  dime <- length(data)
  varcov <- matrix$covmatrix
  if (nrow(varcov) != length(data)) {
    stop("The dimension of the covariance matrix and/or the vector data are not correct\n")
  }
  dimnames(varcov) <- NULL
  data <- as.numeric(data)

  ## Return Q b and diag(Q), where Q is the inverse covariance matrix.  The
  ## Cholesky and SVD paths deliberately use their corresponding decomposition;
  ## in particular an svd object must never be passed to forwardsolve().
  getInv2 <- function(covmatrix, b) {
    if (!covmatrix$sparse) {
      dec <- MatDecomp(covmatrix$covmatrix, method)
      if (is.logical(dec)) {
        stop("Covariance matrix is not positive definite")
      }

      if (identical(method, "cholesky")) {
        vec <- forwardsolve(dec, b)
        Invc <- forwardsolve(dec, vec, transpose = TRUE)
        ## Q = U^{-T} U^{-1}.  GeoDoScores needs only diag(Q).
        I <- diag(nrow(covmatrix$covmatrix))
        Uinv <- forwardsolve(dec, I)
        diagInv <- colSums(Uinv^2)
      } else {
        ## For A = U D V^T, A^{-1} = V D^{-1} U^T.
        invd <- 1 / dec$d
        Invc <- as.numeric(dec$v %*% (invd * crossprod(dec$u, b)))
        diagInv <- rowSums(sweep(dec$v * dec$u, 2L, dec$d, "/"))
      }
    } else {
      if (!identical(method, "cholesky")) {
        stop("method='svd' is not available for sparse GeoCovmatrix objects", call. = FALSE)
      }

      cc <- covmatrix$covmatrix
      if (!spam::is.spam(cc)) cc <- spam::as.spam(cc)

      U <- tryCatch(
        spam::chol.spam(cc),
        error = function(e) stop("Covariance matrix is not positive definite")
      )
      vec <- spam::forwardsolve(U, b)
      Invc <- spam::backsolve(U, vec)
      Inv <- spam::solve.spam(U)
      diagInv <- diag(Inv)
    }

    diagInv <- as.numeric(diagInv)
    if (length(diagInv) != length(b) || any(!is.finite(diagInv)) || any(diagInv <= 0)) {
      stop("Invalid diagonal of the inverse covariance matrix", call. = FALSE)
    }
    list(a = as.numeric(Invc), c = diagInv)
  }

  param <- matrix$param
  external_mean <- matrix$mean_external
  if (is.null(external_mean)) external_mean <- .GeoMean_external(param, n = dime)
  if (!is.null(external_mean)) {
    external_mean <- as.numeric(external_mean)
    if (length(external_mean) != dime || any(!is.finite(external_mean)))
      stop("The external mean stored in the covariance object is invalid.")
    if (!is.null(matrix$X))
      stop("A fixed observation mean and X cannot be used together.")
    MM <- external_mean
  } else {
    matrix$X <- .GeoMean_design(matrix$X, n = dime, name = "X")
    beta <- .GeoMean_beta(param, ncol(matrix$X))
    MM <- as.numeric(matrix$X %*% beta)
  }

  namesnuis <- matrix$namesnuis
  nuisance <- param[namesnuis]
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

  ## LOO identities for Q = Sigma^{-1}:
  ##   e_i       = (Q r)_i / Q_ii
  ##   Var_i_LOO = 1 / Q_ii
  ##   z_i       = e_i / sqrt(Var_i_LOO) = (Q r)_i / sqrt(Q_ii)
  cc <- getInv2(matrix, data)
  temp <- cc$a
  vv <- cc$c
  sd_loo <- 1 / sqrt(vv)
  residual_loo <- temp / vv
  zstd <- temp / sqrt(vv)

  MAD <- median(abs(residual_loo))
  RMSE <- sqrt(mean(residual_loo^2))
  MAE <- mean(abs(residual_loo))
  LSCORE <- mean(0.5 * (log(2 * pi / vv) + zstd^2))
  CRPS <- mean(
    sd_loo * (zstd * (2 * pnorm(zstd) - 1) +
                2 * dnorm(zstd) - 1 / sqrt(pi))
  )

  list(
    RMSE   = RMSE,
    LSCORE = LSCORE,
    MAD    = MAD,
    CRPS   = CRPS,
    MAE    = MAE
  )
}
