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
  model <- matrix$model
  if (isTRUE(matrix$bivariate)) {
    if (!identical(as.integer(model), 1L))
      stop("GeoDoScores supports bivariate covariance objects only for Gaussian margins.",
           call. = FALSE)
    MM <- .GeoBivariateObservedMean(param, matrix$X, matrix$ns,
                                    context = "GeoDoScores")
    marginal_mean <- MM
  } else {
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
    marginal_mean <- .GeoMarginalMean(model, param, MM, n = matrix$n,
                                      copula = matrix$copula)
  }
  if (length(marginal_mean) != dime || any(!is.finite(marginal_mean)))
    stop("The model does not have a finite marginal mean at every observation.",
         call. = FALSE)
  data <- data - marginal_mean

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
