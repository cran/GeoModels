############################################################
### Circulant-embedding helpers for spatial and separable
### spatio-temporal Gaussian random fields.
############################################################

.CEValidateScalar <- function(x, name, lower = -Inf, upper = Inf,
                              lower_open = FALSE, upper_open = FALSE) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x)) {
    stop(name, " must be a finite numeric scalar.", call. = FALSE)
  }
  bad_lower <- if (lower_open) x <= lower else x < lower
  bad_upper <- if (upper_open) x >= upper else x > upper
  if (bad_lower || bad_upper) {
    left <- if (lower_open) "(" else "["
    right <- if (upper_open) ")" else "]"
    stop(name, " must belong to ", left, lower, ", ", upper, right, ".",
         call. = FALSE)
  }
  as.numeric(x)
}

.CERegularAxis <- function(x, n, name) {
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop(name, " must be a finite numeric vector.", call. = FALSE)
  }
  if (length(n) != 1L || !is.numeric(n) || !is.finite(n) ||
      n != as.integer(n) || n < 2L) {
    stop("The dimension associated with ", name,
         " must be an integer greater than or equal to two.", call. = FALSE)
  }
  n <- as.integer(n)
  if (length(x) != n) {
    stop(name, " must contain exactly ", n, " coordinates.", call. = FALSE)
  }
  dx <- diff(x)
  if (any(dx <= 0)) {
    stop(name, " must be strictly increasing.", call. = FALSE)
  }
  tol <- sqrt(.Machine$double.eps) * max(1, abs(dx[1L]), max(abs(x)))
  if (max(abs(dx - dx[1L])) > tol) {
    stop(name, " must be equally spaced for circulant embedding.",
         call. = FALSE)
  }
  unname(as.numeric(x))
}

.CEValidateMaxExt <- function(max.ext) {
  if (length(max.ext) != 1L || !is.numeric(max.ext) || !is.finite(max.ext) ||
      max.ext != as.integer(max.ext) || max.ext < 1L) {
    stop("max.ext must be a positive integer.", call. = FALSE)
  }
  as.integer(max.ext)
}

.CEVarianceParameters <- function(param) {
  sill <- if (is.null(param$sill)) 1 else param$sill
  nugget <- if (is.null(param$nugget)) 0 else param$nugget
  sill <- .CEValidateScalar(sill, "param$sill", lower = 0,
                            lower_open = TRUE)
  nugget <- .CEValidateScalar(nugget, "param$nugget", lower = 0,
                              upper = 1)
  list(sill = sill, nugget = nugget)
}

.CEStandardCorrelationParam <- function(param) {
  out <- param
  out$mean <- 0
  out$sill <- 1
  out$nugget <- 0
  out
}

.CEEigenInfo <- function(lambda) {
  if (!is.numeric(lambda) || any(!is.finite(lambda))) {
    stop("The circulant embedding produced non-finite eigenvalues.",
         call. = FALSE)
  }
  scale <- max(1, max(abs(lambda)))
  tol <- sqrt(.Machine$double.eps) * scale
  list(valid = min(lambda) >= -tol, tol = tol,
       values = pmax(lambda, 0), minimum = min(lambda))
}

.CESeparableModel <- function(corrmodel, param) {
  make_param <- function(prefix = "") {
    suffix <- if (nzchar(prefix)) paste0("_", prefix) else ""
    out <- list(nugget = 0, sill = 1)
    for (nm in c("scale", "smooth", "power2")) {
      value <- param[[paste0(nm, suffix)]]
      if (!is.null(value)) out[[nm]] <- value
    }
    out
  }

  if (identical(corrmodel, "Matern_Matern")) {
    return(list(model = "Matern", spatial = make_param("s"),
                temporal = make_param("t")))
  }
  if (identical(corrmodel, "GenWend_GenWend")) {
    return(list(model = "GenWend", spatial = make_param("s"),
                temporal = make_param("t")))
  }
  if (identical(corrmodel, "GenWend_Matern_GenWend_Matern")) {
    return(list(model = "GenWend_Matern", spatial = make_param("s"),
                temporal = make_param("t")))
  }
  if (identical(corrmodel, "Kummer_Kummer")) {
    return(list(model = "Kummer", spatial = make_param("s"),
                temporal = make_param("t")))
  }
  if (identical(corrmodel, "Kummer_Matern_Kummer_Matern")) {
    return(list(model = "Kummer_Matern", spatial = make_param("s"),
                temporal = make_param("t")))
  }

  stop("CE spatio-temporal simulation is implemented only for supported separable correlation models.",
       call. = FALSE)
}

.CETemporalEmbedding <- function(time.seq, corrmodel, param, max.ext) {
  nt <- length(time.seq)
  if (nt < 2L) {
    stop("At least two times are required for a temporal embedding.",
         call. = FALSE)
  }
  dt <- diff(time.seq)
  tol_time <- sqrt(.Machine$double.eps) *
    max(1, abs(dt[1L]), max(abs(time.seq)))
  if (any(dt <= 0)) {
    stop("coordt must be strictly increasing for CE simulation.",
         call. = FALSE)
  }
  if (max(abs(dt - dt[1L])) > tol_time) {
    stop("coordt must be equally spaced for CE spatio-temporal simulation.",
         call. = FALSE)
  }

  base_size <- 2L * (nt - 1L)
  last <- NULL
  for (attempt in seq_len(max.ext)) {
    emb_size <- as.integer(base_size * 2^(attempt - 1L))
    half <- emb_size %/% 2L
    lags <- (0:half) * dt[1L]
    rho <- as.numeric(GeoCorrFct(
      x = lags, corrmodel = corrmodel, model = "Gaussian",
      param = param
    )$corr)
    second_half <- if (half > 1L) rev(rho[2:half]) else numeric(0)
    first_row <- c(rho, second_half)
    lambda <- Re(fft(first_row, inverse = TRUE))
    info <- .CEEigenInfo(lambda)
    last <- list(size = emb_size, lambda = info$values,
                 valid = info$valid, minimum = info$minimum,
                 attempt = attempt)
    if (info$valid) return(last)
  }

  warning("The temporal circulant embedding has negative eigenvalues at the maximum extension; negative eigenvalues were truncated to zero.",
          call. = FALSE)
  last$truncated <- TRUE
  last
}

######################## space-time case ########################
## Separable simulation using Cholesky in space and a circulant
## embedding/FFT in time. Spatial locations must be fixed.
CE_Space_Time <- function(coords, time.seq, param, corrmodel, distance,
                          max.ext = 1L) {
  coords <- as.matrix(coords)
  if (!is.numeric(coords) || nrow(coords) < 1L ||
      !(ncol(coords) %in% c(2L, 3L)) || any(!is.finite(coords))) {
    stop("coords must be a finite numeric matrix with two or three columns.",
         call. = FALSE)
  }
  if (!is.numeric(time.seq) || !length(time.seq) ||
      any(!is.finite(time.seq))) {
    stop("coordt must be a non-empty finite numeric vector.",
         call. = FALSE)
  }
  max.ext <- .CEValidateMaxExt(max.ext)
  variance <- .CEVarianceParameters(param)
  sep <- .CESeparableModel(corrmodel, param)
  sep$spatial$mean <- 0
  sep$temporal$mean <- 0

  ns <- nrow(coords)
  nt <- length(time.seq)
  normalized_time <- as.numeric(time.seq - time.seq[1L])
  prepared <- .GeoSimulationMemoise(
    "ce_spacetime",
    list(coords = coords, time = normalized_time, sep = sep,
         distance = distance, max.ext = max.ext),
    function() {
      spatial_cov <- GeoCovmatrix(
        coordx = coords, corrmodel = sep$model, model = "Gaussian",
        param = sep$spatial, distance = distance
      )$covmatrix
      spatial_cov <- as.matrix(spatial_cov)
      diag(spatial_cov) <- 1
      chol_space <- tryCatch(
        t(chol(spatial_cov)),
        error = function(e) {
          stop("The spatial correlation matrix is not positive definite in CE_Space_Time: ",
               conditionMessage(e), call. = FALSE)
        }
      )
      embedding <- if (nt > 1L) {
        .CETemporalEmbedding(normalized_time, sep$model, sep$temporal, max.ext)
      } else NULL
      list(chol_space = chol_space, embedding = embedding)
    }
  )
  chol_space <- prepared$chol_space
  if (nt == 1L) {
    base <- matrix(as.numeric(chol_space %*% rnorm(ns)), nrow = 1L)
  } else {
    embedding <- prepared$embedding
    emb_size <- embedding$size
    spectrum <- sqrt(embedding$lambda)
    freq <- matrix(complex(real = 0, imaginary = 0),
                   nrow = emb_size, ncol = ns)
    for (j in seq_len(emb_size)) {
      z <- complex(real = rnorm(ns), imaginary = rnorm(ns))
      freq[j, ] <- spectrum[j] * as.vector(chol_space %*% z)
    }
    temporal <- matrix(complex(real = 0, imaginary = 0),
                       nrow = emb_size, ncol = ns)
    for (j in seq_len(ns)) {
      temporal[, j] <- fft(freq[, j], inverse = TRUE) / sqrt(emb_size)
    }
    base <- Re(temporal[seq_len(nt), , drop = FALSE])
  }

  if (variance$nugget > 0) {
    noise <- matrix(rnorm(nt * ns), nrow = nt, ncol = ns)
    field <- sqrt(variance$sill) *
      (sqrt(1 - variance$nugget) * base +
         sqrt(variance$nugget) * noise)
  } else {
    field <- sqrt(variance$sill) * base
  }

  ## Temporal blocks: all spatial sites at time 1, then time 2, etc.
  c(t(field))
}

############################################################
SimCE <- function(M, N, x, y, z = NULL, corrmodel, param,
                  mean.val = 0, max.ext = 1L) {
  if (!is.null(z)) {
    stop("Spatial CE simulation is implemented only for two-dimensional grids.",
         call. = FALSE)
  }
  x <- .CERegularAxis(x, M, "x")
  y <- .CERegularAxis(y, N, "y")
  max.ext <- .CEValidateMaxExt(max.ext)
  variance <- .CEVarianceParameters(param)

  out <- GRF.CE.Sim(
    nsim = 1L, x = x, y = y, M = length(x), N = length(y),
    corrmodel = corrmodel, param = .CEStandardCorrelationParam(param),
    max.ext = max.ext
  )

  mean.val <- as.numeric(mean.val)
  ngrid <- length(x) * length(y)
  if (any(!is.finite(mean.val)) || !(length(mean.val) %in% c(1L, ngrid))) {
    stop("mean.val must be a finite scalar or contain one value per grid point.",
         call. = FALSE)
  }
  if (length(mean.val) == 1L) mean.val <- rep(mean.val, ngrid)

  latent <- as.matrix(out$X)
  if (variance$nugget > 0) {
    nugget_noise <- matrix(rnorm(length(latent)), nrow = nrow(latent),
                           ncol = ncol(latent))
    field <- sqrt(variance$sill) *
      (sqrt(1 - variance$nugget) * latent +
         sqrt(variance$nugget) * nugget_noise)
  } else {
    field <- sqrt(variance$sill) * latent
  }
  out$X <- sweep(field, 1L, mean.val, "+")
  out
}

##################
GRF.CE.Sim <- function(nsim, x, y, M, N, corrmodel, param, max.ext) {
  if (length(nsim) != 1L || !is.numeric(nsim) || !is.finite(nsim) ||
      nsim != as.integer(nsim) || nsim < 1L) {
    stop("nsim must be a positive integer.", call. = FALSE)
  }
  nsim <- as.integer(nsim)
  max.ext <- .CEValidateMaxExt(max.ext)
  prepared <- .GeoSimulationMemoise(
    "ce_spatial",
    list(x = x, y = y, M = M, N = N, corrmodel = corrmodel,
         param = param, max.ext = max.ext),
    function() {
      decomposition <- NULL
      selected_attempt <- max.ext
      for (attempt in seq_len(max.ext)) {
        ext <- 2^attempt
        candidate <- Cov_Mat_CE(
          x = x, y = y, M = M, N = N, corrmodel = corrmodel,
          param = param, ext = ext, approx = FALSE
        )
        decomposition <- candidate
        selected_attempt <- attempt
        if (isTRUE(candidate$valid)) break
      }
      approx.flag <- !isTRUE(decomposition$valid)
      if (approx.flag) {
        warning("The spatial circulant embedding has negative eigenvalues at the maximum extension; negative eigenvalues were truncated to zero.",
                call. = FALSE)
        decomposition$prefix <- sqrt(pmax(decomposition$eigenvalues, 0))
      }
      list(decomposition = decomposition, selected_attempt = selected_attempt,
           approx.flag = approx.flag)
    }
  )
  decomposition <- prepared$decomposition
  selected_attempt <- prepared$selected_attempt
  approx.flag <- prepared$approx.flag

  simulation <- Cov_Mat_CE.sim(decomposition, nsim = nsim)
  list(
    X = simulation,
    grid.points = as.matrix(decomposition$grid),
    k0 = selected_attempt,
    approx.flag = approx.flag,
    min.eigenvalue = decomposition$min.eigenvalue
  )
}

##########################################
Cov_Mat_CE <- function(x, y, M, N, corrmodel, param, ext = 2,
                       approx = FALSE) {
  mygrid <- grid.prep(x, y, M = M, N = N, ext = ext)
  rx <- mygrid$M.ext * mygrid$cell.width
  ry <- mygrid$N.ext * mygrid$cell.height

  dx_abs <- abs(mygrid$x.ext[1L] - mygrid$x.ext)
  dx <- pmin(dx_abs, rx - dx_abs)
  dy_abs <- abs(mygrid$y.ext[1L] - mygrid$y.ext)
  dy <- pmin(dy_abs, ry - dy_abs)
  torus <- expand.grid(dx = dx, dy = dy)
  distances <- matrix(
    sqrt(torus[, 1L]^2 + torus[, 2L]^2),
    nrow = mygrid$M.ext, ncol = mygrid$N.ext
  )

  sigma_row <- my.GeoCorrFct(
    x = c(distances), corrmodel = corrmodel,
    param = .CEStandardCorrelationParam(param)
  )
  sigma_row <- matrix(sigma_row, nrow = mygrid$M.ext,
                      ncol = mygrid$N.ext)
  lambda <- Re(fft(sigma_row, inverse = TRUE))
  info <- .CEEigenInfo(lambda)
  if (isTRUE(approx)) info$values <- pmax(lambda, 0)

  list(
    dp = length(lambda), d = dim(lambda), sdp = sqrt(length(lambda)),
    M = M, N = N, prefix = sqrt(info$values),
    eigenvalues = lambda, valid = info$valid || isTRUE(approx),
    min.eigenvalue = info$minimum,
    grid = expand.grid(x = mygrid$x, y = mygrid$y)
  )
}

##################
grid.prep <- function(x, y, M, N, ext = 2) {
  x <- .CERegularAxis(x, M, "x")
  y <- .CERegularAxis(y, N, "y")
  if (length(ext) != 1L || !is.numeric(ext) || !is.finite(ext) ||
      ext != as.integer(ext) || ext < 2L) {
    stop("The embedding extension must be an integer greater than or equal to two.",
         call. = FALSE)
  }
  ext <- as.integer(ext)
  M.ext <- ext * M
  N.ext <- ext * N
  cell.width <- diff(x)[1L]
  cell.height <- diff(y)[1L]
  x.ext <- x[1L] + (0:(M.ext - 1L)) * cell.width
  y.ext <- y[1L] + (0:(N.ext - 1L)) * cell.height

  list(
    M = M, N = N, x = x, y = y,
    cell.width = cell.width, cell.height = cell.height,
    M.ext = M.ext, N.ext = N.ext,
    x.ext = x.ext, y.ext = y.ext
  )
}

##################
my.GeoCorrFct <- function(x, corrmodel, param) {
  if (identical(corrmodel, "none")) {
    result <- rep(0, length(x))
    result[x == 0] <- 1
    return(result)
  }
  as.numeric(GeoCorrFct(
    x = x, corrmodel = corrmodel, model = "Gaussian",
    param = .CEStandardCorrelationParam(param)
  )$corr)
}

##################
Cov_Mat_CE.sim <- function(ex1, nsim = 2L) {
  dp <- ex1$dp
  d <- ex1$d
  sdp <- ex1$sdp
  M <- ex1$M
  N <- ex1$N
  prefix <- ex1$prefix
  if (length(prefix) != dp || any(!is.finite(prefix))) {
    stop("Invalid spectral decomposition in Cov_Mat_CE.sim.",
         call. = FALSE)
  }

  X <- matrix(0, nrow = M * N, ncol = nsim)
  for (i in seq_len(nsim)) {
    standard <- rnorm(dp)
    spectral <- prefix * (fft(matrix(standard, d[1L], d[2L])) / sdp)
    draw <- Re(fft(spectral, inverse = TRUE) / sdp)
    X[, i] <- as.vector(draw[seq_len(M), seq_len(N), drop = FALSE])
  }
  X
}
