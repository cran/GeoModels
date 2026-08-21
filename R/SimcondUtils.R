# Internal engines used by GeoSimcond().
# Keeping them outside the public function avoids rebuilding large closures on
# every call and keeps the user-facing implementation focused on dispatch.

.GeoSimcondNormalizeSimulations <- function(x, expected, expected_length = NULL) {
  expected <- as.integer(expected)[1L]
  if (is.list(x)) {
    if (length(x) != expected) stop("Unexpected number of simulations.", call. = FALSE)
    out <- lapply(x, as.numeric)
  } else if (is.matrix(x)) {
    if (!is.null(expected_length) &&
        ncol(x) == expected && nrow(x) == expected_length) {
      out <- lapply(seq_len(expected), function(j) as.numeric(x[, j]))
    } else if (!is.null(expected_length) &&
               nrow(x) == expected && ncol(x) == expected_length) {
      out <- lapply(seq_len(expected), function(j) as.numeric(x[j, ]))
    } else if (is.null(expected_length) && ncol(x) == expected) {
      out <- lapply(seq_len(expected), function(j) as.numeric(x[, j]))
    } else if (is.null(expected_length) && nrow(x) == expected) {
      out <- lapply(seq_len(expected), function(j) as.numeric(x[j, ]))
    } else if (expected == 1L &&
               (is.null(expected_length) || length(x) == expected_length)) {
      out <- list(as.numeric(x))
    } else {
      stop("Unexpected matrix format for simulations.", call. = FALSE)
    }
  } else if (expected == 1L) {
    out <- list(as.numeric(x))
  } else {
    stop("Unexpected format of simulations.", call. = FALSE)
  }

  if (!is.null(expected_length) &&
      any(vapply(out, length, integer(1L)) != expected_length)) {
    stop("A simulation has an unexpected length.", call. = FALSE)
  }
  out
}

.GeoSimcondSkewGaussian <- function(
    coord_obs, loc, data, param, corrmodel,
    nrep = 1, method = "Cholesky", local = FALSE,
    neighb = NULL, L = NULL,
    distance = "Eucl", radius = 1, anisopars = NULL,
    parallel = FALSE, ncores = 1, progress = FALSE,
    n_iter = 1000L, batch_size = NULL, max_batch_mb = 256,
    mean_obs = NULL, mean_loc = NULL
) {

  ##########################################################################
  ## Parameters and checks
  ##########################################################################

  sill       <- as.numeric(param$sill)
  skew       <- as.numeric(param$skew)
  sigma_sqrt <- sqrt(sill)

  if (length(sill) != 1L || !is.finite(sill) || sill <= 0) {
    stop("param$sill must be strictly positive.")
  }

  if (length(skew) != 1L || !is.finite(skew)) {
    stop("param$skew must be finite.")
  }

  if (length(nrep) != 1L || !is.finite(nrep) || nrep < 1) {
    stop("nrep must be at least 1.")
  }

  if (length(n_iter) != 1L || !is.finite(n_iter) || n_iter < 1) {
    stop("n_iter must be at least 1.")
  }

  if (local) {
    stop(
      "Local prediction cannot be used with the ",
      "skew-Gaussian Gibbs sampler."
    )
  }

  if (!(method %in% c("Cholesky", "TB", "CE"))) {
    stop("method must be one of 'Cholesky', 'TB', or 'CE'.")
  }

  nrep  <- as.integer(nrep)
  n_iter <- as.integer(n_iter)

  coord_obs <- as.matrix(coord_obs)
  loc       <- as.matrix(loc)
  data      <- as.numeric(data)

  M0 <- nrow(coord_obs)
  N0 <- nrow(loc)

  if(is.null(mean_obs)) mean_obs <- rep(as.numeric(param$mean)[1L], M0)
  if(is.null(mean_loc)) mean_loc <- rep(as.numeric(param$mean)[1L], N0)
  mean_obs <- as.numeric(mean_obs)
  mean_loc <- as.numeric(mean_loc)
  if(length(mean_obs) != M0 || any(!is.finite(mean_obs)))
    stop("mean_obs must contain one finite value per observation.")
  if(length(mean_loc) != N0 || any(!is.finite(mean_loc)))
    stop("mean_loc must contain one finite value per prediction point.")
  ## Existing worker code uses mm for the prediction-side location vector.
  mm <- mean_loc

  if (length(data) != M0) {
    stop(
      "The length of data must match the number ",
      "of observed locations."
    )
  }

  if (N0 < 1L) {
    stop("loc must contain at least one prediction location.")
  }

  ##########################################################################
  ## Standardized latent representation
  ##
  ## Y(s) = mm + skew |G1(s)| + sqrt(sill) G2(s)
  ##
  ## Z(s) = gamma |G1(s)| + G2(s),
  ## gamma = skew / sqrt(sill)
  ##########################################################################

  gamma <- skew / sigma_sqrt
  z_obs <- (data - mean_obs) / sigma_sqrt

  if (any(!is.finite(z_obs))) {
    stop("The standardized observations contain non-finite values.")
  }

  gaussian_case <- isTRUE(gamma == 0)

  keep       <- CorrParam(corrmodel)
  corr_param <- param[keep]

  latent_nugget <- if (is.null(param$nugget)) 0 else as.numeric(param$nugget)
  if (length(latent_nugget) != 1L || !is.finite(latent_nugget) ||
      latent_nugget < 0 || latent_nugget >= 1) {
    stop("param$nugget must be a finite scalar in [0, 1).")
  }

  latent_param <- c(
    corr_param,
    mean = 0,
    nugget = latent_nugget,
    sill = 1
  )

  ##########################################################################
  ## Gaussian covariance matrix and kriging weights
  ##########################################################################

  cat("Computing Kriging weights ...\n")

  GeoW <- GeoKrigWeights(
    coordx = coord_obs,
    corrmodel = corrmodel,
    loc = loc,
    model = "Gaussian",
    distance = distance,
    radius = radius,
    anisopars = anisopars,
    param = latent_param
  )

  Wglob <- GeoW$weights

  if (is.null(Wglob)) {
    stop("GeoKrigWeights did not return the kriging weights.")
  }

  if (is.null(dim(Wglob))) {
    Wglob <- matrix(Wglob, nrow = M0)
  }

  weights_orientation <- GeoW$weights_orientation
  if (is.null(weights_orientation) ||
      !(weights_orientation %in% c("nobs_by_nloc", "nloc_by_nobs"))) {
    stop("GeoKrigWeights returned an unknown weights orientation.")
  }
  if (weights_orientation == "nobs_by_nloc" &&
      (nrow(Wglob) != M0 || ncol(Wglob) != N0)) {
    stop("The nobs_by_nloc kriging weights have incompatible dimensions.")
  }
  if (weights_orientation == "nloc_by_nobs" &&
      (nrow(Wglob) != N0 || ncol(Wglob) != M0)) {
    stop("The nloc_by_nobs kriging weights have incompatible dimensions.")
  }

  # The precision matrix is unnecessary when gamma = 0.
  if (gaussian_case) {

    Sigma.mat.inv <- numeric(0L)

  } else {

    if (is.null(GeoW$covmatrix)) {
      stop(
        "GeoKrigWeights did not return the ",
        "latent covariance matrix."
      )
    }

    Sigma.mat.inv <- MatInv(mtx = GeoW$covmatrix)

    if (!is.matrix(Sigma.mat.inv) ||
        !identical(dim(Sigma.mat.inv), c(M0, M0)) ||
        any(!is.finite(Sigma.mat.inv))) {
      stop(
        "The latent covariance matrix could not be inverted reliably."
      )
    }
  }

  all_coord <- rbind(coord_obs, loc)

  ##########################################################################
  ## Batch size
  ##########################################################################

  fields_per_rep <- if (gaussian_case) 1L else 2L

  if (is.null(batch_size)) {

    if (length(max_batch_mb) != 1L ||
        !is.finite(max_batch_mb) ||
        max_batch_mb <= 0) {
      stop("max_batch_mb must be strictly positive.")
    }

    bytes_per_rep <- fields_per_rep * (M0 + N0) * 8

    batch_size <- floor(
      max_batch_mb * 1024^2 / max(bytes_per_rep, 1)
    )

    batch_size <- max(
      1L,
      min(nrep, as.integer(batch_size))
    )

  } else {

    if (length(batch_size) != 1L ||
        !is.finite(batch_size) ||
        batch_size < 1) {
      stop("batch_size must be at least 1.")
    }

    batch_size <- max(
      1L,
      min(nrep, as.integer(batch_size))
    )
  }

  ##########################################################################
  ## Normalize the output of GeoSim or GeoSimapprox
  ##########################################################################

  ##########################################################################
  ## Generate one batch of unconditional Gaussian simulations
  ##########################################################################

  generate_batch <- function(current_nrep) {

    current_fields <- fields_per_rep * current_nrep

    if (method %in% c("TB", "CE")) {

      sim_obj <- GeoSimapprox(
        coordx = all_coord,
        corrmodel = corrmodel,
        method = method,
        model = "Gaussian",
        param = latent_param,
        L = L,
        nrep = current_fields,
        progress = FALSE,
        parallel = FALSE,
        distance = distance,
        radius = radius,
        anisopars = anisopars
      )

    } else {

      sim_obj <- GeoSim(
        coordx = all_coord,
        corrmodel = corrmodel,
        model = "Gaussian",
        param = latent_param,
        nrep = current_fields,
        progress = FALSE,
        distance = distance,
        radius = radius,
        anisopars = anisopars
      )
    }

    sim_data <- .GeoSimcondNormalizeSimulations(
      sim_obj$data,
      current_fields
    )

    if (gaussian_case) {

      tasks <- lapply(
        seq_len(current_nrep),
        function(j) {
          list(
            g2 = as.numeric(sim_data[[j]])
          )
        }
      )

    } else {

      tasks <- lapply(
        seq_len(current_nrep),
        function(j) {
          list(
            g1 = as.numeric(sim_data[[2L * j - 1L]]),
            g2 = as.numeric(sim_data[[2L * j]])
          )
        }
      )
    }

    rm(sim_obj, sim_data)

    tasks
  }

  ##########################################################################
  ## Self-contained worker
  ##
  ## Only the Gaussian simulations required by the current task are passed
  ## to the worker. The complete batch is not captured as a global object.
  ##########################################################################

  worker_fun <- function(
      task, z_obs, M0, N0, gamma, mm, skew, sigma_sqrt,
      Sigma_mat_inv, Wglob, weights_orientation,
      n_iter, gaussian_case
  ) {

    tryCatch({

      apply_weights <- function(delta) {

        if (weights_orientation == "nobs_by_nloc") {
          as.numeric(crossprod(Wglob, delta))
        } else {
          as.numeric(Wglob %*% delta)
        }
      }

      conditional_field <- function(
          simulation,
          latent_obs
      ) {

        if (length(simulation) != M0 + N0) {
          stop(
            "An unconditional Gaussian simulation ",
            "has an incorrect length."
          )
        }

        sim_obs <- simulation[seq_len(M0)]
        sim_loc <- simulation[M0 + seq_len(N0)]

        sim_loc + apply_weights(
          latent_obs - sim_obs
        )
      }

      ######################################################################
      ## Gaussian special case
      ######################################################################

      if (gaussian_case) {

        g2_loc <- conditional_field(
          task$g2,
          z_obs
        )

        return(
          list(
            ok = TRUE,
            value = mm + sigma_sqrt * g2_loc
          )
        )
      }

      ######################################################################
      ## Initial state satisfying
      ##
      ## z_i = gamma |G1_i| + G2_i
      ######################################################################
      g1_init <- stats::rnorm(M0)
      g2_init <- z_obs - gamma * abs(g1_init)
      ######################################################################
      ## Exact Gibbs kernel implemented in C
      ######################################################################
      gibbs_result <- dotCall64::.C64(
        "skew_gaussian_gibbs_sampler",
        SIGNATURE = c(
          "double", "integer", "double", "double",
          "integer", "double", "double"
        ),
        data_obs = as.double(z_obs),
        n = as.integer(M0),
        Sigma_mat_inv = as.double(Sigma_mat_inv),
        eta = as.double(c(gamma, 0)),
        n_iter = as.integer(n_iter),
        data_x = as.double(g1_init),
        data_y = as.double(g2_init),
        INTENT = c(
          "r", "r", "r", "r",
          "r", "rw", "rw"
        ),
        VERBOSE = 0,
        NAOK = FALSE,
        PACKAGE = "GeoModels"
      )

      g1_obs <- as.numeric(gibbs_result$data_x)
      g2_obs <- as.numeric(gibbs_result$data_y)

      if (length(g1_obs) != M0 ||
          length(g2_obs) != M0 ||
          any(!is.finite(g1_obs)) ||
          any(!is.finite(g2_obs))) {
        stop(
          "The Gibbs sampler returned invalid ",
          "latent values."
        )
      }

      ######################################################################
      ## Check the nonlinear restriction
      ######################################################################
      constraint_error <- max(
        abs(
          z_obs -
            gamma * abs(g1_obs) -
            g2_obs
        )
      )
      constraint_tolerance <- 1e-7 * (
        1 + max(abs(z_obs))
      )

      if (!is.finite(constraint_error) ||
          constraint_error > constraint_tolerance) {
        stop(
          "The Gibbs output does not satisfy the ",
          "skew-Gaussian constraint."
        )
      }
      ######################################################################
      ## Gaussian substitution conditional simulations
      ######################################################################

      g1_loc <- conditional_field(
        task$g1,
        g1_obs
      )

      g2_loc <- conditional_field(
        task$g2,
        g2_obs
      )

      ######################################################################
      ## Back transformation
      ##
      ## Y(s) = mm + skew |G1(s)| + sqrt(sill) G2(s)
      ######################################################################

      list(
        ok = TRUE,
        value =
          mm +
          skew * abs(g1_loc) +
          sigma_sqrt * g2_loc
      )

    }, error = function(e) {

      list(
        ok = FALSE,
        error = conditionMessage(e)
      )
    })
  }

  ##########################################################################
  ## Parallel configuration
  ##########################################################################

  ncores <- .GeoResolveWorkers(parallel, ncores, n_jobs = nrep)
  parallel <- isTRUE(parallel) && nrep > 1L && ncores > 1L
  if (parallel) {
    old_plan <- future::plan()

    on.exit(
      future::plan(old_plan),
      add = TRUE
    )

    old_limit <- getOption(
      "future.globals.maxSize"
    )

    static_size <- sum(
      as.numeric(object.size(Sigma.mat.inv)),
      as.numeric(object.size(Wglob)),
      as.numeric(object.size(z_obs))
    )

    required_limit <- max(
      500 * 1024^2,
      2 * static_size + max_batch_mb * 1024^2
    )

    if (!is.null(old_limit)) {
      required_limit <- max(
        required_limit,
        old_limit
      )
    }

    options(
      future.globals.maxSize = required_limit
    )
    on.exit(
      options(
        future.globals.maxSize = old_limit
      ),
      add = TRUE
    )
    future::plan(
      future::multisession,
      workers = ncores
    )
  }
  ##########################################################################
  ## Batch execution
  ##########################################################################

  run_batches <- function(progressor) {

    res <- vector("list", nrep)
    failure_messages <- character(0L)
    batch_starts <- seq.int(
      1L,
      nrep,
      by = batch_size
    )
    for (batch_start in batch_starts) {

      batch_end <- min(
        nrep,
        batch_start + batch_size - 1L
      )
      current_indices <- batch_start:batch_end
      current_nrep <- length(current_indices)
      tasks <- generate_batch(current_nrep)
      if (parallel && current_nrep > 1L) {

        batch_result <- future.apply::future_lapply(
          tasks,
          worker_fun,
          z_obs = z_obs,
          M0 = M0,
          N0 = N0,
          gamma = gamma,
          mm = mm,
          skew = skew,
          sigma_sqrt = sigma_sqrt,
          Sigma_mat_inv = Sigma.mat.inv,
          Wglob = Wglob,
          weights_orientation = weights_orientation,
          n_iter = n_iter,
          gaussian_case = gaussian_case,
          future.seed = TRUE,

          # One future chunk per worker. This avoids repeatedly
          # serializing Sigma.mat.inv and Wglob for every simulation.
          future.scheduling = 1,

          # worker_fun is self-contained and all data are passed
          # explicitly.
          future.globals = FALSE,

          future.packages = c(
            "GeoModels",
            "dotCall64"
          )
        )

      } else {

        batch_result <- lapply(
          tasks,
          worker_fun,
          z_obs = z_obs,
          M0 = M0,
          N0 = N0,
          gamma = gamma,
          mm = mm,
          skew = skew,
          sigma_sqrt = sigma_sqrt,
          Sigma_mat_inv = Sigma.mat.inv,
          Wglob = Wglob,
          weights_orientation = weights_orientation,
          n_iter = n_iter,
          gaussian_case = gaussian_case
        )
      }

      for (j in seq_along(batch_result)) {

        global_index <- current_indices[j]

        if (isTRUE(batch_result[[j]]$ok)) {

          res[[global_index]] <-
            batch_result[[j]]$value

        } else {

          failure_messages <- c(
            failure_messages,
            sprintf(
              "simulation %d: %s",
              global_index,
              batch_result[[j]]$error
            )
          )
        }

        progressor()
      }

      rm(tasks, batch_result)
      gc(verbose = FALSE)
    }

    ########################################################################
    ## Remove failed simulations
    ########################################################################

    failed <- vapply(
      res,
      is.null,
      logical(1L)
    )
    if (any(failed)) {
      details <- if (length(failure_messages)) {
        paste(
          utils::head(failure_messages, 3L),
          collapse = "; "
        )
      } else {
        "unknown error"
      }
      warning(
        sprintf(
          "%d simulations out of %d failed. %s",
          sum(failed),
          nrep,
          details
        ),
        call. = FALSE
      )
      res <- res[!failed]
    }
    if (!length(res)) {
      stop("All conditional simulations failed.")
    }

    res
  }
  ##########################################################################
  ## Run
  ##########################################################################
  cat(
    "Performing",
    nrep,
    "conditional simulations",
    if (parallel) {
      paste0(" using ", ncores, " cores")
    } else {
      ""
    },
    "...\n"
  )
  if (progress) {
    restore_progress_handlers <- .GeoProgressHandlersPush(
      handler = "txtprogressbar", global = TRUE
    )
    on.exit(restore_progress_handlers(), add = TRUE)
    res <- progressr::with_progress({
      p <- progressr::progressor(
        steps = nrep
      )
      run_batches(p)
    })
  } else {
    res <- run_batches(
      function(...) invisible(NULL)
    )
  }
  res
}
#################################################################################
### Helper function for conditional Clayton-like copula simulation
#################################################################################
.GeoSimcondClaytonCopula <- function(
    coord_obs, loc, u_obs, nu, latent_param, corrmodel,
    nrep = 1, method = "Cholesky", L = 1000,
    distance = "Eucl", radius = 1, anisopars = NULL,
    parallel = FALSE, ncores = 1, progress = FALSE,
    n_iter = 1000L, max_batch_mb = 256
) {
  coord_obs <- as.matrix(coord_obs)
  loc <- as.matrix(loc)
  u_obs <- as.numeric(u_obs)
  M0 <- nrow(coord_obs)
  N0 <- nrow(loc)

  nu_raw <- as.numeric(nu)
  if (length(nu_raw) != 1L || !is.finite(nu_raw) || nu_raw <= 0 ||
      abs(nu_raw - round(nu_raw)) > sqrt(.Machine$double.eps)) {
    stop("For copula='Clayton', param$nu must be a positive integer.", call. = FALSE)
  }
  nu <- as.integer(round(nu_raw))
  d_latent <- nu + 2L

  if (length(u_obs) != M0 || any(!is.finite(u_obs)) ||
      any(u_obs <= 0 | u_obs >= 1)) {
    stop("Clayton latent uniforms must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (!(method %in% c("Cholesky", "TB"))) {
    stop("Clayton copula conditional simulation supports method='Cholesky' or method='TB'.",
         call. = FALSE)
  }

  ## U_nu = R^(nu/2), where R=A/(A+B).
  ratio_obs <- exp((2 / nu) * log(u_obs))
  eps_r <- .Machine$double.eps
  ratio_obs <- pmin(pmax(ratio_obs, eps_r), 1 - eps_r)

  cat("Computing Clayton latent Kriging weights ...\n")
  GeoW <- GeoKrigWeights(
    coordx = coord_obs, corrmodel = corrmodel, loc = loc,
    model = "Gaussian", distance = distance, radius = radius,
    anisopars = anisopars, param = latent_param
  )
  if (is.null(GeoW$covmatrix) || is.null(GeoW$weights)) {
    stop("GeoKrigWeights did not return the covariance matrix and weights required by the Clayton sampler.",
         call. = FALSE)
  }

  Sigma_inv <- MatInv(mtx = GeoW$covmatrix)
  if (!is.matrix(Sigma_inv) || !identical(dim(Sigma_inv), c(M0, M0)) ||
      any(!is.finite(Sigma_inv))) {
    stop("The Clayton latent covariance matrix could not be inverted reliably.", call. = FALSE)
  }

  Wglob <- GeoW$weights
  if (is.null(dim(Wglob))) Wglob <- matrix(Wglob, nrow = M0)
  weights_orientation <- GeoW$weights_orientation
  if (is.null(weights_orientation) ||
      !(weights_orientation %in% c("nobs_by_nloc", "nloc_by_nobs"))) {
    stop("GeoKrigWeights returned an unknown weights orientation.", call. = FALSE)
  }
  if (weights_orientation == "nobs_by_nloc" &&
      (nrow(Wglob) != M0 || ncol(Wglob) != N0)) {
    stop("The nobs_by_nloc Clayton latent weights have incompatible dimensions.",
         call. = FALSE)
  }
  if (weights_orientation == "nloc_by_nobs" &&
      (nrow(Wglob) != N0 || ncol(Wglob) != M0)) {
    stop("The nloc_by_nobs Clayton latent weights have incompatible dimensions.",
         call. = FALSE)
  }

  apply_weights <- function(delta) {
    if (weights_orientation == "nobs_by_nloc") {
      as.numeric(crossprod(Wglob, delta))
    } else {
      as.numeric(Wglob %*% delta)
    }
  }

  initialize_latent <- function() {
    U <- matrix(stats::rnorm(M0 * d_latent), nrow = M0, ncol = d_latent)
    for (ii in seq_len(M0)) {
      q <- sqrt(sum(U[ii, ]^2))
      if (!is.finite(q) || q <= 0) q <- sqrt(d_latent)

      za <- U[ii, seq_len(nu)]
      nz <- sqrt(sum(za^2))
      if (!is.finite(nz) || nz <= 0) {
        za <- rep(0, nu); za[1L] <- 1; nz <- 1
      }
      wb <- U[ii, nu + seq_len(2L)]
      nw <- sqrt(sum(wb^2))
      if (!is.finite(nw) || nw <= 0) {
        wb <- c(1, 0); nw <- 1
      }
      U[ii, seq_len(nu)] <- q * sqrt(ratio_obs[ii]) * za / nz
      U[ii, nu + seq_len(2L)] <- q * sqrt(1 - ratio_obs[ii]) * wb / nw
    }
    U
  }

  run_one <- function(field_list) {
    if (length(field_list) != d_latent)
      stop("Internal Clayton latent-field count mismatch.", call. = FALSE)

    U0 <- initialize_latent()
    gibbs <- dotCall64::.C64(
      "clayton_gibbs_sampler",
      SIGNATURE = c("double", "double", "integer", "integer", "double", "integer"),
      SigmaInv = as.double(Sigma_inv),
      ratio = as.double(ratio_obs),
      n = as.integer(M0),
      nu = as.integer(nu),
      U = as.double(U0),
      nIter = as.integer(n_iter),
      INTENT = c("r", "r", "r", "r", "rw", "r"),
      VERBOSE = 0, NAOK = FALSE, PACKAGE = "GeoModels"
    )
    Uobs <- matrix(gibbs$U, nrow = M0, ncol = d_latent)
    if (any(!is.finite(Uobs)))
      stop("The Clayton latent Gibbs sampler returned non-finite values.", call. = FALSE)

    ## The C sampler enforces A/(A+B)=ratio_obs at every Gibbs update and
    ## performs a final round-off projection onto that manifold.  Recompute the
    ## ratio with scaled row norms rather than raw sums of squares: this is more
    ## stable when one of the two latent radii is very small.
    Zobs <- Uobs[, seq_len(nu), drop = FALSE]
    Wobs <- Uobs[, nu + seq_len(2L), drop = FALSE]
    znorm <- sqrt(rowSums(Zobs^2))
    wnorm <- sqrt(rowSums(Wobs^2))
    qnorm <- sqrt(znorm^2 + wnorm^2)
    ratio_check <- (znorm / qnorm)^2
    ratio_err <- abs(ratio_check - ratio_obs)
    if (any(!is.finite(ratio_check)) || any(!is.finite(ratio_err))) {
      stop("The Clayton latent Gibbs sampler returned an invalid observed ratio constraint.",
           call. = FALSE)
    }
    max_ratio_err <- max(ratio_err)
    ## A discrepancy above this level is no longer attributable to ordinary
    ## floating-point round-off and should still fail loudly.
    if (max_ratio_err > 1e-7) {
      stop(sprintf(
        "The Clayton latent Gibbs sampler violated the observed ratio constraint (max abs error %.3g).",
        max_ratio_err
      ), call. = FALSE)
    }

    Upred <- matrix(0, nrow = N0, ncol = d_latent)
    for (jj in seq_len(d_latent)) {
      simj <- as.numeric(field_list[[jj]])
      if (length(simj) != M0 + N0)
        stop("A Clayton latent Gaussian simulation has an invalid length.", call. = FALSE)
      sim_obs <- simj[seq_len(M0)]
      sim_loc <- simj[M0 + seq_len(N0)]
      Upred[, jj] <- sim_loc + apply_weights(Uobs[, jj] - sim_obs)
    }

    Apred <- rowSums(Upred[, seq_len(nu), drop = FALSE]^2)
    Bpred <- rowSums(Upred[, nu + seq_len(2L), drop = FALSE]^2)
    den <- Apred + Bpred
    if (any(!is.finite(den)) || any(den <= 0))
      stop("Invalid Clayton latent radius at a prediction location.", call. = FALSE)
    beta_pred <- pmin(pmax(Apred / den, 0), 1)
    u_pred <- beta_pred^(nu / 2)
    pmin(pmax(u_pred, .Machine$double.eps), 1 - .Machine$double.eps)
  }

  all_coord <- rbind(coord_obs, loc)
  bytes_per_rep <- d_latent * (M0 + N0) * 8
  batch_size <- floor(max_batch_mb * 1024^2 / max(bytes_per_rep, 1))
  batch_size <- max(1L, min(nrep, as.integer(batch_size)))
  n_batches <- ceiling(nrep / batch_size)
  res <- vector("list", nrep)

  if (isTRUE(parallel) && nrep > 1L) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = ncores)
  }

  cat("Performing", nrep, "Clayton conditional simulations",
      if (isTRUE(parallel) && nrep > 1L) paste0(" using ", ncores, " cores") else "",
      "...\n")

  out_pos <- 1L
  for (bb in seq_len(n_batches)) {
    nb <- min(batch_size, nrep - out_pos + 1L)
    nfields <- d_latent * nb

    if (method == "TB") {
      gs <- GeoSimapprox(
        coordx = all_coord, corrmodel = corrmodel, method = "TB",
        model = "Gaussian", param = latent_param, L = L,
        nrep = nfields, progress = FALSE, parallel = FALSE,
        distance = distance, radius = radius, anisopars = anisopars
      )
    } else {
      gs <- GeoSim(
        coordx = all_coord, corrmodel = corrmodel,
        model = "Gaussian", param = latent_param,
        nrep = nfields, progress = FALSE,
        distance = distance, radius = radius, anisopars = anisopars
      )
    }
    gdata <- .GeoSimcondNormalizeSimulations(gs$data, nfields, M0 + N0)
    groups <- split(seq_len(nfields), ceiling(seq_len(nfields) / d_latent))

    one_batch <- function(ii) run_one(gdata[groups[[ii]]])
    if (isTRUE(parallel) && nb > 1L) {
      batch_res <- future.apply::future_lapply(
        seq_len(nb), one_batch, future.seed = TRUE,
        future.packages = c("GeoModels", "dotCall64")
      )
    } else {
      batch_res <- lapply(seq_len(nb), one_batch)
    }

    res[out_pos:(out_pos + nb - 1L)] <- batch_res
    out_pos <- out_pos + nb
    if (isTRUE(progress))
      cat(sprintf("  Clayton batch %d/%d completed\n", bb, n_batches))
    rm(gs, gdata, batch_res)
  }

  res
}

################################################################################# 
### Helper function for conditional Gamma simulation 
################################################################################# 
################################################################################# 
.GeoSimcondGamma <- function(coord_obs, loc, data, param, corrmodel,
                         nrep = 1, method = "Cholesky", local = FALSE, neighb = NULL,
                         L = NULL, distance = "Eucl", radius = 1, anisopars = NULL,
                         parallel = FALSE, ncores = 1, progress = FALSE,
                         mean_obs = NULL, mean_loc = NULL) {
  shape_raw <- as.numeric(param$shape)
  if (length(shape_raw) != 1L || !is.finite(shape_raw) || shape_raw <= 0 ||
      abs(shape_raw - round(shape_raw)) > sqrt(.Machine$double.eps)) {
    stop("Direct Gamma conditional simulation requires param$shape to be a positive integer.",
         call. = FALSE)
  }
  shape <- as.integer(round(shape_raw))

  keep <- CorrParam(corrmodel)
  corr_param <- param[keep]
  latent_nugget <- if (is.null(param$nugget)) 0 else as.numeric(param$nugget)
  if (length(latent_nugget) != 1L || !is.finite(latent_nugget) ||
      latent_nugget < 0 || latent_nugget >= 1) {
    stop("param$nugget must be a finite scalar in [0, 1).")
  }
  latent_param <- c(corr_param, list(mean = 0, nugget = latent_nugget, sill = 1))
  
  cat("Computing Kriging weights ...\n")
  if (!local) 
    GeoW <- GeoKrigWeights(coordx = coord_obs, corrmodel = corrmodel,
                           loc = loc, model = "Gaussian",
                           distance = distance, radius = radius,
                           anisopars = anisopars,
                           param = latent_param)
  else stop("Local option cannot be used for Gibbs sampling.\n")
  
  # Static quantities used by all Gibbs iterations.
  Sigma.mat.inv <- MatInv(mtx = GeoW$covmatrix)
  M0 <- nrow(coord_obs)
  N0 <- nrow(loc)
  Wglob <- GeoW$weights
  weights_orientation <- GeoW$weights_orientation
  if (is.null(weights_orientation) ||
      !(weights_orientation %in% c("nobs_by_nloc", "nloc_by_nobs"))) {
    stop("GeoKrigWeights returned an unknown weights orientation.")
  }
  if (!is.matrix(Sigma.mat.inv) ||
      !identical(dim(Sigma.mat.inv), c(M0, M0)) ||
      any(!is.finite(Sigma.mat.inv))) {
    stop("The latent covariance matrix could not be inverted reliably.")
  }
  if (weights_orientation == "nobs_by_nloc" &&
      (nrow(Wglob) != M0 || ncol(Wglob) != N0)) {
    stop("The nobs_by_nloc kriging weights have incompatible dimensions.")
  }
  if (weights_orientation == "nloc_by_nobs" &&
      (nrow(Wglob) != N0 || ncol(Wglob) != M0)) {
    stop("The nloc_by_nobs kriging weights have incompatible dimensions.")
  }
  if(is.null(mean_obs)) mean_obs <- rep(as.numeric(param$mean)[1L], M0)
  if(is.null(mean_loc)) mean_loc <- rep(as.numeric(param$mean)[1L], N0)
  mean_obs <- as.numeric(mean_obs)
  mean_loc <- as.numeric(mean_loc)
  if(length(mean_obs) != M0 || any(!is.finite(mean_obs)))
    stop("mean_obs must contain one finite value per observation.")
  if(length(mean_loc) != N0 || any(!is.finite(mean_loc)))
    stop("mean_loc must contain one finite value per prediction point.")
  data_transformed <- exp(-mean_obs) * data
  if (any(!is.finite(data_transformed)) || any(data_transformed < 0)) {
    stop("Gamma conditional simulation requires finite non-negative observations.")
  }
  # gamma_gibbs_sampler constrains ||U_i||^2 = 2 * y_i.  For the
  # GeoSim Gamma construction, ||U_i||^2 must equal shape * Y_i / exp(eta_i).
  gamma_constraint <- 0.5 * shape * data_transformed
  exp_mean_over_shape <- exp(mean_loc) / shape

  # STRATEGIA MEMORY-EFFICIENT: Genera simulazioni in batch
  if (parallel && nrep > 1) {
    
    # Calcola memoria necessaria per un batch
    coord_sim <- rbind(coord_obs, loc)
    n_total <- nrow(coord_sim)
    
    # Stima: ogni simulazione ~ 8 bytes * n_total
    memory_per_sim_mb <- (8 * n_total) / (1024^2)
    memory_per_batch_mb <- memory_per_sim_mb * shape
    
    # Determina batch size: massimo 200 MB per batch
    max_batch_memory_mb <- 200
    sims_per_batch <- max(1, floor(max_batch_memory_mb / memory_per_batch_mb))
    
    #cat(sprintf("Memory-efficient mode: processing %d simulations per batch\n", sims_per_batch))
    #cat(sprintf("Total batches needed: %d\n", ceiling(nrep / sims_per_batch)))
    
    # Setup progress
    if (progress) {
      restore_progress_handlers <- .GeoProgressHandlersPush(
        handler = "txtprogressbar", global = TRUE
      )
      on.exit(restore_progress_handlers(), add = TRUE)
      pb <- progressr::progressor(
        along = seq_len(nrep), on_exit = FALSE, auto_finish = FALSE
      )
      ## Finish before restoring the temporary global handler, including when
      ## a batch or Gibbs simulation exits early.
      on.exit(pb(type = "finish"), add = TRUE, after = FALSE)
    }
    
    # Parallelization by batch. Preserve the user's future plan.
    old_plan_batch <- future::plan()
    on.exit(try(future::plan(old_plan_batch), silent = TRUE), add = TRUE)
    future::plan(future::multisession, workers = ncores)

    res <- vector("list", nrep)
    n_batches <- ceiling(nrep / sims_per_batch)
    
    for (batch_idx in seq_len(n_batches)) {
      start_idx <- (batch_idx - 1) * sims_per_batch + 1
      end_idx <- min(batch_idx * sims_per_batch, nrep)
      batch_nrep <- end_idx - start_idx + 1
      
      cat(sprintf("\n=== Batch %d/%d: simulations %d to %d ===\n", 
                  batch_idx, n_batches, start_idx, end_idx))
      
      # Genera simulazioni Gaussiane SOLO per questo batch
      if (method == "TB" || method == "CE") {
        Gauss_batch <- GeoSimapprox(
          coordx = coord_sim, 
          corrmodel = corrmodel, 
          method = method, 
          model = "Gaussian", 
          param = latent_param, 
          L = L, 
          nrep = shape * batch_nrep, 
          progress = FALSE, 
          parallel = FALSE,
          distance = distance,
          radius = radius,
          anisopars = anisopars
        )
      } else {
        Gauss_batch <- GeoSim(
          coordx = coord_sim, 
          corrmodel = corrmodel, 
          model = "Gaussian", 
          param = latent_param, 
          nrep = shape * batch_nrep, 
          progress = FALSE,
          distance = distance,
          radius = radius,
          anisopars = anisopars
        )
      }
      
      pair_idx_batch <- split(seq_len(shape * batch_nrep), 
                             ceiling(seq_len(shape * batch_nrep) / shape))
      
      # Normalize also the single-field case returned as a vector.
      Gauss_data_batch <- .GeoSimcondNormalizeSimulations(
        Gauss_batch$data, shape * batch_nrep, M0 + N0
      )
      
      # Variabili essenziali per questo batch
      essential_vars <- list(
        shape = shape,
        data_transformed = data_transformed,
        gamma_constraint = gamma_constraint,
        M0 = M0,
        N0 = N0,
        Sigma.mat.inv = Sigma.mat.inv,
        Wglob = Wglob,
        weights_orientation = weights_orientation,
        exp_mean_over_shape = exp_mean_over_shape,
        pair_idx = pair_idx_batch
      )
      
      # Worker function ottimizzata
      run_simulation_worker <- function(local_i, vars, Gauss_data) {
        shape <- vars$shape
        data_transformed <- vars$data_transformed
        gamma_constraint <- vars$gamma_constraint
        M0 <- vars$M0
        N0 <- vars$N0
        Sigma.mat.inv <- vars$Sigma.mat.inv
        Wglob <- vars$Wglob
        weights_orientation <- vars$weights_orientation
        exp_mean_over_shape <- vars$exp_mean_over_shape
        pair_idx <- vars$pair_idx
        
        # Gibbs sampling
        random_sign <- matrix(ifelse(runif(M0 * shape) <= 0.5, 1, -1), 
                             nrow = M0, ncol = shape)
        ini_vals <- random_sign * matrix(rep(sqrt(data_transformed), shape),
                                         nrow = M0, ncol = shape)
        
        gibbs.result <- tryCatch({
          dotCall64::.C64("gamma_gibbs_sampler",
                          SIGNATURE = c("double", "double", "int", "int", 
                                       "double", "int", "int"),
                          SigmaInv = as.double(Sigma.mat.inv),
                          y        = as.double(gamma_constraint),
                          n        = as.integer(M0),
                          v        = as.integer(shape),
                          U        = as.double(ini_vals),
                          nIte     = as.integer(1000),
                          nRep     = as.integer(150),
                          INTENT   = c("r", "r", "r", "r", "rw", "r", "r"), 
                          VERBOSE = 0, 
                          NAOK = FALSE, 
                          PACKAGE = "GeoModels")
        }, error = function(e) {
          return(list(U = NULL))
        })
        
        if (is.null(gibbs.result$U)) return(NULL)
        
        gibbs_U <- matrix(gibbs.result$U, nrow = M0, ncol = shape)
        
        # Processamento
        idx_xy <- pair_idx[[local_i]]
        res_i <- matrix(0, nrow = N0, ncol = shape)
        
        for (kk in seq_len(shape)) {
          simX <- Gauss_data[[idx_xy[kk]]]
          sim_train_x <- simX[seq_len(M0)]
          sim_valid_x <- simX[(M0 + 1):(N0 + M0)]
          delta <- gibbs_U[, kk] - sim_train_x
          sk_pred_x <- if (weights_orientation == "nobs_by_nloc") {
            as.numeric(crossprod(Wglob, delta))
          } else {
            as.numeric(Wglob %*% delta)
          }
          res_i[, kk] <- sim_valid_x + sk_pred_x
        }
        
        return(exp_mean_over_shape * rowSums(res_i^2))
      }
      
      # Parallelizzazione del batch corrente
      old_limit <- options(future.globals.maxSize = 2000 * 1024^2)
      on.exit(options(old_limit), add = TRUE)
      
      batch_results <- future.apply::future_lapply(
        seq_len(batch_nrep),
        function(local_i) {
          result <- run_simulation_worker(local_i, essential_vars, Gauss_data_batch)
          return(result)
        },
        future.seed = TRUE,
        future.packages = c("dotCall64", "GeoModels"),
        future.globals = structure(TRUE, add = c("essential_vars", "Gauss_data_batch"))
      )
      
      # Salva risultati del batch
      for (local_i in seq_len(batch_nrep)) {
        res[[start_idx + local_i - 1]] <- batch_results[[local_i]]
        if (progress) {
          pb(sprintf("Simulation %d/%d completed", start_idx + local_i - 1, nrep))
        }
      }
      
      # Pulizia memoria dopo ogni batch
      rm(Gauss_batch, Gauss_data_batch, batch_results, essential_vars)
      gc()
    }
    
  } else {
    # Modalità sequenziale (codice originale)
    if (method == "TB" || method == "CE") {
      Gauss.all <- GeoSimapprox(
        coordx = rbind(coord_obs, loc), 
        corrmodel = corrmodel, 
        method = method, 
        model = "Gaussian", 
        param = latent_param, 
        L = L, 
        nrep = shape * nrep, 
        progress = FALSE, 
        parallel = FALSE,
        distance = distance,
        radius = radius,
        anisopars = anisopars
      )
    } else {
      Gauss.all <- GeoSim(
        coordx = rbind(coord_obs, loc), 
        corrmodel = corrmodel, 
        model = "Gaussian", 
        param = latent_param, 
        nrep = shape * nrep, 
        progress = FALSE,
        distance = distance,
        radius = radius,
        anisopars = anisopars
      )
    }
    
    Gauss_data <- .GeoSimcondNormalizeSimulations(
      Gauss.all$data, shape * nrep, M0 + N0
    )
    pair_idx <- split(seq_len(shape * nrep), ceiling(seq_len(shape * nrep) / shape))
    
    cat("Performing", nrep, "conditional simulations ...\n")
    if (progress) {
      restore_progress_handlers <- .GeoProgressHandlersPush(
        handler = "txtprogressbar", global = TRUE
      )
      on.exit(restore_progress_handlers(), add = TRUE)
      pb <- progressr::progressor(
        along = seq_len(nrep), on_exit = FALSE, auto_finish = FALSE
      )
      ## Finish before restoring the temporary global handler, including when
      ## a batch or Gibbs simulation exits early.
      on.exit(pb(type = "finish"), add = TRUE, after = FALSE)
    }
    
    res <- vector("list", nrep)
    
    for (i in seq_len(nrep)) {
      random_sign <- matrix(ifelse(runif(M0 * shape) <= 0.5, 1, -1), 
                           nrow = M0, ncol = shape)
      ini_vals <- random_sign * matrix(rep(sqrt(data_transformed), shape),
                                         nrow = M0, ncol = shape)
      
      gibbs.result <- dotCall64::.C64("gamma_gibbs_sampler",
                                      SIGNATURE = c("double", "double", "int", "int", 
                                                   "double", "int", "int"),
                                      SigmaInv = as.double(Sigma.mat.inv),
                                      y        = as.double(gamma_constraint),
                                      n        = as.integer(M0),
                                      v        = as.integer(shape),
                                      U        = as.double(ini_vals),
                                      nIte     = as.integer(1000),
                                      nRep     = as.integer(150),
                                      INTENT   = c("r", "r", "r", "r", "rw", "r", "r"), 
                                      VERBOSE = 0, 
                                      NAOK = FALSE, 
                                      PACKAGE = "GeoModels")
      
      if (is.null(gibbs.result$U)) {
        warning(paste("Gibbs sampler failed at simulation", i))
        next
      }
      
      gibbs_U <- matrix(gibbs.result$U, nrow = M0, ncol = shape)
      
      idx_xy <- pair_idx[[i]]
      res_i <- matrix(0, nrow = N0, ncol = shape)
      
      for (kk in seq_len(shape)) {
        simX <- Gauss_data[[idx_xy[kk]]]
        sim_train_x <- simX[seq_len(M0)]
        sim_valid_x <- simX[(M0 + 1):(N0 + M0)]
        delta <- gibbs_U[, kk] - sim_train_x
        sk_pred_x <- if (weights_orientation == "nobs_by_nloc") {
          as.numeric(crossprod(Wglob, delta))
        } else {
          as.numeric(Wglob %*% delta)
        }
        res_i[, kk] <- sim_valid_x + sk_pred_x
      }
      
      res[[i]] <- exp_mean_over_shape * rowSums(res_i^2)
      
      if (progress) {
        pb(sprintf("Simulation %d/%d completed", i, nrep))
      }
    }
  }
  
  # Rimuove risultati nulli
  res <- res[!sapply(res, is.null)]
  failed_sims <- nrep - length(res)
  if (failed_sims > 0) {
    warning(paste(failed_sims, "simulazioni su", nrep, "sono fallite"))
  }
  
  return(res)
}
################################################################################# 
### internal function: conditional simulation of Gaussian RF 
################################################################################# 
.GeoSimcondPrepareLocalWeightMap <- function(weights_list, n_obs, n_loc) {
  if (!is.list(weights_list) || length(weights_list) != n_loc) {
    stop("GeoKriglocWeights returned an incompatible local weights list.")
  }

  lens <- integer(n_loc)
  idx_parts <- vector("list", n_loc)
  weight_parts <- vector("list", n_loc)

  for (j in seq_len(n_loc)) {
    wj <- weights_list[[j]]
    if (is.null(wj)) next

    idx <- as.integer(wj$neighbor_indices)
    ww <- as.numeric(wj$weights)
    if (length(idx) != length(ww)) {
      stop("A local kriging system returned inconsistent indices and weights.")
    }
    if (!length(idx)) next
    if (anyNA(idx) || any(idx < 1L | idx > n_obs)) {
      stop("A local kriging system returned invalid observation indices.")
    }
    if (any(!is.finite(ww))) {
      stop("A local kriging system returned non-finite weights.")
    }

    lens[j] <- length(idx)
    idx_parts[[j]] <- idx
    weight_parts[[j]] <- ww
  }

  active <- which(lens > 0L)
  if (!length(active)) {
    return(list(
      index = integer(0L), weight = numeric(0L), target = integer(0L),
      active = integer(0L), nnz = 0L, sparse_weights = NULL
    ))
  }

  out <- list(
    index = as.integer(unlist(idx_parts[active], use.names = FALSE)),
    weight = as.numeric(unlist(weight_parts[active], use.names = FALSE)),
    target = rep(as.integer(active), times = lens[active]),
    active = as.integer(active),
    nnz = as.integer(sum(lens)),
    sparse_weights = NULL
  )

  ## spam is already a GeoModels dependency.  Build the n_loc x n_obs local
  ## weight operator once so that multiple conditional Gaussian replicates can
  ## share compiled sparse matrix multiplication.  Keep the compact triplet map
  ## as a fallback for compatibility with older spam installations.
  out$sparse_weights <- tryCatch({
    sp <- spam::spam(
      list(i = out$target, j = out$index, values = out$weight),
      nrow = n_loc, ncol = n_obs
    )
    ## Validate the sparse multiplication method before discarding the triplet
    ## representation.  If this probe fails, retain the compact base-R fallback.
    probe <- as.numeric(sp %*% numeric(n_obs))
    if (length(probe) != n_loc || any(!is.finite(probe))) NULL else sp
  }, error = function(e) NULL)

  if (!is.null(out$sparse_weights)) {
    ## Avoid storing both sparse and triplet representations for massive local
    ## graphs.  nnz is retained for the empty/non-empty decision.
    out$index <- integer(0L)
    out$weight <- numeric(0L)
    out$target <- integer(0L)
    out$active <- integer(0L)
  }

  out
}

.GeoSimcondLocalWeightProduct <- function(x, weight_map, n_loc) {
  out <- numeric(n_loc)
  if (is.null(weight_map$nnz) || weight_map$nnz == 0L) return(out)

  if (!is.null(weight_map$sparse_weights)) {
    ans <- as.numeric(weight_map$sparse_weights %*% as.numeric(x))
    if (length(ans) != n_loc || any(!is.finite(ans))) {
      stop("Sparse local conditioning returned invalid values.")
    }
    return(ans)
  }

  weighted <- weight_map$weight * as.numeric(x)[weight_map$index]
  dim(weighted) <- c(length(weighted), 1L)
  sums <- rowsum(weighted, group = weight_map$target, reorder = FALSE)
  out[weight_map$active] <- as.numeric(sums[, 1L])
  out
}

.GeoSimcondLocalConditionalBatch <- function(sim_data_list, krig_pred_vec, weight_map,
                                    n_obs, n_loc, progressor = NULL,
                                    target_temp_mb = 128) {
  nrep <- length(sim_data_list)
  out <- vector("list", nrep)
  if (!nrep) return(out)
  if (is.null(progressor)) progressor <- function(...) NULL

  ## For one realization or an empty local graph the compact vector path is
  ## cheaper and preserves the maxdist-only no-neighbor semantics exactly.
  if (is.null(weight_map$nnz) || weight_map$nnz == 0L ||
      is.null(weight_map$sparse_weights) || nrep == 1L) {
    for (i in seq_len(nrep)) {
      out[[i]] <- .GeoSimcondOne(
        i, sim_data_list, krig_pred_vec, TRUE,
        weight_map, NA_character_, n_obs, n_loc
      )
      progressor()
    }
    return(out)
  }

  ## Bound only the temporary dense matrices created by the sparse batch.  The
  ## unconditional simulations are already retained by the simulator.  This
  ## avoids a large nrep * nnz temporary object when many realizations are asked
  ## for on a massive data set.
  target_bytes <- max(1, as.numeric(target_temp_mb)) * 1024^2
  bytes_per_rep <- 8 * (n_obs + 4 * n_loc)
  block_size <- max(
    1L,
    min(nrep, as.integer(floor(target_bytes / max(bytes_per_rep, 1))))
  )
  blocks <- split(seq_len(nrep), ceiling(seq_len(nrep) / block_size))

  sparse_ok <- TRUE
  for (ids in blocks) {
    nb <- length(ids)
    obs_mat <- matrix(
      unlist(lapply(ids, function(i) sim_data_list[[i]][seq_len(n_obs)]),
             use.names = FALSE),
      nrow = n_obs, ncol = nb
    )
    loc_mat <- matrix(
      unlist(lapply(ids, function(i) sim_data_list[[i]][n_obs + seq_len(n_loc)]),
             use.names = FALSE),
      nrow = n_loc, ncol = nb
    )

    obs_pred <- if (sparse_ok) {
      tryCatch(
        as.matrix(weight_map$sparse_weights %*% obs_mat),
        error = function(e) NULL
      )
    } else NULL

    if (is.null(obs_pred) || !identical(dim(obs_pred), c(n_loc, nb)) ||
        any(!is.finite(obs_pred))) {
      ## Defensive fallback: do not make local conditional simulation depend on
      ## sparse-matrix method details of a particular spam version.
      sparse_ok <- FALSE
      for (i in ids) {
        out[[i]] <- .GeoSimcondOne(
          i, sim_data_list, krig_pred_vec, TRUE,
          weight_map, NA_character_, n_obs, n_loc
        )
        progressor()
      }
      next
    }

    cond_mat <- loc_mat - obs_pred
    cond_mat <- cond_mat + krig_pred_vec
    out[ids] <- lapply(seq_len(nb), function(k) as.numeric(cond_mat[, k]))
    for (unused in ids) progressor()
  }

  out
}

.GeoSimcondGaussianMeans <- function(param, X, Xloc, Mloc, n_obs, n_loc) {
  external <- .GeoMean_external(param, n = n_obs)
  mu_obs <- .GeoMean_vector(
    param, X = X, n = n_obs, external = external, xname = "X"
  )

  if (!is.null(Mloc)) {
    mu_loc <- as.numeric(Mloc)
    if (length(mu_loc) == 1L) mu_loc <- rep(mu_loc, n_loc)
    if (length(mu_loc) != n_loc || any(!is.finite(mu_loc))) {
      stop("Mloc must contain one finite value per prediction point.")
    }
  } else if (!is.null(Xloc)) {
    Xloc_mat <- .GeoMean_design(Xloc, n = n_loc, name = "Xloc")
    beta <- .GeoMean_beta(param, ncol(Xloc_mat))
    mu_loc <- as.numeric(Xloc_mat %*% beta)
  } else {
    if (!is.null(X)) {
      stop("Xloc or Mloc is required when X is supplied.")
    }
    mu_loc <- .GeoMean_vector(param, X = NULL, n = n_loc, xname = "Xloc")
  }

  if (length(mu_obs) != n_obs || length(mu_loc) != n_loc ||
      any(!is.finite(c(mu_obs, mu_loc)))) {
    stop("Could not construct finite Gaussian means for local conditioning.")
  }

  list(obs = as.numeric(mu_obs), loc = as.numeric(mu_loc))
}

.GeoSimcondOne <- function(i, sim_data_list, krig_pred_vec, local_flag,
                             weights_obj, weights_orientation, n_obs, n_loc) {
  sim_i <- sim_data_list[[i]]
  sim_nc_obs_data <- sim_i[seq_len(n_obs)]
  sim_nc_loc_data <- sim_i[n_obs + seq_len(n_loc)]

  if (local_flag) {
    sim_nc_obs_pred <- .GeoSimcondLocalWeightProduct(sim_nc_obs_data, weights_obj, n_loc)
  } else {
    W <- weights_obj
    xvec <- as.vector(sim_nc_obs_data)

    if (identical(weights_orientation, "nobs_by_nloc")) {
      sim_nc_obs_pred <- as.vector(crossprod(W, xvec))
    } else if (identical(weights_orientation, "nloc_by_nobs")) {
      sim_nc_obs_pred <- as.vector(W %*% xvec)
    } else {
      stop("Unknown kriging weights orientation.")
    }
  }

  krig_pred_vec + sim_nc_loc_data - sim_nc_obs_pred
}

#################################################################################
### main function: conditional Gaussian simulations
#################################################################################
.GeoSimcondGaussian <- function(data, corrmodel, nrep, method, L,
                     param,
                     coord_obs, loc, coordt_use, time, X, Xloc, Mloc, distance, radius, anisopars,
                     local, neighb, maxdist, maxtime,
                     space, spacetime, bivariate, parallel, ncores, progress) {

  #############################################
  # Input validation
  #############################################
  if (is.null(loc) || nrow(loc) == 0) {
    stop("loc must be a non-empty matrix")
  }
  if (nrep < 1) {
    stop("nrep must be at least 1")
  }
  if (local && is.null(neighb) && is.null(maxdist)) {
    stop("For local conditional simulation, specify neighb or maxdist.",
         call. = FALSE)
  }

  #############################################
  # Build combined coordinates
  #############################################
  coord_sim <- rbind(coord_obs, loc)
  time_sim  <- c(coordt_use, time)
  n_obs <- length(data)
  n_loc <- if(!is.null(Mloc)) length(Mloc) else if(!is.null(Xloc)) nrow(Xloc) else nrow(loc)

  ## The substitution residual must be a zero-mean Gaussian field:
  ## Z_c = E(Z | data) + Z*loc - W Z*obs.
  param_sim <- .GeoMean_intercept_param(param, 0)
  X_sim <- NULL

  #############################################
  # Unconditional simulations
  #############################################
  if (method %in% c("TB", "CE") && isTRUE(parallel)) {
    cat("Performing", nrep, "unconditional simulations using", method,
        "with", ncores, "cores ...\n")
  } else {
    cat("Performing", nrep, "unconditional simulations using", method, "...\n")
  }

  if (method == "Cholesky") {
    sim_args <- list(
      coordx = coord_sim, coordt = time_sim, corrmodel = corrmodel,
      progress = FALSE, X = X_sim, nrep = nrep, distance = distance,
      radius = radius, anisopars = anisopars
    )
    sim_args <- c(sim_args, list(param = param_sim))
    ## GeoSim() prints its own Cholesky banner when nrep > 1.  Suppress that
    ## internal banner here so GeoSimcond reports this phase exactly once.
    invisible(utils::capture.output(
      sim_nc <- do.call(GeoSim, sim_args),
      type = "output"
    ))
  } else if (method == "TB" || method == "CE") {
    sim_args_approx <- list(
      coordx = coord_sim, coordt = time_sim, corrmodel = corrmodel,
      progress = FALSE, method = method, L = L, parallel = parallel,
      ncores = ncores,
      X = X_sim, nrep = nrep, distance = distance, radius = radius,
      anisopars = anisopars
    )
    sim_args_approx <- c(sim_args_approx, list(param = param_sim))
    sim_nc <- do.call(GeoSimapprox, sim_args_approx)
  } else {
    stop("Unsupported unconditional simulation method.")
  }

  sim_data_list <- .GeoSimcondNormalizeSimulations(
    sim_nc$data, nrep, n_obs + n_loc
  )
  rm(sim_nc)

  #############################################
  # Prediction and kriging weights
  #############################################
  weights_orientation <- NA_character_

  if (local) {
    local_desc <- character(0L)
    if (!is.null(neighb)) {
      local_desc <- c(local_desc, paste0(neighb, " nearest neighbours"))
    }
    if (!is.null(maxdist)) {
      local_desc <- c(local_desc, paste0("maxdist = ", format(maxdist)))
    }
    if (!is.null(maxtime)) {
      local_desc <- c(local_desc, paste0("maxtime = ", format(maxtime)))
    }
    local_desc <- if (length(local_desc)) paste(local_desc, collapse = "; ") else "local"
    weight_parallel_note <- if (isTRUE(parallel)) paste0("; ", ncores, " cores") else ""
    cat("Computing compact local kriging weights (", local_desc,
        weight_parallel_note, ") ...\n", sep = "")

    weights_args <- list(
      coordx = coord_obs, coordt = coordt_use,
      corrmodel = corrmodel, loc = loc,
      X = X, Xloc = Xloc, Mloc = Mloc, time = time,
      distance = distance, radius = radius, anisopars = anisopars,
      param = param, neighb = neighb, maxdist = maxdist,
      maxtime = maxtime, parallel = parallel, ncores = ncores,
      compact = TRUE
    )
    krig_weights_obj <- do.call(GeoKriglocWeights, weights_args)
    local_weights <- .GeoSimcondPrepareLocalWeightMap(
      krig_weights_obj$weights, n_obs = n_obs, n_loc = n_loc
    )
    rm(krig_weights_obj)

    means <- .GeoSimcondGaussianMeans(
      param = param, X = X, Xloc = Xloc, Mloc = Mloc,
      n_obs = n_obs, n_loc = n_loc
    )
    krig_pred_vec <- means$loc + .GeoSimcondLocalWeightProduct(
      data - means$obs, local_weights, n_loc
    )
    rm(means)

    weights_obj <- local_weights
    rm(local_weights)
  } else {
    cat("Computing global kriging predictor and weights ...\n")
    krig_sim_args <- list(
      coordx = coord_obs, coordt = coordt_use,
      data = data, corrmodel = corrmodel,
      loc = loc, X = X, Xloc = Xloc, Mloc = Mloc, time = time,
      distance = distance, radius = radius, anisopars = anisopars,
      param = param
    )
    krig_sim <- do.call(GeoKrig, krig_sim_args)
    krig_pred_vec <- as.numeric(krig_sim$pred)
    rm(krig_sim)

    weights_args <- list(
      coordx = coord_obs, coordt = coordt_use,
      corrmodel = corrmodel, loc = loc,
      X = X, Xloc = Xloc, Mloc = Mloc, time = time,
      distance = distance, radius = radius, anisopars = anisopars,
      param = param
    )
    krig_weights_obj <- do.call(GeoKrigWeights, weights_args)
    W <- krig_weights_obj$weights
    weights_orientation <- krig_weights_obj$weights_orientation
    if (is.null(weights_orientation) ||
        !(weights_orientation %in% c("nobs_by_nloc", "nloc_by_nobs"))) {
      stop("GeoKrigWeights returned an unknown weights orientation.")
    }
    if (weights_orientation == "nobs_by_nloc" &&
        (nrow(W) != n_obs || ncol(W) != n_loc)) {
      stop("The nobs_by_nloc kriging weights have incompatible dimensions.")
    }
    if (weights_orientation == "nloc_by_nobs" &&
        (nrow(W) != n_loc || ncol(W) != n_obs)) {
      stop("The nloc_by_nobs kriging weights have incompatible dimensions.")
    }
    rm(krig_weights_obj)
    weights_obj <- W
    rm(W)
  }

  if (length(krig_pred_vec) != n_loc || any(!is.finite(krig_pred_vec))) {
    stop("The kriging predictor returned invalid conditional means.")
  }

  #############################################
  # Conditional simulations with progress bar
  #############################################
  if (local) {
    cat("Applying local conditional correction to", nrep,
        "simulations ...\n")
  } else if (isTRUE(parallel) && nrep > 1L) {
    cat("Performing", nrep, "global conditional simulations in parallel using",
        ncores, "cores ...\n")
  } else {
    cat("Performing", nrep, "global conditional simulations ...\n")
  }

  if (progress) {
    restore_progress_handlers <- .GeoProgressHandlersPush(
      handler = "txtprogressbar", global = TRUE
    )
    on.exit(restore_progress_handlers(), add = TRUE)
    p <- progressr::progressor(
      along = seq_len(nrep), on_exit = FALSE, auto_finish = FALSE
    )
    ## Finish before restoring the temporary global handler on every exit path.
    on.exit(p(type = "finish"), add = TRUE, after = FALSE)
  } else {
    p <- function(...) NULL
  }

  if (space) {
    ## For local conditioning the expensive stages (TB/CE simulation and local
    ## kriging-weight construction) can already use the requested workers.
    ## Keep the final substitution in the main R process: with compact local
    ## weights it is a cheap vectorized operation, while a PSOCK cluster would
    ## duplicate sim_data_list and the local weight map on every worker.
    if (local) {
      ## Fast path shared by Gaussian fields, one-Gaussian monotone transforms
      ## and Gaussian-copula margins.  The local graph is identical for every
      ## realization, so apply the same sparse operator in bounded batches.
      sim_cond <- .GeoSimcondLocalConditionalBatch(
        sim_data_list = sim_data_list,
        krig_pred_vec = krig_pred_vec,
        weight_map = weights_obj,
        n_obs = n_obs,
        n_loc = n_loc,
        progressor = p
      )
    } else if (parallel && !is.null(ncores) && nrep > 1) {
      memory_needed <- object.size(sim_data_list) +
                       object.size(weights_obj) +
                       object.size(krig_pred_vec)
      memory_gb <- as.numeric(memory_needed) / 1024^3

      if (memory_gb > 10) {
        chunk_size <- max(1, min(20, floor(nrep / 4)))
        chunks <- split(seq_len(nrep), ceiling(seq_len(nrep) / chunk_size))
        sim_cond <- vector("list", nrep)

        for (chunk_idx in seq_along(chunks)) {
          chunk_indices <- chunks[[chunk_idx]]
          cl <- parallel::makeCluster(min(ncores, length(chunk_indices)))
          chunk_sim_data_only <- sim_data_list[chunk_indices]
          krig_pred_only <- krig_pred_vec
          weights_only <- weights_obj
          orientation_only <- weights_orientation

          compute_one <- .GeoSimcondOne
          chunk_results <- tryCatch({
            parallel::clusterExport(cl, "compute_one", envir = environment())
            parallel::parLapply(cl, seq_along(chunk_indices), function(local_i) {
              compute_one(
                local_i, chunk_sim_data_only, krig_pred_only, FALSE,
                weights_only, orientation_only, n_obs, n_loc
              )
            })
          }, finally = {
            parallel::stopCluster(cl)
          })

          for (local_i in seq_along(chunk_indices)) {
            sim_cond[[chunk_indices[local_i]]] <- chunk_results[[local_i]]
            p()
          }
          rm(chunk_sim_data_only, krig_pred_only, weights_only,
             orientation_only, chunk_results)
          if (chunk_idx %% 5 == 0) gc()
        }
      } else {
        cl <- parallel::makeCluster(ncores)
        compute_one <- .GeoSimcondOne
        sim_cond <- tryCatch({
          parallel::clusterExport(cl, "compute_one", envir = environment())
          parallel::parLapply(cl, seq_len(nrep), function(i) {
            compute_one(
              i, sim_data_list, krig_pred_vec, FALSE,
              weights_obj, weights_orientation, n_obs, n_loc
            )
          })
        }, finally = {
          parallel::stopCluster(cl)
        })

        for (i in seq_len(nrep)) p()
      }
    } else {
      sim_cond <- vector("list", nrep)
      for (i in seq_len(nrep)) {
        sim_cond[[i]] <- tryCatch({
          .GeoSimcondOne(
            i, sim_data_list, krig_pred_vec, FALSE,
            weights_obj, weights_orientation, n_obs, n_loc
          )
        }, error = function(e) {
          warning(paste("Simulation", i, "failed:", e$message))
          NULL
        })
        p()
      }

      sim_cond <- sim_cond[!vapply(sim_cond, is.null, logical(1L))]
      failed_sims <- nrep - length(sim_cond)
      if (failed_sims > 0) {
        warning(paste(failed_sims, "simulations out of", nrep, "failed"))
      }
    }
  } else {
    stop("Gauss_cd currently not implemented for spacetime models")
  }

  sim_cond
}

#################################################################################
### Direct Binomial / Negative-Binomial count-constrained conditional sampler
#################################################################################
.GeoSimcondCountRF <- function(
    model, coord_obs, loc, data, n, nloc, param, corrmodel,
    mean_obs, mean_loc, distance = "Eucl", radius = 1,
    anisopars = NULL, nrep = 1L, n_iter = 25L, mcmc_thin = 1L,
    method = "Cholesky", local = FALSE, parallel = FALSE, progress = FALSE
) {
  if (!(model %in% c("Binomial", "BinomialNeg"))) {
    stop("Internal error: unsupported count model in .GeoSimcondCountRF().",
         call. = FALSE)
  }
  if (!identical(method, "Cholesky")) {
    stop(
      paste0(
        "Direct ", model, " conditional simulation currently uses the exact dense ",
        "latent Gaussian covariance and requires method = 'Cholesky'."
      ),
      call. = FALSE
    )
  }
  if (isTRUE(local)) {
    stop(
      paste0(
        "Direct ", model, " conditional simulation currently requires local = FALSE; ",
        "the count-constrained latent Gibbs sampler conditions on the full observed field."
      ),
      call. = FALSE
    )
  }
  if (!is.numeric(mcmc_thin) || length(mcmc_thin) != 1L ||
      !is.finite(mcmc_thin) || mcmc_thin < 1 ||
      mcmc_thin != as.integer(mcmc_thin)) {
    stop("mcmc_thin must be a positive integer.", call. = FALSE)
  }
  mcmc_thin <- as.integer(mcmc_thin)

  coord_obs <- as.matrix(coord_obs)
  loc <- as.matrix(loc)
  data <- as.numeric(data)
  mean_obs <- as.numeric(mean_obs)
  mean_loc <- as.numeric(mean_loc)
  m <- nrow(coord_obs)
  p <- nrow(loc)

  if (length(mean_obs) != m || length(mean_loc) != p ||
      any(!is.finite(c(mean_obs, mean_loc)))) {
    stop("Invalid latent probit means for count conditional simulation.", call. = FALSE)
  }
  tol <- sqrt(.Machine$double.eps)
  if (any(data < 0) || any(abs(data - round(data)) > tol)) {
    stop(sprintf("%s observations must be non-negative integers.", model),
         call. = FALSE)
  }
  y <- as.integer(round(data))

  sizes <- .GeoKrigCountSizes(
    model = model, n = n, nloc = nloc,
    nobs = m, npred = p, context = "GeoSimcond"
  )

  if (model == "Binomial") {
    if (any(y > sizes$obs)) {
      stop("Binomial observations cannot exceed the corresponding number of trials n.",
           call. = FALSE)
    }
    free_len <- as.integer(sizes$obs)
    success_count <- y
    forced_success <- integer(m)
    Lobs <- max(free_len)
    model_code <- 0L
  } else {
    r <- as.integer(sizes$obs[1L])
    free_len <- as.integer(y + r - 1L)
    success_count <- rep.int(r - 1L, m)
    forced_success <- as.integer(y + r)
    Lobs <- max(forced_success)
    model_code <- 1L
  }

  latent_cells <- as.double(Lobs) * as.double(m)
  if (!is.finite(latent_cells) || latent_cells > .Machine$integer.max) {
    stop("The latent count-constraint array is too large for the current implementation.",
         call. = FALSE)
  }
  if (latent_cells * 16 / 1024^2 > 1024) {
    warning(
      "The count-constrained latent Gibbs state exceeds approximately 1 GB; conditional simulation may be memory intensive.",
      call. = FALSE
    )
  }

  ## The direct count constructions threshold unit-variance Gaussian fields.
  ## Build exactly that latent covariance, retaining the fitted nugget and
  ## correlation parameters but removing marginal count parameters.
  corr_names <- CorrParam(corrmodel)
  latent_param <- as.list(param[corr_names])
  latent_param$mean <- 0
  latent_param$nugget <- if (is.null(param$nugget)) 0 else as.numeric(param$nugget)[1L]
  latent_param$sill <- 1

  all_coord <- rbind(coord_obs, loc)
  latent_cov <- GeoCovmatrix(
    coordx = all_coord, corrmodel = corrmodel, distance = distance,
    grid = FALSE, model = "Gaussian", n = 1,
    param = latent_param, anisopars = anisopars, radius = radius,
    sparse = FALSE, copula = NULL, X = NULL,
    check.duplicates = FALSE
  )$covmatrix
  Rall <- as.matrix(latent_cov)
  rm(latent_cov)

  expected_dim <- m + p
  if (!identical(dim(Rall), c(expected_dim, expected_dim)) ||
      any(!is.finite(Rall))) {
    stop("Could not construct a finite latent Gaussian covariance matrix.",
         call. = FALSE)
  }
  Rall <- 0.5 * (Rall + t(Rall))
  io <- seq_len(m)
  il <- m + seq_len(p)
  Roo <- Rall[io, io, drop = FALSE]
  Rol <- Rall[io, il, drop = FALSE]
  Rll <- Rall[il, il, drop = FALSE]
  rm(Rall)

  chol_obs <- tryCatch(chol(Roo), error = function(e) NULL)
  if (is.null(chol_obs)) {
    stop("The latent observation covariance matrix is not positive definite.",
         call. = FALSE)
  }
  Q <- chol2inv(chol_obs)
  Q <- 0.5 * (Q + t(Q))
  W <- Q %*% Rol

  cond_cov <- Rll - crossprod(Rol, W)
  cond_cov <- 0.5 * (cond_cov + t(cond_cov))
  chol_cond <- tryCatch(chol(cond_cov), error = function(e) NULL)
  if (is.null(chol_cond)) {
    stop(
      paste0(
        "The latent conditional covariance at loc is singular or not positive definite. ",
        "For direct count conditional simulation, prediction locations must not duplicate ",
        "observed locations or create an exactly singular prediction covariance."
      ),
      call. = FALSE
    )
  }
  chol_loc <- tryCatch(chol(Rll), error = function(e) NULL)
  if (is.null(chol_loc)) {
    stop("The latent covariance among prediction locations is not positive definite.",
         call. = FALSE)
  }

  if (isTRUE(parallel)) {
    warning(
      "parallel=TRUE is currently ignored by the native single-chain count-constrained Gibbs sampler.",
      call. = FALSE
    )
  }
  if (isTRUE(progress)) {
    message(
      sprintf(
        "Running count-constrained latent Gibbs sampler: %d burn-in sweeps, %d retained simulations, thinning %d.",
        as.integer(n_iter), as.integer(nrep), mcmc_thin
      )
    )
  }

  ans <- .Call(
    "GeoCountCondSim",
    Q,
    as.numeric(mean_obs),
    as.integer(free_len),
    as.integer(success_count),
    as.integer(forced_success),
    W,
    as.numeric(mean_loc),
    chol_cond,
    chol_loc,
    as.integer(sizes$pred),
    as.integer(model_code),
    as.integer(n_iter),
    as.integer(mcmc_thin),
    as.integer(nrep),
    as.integer(Lobs),
    PACKAGE = "GeoModels"
  )

  if (!is.matrix(ans) || !identical(dim(ans), c(p, as.integer(nrep))) ||
      any(!is.finite(ans))) {
    stop("The native count conditional simulator returned an invalid result.",
         call. = FALSE)
  }

  lapply(seq_len(nrep), function(j) as.numeric(ans[, j]))
}
