# GeoVarest: parametric score bootstrap for GeoModels
# Estimates J = Var{score(theta_hat)} using GeoFit(score=TRUE),
# then G^{-1} = H^{-1} J H^{-1}.

GeoVarest <- function(fit, K = 100, sparse = FALSE,
                      method = c("cholesky", "TB", "CE"),
                      alpha = 0.95, L = 10000,
                      parallel = TRUE, ncores = NULL, progress = TRUE,
                      seed = NULL) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  method <- match.arg(method)

  ## ---------------------------------------------------------------
  ## 0. Checks
  ## ---------------------------------------------------------------
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) || K < 2) {
    stop("K must be an integer >= 2", call. = FALSE)
  }
  K <- as.integer(K)

  if (!is.logical(progress) || length(progress) != 1L) {
    stop("progress must be logical", call. = FALSE)
  }

  if (!(is.numeric(alpha) && length(alpha) == 1L && alpha > 0 && alpha < 1)) {
    stop("alpha must be a single numeric in (0,1)", call. = FALSE)
  }


  if (is.null(fit$sensmat)) {
    stop("Sensitivity matrix is missing: use sensitivity = TRUE in GeoFit", call. = FALSE)
  }

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }

    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)

    set.seed(seed)
  }

  ## ---------------------------------------------------------------
  ## Progress handlers
  ## ---------------------------------------------------------------
  use_progressr <- isTRUE(progress) && requireNamespace("progressr", quietly = TRUE)

  if (isTRUE(progress) && !use_progressr) {
    warning(
      "progress=TRUE but progressr is not available; using text messages only.",
      call. = FALSE
    )
  }

  if (use_progressr) {
    old_handlers <- progressr::handlers()

    on.exit({
      if (length(old_handlers) > 0L) {
        progressr::handlers(old_handlers)
      } else {
        progressr::handlers("default")
      }
    }, add = TRUE)

    progressr::handlers(global = TRUE)
    progressr::handlers("txtprogressbar")
  }

  gm_fun <- function(name) {
    if (requireNamespace("GeoModels", quietly = TRUE) &&
        exists(name, envir = asNamespace("GeoModels"), inherits = FALSE)) {
      return(getFromNamespace(name, "GeoModels"))
    }
    if (exists(name, mode = "function", inherits = TRUE)) {
      return(get(name, mode = "function", inherits = TRUE))
    }
    stop("Function ", name, " is not available. Load GeoModels first.", call. = FALSE)
  }


  set_proj_env_for_workers <- function() {
    get_proj_path <- function() {
      path <- Sys.getenv("PROJ_DATA", unset = "")

      if (!nzchar(path)) {
        path <- Sys.getenv("PROJ_LIB", unset = "")
      }

     # if (!nzchar(path) && requireNamespace("sf", quietly = TRUE)) {
     #   path <- tryCatch(sf::sf_proj_info("path"), error = function(e) "")
     # }

      if (length(path) != 1L || !nzchar(path)) {
        return("")
      }

      normalizePath(path, winslash = "/", mustWork = FALSE)
    }

    proj_path <- get_proj_path()

    if (!nzchar(proj_path)) {
      return(character(0))
    }

    Sys.setenv(PROJ_DATA = proj_path)
    Sys.setenv(PROJ_LIB  = proj_path)

    c(PROJ_DATA = proj_path, PROJ_LIB = proj_path)
  }

  GeoFit_fun <- gm_fun("GeoFit")
  GeoSim_fun <- gm_fun("GeoSim")

  if (length(fit$coordt) == 1L) fit$coordt <- NULL

  fit$thin_method <- fit$thin_method %||% "bernoulli"
  fit$p_neighb    <- fit$p_neighb %||% 1

  optimizer <- fit$optimizer %||% "Nelder-Mead"
  lower     <- fit$lower
  upper     <- fit$upper

  theta_hat <- as.numeric(unlist(fit$param, use.names = TRUE))
  names(theta_hat) <- names(unlist(fit$param, use.names = TRUE))

  num_params <- length(theta_hat)
  if (num_params < 1L) {
    stop("fit$param is empty", call. = FALSE)
  }

  H <- as.matrix(fit$sensmat)
  storage.mode(H) <- "double"

  if (!all(dim(H) == c(num_params, num_params))) {
    stop("fit$sensmat has incompatible dimension", call. = FALSE)
  }

  H <- (H + t(H)) / 2
  rownames(H) <- colnames(H) <- names(theta_hat)

  safe_inverse <- function(A, label) {
    A <- as.matrix(A)
    storage.mode(A) <- "double"
    A <- (A + t(A)) / 2

    out <- try(solve(A), silent = TRUE)

    if (!inherits(out, "try-error") && all(is.finite(out))) {
      return((out + t(out)) / 2)
    }

    ee <- eigen(A, symmetric = TRUE)
    tol <- max(dim(A)) * max(abs(ee$values)) * .Machine$double.eps
    keep <- abs(ee$values) > tol

    if (!any(keep)) {
      stop("Cannot invert ", label, call. = FALSE)
    }

    out <- ee$vectors[, keep, drop = FALSE] %*%
      diag(1 / ee$values[keep], nrow = sum(keep)) %*%
      t(ee$vectors[, keep, drop = FALSE])

    warning("Used spectral generalized inverse for ", label, call. = FALSE)
    (out + t(out)) / 2
  }

  ## Dimension used in CLBIC penalty.
  if (!is.null(fit$coordx_dyn)) {
    dimat <- sum(fit$ns)
  } else {
    dimat <- fit$numtime * fit$numcoord
  }

  ## Same X handling as GeoVarestbootstrap.
  if (!is.null(fit$X) && !is.null(dim(fit$X)) && ncol(fit$X) == 1L) {
    ncheck <- min(NROW(fit$X), dimat)
    X_use <- if (ncheck > 0 && all(fit$X[seq_len(ncheck), 1] == 1)) NULL else fit$X
  } else {
    X_use <- fit$X
  }

  coords <- if (!is.null(fit$coordz)) {
    cbind(fit$coordx, fit$coordy, fit$coordz)
  } else {
    cbind(fit$coordx, fit$coordy)
  }

  if (isTRUE(fit$bivariate) && is.null(fit$coordx_dyn)) {
    if (nrow(coords) %% 2L != 0L) {
      stop("bivariate=TRUE but odd number of coordinates", call. = FALSE)
    }
    coords <- coords[seq_len(nrow(coords) / 2L), , drop = FALSE]
  }

  ## Misspecification mapping.
  model_sim <- fit$model
  model_est <- fit$model

  if (isTRUE(fit$missp)) {
    model_map <- c(
      StudentT     = "Gaussian_misp_StudentT",
      Poisson      = "Gaussian_misp_Poisson",
      PoissonZIP   = "Gaussian_misp_PoissonZIP",
      SkewStudentT = "Gaussian_misp_SkewStudentT",
      Tukeygh      = "Gaussian_misp_Tukeygh"
    )

    if (as.character(fit$model) %in% names(model_map)) {
      model_est <- model_map[[as.character(fit$model)]]
    }
  }

  estimated_size_mb <- (K * fit$numtime * fit$numcoord * 8) / (1024^2)

  if (estimated_size_mb > 500) {
    warning(
      sprintf("Estimated simulated dataset size: %.1f MB", estimated_size_mb),
      call. = FALSE
    )
  }

  ## ---------------------------------------------------------------
  ## 1. Simulate K datasets
  ## ---------------------------------------------------------------
  sim_args <- list(
    coordx     = coords,
    coordt     = fit$coordt,
    coordx_dyn = fit$coordx_dyn,
    anisopars  = fit$anisopars,
    corrmodel  = fit$corrmodel,
    model      = model_sim,
    param      = append(fit$param, fit$fixed),
    grid       = fit$grid,
    X          = fit$X,
    n          = fit$n,
    distance   = fit$distance,
    radius     = fit$radius,
    nrep       = K,
    progress   = progress
  )

  if (is.null(fit$copula)) {
    if (method == "cholesky") {
      data_sim_full <- do.call(
        GeoSim_fun,
        c(sim_args, list(sparse = sparse, method = method))
      )
    } else {
      GeoSimapprox_fun <- gm_fun("GeoSimapprox")

      ## GeoSimapprox may parallelize over nrep = K. This simulation phase
      ## is completed before score evaluation starts, so the score
      ## parallelism below is not nested. We defensively restore the current
      ## future plan before entering the score phase.
      old_plan_before_sim <- NULL
      if (isTRUE(parallel) && requireNamespace("future", quietly = TRUE)) {
        old_plan_before_sim <- future::plan()
      }

      data_sim_full <- do.call(
        GeoSimapprox_fun,
        c(
          sim_args,
          list(
            method   = method,
            L        = L,
            parallel = parallel,
            ncores   = ncores
          )
        )
      )

      if (!is.null(old_plan_before_sim)) {
        try(future::plan(old_plan_before_sim), silent = TRUE)
      }
    }
  } else {
    if (method != "cholesky") {
      stop("Unsupported method for copula simulation", call. = FALSE)
    }

    GeoSimCopula_fun <- gm_fun("GeoSimCopula")

    data_sim_full <- do.call(
      GeoSimCopula_fun,
      c(sim_args, list(copula = fit$copula, sparse = sparse, method = method))
    )
  }

  data_sim <- data_sim_full$data
  rm(data_sim_full)
  gc(verbose = FALSE, full = TRUE)

  ## ---------------------------------------------------------------
  ## 2. Evaluate composite likelihood and score through GeoFit
  ## ---------------------------------------------------------------
 is_stochastic_thinning <- function() {
  tm <- tolower(as.character(fit$thin_method %||% ""))
  tm %in% c(
    "bernoulli",
    "fixedbudget",
    "targetbalanced",
    "match"
  ) && isTRUE(fit$p_neighb < 1)
}

  with_seed_safe <- function(seed_value, expr) {
    if (is.null(seed_value) || !is.finite(seed_value)) {
      return(force(expr))
    }

    if (requireNamespace("withr", quietly = TRUE)) {
      return(withr::with_seed(as.integer(seed_value), force(expr)))
    }

    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }

    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)

    set.seed(as.integer(seed_value))
    force(expr)
  }

  make_start <- function(par) {
    par <- as.numeric(par)
    names(par) <- names(theta_hat)
    as.list(par)
  }

  geofit_onlyvar <- function(current_data, par, want_score, seed_value) {
    call_it <- function() {
      GeoFit_fun(
        data        = current_data,
        start       = make_start(par),
        fixed       = fit$fixed,
        coordx      = coords,
        coordt      = fit$coordt,
        coordx_dyn  = fit$coordx_dyn,
        copula      = fit$copula,
        anisopars   = fit$anisopars,
        est.aniso   = fit$est.aniso,
        thin_method = fit$thin_method,
        lower       = lower,
        upper       = upper,
        neighb      = fit$neighb,
        p_neighb    = fit$p_neighb,
        corrmodel   = fit$corrmodel,
        model       = model_est,
        sparse      = FALSE,
        n           = fit$n,
        maxdist     = fit$maxdist,
        maxtime     = fit$maxtime,
        memdist     = fit$memdist %||% TRUE,
        optimizer   = optimizer,
        grid        = fit$grid,
        likelihood  = fit$likelihood,
        type        = fit$type,
        X           = X_use,
        distance    = fit$distance,
        radius      = fit$radius,
        onlyvar     = TRUE,
        score       = want_score,
        sensitivity = FALSE,
        varest      = FALSE,
        weighted    = fit$weighted %||% FALSE
      )
    }

    if (is_stochastic_thinning()) {
      with_seed_safe(seed_value, call_it())
    } else {
      call_it()
    }
  }

  geofit_score <- function(current_data, par, seed_value) {
    capture.output({
      ff <- geofit_onlyvar(
        current_data = current_data,
        par          = par,
        want_score   = TRUE,
        seed_value   = seed_value
      )
    }, file = nullfile())

    sc <- as.numeric(unlist(ff$score, use.names = TRUE))

    if (length(sc) != num_params || any(!is.finite(sc))) {
      stop("GeoFit returned an invalid score", call. = FALSE)
    }

    ## In CompLik2 the numerical score is the gradient of the minimized
    ## objective, i.e. -CL/npairs when the normalized CompLik2 is used.
    ## Change sign so J is based on the composite-likelihood score.
    sc <- -sc
    names(sc) <- names(theta_hat)

    val <- if (!is.null(ff$logCompLik)) ff$logCompLik else ff$logLik
    val <- as.numeric(val)

    if (length(val) != 1L || !is.finite(val)) {
      stop("non-finite log composite likelihood", call. = FALSE)
    }

    list(score = sc, logCompLik = val)
  }

  make_score <- function(k, current_data) {
    tryCatch({
      seed_thin <- as.integer(1234567L + k * 9999L)

      out <- geofit_score(current_data, theta_hat, seed_thin)

      c(as.numeric(out$score), logCompLik = as.numeric(out$logCompLik))

    }, error = function(e) {
      z <- rep(NA_real_, num_params + 1L)
      attr(z, "error") <- conditionMessage(e)
      z
    })
  }

  ## ---------------------------------------------------------------
  ## 3. Score loop. Parallelization is only here.
  ## ---------------------------------------------------------------
  if (progress) {
    cat("Computing", K, "scores ...\n")
  }

  use_parallel <- isTRUE(parallel) && K > 1L &&
    requireNamespace("future", quietly = TRUE) &&
    requireNamespace("future.apply", quietly = TRUE)

  if (isTRUE(parallel) && !use_parallel) {
    warning(
      "Parallel evaluation requested but future/future.apply is unavailable; using sequential evaluation.",
      call. = FALSE
    )
  }

  if (use_parallel) {
    coremax <- parallel::detectCores()

    if (is.na(coremax) || coremax <= 1L) {
      use_parallel <- FALSE
      ncores <- 1L
    } else {
      ncores <- max(
        1L,
        min(if (is.null(ncores)) getOption("mc.cores", 4L) else as.integer(ncores), K)
      )
    }
  }

  if (use_parallel) {
    old_plan <- future::plan()
    on.exit(try(future::plan(old_plan), silent = TRUE), add = TRUE)

    ## Il payload esportato ai worker (fit, coords, neighb, ecc. catturati
    ## dalle closure geofit_onlyvar/geofit_score/make_score) e'
    ## necessario per rifare la verosimiglianza composita su ogni dataset
    ## simulato: con dataset grandi puo' legittimamente superare il default
    ## di 500 MiB di future.globals.maxSize. Alziamo il limite in modo
    ## dinamico, proporzionato alla dimensione reale, e lo ripristiniamo a
    ## fine chiamata per non alterare in modo permanente le opzioni globali
    ## dell'utente (rilevante soprattutto trattandosi di una funzione di
    ## pacchetto).
    needed_size <- tryCatch({
      max(vapply(
        list(make_score, geofit_onlyvar, geofit_score),
        function(f) as.numeric(utils::object.size(f)),
        numeric(1)
      ))
    }, error = function(e) NA_real_)

    old_max_size <- getOption("future.globals.maxSize")
    on.exit(options(future.globals.maxSize = old_max_size), add = TRUE)

    safety_margin <- 1.5
    default_cap   <- 500 * 1024^2

    options(future.globals.maxSize = if (is.na(needed_size)) {
      max(default_cap, 4 * 1024^3)
    } else {
      max(default_cap, needed_size * safety_margin)
    })

    old_proj_data <- Sys.getenv("PROJ_DATA", unset = NA)
    old_proj_lib  <- Sys.getenv("PROJ_LIB",  unset = NA)

    on.exit({
      if (is.na(old_proj_data)) {
        Sys.unsetenv("PROJ_DATA")
      } else {
        Sys.setenv(PROJ_DATA = old_proj_data)
      }
      if (is.na(old_proj_lib)) {
        Sys.unsetenv("PROJ_LIB")
      } else {
        Sys.setenv(PROJ_LIB = old_proj_lib)
      }
    }, add = TRUE)

    rscript_envs <- set_proj_env_for_workers()

    future::plan(future::multisession, workers = ncores)

    temp_dir <- tempdir()
    tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

    data_files <- vapply(seq_len(K), function(k) {
      fp <- file.path(temp_dir, sprintf("scoreJ_%s_%04d.rds", tag, k))
      saveRDS(data_sim[[k]], fp, compress = FALSE)
      fp
    }, character(1))

    on.exit(unlink(data_files), add = TRUE)

    rm(data_sim)
    gc(verbose = FALSE, full = TRUE)

    if (use_progressr) {
      score_list <- progressr::with_progress({
        pb <- progressr::progressor(along = seq_len(K))

        future.apply::future_lapply(seq_len(K), function(k) {
          if (requireNamespace("GeoModels", quietly = TRUE)) invisible(NULL)
          out <- make_score(k, readRDS(data_files[[k]]))
          pb(sprintf("score %d", k))
          out
        }, future.seed = TRUE)
      })
    } else {
      score_list <- future.apply::future_lapply(seq_len(K), function(k) {
        if (requireNamespace("GeoModels", quietly = TRUE)) invisible(NULL)
        make_score(k, readRDS(data_files[[k]]))
      }, future.seed = TRUE)
    }

  } else {
    if (use_progressr) {
      score_list <- progressr::with_progress({
        pb <- progressr::progressor(along = seq_len(K))
        out <- vector("list", K)

        for (k in seq_len(K)) {
          out[[k]] <- make_score(k, data_sim[[k]])
          pb(sprintf("score %d", k))
        }

        out
      })
    } else {
      score_list <- vector("list", K)

      for (k in seq_len(K)) {
        if (progress && (k == 1L || k %% 10L == 0L || k == K)) {
          cat("score", k, "of", K, "\n")
        }

        score_list[[k]] <- make_score(k, data_sim[[k]])
      }
    }

    rm(data_sim)
    gc(verbose = FALSE, full = TRUE)
  }

  errors <- vapply(
    score_list,
    function(x) attr(x, "error") %||% "",
    character(1)
  )

  score_raw <- do.call(rbind, score_list)
  colnames(score_raw) <- c(names(theta_hat), "logCompLik")

  valid <- apply(
    score_raw[, seq_len(num_params), drop = FALSE],
    1L,
    function(z) all(is.finite(z))
  ) & is.finite(score_raw[, num_params + 1L])

  score_out <- score_raw[valid, , drop = FALSE]
  n_successful <- nrow(score_out)

  if (n_successful < 2L) {
    shown <- utils::head(errors[nzchar(errors)], 8L)

    msg <- paste0(
      "Insufficient successful score evaluations: ",
      n_successful,
      " out of ",
      K
    )

    if (length(shown) > 0L) {
      msg <- paste0(
        msg,
        "\nFirst errors:\n- ",
        paste(shown, collapse = "\n- ")
      )
    }

    stop(msg, call. = FALSE)
  }

  if (progress) {
    cat("Successful score evaluations:", n_successful, "out of", K, "\n")
  }

  ## ---------------------------------------------------------------
  ## 4. Sandwich/Godambe and criteria
  ## ---------------------------------------------------------------
  score_mat <- score_out[, seq_len(num_params), drop = FALSE]

  J <- stats::var(score_mat)
  J <- (J + t(J)) / 2
  rownames(J) <- colnames(J) <- names(theta_hat)

  Hinv <- safe_inverse(H, "H/sensmat")

  Ginv <- Hinv %*% J %*% Hinv
  Ginv <- (Ginv + t(Ginv)) / 2
  rownames(Ginv) <- colnames(Ginv) <- names(theta_hat)

  G <- safe_inverse(Ginv, "Ginv/varcov")
  rownames(G) <- colnames(G) <- names(theta_hat)

  stderr <- sqrt(pmax(0, diag(Ginv)))
  names(stderr) <- names(theta_hat)

  penalty <- sum(diag(Hinv %*% J))
  penalty_alt <- sum(diag(H %*% Ginv))

  fit_loglik <- if (!is.null(fit$logCompLik)) fit$logCompLik else fit$logLik
  fit_loglik <- as.numeric(fit_loglik)

  if (length(fit_loglik) != 1L || !is.finite(fit_loglik)) {
    stop("fit does not contain a finite logCompLik/logLik value", call. = FALSE)
  }

  lik <- as.character(fit$likelihood)
  typ <- as.character(fit$type)

  if ((lik == "Marginal" && typ %in% c("Independence", "Pairwise")) ||
      (lik == "Conditional" && typ == "Pairwise")) {

    claic <- -2 * fit_loglik + 2 * penalty
    clbic <- -2 * fit_loglik + log(dimat) * penalty
    fit$varimat <- H %*% Ginv %*% H

  } else if (lik == "Full" && typ == "Standard") {

    claic <- -2 * fit_loglik + 2 * num_params
    clbic <- -2 * fit_loglik + log(dimat) * num_params

  } else {
    claic <- clbic <- NA_real_
  }

  z_alpha <- stats::qnorm(1 - (1 - alpha) / 2)

  fit$stderr <- stderr
  fit$varcov <- Ginv
  fit$godambe <- G
  fit$Jmat <- J
  fit$Hinv <- Hinv

  fit$claic <- claic
  fit$clic <- claic
  fit$clbic <- clbic
  fit$clic_penalty <- penalty
  fit$clic_penalty_alt <- penalty_alt

  fit$conf.int <- rbind(
    theta_hat - z_alpha * stderr,
    theta_hat + z_alpha * stderr
  )
  colnames(fit$conf.int) <- names(theta_hat)
  rownames(fit$conf.int) <- c("Lower", "Upper")

  fit$pvalues <- 2 * stats::pnorm(-abs(theta_hat / stderr))
  names(fit$pvalues) <- names(theta_hat)

  fit$scores <- score_mat
  fit$score_logCompLik <- score_out[, "logCompLik"]

  fit$score_failures <- data.frame(
    iteration = which(!valid),
    error = errors[!valid],
    stringsAsFactors = FALSE
  )

  fit$estimates <- score_out
  attr(fit$estimates, "content") <- "scores at theta_hat, not refitted parameter estimates"

  fit$conf.int.bootstrap_percentile <- NULL
  fit$pvalues_bootstrap <- NULL

  fit$bootstrap_type <- "parametric_score_J_geofit_score"
  fit$bootstrap_K <- K
  fit$bootstrap_successful <- n_successful
  fit$bootstrap_success_rate <- n_successful / K
  fit$bootstrap_thin_method <- fit$thin_method
  fit$bootstrap_p_neighb <- fit$p_neighb
  fit$bootstrap_thinning_seed_rule <- "1234567 + k * 9999 for stochastic thinning"

  fit
}
