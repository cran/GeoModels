## Internal utilities shared by GeoVarest and GeoVarestbootstrap.
## Keep bootstrap plumbing in one place while the two public functions retain
## their distinct statistical post-processing.

.GeoVarest_preserve_seed <- function(seed) {
  if (is.null(seed)) return(function() invisible(NULL))

  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }

  set.seed(seed)

  function() {
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    invisible(NULL)
  }
}

.GeoVarest_prepare_geometry <- function(fit, coordt_use = fit$coordt) {
  dimat <- if (!is.null(fit$coordx_dyn)) sum(fit$ns) else fit$numtime * fit$numcoord

  if (!is.null(fit$X) && !is.null(dim(fit$X)) && ncol(fit$X) == 1L) {
    ncheck <- min(NROW(fit$X), dimat)
    X_use <- if (ncheck > 0L && all(fit$X[seq_len(ncheck), 1L] == 1)) NULL else fit$X
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

  list(dimat = dimat, X = X_use, coords = coords, coordt = coordt_use)
}

.GeoVarest_make_fit_spec <- function(fit, boot_par, model_est, X_use,
                                     optimizer, lower, upper, pair_index,
                                     coordt_use, memdist = TRUE,
                                     weighted = FALSE, include_start = FALSE) {
  list(
    param       = if (isTRUE(include_start)) boot_par$start else NULL,
    fixed       = boot_par$fixed,
    coordt      = coordt_use,
    coordx_dyn  = fit$coordx_dyn,
    copula      = fit$copula,
    anisopars   = boot_par$anisopars,
    est.aniso   = fit$est.aniso,
    aniso_names = if (is.null(boot_par$anisopars)) character(0) else
      intersect(c("angle", "ratio"), names(unlist(boot_par$anisopars, use.names = TRUE))),
    thin_method = fit$thin_method,
    lower       = lower,
    upper       = upper,
    neighb      = fit$neighb,
    p_neighb    = fit$p_neighb,
    corrmodel   = fit$corrmodel,
    model       = model_est,
    n           = fit$n,
    maxdist     = fit$maxdist,
    maxtime     = fit$maxtime,
    memdist     = memdist,
    optimizer   = optimizer,
    grid        = fit$grid,
    likelihood  = fit$likelihood,
    type        = fit$type,
    X           = X_use,
    distance    = fit$distance,
    radius      = fit$radius,
    weighted    = weighted,
    bivariate   = fit$bivariate,
    spacetime   = fit$spacetime,
    pair_index  = pair_index
  )
}

.GeoVarest_start_from_theta <- function(theta_hat, aniso_names = character(0)) {
  out <- as.list(as.numeric(theta_hat))
  names(out) <- names(theta_hat)
  if (length(aniso_names)) out[intersect(names(out), aniso_names)] <- NULL
  out
}

.GeoVarest_fit_args <- function(data, coords, fit_spec, start,
                                score = FALSE) {
  args <- list(
    data        = data,
    start       = start,
    fixed       = fit_spec$fixed,
    coordx      = coords,
    coordt      = fit_spec$coordt,
    coordx_dyn  = fit_spec$coordx_dyn,
    copula      = fit_spec$copula,
    anisopars   = fit_spec$anisopars,
    est.aniso   = fit_spec$est.aniso,
    thin_method = fit_spec$thin_method,
    lower       = fit_spec$lower,
    upper       = fit_spec$upper,
    neighb      = fit_spec$neighb,
    p_neighb    = fit_spec$p_neighb,
    corrmodel   = fit_spec$corrmodel,
    model       = fit_spec$model,
    sparse      = FALSE,
    n           = fit_spec$n,
    maxdist     = fit_spec$maxdist,
    maxtime     = fit_spec$maxtime,
    memdist     = fit_spec$memdist,
    optimizer   = fit_spec$optimizer,
    grid        = fit_spec$grid,
    likelihood  = fit_spec$likelihood,
    type        = fit_spec$type,
    X           = fit_spec$X,
    distance    = fit_spec$distance,
    radius      = fit_spec$radius,
    weighted    = fit_spec$weighted
  )

  if (isTRUE(score)) {
    args$onlyvar <- TRUE
    args$score <- TRUE
    args$sensitivity <- FALSE
    args$varest <- FALSE
  }
  args
}

.GeoVarest_worker_context <- function(fit_spec, coords) {
  list(
    pair_info = .GeoVarest_prepare_pair_info(fit_spec$pair_index, coords, fit_spec),
    fit_context = if (as.character(fit_spec$type)[1L] %in% c("Pairwise", "Independence")) {
      new.env(parent = emptyenv())
    } else {
      NULL
    }
  )
}

.GeoVarest_run_fit <- function(data, coords, fit_spec, start,
                               pair_info, fit_context, score = FALSE) {
  GeoFit_fun <- getFromNamespace("GeoFit", "GeoModels")
  fit_args <- .GeoVarest_fit_args(data, coords, fit_spec, start, score = score)
  out <- NULL
  capture.output({
    out <- .GeoVarest_with_selected_pairs(
      fit_spec$pair_index,
      do.call(GeoFit_fun, fit_args),
      pair_info = pair_info,
      fit_context = fit_context
    )
  }, file = nullfile())
  out
}

.GeoVarest_extract_score <- function(fit_result, target_names) {
  num_params <- length(target_names)
  sc <- unlist(fit_result$score, use.names = TRUE)
  if (length(sc) != num_params || any(!is.finite(sc))) {
    stop("GeoFit returned an invalid score", call. = FALSE)
  }
  if (!is.null(names(sc)) && all(target_names %in% names(sc))) sc <- sc[target_names]

  ll <- if (!is.null(fit_result$logCompLik)) fit_result$logCompLik else fit_result$logLik
  ll <- as.numeric(ll)
  if (length(ll) != 1L || !is.finite(ll)) {
    stop("non-finite log composite likelihood", call. = FALSE)
  }

  ans <- c(-as.numeric(sc), logCompLik = ll)
  names(ans) <- c(target_names, "logCompLik")
  ans
}

.GeoVarest_extract_refit <- function(fit_result, target_names) {
  bad <- function() {
    z <- rep(NA_real_, length(target_names) + 1L)
    names(z) <- c(target_names, "logCompLik")
    z
  }

  if (is.null(fit_result) || !identical(fit_result$convergence, "Successful")) return(bad())
  ll <- if (!is.null(fit_result$logCompLik)) fit_result$logCompLik else fit_result$logLik
  if (length(ll) != 1L || !is.finite(ll)) return(bad())

  th <- unlist(fit_result$param, use.names = TRUE)
  if (!is.null(names(th)) && all(target_names %in% names(th))) th <- th[target_names]
  if (length(th) != length(target_names) || any(!is.finite(th))) return(bad())

  ans <- c(as.numeric(th), logCompLik = as.numeric(ll))
  names(ans) <- c(target_names, "logCompLik")
  ans
}

.GeoVarest_failed_result <- function(target_names, error = NULL) {
  z <- rep(NA_real_, length(target_names) + 1L)
  names(z) <- c(target_names, "logCompLik")
  if (!is.null(error)) attr(z, "error") <- conditionMessage(error)
  z
}
