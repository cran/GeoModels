GeoNeighbSelect <- function(
  data, coordx, coordy=NULL, coordz=NULL, coordt=NULL, coordx_dyn=NULL,
  copula=NULL, corrmodel=NULL, distance="Eucl", fixed=NULL, anisopars=NULL,
  est.aniso=c(FALSE,FALSE), grid=FALSE, likelihood='Marginal',
  lower=NULL, neighb=c(1,2,3,4,5), p_neighb=1,
  maxtime=Inf, memdist=TRUE, model='Gaussian', n=1, ncores=NULL,
  optimizer='Nelder-Mead', parallel=FALSE, bivariate=FALSE, radius=1, start=NULL,
  type='Pairwise', upper=NULL, weighted=FALSE,
  X=NULL, nosym=FALSE, spobj=NULL, spdata=NULL, vario=NULL, progress=TRUE
)
{
  # --- Input validation (FIX: no coercion warning) ---
  if (!is.numeric(neighb) || length(neighb) == 0L)
    stop("neighb must be a non-empty numeric vector")

  if (anyNA(neighb) || any(!is.finite(neighb)))
    stop("neighb must contain only finite, non-missing values")

  if (any(neighb <= 0))
    stop("neighb must contain only positive values")

  tol <- sqrt(.Machine$double.eps)
  if (any(abs(neighb - round(neighb)) > tol))
    stop("neighb must be a positive integer numeric vector")

  # normalize to integers (recommended)
  neighb <- as.integer(round(neighb))

  if (!is.numeric(p_neighb) || length(p_neighb) != 1L || !is.finite(p_neighb))
    stop("p_neighb must be a finite numeric scalar.")
  if (p_neighb <= 0 || p_neighb > 1)
    stop("p_neighb must be in (0,1].")

  # --- Plateau tolerance ---
  tau <- 0.15

  # Validate correlation model
  corrmodel_valid <- try(CkCorrModel(corrmodel), silent = TRUE)
  if(inherits(corrmodel_valid, "try-error")) stop("correlation model is not valid\n")

  # NOTE: we overwrite 'bivariate' to match the model type (as in your original logic)
  bivariate <- CheckBiv(corrmodel_valid)
  spacetime <- CheckST(corrmodel_valid)
  space <- !spacetime && !bivariate

  # Validate variogram object
  if(!is.null(vario)) {
    if(!inherits(vario, "GeoVariogram")) stop("A GeoVariogram object is needed as input for vario\n")
    ok_type <- (vario$bivariate && bivariate) ||
      (!is.null(vario$bint) && spacetime) ||
      (is.null(vario$bint) && !bivariate)
    if(!ok_type) stop("The GeoVariogram object is not of the same type of the correlation model\n")
    semiv <- vario
  } else stop("A GeoVariogram object is needed as input for vario\n")

  # Validate start parameter
  if(is.null(start) || length(start) == 0) stop("start parameter must be provided\n")

  # Pre-compute common values
  coremax <- parallel::detectCores()
  if(is.na(coremax) || coremax <= 1) parallel <- FALSE
  K <- length(neighb)
  P <- NULL

  # Pre-allocation
  if(spacetime) {
    if(!is.numeric(maxtime)) stop("maxtime must be a numeric vector")
    P <- length(maxtime)
    estimates <- matrix(NA_real_, nrow = K * P, ncol = length(start))
    res <- rep(NA_real_, K * P)
  } else {
    estimates <- matrix(NA_real_, nrow = K, ncol = length(start))
    res <- rep(NA_real_, K)
  }

  # Common parameters for GeoFit
  # NOTE: p_neighb is passed here -> thinning happens inside
  common_params <- list(
    data = data, coordx = coordx, coordy = coordy, coordz = coordz, coordt = coordt,
    coordx_dyn = coordx_dyn, copula = copula, corrmodel = corrmodel, distance = distance,
    fixed = fixed, anisopars = anisopars, est.aniso = est.aniso, grid = grid,
    likelihood = likelihood, lower = lower, memdist = memdist, model = model, n = n,
    optimizer = optimizer, radius = radius, start = start, type = type, upper = upper,
    weighted = weighted, X = X, nosym = nosym, spobj = spobj, spdata = spdata,
    p_neighb = p_neighb
  )

  # Core fitting function
  compute_fit <- function(neighb_val, maxtime_val = Inf) {
    params <- common_params
    params$neighb  <- neighb_val
    params$maxtime <- maxtime_val

    aa <- try(suppressWarnings(do.call(GeoFit, params)), silent = TRUE)
    if (inherits(aa, "try-error") || is.null(aa)) {
      return(list(estimates = rep(NA_real_, length(start)), res = Inf))
    }

    estimates_val <- unlist(aa$param)

    if(!(aa$logCompLik == 0 || aa$logCompLik > 1e+13)) {
      # IMPORTANT: keep your "variogram vs corr" choice unchanged
      cc <- suppressWarnings(
        GeoCorrFct(
          semiv$centers, t = semiv$centert, corrmodel = corrmodel,
          model = model, distance = distance,
          param = c(aa$param, aa$fixed), radius = radius, n = n,
          covariance = TRUE, variogram = TRUE
        )$corr
      )
      res_val <- sum((semiv$variograms - cc)^2)
    } else {
      res_val <- Inf
    }

    list(estimates = estimates_val, res = res_val)
  }

  #################### SPATIAL / BIVARIATE ####################
  if(space || bivariate) {
    if(!parallel) {
      if(progress) {
        progressr::handlers(global = TRUE)
        progressr::handlers("txtprogressbar")
        pb <- progressr::progressor(along = 1:K)
      }
      for(M in seq_len(K)) {
        if(progress) pb(sprintf("k=%g", M))
        out <- compute_fit(neighb[M])
        estimates[M, ] <- out$estimates
        res[M] <- out$res
      }
    } else {
      n.cores <- if(is.null(ncores)) {
        min(coremax - 1L, K)
      } else {
        if(!is.numeric(ncores) || ncores > coremax || ncores < 1) stop("number of cores not valid\n")
        min(ncores, K)
      }

      cat("Performing", K, "estimations using", n.cores, "cores...\n")

      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(future::multisession, workers = n.cores)

      n_params <- length(start)

      if(progress) {
        progressr::handlers(global = TRUE)
        progressr::handlers("txtprogressbar")

        results <- progressr::with_progress({
          pb <- progressr::progressor(along = 1:K)

          single_fit <- function(M) {
            out <- compute_fit(neighb[M])
            pb(sprintf("k=%g", M))
            c(out$estimates, res = out$res)
          }

          future.apply::future_lapply(
            seq_len(K), single_fit,
            future.seed = TRUE,
            future.stdout = FALSE,
            future.conditions = "none"
          )
        })
      } else {
        single_fit <- function(M) {
          out <- compute_fit(neighb[M])
          c(out$estimates, res = out$res)
        }

        results <- future.apply::future_lapply(
          seq_len(K), single_fit,
          future.seed = TRUE,
          future.stdout = FALSE,
          future.conditions = "none"
        )
      }

      estimates <- do.call(rbind, lapply(results, function(x) x[seq_len(n_params)]))
      res <- vapply(results, function(x) x[n_params + 1L], numeric(1))
    }
  }

  #################### SPATIO-TEMPORAL ####################
  if(spacetime) {
    total_iter <- K * P

    # FIX: canonical grid consistent with index -> (neighb, maxtime)
    grid_result <- expand.grid(neighb = neighb, maxtime = maxtime)

    if(!parallel) {
      if(progress) {
        progressr::handlers(global = TRUE)
        progressr::handlers("txtprogressbar")
        pb <- progressr::progressor(along = 1:total_iter)
      }

      idx <- 1L
      for(L in seq_len(P)) {
        for(M in seq_len(K)) {
          if(progress) pb(sprintf("k=%g", idx))
          out <- compute_fit(neighb[M], maxtime[L])
          estimates[idx, ] <- out$estimates
          res[idx] <- out$res
          idx <- idx + 1L
        }
      }
    } else {
      n.cores <- if(is.null(ncores)) {
        min(coremax - 1L, total_iter)
      } else {
        if(!is.numeric(ncores) || ncores > coremax || ncores < 1) stop("number of cores not valid\n")
        min(ncores, total_iter)
      }

      cat("Performing", total_iter, "estimations using", n.cores, "cores...\n")

      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(future::multisession, workers = n.cores)

      n_params <- length(start)

      if(progress) {
        progressr::handlers(global = TRUE)
        progressr::handlers("txtprogressbar")

        results <- progressr::with_progress({
          pb <- progressr::progressor(along = 1:total_iter)

          single_fit_st <- function(i) {
            out <- compute_fit(grid_result$neighb[i], grid_result$maxtime[i])
            pb(sprintf("k=%g", i))
            c(out$estimates, res = out$res)
          }

          future.apply::future_lapply(
            seq_len(total_iter), single_fit_st,
            future.seed = TRUE,
            future.stdout = FALSE,
            future.conditions = "none"
          )
        })
      } else {
        single_fit_st <- function(i) {
          out <- compute_fit(grid_result$neighb[i], grid_result$maxtime[i])
          c(out$estimates, res = out$res)
        }

        results <- future.apply::future_lapply(
          seq_len(total_iter), single_fit_st,
          future.seed = TRUE,
          future.stdout = FALSE,
          future.conditions = "none"
        )
      }

      estimates <- do.call(rbind, lapply(results, function(x) x[seq_len(n_params)]))
      res <- vapply(results, function(x) x[n_params + 1L], numeric(1))
    }
  }

  # ---- BEST: absolute minimum ----
  indexmin <- which.min(res)

  if (space || bivariate) {
    bestK  <- neighb[indexmin]
    best_T <- NULL
  } else if (spacetime) {
    # FIX: map using the same grid used for computation
    grid_result <- expand.grid(neighb = neighb, maxtime = maxtime)
    bestK  <- grid_result$neighb[indexmin]
    best_T <- grid_result$maxtime[indexmin]
  } else {
    bestK  <- neighb[indexmin]
    best_T <- NULL
  }

  best_estimates <- estimates[indexmin, ]

  # ---- SUGGESTED: plateau rule (tau) ----
  finite_ok <- is.finite(res)
  if (!any(finite_ok)) {
    sugg_neighb <- bestK
    sugg_time <- best_T
  } else {
    minres <- min(res[finite_ok], na.rm = TRUE)
    idx_ok <- which(res <= (1 + tau) * minres)

    if (space || bivariate) {
      sugg_neighb <- min(neighb[idx_ok], na.rm = TRUE)
      sugg_neighb <- min(sugg_neighb, bestK)
      sugg_time <- NULL
    } else if (spacetime) {
      grid_result <- expand.grid(neighb = neighb, maxtime = maxtime)

      sugg_neighb <- min(grid_result$neighb[idx_ok], na.rm = TRUE)
      sugg_neighb <- min(sugg_neighb, bestK)

      idx_m <- idx_ok[grid_result$neighb[idx_ok] == sugg_neighb]
      if (length(idx_m) == 0L) {
        idx_m_all <- which(grid_result$neighb == sugg_neighb)
        idx_s <- idx_m_all[which.min(res[idx_m_all])]
      } else {
        idx_s <- idx_m[which.min(res[idx_m])]
      }
      sugg_time <- grid_result$maxtime[idx_s]
    } else {
      sugg_neighb <- bestK
      sugg_time <- best_T
    }
  }

  list(
    best_neighb    = as.numeric(bestK),
    best_maxtime   = best_T,
    best_estimates = best_estimates,
    res            = unname(res),
    estimates      = estimates,
    sugg_neighb    = as.numeric(sugg_neighb),
    sugg_time      = sugg_time
  )
}
