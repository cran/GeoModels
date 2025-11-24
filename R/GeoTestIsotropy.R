########################################
##parametric boostrap test for isotropy
########################################
GeoTestIsotropy <- function(data, coordx,
                            start, fixed,
                            optimizer = "bobyqa",
                            model = "Gaussian",
                            corrmodel = "Matern",
                            lower = NULL, upper = NULL,
                            B = 1000,
                            likelihood = NULL,
                            type = NULL,
                            copula = NULL,
                            neighb = 5,
                            parallel = TRUE,
                            ncores = NULL,
                            progress = TRUE) {

  ## ====== VALIDAZIONE INPUT INIZIALE ======
  if (!is.matrix(coordx)) coordx <- as.matrix(coordx)
  if (!is.numeric(data))  data   <- as.numeric(data)
  
  # Controlli dimensioni e valori
  if (nrow(coordx) != length(data)) 
    stop("coordx rows (", nrow(coordx), ") must match data length (", length(data), ")")
  if (B < 1) 
    stop("B must be at least 1")
  if (neighb < 1) 
    stop("neighb must be positive")

if(is.null(likelihood)||is.null(type)){
  ## ====== AUTOMATIC LIKELIHOOD SELECTION ======
  if (NN > 2000) {
    likelihood <- "Marginal"
    type <- "Pairwise"
  } else {
    likelihood <- "Full"
    type <- "Standard"
    neighb <- NULL
  }
 }
  
  NN <- nrow(coordx)

  # --- Setup progress handlers UNICA VOLTA ---
  if (progress) {
    progressr::handlers(global = TRUE)
    progressr::handlers(progressr::handler_txtprogressbar(clear = TRUE))
  } else {
    progressr::handlers("void")
  }

  method <- if (NN > 10000) "TB" else "Cholesky"
  corrmodel_used <- corrmodel
  model_used <- model

  ## ====== VALIDAZIONI PARAMETRI ======
  corr_names <- CorrParam(corrmodel_used)
  nuis_names <- NuisParam(model_used)
  expected <- sort(unique(c(corr_names, nuis_names)))
  provided <- sort(unique(c(names(start), names(fixed))))
  if (!identical(expected, provided)) {
    missing <- setdiff(expected, provided)
    extra   <- setdiff(provided, expected)
    msg <- c()
    if (length(missing)) msg <- c(msg, paste0("missing: ", paste(missing, collapse = ", ")))
    if (length(extra))   msg <- c(msg, paste0("extra: ",   paste(extra,   collapse = ", ")))
    stop("Parameter names in start+fixed must match CorrParam(corrmodel)+NuisParam(model). ",
         paste(msg, collapse = " | "))
  }

  # Bounds SOLO per 'start'
  start_names <- sort(names(start))
  norm_bounds <- function(b, which = c("lower","upper")) {
    which <- match.arg(which)
    if (is.null(b)) b <- list()
    b_names <- names(b)
    if (length(b) && is.null(b_names))
      stop(which, " must be a *named* list with names matching start")
    if (length(setdiff(b_names, start_names)))
      stop(which, " contains names not in start: ",
           paste(setdiff(b_names, start_names), collapse = ", "))
    out <- setNames(rep(if (which == "lower") -Inf else Inf, length(start_names)), start_names)
    if (length(b_names)) out[b_names] <- unlist(b, use.names = FALSE)
    as.list(out)
  }
  lower <- norm_bounds(lower, "lower")
  upper <- norm_bounds(upper, "upper")

  ## ====== ANISOTROPIA ======
  aniso_H0    <- list(angle = 0, ratio = 1)
  aniso_start <- list(angle = 0, ratio = 1.2)

  ## ====== HELPERS DI FIT ======
  fit_H0 <- function(Z) {
    GeoFit(
      data = Z, coordx = coordx, corrmodel = corrmodel_used,
      likelihood = likelihood, optimizer = optimizer,
      model = model_used, type = type, neighb = neighb,
      start = start, fixed = fixed,
      copula = copula,
      lower = lower, upper = upper,
      anisopars = aniso_H0
    )
  }
  fit_H1 <- function(Z) {
    GeoFit(
      data = Z, coordx = coordx, corrmodel = corrmodel_used,
      likelihood = likelihood, optimizer = optimizer,
      model = model_used, type = type, neighb = neighb,
      start = start, fixed = fixed,
      copula = copula,
      lower = lower, upper = upper,
      anisopars = aniso_start,
      est.aniso = c(TRUE, TRUE)
    )
  }

  ## ====== FIT OSSERVATO ======
  fit0 <- fit_H0(data)
  fit1 <- fit_H1(data)
  LL0 <- fit0$logCompLik
  LL1 <- fit1$logCompLik
  
  # Controllo e warning per LRT negativo
  if (LL1 < LL0) {
    message("Note: H1 log-likelihood (", round(LL1, 3), 
            ") < H0 log-likelihood (", round(LL0, 3), "), LRT set to 0")
  }
  LRT_obs <- max(0, 2 * (LL1 - LL0))

  ## ====== SIMULAZIONI H0 ======
  par0 <- c(fit0$param, fit0$fixed)
  param_H0 <- as.list(par0[expected])

  if (is.null(copula)) {
    if (method == "Cholesky") {
      sim <- GeoSim(coordx = coordx, corrmodel = corrmodel_used,
                    param = param_H0, model = model_used,
                    nrep = B, anisopars = aniso_H0, progress = FALSE)
    } else if (method == "TB") {
      sim <- GeoSimapprox(coordx = coordx, corrmodel = corrmodel_used, model = model_used, method = "TB",
                          param = param_H0, nrep = B, anisopars = aniso_H0, progress=FALSE)
    } else stop("method must be 'Cholesky' or 'TB'.")
  } else {
    sim <- GeoSimCopula(coordx = coordx, corrmodel = corrmodel_used,
                        method = "Cholesky",
                        param = param_H0,  model = model_used, copula = copula,
                        nrep = B, anisopars = aniso_H0)
  }

  get_rep <- function(simobj, b) {
    X <- simobj$data
    if (is.list(X)) X[[b]] else X[, b]
  }

  ## ====== GESTIONE CORES ROBUSTA ======
  coremax <- parallel::detectCores()
  if (is.na(coremax) || coremax <= 1) {
    parallel <- FALSE
    ncores   <- 1
  } else {
    # Determina ncores ottimale
    if (is.null(ncores)) {
      ncores <- min(coremax - 1, B)
    }
    # Limita a valori sensati
    ncores <- max(1, min(ncores, B, coremax))
    
    # Se ncores = 1, disabilita parallelizzazione
    if (ncores == 1) parallel <- FALSE
  }

  ## ====== BOOTSTRAP ======
  if (!parallel) {
    message(sprintf("Performing %d estimations under H0 and H1 sequentially...", B))

    if (progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
    }
    pb <- progressr::progressor(along = seq_len(B))

    LRT_boot <- vector("numeric", B)
    for (b in seq_len(B)) {
      if (progress) pb(sprintf("b=%d", b))
      
      Zb  <- get_rep(sim, b)
      f0b <- suppressWarnings(try(fit_H0(Zb), silent = TRUE))
      f1b <- suppressWarnings(try(fit_H1(Zb), silent = TRUE))
      
      if (inherits(f0b, "try-error") || inherits(f1b, "try-error") ||
          !is.finite(f0b$logCompLik) || !is.finite(f1b$logCompLik)) {
        LRT_boot[b] <- NA_real_
      } else {
        LRT_boot[b] <- max(0, 2 * (f1b$logCompLik - f0b$logCompLik))
      }
    }
    # Rimuovi solo NA e valori non finiti (max(0,...) già garantisce >= 0)
    LRT_boot <- LRT_boot[is.finite(LRT_boot)]

  } else {
    message(sprintf("Performing %d estimations under H0 and H1 using %d cores...", B, ncores))
    
    future::plan(future::multisession, workers = ncores)
    on.exit(future::plan(future::sequential), add = TRUE)

    if (progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
      pb <- progressr::progressor(along = seq_len(B))
    }

    xx <- foreach::foreach(
      b = seq_len(B),
      .combine = c,
      .options.future = list(
        packages   = "GeoModels",
        seed       = TRUE,
        stdout     = FALSE,
        conditions = "none"
      )
    ) %dofuture% {
      if (progress) pb(sprintf("b=%d", b))
      
      Zb  <- get_rep(sim, b)
      f0b <- suppressWarnings(try(fit_H0(Zb), silent = TRUE))
      f1b <- suppressWarnings(try(fit_H1(Zb), silent = TRUE))

      if (inherits(f0b, "try-error") || inherits(f1b, "try-error") ||
          !is.finite(f0b$logCompLik) || !is.finite(f1b$logCompLik)) {
        as.numeric(NA)
      } else {
        as.numeric(max(0, 2 * (f1b$logCompLik - f0b$logCompLik)))
      }
    }

    # Rimuovi solo NA e valori non finiti
    LRT_boot <- xx[is.finite(xx)]
  }

  ## ====== RISULTATI CON WARNING INFORMATIVO ======
  failed <- B - length(LRT_boot)
  if (failed > 0) {
    warning(sprintf("%d/%d bootstrap replications failed (%.1f%%)", 
                    failed, B, 100 * failed / B))
  }
  
  if (!length(LRT_boot)) {
    warning("No valid bootstrap replications; p-value set to NA")
    pval <- NA_real_
  } else {
    pval <- (1 + sum(LRT_boot >= LRT_obs)) / (length(LRT_boot) + 1)
  }

  ratio_hat <- fit1$anisopars$ratio
  angle_hat <- fit1$anisopars$angle

  cat("\n=== FINAL RESULTS ===\n")
  cat("Observed LRT statistic   =", round(LRT_obs, 5), "\n")
  cat("Valid bootstrap reps     =", length(LRT_boot), "/", B, "\n")
  cat("Estimated anisotropy: ratio =", round(ratio_hat, 4),
      "| angle =", round(angle_hat, 4), "\n")
  cat("p-value                  =", round(pval, 5), "\n\n")

  concl <- if (!is.na(pval) && pval < 0.05)
    "Reject H0 => significant anisotropy detected"
  else
    "Do not reject H0 => data consistent with isotropy"
  cat("Conclusion:", concl, "\n\n")

  invisible(list(statistic = LRT_obs, pvalue = pval,
                 ratio_hat = ratio_hat, angle_hat = angle_hat,
                 corrmodel = corrmodel, model =  model_used,
                 fit_H0 = fit0, fit_H1 = fit1, B_rep = LRT_boot))
}