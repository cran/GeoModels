## Internal utilities shared by the parametric-bootstrap test functions.
## They centralize runtime plumbing only; each public test keeps its own
## statistical model, fitting strategy, and result construction.

.GeoTestPreserveRuntime <- function(gc_on_exit = FALSE) {
  old_plan <- if (requireNamespace("future", quietly = TRUE)) future::plan() else NULL
  old_handlers <- if (requireNamespace("progressr", quietly = TRUE)) progressr::handlers() else NULL
  old_progress_global <- if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::handlers(global = NA)
  } else {
    NULL
  }

  function() {
    if (!is.null(old_plan)) try(future::plan(old_plan), silent = TRUE)
    if (requireNamespace("progressr", quietly = TRUE)) {
      try(progressr::handlers(old_handlers), silent = TRUE)
      if (!is.null(old_progress_global)) {
        ## Do not wrap global handler changes in try()/tryCatch(); progressr
        ## must call base::globalCallingHandlers() without handlers on stack.
        progressr::handlers(global = old_progress_global)
      }
    }
    if (isTRUE(gc_on_exit)) gc(verbose = FALSE, full = TRUE)
    invisible(NULL)
  }
}

.GeoTestSetSeed <- function(seed) {
  if (is.null(seed)) return(function() invisible(NULL))
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    stop("seed must be NULL or a single finite numeric value", call. = FALSE)
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  set.seed(as.integer(seed))

  function() {
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    invisible(NULL)
  }
}

.GeoTestConfigureProgress <- function(progress, warn_missing = TRUE) {
  use_progressr <- isTRUE(progress) && requireNamespace("progressr", quietly = TRUE)
  if (isTRUE(progress) && !use_progressr && isTRUE(warn_missing)) {
    warning(
      "progress=TRUE but progressr is not available; progress bars disabled.",
      call. = FALSE
    )
  }
  if (use_progressr) {
    progressr::handlers(global = TRUE)
    progressr::handlers(progressr::handler_txtprogressbar(clear = TRUE))
  } else if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::handlers("void")
  }
  use_progressr
}

.GeoTestCoreConfig <- function(parallel, ncores, B, strict_positive = TRUE) {
  ## strict_positive is retained for internal-call compatibility; all public
  ## ncores values now share the same positive-integer validation.
  ncores <- .GeoResolveWorkers(parallel, ncores, n_jobs = B)
  list(parallel = isTRUE(parallel) && ncores > 1L, ncores = ncores)
}

.GeoTestNormalizeBounds <- function(b, which = c("lower", "upper"),
                                    start_names, require_list = TRUE) {
  which <- match.arg(which)
  if (is.null(b)) b <- list()
  if (isTRUE(require_list) && !is.list(b)) {
    stop(which, " must be a named list", call. = FALSE)
  }

  b_names <- names(b)
  if (length(b) > 0L && (is.null(b_names) || any(b_names == ""))) {
    stop(which, " must be a named list with names matching start", call. = FALSE)
  }
  unknown <- setdiff(b_names, start_names)
  if (length(unknown) > 0L) {
    stop(which, " contains names not in start: ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  base_val <- if (which == "lower") -Inf else Inf
  out <- setNames(rep(base_val, length(start_names)), start_names)
  if (length(b) > 0L) out[b_names] <- unlist(b, use.names = FALSE)
  as.list(out)
}

.GeoTestLogLikValue <- function(fit_obj) {
  if (inherits(fit_obj, "try-error") || is.null(fit_obj)) return(NA_real_)
  val <- if (!is.null(fit_obj$logCompLik)) fit_obj$logCompLik else fit_obj$logLik
  val <- as.numeric(val)
  if (length(val) != 1L || !is.finite(val)) NA_real_ else val
}

.GeoTestSamePairs <- function(fit0, fit1, type, require_indices = FALSE) {
  if (!identical(type, "Pairwise")) return(TRUE)
  if (isTRUE(require_indices) &&
      (is.null(fit0$rowidx) || is.null(fit0$colidx) ||
       is.null(fit1$rowidx) || is.null(fit1$colidx))) {
    return(FALSE)
  }
  identical(fit0$rowidx, fit1$rowidx) && identical(fit0$colidx, fit1$colidx)
}

.GeoTestGetRep <- function(sim, b) {
  x <- sim$data
  if (is.list(x)) return(x[[b]])
  if (is.null(dim(x))) return(x)
  x[, b]
}

.GeoTestSequentialApply <- function(sim, B, FUN, progress = FALSE) {
  run <- function() {
    pb <- if (isTRUE(progress)) progressr::progressor(along = seq_len(B)) else NULL
    out <- vector("list", B)
    for (b in seq_len(B)) {
      out[[b]] <- FUN(.GeoTestGetRep(sim, b))
      if (!is.null(pb)) pb(sprintf("b=%d", b))
    }
    out
  }
  if (isTRUE(progress)) progressr::with_progress(run()) else run()
}

.GeoTestSaveReplicates <- function(sim, B, prefix) {
  tag <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
  vapply(seq_len(B), function(b) {
    path <- file.path(tempdir(), sprintf("%s_%s_%05d.rds", prefix, tag, b))
    saveRDS(.GeoTestGetRep(sim, b), path, compress = FALSE)
    path
  }, character(1))
}

.GeoTestParallelFileApply <- function(data_files, FUN, ncores,
                                      progressor = NULL,
                                      future_packages = NULL,
                                      globals = NULL,
                                      progress_every = FALSE,
                                      error_value = NULL) {
  old_plan <- future::plan()
  on.exit(try(future::plan(old_plan), silent = TRUE), add = TRUE)
  future::plan(future::multisession, workers = ncores)

  if (!is.null(globals)) {
    old_max_size <- getOption("future.globals.maxSize")
    on.exit(options(future.globals.maxSize = old_max_size), add = TRUE)
    required_size <- sum(vapply(
      globals, function(x) as.numeric(utils::object.size(x)), numeric(1)
    ))
    options(future.globals.maxSize = max(500 * 1024^2, 1.5 * required_size))
  }

  ## Use a bounded queue of individual futures instead of future_lapply().
  ## This lets the main R process detect each completed bootstrap refit and
  ## advance the package-standard progress bar immediately.  In particular,
  ## progress updates are not buffered until an entire future.apply chunk or
  ## bootstrap batch has finished.
  n_jobs <- length(data_files)
  if (n_jobs == 0L) return(vector("list", 0L))
  workers <- max(1L, min(as.integer(ncores), n_jobs))
  ans <- vector("list", n_jobs)
  active <- list()
  next_job <- 1L

  launch_one <- function(b) {
    path <- data_files[[b]]
    future::future(
      {
        FUN(readRDS(path))
      },
      globals = list(FUN = FUN, path = path),
      packages = future_packages,
      seed = TRUE
    )
  }

  launch_available <- function() {
    while (length(active) < workers && next_job <= n_jobs) {
      b <- next_job
      next_job <<- next_job + 1L
      active[[as.character(b)]] <<- launch_one(b)
    }
    invisible(NULL)
  }

  launch_available()
  while (length(active) > 0L) {
    ready <- vapply(active, future::resolved, logical(1))
    if (!any(ready)) {
      Sys.sleep(0.05)
      next
    }

    ids <- names(active)[ready]
    for (id in ids) {
      b <- as.integer(id)
      out <- tryCatch(
        future::value(active[[id]]),
        error = function(e) {
          if (!is.null(error_value)) return(error_value)
          list(
            valid = FALSE,
            lrt = NA_real_,
            ratio = NA_real_,
            angle = NA_real_,
            boundary_selected = FALSE,
            h1_attempts = 0L,
            h1_valid_attempts = 0L,
            error = conditionMessage(e)
          )
        }
      )
      ans[[b]] <- out
      active[[id]] <- NULL

      if (!is.null(progressor)) {
        ## Most bootstrap tests return a structured list with valid/lrt/ratio,
        ## but GeoTestIndependence() and GeoTestsupp_space() return named
        ## numeric vectors.  Never use $ on an arbitrary worker result.
        progress_ok <- isTRUE(progress_every)
        if (!progress_ok && is.list(out)) {
          progress_ok <- isTRUE(out$valid) &&
            length(out$lrt) == 1L && is.finite(out$lrt) &&
            length(out$ratio) == 1L && is.finite(out$ratio)
        }
        if (progress_ok) progressor(sprintf("b=%d", b))
      }
      launch_available()
    }
  }

  ans
}
