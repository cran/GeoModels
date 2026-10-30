## Internal parallel-worker utilities shared across GeoModels.

.GeoValidateNcores <- function(ncores, context = "ncores") {
  if (is.null(ncores)) return(NULL)
  tol <- sqrt(.Machine$double.eps)
  if (!is.numeric(ncores) || length(ncores) != 1L || !is.finite(ncores) ||
      ncores < 1 || abs(ncores - round(ncores)) > tol) {
    stop(context, " must be NULL or a positive integer", call. = FALSE)
  }
  as.integer(round(ncores))
}

.GeoResolveWorkers <- function(parallel = TRUE, ncores = NULL, n_jobs = Inf,
                               auto_max = 6L, context = "ncores") {
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel))
    stop("parallel must be TRUE or FALSE", call. = FALSE)

  ncores <- .GeoValidateNcores(ncores, context)
  ## Respect scheduler/cgroup quotas and R's check-time core limits, rather
  ## than treating all hardware CPUs as available to this R process.
  coremax <- suppressWarnings(parallelly::availableCores())
  if (is.na(coremax) || coremax < 1L) coremax <- 1L
  coremax <- as.integer(coremax)

  jobs <- suppressWarnings(as.numeric(n_jobs)[1L])
  if (is.na(jobs) || jobs <= 0) jobs <- 1
  if (!is.finite(jobs)) jobs <- Inf

  if (!isTRUE(parallel) || coremax <= 1L || jobs <= 1) return(1L)

  if (is.null(ncores)) {
    auto_max <- suppressWarnings(as.numeric(auto_max)[1L])
    if (!is.finite(auto_max) || auto_max < 1) auto_max <- 6L
    workers <- min(as.integer(auto_max), max(1L, coremax - 1L))
  } else {
    ## An explicit request may use all detected cores; the automatic path
    ## leaves one core free and is capped at six workers package-wide.
    workers <- min(ncores, coremax)
  }

  if (is.finite(jobs)) workers <- min(workers, as.integer(floor(jobs)))
  as.integer(max(1L, workers))
}
