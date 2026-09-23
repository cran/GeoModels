#' Configure native OpenMP acceleration
#'
#' Configure or inspect the optional OpenMP acceleration used by the
#' univariate Gaussian spatial pairwise marginal and conditional likelihoods.
#' The parallel native path covers the Matérn correlation model for any
#' positive smoothness. Exact half-integer smoothness uses elementary formulas;
#' general smoothness uses an internal allocation-free Bessel-K kernel that
#' does not call the R API from OpenMP workers. Other correlation/likelihood
#' combinations use the serial implementation.
#'
#' The default is one thread. This opt-in default avoids oversubscription when
#' GeoModels is called inside outer parallel workflows such as future-based
#' simulation, bootstrap, or cross-validation. Use either outer R-level
#' parallelism or native OpenMP parallelism, rather than both simultaneously.
#'
#' @param threads Optional integer. `1` selects the serial implementation;
#'   values greater than one request that many OpenMP threads; `0` requests
#'   the OpenMP runtime default. If omitted, the current setting is retained.
#' @param min_pairs Optional positive integer giving the minimum number of
#'   retained pairs required before the OpenMP path is used. The default is
#'   50000.
#'
#' @return Returns a list with components `compiled`, `threads`,
#'   `max_threads`, and `min_pairs`.
#' @export
#'
#' @examples
#' GeoOMP()
#' GeoOMP(threads = 1L)
#' \dontrun{GeoOMP(threads = 8L)}
GeoOMP <- function(threads = NULL, min_pairs = NULL) {

  if (!is.null(threads)) {
    if (length(threads) != 1L || is.na(threads) || !is.finite(threads) ||
        threads < 0 || threads != as.integer(threads)) {
      stop("threads must be one non-negative integer", call. = FALSE)
    }
    .C(
      "GeoModels_set_omp_threads",
      n = as.integer(threads),
      PACKAGE = "GeoModels"
    )
  }

  if (!is.null(min_pairs)) {
    if (length(min_pairs) != 1L || is.na(min_pairs) || !is.finite(min_pairs) ||
        min_pairs < 1 || min_pairs != as.integer(min_pairs)) {
      stop("min_pairs must be one positive integer", call. = FALSE)
    }
    .C(
      "GeoModels_set_omp_min_pairs",
      n = as.integer(min_pairs),
      PACKAGE = "GeoModels"
    )
  }

  ans <- .C(
    "GeoModels_openmp_info",
    compiled = integer(1L),
    threads = integer(1L),
    max_threads = integer(1L),
    min_pairs = integer(1L),
    PACKAGE = "GeoModels"
  )

  out <- list(
    compiled = as.logical(ans$compiled),
    threads = as.integer(ans$threads),
    max_threads = as.integer(ans$max_threads),
    min_pairs = as.integer(ans$min_pairs)
  )

  if (!out$compiled && !is.null(threads) && threads > 1L) {
    warning(
      "GeoModels was compiled without OpenMP support; one thread will be used.",
      call. = FALSE
    )
  }

  out
}
