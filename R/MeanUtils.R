####################################################
### Internal utilities for mean specifications
####################################################

.GeoMean_names <- function(p, prefix = "mean") {
  if (length(p) != 1L || is.na(p) || p < 1L || p != as.integer(p)) {
    stop("The number of mean coefficients must be a positive integer.", call. = FALSE)
  }
  p <- as.integer(p)
  if (p == 1L) return(prefix)
  c(prefix, paste0(prefix, seq_len(p - 1L)))
}


.GeoMean_supplied_names <- function(param, prefix = "mean") {
  if (is.null(names(param))) return(character(0))
  pat <- paste0("^", prefix, "[0-9]*$")
  unique(names(param)[grepl(pat, names(param))])
}

.GeoMean_zero <- function(param, p, prefix = "mean") {
  for (nm in .GeoMean_names(p, prefix)) param[[nm]] <- 0
  param
}

.GeoMean_intercept_param <- function(param, value = 0, prefix = "mean") {
  for (nm in .GeoMean_supplied_names(param, prefix)) param[[nm]] <- NULL
  param[[prefix]] <- value
  param
}

.GeoMean_external <- function(param, n = NULL, name = "mean") {
  ## `fixed` and other parameter containers are not guaranteed to remain
  ## lists.  Several fitted objects store them as named atomic vectors, for
  ## which param[[name]] throws "subscript out of bounds" when the requested
  ## name is absent.  Absence of `mean` means that no external mean was
  ## supplied and must therefore return NULL.
  if (is.null(param) || !length(param) || is.null(names(param)) ||
      !(name %in% names(param))) {
    return(NULL)
  }

  value <- if (is.list(param)) {
    param[[name]]
  } else {
    param[names(param) == name]
  }
  if (is.null(value) || !length(value)) return(NULL)

  out <- as.numeric(unlist(value, use.names = FALSE))
  if (length(out) <= 1L) return(NULL)
  if (any(!is.finite(out))) {
    stop("The external mean vector contains non-finite values.", call. = FALSE)
  }
  if (!is.null(n) && length(out) != n) {
    stop("The external mean vector must have one value per observation.", call. = FALSE)
  }
  out
}

.GeoMean_beta <- function(param, p, prefix = "mean") {
  nm <- .GeoMean_names(p, prefix)
  if (is.null(names(param))) {
    stop("Mean coefficients must be supplied as named parameters.", call. = FALSE)
  }
  missing_nm <- nm[!nm %in% names(param)]
  if (length(missing_nm)) {
    stop("Missing mean coefficient(s): ", paste(missing_nm, collapse = ", "), call. = FALSE)
  }
  vals <- lapply(nm, function(z) param[[z]])
  bad <- vapply(vals, function(z) length(z) != 1L || !is.numeric(z) || !is.finite(z), logical(1L))
  if (any(bad)) {
    stop("Each mean coefficient must be a finite numeric scalar.", call. = FALSE)
  }
  out <- as.numeric(unlist(vals, use.names = FALSE))
  names(out) <- nm
  out
}

.GeoMean_design <- function(X, n, p = NULL, name = "X", allow_null = TRUE) {
  if (length(n) != 1L || is.na(n) || n < 1L || n != as.integer(n)) {
    stop("Invalid number of rows requested for ", name, ".", call. = FALSE)
  }
  n <- as.integer(n)

  if (is.null(X)) {
    if (!allow_null) stop(name, " is required.", call. = FALSE)
    if (!is.null(p) && p != 1L) {
      stop(name, " is required when more than one mean coefficient is used.", call. = FALSE)
    }
    return(matrix(1, nrow = n, ncol = 1L))
  }

  if (is.data.frame(X)) X <- as.matrix(X)
  if (is.list(X) && !is.matrix(X)) X <- do.call(rbind, lapply(X, as.matrix))

  if (is.vector(X) && !is.list(X)) {
    if (!is.null(p) && p == 1L && length(X) == n) {
      X <- matrix(X, nrow = n, ncol = 1L)
    } else if (n == 1L && (is.null(p) || length(X) == p)) {
      X <- matrix(X, nrow = 1L)
    } else {
      stop(name, " supplied as a vector is ambiguous; use an explicit matrix.", call. = FALSE)
    }
  } else {
    X <- as.matrix(X)
  }

  if (!is.numeric(X) || any(!is.finite(X))) {
    stop(name, " must be a finite numeric matrix.", call. = FALSE)
  }
  if (nrow(X) != n) {
    stop(name, " must have ", n, " rows.", call. = FALSE)
  }
  if (!is.null(p) && ncol(X) != p) {
    stop(name, " must have ", p, " columns, one for each mean coefficient.", call. = FALSE)
  }
  unname(X)
}

.GeoMean_vector <- function(param, X = NULL, n, external = NULL,
                            prefix = "mean", xname = "X") {
  if (!is.null(external)) {
    external <- as.numeric(external)
    if (length(external) == 1L) external <- rep(external, n)
    if (length(external) != n || any(!is.finite(external))) {
      stop("The fixed mean must contain one finite value per location.", call. = FALSE)
    }
    return(external)
  }

  if (is.null(X)) {
    X <- matrix(1, nrow = n, ncol = 1L)
    p <- 1L
  } else {
    X <- .GeoMean_design(X, n = n, name = xname)
    p <- ncol(X)
  }
  beta <- .GeoMean_beta(param, p = p, prefix = prefix)
  as.numeric(X %*% beta)
}

.GeoMean_validate_fit <- function(start, fixed, X, nobs) {
  start_mean <- .GeoMean_supplied_names(start)
  fixed_mean <- .GeoMean_supplied_names(fixed)
  duplicated_mean <- intersect(start_mean, fixed_mean)
  if (length(duplicated_mean)) {
    stop(
      "A mean coefficient cannot be supplied in both start and fixed: ",
      paste(duplicated_mean, collapse = ", "), ".",
      call. = FALSE
    )
  }

  external <- .GeoMean_external(fixed, n = nobs)
  if (!is.null(external)) {
    if (!is.null(X)) {
      stop("A site-specific fixed mean and X cannot be used together.", call. = FALSE)
    }
    if (length(start_mean)) {
      stop(
        "A site-specific fixed mean cannot be combined with estimated mean coefficients in start.",
        call. = FALSE
      )
    }
    extra_fixed <- setdiff(fixed_mean, "mean")
    if (length(extra_fixed)) {
      stop(
        "A site-specific fixed mean cannot be combined with additional mean coefficients: ",
        paste(extra_fixed, collapse = ", "), ".",
        call. = FALSE
      )
    }
    return(list(X = NULL, external = external))
  }

  supplied <- unique(c(start_mean, fixed_mean))

  if (!is.null(X)) {
    X <- .GeoMean_design(X, n = nobs, name = "X")
    expected <- .GeoMean_names(ncol(X))
    unexpected <- setdiff(supplied, expected)
    if (length(unexpected)) {
      stop(
        "For ncol(X) = ", ncol(X), ", regression coefficients can only be named: ",
        paste(expected, collapse = ", "), ". Unexpected: ",
        paste(unexpected, collapse = ", "), ".",
        call. = FALSE
      )
    }
    if (length(supplied)) {
      all_values <- c(start, fixed)
      vals <- lapply(supplied, function(nm) all_values[[nm]])
      bad <- vapply(
        vals,
        function(z) length(z) != 1L || !is.numeric(z) || !is.finite(z),
        logical(1L)
      )
      if (any(bad)) {
        stop("Supplied mean coefficients must be finite numeric scalars.", call. = FALSE)
      }
    }
  } else {
    unexpected <- setdiff(supplied, "mean")
    if (length(unexpected)) {
      stop(
        "Multiple mean coefficients require an explicit design matrix X; ",
        "without X the model is intercept-only and uses parameter 'mean'.",
        call. = FALSE
      )
    }
    if (length(supplied)) {
      all_values <- c(start, fixed)
      value <- all_values[["mean"]]
      if (length(value) != 1L || !is.numeric(value) || !is.finite(value)) {
        stop("The intercept parameter 'mean' must be a finite numeric scalar.", call. = FALSE)
      }
    }
  }

  list(X = X, external = NULL)
}
