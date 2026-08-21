####################################################
### Internal utilities for bivariate Gaussian paths
####################################################

.GeoBivariateCoordinateInfo <- function(coordx = NULL, coordy = NULL,
                                        coordz = NULL, coordx_dyn = NULL,
                                        grid = FALSE,
                                        context = "Bivariate model") {
  fail <- function(...) stop(context, ": ", ..., call. = FALSE)

  normalize_block <- function(z, label) {
    if (is.data.frame(z)) z <- as.matrix(z)
    z <- as.matrix(z)
    if (!is.numeric(z) || any(!is.finite(z))) {
      fail(label, " must be a finite numeric coordinate matrix.")
    }
    if (!nrow(z)) fail(label, " must contain at least one location.")
    if (!(ncol(z) %in% c(2L, 3L))) {
      fail(label, " must have two or three coordinate columns.")
    }
    unname(z)
  }

  if (!is.null(coordx_dyn)) {
    if (!is.list(coordx_dyn) || length(coordx_dyn) != 2L) {
      fail("coordx_dyn must be list(coords_variable_1, coords_variable_2).")
    }
    blocks <- list(
      normalize_block(coordx_dyn[[1L]], "coordx_dyn[[1]]"),
      normalize_block(coordx_dyn[[2L]], "coordx_dyn[[2]]")
    )
    if (ncol(blocks[[1L]]) != ncol(blocks[[2L]])) {
      fail("the two coordinate blocks must have the same spatial dimension.")
    }
    ns <- as.integer(vapply(blocks, nrow, integer(1L)))
    return(list(dynamic = TRUE, blocks = blocks,
                coords = do.call(rbind, blocks), ns = ns,
                dimension = ncol(blocks[[1L]])))
  }

  if (isTRUE(grid)) {
    if (is.null(coordx) || is.null(coordy)) {
      fail("coordx and coordy are required when grid = TRUE.")
    }
    coords <- if (is.null(coordz)) {
      as.matrix(expand.grid(coordx, coordy))
    } else {
      as.matrix(expand.grid(coordx, coordy, coordz))
    }
  } else if (!is.null(coordy)) {
    coords <- as.matrix(cbind(coordx, coordy, coordz))
  } else {
    coords <- as.matrix(coordx)
  }

  coords <- normalize_block(coords, "coordinates")
  ns <- rep.int(as.integer(nrow(coords)), 2L)
  list(dynamic = FALSE, blocks = list(coords, coords), coords = coords,
       ns = ns, dimension = ncol(coords))
}


.GeoBivariateDataVector <- function(data, ns, context = "Bivariate model",
                                      allow_na = FALSE) {
  fail <- function(...) stop(context, ": ", ..., call. = FALSE)
  ns <- as.integer(ns)
  if (length(ns) != 2L || any(is.na(ns)) || any(ns < 1L)) {
    fail("invalid bivariate block sizes.")
  }

  if (is.list(data) && !is.data.frame(data)) {
    if (length(data) != 2L) fail("data must be a list of length two.")
    out1 <- as.numeric(data[[1L]])
    out2 <- as.numeric(data[[2L]])
  } else if (is.matrix(data) || is.data.frame(data)) {
    dm <- as.matrix(data)
    if (nrow(dm) == 2L && ns[1L] == ns[2L] && ncol(dm) == ns[1L]) {
      out1 <- as.numeric(dm[1L, ])
      out2 <- as.numeric(dm[2L, ])
    } else if (ncol(dm) == 2L && ns[1L] == ns[2L] && nrow(dm) == ns[1L]) {
      out1 <- as.numeric(dm[, 1L])
      out2 <- as.numeric(dm[, 2L])
    } else if (length(dm) == sum(ns)) {
      vv <- as.numeric(dm)
      out1 <- vv[seq_len(ns[1L])]
      out2 <- vv[ns[1L] + seq_len(ns[2L])]
    } else {
      fail("data dimensions do not match the two coordinate blocks.")
    }
  } else {
    vv <- as.numeric(data)
    if (length(vv) != sum(ns)) {
      fail("data must contain ", sum(ns), " observations in variable-block order.")
    }
    out1 <- vv[seq_len(ns[1L])]
    out2 <- vv[ns[1L] + seq_len(ns[2L])]
  }

  if (length(out1) != ns[1L] || length(out2) != ns[2L]) {
    fail("data lengths do not match the two coordinate blocks.")
  }
  out <- c(out1, out2)
  if (!is.numeric(out) ||
      (!isTRUE(allow_na) && any(!is.finite(out))) ||
      (isTRUE(allow_na) && any(is.infinite(out)))) {
    fail(if (isTRUE(allow_na))
           "data must be numeric and cannot contain infinite values."
         else "data must be finite and numeric.")
  }
  unname(out)
}


.GeoBivariateDesign <- function(X, ns, context = "Bivariate model",
                                allow_null = TRUE) {
  fail <- function(...) stop(context, ": ", ..., call. = FALSE)
  ns <- as.integer(ns)
  if (length(ns) != 2L || any(is.na(ns)) || any(ns < 1L)) {
    fail("invalid bivariate block sizes.")
  }

  if (is.null(X)) {
    if (!allow_null) fail("X is required.")
    return(matrix(1, nrow = sum(ns), ncol = 1L))
  }

  if (is.list(X) && !is.data.frame(X) && !is.matrix(X)) {
    if (length(X) != 2L) fail("X must be a list of length two.")
    X1 <- as.matrix(X[[1L]])
    X2 <- as.matrix(X[[2L]])
    if (nrow(X1) != ns[1L] || nrow(X2) != ns[2L]) {
      fail("the rows of X[[1]] and X[[2]] must match their coordinate blocks.")
    }
    if (ncol(X1) != ncol(X2)) {
      fail("X[[1]] and X[[2]] must have the same number of columns.")
    }
    X <- rbind(X1, X2)
  } else {
    X <- as.matrix(X)
    if (nrow(X) == sum(ns)) {
      # already stacked by variable
    } else if (ns[1L] == ns[2L] && nrow(X) == ns[1L]) {
      X <- rbind(X, X)
    } else {
      fail("X must have ", sum(ns), " rows, or one shared block of ",
           ns[1L], " rows when both variables use the same locations.")
    }
  }

  if (!is.numeric(X) || any(!is.finite(X)) || ncol(X) < 1L) {
    fail("X must be a finite numeric matrix.")
  }
  unname(X)
}


.GeoBivariateValidateFitMeans <- function(start, fixed, X, ns,
                                           context = "Bivariate model") {
  p <- if(is.null(X)) 1L else ncol(.GeoBivariateDesign(X, ns, context=context))
  expected <- c(.GeoMean_names(p, "mean_1"),
                .GeoMean_names(p, "mean_2"))
  start_names <- unique(c(.GeoMean_supplied_names(start, "mean_1"),
                          .GeoMean_supplied_names(start, "mean_2")))
  fixed_names <- unique(c(.GeoMean_supplied_names(fixed, "mean_1"),
                          .GeoMean_supplied_names(fixed, "mean_2")))

  overlap <- intersect(start_names, fixed_names)
  if(length(overlap)) {
    stop(context, ": a mean coefficient cannot be supplied in both start and fixed: ",
         paste(overlap, collapse=", "), ".", call.=FALSE)
  }
  supplied <- unique(c(start_names, fixed_names))
  unexpected <- setdiff(supplied, expected)
  if(length(unexpected)) {
    stop(context, ": for ", p, " design column(s), mean coefficients can only be named ",
         paste(expected, collapse=", "), ". Unexpected: ",
         paste(unexpected, collapse=", "), ".", call.=FALSE)
  }

  values <- c(as.list(start), as.list(fixed))
  for(nm in supplied) {
    value <- values[[nm]]
    if(length(value) != 1L || !is.numeric(value) || !is.finite(value)) {
      stop(context, ": mean coefficient '", nm,
           "' must be a finite numeric scalar.", call.=FALSE)
    }
  }
  invisible(expected)
}


.GeoBivariateBetas <- function(param, p, variable,
                               context = "Bivariate model") {
  if (!(variable %in% c(1L, 2L))) {
    stop(context, ": variable must be 1 or 2.", call. = FALSE)
  }
  .GeoMean_beta(param, p = p, prefix = paste0("mean_", variable))
}


.GeoBivariateValidateMeanParameters <- function(param, p,
                                                   context = "Bivariate model") {
  expected <- c(.GeoMean_names(p, "mean_1"),
                .GeoMean_names(p, "mean_2"))
  supplied <- unique(c(
    .GeoMean_supplied_names(param, "mean_1"),
    .GeoMean_supplied_names(param, "mean_2")
  ))
  unexpected <- setdiff(supplied, expected)

  ## Intercept-only bivariate models historically default both means to zero.
  ## Preserve that convention, but require every coefficient when X has more
  ## than one column.
  param_use <- as.list(param)
  if (p == 1L) {
    if (!("mean_1" %in% names(param_use))) param_use$mean_1 <- 0
    if (!("mean_2" %in% names(param_use))) param_use$mean_2 <- 0
  }
  missing <- setdiff(expected, names(param_use))

  if (length(missing) || length(unexpected)) {
    msg <- paste0(context, ": the mean specification must contain ",
                  paste(expected, collapse = ", "), ".")
    if (length(missing))
      msg <- paste0(msg, " Missing: ", paste(missing, collapse = ", "), ".")
    if (length(unexpected))
      msg <- paste0(msg, " Unexpected: ", paste(unexpected, collapse = ", "), ".")
    stop(msg, call. = FALSE)
  }
  list(
    beta1 = .GeoBivariateBetas(param_use, p, 1L, context),
    beta2 = .GeoBivariateBetas(param_use, p, 2L, context)
  )
}


.GeoBivariateObservedMean <- function(param, X, ns,
                                      context = "Bivariate model") {
  X <- .GeoBivariateDesign(X, ns, context = context)
  p <- ncol(X)
  beta <- .GeoBivariateValidateMeanParameters(param, p, context)
  b1 <- beta$beta1
  b2 <- beta$beta2
  i1 <- seq_len(ns[1L])
  i2 <- ns[1L] + seq_len(ns[2L])
  as.numeric(c(X[i1, , drop = FALSE] %*% b1,
               X[i2, , drop = FALSE] %*% b2))
}


.GeoBivariatePredictionMean <- function(param, Xloc, nloc, which,
                                        p, Mloc = NULL,
                                        context = "Bivariate model") {
  if (!is.null(Mloc)) {
    Mloc <- as.numeric(Mloc)
    if (length(Mloc) == 1L) Mloc <- rep(Mloc, nloc)
    if (length(Mloc) != nloc || any(!is.finite(Mloc))) {
      stop(context, ": Mloc must contain one finite value per prediction location.",
           call. = FALSE)
    }
    return(Mloc)
  }
  Xloc <- .GeoMean_design(Xloc, n = nloc, p = p, name = "Xloc")
  beta <- .GeoBivariateValidateMeanParameters(param, p, context)
  beta_target <- if (as.integer(which) == 1L) beta$beta1 else beta$beta2
  as.numeric(Xloc %*% beta_target)
}


.GeoBivariateBlockDesign <- function(X, ns) {
  X <- .GeoBivariateDesign(X, ns)
  p <- ncol(X)
  out <- matrix(0, nrow = sum(ns), ncol = 2L * p)
  i1 <- seq_len(ns[1L])
  i2 <- ns[1L] + seq_len(ns[2L])
  out[i1, seq_len(p)] <- X[i1, , drop = FALSE]
  out[i2, p + seq_len(p)] <- X[i2, , drop = FALSE]
  colnames(out) <- c(.GeoMean_names(p, "mean_1"),
                     .GeoMean_names(p, "mean_2"))
  out
}


.GeoBivariateTargetBlockDesign <- function(Xloc, nloc, p, which) {
  Xloc <- .GeoMean_design(Xloc, n = nloc, p = p, name = "Xloc")
  out <- matrix(0, nrow = nloc, ncol = 2L * p)
  if (which == 1L) out[, seq_len(p)] <- Xloc
  else out[, p + seq_len(p)] <- Xloc
  colnames(out) <- c(.GeoMean_names(p, "mean_1"),
                     .GeoMean_names(p, "mean_2"))
  out
}
