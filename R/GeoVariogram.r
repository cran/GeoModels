####################################################
### File name: GeoVariogram.r
####################################################

.GeoVariogram_bin_values <- function(lags, values, breaks) {
  nbin <- length(breaks) - 1L
  out_count <- integer(max(0L, nbin))
  out_sum <- numeric(max(0L, nbin))
  if (nbin <= 0L || !length(lags)) {
    return(list(count = out_count, sum = out_sum,
                mean = rep(NA_real_, max(0L, nbin))))
  }

  lags <- as.numeric(lags)
  values <- as.numeric(values)
  ok <- is.finite(lags) & is.finite(values)
  if (!any(ok)) {
    return(list(count = out_count, sum = out_sum,
                mean = rep(NA_real_, nbin)))
  }
  lags <- lags[ok]
  values <- values[ok]

  ib <- findInterval(lags, breaks, rightmost.closed = TRUE)
  ok <- ib >= 1L & ib <= nbin
  if (!any(ok)) {
    return(list(count = out_count, sum = out_sum,
                mean = rep(NA_real_, nbin)))
  }
  ib <- ib[ok]
  values <- values[ok]

  out_count <- tabulate(ib, nbins = nbin)
  rr <- rowsum(values, group = ib, reorder = FALSE, na.rm = TRUE)
  ii <- suppressWarnings(as.integer(rownames(rr)))
  out_sum[ii] <- rr[, 1L]
  out_mean <- rep(NA_real_, nbin)
  nz <- out_count > 0L
  out_mean[nz] <- out_sum[nz] / out_count[nz]

  list(count = out_count, sum = out_sum, mean = out_mean)
}

.GeoVariogram_bin_surface <- function(lags, lagt, values,
                                      spatial_breaks, temporal_breaks) {
  nh <- length(spatial_breaks) - 1L
  nt <- length(temporal_breaks) - 1L
  count <- matrix(0L, nrow = max(0L, nh), ncol = max(0L, nt))
  sums <- matrix(0, nrow = max(0L, nh), ncol = max(0L, nt))
  means <- matrix(NA_real_, nrow = max(0L, nh), ncol = max(0L, nt))

  if (nh <= 0L || nt <= 0L || !length(lags)) {
    return(list(count = count, sum = sums, mean = means))
  }

  lags <- as.numeric(lags)
  lagt <- as.numeric(lagt)
  values <- as.numeric(values)
  ok <- is.finite(lags) & is.finite(lagt) & is.finite(values)
  if (!any(ok)) return(list(count = count, sum = sums, mean = means))

  hb <- findInterval(lags[ok], spatial_breaks, rightmost.closed = TRUE)
  tb <- findInterval(lagt[ok], temporal_breaks, rightmost.closed = TRUE)
  values <- values[ok]
  ok2 <- hb >= 1L & hb <= nh & tb >= 1L & tb <= nt
  if (!any(ok2)) return(list(count = count, sum = sums, mean = means))

  hb <- hb[ok2]
  tb <- tb[ok2]
  values <- values[ok2]
  grp <- (hb - 1L) * nt + tb

  cc <- tabulate(grp, nbins = nh * nt)
  ss <- numeric(nh * nt)
  rr <- rowsum(values, group = grp, reorder = FALSE, na.rm = TRUE)
  ii <- suppressWarnings(as.integer(rownames(rr)))
  ss[ii] <- rr[, 1L]

  count <- matrix(cc, nrow = nh, ncol = nt, byrow = TRUE)
  sums <- matrix(ss, nrow = nh, ncol = nt, byrow = TRUE)
  nz <- count > 0L
  means[nz] <- sums[nz] / count[nz]
  list(count = count, sum = sums, mean = means)
}

.GeoVariogram_make_breaks <- function(lags, numbins) {
  lags <- as.numeric(lags)
  lags <- lags[is.finite(lags)]
  if (!length(lags)) return(numeric(0L))

  lo <- min(lags)
  hi <- max(lags)
  if (hi <= lo) {
    delta <- max(1, abs(lo)) * sqrt(.Machine$double.eps)
    lo <- max(0, lo - delta)
    hi <- hi + delta
  }
  seq(lo, hi, length.out = as.integer(numbins))
}

.GeoVariogram_time_bins <- function(coordt, maxtime = NULL, numbins_t = NULL) {
  coordt <- as.numeric(coordt)
  if (length(coordt) < 2L || any(!is.finite(coordt))) {
    stop("At least two finite temporal coordinates are required")
  }

  tt <- sort(coordt)
  dtt <- diff(tt)
  tol_t <- sqrt(.Machine$double.eps) * max(1, max(abs(tt)))
  if (any(dtt <= tol_t)) {
    stop("Temporal coordinates must be distinct (up to numerical tolerance)")
  }

  max_available <- tt[length(tt)] - tt[1L]
  if (is.null(maxtime)) {
    maxlag <- max_available
  } else {
    if (!is.numeric(maxtime) || length(maxtime) != 1L ||
        !is.finite(maxtime) || maxtime <= 0) {
      stop("maxtime must be a positive finite numeric scalar")
    }
    maxlag <- min(as.numeric(maxtime), max_available)
  }
  if (!(maxlag > 0)) {
    stop("maxtime does not include any positive temporal lag")
  }

  dt0 <- stats::median(dtt)
  tol_d <- sqrt(.Machine$double.eps) * max(1, abs(dt0), max(abs(dtt)))
  regular <- all(abs(dtt - dt0) <= tol_d)
  if (maxlag + tol_d < min(dtt)) {
    stop("maxtime does not include any positive temporal lag")
  }

  if (regular) {
    ## For regular time coordinates the attainable positive lags are integer
    ## multiples of the common spacing.  Construct them directly: unlike the
    ## historical outer(coordt, coordt, '-') implementation this is O(T), not
    ## O(T^2), and it does not split numerically equal decimal lags.
    nlag <- floor((maxlag + tol_d) / dt0)
    if (nlag < 1L) {
      stop("maxtime does not include any positive temporal lag")
    }
    centers <- dt0 * seq_len(nlag)
    centers <- centers[centers <= maxlag + tol_d]
    if (!length(centers)) {
      stop("maxtime does not include any positive temporal lag")
    }

    if (length(centers) == 1L) {
      upper <- max(maxlag, centers[1L] + 0.5 * dt0)
      breaks <- c(0, upper)
    } else {
      mids <- (centers[-1L] + centers[-length(centers)]) / 2
      upper <- max(maxlag, centers[length(centers)] + 0.5 * dt0)
      breaks <- c(0, mids, upper)
    }
  } else {
    ## Irregular times can generate O(T^2) distinct lags.  Do not enumerate
    ## them.  Use a controlled number of temporal bins instead.
    if (is.null(numbins_t)) numbins_t <- min(12L, length(coordt) - 1L)
    if (!is.numeric(numbins_t) || length(numbins_t) != 1L ||
        !is.finite(numbins_t) || numbins_t < 1 || numbins_t %% 1 != 0) {
      stop("numbins_t must be a positive integer")
    }
    ntbin <- min(as.integer(numbins_t), length(coordt) - 1L)
    breaks <- seq(0, maxlag, length.out = ntbin + 1L)
    centers <- breaks[-length(breaks)] + diff(breaks) / 2
  }

  list(breaks = as.numeric(breaks), centers = as.numeric(centers),
       regular = regular, maxlag = maxlag)
}

.GeoVariogram_validate_data <- function(data, context = "GeoVariogram") {
  if (is.list(data)) {
    if (!all(vapply(data, is.numeric, logical(1L)))) {
      stop(context, ": data must be numeric", call. = FALSE)
    }
    z <- unlist(data, use.names = FALSE)
  } else {
    if (!is.numeric(data)) stop(context, ": data must be numeric", call. = FALSE)
    z <- as.numeric(data)
  }
  if (any(is.infinite(z))) {
    stop(context, ": data may contain NA/NaN but not Inf or -Inf", call. = FALSE)
  }
  invisible(TRUE)
}

.GeoVariogram_validate_dynamic_data <- function(data, coordx_dyn,
                                                 context = "GeoVariogram") {
  if (!is.list(data) || length(data) != length(coordx_dyn)) {
    stop(context, ": data must be a list with one component per temporal instant when coordx_dyn is used.",
         call. = FALSE)
  }
  sizes <- vapply(coordx_dyn, nrow, integer(1L))
  for (i in seq_along(data)) {
    if (!is.numeric(data[[i]]) || any(is.infinite(data[[i]]))) {
      stop(context, ": data[[", i, "]] must be numeric and may contain NA/NaN but not Inf or -Inf.",
           call. = FALSE)
    }
    if (length(data[[i]]) != sizes[i]) {
      stop(context, ": data[[", i, "]] must contain one value per row of coordx_dyn[[", i, "]].",
           call. = FALSE)
    }
  }
  invisible(sizes)
}

.GeoVariogram_normalize_grid <- function(data, coordx, coordy = NULL, coordz = NULL,
                                         coordt = NULL, bivariate = FALSE) {
  if (is.matrix(coordx)) {
    coords <- as.matrix(coordx)
  } else {
    if (is.null(coordy)) {
      stop("grid=TRUE requires coordy when coordx is not an explicit coordinate matrix")
    }
    axes <- list(coordx, coordy)
    if (!is.null(coordz)) axes[[3L]] <- coordz
    coords <- as.matrix(do.call(expand.grid, axes))
  }
  storage.mode(coords) <- "double"
  nspace <- nrow(coords)

  if (isTRUE(bivariate)) {
    if (is.list(data)) {
      if (length(data) != 2L ||
          any(vapply(data, length, integer(1L)) != nspace)) {
        stop("For grid=TRUE and bivariate=TRUE, data must contain two fields on the full spatial grid")
      }
      data_out <- lapply(data, as.numeric)
    } else if (is.matrix(data) &&
               ((nrow(data) == 2L && ncol(data) == nspace) ||
                (ncol(data) == 2L && nrow(data) == nspace))) {
      data_out <- data
    } else if (!is.null(dim(data)) && tail(dim(data), 1L) == 2L &&
               prod(head(dim(data), -1L)) == nspace) {
      zz <- matrix(as.numeric(data), nrow = nspace, ncol = 2L)
      data_out <- t(zz)
    } else if (length(data) == 2L * nspace) {
      data_out <- as.numeric(data)
    } else {
      stop("Bivariate grid data are not compatible with the supplied spatial axes")
    }
  } else if (!is.null(coordt) && length(coordt) > 1L) {
    nt <- length(coordt)
    if (is.matrix(data) && nrow(data) == nt && ncol(data) == nspace) {
      data_out <- data
    } else if (is.matrix(data) && ncol(data) == nt && nrow(data) == nspace) {
      data_out <- t(data)
    } else if (!is.null(dim(data)) && tail(dim(data), 1L) == nt &&
               prod(head(dim(data), -1L)) == nspace) {
      zz <- matrix(as.numeric(data), nrow = nspace, ncol = nt)
      data_out <- t(zz)
    } else {
      stop("Space-time grid data are not compatible with the supplied spatial axes and coordt")
    }
  } else {
    if (length(data) != nspace) {
      stop("Grid data length does not match the number of spatial grid locations")
    }
    data_out <- as.numeric(data)
  }

  list(data = data_out, coordx = coords, coordy = NULL, coordz = NULL)
}

.GeoVariogram_deduplicate_pairs <- function(idx, directed = FALSE) {
  if (isTRUE(directed) || is.null(idx$rowidx) || !length(idx$rowidx)) return(idx)
  rr <- as.integer(idx$rowidx)
  cc <- as.integer(idx$colidx)
  key <- paste(pmin(rr, cc), pmax(rr, cc), sep = ":")
  keep <- !duplicated(key)
  n <- length(keep)
  for (nm in names(idx)) {
    x <- idx[[nm]]
    if (length(x) == n) idx[[nm]] <- x[keep]
  }
  idx
}

.GeoVariogram_exact_zero_tol <- function(lags) {
  z <- as.numeric(lags)
  z <- z[is.finite(z)]
  scale <- if (length(z)) max(1, max(abs(z))) else 1
  100 * .Machine$double.eps * scale
}

.GeoVariogram_st_data_vector <- function(data, coordt, ns, dynamic = FALSE) {
  if (dynamic) {
    if (!is.list(data) || length(data) != length(coordt)) {
      stop("Dynamic space-time data must be a list with one vector per time")
    }
    out <- unlist(data, use.names = FALSE)
  } else if (is.matrix(data)) {
    nt <- length(coordt)
    if (nrow(data) == nt) {
      out <- c(t(data))
    } else if (ncol(data) == nt) {
      out <- as.numeric(data)
    } else {
      stop("Fixed-site space-time data must have one row (or legacy column) per temporal coordinate")
    }
  } else {
    out <- as.numeric(data)
  }

  if (length(out) != sum(ns)) {
    stop("The space-time data length is not compatible with the spatial/temporal layout")
  }
  out
}

.GeoVariogram_biv_data_vector <- function(data, n1, n2) {
  if (is.list(data)) {
    if (length(data) != 2L || length(data[[1L]]) != n1 || length(data[[2L]]) != n2) {
      stop("Bivariate list data must contain one vector for each variable")
    }
    return(c(as.numeric(data[[1L]]), as.numeric(data[[2L]])))
  }

  if (is.matrix(data)) {
    if (nrow(data) == 2L && ncol(data) == n1 && n1 == n2) {
      return(c(as.numeric(data[1L, ]), as.numeric(data[2L, ])))
    }
    if (ncol(data) == 2L && nrow(data) == n1 && n1 == n2) {
      return(c(as.numeric(data[, 1L]), as.numeric(data[, 2L])))
    }
  }

  out <- as.numeric(data)
  if (length(out) != n1 + n2) {
    stop("The bivariate data length is not compatible with the two spatial supports")
  }
  out
}

.GeoVariogram_neighbor_coords <- function(coords, distance) {
  coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  key <- tolower(gsub("[[:blank:]]", "", as.character(distance)[1L]))
  if (key != "eucl" && ncol(coords) == 3L) coords <- coords[, 1:2, drop = FALSE]
  coords
}

.GeoVariogram_safe_max <- function(x, default = 1) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  max(x)
}

GeoVariogram <- function(data, coordx, coordy=NULL, coordz=NULL, coordt=NULL, coordx_dyn=NULL,
                         cloud=FALSE, distance="Eucl",
                         grid=FALSE, maxdist=NULL, neighb=NULL, maxtime=NULL, numbins=NULL,
                         radius=1, type='variogram', bivariate=FALSE, subsample=1, subsample_t=1,
                         numbins_t=NULL, directed=FALSE)
{
  call <- match.call()
  grid_requested <- isTRUE(grid)
  corrmodel <- 'exponential'

  if (missing(coordx)) {
    if (is.null(coordx_dyn))
      stop("argument 'coordx' is missing: provide 'coordx' or 'coordx_dyn'.")
    if (!is.list(coordx_dyn) || length(coordx_dyn) < 1L)
      stop("'coordx_dyn' must be a non-empty list when 'coordx' is omitted.")
    base <- coordx_dyn[[1L]]
    if (is.data.frame(base)) base <- as.matrix(base)
    if (!is.matrix(base))
      stop("coordx_dyn[[1]] must be a matrix/data.frame (n x d).")
    coordx <- base
  }

  if (is.data.frame(coordx)) coordx <- as.matrix(coordx)

  if (grid_requested && !is.null(coordx_dyn)) {
    stop("grid=TRUE cannot be combined with coordx_dyn")
  }
  if (grid_requested) {
    gg <- .GeoVariogram_normalize_grid(
      data = data, coordx = coordx, coordy = coordy, coordz = coordz,
      coordt = coordt, bivariate = bivariate
    )
    data <- gg$data
    coordx <- gg$coordx
    coordy <- gg$coordy
    coordz <- gg$coordz
    ## All subsequent validation/computation uses the canonical explicit
    ## coordinate representation.  grid_requested is retained in the output.
    grid <- FALSE
  }

  dynamic_st <- !isTRUE(bivariate) && !is.null(coordx_dyn) && !is.null(coordt)
  if (dynamic_st) {
    coordx_dyn <- .GeoNormalizeDynamicCoordinates(
      coordx_dyn, coordt = coordt, context = "GeoVariogram"
    )
    .GeoVariogram_validate_dynamic_data(
      data, coordx_dyn, context = "GeoVariogram"
    )
    coordx <- coordx_dyn[[1L]]
    coordy <- NULL
    coordz <- NULL
  } else {
    .GeoVariogram_validate_data(data, context = "GeoVariogram")
  }

  if (!is.numeric(subsample) || length(subsample) != 1L ||
      !is.finite(subsample) || subsample <= 0 || subsample > 1)
    stop("subsample must be a numeric value in (0, 1]")
  if (!is.numeric(subsample_t) || length(subsample_t) != 1L ||
      !is.finite(subsample_t) || subsample_t <= 0 || subsample_t > 1)
    stop("subsample_t must be a numeric value in (0, 1]")
  if (!is.logical(cloud) || length(cloud) != 1L || is.na(cloud))
    stop("cloud must be TRUE or FALSE")
  if (!is.logical(bivariate) || length(bivariate) != 1L || is.na(bivariate))
    stop("bivariate must be TRUE or FALSE")
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed))
    stop("directed must be TRUE or FALSE")

  if (!is.null(maxdist)) {
    allowed_maxdist_lengths <- if (bivariate) c(1L, 3L) else 1L
    if (!is.numeric(maxdist) || !(length(maxdist) %in% allowed_maxdist_lengths) ||
        any(!is.finite(maxdist)) || any(maxdist <= 0)) {
      stop(if (bivariate)
        "maxdist must be a positive finite numeric scalar or a length-3 vector in bivariate mode"
      else "maxdist must be a positive finite numeric scalar")
    }
    maxdist <- as.numeric(maxdist)
  }

  if (!is.null(numbins_t)) {
    if (!is.numeric(numbins_t) || length(numbins_t) != 1L ||
        !is.finite(numbins_t) || numbins_t < 1 || numbins_t %% 1 != 0) {
      stop("numbins_t must be a positive integer")
    }
    numbins_t <- as.integer(numbins_t)
  }

  if (is.null(type)) type <- "variogram"
  if (!is.character(type) || length(type) != 1L ||
      tolower(type) != "variogram")
    stop("type must be 'variogram'")
  type <- "variogram"
  model <- "Gaussian"

  if (is.null(numbins)) numbins <- 13L
  if (!is.numeric(numbins) || length(numbins) != 1L || !is.finite(numbins) ||
      numbins < 2 || numbins %% 1 != 0)
    stop("numbins must be an integer >= 2")
  numbins <- as.integer(numbins)

  if (!is.numeric(radius) || length(radius) != 1L || !is.finite(radius) || radius <= 0)
    stop("radius must be a positive finite numeric scalar")

  if (bivariate) coordt <- c(0, 1)
  if (!is.null(coordt) && is.numeric(coordt) && length(coordt) > 1L)
    corrmodel <- "gneiting"

  if (bivariate) subsample_t <- 1

  ## Temporal subsampling follows the documented T x N fixed-site ordering,
  ## while accepting the historical N x T orientation for compatibility.
  if (!bivariate && !is.null(coordt) && length(coordt) > 1L && subsample_t < 1) {
    nt0 <- length(coordt)
    mt <- min(nt0, max(2L, floor(subsample_t * nt0)))
    idx_t <- sort(sample.int(nt0, mt, replace = FALSE))

    if (is.list(data) && length(data) == nt0) {
      coordt <- coordt[idx_t]
      data <- data[idx_t]
      if (dynamic_st) {
        coordx_dyn <- coordx_dyn[idx_t]
        coordx <- coordx_dyn[[1L]]
      }
    } else if (is.matrix(data) && nrow(data) == nt0) {
      coordt <- coordt[idx_t]
      data <- data[idx_t, , drop = FALSE]
    } else if (is.matrix(data) && ncol(data) == nt0) {
      coordt <- coordt[idx_t]
      data <- data[, idx_t, drop = FALSE]
    } else if (is.array(data) && length(dim(data)) >= 3L && tail(dim(data), 1L) == nt0) {
      subs <- lapply(dim(data), seq_len)
      subs[[length(subs)]] <- idx_t
      data <- do.call(`[`, c(list(data), subs, list(drop = FALSE)))
      coordt <- coordt[idx_t]
    } else {
      warning("subsample_t was requested but the data layout is not compatible with coordt; ignoring subsample_t")
    }
  }

  .subset_data_spatial <- function(z, idx0, n0) {
    if (is.matrix(z)) {
      if (nrow(z) == n0) return(z[idx0, , drop = FALSE])
      if (ncol(z) == n0) return(z[, idx0, drop = FALSE])
      stop("Spatial subsampling cannot identify the location dimension of data")
    }
    if (is.list(z)) {
      if (!all(vapply(z, function(d) length(d) == n0, logical(1))))
        stop("Spatial subsampling requires one value per location in every data block")
      return(lapply(z, function(d) d[idx0]))
    }
    if (length(z) != n0)
      stop("Spatial subsampling requires one data value per location")
    z[idx0]
  }

  ## For a regular grid, sample actual grid locations rather than coordinate axes.
  if (subsample < 1 && isTRUE(grid) && is.null(coordx_dyn)) {
    if (is.matrix(coordx)) {
      full_coords <- as.matrix(coordx)
    } else {
      if (is.null(coordy)) stop("grid=TRUE requires coordy when coordx is not a coordinate matrix")
      full_coords <- if (is.null(coordz)) {
        as.matrix(expand.grid(coordx, coordy))
      } else {
        as.matrix(expand.grid(coordx, coordy, coordz))
      }
    }
    n0 <- nrow(full_coords)
    m <- min(n0, max(2L, floor(subsample * n0)))
    idx0 <- sort(sample.int(n0, m, replace = FALSE))

    if (bivariate) {
      data <- .subset_data_spatial(data, idx0, n0)
    } else if (!is.null(coordt) && length(coordt) > 1L) {
      nt <- length(coordt)
      if (is.array(data) && length(dim(data)) >= 3L && tail(dim(data), 1L) == nt &&
          prod(head(dim(data), -1L)) == n0) {
        z <- matrix(as.numeric(data), nrow = n0, ncol = nt)
        data <- t(z[idx0, , drop = FALSE])
      } else {
        data <- .subset_data_spatial(data, idx0, n0)
      }
    } else {
      z <- as.numeric(data)
      if (length(z) != n0) stop("Grid data length does not match the number of grid locations")
      data <- z[idx0]
    }

    coordx <- full_coords[idx0, , drop = FALSE]
    coordy <- NULL
    coordz <- NULL
    grid <- FALSE
  } else if (subsample < 1 && dynamic_st) {
    sampled <- lapply(seq_along(coordx_dyn), function(i) {
      n0 <- nrow(coordx_dyn[[i]])
      m <- min(n0, max(1L, floor(subsample * n0)))
      idx0 <- sort(sample.int(n0, m, replace = FALSE))
      list(coords = coordx_dyn[[i]][idx0, , drop = FALSE], data = data[[i]][idx0])
    })
    coordx_dyn <- lapply(sampled, `[[`, "coords")
    data <- lapply(sampled, `[[`, "data")
    coordx <- coordx_dyn[[1L]]
  } else if (subsample < 1 && bivariate && !is.null(coordx_dyn)) {
    if (!is.list(coordx_dyn) || length(coordx_dyn) != 2L ||
        nrow(coordx_dyn[[1L]]) != nrow(coordx_dyn[[2L]])) {
      stop("Bivariate spatial subsampling requires two supports with the same number of locations")
    }
    n0 <- nrow(coordx_dyn[[1L]])
    m <- min(n0, max(2L, floor(subsample * n0)))
    idx0 <- sort(sample.int(n0, m, replace = FALSE))
    coordx_dyn <- lapply(coordx_dyn, function(x) as.matrix(x)[idx0, , drop = FALSE])
    data <- .subset_data_spatial(data, idx0, n0)
    coordx <- coordx_dyn[[1L]]
    coordy <- NULL
    coordz <- NULL
  } else if (subsample < 1) {
    if (is.matrix(coordx)) {
      n0 <- nrow(coordx)
      m <- min(n0, max(2L, floor(subsample * n0)))
      idx0 <- sort(sample.int(n0, m, replace = FALSE))
      coordx <- coordx[idx0, , drop = FALSE]
      if (!is.null(coordy) && length(coordy) == n0) coordy <- coordy[idx0]
      if (!is.null(coordz) && length(coordz) == n0) coordz <- coordz[idx0]
      data <- .subset_data_spatial(data, idx0, n0)
    } else {
      n0 <- length(coordx)
      m <- min(n0, max(2L, floor(subsample * n0)))
      idx0 <- sort(sample.int(n0, m, replace = FALSE))
      coordx <- coordx[idx0]
      if (!is.null(coordy) && length(coordy) == n0) coordy <- coordy[idx0]
      if (!is.null(coordz) && length(coordz) == n0) coordz <- coordz[idx0]
      data <- .subset_data_spatial(data, idx0, n0)
    }
  }

  checkinput <- CkInput(coordx, coordy, coordz, coordt, coordx_dyn, corrmodel, data,
                        distance, "Fitting", NULL, grid, "None", maxdist, maxtime,
                        model, NULL, "Nelder-Mead", NULL, radius, NULL, NULL, NULL,
                        "GeoWLS", FALSE, FALSE, NULL, NULL)
  if (!is.null(checkinput$error)) stop(checkinput$error)

  if (bivariate) corrmodel <- "Bi_matern_sep"

  initparam <- StartParam(coordx, coordy, coordz, coordt, coordx_dyn, corrmodel, data,
                          distance, "Fitting", NULL, grid, "None", maxdist, neighb,
                          maxtime, model, 1, NULL, NULL, FALSE, radius, NULL, NULL,
                          NULL, "GeoWLS", "GeoWLS", FALSE, NULL, NULL, FALSE, FALSE)
  if (!is.null(initparam$error)) stop(initparam$error)

  ## GeoWLS setup always uses SetGlobalVar(), including neighborhood calls.
  ## Register the matching cleanup once so errors after setup cannot leak C state.
  on.exit({
    try(dotCall64::.C64("DeleteGlobalVar", SIGNATURE = character(0),
                        INTENT = character(0), PACKAGE = "GeoModels",
                        VERBOSE = 0, NAOK = TRUE), silent = TRUE)
  }, add = TRUE)

  coordx <- initparam$coordx
  coordy <- initparam$coordy
  coordz <- initparam$coordz
  coordt <- initparam$coordt
  if (is.null(coordz)) coordz <- double(length(coordx))

  spacetime_dyn <- if (isTRUE(initparam$spacetime)) !is.null(coordx_dyn) else FALSE
  bins <- bint <- centers <- centert <- NULL
  lenbins <- lenbinst <- lenbint <- NULL
  variograms <- variogramst <- variogramt <- NULL
  regular <- NULL
  time_breaks <- NULL

  if (cloud && (isTRUE(initparam$spacetime) || isTRUE(initparam$bivariate))) {
    stop("cloud=TRUE is supported for univariate spatial variograms only")
  }

  ##########################################################################
  ## Bivariate
  ##########################################################################
  if (isTRUE(initparam$bivariate)) {
    n_var <- initparam$numtime
    if (n_var != 2L) stop("GeoVariogram currently supports exactly two variables in bivariate mode")
    ns <- as.integer(initparam$ns)
    if (length(ns) < 2L) stop("Invalid bivariate spatial layout")
    n1 <- ns[1L]
    n2 <- ns[2L]

    if (!is.null(coordx_dyn)) {
      if (length(coordx_dyn) != 2L) stop("Bivariate coordx_dyn must contain two spatial supports")
      X1 <- as.matrix(coordx_dyn[[1L]])
      X2 <- as.matrix(coordx_dyn[[2L]])
      scale0 <- suppressWarnings(max(abs(c(X1, X2)), na.rm = TRUE))
      if (!is.finite(scale0) || scale0 < 1) scale0 <- 1
      tol <- 100 * .Machine$double.eps * scale0
      aligned <- identical(dim(X1), dim(X2)) && all(abs(X1 - X2) <= tol)
      if (!aligned) {
        stop("A classical bivariate cross-semivariogram requires aligned spatial supports; the two coordx_dyn grids are not aligned")
      }
    }

    data_biv <- .GeoVariogram_biv_data_vector(data, n1, n2)
    numvario <- numbins - 1L

    if (is.null(neighb)) {
      NS <- c(0L, n1)
      moments_marg <- double(2L * numvario)
      lenbins_marg <- integer(2L * numvario)
      moments_cross <- double(numvario)
      lenbins_cross <- integer(numvario)

      EV <- dotCall64::.C64(
        "Binned_Variogram_biv2",
        bins = double(numbins), coordx = as.double(coordx), coordy = as.double(coordy),
        coordz = as.double(coordz), coordt = as.double(coordt), data = as.double(data_biv),
        lenbins_cross = lenbins_cross, moments_cross = moments_cross,
        numbins = as.integer(numbins), lenbins_marg = lenbins_marg,
        moments_marg = moments_marg, ns = ns, NS = as.integer(NS),
        SIGNATURE = c("double","double","double","double","double","double",
                      "integer","double","integer","integer","double","integer","integer"),
        INTENT = c("rw","r","r","r","r","r","rw","rw","r","rw","rw","r","r"),
        NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0
      )

      bins <- EV$bins
      l11 <- EV$lenbins_marg[seq_len(numvario)]
      l22 <- EV$lenbins_marg[numvario + seq_len(numvario)]
      l12 <- EV$lenbins_cross[seq_len(numvario)]
      m11 <- EV$moments_marg[seq_len(numvario)]
      m22 <- EV$moments_marg[numvario + seq_len(numvario)]
      m12 <- EV$moments_cross[seq_len(numvario)]
    } else {
      if (is.null(coordx_dyn)) {
        base <- cbind(coordx[seq_len(n1)], coordy[seq_len(n1)], coordz[seq_len(n1)])
        pair_coords <- rbind(base, base)
        pair_coords <- .GeoVariogram_neighbor_coords(pair_coords, distance)
        idx <- GeoNeighIndex(coordx = pair_coords, distance = distance,
                             neighb = neighb, maxdist = maxdist, radius = radius,
                             bivariate = TRUE)
      } else {
        idx <- GeoNeighIndex(coordx_dyn = coordx_dyn, distance = distance,
                             neighb = neighb, maxdist = maxdist, radius = radius,
                             bivariate = TRUE)
      }

      breaks <- .GeoVariogram_make_breaks(idx$lags, numbins)
      if (!length(breaks)) {
        bins <- numeric(0L)
        centers <- numeric(0L)
        lenbins <- matrix(integer(0L), nrow = 2L)
        lenbinst <- integer(0L)
        variograms <- matrix(numeric(0L), nrow = 2L)
        variogramst <- numeric(0L)
      } else {
        bins <- breaks
        row <- as.integer(idx$rowidx)
        col <- as.integer(idx$colidx)
        first <- as.integer(idx$first)
        second <- as.integer(idx$second)
        lagv <- as.numeric(idx$lags)

        z1 <- data_biv[seq_len(n1)]
        z2 <- data_biv[n1 + seq_len(n2)]

        sel11 <- which(first == 0L & second == 0L)
        sel22 <- which(first == 1L & second == 1L)
        sel12 <- which(first != second)

        unique_within <- function(ii) {
          if (isTRUE(directed) || !length(ii)) return(ii)
          key <- paste(pmin(row[ii], col[ii]), pmax(row[ii], col[ii]), sep = ":")
          ii[!duplicated(key)]
        }
        sel11 <- unique_within(sel11)
        sel22 <- unique_within(sel22)

        b11 <- .GeoVariogram_bin_values(
          lagv[sel11], 0.5 * (data_biv[row[sel11]] - data_biv[col[sel11]])^2, breaks)
        b22 <- .GeoVariogram_bin_values(
          lagv[sel22], 0.5 * (data_biv[row[sel22]] - data_biv[col[sel22]])^2, breaks)

        cross_lag <- numeric(0L)
        cross_val <- numeric(0L)
        if (length(sel12)) {
          rr <- row[sel12]
          cc <- col[sel12]
          ff <- first[sel12]
          ii <- ifelse(ff == 0L, rr, cc)
          jj <- ifelse(ff == 0L, cc - n1, rr - n1)
          zero_tol <- .GeoVariogram_exact_zero_tol(lagv[sel12])
          ok <- n1 == n2 & ii >= 1L & ii <= n1 & jj >= 1L & jj <= n2 &
            ii != jj & is.finite(lagv[sel12]) & lagv[sel12] > zero_tol
          if (any(ok)) {
            ii <- ii[ok]
            jj <- jj[ok]
            kk <- sel12[ok]
            if (!isTRUE(directed)) {
              key <- paste(pmin(ii, jj), pmax(ii, jj), sep = ":")
              keep <- !duplicated(key)
              ii <- ii[keep]
              jj <- jj[keep]
              kk <- kk[keep]
            }
            cross_lag <- lagv[kk]
            cross_val <- 0.5 * (z1[ii] - z1[jj]) * (z2[ii] - z2[jj])
          }
        }
        b12 <- .GeoVariogram_bin_values(cross_lag, cross_val, breaks)

        l11 <- b11$count; l22 <- b22$count; l12 <- b12$count
        m11 <- b11$sum; m22 <- b22$sum; m12 <- b12$sum
      }
    }

    if (is.null(centers)) {
      centers_all <- bins[-length(bins)] + diff(bins) / 2
      keep <- (l11 + l22 + l12) > 0L
      centers <- centers_all[keep]
      l11 <- l11[keep]; l22 <- l22[keep]; l12 <- l12[keep]
      m11 <- m11[keep]; m22 <- m22[keep]; m12 <- m12[keep]

      v11 <- rep(NA_real_, length(l11)); v22 <- v11; v12 <- v11
      ii <- l11 > 0L; v11[ii] <- m11[ii] / l11[ii]
      ii <- l22 > 0L; v22[ii] <- m22[ii] / l22[ii]
      ii <- l12 > 0L; v12[ii] <- m12[ii] / l12[ii]
      lenbins <- rbind(l11, l22)
      lenbinst <- l12
      variograms <- rbind(v11, v22)
      variogramst <- v12
    }
  }

  ##########################################################################
  ## Space-time
  ##########################################################################
  if (isTRUE(initparam$spacetime)) {
    ns <- as.integer(initparam$ns)
    NS <- as.integer(c(0L, cumsum(ns))[-(length(ns) + 1L)])
    data_st <- .GeoVariogram_st_data_vector(data, coordt, ns, dynamic = spacetime_dyn)
    ti <- .GeoVariogram_time_bins(coordt, maxtime,
                                    numbins_t = if (is.null(numbins_t)) numbins - 1L else numbins_t)
    time_breaks <- ti$breaks
    centert_all <- ti$centers
    regular <- ti$regular
    ntbin <- length(time_breaks) - 1L
    nhbin <- numbins - 1L

    if (is.null(neighb)) {
      EV <- dotCall64::.C64(
        if (spacetime_dyn) "Binned_Variogram_st2_dyn" else "Binned_Variogram_st2",
        bins = double(numbins), bint = as.double(time_breaks),
        coordx = as.double(coordx), coordy = as.double(coordy), coordz = as.double(coordz),
        coordt = as.double(coordt), data = as.double(data_st),
        lenbins = integer(nhbin), lenbinst = integer(nhbin * ntbin),
        lenbint = integer(ntbin), moments = double(nhbin),
        momentst = double(nhbin * ntbin), momentt = double(ntbin),
        numbins = as.integer(numbins), numbint = as.integer(length(time_breaks)),
        ns = ns, NS = NS,
        SIGNATURE = c(rep("double", 7), rep("integer", 3), rep("double", 3), rep("integer", 4)),
        INTENT = c("rw","rw",rep("r",5),rep("rw",3),rep("rw",3),rep("r",4)),
        NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0
      )

      bins <- EV$bins
      spatial_count <- EV$lenbins
      spatial_sum <- EV$moments
      temporal_count <- EV$lenbint
      temporal_sum <- EV$momentt
      st_count <- matrix(EV$lenbinst, nrow = nhbin, ncol = ntbin, byrow = TRUE)
      st_sum <- matrix(EV$momentst, nrow = nhbin, ncol = ntbin, byrow = TRUE)
    } else {
      if (spacetime_dyn) {
        idx <- GeoNeighIndex(coordx_dyn = coordx_dyn, coordt = coordt,
                             distance = distance, neighb = neighb,
                             maxdist = maxdist, maxtime = ti$maxlag, radius = radius)
      } else {
        nsite <- ns[1L]
        base <- cbind(coordx[seq_len(nsite)], coordy[seq_len(nsite)], coordz[seq_len(nsite)])
        base <- .GeoVariogram_neighbor_coords(base, distance)
        idx <- GeoNeighIndex(coordx = base, coordt = coordt,
                             distance = distance, neighb = neighb,
                             maxdist = maxdist, maxtime = ti$maxlag, radius = radius)
      }

      idx <- .GeoVariogram_deduplicate_pairs(idx, directed = directed)
      lagv <- as.numeric(idx$lags)
      lagt <- as.numeric(idx$lagt)
      rr <- as.integer(idx$rowidx)
      cc <- as.integer(idx$colidx)
      ok_idx <- rr >= 1L & rr <= length(data_st) & cc >= 1L & cc <= length(data_st)
      lagv <- lagv[ok_idx]; lagt <- lagt[ok_idx]
      rr <- rr[ok_idx]; cc <- cc[ok_idx]
      values <- 0.5 * (data_st[rr] - data_st[cc])^2

      positive_spatial <- lagv[is.finite(lagv) & lagv > 0]
      bins <- .GeoVariogram_make_breaks(positive_spatial, numbins)
      if (!length(bins)) {
        bins <- seq(0, 1, length.out = numbins)
      }
      nhbin <- length(bins) - 1L

      sel_s <- is.finite(lagt) & abs(lagt) <= sqrt(.Machine$double.eps) & lagv > 0
      zero_tol <- .GeoVariogram_exact_zero_tol(lagv)
      sel_t <- is.finite(lagv) & abs(lagv) <= zero_tol & lagt > 0
      ## The temporal margin is gamma(0,u): collocated pairs are not mixed
      ## into the positive-distance space-time surface.
      sel_st <- lagv > zero_tol & lagt > 0

      bs <- .GeoVariogram_bin_values(lagv[sel_s], values[sel_s], bins)
      bt <- .GeoVariogram_bin_values(lagt[sel_t], values[sel_t], time_breaks)
      bst <- .GeoVariogram_bin_surface(lagv[sel_st], lagt[sel_st], values[sel_st],
                                       bins, time_breaks)
      spatial_count <- bs$count; spatial_sum <- bs$sum
      temporal_count <- bt$count; temporal_sum <- bt$sum
      st_count <- bst$count; st_sum <- bst$sum
    }

    centers_all <- bins[-length(bins)] + diff(bins) / 2
    spatial_mean <- rep(NA_real_, length(spatial_count))
    ii <- spatial_count > 0L
    spatial_mean[ii] <- spatial_sum[ii] / spatial_count[ii]
    temporal_mean <- rep(NA_real_, length(temporal_count))
    ii <- temporal_count > 0L
    temporal_mean[ii] <- temporal_sum[ii] / temporal_count[ii]
    st_mean <- matrix(NA_real_, nrow = nrow(st_count), ncol = ncol(st_count))
    ii <- st_count > 0L
    st_mean[ii] <- st_sum[ii] / st_count[ii]

    keep_s <- spatial_count > 0L | rowSums(st_count) > 0L
    keep_t <- temporal_count > 0L | colSums(st_count) > 0L
    centers <- centers_all[keep_s]
    centert <- centert_all[keep_t]
    spatial_count <- spatial_count[keep_s]
    spatial_mean <- spatial_mean[keep_s]
    temporal_count <- temporal_count[keep_t]
    temporal_mean <- temporal_mean[keep_t]
    st_count <- st_count[keep_s, keep_t, drop = FALSE]
    st_mean <- st_mean[keep_s, keep_t, drop = FALSE]

    lenbins <- spatial_count
    variograms <- spatial_mean
    lenbint <- temporal_count
    variogramt <- temporal_mean
    lenbinst <- as.integer(as.vector(t(st_count)))
    variogramst <- as.numeric(as.vector(t(st_mean)))
    bint <- centert
  }

  ##########################################################################
  ## Purely spatial
  ##########################################################################
  if (!isTRUE(initparam$bivariate) && !isTRUE(initparam$spacetime)) {
    nsite <- length(coordx)
    coords <- cbind(coordx, coordy, coordz)
    data_sp <- as.numeric(data)
    if (length(data_sp) != nsite) stop("The number of spatial data values differs from the number of locations")

    if (cloud) {
      if (!is.null(neighb)) {
        nncoords <- .GeoVariogram_neighbor_coords(coords, distance)
        idx <- GeoNeighIndex(coordx = nncoords, distance = distance, neighb = neighb,
                             maxdist = maxdist, radius = radius)
        idx <- .GeoVariogram_deduplicate_pairs(idx, directed = directed)
        rr <- as.integer(idx$rowidx); cc <- as.integer(idx$colidx)
        lagv <- as.numeric(idx$lags)
        val <- 0.5 * (data_sp[rr] - data_sp[cc])^2
        ok <- is.finite(lagv) & is.finite(val)
        bins <- centers <- lagv[ok]
        variograms <- val[ok]
        lenbins <- rep.int(1L, length(variograms))
      } else if (!is.null(maxdist) && length(maxdist) == 1L && is.finite(maxdist)) {
        pp <- .GeoScatterplot_spatial_pairs_maxdist(
          data_sp, coords, maxdist = maxdist, bins = c(0, maxdist),
          distance = distance, radius = radius
        )
        bins <- centers <- pp$lags
        variograms <- 0.5 * (pp$v1 - pp$v2)^2
        lenbins <- rep.int(1L, length(variograms))
      } else {
        nmax <- as.double(nsite) * (nsite - 1) / 2
        if (!is.finite(nmax) || nmax > .Machine$integer.max) {
          stop("The requested variogram cloud is too large for the native interface; use maxdist, neighb, or subsample")
        }
        nmax <- as.integer(nmax)
        EV <- dotCall64::.C64(
          "Cloud_Variogram2",
          bins = double(nmax), coordx = as.double(coordx), coordy = as.double(coordy),
          coordz = as.double(coordz), coordt = 0, data = data_sp,
          lenbins = integer(nmax), moments = double(nmax), numbins = nmax,
          SIGNATURE = c("double","double","double","double","double","double","integer","double","integer"),
          INTENT = c("rw","r","r","r","r","r","rw","rw","r"),
          NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0
        )
        keep <- EV$lenbins > 0L
        bins <- centers <- EV$bins[keep]
        lenbins <- EV$lenbins[keep]
        variograms <- EV$moments[keep]
      }
    } else if (is.null(neighb)) {
      EV <- dotCall64::.C64(
        "Binned_Variogram2",
        bins = double(numbins), coordx = as.double(coordx), coordy = as.double(coordy),
        coordz = as.double(coordz), coordt = 0, data = data_sp,
        lenbins = integer(numbins - 1L), moments = double(numbins - 1L),
        numbins = as.integer(numbins),
        SIGNATURE = c("double","double","double","double","double","double","integer","double","integer"),
        INTENT = c("rw","r","r","r","r","r","rw","rw","r"),
        NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0
      )
      bins <- EV$bins
      centers_all <- bins[-length(bins)] + diff(bins) / 2
      keep <- EV$lenbins > 0L
      centers <- centers_all[keep]
      lenbins <- EV$lenbins[keep]
      variograms <- EV$moments[keep] / lenbins
    } else {
      nncoords <- .GeoVariogram_neighbor_coords(coords, distance)
      idx <- GeoNeighIndex(coordx = nncoords, distance = distance, neighb = neighb,
                           maxdist = maxdist, radius = radius)
      idx <- .GeoVariogram_deduplicate_pairs(idx, directed = directed)
      lagv <- as.numeric(idx$lags)
      rr <- as.integer(idx$rowidx); cc <- as.integer(idx$colidx)
      breaks <- .GeoVariogram_make_breaks(lagv, numbins)
      if (!length(breaks)) {
        bins <- centers <- variograms <- numeric(0L)
        lenbins <- integer(0L)
      } else {
        bv <- .GeoVariogram_bin_values(
          lagv, 0.5 * (data_sp[rr] - data_sp[cc])^2, breaks
        )
        bins <- breaks
        centers_all <- bins[-length(bins)] + diff(bins) / 2
        keep <- bv$count > 0L
        centers <- centers_all[keep]
        lenbins <- bv$count[keep]
        variograms <- bv$mean[keep]
      }
    }
  }

  ans <- list(
    bins = bins,
    bint = bint,
    bivariate = isTRUE(initparam$bivariate),
    cloud = cloud,
    centers = centers,
    centert = centert,
    distance = distance,
    radius = radius,
    grid = grid_requested,
    neighb = neighb,
    directed = directed,
    lenbins = lenbins,
    lenbinst = lenbinst,
    lenbint = lenbint,
    maxdist = maxdist,
    maxtime = maxtime,
    numbins_t = numbins_t,
    regular = regular,
    time.breaks = time_breaks,
    spacetime_dyn = spacetime_dyn,
    temporal.margin = if (isTRUE(spacetime_dyn)) "collocated" else if (isTRUE(initparam$spacetime)) "same-site" else NULL,
    subsample = subsample,
    subsample_t = subsample_t,
    variograms = variograms,
    variogramst = variogramst,
    variogramt = variogramt,
    type = type
  )

  structure(c(ans, call = call), class = "GeoVariogram")
}

####################################################
### Plot method
####################################################
plot.GeoVariogram <- function(x, ...) {
  if (!inherits(x, "GeoVariogram"))
    stop("Enter an object obtained from the function GeoVariogram")

  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  dots <- list(...)
  xlab <- if (!is.null(dots$xlab)) dots$xlab else "Distance"
  ylab <- if (!is.null(dots$ylab)) dots$ylab else "Semi-Variogram"
  dots$xlab <- NULL
  dots$ylab <- NULL

  plot_with_defaults <- function(defaults) {
    do.call(plot.default, utils::modifyList(defaults, dots))
  }

  safe_ylim <- function(z, include_zero = TRUE, symmetric = FALSE) {
    z <- as.numeric(z)
    z <- z[is.finite(z)]
    if (!length(z)) return(c(0, 1))
    if (symmetric) {
      a <- max(abs(z))
      if (!is.finite(a) || a <= 0) a <- 1
      return(c(-a, a))
    }
    lo <- if (include_zero) min(0, z) else min(z)
    hi <- max(z)
    if (!is.finite(hi) || hi <= lo) hi <- lo + 1
    c(lo, hi)
  }

  ispatim <- !is.null(x$centert) || !is.null(x$bint)
  bivariate <- isTRUE(x$bivariate)

  if (bivariate) {
    par(mfrow = c(2, 2))
    xmax <- .GeoVariogram_safe_max(x$centers)

    plot_with_defaults(list(
      x = x$centers, y = x$variograms[1L, ], main = "First semi-variogram",
      ylim = safe_ylim(x$variograms[1L, ]), xlim = c(0, xmax),
      xlab = xlab, ylab = ylab))

    plot_with_defaults(list(
      x = x$centers, y = x$variogramst, main = "Cross semi-variogram",
      ylim = safe_ylim(x$variogramst, symmetric = TRUE), xlim = c(0, xmax),
      xlab = xlab, ylab = ylab))

    plot_with_defaults(list(
      x = x$centers, y = x$variograms[2L, ], main = "Second semi-variogram",
      ylim = safe_ylim(x$variograms[2L, ]), xlim = c(0, xmax),
      xlab = xlab, ylab = ylab))

    plot.new()
    return(invisible())
  }

  if (ispatim) {
    par(mfrow = c(2, 2))
    temporal_lags <- if (!is.null(x$centert)) x$centert else x$bint
    nx <- length(x$centers)
    nt <- length(temporal_lags)
    if (length(x$variogramst) != nx * nt)
      stop("The space-time variogram surface is not rectangular")

    cross <- matrix(x$variogramst, nrow = nx, ncol = nt, byrow = TRUE)
    evario <- rbind(c(0, x$variogramt), cbind(x$variograms, cross))
    gx <- c(0, x$centers)
    gt <- c(0, temporal_lags)
    egrid <- expand.grid(gx, gt)
    zz <- as.numeric(evario)
    ok <- is.finite(zz)

    if (any(ok) && requireNamespace("scatterplot3d", quietly = TRUE)) {
      scatterplot3d::scatterplot3d(
        egrid[ok, 1L], egrid[ok, 2L], zz[ok], type = "h", highlight.3d = TRUE,
        cex.axis = .7, cex.lab = .7, main = "Empirical space-time semi-variogram",
        xlab = "Distance", ylab = "Time", zlab = ylab,
        mar = c(2, 2, 2, 2), mgp = c(0, 0, 0)
      )
    } else {
      plot.new(); title("Empirical space-time semi-variogram")
      if (any(ok)) {
        graphics::text(0.5, 0.5, "Install 'scatterplot3d' for the 3D empirical panel")
      }
    }

    par(mai = c(.2, .2, .2, .2))
    persp(x = gx, y = gt, z = evario,
          xlab = "h", ylab = "u", zlab = expression(gamma(h, u)),
          ltheta = 90, shade = 0.75, ticktype = "detailed", phi = 30,
          theta = 30, main = "Smoothed space-time semi-variogram",
          cex.axis = .8, cex.lab = .8)

    par(mai = c(.5, .5, .5, .5), mgp = c(1.6, .6, 0))
    if (any(is.finite(x$variogramt))) {
      plot_with_defaults(list(
        x = temporal_lags, y = x$variogramt,
        xlab = "t", ylab = expression(gamma(t)),
        ylim = safe_ylim(x$variogramt),
        xlim = c(0, .GeoVariogram_safe_max(temporal_lags)),
        main = "Marginal temporal semi-variogram"))
    } else {
      plot.new()
      title(main = "Marginal temporal semi-variogram")
      if (isTRUE(x$spacetime_dyn)) {
        text(0.5, 0.56, "Temporal marginal unavailable", font = 2, cex = 0.9)
        text(0.5, 0.43, "(no repeated spatial locations)", cex = 0.72)
      } else {
        graphics::text(0.5, 0.5, "No temporal marginal pairs available", cex = 0.85)
      }
      box()
    }

    plot_with_defaults(list(
      x = x$centers, y = x$variograms,
      xlab = "h", ylab = expression(gamma(h)),
      ylim = safe_ylim(x$variograms),
      xlim = c(0, .GeoVariogram_safe_max(x$centers)),
      main = "Marginal spatial semi-variogram"))
    return(invisible())
  }

  plot_with_defaults(list(x = x$centers, y = x$variograms,
                          xlab = xlab, ylab = ylab))
  invisible()
}
