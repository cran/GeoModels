####################################################
### File name: GeoNeighborhood.R
####################################################

GeoNeighborhood <- function(data = NULL, coordx = NULL, coordy = NULL,
                            coordz = NULL, coordt = NULL,
                            coordx_dyn = NULL, bivariate = FALSE,
                            distance = "Eucl", grid = FALSE, loc,
                            neighb = NULL, maxdist = NULL,
                            maxtime = NULL, radius = 1, time = NULL,
                            X = NULL, M = NULL, spobj = NULL,
                            spdata = NULL, parallel = FALSE,
                            ncores = 6, check.duplicates = FALSE) {

  ##########################################################################
  ## Small validation helpers
  ##########################################################################

  .is_scalar_logical <- function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
  }

  .check_nonnegative_scalar <- function(x, name, allow_null = TRUE) {
    if (is.null(x) && allow_null) return(NULL)
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
      stop(name, " must be a finite non-negative scalar.", call. = FALSE)
    }
    as.numeric(x)
  }

  .check_numeric_matrix <- function(x, name, ncol_allowed = c(2L, 3L)) {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (!is.matrix(x) || !is.numeric(x)) {
      stop(name, " must be a numeric matrix.", call. = FALSE)
    }
    if (nrow(x) < 1L || !(ncol(x) %in% ncol_allowed)) {
      stop(name, " must have at least one row and two or three columns.",
           call. = FALSE)
    }
    if (any(!is.finite(x))) {
      stop(name, " contains non-finite values.", call. = FALSE)
    }
    unname(x)
  }

  .check_spherical_coords <- function(x, name) {
    if (ncol(x) != 2L) {
      stop(name, " must have exactly two columns (longitude, latitude) ",
           "for Geod or Chor distance.", call. = FALSE)
    }
    if (any(x[, 2L] < -90 | x[, 2L] > 90)) {
      stop(name, " contains latitudes outside [-90, 90].", call. = FALSE)
    }
    invisible(TRUE)
  }

  .normalise_neighb <- function(x) {
    if (is.null(x)) return(NULL)
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1) {
      stop("neighb must be a finite positive integer.", call. = FALSE)
    }
    xr <- round(x)
    if (abs(x - xr) > sqrt(.Machine$double.eps)) {
      stop("neighb must be an integer.", call. = FALSE)
    }
    as.integer(xr)
  }

  .build_coords <- function(coordx, coordy, coordz, grid) {
    if (is.null(coordx)) return(NULL)

    if ((is.matrix(coordx) || is.data.frame(coordx)) &&
        is.null(coordy) && is.null(coordz)) {
      return(.check_numeric_matrix(as.matrix(coordx), "coordx"))
    }

    if (is.null(coordy)) {
      stop("coordy is required when coordx is not a coordinate matrix.",
           call. = FALSE)
    }

    if (!is.numeric(coordx) || !is.numeric(coordy) ||
        any(!is.finite(coordx)) || any(!is.finite(coordy))) {
      stop("coordx and coordy must be finite numeric vectors.", call. = FALSE)
    }

    if (!is.null(coordz) &&
        (!is.numeric(coordz) || any(!is.finite(coordz)))) {
      stop("coordz must be a finite numeric vector.", call. = FALSE)
    }

    if (isTRUE(grid)) {
      out <- if (is.null(coordz)) {
        as.matrix(expand.grid(coordx, coordy))
      } else {
        as.matrix(expand.grid(coordx, coordy, coordz))
      }
    } else {
      if (length(coordx) != length(coordy)) {
        stop("coordx and coordy must have the same length.", call. = FALSE)
      }
      if (!is.null(coordz) && length(coordz) != length(coordx)) {
        stop("coordz must have the same length as coordx and coordy.",
             call. = FALSE)
      }
      out <- if (is.null(coordz)) {
        cbind(coordx, coordy)
      } else {
        cbind(coordx, coordy, coordz)
      }
    }

    .check_numeric_matrix(out, "coordinates")
  }

  .map <- function(index, fun) {
    if (!isTRUE(parallel) || length(index) < 2L) {
      return(lapply(index, fun))
    }

    workers <- .GeoResolveWorkers(TRUE, ncores, n_jobs = length(index))
    if (workers <= 1L) return(lapply(index, fun))

    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop("parallel=TRUE requires the future and future.apply packages.",
           call. = FALSE)
    }

    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = workers)

    future.apply::future_lapply(
      index,
      fun,
      future.seed = TRUE,
      future.globals = TRUE,
      future.packages = "nabor"
    )
  }

  ##########################################################################
  ## Exact spatial distances used to verify/order the nabor candidates
  ##########################################################################

  .central_angle <- function(ref, query) {
    deg_to_rad <- pi / 180
    phi_q <- query[2L] * deg_to_rad
    phi_r <- ref[, 2L] * deg_to_rad
    delta_lon <- (ref[, 1L] - query[1L]) * deg_to_rad

    sin_q <- sin(phi_q)
    cos_q <- cos(phi_q)
    sin_r <- sin(phi_r)
    cos_r <- cos(phi_r)

    x <- cos_r * sin(delta_lon)
    y <- cos_q * sin_r - sin_q * cos_r * cos(delta_lon)
    z <- sin_q * sin_r + cos_q * cos_r * cos(delta_lon)

    atan2(sqrt(x * x + y * y), z)
  }

  .exact_distance <- function(ref, query, distance, radius) {
    if (!nrow(ref)) return(numeric(0L))

    if (identical(distance, "Eucl")) {
      delta <- sweep(ref, 2L, query, FUN = "-")
      return(sqrt(rowSums(delta * delta)))
    }

    angle <- .central_angle(ref, query)
    if (identical(distance, "Geod")) {
      radius * angle
    } else {
      2 * radius * sin(angle / 2)
    }
  }

  .to_xyz <- function(x, radius) {
    lat <- x[, 2L] * pi / 180
    lon <- x[, 1L] * pi / 180
    radius * cbind(
      cos(lat) * cos(lon),
      cos(lat) * sin(lon),
      sin(lat)
    )
  }

  .search_radius <- function(maxdist, distance, radius) {
    if (is.null(maxdist)) return(NULL)
    if (maxdist == 0) return(NULL)

    if (identical(distance, "Geod")) {
      if (maxdist >= pi * radius) return(NULL)
      return(2 * radius * sin(maxdist / (2 * radius)))
    }

    if (identical(distance, "Chor")) {
      if (maxdist >= 2 * radius) return(NULL)
      return(maxdist)
    }

    maxdist
  }

  .find_neighbors <- function(ref, query, distance, maxdist, neighb, radius) {
    ref <- .check_numeric_matrix(ref, "reference coordinates")
    query <- .check_numeric_matrix(query, "prediction locations")

    if (ncol(ref) != ncol(query)) {
      stop("Reference coordinates and prediction locations must have ",
           "the same number of columns.", call. = FALSE)
    }

    if (distance %in% c("Geod", "Chor")) {
      .check_spherical_coords(ref, "reference coordinates")
      .check_spherical_coords(query, "prediction locations")
      ref_search <- .to_xyz(ref, radius)
      query_search <- .to_xyz(query, radius)
    } else {
      ref_search <- ref
      query_search <- query
    }

    nref <- nrow(ref)
    k <- if (is.null(neighb)) nref else min(neighb, nref)
    if (k < 1L) stop("No reference coordinates are available.", call. = FALSE)

    search_radius <- .search_radius(maxdist, distance, radius)
    cutoff_tol <- if (is.null(maxdist)) {
      0
    } else {
      64 * .Machine$double.eps * max(1, abs(maxdist))
    }

    finish_one <- function(idx, dist_hint, qrow) {
      ok <- !is.na(idx) & idx > 0L & idx <= nref & is.finite(dist_hint)
      idx <- as.integer(idx[ok])
      if (!length(idx)) {
        return(list(indices = integer(0L), distances = numeric(0L)))
      }

      idx <- unique(idx)
      exact <- .exact_distance(ref[idx, , drop = FALSE], qrow,
                               distance, radius)

      if (!is.null(maxdist)) {
        keep <- exact <= maxdist + cutoff_tol
        idx <- idx[keep]
        exact <- exact[keep]
      }

      if (!length(idx)) {
        return(list(indices = integer(0L), distances = numeric(0L)))
      }

      ord <- order(exact, idx)
      idx <- idx[ord]
      exact <- exact[ord]

      if (!is.null(neighb) && length(idx) > neighb) {
        idx <- idx[seq_len(neighb)]
        exact <- exact[seq_len(neighb)]
      }

      list(indices = idx, distances = exact)
    }

    ## Radius-only neighborhoods may contain all reference points. Processing
    ## one query at a time avoids allocating a nquery x nref result matrix.
    radius_only <- is.null(neighb) && !is.null(maxdist)

    if (radius_only) {
      return(.map(seq_len(nrow(query_search)), function(i) {
        args <- list(
          data = ref_search,
          query = matrix(query_search[i, ], nrow = 1L),
          k = k
        )
        if (!is.null(search_radius)) args$radius <- search_radius
        ans <- do.call(nabor::knn, args)
        finish_one(ans$nn.idx[1L, ], ans$nn.dists[1L, ],
                   query[i, ])
      }))
    }

    args <- list(data = ref_search, query = query_search, k = k)
    if (!is.null(search_radius)) args$radius <- search_radius
    ans <- do.call(nabor::knn, args)

    lapply(seq_len(nrow(query)), function(i) {
      finish_one(ans$nn.idx[i, ], ans$nn.dists[i, ], query[i, ])
    })
  }

  .select_times <- function(target_time, coordt, maxtime) {
    delta <- abs(coordt - target_time)
    idx <- if (is.null(maxtime)) {
      seq_along(coordt)
    } else {
      which(delta <= maxtime +
              64 * .Machine$double.eps * max(1, abs(maxtime)))
    }
    if (!length(idx)) return(integer(0L))
    idx[order(coordt[idx], idx)]
  }

  ##########################################################################
  ## Handle sp/spacetime objects before determining the data structure
  ##########################################################################

  if (!is.null(spobj)) {
    spatial_info <- sp2Geo(spobj, spdata)
    coordx <- spatial_info$coords
    if (!is.null(spatial_info$coordt)) coordt <- spatial_info$coordt
    if (!isTRUE(spatial_info$pj) && !identical(distance, "Chor")) {
      distance <- "Geod"
    }
    if (!is.null(spatial_info$Y)) data <- spatial_info$Y
    if (!is.null(spatial_info$X)) X <- spatial_info$X
  }

  ##########################################################################
  ## Main input validation
  ##########################################################################

  if (missing(loc)) stop("Parameter 'loc' is required.", call. = FALSE)
  if (!.is_scalar_logical(bivariate)) {
    stop("bivariate must be TRUE or FALSE.", call. = FALSE)
  }
  if (!.is_scalar_logical(grid)) {
    stop("grid must be TRUE or FALSE.", call. = FALSE)
  }
  if (!.is_scalar_logical(parallel)) {
    stop("parallel must be TRUE or FALSE.", call. = FALSE)
  }

  distance <- match.arg(distance, c("Eucl", "Chor", "Geod"))
  neighb <- .normalise_neighb(neighb)
  maxdist <- .check_nonnegative_scalar(maxdist, "maxdist")
  maxtime <- .check_nonnegative_scalar(maxtime, "maxtime")

  if (!is.numeric(radius) || length(radius) != 1L ||
      !is.finite(radius) || radius <= 0) {
    stop("radius must be a finite positive scalar.", call. = FALSE)
  }
  radius <- as.numeric(radius)

  if (is.vector(loc) && !is.list(loc)) loc <- matrix(loc, nrow = 1L)
  loc <- .check_numeric_matrix(as.matrix(loc), "loc")

  coords <- .build_coords(coordx, coordy, coordz, grid)

  .GeoValidate3D(
    coordx = coords, coordx_dyn = coordx_dyn,
    distance = distance, grid = grid,
    context = "GeoNeighborhood"
  )
  .GeoValidate3D(
    coordx = loc, distance = distance, grid = FALSE,
    context = "GeoNeighborhood prediction locations"
  )
  spacetime <- !is.null(coordt)
  is_dynamic <- !is.null(coordx_dyn)

  ## A two-element dynamic coordinate list without temporal coordinates is the
  ## standard representation of a bivariate field on different locations.
  if (!isTRUE(bivariate) && !spacetime && is_dynamic &&
      is.list(coordx_dyn) && length(coordx_dyn) == 2L) {
    bivariate <- TRUE
  }

  if (spacetime && isTRUE(bivariate)) {
    stop("Simultaneous bivariate and spatiotemporal neighborhood selection ",
         "is not implemented.", call. = FALSE)
  }

  if (!spacetime && !isTRUE(bivariate) &&
      is.null(neighb) && is.null(maxdist)) {
    stop("For spatial data, specify neighb or maxdist.", call. = FALSE)
  }

  if (spacetime) {
    if (!is.numeric(coordt) || !length(coordt) || any(!is.finite(coordt))) {
      stop("coordt must be a non-empty finite numeric vector.", call. = FALSE)
    }
    if (is.null(time) || !is.numeric(time) || !length(time) ||
        any(!is.finite(time))) {
      stop("time must be a non-empty finite numeric vector for ",
           "spatiotemporal neighborhoods.", call. = FALSE)
    }
    coordt <- as.numeric(coordt)
    time <- as.numeric(time)
  }

  if (distance %in% c("Geod", "Chor")) {
    .check_spherical_coords(loc, "loc")
  }

  if (!is_dynamic && is.null(coords)) {
    stop("Observed coordinates are missing.", call. = FALSE)
  }

  if (!is.null(coords) && distance %in% c("Geod", "Chor")) {
    .check_spherical_coords(coords, "observed coordinates")
  }

  if (spacetime && is_dynamic) {
    if (!is.list(coordx_dyn) || length(coordx_dyn) != length(coordt)) {
      stop("For dynamic spatiotemporal data, coordx_dyn must be a list ",
           "with one coordinate matrix for each element of coordt.",
           call. = FALSE)
    }
    coordx_dyn <- lapply(seq_along(coordx_dyn), function(i) {
      z <- .check_numeric_matrix(coordx_dyn[[i]],
                                 paste0("coordx_dyn[[", i, "]]"))
      if (distance %in% c("Geod", "Chor")) {
        .check_spherical_coords(z, paste0("coordx_dyn[[", i, "]]"))
      }
      z
    })
  }

  if (isTRUE(bivariate) && is_dynamic) {
    if (!is.list(coordx_dyn) || length(coordx_dyn) != 2L) {
      stop("For a bivariate field on different locations, coordx_dyn ",
           "must contain two coordinate matrices.", call. = FALSE)
    }
    coordx_dyn <- lapply(seq_len(2L), function(i) {
      z <- .check_numeric_matrix(coordx_dyn[[i]],
                                 paste0("coordx_dyn[[", i, "]]"))
      if (distance %in% c("Geod", "Chor")) {
        .check_spherical_coords(z, paste0("coordx_dyn[[", i, "]]"))
      }
      z
    })
  }

  .GeoCheckDuplicateLocations(
    coordx = coords, coordt = coordt, coordx_dyn = coordx_dyn,
    spacetime = spacetime, bivariate = bivariate, grid = FALSE,
    context = "GeoNeighborhood", check.duplicates = check.duplicates
  )

  ##########################################################################
  ## Component subsetting helpers
  ##########################################################################

  .subset_rows <- function(x, idx, n_total, name) {
    if (is.null(x)) return(NULL)
    if (is.vector(x) && !is.list(x)) {
      if (length(x) != n_total) {
        stop(name, " has an incompatible length.", call. = FALSE)
      }
      return(x[idx])
    }
    x <- as.matrix(x)
    if (nrow(x) != n_total) {
      stop(name, " has an incompatible number of rows.", call. = FALSE)
    }
    x[idx, , drop = FALSE]
  }

  .subset_static_data <- function(x, flat_idx, time_idx, spatial_idx,
                                  nt, ns) {
    if (is.null(x)) return(NULL)
    if (is.matrix(x)) {
      if (identical(dim(x), c(nt, ns))) {
        return(as.numeric(t(x[time_idx, spatial_idx, drop = FALSE])))
      }
      if (identical(dim(x), c(ns, nt))) {
        return(as.numeric(x[spatial_idx, time_idx, drop = FALSE]))
      }
    }
    if (is.vector(x) && !is.list(x) && length(x) == nt * ns) {
      return(x[flat_idx])
    }
    stop("data has an incompatible format for static spatiotemporal data.",
         call. = FALSE)
  }

  .dynamic_piece <- function(x, time_idx, local_idx, offsets, name,
                             matrix_result = FALSE) {
    if (is.null(x)) return(NULL)

    if (is.list(x)) {
      if (length(x) < time_idx) {
        stop(name, " has too few temporal components.", call. = FALSE)
      }
      piece <- x[[time_idx]]
      if (matrix_result) {
        piece <- as.matrix(piece)
        return(piece[local_idx, , drop = FALSE])
      }
      return(piece[local_idx])
    }

    global_idx <- offsets[time_idx] + local_idx
    if (matrix_result) {
      x <- as.matrix(x)
      if (nrow(x) < max(global_idx, 0L)) {
        stop(name, " has an incompatible number of rows.", call. = FALSE)
      }
      return(x[global_idx, , drop = FALSE])
    }

    if (length(x) < max(global_idx, 0L)) {
      stop(name, " has an incompatible length.", call. = FALSE)
    }
    x[global_idx]
  }

  ##########################################################################
  ## Pure spatial processor
  ##########################################################################

  .process_spatial <- function() {
    n_obs <- nrow(coords)
    neighbors <- .find_neighbors(coords, loc, distance, maxdist,
                                 neighb, radius)

    pieces <- .map(seq_len(nrow(loc)), function(i) {
      idx <- neighbors[[i]]$indices
      list(
        coordx = coords[idx, , drop = FALSE],
        indices = idx,
        distances = neighbors[[i]]$distances,
        data = .subset_rows(data, idx, n_obs, "data"),
        X = .subset_rows(X, idx, n_obs, "X"),
        M = .subset_rows(M, idx, n_obs, "M"),
        numpoints = length(idx),
        numcoord = length(idx),
        numtime = 1L,
        block_sizes = length(idx)
      )
    })

    list(
      coordx = lapply(pieces, `[[`, "coordx"),
      coordx_dyn = NULL,
      coordt = NULL,
      coordt_obs = NULL,
      indices = lapply(pieces, `[[`, "indices"),
      indices_by_time = NULL,
      distances = lapply(pieces, `[[`, "distances"),
      data = if (is.null(data)) NULL else lapply(pieces, `[[`, "data"),
      X = if (is.null(X)) NULL else lapply(pieces, `[[`, "X"),
      M = if (is.null(M)) NULL else lapply(pieces, `[[`, "M"),
      numpoints = vapply(pieces, `[[`, integer(1L), "numpoints"),
      numcoord = vapply(pieces, `[[`, integer(1L), "numcoord"),
      numtime = vapply(pieces, `[[`, integer(1L), "numtime"),
      block_sizes = lapply(pieces, `[[`, "block_sizes")
    )
  }

  ##########################################################################
  ## Static spatiotemporal processor
  ##########################################################################

  .process_spacetime_static <- function() {
    ns <- nrow(coords)
    nt <- length(coordt)
    n_total <- ns * nt

    spatial_neighbors <- .find_neighbors(coords, loc, distance, maxdist,
                                         neighb, radius)
    temporal_neighbors <- lapply(time, .select_times,
                                 coordt = coordt, maxtime = maxtime)

    ## Local prediction tasks use location-major order throughout GeoKrigloc:
    ## all requested times for location 1, then all times for location 2, etc.
    ## Do not use expand.grid(location, target_time) here because its first
    ## argument varies fastest and would therefore create time-major blocks.
    task_grid <- data.frame(
      location = rep(seq_len(nrow(loc)), each = length(time)),
      target_time = rep(seq_along(time), times = nrow(loc))
    )

    pieces <- .map(seq_len(nrow(task_grid)), function(k) {
      i <- task_grid$location[k]
      j <- task_grid$target_time[k]
      sidx <- spatial_neighbors[[i]]$indices
      tidx <- temporal_neighbors[[j]]

      flat_idx <- if (length(sidx) && length(tidx)) {
        as.integer(unlist(lapply(tidx, function(tt) sidx + ns * (tt - 1L)),
                          use.names = FALSE))
      } else {
        integer(0L)
      }

      list(
        coordx = coords[sidx, , drop = FALSE],
        coordt = coordt[tidx],
        coordt_obs = rep(coordt[tidx], each = length(sidx)),
        indices = flat_idx,
        indices_by_time = rep(list(sidx), length(tidx)),
        distances = spatial_neighbors[[i]]$distances,
        data = .subset_static_data(data, flat_idx, tidx, sidx, nt, ns),
        X = .subset_rows(X, flat_idx, n_total, "X"),
        M = .subset_rows(M, flat_idx, n_total, "M"),
        numpoints = length(flat_idx),
        numcoord = length(sidx),
        numtime = length(tidx),
        block_sizes = rep(length(sidx), length(tidx))
      )
    })

    list(
      coordx = lapply(pieces, `[[`, "coordx"),
      coordx_dyn = NULL,
      coordt = lapply(pieces, `[[`, "coordt"),
      coordt_obs = lapply(pieces, `[[`, "coordt_obs"),
      indices = lapply(pieces, `[[`, "indices"),
      indices_by_time = lapply(pieces, `[[`, "indices_by_time"),
      distances = lapply(pieces, `[[`, "distances"),
      data = if (is.null(data)) NULL else lapply(pieces, `[[`, "data"),
      X = if (is.null(X)) NULL else lapply(pieces, `[[`, "X"),
      M = if (is.null(M)) NULL else lapply(pieces, `[[`, "M"),
      numpoints = vapply(pieces, `[[`, integer(1L), "numpoints"),
      numcoord = vapply(pieces, `[[`, integer(1L), "numcoord"),
      numtime = vapply(pieces, `[[`, integer(1L), "numtime"),
      block_sizes = lapply(pieces, `[[`, "block_sizes")
    )
  }

  ##########################################################################
  ## Dynamic spatiotemporal processor
  ##########################################################################

  .process_spacetime_dynamic <- function() {
    nt <- length(coordt)
    sizes <- vapply(coordx_dyn, nrow, integer(1L))
    offsets <- c(0L, cumsum(sizes))[seq_len(nt)]
    temporal_neighbors <- lapply(time, .select_times,
                                 coordt = coordt, maxtime = maxtime)

    ## Local prediction tasks use location-major order throughout GeoKrigloc:
    ## all requested times for location 1, then all times for location 2, etc.
    ## Do not use expand.grid(location, target_time) here because its first
    ## argument varies fastest and would therefore create time-major blocks.
    task_grid <- data.frame(
      location = rep(seq_len(nrow(loc)), each = length(time)),
      target_time = rep(seq_along(time), times = nrow(loc))
    )

    pieces <- .map(seq_len(nrow(task_grid)), function(k) {
      i <- task_grid$location[k]
      j <- task_grid$target_time[k]
      tidx <- temporal_neighbors[[j]]

      if (!length(tidx)) {
        return(list(
          coordx = matrix(numeric(0L), nrow = 0L, ncol = ncol(loc)),
          coordx_dyn = list(), coordt = numeric(0L),
          coordt_obs = numeric(0L), indices = integer(0L),
          indices_by_time = list(), distances = list(), data = NULL,
          X = NULL, M = NULL, numpoints = 0L, numcoord = 0L,
          numtime = 0L, block_sizes = integer(0L)
        ))
      }

      blocks <- lapply(tidx, function(tt) {
        nb <- .find_neighbors(
          coordx_dyn[[tt]],
          matrix(loc[i, ], nrow = 1L),
          distance, maxdist, neighb, radius
        )[[1L]]

        idx <- nb$indices
        list(
          coords = coordx_dyn[[tt]][idx, , drop = FALSE],
          local_indices = idx,
          global_indices = offsets[tt] + idx,
          distances = nb$distances,
          data = .dynamic_piece(data, tt, idx, offsets, "data"),
          X = .dynamic_piece(X, tt, idx, offsets, "X", TRUE),
          M = .dynamic_piece(M, tt, idx, offsets, "M")
        )
      })

      block_sizes <- vapply(blocks, function(z) nrow(z$coords), integer(1L))
      keep_blocks <- block_sizes > 0L
      if (!any(keep_blocks)) {
        return(list(
          coordx = matrix(numeric(0L), nrow = 0L, ncol = ncol(loc)),
          coordx_dyn = list(), coordt = numeric(0L),
          coordt_obs = numeric(0L), indices = integer(0L),
          indices_by_time = list(), distances = list(), data = NULL,
          X = NULL, M = NULL, numpoints = 0L, numcoord = 0L,
          numtime = 0L, block_sizes = integer(0L)
        ))
      }
      blocks <- blocks[keep_blocks]
      tidx <- tidx[keep_blocks]
      block_sizes <- block_sizes[keep_blocks]
      coord_blocks <- lapply(blocks, `[[`, "coords")
      flat_coords <- do.call(rbind, coord_blocks)

      list(
        coordx = flat_coords,
        coordx_dyn = coord_blocks,
        coordt = coordt[tidx],
        coordt_obs = rep(coordt[tidx], block_sizes),
        indices = as.integer(unlist(lapply(blocks, `[[`, "global_indices"),
                                    use.names = FALSE)),
        indices_by_time = lapply(blocks, `[[`, "local_indices"),
        distances = lapply(blocks, `[[`, "distances"),
        data = if (is.null(data)) NULL else
          unlist(lapply(blocks, `[[`, "data"), use.names = FALSE),
        X = if (is.null(X)) NULL else do.call(rbind, lapply(blocks, `[[`, "X")),
        M = if (is.null(M)) NULL else
          unlist(lapply(blocks, `[[`, "M"), use.names = FALSE),
        numpoints = sum(block_sizes),
        numcoord = sum(block_sizes),
        numtime = length(tidx),
        block_sizes = block_sizes
      )
    })

    list(
      coordx = lapply(pieces, `[[`, "coordx"),
      coordx_dyn = lapply(pieces, `[[`, "coordx_dyn"),
      coordt = lapply(pieces, `[[`, "coordt"),
      coordt_obs = lapply(pieces, `[[`, "coordt_obs"),
      indices = lapply(pieces, `[[`, "indices"),
      indices_by_time = lapply(pieces, `[[`, "indices_by_time"),
      distances = lapply(pieces, `[[`, "distances"),
      data = if (is.null(data)) NULL else lapply(pieces, `[[`, "data"),
      X = if (is.null(X)) NULL else lapply(pieces, `[[`, "X"),
      M = if (is.null(M)) NULL else lapply(pieces, `[[`, "M"),
      numpoints = vapply(pieces, `[[`, integer(1L), "numpoints"),
      numcoord = vapply(pieces, `[[`, integer(1L), "numcoord"),
      numtime = vapply(pieces, `[[`, integer(1L), "numtime"),
      block_sizes = lapply(pieces, `[[`, "block_sizes")
    )
  }

  ##########################################################################
  ## Bivariate processor
  ##########################################################################

  .process_bivariate <- function() {
    if (!is_dynamic) {
      n1 <- nrow(coords)
      n2 <- n1
      neighbors <- .find_neighbors(coords, loc, distance, maxdist,
                                   neighb, radius)

      pieces <- .map(seq_len(nrow(loc)), function(i) {
        idx <- neighbors[[i]]$indices
        global_idx <- c(idx, n1 + idx)

        data_i <- NULL
        if (!is.null(data)) {
          if (is.matrix(data) && nrow(data) == 2L && ncol(data) == n1) {
            data_i <- c(data[1L, idx], data[2L, idx])
          } else if (is.list(data) && length(data) == 2L) {
            data_i <- c(data[[1L]][idx], data[[2L]][idx])
          } else if (is.vector(data) && length(data) == n1 + n2) {
            data_i <- data[global_idx]
          } else {
            stop("data has an incompatible format for bivariate data.",
                 call. = FALSE)
          }
        }

        M_i <- NULL
        if (!is.null(M)) {
          if (is.list(M) && length(M) == 2L) {
            M_i <- c(M[[1L]][idx], M[[2L]][idx])
          } else if (is.vector(M) && length(M) == n1 + n2) {
            M_i <- M[global_idx]
          } else if (is.vector(M) && length(M) == 2L) {
            M_i <- M
          } else {
            stop("M has an incompatible format for bivariate data.",
                 call. = FALSE)
          }
        }

        X_i <- NULL
        if (!is.null(X)) {
          if (is.list(X) && length(X) == 2L) {
            X_i <- rbind(as.matrix(X[[1L]])[idx, , drop = FALSE],
                         as.matrix(X[[2L]])[idx, , drop = FALSE])
          } else {
            Xm <- as.matrix(X)
            if (nrow(Xm) == n1 + n2) {
              X_i <- Xm[global_idx, , drop = FALSE]
            } else if (nrow(Xm) == n1) {
              X_i <- rbind(Xm[idx, , drop = FALSE],
                           Xm[idx, , drop = FALSE])
            } else {
              stop("X has an incompatible number of rows for bivariate data.",
                   call. = FALSE)
            }
          }
        }

        list(
          coordx = coords[idx, , drop = FALSE],
          coordx_dyn = list(coords[idx, , drop = FALSE],
                            coords[idx, , drop = FALSE]),
          indices = global_idx,
          indices_by_variable = list(idx, idx),
          distances = list(neighbors[[i]]$distances,
                           neighbors[[i]]$distances),
          data = data_i, X = X_i, M = M_i,
          numpoints = length(global_idx),
          numcoord = length(idx), numtime = 1L,
          block_sizes = c(length(idx), length(idx))
        )
      })
    } else {
      coords1 <- coordx_dyn[[1L]]
      coords2 <- coordx_dyn[[2L]]
      n1 <- nrow(coords1)
      n2 <- nrow(coords2)
      neighbors1 <- .find_neighbors(coords1, loc, distance, maxdist,
                                    neighb, radius)
      neighbors2 <- .find_neighbors(coords2, loc, distance, maxdist,
                                    neighb, radius)

      pieces <- .map(seq_len(nrow(loc)), function(i) {
        idx1 <- neighbors1[[i]]$indices
        idx2 <- neighbors2[[i]]$indices
        global_idx <- c(idx1, n1 + idx2)
        block1 <- coords1[idx1, , drop = FALSE]
        block2 <- coords2[idx2, , drop = FALSE]

        data_i <- NULL
        if (!is.null(data)) {
          if (is.list(data) && length(data) == 2L) {
            data_i <- c(data[[1L]][idx1], data[[2L]][idx2])
          } else if (is.vector(data) && length(data) == n1 + n2) {
            data_i <- data[global_idx]
          } else {
            stop("data has an incompatible format for bivariate data ",
                 "on different locations.", call. = FALSE)
          }
        }

        M_i <- NULL
        if (!is.null(M)) {
          if (is.list(M) && length(M) == 2L) {
            M_i <- c(M[[1L]][idx1], M[[2L]][idx2])
          } else if (is.vector(M) && length(M) == n1 + n2) {
            M_i <- M[global_idx]
          } else if (is.vector(M) && length(M) == 2L) {
            M_i <- M
          } else {
            stop("M has an incompatible format for bivariate data ",
                 "on different locations.", call. = FALSE)
          }
        }

        X_i <- NULL
        if (!is.null(X)) {
          if (is.list(X) && length(X) == 2L) {
            X_i <- rbind(as.matrix(X[[1L]])[idx1, , drop = FALSE],
                         as.matrix(X[[2L]])[idx2, , drop = FALSE])
          } else {
            Xm <- as.matrix(X)
            if (nrow(Xm) != n1 + n2) {
              stop("X has an incompatible number of rows for bivariate data ",
                   "on different locations.", call. = FALSE)
            }
            X_i <- Xm[global_idx, , drop = FALSE]
          }
        }

        list(
          coordx = rbind(block1, block2),
          coordx_dyn = list(block1, block2),
          indices = global_idx,
          indices_by_variable = list(idx1, idx2),
          distances = list(neighbors1[[i]]$distances,
                           neighbors2[[i]]$distances),
          data = data_i, X = X_i, M = M_i,
          numpoints = length(global_idx),
          numcoord = length(global_idx), numtime = 1L,
          block_sizes = c(length(idx1), length(idx2))
        )
      })
    }

    list(
      coordx = lapply(pieces, `[[`, "coordx"),
      coordx_dyn = lapply(pieces, `[[`, "coordx_dyn"),
      coordt = NULL,
      coordt_obs = NULL,
      indices = lapply(pieces, `[[`, "indices"),
      indices_by_time = lapply(pieces, `[[`, "indices_by_variable"),
      distances = lapply(pieces, `[[`, "distances"),
      data = if (is.null(data)) NULL else lapply(pieces, `[[`, "data"),
      X = if (is.null(X)) NULL else lapply(pieces, `[[`, "X"),
      M = if (is.null(M)) NULL else lapply(pieces, `[[`, "M"),
      numpoints = vapply(pieces, `[[`, integer(1L), "numpoints"),
      numcoord = vapply(pieces, `[[`, integer(1L), "numcoord"),
      numtime = vapply(pieces, `[[`, integer(1L), "numtime"),
      block_sizes = lapply(pieces, `[[`, "block_sizes")
    )
  }

  ##########################################################################
  ## Dispatch and output
  ##########################################################################

  results <- if (isTRUE(bivariate)) {
    .process_bivariate()
  } else if (spacetime && is_dynamic) {
    .process_spacetime_dynamic()
  } else if (spacetime) {
    .process_spacetime_static()
  } else {
    .process_spatial()
  }

  list(
    data = results$data,
    coordx = results$coordx,
    coordx_dyn = results$coordx_dyn,
    coordt = results$coordt,
    coordt_obs = results$coordt_obs,
    indices = results$indices,
    indices_by_time = results$indices_by_time,
    distances = results$distances,
    distance = distance,
    numpoints = results$numpoints,
    numcoord = results$numcoord,
    numtime = results$numtime,
    block_sizes = results$block_sizes,
    radius = radius,
    spacetime = spacetime,
    bivariate = bivariate,
    dynamic = is_dynamic,
    X = results$X,
    M = results$M
  )
}
