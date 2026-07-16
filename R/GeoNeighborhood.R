GeoNeighborhood <- function(data = NULL, coordx = NULL, coordy = NULL, coordz = NULL,
                            coordt = NULL, coordx_dyn = NULL, bivariate = FALSE,
                            distance = "Eucl", grid = FALSE, loc, neighb = NULL,
                            maxdist = NULL, maxtime = NULL, radius = 1, time = NULL,
                            X = NULL, M = NULL, spobj = NULL, spdata = NULL,
                            parallel = FALSE, ncores = NULL) {

  .valid_idx <- function(idx) {
    idx <- as.integer(idx)
    idx[idx > 0L & !is.na(idx)]
  }

  .knn_radius <- function(query, ref, k, rad = NULL) {
    if (is.null(rad) || is.infinite(rad)) {
      nabor::knn(ref, query, k = k)
    } else {
      nabor::knn(ref, query, k = k, radius = rad)
    }
  }

  find_neighbors <- function(coords, locations, distance, maxdist, neighb, radius) {
    coords <- as.matrix(coords)
    locations <- as.matrix(locations)
    if (is.null(neighb)) neighb <- min(100L, nrow(coords))
    neighb <- min(as.integer(round(neighb)), nrow(coords))
    if (neighb < 1L) stop("'neighb' must be a positive integer")

    if (distance %in% c("Geod", "Chor")) {
      ## Approximate neighbor search for lon/lat using 3D Cartesian coordinates.
      ## Exact model distances are still handled downstream by GeoModels.
      to_xyz <- function(mat) {
        mat <- as.matrix(mat)
        lat <- mat[, 2L] * pi / 180
        lon <- mat[, 1L] * pi / 180
        cbind(cos(lat) * cos(lon),
              cos(lat) * sin(lon),
              sin(lat)) * radius
      }
      coords_3d <- to_xyz(coords)
      loc_3d <- to_xyz(locations)
      result <- .knn_radius(loc_3d, coords_3d, k = neighb, rad = maxdist)
    } else {
      result <- .knn_radius(locations, coords, k = neighb, rad = maxdist)
    }
    result
  }

  #################
  # Input validation
  #################
  if (missing(loc)) stop("Parameter 'loc' is required")
  if (!is.logical(bivariate)) stop("Parameter 'bivariate' must be logical")
  if (is.null(neighb) && is.null(maxdist)) {
    stop("Either 'maxdist' (maxtime) or 'neighb' must be specified")
  }
  if (!distance %in% c("Eucl", "Geod", "Chor")) {
    stop("Distance must be one of: 'Eucl', 'Geod', 'Chor'")
  }

  if (is.vector(loc)) loc <- t(as.matrix(loc))
  if (!is.matrix(loc)) stop("Parameter 'loc' must be a matrix")
  colnames(loc) <- NULL

  spacetime <- !is.null(coordt)
  if (spacetime && is.null(time)) stop("Parameter 'time' is required for spacetime analysis")
  if (spacetime && !is.vector(time)) stop("Parameter 'time' must be a vector")

  is_dynamic <- !is.null(coordx_dyn)
  if (is_dynamic && spacetime && !is.null(data) && length(coordx_dyn) != length(data)) {
    stop("Length of 'coordx_dyn' and 'data' must match for dynamic spatiotemporal analysis")
  }
  if (!is.null(neighb)) neighb <- round(neighb)

  #################
  # Handle Spatial Objects
  #################
  if (!is.null(spobj)) {
    spatial_info <- sp2Geo(spobj, spdata)
    coordx <- spatial_info$coords
    if (spacetime) coordt <- spatial_info$coordt
    if (!spatial_info$pj && distance != "Chor") distance <- "Geod"

    if (!is.null(spatial_info$Y) && !is.null(spatial_info$X)) {
      data <- spatial_info$Y
      X <- spatial_info$X
    }
  }

  #################
  # Prepare Coordinates
  #################
  coords <- NULL
  if (!is.null(coordx)) {
    if (is.null(coordy) && is.null(coordz)) {
      coords <- as.matrix(coordx)
    } else {
      coords <- if (grid) as.matrix(expand.grid(coordx, coordy, coordz)) else as.matrix(cbind(coordx, coordy, coordz))
    }
    colnames(coords) <- NULL
  }

  #################
  # PROCESSORS
  #################

  ## ---- Pure spatial ----
  process_spatial <- function() {
    num_locations <- nrow(loc)
    num_coords <- nrow(coords)
    if (is.null(neighb)) neighb <- min(100L, num_coords)

    neighbors <- find_neighbors(coords, loc, distance, maxdist, neighb, radius)
    neighbor_indices <- neighbors$nn.idx

    selected_coords <- vector("list", num_locations)
    selected_data <- if (!is.null(data)) vector("list", num_locations) else NULL
    selected_X <- if (!is.null(X)) vector("list", num_locations) else NULL
    selected_M <- if (!is.null(M)) vector("list", num_locations) else NULL
    num_points <- integer(num_locations)

    for (i in seq_len(num_locations)) {
      indices <- .valid_idx(neighbor_indices[i, ])
      selected_coords[[i]] <- coords[indices, , drop = FALSE]
      num_points[i] <- length(indices)
      if (!is.null(data)) selected_data[[i]] <- data[indices]
      if (!is.null(M)) selected_M[[i]] <- M[indices]
      if (!is.null(X)) selected_X[[i]] <- X[indices, , drop = FALSE]
    }
    list(coords = selected_coords, times = NULL, data = selected_data,
         X = selected_X, M = selected_M, num_points = num_points, num_time = 1)
  }

  ## ---- Spatiotemporal static ----
  process_spatiotemporal_static <- function() {
    num_locations <- nrow(loc)
    num_times <- length(time)
    num_coords <- nrow(coords)
    num_time_coords <- length(coordt)
    if (is.null(neighb)) neighb <- min(100L, num_coords)

    if (!is.null(M)) {
      M_matrix <- matrix(M, nrow = num_time_coords, ncol = num_coords, byrow = TRUE)
    }
    if (!is.null(X)) {
      X_arrays <- lapply(seq_len(ncol(X)), function(j) {
        matrix(X[, j], nrow = num_time_coords, ncol = num_coords, byrow = TRUE)
      })
    }

    spatial_neighbors <- find_neighbors(coords, loc, distance, maxdist, neighb, radius)
    max_temporal_k <- if (is.null(neighb)) length(coordt) else min(neighb, length(coordt))

    temporal_results <- lapply(time, function(t) {
      temp_neighbors <- .knn_radius(matrix(t, ncol = 1L), matrix(coordt, ncol = 1L),
                                    k = max_temporal_k, rad = maxtime)
      indices <- .valid_idx(temp_neighbors$nn.idx[1, ])
      if (length(indices) == 0L) {
        closest <- nabor::knn(matrix(coordt, ncol = 1L), matrix(t, ncol = 1L), k = 1L)
        indices <- closest$nn.idx[1, 1]
      }
      indices
    })

    total_points <- num_locations * num_times
    selected_coords <- vector("list", total_points)
    selected_times <- vector("list", total_points)
    selected_data <- if (!is.null(data)) vector("list", total_points) else NULL
    selected_X <- if (!is.null(X)) vector("list", total_points) else NULL
    selected_M <- if (!is.null(M)) vector("list", total_points) else NULL
    num_points <- integer(total_points)

    k <- 1L
    for (i in seq_len(num_locations)) {
      spatial_indices <- .valid_idx(spatial_neighbors$nn.idx[i, ])
      selected_spatial_coords <- coords[spatial_indices, , drop = FALSE]

      for (j in seq_len(num_times)) {
        current_temporal_indices <- temporal_results[[j]]
        temporal_order <- order(coordt[current_temporal_indices])
        ordered_temporal_indices <- current_temporal_indices[temporal_order]

        selected_coords[[k]] <- selected_spatial_coords
        selected_times[[k]] <- coordt[ordered_temporal_indices]
        num_points[k] <- nrow(selected_spatial_coords)

        if (!is.null(data)) {
          selected_data[[k]] <- data[ordered_temporal_indices, spatial_indices, drop = FALSE]
        }
        if (!is.null(M)) {
          M_subset <- M_matrix[ordered_temporal_indices, spatial_indices, drop = FALSE]
          selected_M[[k]] <- as.numeric(t(M_subset))
        }
        if (!is.null(X)) {
          X_combined <- do.call(cbind, lapply(X_arrays, function(x_array) {
            X_subset <- x_array[ordered_temporal_indices, spatial_indices, drop = FALSE]
            as.numeric(t(X_subset))
          }))
          selected_X[[k]] <- X_combined
        }
        k <- k + 1L
      }
    }

    list(coords = selected_coords, times = selected_times, data = selected_data,
         X = selected_X, M = selected_M, num_points = num_points, num_time = length(time))
  }

  ## ---- Spatiotemporal dynamic ----
  process_spatiotemporal_dynamic <- function() {
    num_locations <- nrow(loc)
    num_times <- length(time)
    max_temporal_k <- if (is.null(neighb)) length(coordt) else min(neighb, length(coordt))

    temporal_neighbors <- .knn_radius(matrix(time, ncol = 1L), matrix(coordt, ncol = 1L),
                                      k = max_temporal_k, rad = maxtime)

    total_points <- num_locations * num_times
    selected_coords <- vector("list", total_points)
    selected_times <- vector("list", total_points)
    selected_data <- if (!is.null(data)) vector("list", total_points) else NULL
    selected_X <- if (!is.null(X)) vector("list", total_points) else NULL
    selected_M <- if (!is.null(M)) vector("list", total_points) else NULL
    num_points <- integer(total_points)

    k <- 1L
    for (i in seq_len(num_locations)) {
      for (j in seq_len(num_times)) {
        temporal_indices <- .valid_idx(temporal_neighbors$nn.idx[j, ])
        if (length(temporal_indices) == 0L) {
          closest <- nabor::knn(matrix(coordt, ncol = 1L), matrix(time[j], ncol = 1L), k = 1L)
          temporal_indices <- closest$nn.idx[1, 1]
        }
        temporal_order <- order(coordt[temporal_indices])
        ordered_temporal_indices <- temporal_indices[temporal_order]
        ordered_times <- coordt[ordered_temporal_indices]

        spatial_results <- lapply(ordered_temporal_indices, function(t_idx) {
          current_coords <- as.matrix(coordx_dyn[[t_idx]])
          current_location <- matrix(loc[i, ], nrow = 1L)
          spatial_neighbors <- find_neighbors(current_coords, current_location, distance, maxdist, neighb, radius)
          spatial_indices <- .valid_idx(spatial_neighbors$nn.idx[1, ])

          list(
            coords = current_coords[spatial_indices, , drop = FALSE],
            data = if (!is.null(data)) data[[t_idx]][spatial_indices] else NULL,
            M = if (!is.null(M)) M[[t_idx]][spatial_indices] else NULL,
            X = if (!is.null(X)) X[[t_idx]][spatial_indices, , drop = FALSE] else NULL
          )
        })

        selected_coords[[k]] <- spatial_results[[1L]]$coords
        selected_times[[k]] <- ordered_times
        num_points[k] <- nrow(spatial_results[[1L]]$coords)

        if (!is.null(data)) {
          selected_data[[k]] <- do.call(rbind, lapply(spatial_results, function(x) x$data))
        }
        if (!is.null(M)) {
          selected_M[[k]] <- unlist(lapply(spatial_results, function(x) x$M))
        }
        if (!is.null(X)) {
          selected_X[[k]] <- do.call(rbind, lapply(spatial_results, function(x) x$X))
        }

        k <- k + 1L
      }
    }

    list(coords = selected_coords, times = selected_times, data = selected_data,
         X = selected_X, M = selected_M, num_points = num_points, num_time = length(time))
  }

  ## ---- Bivariate ----
  process_bivariate <- function() {
    num_locations <- nrow(loc)
    selected_coords <- vector("list", num_locations)
    selected_data <- if (!is.null(data)) vector("list", num_locations) else NULL
    selected_X <- if (!is.null(X)) vector("list", num_locations) else NULL
    selected_M <- if (!is.null(M)) vector("list", num_locations) else NULL
    num_points <- integer(num_locations)

    if (is_dynamic) {
      coords1 <- as.matrix(coordx_dyn[[1L]])
      coords2 <- as.matrix(coordx_dyn[[2L]])
      neighbors1 <- find_neighbors(coords1, loc, distance, maxdist, neighb, radius)
      neighbors2 <- find_neighbors(coords2, loc, distance, maxdist, neighb, radius)

      for (i in seq_len(num_locations)) {
        indices1 <- .valid_idx(neighbors1$nn.idx[i, ])
        indices2 <- .valid_idx(neighbors2$nn.idx[i, ])

        selected_coords[[i]] <- rbind(coords1[indices1, , drop = FALSE],
                                      coords2[indices2, , drop = FALSE])
        num_points[i] <- length(indices1) + length(indices2)

        if (!is.null(data)) {
          selected_data[[i]] <- matrix(c(data[1, indices1], data[2, indices2]), nrow = 2L)
        }
        if (!is.null(M)) {
          selected_M[[i]] <- list(M[[1L]][indices1], M[[2L]][indices2])
        }
        if (!is.null(X)) {
          selected_X[[i]] <- list(X[[1L]][indices1, , drop = FALSE],
                                  X[[2L]][indices2, , drop = FALSE])
        }
      }
    } else {
      num_coords <- nrow(coords)
      neighbors <- find_neighbors(coords, loc, distance, maxdist, neighb, radius)

      for (i in seq_len(num_locations)) {
        indices <- .valid_idx(neighbors$nn.idx[i, ])
        selected_coords[[i]] <- coords[indices, , drop = FALSE]
        num_points[i] <- length(indices)

        if (!is.null(data)) {
          selected_data[[i]] <- matrix(data[, indices], nrow = 2L)
        }
        if (!is.null(M)) {
          selected_M[[i]] <- c(M[indices], M[num_coords + indices])
        }
        if (!is.null(X)) {
          selected_X[[i]] <- rbind(X[indices, , drop = FALSE],
                                   X[num_coords + indices, , drop = FALSE])
        }
      }
    }
    list(coords = selected_coords, times = NULL, data = selected_data,
         X = selected_X, M = selected_M, num_points = num_points, num_time = 1)
  }

  #################
  # MAIN LOGIC
  #################
  if (!spacetime && !bivariate) {
    results <- process_spatial()
  } else if (spacetime && !is_dynamic) {
    results <- process_spatiotemporal_static()
  } else if (spacetime && is_dynamic) {
    results <- process_spatiotemporal_dynamic()
  } else if (bivariate) {
    results <- process_bivariate()
  }

  #################
  # Return results
  #################
  list(
    data = results$data,
    coordx = results$coords,
    coordt = results$times,
    distance = distance,
    numpoints = results$num_points,
    numtime = results$num_time,
    radius = radius,
    spacetime = spacetime,
    X = results$X,
    M = results$M
  )
}