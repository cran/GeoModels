
GeoNeighborhood <- function(data = NULL, coordx = NULL, coordy = NULL, coordz = NULL,
                           coordt = NULL, coordx_dyn = NULL, bivariate = FALSE,
                           distance = "Eucl", grid = FALSE, loc, neighb = NULL,
                           maxdist = NULL, maxtime = NULL, radius = 1, time = NULL,
                           X = NULL, M = NULL, spobj = NULL, spdata = NULL,
                           parallel = FALSE, ncores = NULL) {
  
  # =============================================================================
  # Input Validation and Preprocessing
  # =============================================================================
  # Basic validation
  if (missing(loc)) stop("Parameter 'loc' is required")
  if (!is.logical(bivariate)) stop("Parameter 'bivariate' must be logical")
  if (is.null(neighb) && is.null(maxdist)) {
    stop("Either 'maxdist' (maxtime) or 'neighb' must be specified")
  }
  if (!distance %in% c("Eucl", "Geod", "Chor")) {
    stop("Distance must be one of: 'Eucl', 'Geod', 'Chor'")
  }
  # Standardize locations
  if (is.vector(loc)) loc <- t(as.matrix(loc))
  if (!is.matrix(loc)) stop("Parameter 'loc' must be a matrix")
  colnames(loc) <- NULL
  # Determine analysis type
  spacetime <- !is.null(coordt)
  if (spacetime && is.null(time)) stop("Parameter 'time' is required for spacetime analysis")
  if (spacetime && !is.vector(time)) stop("Parameter 'time' must be a vector")
  
  # Check for dynamic coordinates
  is_dynamic <- !is.null(coordx_dyn)
  if (is_dynamic && spacetime && length(coordx_dyn) != length(data)) {
    stop("Length of 'coordx_dyn' and 'data' must match for dynamic spatiotemporal analysis")
  }
  # Set default neighbor count
  if (!is.null(neighb)) neighb <- round(neighb)
  # =============================================================================
  # Handle Spatial Objects
  # =============================================================================
  if (!is.null(spobj)) {
    if (spacetime) {
      spatial_info <- sp2Geo(spobj, spdata)
      coordx <- spatial_info$coords
      coordt <- spatial_info$coordt
      if (!spatial_info$pj && distance != "Chor") distance <- "Geod"
    } else {
      spatial_info <- sp2Geo(spobj, spdata)
      coordx <- spatial_info$coords
      if (!spatial_info$pj && distance != "Chor") distance <- "Geod"
    }
    
    if (!is.null(spatial_info$Y) && !is.null(spatial_info$X)) {
      data <- spatial_info$Y
      X <- spatial_info$X
    }
  }
  # =============================================================================
  # Prepare Coordinates
  # =============================================================================
  coords <- NULL
  if (!is.null(coordx)) {
    if (is.null(coordy) && is.null(coordz)) {
      coords <- as.matrix(coordx)
    } else {
      if (grid) {
        coords <- as.matrix(expand.grid(coordx, coordy, coordz))
      } else {
        coords <- as.matrix(cbind(coordx, coordy, coordz))
      }
    }
    colnames(coords) <- NULL
  }
  
  # =============================================================================
  # Main Processing Logic
  # =============================================================================
  
  num_locations <- nrow(loc)
  
  # Initialize result variables
  selected_coords <- NULL
  selected_times <- NULL
  selected_data <- NULL
  selected_X <- NULL
  selected_M <- NULL
  num_points <- NULL
  num_time <- 1
  
  # =============================================================================
  # SPATIAL ANALYSIS
  # =============================================================================
  
  if (!spacetime && !bivariate) {
    
    num_coords <- nrow(coords)
    if (is.null(neighb)) neighb <- min(100, num_coords)
    
    # Find spatial neighbors
    neighbors <- find_neighbors(coords, loc, distance, maxdist, neighb, radius)
    neighbor_indices <- neighbors$nn.idx
    
    # Pre-allocate result lists
    selected_coords <- vector("list", num_locations)
    selected_data <- if (!is.null(data)) vector("list", num_locations) else NULL
    selected_X <- if (!is.null(X)) vector("list", num_locations) else NULL
    selected_M <- if (!is.null(M)) vector("list", num_locations) else NULL
    num_points <- integer(num_locations)
    
    # Process each location
    for (i in seq_len(num_locations)) {
      indices <- neighbor_indices[i, ]
      
      selected_coords[[i]] <- coords[indices, , drop = FALSE]
      num_points[i] <- length(indices)
      
      if (!is.null(data)) selected_data[[i]] <- data[indices]
      if (!is.null(M)) selected_M[[i]] <- M[indices]
      if (!is.null(X)) selected_X[[i]] <- X[indices, , drop = FALSE]
    }
  }
  # =============================================================================
  # SPATIOTEMPORAL ANALYSIS
  # =============================================================================
  if (spacetime) {
    
    if (!is_dynamic) {
      # Static spatiotemporal case
      num_times <- length(time)
      num_coords <- nrow(coords)
      num_time_coords <- length(coordt)
      
      if (is.null(neighb)) neighb <- min(100, num_coords)
      
      # Pre-allocate matrices for efficient access
      if (!is.null(M)) {
        M_matrix <- matrix(M, nrow = num_time_coords, ncol = num_coords, byrow = TRUE)
      }
      if (!is.null(X)) {
        X_arrays <- lapply(seq_len(ncol(X)), function(j) {
          matrix(X[, j], nrow = num_time_coords, ncol = num_coords, byrow = TRUE)
        })
      }
      
      # Find spatial neighbors
      spatial_neighbors <- find_neighbors(coords, loc, distance, maxdist, neighb, radius)
      
      # Find temporal neighbors
      temporal_neighbors <- nabor::knn(coordt, time, k = length(coordt), radius = maxtime)
      
      # Process temporal indices
      temporal_query <- temporal_neighbors$nn.idx[1, ]
      valid_temporal <- which(temporal_query != 0)
      temporal_indices_with_times <- cbind(
        as.vector(temporal_neighbors$nn.idx[, valid_temporal]), 
        time
      )
      
      # Pre-allocate results
      total_points <- num_locations * num_times
      selected_coords <- vector("list", total_points)
      selected_times <- vector("list", total_points)
      selected_data <- if (!is.null(data)) vector("list", total_points) else NULL
      selected_X <- if (!is.null(X)) vector("list", total_points) else NULL
      selected_M <- if (!is.null(M)) vector("list", total_points) else NULL
      num_points <- integer(total_points)
      
      # Process each space-time combination
      k <- 1
      for (i in seq_len(num_locations)) {
        spatial_indices <- spatial_neighbors$nn.idx[i, ]
        selected_spatial_coords <- matrix(coords[spatial_indices, ], ncol = ncol(coords))
        
        for (j in seq_len(num_times)) {
          # Find temporal neighbors for this time point
          time_mask <- temporal_indices_with_times[, 2] == time[j]
          
          if (any(time_mask)) {
            current_temporal_indices <- temporal_indices_with_times[time_mask, 1]
            temporal_order <- order(coordt[current_temporal_indices])
            ordered_temporal_indices <- current_temporal_indices[temporal_order]
            
            selected_coords[[k]] <- selected_spatial_coords
            selected_times[[k]] <- coordt[ordered_temporal_indices]
            num_points[k] <- nrow(selected_spatial_coords)
            
            if (!is.null(data)) {
              selected_data[[k]] <- data[ordered_temporal_indices, spatial_indices, drop = FALSE]
              selected_data[[k]] <- selected_data[[k]][order(coordt[ordered_temporal_indices]), , drop = FALSE]
            }
            if (!is.null(M)) {
              M_subset <- M_matrix[ordered_temporal_indices, spatial_indices, drop = FALSE]
              M_subset <- M_subset[order(coordt[ordered_temporal_indices]), , drop = FALSE]
              selected_M[[k]] <- as.numeric(t(M_subset))
            }
            if (!is.null(X)) {
              X_combined <- do.call(cbind, lapply(X_arrays, function(x_array) {
                X_subset <- x_array[ordered_temporal_indices, spatial_indices, drop = FALSE]
                X_subset <- X_subset[order(coordt[ordered_temporal_indices]), , drop = FALSE]
                as.numeric(t(X_subset))
              }))
              selected_X[[k]] <- X_combined
            }
          } else {
            # No temporal neighbors for this time point
            selected_coords[[k]] <- selected_spatial_coords
            selected_times[[k]] <- NULL
            num_points[k] <- nrow(selected_spatial_coords)
          }
          k <- k + 1
        }
      }
      
    } else {
      # Dynamic spatiotemporal case
      num_times <- length(time)
      
      # Find temporal neighbors
      temporal_neighbors <- nabor::knn(coordt, time, k = length(coordt), radius = maxtime)
      temporal_indices <- temporal_neighbors$nn.idx[1, ]
      temporal_indices <- temporal_indices[temporal_indices > 0]
      temporal_order <- order(coordt[temporal_indices])
      ordered_temporal_indices <- temporal_indices[temporal_order]
      ordered_times <- coordt[ordered_temporal_indices]
      
      # Pre-allocate results
      total_points <- num_locations * num_times
      selected_coords <- vector("list", total_points)
      selected_times <- vector("list", total_points)
      selected_data <- if (!is.null(data)) vector("list", total_points) else NULL
      selected_X <- if (!is.null(X)) vector("list", total_points) else NULL
      selected_M <- if (!is.null(M)) vector("list", total_points) else NULL
      num_points <- integer(total_points)
      
      # Process each space-time combination
      k <- 1
      for (i in seq_len(num_locations)) {
        for (j in seq_len(num_times)) {
          # For each temporal index, find spatial neighbors using dynamic coordinates
          spatial_results <- lapply(ordered_temporal_indices, function(t_idx) {
            current_coords <- coordx_dyn[[t_idx]]
            current_location <- matrix(loc[i, ], nrow = 1)
            
            spatial_neighbors <- find_neighbors(current_coords, current_location, 
                                              distance, maxdist, neighb, radius)
            spatial_indices <- spatial_neighbors$nn.idx[1, ]
            
            list(
              coords = current_coords[spatial_indices, , drop = FALSE],
              data = if (!is.null(data)) data[[t_idx]][spatial_indices] else NULL,
              M = if (!is.null(M)) M[[t_idx]][spatial_indices] else NULL,
              X = if (!is.null(X)) X[[t_idx]][spatial_indices, , drop = FALSE] else NULL
            )
          })
          
          # Combine results
          selected_coords[[k]] <- spatial_results[[1]]$coords
          selected_times[[k]] <- ordered_times
          num_points[k] <- nrow(spatial_results[[1]]$coords)
          
          if (!is.null(data)) {
            selected_data[[k]] <- do.call(rbind, lapply(spatial_results, function(x) x$data))
          }
          if (!is.null(M)) {
            selected_M[[k]] <- unlist(lapply(spatial_results, function(x) x$M))
          }
          if (!is.null(X)) {
            selected_X[[k]] <- do.call(rbind, lapply(spatial_results, function(x) x$X))
          }
          
          k <- k + 1
        }
      }
    }
    
    num_time <- length(time)
  }
  # =============================================================================
  # BIVARIATE ANALYSIS
  # =============================================================================
  if (bivariate) {
    
    # Pre-allocate results
    selected_coords <- vector("list", num_locations)
    selected_data <- if (!is.null(data)) vector("list", num_locations) else NULL
    selected_X <- if (!is.null(X)) vector("list", num_locations) else NULL
    selected_M <- if (!is.null(M)) vector("list", num_locations) else NULL
    num_points <- integer(num_locations)
    
    if (is_dynamic) {
      # Dynamic bivariate case
      coords1 <- coords[[1]]
      coords2 <- coords[[2]]
      num_coords1 <- nrow(coords1)
      
      neighbors1 <- find_neighbors(coords1, loc, distance, maxdist, neighb, radius)
      neighbors2 <- find_neighbors(coords2, loc, distance, maxdist, neighb, radius)
      
      for (i in seq_len(num_locations)) {
        indices1 <- neighbors1$nn.idx[i, ]
        indices2 <- neighbors2$nn.idx[i, ]
        
        selected_coords[[i]] <- rbind(coords1[indices1, ], coords2[indices2, ])
        num_points[i] <- length(indices1) + length(indices2)
        
        if (!is.null(data)) {
          selected_data[[i]] <- matrix(c(data[1, indices1], data[2, indices2]), nrow = 2)
        }
        if (!is.null(M)) {
          selected_M[[i]] <- list(M[[1]][indices1], M[[2]][indices2])
        }
        if (!is.null(X)) {
          selected_X[[i]] <- list(X[[1]][indices1, , drop = FALSE], 
                                 X[[2]][indices2, , drop = FALSE])
        }
      }
    } else {
      # Static bivariate case
      num_coords <- nrow(coords)
      neighbors <- find_neighbors(coords, loc, distance, maxdist, neighb, radius)
      
      for (i in seq_len(num_locations)) {
        indices <- neighbors$nn.idx[i, ]
        
        selected_coords[[i]] <- coords[indices, ]
        num_points[i] <- length(indices)
        
        if (!is.null(data)) {
          selected_data[[i]] <- matrix(data[, indices], nrow = 2)
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
  }
  
  # =============================================================================
  # Clean up empty lists
  # =============================================================================
  
  if (length(selected_M) == 0) selected_M <- NULL
  if (length(selected_X) == 0) selected_X <- NULL
  if (length(selected_data) == 0) selected_data <- NULL
  
  # =============================================================================
  # Return Results
  # =============================================================================
  
  return(list(
    data = selected_data,
    coordx = selected_coords,
    coordt = selected_times,
    distance = distance,
    numpoints = num_points,
    numtime = num_time,
    radius = radius,
    spacetime = spacetime,
    X = selected_X,
    M = selected_M
  ))
}

# =============================================================================
# Helper Functions
# =============================================================================
find_neighbors <- function(coords, locations, distance, maxdist, neighb, radius) {
  
  # Set default search parameters
  if (is.null(maxdist)) maxdist <- 0
  if (is.null(neighb)) neighb <- min(100, nrow(coords))
  
  # Handle geodetic and chord distances with projection
  if (distance %in% c("Geod", "Chor")) {
    # Project coordinates using sinusoidal projection
    coords_proj <- mapproj::mapproject(coords[, 1], coords[, 2], projection = "sinusoidal")
    coords_projected <- radius * cbind(coords_proj$x, coords_proj$y)
    
    loc_proj <- mapproj::mapproject(locations[, 1], locations[, 2], projection = "sinusoidal")
    locations_projected <- radius * cbind(loc_proj$x, loc_proj$y)
    
    result <- nabor::knn(coords_projected, locations_projected, 
                        k = neighb, radius = maxdist)
  } else {
    # Euclidean distance
    result <- nabor::knn(coords, locations, k = neighb, radius = maxdist)
  }
  
  return(result)
}


sp2Geo <- function(spobj, spdata) {
  # This is a placeholder - the actual implementation would depend on the sp package
  # and the specific spatial object types being used
  warning("sp2Geo function is not implemented - placeholder used")
  return(list(
    coords = NULL,
    coordt = NULL,
    Y = NULL,
    X = NULL,
    pj = TRUE  # Assuming projected coordinates
  ))
}