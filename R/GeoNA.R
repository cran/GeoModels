
GeoNA <- function(data, coordx = NULL, coordy = NULL, coordz = NULL, coordt = NULL, 
                  coordx_dyn = NULL, grid = FALSE, X = NULL, setting = "spatial") {
  

#' Validate function inputs
validate_inputs <- function(data, coordx, coordy, coordz, coordt, coordx_dyn, setting) {
  
  # Check setting
  if (!setting %in% c("spatial", "spacetime", "bivariate")) {
    stop("'setting' must be one of: 'spatial', 'spacetime', 'bivariate'")
  }
  
  # Check required inputs for spacetime
  if (setting == "spacetime" && is.null(coordt)) {
    stop("'coordt' is required for spacetime setting")
  }
  
  # Check coordx is provided
  if (is.null(coordx) && is.null(coordx_dyn)) {
    stop("Either 'coordx' or 'coordx_dyn' must be provided")
  }
  
  # Check bivariate specific requirements
  if (setting == "bivariate") {
    if (!is.null(coordx_dyn)) {
      if (!is.list(coordx_dyn) || length(coordx_dyn) != 2) {
        stop("For bivariate with coordx_dyn, coordx_dyn must be a list of length 2")
      }
      if (!is.list(data) || length(data) != 2) {
        stop("For bivariate with coordx_dyn, data must be a list of length 2")
      }
    } else {
      if (!is.matrix(data) || nrow(data) != 2) {
        stop("For bivariate with fixed coordinates, data must be a 2 x N matrix")
      }
    }
  }
}

#' Build coordinate matrix from input vectors
build_coordinate_matrix <- function(coordx, coordy, coordz, grid, coordx_dyn) {
  
  # Skip if using dynamic coordinates
  if (!is.null(coordx_dyn)) {
    return(NULL)
  }
  
  # Handle 1D case
  if (is.null(coordy)) {
    return(as.matrix(coordx))
  }
  
  # Handle 2D/3D cases
  if (grid) {
    # Create grid from coordinate vectors
    if (is.null(coordz)) {
      coords <- as.matrix(expand.grid(coordx, coordy))
    } else {
      coords <- as.matrix(expand.grid(coordx, coordy, coordz))
    }
  } else {
    # Use coordinates as-is
    if (is.null(coordz)) {
      coords <- cbind(coordx, coordy)
    } else {
      coords <- cbind(coordx, coordy, coordz)
    }
  }
  
  return(coords)
}

#' Identify missing/invalid values
identify_missing <- function(x) {
  is.na(x) | is.nan(x) | is.infinite(x)
}

#' Process spatial data
process_spatial <- function(data, coords, X) {
  
  # Find missing values
  missing_mask <- identify_missing(data)
  valid_mask <- !missing_mask
  
  # Calculate removal percentage
  perc <- sum(missing_mask) / length(data)
  
  # Remove missing values if any exist
  if (perc > 0) {
    data <- data[valid_mask]
    coords <- coords[valid_mask, , drop = FALSE]
    if (!is.null(X)) {
      X <- X[valid_mask, , drop = FALSE]
    }
  }
  
  return(list(
    data = data,
    coords = coords,
    coordx_dyn = NULL,
    X = X,
    perc = perc
  ))
}

#' Process spacetime data
process_spacetime <- function(data, coords, coordt, coordx_dyn, X) {
  
  if (!is.null(coordx_dyn)) {
    # Dynamic coordinates case
    return(process_spacetime_dynamic(data, coordx_dyn, X))
  } else {
    # Fixed coordinates case
    return(process_spacetime_fixed(data, coords, X))
  }
}

#' Process spacetime data with dynamic coordinates
process_spacetime_dynamic <- function(data, coordx_dyn, X) {
  
  if (!is.list(data)) {
    stop("For spacetime with dynamic coordinates, data must be a list")
  }
  
  cleaned_data <- list()
  cleaned_coords <- list()
  cleaned_X <- if (!is.null(X)) list() else NULL
  total_removed <- 0
  total_points <- 0
  
  for (t in seq_along(data)) {
    # Find missing values for time t
    missing_mask <- identify_missing(data[[t]])
    valid_mask <- !missing_mask
    
    # Update totals for percentage calculation
    total_removed <- total_removed + sum(missing_mask)
    total_points <- total_points + length(data[[t]])
    
    # Clean data for time t
    cleaned_data[[t]] <- data[[t]][valid_mask]
    cleaned_coords[[t]] <- coordx_dyn[[t]][valid_mask, , drop = FALSE]
    
    if (!is.null(X)) {
      cleaned_X[[t]] <- X[[t]][valid_mask, , drop = FALSE]
    }
  }
  
  perc <- total_removed / total_points
  
  return(list(
    data = cleaned_data,
    coords = NULL,
    coordx_dyn = cleaned_coords,
    X = cleaned_X,
    perc = perc
  ))
}

#' Process spacetime data with fixed coordinates
process_spacetime_fixed <- function(data, coords, X) {
  
  if (!is.matrix(data)) {
    stop("For spacetime with fixed coordinates, data must be a matrix (T x N)")
  }
  
  # Find locations with any missing values across time
  missing_mask <- apply(data, 2, function(col) any(identify_missing(col)))
  valid_mask <- !missing_mask
  
  perc <- sum(missing_mask) / ncol(data)
  
  # Remove locations with missing values
  if (perc > 0) {
    data <- data[, valid_mask, drop = FALSE]
    coords <- coords[valid_mask, , drop = FALSE]
    if (!is.null(X)) {
      X <- X[, valid_mask, drop = FALSE]
    }
  }
  
  return(list(
    data = data,
    coords = coords,
    coordx_dyn = NULL,
    X = X,
    perc = perc
  ))
}

#' Process bivariate data
process_bivariate <- function(data, coords, coordx_dyn, X) {
  
  if (!is.null(coordx_dyn)) {
    # Dynamic coordinates case
    return(process_bivariate_dynamic(data, coordx_dyn, X))
  } else {
    # Fixed coordinates case
    return(process_bivariate_fixed(data, coords, X))
  }
}

#' Process bivariate data with dynamic coordinates
process_bivariate_dynamic <- function(data, coordx_dyn, X) {
  
  cleaned_data <- list()
  cleaned_coords <- list()
  cleaned_X <- if (!is.null(X)) list() else NULL
  perc <- numeric(2)
  
  for (i in 1:2) {
    # Validate dimensions
    if (length(data[[i]]) != nrow(coordx_dyn[[i]])) {
      stop(paste0("Length mismatch between data[[", i, "]] and coordx_dyn[[", i, "]]"))
    }
    
    # Find valid observations
    missing_mask <- identify_missing(data[[i]])
    valid_mask <- !missing_mask
    
    perc[i] <- sum(missing_mask) / length(data[[i]])
    
    # Clean data for variable i
    cleaned_data[[i]] <- data[[i]][valid_mask]
    cleaned_coords[[i]] <- coordx_dyn[[i]][valid_mask, , drop = FALSE]
    
    if (!is.null(X)) {
      if (!is.matrix(X[[i]]) || nrow(X[[i]]) != length(data[[i]])) {
        stop(paste0("X[[", i, "]] dimension mismatch with data[[", i, "]]"))
      }
      cleaned_X[[i]] <- X[[i]][valid_mask, , drop = FALSE]
    }
  }
  
  return(list(
    data = cleaned_data,
    coords = NULL,
    coordx_dyn = cleaned_coords,
    X = cleaned_X,
    perc = perc
  ))
}

#' Process bivariate data with fixed coordinates
process_bivariate_fixed <- function(data, coords, X) {
  
  # Find locations valid for both variables
  valid_mask1 <- !identify_missing(data[1, ])
  valid_mask2 <- !identify_missing(data[2, ])
  valid_mask <- valid_mask1 & valid_mask2
  
  perc <- 1 - sum(valid_mask) / length(valid_mask)
  
  # Remove locations with missing values in either variable
  if (perc > 0) {
    data <- data[, valid_mask, drop = FALSE]
    coords <- coords[valid_mask, , drop = FALSE]
    if (!is.null(X)) {
      X <- X[valid_mask, , drop = FALSE]
    }
  }
  
  return(list(
    data = data,
    coords = coords,
    coordx_dyn = NULL,
    X = X,
    perc = perc
  ))
}

# Example usage and tests
if (FALSE) {
  # Test spatial case
  set.seed(123)
  data_spatial <- c(1, 2, NA, 4, Inf, 6)
  coords_x <- 1:6
  coords_y <- 1:6
  
  result_spatial <- GeoNA(data_spatial, coords_x, coords_y, setting = "spatial")
  cat("Spatial test - Removed percentage:", result_spatial$perc, "\n")
  
  # Test spacetime case with fixed coordinates
  data_st <- matrix(c(1, 2, NA, 4, 5, 6, 7, 8), nrow = 2)
  coords_st <- cbind(1:4, 1:4)
  
  result_st <- GeoNA(data_st, coords_st[,1], coords_st[,2], 
                     coordt = 1:2, setting = "spacetime")
  cat("Spacetime test - Removed percentage:", result_st$perc, "\n")
  
  # Test bivariate case
  data_biv <- rbind(c(1, 2, NA, 4), c(5, 6, 7, Inf))
  coords_biv <- cbind(1:4, 1:4)
  
  result_biv <- GeoNA(data_biv, coords_biv[,1], coords_biv[,2], 
                      setting = "bivariate")
  cat("Bivariate test - Removed percentage:", result_biv$perc, "\n")
}





  # Input validation
  validate_inputs(data, coordx, coordy, coordz, coordt, coordx_dyn, setting)
  
  # Build coordinate matrix for static coordinates
  coords <- build_coordinate_matrix(coordx, coordy, coordz, grid, coordx_dyn)
  
  # Process data based on setting
  result <- switch(setting,
    "spatial" = process_spatial(data, coords, X),
    "spacetime" = process_spacetime(data, coords, coordt, coordx_dyn, X),
    "bivariate" = process_bivariate(data, coords, coordx_dyn, X),
    stop("Invalid setting. Must be 'spatial', 'spacetime', or 'bivariate'")
  )
  
  # Return standardized output
  return(list(
    coordx = if (is.null(coordx_dyn)) result$coords else coordx,
    coordy = coordy,
    coordz = coordz,
    coordt = coordt,
    coordx_dyn = result$coordx_dyn,
    data = result$data,
    grid = grid,
    perc = result$perc,
    setting = setting,
    X = result$X
  ))
}


