GeoNA <- function(data, coordx = NULL, coordy = NULL, coordz = NULL, coordt = NULL, 
                  coordx_dyn = NULL, grid = FALSE, X = NULL, setting = "spatial") {
  
  if (is.null(coordx_dyn)) {
    if (is.null(coordy)) {
      coords <- as.matrix(coordx)
    } else {
      if (grid) {
        if (is.null(coordz)) {
          coords <- as.matrix(expand.grid(coordx, coordy))
        } else {
          coords <- as.matrix(expand.grid(coordx, coordy, coordz))
        }
      } else {
        if (is.null(coordz)) {
          coords <- cbind(coordx, coordy)
        } else {
          coords <- cbind(coordx, coordy, coordz)
        }
      }
    }
  }

  perc <- NULL

  ### SPATIAL case
  if (setting == "spatial") {
    nasel <- is.nan(data) | is.infinite(data) | is.na(data)
    perc <- sum(nasel) / length(data)
    sel <- !nasel
    
    if (perc > 0) {
      data <- data[sel]
      coords <- coords[sel, , drop = FALSE]
      if (!is.null(X)) X <- X[sel, , drop = FALSE]
    }
  }

  ### SPACETIME case
  if (setting == "spacetime") {

    if (is.null(coordt)) {
      stop("In 'spacetime' setting, 'coordt' cannot be NULL. It represents temporal coordinates.")
    }

    if (!is.null(coordx_dyn)) {
      # Dynamic coordinates: data is a list of matrices (one for each time)
      new_data <- list()
      for (t in seq_along(data)) {
        nasel <- is.na(data[[t]]) | is.nan(data[[t]]) | is.infinite(data[[t]])
        perc <- sum(nasel) / length(data[[t]])
        sel <- !nasel

        new_data[[t]] <- data[[t]][sel]
        coordx_dyn[[t]] <- coordx_dyn[[t]][sel, , drop = FALSE]
        if (!is.null(X)) X[[t]] <- X[[t]][sel, , drop = FALSE]
      }
      data <- new_data

    } else {
      # Fixed coordinates: data is a T x N matrix, coords is N x d
      if (!is.matrix(data)) stop("In 'spacetime' with fixed coordinates, 'data' must be a matrix with time in rows.")
      nasel <- apply(data, 2, function(col) any(is.na(col) | is.nan(col) | is.infinite(col)))
      perc <- sum(nasel) / ncol(data)
      sel <- !nasel

      data <- data[, sel, drop = FALSE]
      coords <- coords[sel, , drop = FALSE]
      if (!is.null(X)) X <- X[, sel, drop = FALSE]
    }
  }

  ### BIVARIATE case
  if (setting == "bivariate") {
    if (!is.null(coordx_dyn)) {
      # Dynamic coordinates: data, coordx_dyn, X are lists of length two
      if (!is.list(data) || length(data) != 2 ||
          !is.list(coordx_dyn) || length(coordx_dyn) != 2) {
        stop("In 'bivariate' with coordx_dyn, 'data' and 'coordx_dyn' must be lists of length 2.")
      }
      if (!is.null(X) && (!is.list(X) || length(X) != 2)) {
        stop("If coordx_dyn is not NULL, X must be a list of two matrices.")
      }

      perc <- numeric(2)

      for (i in 1:2) {
        d <- data[[i]]
        c <- coordx_dyn[[i]]
        if (length(d) != nrow(c)) {
          stop(paste0("Length of data[[", i, "]] and number of rows in coordx_dyn[[", i, "]] do not match."))
        }
        na_mask <- !(is.na(d) | is.nan(d) | is.infinite(d))
        perc[i] <- sum(!na_mask) / length(d)
        data[[i]] <- d[na_mask]
        coordx_dyn[[i]] <- c[na_mask, , drop = FALSE]
        if (!is.null(X)) {
          Xi <- X[[i]]
          if (!is.matrix(Xi) || nrow(Xi) != length(d)) {
            stop(paste0("X[[", i, "]] must be a matrix with the same number of rows as data[[", i, "]]."))
          }
          X[[i]] <- Xi[na_mask, , drop = FALSE]
        }
      }

    } else {
      # Fixed coordinates: data is a 2 x N matrix
      if (!is.matrix(data) || nrow(data) != 2) {
        stop("In 'bivariate' mode without coordx_dyn, 'data' must be a 2 x N matrix.")
      }

      na_mask1 <- !(is.na(data[1, ]) | is.nan(data[1, ]) | is.infinite(data[1, ]))
      na_mask2 <- !(is.na(data[2, ]) | is.nan(data[2, ]) | is.infinite(data[2, ]))
      sel <- na_mask1 & na_mask2
      perc <- 1 - sum(sel) / length(sel)

      data <- data[, sel, drop = FALSE]
      coords <- coords[sel, , drop = FALSE]
      if (!is.null(X)) X <- X[sel, , drop = FALSE]
    }
  }

  return(list(
    coordx = if (is.null(coordx_dyn)) coords else coordx,
    coordy = coordy,
    coordz = coordz,
    coordt = coordt,
    coordx_dyn = coordx_dyn,
    data = data,
    grid = grid,
    perc = perc,
    setting = setting,
    X = X
  ))
}
