####################################################
### File name: GeoSpoutlier.r
####################################################

GeoOutlier <- function(data, coordx, coordy = NULL, coordz = NULL,
                       coordt = NULL, coordx_dyn = NULL,
                       distance = "Eucl", grid = FALSE, neighb = 10,
                       alpha = 0.001, method = "Z-Median", radius = 1,
                       bivariate = FALSE, X = NULL) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1", call. = FALSE)
  }
  if (!is.logical(bivariate) || length(bivariate) != 1L || is.na(bivariate))
    stop("bivariate must be logical", call. = FALSE)
  if (!is.numeric(neighb) || length(neighb) != 1L || !is.finite(neighb) ||
      neighb < 2 || neighb != as.integer(neighb)) {
    stop("neighb must be an integer >= 2", call. = FALSE)
  }
  neighb <- as.integer(neighb)
  if (!identical(method, "Z-Median")) stop("method is not valid", call. = FALSE)
  if ((is.numeric(coordt) && length(coordt) > 1L) || !is.null(coordx_dyn)) {
    stop("GeoOutlier is currently implemented only for purely spatial data.",
         call. = FALSE)
  }
  if (isTRUE(bivariate)) {
    stop("GeoOutlier is currently implemented only for univariate data.",
         call. = FALSE)
  }

  if ((is.matrix(coordx) || is.data.frame(coordx)) &&
      is.null(coordy) && is.null(coordz)) {
    coords <- as.matrix(coordx)
  } else {
    if (is.null(coordy)) {
      stop("coordy is required when coordx is not a coordinate matrix.",
           call. = FALSE)
    }
    coords <- if (is.null(coordz)) {
      if (isTRUE(grid)) as.matrix(expand.grid(coordx, coordy))
      else cbind(coordx, coordy)
    } else {
      if (isTRUE(grid)) as.matrix(expand.grid(coordx, coordy, coordz))
      else cbind(coordx, coordy, coordz)
    }
  }
  coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  .GeoValidate3D(coords, distance = distance, grid = grid,
                 context = "GeoOutlier")

  data <- as.numeric(c(t(data)))
  if (length(data) != nrow(coords)) {
    stop("The number of observations must match the number of coordinate rows.",
         call. = FALSE)
  }
  if (neighb > nrow(coords)) {
    stop("neighb cannot exceed the number of observations.", call. = FALSE)
  }
  if (!is.null(X) && nrow(as.matrix(X)) != length(data)) {
    stop("X must have one row per observation.", call. = FALSE)
  }

  distance <- gsub("[[:blank:]]", "", distance)
  coords_original <- coords
  if (distance %in% c("Geod", "Chor")) {
    if (!requireNamespace("mapproj", quietly = TRUE)) {
      stop("Package 'mapproj' is required by GeoOutlier() for geodesic or chordal distance.",
           call. = FALSE)
    }
    prj <- mapproj::mapproject(coords[, 1L], coords[, 2L],
                              projection = "sinusoidal")
    coords <- radius * cbind(prj$x, prj$y)
  } else if (!identical(distance, "Eucl")) {
    stop("distance must be one of 'Eucl', 'Geod', or 'Chor'.", call. = FALSE)
  }

  aa <- nabor::knn(coords, coords, k = neighb)
  local_values <- lapply(seq_len(nrow(coords)), function(i) {
    data[aa$nn.idx[i, -1L]]
  })
  res <- matrix(unlist(local_values, use.names = FALSE),
                ncol = neighb - 1L, byrow = TRUE)
  local_median <- apply(res, 1L, stats::median)
  deviation <- data - local_median
  sd_deviation <- stats::sd(deviation)
  if (!is.finite(sd_deviation) || sd_deviation <= 0) {
    selected <- rep(FALSE, length(deviation))
  } else {
    zscore <- (deviation - mean(deviation)) / sd_deviation
    selected <- zscore > stats::qnorm(1 - alpha / 2) |
      zscore < stats::qnorm(alpha / 2)
  }

  output_coords <- if (distance %in% c("Geod", "Chor")) coords_original else coords
  if (is.null(X)) {
    cbind(output_coords, data = data)[selected, , drop = FALSE]
  } else {
    cbind(output_coords, data = data, as.matrix(X))[selected, , drop = FALSE]
  }
}
