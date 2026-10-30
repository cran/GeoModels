## Geographic nearest-neighbor searches use a unit-sphere embedding.  Chord
## distance is monotone in great-circle distance on [0, pi], including across
## the antimeridian and near the poles.
.GeoVariogramDirSphere <- function(coords) {
  lon <- coords[, 1L] * pi / 180
  lat <- coords[, 2L] * pi / 180
  cbind(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
}

.GeoVariogramDirSphereAngle <- function(xyz, i, j) {
  midpoint <- xyz[i, , drop = FALSE] + xyz[j, , drop = FALSE]
  tangent <- xyz[j, , drop = FALSE] - xyz[i, , drop = FALSE]
  norm_mid <- sqrt(rowSums(midpoint^2))
  horizontal <- sqrt(rowSums(midpoint[, 1:2, drop = FALSE]^2))
  lon <- atan2(midpoint[, 2L], midpoint[, 1L])
  lat <- atan2(midpoint[, 3L], horizontal)
  east <- cbind(-sin(lon), cos(lon), 0)
  north <- cbind(-sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat))
  angle <- (atan2(rowSums(tangent * north), rowSums(tangent * east)) * 180 / pi) %% 180
  ## Antipodal endpoints, coincident endpoints and a polar midpoint do not
  ## define this local geographic azimuth.  They are included only when the
  ## tolerance is 90 degrees (all directions).
  undefined <- norm_mid < 64 * .Machine$double.eps |
    horizontal < 64 * .Machine$double.eps |
    rowSums(tangent^2) < (64 * .Machine$double.eps)^2
  angle[undefined] <- NA_real_
  angle
}

GeoVariogramDir <- function(data, coordx, coordy = NULL, coordz = NULL,
                            directions = c(0, 45, 90, 135),
                            tolerance = 22.5, numbins = 13,
                            maxdist = NULL, neighb = NULL,
                            distance = "Eucl", subsample = 1) {
  if (!is.numeric(subsample) || length(subsample) != 1L || !is.finite(subsample) ||
      subsample <= 0 || subsample > 1)
    stop("subsample must be a single numeric value in (0,1].")
  if (!is.numeric(tolerance) || length(tolerance) != 1L || !is.finite(tolerance) ||
      tolerance <= 0 || tolerance > 90)
    stop("tolerance must be in (0,90].")
  if (!is.numeric(numbins) || length(numbins) != 1L || !is.finite(numbins) ||
      numbins < 2 || numbins %% 1 != 0 || numbins >= .Machine$integer.max)
    stop("numbins must be an integer >= 2.")
  if (!is.numeric(directions) || !length(directions) || any(!is.finite(directions)))
    stop("directions must be a non-empty finite numeric vector.")

  if (is.matrix(coordx) || is.data.frame(coordx)) {
    coords <- as.matrix(coordx)
  } else {
    if (is.null(coordy)) stop("If coordx is a vector, coordy must be provided.")
    if (length(coordx) != length(coordy) ||
        (!is.null(coordz) && length(coordz) != length(coordx)))
      stop("Coordinate vectors must have equal lengths.")
    coords <- if (is.null(coordz)) cbind(coordx, coordy) else cbind(coordx, coordy, coordz)
  }
  if (!is.numeric(coords) || !(ncol(coords) %in% c(2L, 3L)) ||
      nrow(coords) < 2L || any(!is.finite(coords)))
    stop("Coordinates must be a finite numeric matrix with at least two rows and two or three columns.")
  storage.mode(coords) <- "double"
  distance <- match.arg(distance, c("Eucl", "Chor", "Geod"))
  .GeoValidate3D(coords, distance = distance, context = "GeoVariogramDir")
  use_geo <- distance %in% c("Geod", "Chor")
  if (use_geo && any(abs(coords[, 2L]) > 90))
    stop("Latitude must lie between -90 and 90 degrees.")
  if (!is.numeric(data) || length(data) != nrow(coords))
    stop("data must be numeric with one value per coordinate row.")

  if (subsample < 1) {
    m <- min(nrow(coords), max(2L, floor(subsample * nrow(coords))))
    idx0 <- sort(sample.int(nrow(coords), m, replace = FALSE))
    coords <- coords[idx0, , drop = FALSE]
    data <- data[idx0]
  }
  n <- length(data)
  if (is.null(maxdist)) maxdist <- Inf
  if (!is.numeric(maxdist) || length(maxdist) != 1L || is.na(maxdist) || maxdist <= 0)
    stop("maxdist must be a positive number or Inf.")

  coords_knn <- if (use_geo) .GeoVariogramDirSphere(coords) else coords
  metric_distance <- function(chord) {
    if (identical(distance, "Geod")) 2 * asin(pmin(1, pmax(0, chord / 2))) else chord
  }

  if (is.null(neighb)) {
    if (!is.finite(maxdist)) {
      neighb <- n - 1L
    } else {
      if (use_geo) {
        lon <- sort((coords[, 1L] * pi / 180) %% (2 * pi))
        lon_width <- 2 * pi - max(diff(c(lon, lon[1L] + 2 * pi)))
        area <- lon_width * diff(range(sin(coords[, 2L] * pi / 180)))
        theta <- if (distance == "Geod") min(pi, maxdist) else
          2 * asin(min(1, maxdist / 2))
        neighborhood_volume <- 4 * pi * sin(theta / 2)^2
        extent <- area
      } else {
        spans <- apply(coords, 2L, function(z) diff(range(z)))
        extent <- prod(spans)
        neighborhood_volume <- if (ncol(coords) == 2L) pi * maxdist^2 else
          (4 / 3) * pi * maxdist^3
      }
      ## This only chooses the automatic neighbor cap, as in the historical
      ## API.  Degenerate bounding rectangles use all available neighbors.
      expected <- if (is.finite(extent) && extent > 0)
        2 * n * neighborhood_volume / extent else n - 1L
      neighb <- as.integer(min(n - 1L, max(100, expected)))
    }
  } else {
    if (!is.numeric(neighb) || length(neighb) != 1L || !is.finite(neighb) ||
        neighb < 1 || neighb %% 1 != 0)
      stop("neighb must be a positive integer.")
    neighb <- as.integer(min(neighb, n - 1L))
  }

  ## Bound temporary neighbor matrices independently of the number of sites.
  chunk_size <- max(1L, min(n, floor(64 * 1024^2 / (24 * (neighb + 1)))))
  if (is.finite(maxdist)) {
    maxdist_bin <- maxdist
  } else if (n > 10000L) {
    samp <- sample.int(n, min(1000L, n))
    sampled <- as.numeric(dist(coords_knn[samp, , drop = FALSE]))
    if (use_geo) sampled <- metric_distance(sampled)
    maxdist_bin <- as.numeric(quantile(sampled, 0.95, na.rm = TRUE))
  } else {
    spans <- apply(coords_knn, 2L, function(z) diff(range(z)))
    bound <- sqrt(sum(spans^2))
    maxdist_bin <- if (use_geo) metric_distance(min(2, bound)) else bound
  }
  if (!is.finite(maxdist_bin) || maxdist_bin <= 0)
    stop("No positive distance range is available for variogram bins.")

  breaks <- seq(0, maxdist_bin, length.out = numbins + 1L)
  centers <- breaks[-1L] - diff(breaks) / 2
  mm <- c(0, maxdist_bin)
  ndir <- length(directions)
  ang_min <- (directions - tolerance) %% 180
  ang_max <- (directions + tolerance) %% 180
  sums_dir <- replicate(ndir, numeric(numbins), simplify = FALSE)
  cnts_dir <- replicate(ndir, numeric(numbins), simplify = FALSE)

  radius_knn <- if (!is.finite(maxdist)) Inf else if (distance == "Geod")
    2 * sin(min(pi, maxdist) / 2) else if (distance == "Chor") min(2, maxdist) else maxdist
  ## Radius searches use a strict floating-point boundary in some backends.
  if (is.finite(radius_knn))
    radius_knn <- radius_knn * (1 + 64 * .Machine$double.eps)
  n_chunks <- ceiling(n / chunk_size)
  for (chunk_id in seq_len(n_chunks)) {
    chunk_indices <- seq.int((chunk_id - 1L) * chunk_size + 1L,
                             min(chunk_id * chunk_size, n))
    knn_args <- list(data = coords_knn, query = coords_knn[chunk_indices, , drop = FALSE],
                     k = neighb + 1L)
    if (is.finite(radius_knn)) knn_args$radius <- radius_knn
    knn_res <- do.call(nabor::knn, knn_args)
    idx_mat <- knn_res$nn.idx
    dist_mat <- if (use_geo) metric_distance(knn_res$nn.dists) else knn_res$nn.dists
    ## pmin/pmax drop matrix attributes on some R versions.
    dim(dist_mat) <- dim(idx_mat)
    keep <- idx_mat > 0 & is.finite(dist_mat)
    if (is.finite(maxdist)) keep <- keep & dist_mat <= maxdist * (1 + 64 * .Machine$double.eps)
    if (!any(keep)) next
    i <- rep(chunk_indices, times = ncol(idx_mat))[keep]
    j <- idx_mat[keep]
    lag <- dist_mat[keep]
    valid <- i != j & is.finite(data[i]) & is.finite(data[j])
    i <- i[valid]; j <- j[valid]; lag <- lag[valid]
    if (!length(i)) next
    if (is.finite(maxdist)) lag <- pmin(lag, maxdist)
    angle <- if (use_geo) .GeoVariogramDirSphereAngle(coords_knn, i, j) else
      (atan2(coords[j, 2L] - coords[i, 2L], coords[j, 1L] - coords[i, 1L]) * 180 / pi) %% 180

    ## One binning implementation is used for every chunk size.  This avoids
    ## different edge or missing-data behavior in R and C branches.
    for (k in seq_along(directions)) {
      mask <- if (tolerance == 90) rep(TRUE, length(lag)) else if (ang_min[k] < ang_max[k])
        is.finite(angle) & angle >= ang_min[k] & angle <= ang_max[k] else
        is.finite(angle) & (angle >= ang_min[k] | angle <= ang_max[k])
      if (!any(mask)) next
      ans <- dotCall64::.C64(
        "Binned_Variogram2new",
        SIGNATURE = c("double", "integer", "double", "double", "double", "integer", "double", "integer", "double"),
        INTENT = c("w", "r", "r", "r", "r", "rw", "rw", "r", "r"),
        bins = double(numbins + 1L), np = as.integer(sum(mask)),
        data1 = as.double(data[i[mask]]), data2 = as.double(data[j[mask]]),
        vdist = as.double(lag[mask]), lbins = integer(numbins), moms = double(numbins),
        nbins = as.integer(numbins + 1L), mm = mm,
        NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0)
      sums_dir[[k]] <- sums_dir[[k]] + ans$moms
      cnts_dir[[k]] <- cnts_dir[[k]] + ans$lbins
    }
  }

  result <- lapply(seq_along(directions), function(k) {
    cnt <- cnts_dir[[k]]
    list(centers = centers, gamma = as.numeric(ifelse(cnt > 0, sums_dir[[k]] / cnt, NA_real_)),
         npairs = as.numeric(cnt))
  })
  names(result) <- paste0("dir_", directions, "deg")
  attr(result, "breaks") <- breaks
  attr(result, "maxdist") <- maxdist_bin
  attr(result, "neighb") <- neighb
  attr(result, "n_chunks") <- n_chunks
  attr(result, "subsample") <- subsample
  attr(result, "direction_definition") <- if (use_geo) {
    "Unoriented great-circle tangent at the pair midpoint, measured from local east toward north."
  } else if (ncol(coords) == 3L) {
    "Azimuth in the xy plane; Euclidean lag distances use x, y, and z."
  } else "Azimuth in the xy plane."
  structure(result, class = "GeoVariogramDir")
}

################################################################################
plot.GeoVariogramDir <- function(x, ...,
                                 main = "Directional Empirical Semivariograms",
                                 pch = 20,
                                 lwd = 1,
                                 col = 1:8,
                                 ylab = "Semivariogram",
                                 xlab = "Distance") {
  directions <- names(x)
  ndir <- length(directions)
  col <- rep(col, length.out = ndir)

  dots <- list(...)
  legend_cex   <- if (!is.null(dots$legend.cex))   dots$legend.cex   else 0.8
  legend_inset <- if (!is.null(dots$legend.inset)) dots$legend.inset else 0.02
  dots$legend.cex <- dots$legend.inset <- NULL

  # limiti default
  all_centers <- unlist(lapply(x, function(vg) vg$centers))
  all_gamma   <- unlist(lapply(x, function(vg) vg$gamma))
  xlim_default <- range(all_centers, na.rm = TRUE)
  ylim_default <- range(all_gamma,   na.rm = TRUE)
  if (!is.finite(ylim_default[1]) || !is.finite(ylim_default[2])) ylim_default <- c(0, 1)

  xlim <- if (!is.null(dots$xlim)) dots$xlim else xlim_default
  ylim <- if (!is.null(dots$ylim)) dots$ylim else ylim_default
  main <- if (!is.null(dots$main)) dots$main else main
  xlab <- if (!is.null(dots$xlab)) dots$xlab else xlab
  ylab <- if (!is.null(dots$ylab)) dots$ylab else ylab
  dots$xlim <- dots$ylim <- dots$main <- dots$xlab <- dots$ylab <- NULL

  # setup margini + xpd per disegnare fuori
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  par(mar = par("mar") + c(0, 0, 0, 4), xpd = NA)  # spazio extra a destra

  do.call(plot, c(list(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main), dots))
  for (i in seq_along(directions)) {
    vg <- x[[i]]
    xvals <- vg$centers
    yvals <- vg$gamma
    if (length(yvals) == 0 || all(is.na(yvals))) next
    lines(xvals, yvals, type = "b", lwd = lwd, pch = pch, col = col[i])
  }

  degs <- as.integer(gsub("dir_(\\d+)deg", "\\1", directions))
  legend_labels <- lapply(degs, function(d) bquote(.(d)*degree))

  usr <- par("usr")
  x_off <- usr[2] + legend_inset * diff(usr[1:2])  # un po' oltre il bordo destro
  y_top <- usr[4]

  legend(x = x_off, y = y_top,
         legend = legend_labels, col = col, lwd = lwd, pch = pch,
         xjust = 0, yjust = 1, bty = "n", cex = legend_cex)
}
