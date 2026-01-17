GeoVariogramDir <- function(data, coordx, coordy = NULL, coordz = NULL,
                            directions = c(0, 45, 90, 135),
                            tolerance = 22.5,
                            numbins = 13,
                            maxdist = NULL,
                            neighb = NULL,
                            distance = "Eucl",
                            subsample = 1) {
  # ---------------- Checks ----------------
  if (!is.numeric(subsample) || length(subsample) != 1L || !is.finite(subsample) ||
      subsample <= 0 || subsample > 1) {
    stop("subsample must be a single numeric value in (0,1].")
  }

  if (!is.numeric(tolerance) || length(tolerance) != 1L || tolerance <= 0 || tolerance > 90)
    stop("tolerance must be in (0,90].")
  if (!is.numeric(numbins) || length(numbins) != 1L || numbins < 2 || numbins %% 1 != 0)
    stop("numbins must be an integer >= 2.")

  # ---------------- Build coords ----------------
  is_matrix <- is.matrix(coordx) || is.data.frame(coordx)

  if (is_matrix) {
    coords <- as.matrix(coordx)
    if (!(ncol(coords) %in% c(2,3))) stop("coordx must have 2 or 3 columns if matrix/data.frame.")
  } else {
    if (is.null(coordy)) stop("If coordx is a vector, coordy must be provided.")
    coords <- if (is.null(coordz)) cbind(coordx, coordy) else cbind(coordx, coordy, coordz)
  }

  n0 <- nrow(coords)
  if (length(data) != n0) stop("length(data) must match number of coordinate rows.")

  # ---------------- EARLY spatial subsample ----------------
  if (subsample < 1) {
    m <- max(2L, floor(subsample * n0))
    m <- min(m, n0)
    idx0 <- sort(sample.int(n0, m, replace = FALSE))

    coords <- coords[idx0, , drop = FALSE]
    data   <- data[idx0]
  }

  # update n after subsample
  n <- length(data)

  # ---------------- Defaults ----------------
  if (is.null(maxdist)) maxdist <- Inf

  # Smart neighb selection based on n and maxdist
  if (is.null(neighb)) {
    if (is.finite(maxdist)) {
      if (ncol(coords) == 2) {
        area_total <- diff(range(coords[,1])) * diff(range(coords[,2]))
        if (!is.finite(area_total) || area_total <= 0) stop("Degenerate coordinate extent after subsample.")
        density <- n / area_total
        area_circle <- pi * maxdist^2
        expected_neighb <- density * area_circle
        neighb <- min(n - 1L, max(100L, as.integer(expected_neighb * 2)))
      } else {
        volume_total <- diff(range(coords[,1])) * diff(range(coords[,2])) * diff(range(coords[,3]))
        if (!is.finite(volume_total) || volume_total <= 0) stop("Degenerate coordinate extent after subsample.")
        density <- n / volume_total
        volume_sphere <- (4/3) * pi * maxdist^3
        expected_neighb <- density * volume_sphere
        neighb <- min(n - 1L, max(100L, as.integer(expected_neighb * 2)))
      }
    } else {
      neighb <- n - 1L
    }
  } else {
    if (!is.numeric(neighb) || length(neighb) != 1L || neighb < 1 || neighb %% 1 != 0)
      stop("neighb must be a positive integer.")
    neighb <- min(as.integer(neighb), n - 1L)
  }

  # Adaptive chunking based on neighb size
  if (neighb < 500) {
    chunk_size <- n
  } else if (neighb < 5000) {
    chunk_size <- max(500L, min(n, as.integer(1000e6 / (neighb * 24))))
  } else {
    chunk_size <- max(100L, min(n, as.integer(500e6 / (neighb * 24))))
  }

  # ---------------- Distance bins ----------------
  if (is.finite(maxdist)) {
    maxdist_bin <- maxdist
  } else {
    if (n > 10000) {
      samp <- sample(n, min(1000, n))
      dists_samp <- as.matrix(dist(coords[samp, ]))
      maxdist_bin <- as.numeric(quantile(dists_samp[upper.tri(dists_samp)], 0.95, na.rm = TRUE))
    } else {
      xr <- range(coords[,1], na.rm = TRUE)
      yr <- range(coords[,2], na.rm = TRUE)
      maxdist_bin <- sqrt(diff(xr)^2 + diff(yr)^2)
    }
  }

  breaks   <- seq(0, maxdist_bin, length.out = numbins + 1)
  widths   <- diff(breaks)
  centers  <- breaks[-1] - widths/2
  numvario <- numbins
  mm       <- c(breaks[1L], breaks[length(breaks)])

  # ---------------- Directional ranges ----------------
  ang_min <- (directions - tolerance) %% 180
  ang_max <- (directions + tolerance) %% 180
  ndir <- length(directions)

  sums_dir <- replicate(ndir, numeric(numvario), simplify = FALSE)
  cnts_dir <- replicate(ndir, integer(numvario), simplify = FALSE)

  # Helper function: C binning
  .bin_pairs_c <- function(dist_sub, i_sub, j_sub) {
    np <- length(dist_sub)
    if (np == 0L) return(list(moms = numeric(numvario), lbins = integer(numvario)))

    zi <- data[i_sub]
    zj <- data[j_sub]

    binsC  <- double(numbins + 1)
    lbinsC <- integer(numvario)
    momsC  <- double(numvario)

    EV <- dotCall64::.C64(
      "Binned_Variogram2new",
      SIGNATURE = c("double","integer","double","double","double","integer","double","integer","double"),
      INTENT    = c("rw","r","r","r","r","rw","rw","r","r"),
      bins  = binsC,
      np    = as.integer(np),
      data1 = as.double(zi),
      data2 = as.double(zj),
      vdist = as.double(dist_sub),
      lbins = lbinsC,
      moms  = momsC,
      nbins = as.integer(numbins + 1),
      mm    = as.double(mm),
      NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0
    )
    list(moms = EV$moms, lbins = EV$lbins)
  }

  # ---------------- Distance type ----------------
  use_geo <- distance %in% c("Geo", "Chor")

  proj_sinusoidal <- function(M) {
    lon <- M[,1] * pi / 180
    lat <- M[,2] * pi / 180
    cbind(lon * cos(lat), lat)
  }

  chord_from_arc <- function(theta) 2 * sin(theta / 2)

  coords_knn <- if (use_geo) proj_sinusoidal(coords) else coords

  # ---------------- Chunk processing ----------------
  n_chunks <- ceiling(n / chunk_size)

  for (chunk_id in seq_len(n_chunks)) {
    start_idx <- (chunk_id - 1) * chunk_size + 1
    end_idx <- min(chunk_id * chunk_size, n)
    chunk_indices <- start_idx:end_idx
    n_chunk <- length(chunk_indices)

    query_coords <- coords_knn[chunk_indices, , drop = FALSE]

    if (is.finite(maxdist)) {
      if (use_geo) {
        radius_proj <- maxdist * 1.1
        knn_res <- nabor::knn(coords_knn, query_coords, k = neighb + 1, radius = radius_proj)
      } else {
        knn_res <- nabor::knn(coords_knn, query_coords, k = neighb + 1, radius = maxdist)
      }
    } else {
      knn_res <- nabor::knn(coords_knn, query_coords, k = neighb + 1)
    }

    idx_mat <- knn_res$nn.idx
    dists_knn <- knn_res$nn.dists

    # Compute true distances for geodesic cases
    if (use_geo) {
      n_query <- nrow(query_coords)
      K <- ncol(idx_mat)
      dist_true <- matrix(NA_real_, n_query, K)

      for (r in seq_len(n_query)) {
        idx <- idx_mat[r, ]
        valid <- idx > 0
        if (!any(valid)) next

        global_idx <- chunk_indices[r]
        neighbor_idx <- idx[valid]

        arc <- fields::rdist.earth.vec(
          x1 = matrix(coords[global_idx, 1:2], nrow = 1),
          x2 = coords[neighbor_idx, 1:2, drop = FALSE],
          miles = FALSE, R = 1
        )

        if (distance == "Chor") {
          dist_true[r, which(valid)] <- chord_from_arc(arc)
        } else {
          dist_true[r, which(valid)] <- arc
        }
      }
      dist_mat <- dist_true
    } else {
      dist_mat <- dists_knn
    }

    keep <- idx_mat > 0
    if (is.finite(maxdist)) keep <- keep & (dist_mat <= maxdist)
    if (!any(keep)) next

    row_indices <- matrix(rep(chunk_indices, times = ncol(idx_mat)),
                          nrow = n_chunk, byrow = FALSE)
    i <- as.integer(row_indices[keep])
    j <- as.integer(idx_mat[keep])
    dist <- as.double(dist_mat[keep])

    # Remove self-pairs
    self_pairs <- i == j
    if (any(self_pairs)) {
      i <- i[!self_pairs]
      j <- j[!self_pairs]
      dist <- dist[!self_pairs]
    }
    if (length(i) == 0) next

    dx <- coords[j, 1] - coords[i, 1]
    dy <- coords[j, 2] - coords[i, 2]
    angle <- (atan2(dy, dx) * 180 / pi) %% 180

    mask_dir_mat <- matrix(FALSE, nrow = length(angle), ncol = ndir)
    for (k in seq_along(directions)) {
      if (ang_min[k] < ang_max[k]) {
        mask_dir_mat[,k] <- (angle >= ang_min[k]) & (angle <= ang_max[k])
      } else {
        mask_dir_mat[,k] <- (angle >= ang_min[k]) | (angle <= ang_max[k])
      }
    }

    n_pairs_chunk <- length(i)
    use_c_chunk <- n_pairs_chunk > 5000

    if (use_c_chunk) {
      for (k in seq_along(directions)) {
        mask <- mask_dir_mat[,k]
        if (any(mask)) {
          bp <- .bin_pairs_c(dist[mask], i[mask], j[mask])
          sums_dir[[k]] <- sums_dir[[k]] + bp$moms
          cnts_dir[[k]] <- cnts_dir[[k]] + bp$lbins
        }
      }
    } else {
      diff2 <- (data[j] - data[i])^2 / 2
      bin_idx_all <- cut(dist, breaks, labels = FALSE, include.lowest = TRUE)

      for (k in seq_along(directions)) {
        mask <- mask_dir_mat[,k]
        if (any(mask)) {
          for (b in seq_len(numvario)) {
            bin_mask <- !is.na(bin_idx_all[mask]) & bin_idx_all[mask] == b
            if (any(bin_mask)) {
              vals <- diff2[mask][bin_mask]
              sums_dir[[k]][b] <- sums_dir[[k]][b] + sum(vals, na.rm = TRUE)
              cnts_dir[[k]][b] <- cnts_dir[[k]][b] + length(vals)
            }
          }
        }
      }
    }
  }

  # ---------------- Results ----------------
  result <- vector("list", ndir)
  names(result) <- paste0("dir_", directions, "deg")

  for (k in seq_len(ndir)) {
    cnt <- cnts_dir[[k]]
    gamma <- ifelse(cnt > 0, sums_dir[[k]] / cnt, NA_real_)
    result[[k]] <- list(
      centers = centers,
      gamma = as.numeric(gamma),
      npairs = as.numeric(cnt)
    )
  }

  attr(result, "breaks") <- breaks
  attr(result, "maxdist") <- maxdist_bin
  attr(result, "neighb") <- neighb
  attr(result, "n_chunks") <- n_chunks
  attr(result, "subsample") <- subsample
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
