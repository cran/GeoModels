GeoVariogramDir <- function(data, coordx, coordy = NULL, coordz = NULL,
                            directions = c(0, 45, 90, 135),
                            tolerance = 22.5,
                            numbins = 13,
                            maxdist = NULL,
                            neighb = NULL,
                            distance = "Eucl") {
  # Check input and set defaults
  n <- length(data)
  is_matrix <- is.matrix(coordx) || is.data.frame(coordx)
  if (is_matrix) {
    if (!(ncol(coordx) %in% c(2,3))) stop("coordx must have 2 or 3 columns if matrix/data.frame.")
  } else {
    if (is.null(coordy)) stop("If coordx is a vector, coordy must be provided.")
  }
  if (is.null(maxdist)) maxdist <- Inf
  if (is.null(neighb)) neighb <- n - 1

  # Find neighbor pairs using GeoNeighIndex
  neigh <- GeoNeighIndex(coordx, coordy, coordz, distance = distance, neighb = neighb, maxdist = maxdist)
  i <- neigh$rowidx
  j <- neigh$colidx
  dist <- neigh$lags

  # Compute coordinate differences for angle calculation
  if (is_matrix) {
    dx <- coordx[j,1] - coordx[i,1]
    dy <- coordx[j,2] - coordx[i,2]
  } else {
    dx <- coordx[j] - coordx[i]
    dy <- coordy[j] - coordy[i]
  }
  angle <- (atan2(dy, dx) * 180 / pi) %% 180
  diff2 <- (data[j] - data[i])^2 / 2

  # Precompute distance bins
  maxdist_bin <- if (is.infinite(maxdist)) max(dist) else maxdist
  breaks <- seq(0, maxdist_bin, length.out = numbins + 1)
  centers <- breaks[-1] - diff(breaks)/2
  bin_idx_all <- cut(dist, breaks, labels = FALSE, include.lowest = TRUE)

  # Precompute directional masks for all pairs and directions
  ang_min <- (directions - tolerance) %% 180
  ang_max <- (directions + tolerance) %% 180
  ndir <- length(directions)
  mask_dir_mat <- matrix(FALSE, nrow = length(angle), ncol = ndir)
  for (k in seq_along(directions)) {
    if (ang_min[k] < ang_max[k]) {
      mask_dir_mat[,k] <- (angle >= ang_min[k]) & (angle <= ang_max[k])
    } else {
      mask_dir_mat[,k] <- (angle >= ang_min[k]) | (angle <= ang_max[k])
    }
  }

  # Compute directional empirical semivariograms
  result <- vector("list", ndir)
  names(result) <- paste0("dir_", directions, "deg")
  for (k in seq_along(directions)) {
    mask <- mask_dir_mat[,k]
    if(any(mask)) {
      gamma <- tapply(diff2[mask], bin_idx_all[mask], mean)
      npairs <- tapply(diff2[mask], bin_idx_all[mask], length)
    } else {
      gamma <- npairs <- numeric(0)
    }
    result[[k]] <- list(
      centers = centers,
      gamma = as.numeric(gamma),
      npairs = as.numeric(npairs)
    )
  }
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
  
  # Raccogli tutti gli argomenti grafici passati in ...
  dots <- list(...)
  
  # Calcola i limiti di default per x e y
  all_centers <- unlist(lapply(x, function(vg) vg$centers))
  all_gamma <- unlist(lapply(x, function(vg) vg$gamma))
  xlim_default <- range(all_centers, na.rm = TRUE)
  ylim_default <- range(all_gamma, na.rm = TRUE)
  if (!is.finite(ylim_default[1]) || !is.finite(ylim_default[2])) ylim_default <- c(0, 1)
  
  # Usa i limiti passati dall'utente se presenti, altrimenti quelli di default
  xlim <- if (!is.null(dots$xlim)) dots$xlim else xlim_default
  ylim <- if (!is.null(dots$ylim)) dots$ylim else ylim_default
  main <- if (!is.null(dots$main)) dots$main else main
  xlab <- if (!is.null(dots$xlab)) dots$xlab else xlab
  ylab <- if (!is.null(dots$ylab)) dots$ylab else ylab
  
  # Rimuovi questi argomenti da dots per non passarli due volte a plot()
  dots$xlim <- dots$ylim <- dots$main <- dots$xlab <- dots$ylab <- NULL
  
  # Disegna il plot vuoto con tutti i parametri (default o utente)
  do.call(plot, c(list(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main), dots))
  
  for (i in seq_along(directions)) {
    vg <- x[[i]]
    xvals <- vg$centers
    yvals <- vg$gamma
    if (length(yvals) == 0 || all(is.na(yvals))) next
    lines(xvals, yvals, type = "b", lwd = lwd, pch = pch, col = col[i])
  }
  
  # Legenda con simbolo di grado (ASCII safe)
  degs <- as.integer(gsub("dir_(\\d+)deg", "\\1", directions))
  legend_labels <- lapply(degs, function(d) bquote(.(d)*degree))
  legend("topleft", legend = legend_labels, col = col, lwd = lwd, pch = pch, bty = "n")
}
