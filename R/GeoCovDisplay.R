####################################################
### File name: GeoCovDisplay.R
####################################################

GeoCovDisplay <- function(covmatrix, limits = FALSE, pch = 2) {

  ## Internal display controls.
  tol <- 1e-150
  max_points <- 250000L
  max_raster_cells <- 1000000L
  main <- NULL

  ##########################################################################
  ## Input checks
  ##########################################################################

  if (!inherits(covmatrix, "GeoCovmatrix")) {
    stop("A GeoCovmatrix object is needed as input.")
  }

  if (is.null(covmatrix$covmatrix)) {
    stop("The GeoCovmatrix object does not contain a covariance matrix.")
  }

  if (length(limits) != 1L || is.na(limits)) {
    stop("limits must be TRUE or FALSE.")
  }

  limits <- isTRUE(limits)

  if (length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("tol must be a finite non-negative scalar.")
  }

  if (length(max_points) != 1L ||
      !is.finite(max_points) || max_points < 1) {
    stop("max_points must be a positive finite scalar.")
  }

  if (length(max_raster_cells) != 1L ||
      !is.finite(max_raster_cells) || max_raster_cells < 4) {
    stop("max_raster_cells must be a finite scalar not smaller than 4.")
  }

  max_points <- as.integer(max_points)
  max_raster_cells <- as.integer(max_raster_cells)

  covmat <- covmatrix$covmatrix
  is_spam <- inherits(covmat, "spam")
  is_dense_matrix <- is.matrix(covmat)

  if (!is_spam && !is_dense_matrix) {
    stop(
      "covmatrix$covmatrix must be either a base R matrix or a spam matrix."
    )
  }

  ##########################################################################
  ## Utilities
  ##########################################################################

  make_raster_dimensions <- function(nr, nc, max_cells) {
    ratio <- nc / nr

    nx <- floor(sqrt(max_cells * ratio))
    nx <- max(1L, min(nc, as.integer(nx)))

    ny <- floor(max_cells / nx)
    ny <- max(1L, min(nr, as.integer(ny)))

    ## If clipping ny permits a wider raster without exceeding the budget,
    ## increase nx accordingly.
    nx <- max(1L, min(nc, as.integer(floor(max_cells / ny))))

    c(nx = nx, ny = ny)
  }

  add_to_raster <- function(z, rows, cols, nr, nc) {
    if (!length(rows)) {
      return(z)
    }

    nx <- nrow(z)
    ny <- ncol(z)

    bx <- pmin(
      nx,
      as.integer(floor((cols - 1) * nx / nc) + 1L)
    )

    ## Matrix row 1 is displayed at the top of the graph.
    by_top <- pmin(
      ny,
      as.integer(floor((rows - 1) * ny / nr) + 1L)
    )
    by <- ny - by_top + 1L

    z[cbind(bx, by)] <- 1L
    z
  }

  draw_empty_frame <- function() {
    graphics::plot.new()
    graphics::plot.window(
      xlim = c(0, 1),
      ylim = c(0, 1),
      xaxs = "i",
      yaxs = "i"
    )

    if (!is.null(main)) {
      graphics::title(main = main)
    }
  }

  draw_points <- function(rows, cols, nr, nc) {
    draw_empty_frame()

    if (length(rows)) {
      graphics::points(
        x = (cols - 0.5) / nc,
        y = 1 - (rows - 0.5) / nr,
        pch = pch
      )
    }
  }

  draw_raster <- function(z) {
    nx <- nrow(z)
    ny <- ncol(z)

    ## image() needs at least two coordinates in each direction. Degenerate
    ## one-row or one-column cases are therefore displayed as points.
    if (nx < 2L || ny < 2L) {
      idx <- which(z != 0, arr.ind = TRUE)
      draw_empty_frame()

      if (nrow(idx)) {
        graphics::points(
          x = (idx[, 1L] - 0.5) / nx,
          y = (idx[, 2L] - 0.5) / ny,
          pch = pch
        )
      }
      return(invisible(NULL))
    }

    graphics::image(
      x = (seq_len(nx) - 0.5) / nx,
      y = (seq_len(ny) - 0.5) / ny,
      z = z,
      zlim = c(0, 1),
      col = grDevices::gray.colors(2L, start = 1, end = 0),
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      xaxs = "i",
      yaxs = "i",
      main = main,
      useRaster = TRUE
    )
  }

  draw_limits <- function(nr, nc) {
    if (!limits) {
      return(invisible(NULL))
    }

    horizontal_cuts <- numeric(0L)
    vertical_cuts <- numeric(0L)

    if (isTRUE(covmatrix$bivariate)) {
      ns <- covmatrix$ns

      if (!is.null(ns) &&
          length(ns) >= 2L &&
          all(is.finite(ns)) &&
          all(ns > 0) &&
          sum(ns) == nr &&
          nr == nc) {
        cuts <- cumsum(ns) / sum(ns)
        cuts <- cuts[-length(cuts)]
        horizontal_cuts <- c(horizontal_cuts, cuts)
        vertical_cuts <- c(vertical_cuts, cuts)
      } else {
        ## Backward-compatible fallback when block sizes are unavailable.
        horizontal_cuts <- c(horizontal_cuts, 0.5)
        vertical_cuts <- c(vertical_cuts, 0.5)
      }
    }

    if (isTRUE(covmatrix$spacetime)) {
      nt <- suppressWarnings(as.integer(covmatrix$numtime))
      ns <- covmatrix$ns

      if (!is.null(ns) &&
          length(ns) > 1L &&
          all(is.finite(ns)) &&
          all(ns > 0) &&
          sum(ns) == nr &&
          nr == nc) {
        cuts <- cumsum(ns) / sum(ns)
        cuts <- cuts[-length(cuts)]
      } else if (length(nt) == 1L && is.finite(nt) && nt > 1L) {
        cuts <- seq_len(nt - 1L) / nt
      } else {
        cuts <- numeric(0L)
      }

      horizontal_cuts <- c(horizontal_cuts, cuts)
      vertical_cuts <- c(vertical_cuts, cuts)
    }

    horizontal_cuts <- unique(horizontal_cuts)
    vertical_cuts <- unique(vertical_cuts)

    if (length(horizontal_cuts)) {
      graphics::abline(h = horizontal_cuts)
    }
    if (length(vertical_cuts)) {
      graphics::abline(v = vertical_cuts)
    }

    invisible(NULL)
  }

  ##########################################################################
  ## Extract and display a spam sparsity pattern without densification
  ##########################################################################

  display_spam <- function(x) {
    entries <- methods::slot(x, "entries")
    colindices <- methods::slot(x, "colindices")
    rowpointers <- methods::slot(x, "rowpointers")
    dims <- methods::slot(x, "dimension")

    nr <- as.integer(dims[1L])
    nc <- as.integer(dims[2L])

    if (length(rowpointers) != nr + 1L) {
      stop("Invalid spam matrix: inconsistent row pointers.")
    }

    if (length(entries) != length(colindices)) {
      stop("Invalid spam matrix: entries and column indices differ in length.")
    }

    if (any(!is.finite(entries))) {
      stop("The covariance matrix contains non-finite stored entries.")
    }

    keep <- abs(entries) > tol
    nnz <- sum(keep)
    total_cells <- as.double(nr) * as.double(nc)

    if (as.double(nnz) == total_cells) {
      message("No image plot: the covariance matrix is dense.")
      return(invisible(list(
        nrow = nr,
        ncol = nc,
        nonzero = nnz,
        mode = "dense"
      )))
    }

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)

    ## For a moderate number of stored non-zero entries, plotting their
    ## coordinates directly is exact and requires O(nnz) rather than O(nr*nc)
    ## memory.
    if (nnz <= max_points) {
      row_counts <- diff(rowpointers)
      rows <- rep.int(seq_len(nr), row_counts)

      rows <- rows[keep]
      cols <- colindices[keep]

      draw_points(rows, cols, nr, nc)
      draw_limits(nr, nc)

      return(invisible(list(
        nrow = nr,
        ncol = nc,
        nonzero = nnz,
        mode = "points"
      )))
    }

    ## For very large sparse patterns, aggregate non-zero entries into a
    ## bounded raster. The full nr x nc logical matrix is never allocated.
    raster_dims <- make_raster_dimensions(
      nr,
      nc,
      max_raster_cells
    )

    z <- matrix(
      0L,
      nrow = raster_dims[["nx"]],
      ncol = raster_dims[["ny"]]
    )

    ## Process row blocks to avoid allocating a row-index vector of length nnz.
    target_entries_per_chunk <- 1000000L
    stored_nnz <- length(entries)
    mean_per_row <- max(1, ceiling(stored_nnz / max(nr, 1L)))
    rows_per_chunk <- max(
      1L,
      as.integer(floor(target_entries_per_chunk / mean_per_row))
    )

    row_starts <- seq.int(1L, nr, by = rows_per_chunk)

    for (r1 in row_starts) {
      r2 <- min(nr, r1 + rows_per_chunk - 1L)
      p1 <- rowpointers[r1]
      p2 <- rowpointers[r2 + 1L] - 1L

      if (p1 > p2) {
        next
      }

      local_counts <- diff(rowpointers[r1:(r2 + 1L)])
      local_rows <- rep.int(r1:r2, local_counts)
      local_keep <- keep[p1:p2]

      if (!any(local_keep)) {
        next
      }

      local_rows <- local_rows[local_keep]
      local_cols <- colindices[p1:p2][local_keep]

      z <- add_to_raster(
        z,
        local_rows,
        local_cols,
        nr,
        nc
      )
    }

    draw_raster(z)
    draw_limits(nr, nc)

    invisible(list(
      nrow = nr,
      ncol = nc,
      nonzero = nnz,
      mode = "raster",
      raster_dimension = dim(z)
    ))
  }

  ##########################################################################
  ## Display an ordinary dense matrix
  ##########################################################################

  display_dense_matrix <- function(x) {
    nr <- nrow(x)
    nc <- ncol(x)

    if (any(!is.finite(x))) {
      stop("The covariance matrix contains non-finite values.")
    }

    total_cells <- as.double(nr) * as.double(nc)

    ## Direct image for moderate matrices.
    if (total_cells <= max_raster_cells) {
      nonzero <- abs(x) > tol
      nnz <- sum(nonzero)

      if (as.double(nnz) == total_cells) {
        message("No image plot: the covariance matrix is dense.")
        return(invisible(list(
          nrow = nr,
          ncol = nc,
          nonzero = nnz,
          mode = "dense"
        )))
      }

      old_par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old_par), add = TRUE)

      ## Transpose after reversing rows to preserve the conventional matrix
      ## display with row 1 at the top, as in the original implementation.
      z <- t(nonzero[nr:1L, , drop = FALSE])
      storage.mode(z) <- "integer"

      draw_raster(z)
      draw_limits(nr, nc)

      return(invisible(list(
        nrow = nr,
        ncol = nc,
        nonzero = nnz,
        mode = "direct"
      )))
    }

    ## A large dense input is already resident in memory. To avoid allocating
    ## an additional nr x nc logical matrix, scan it in row blocks and build a
    ## bounded display raster only.
    raster_dims <- make_raster_dimensions(
      nr,
      nc,
      max_raster_cells
    )

    z <- matrix(
      0L,
      nrow = raster_dims[["nx"]],
      ncol = raster_dims[["ny"]]
    )

    nnz <- 0
    rows_per_chunk <- max(1L, min(nr, 1000L))
    row_starts <- seq.int(1L, nr, by = rows_per_chunk)

    for (r1 in row_starts) {
      r2 <- min(nr, r1 + rows_per_chunk - 1L)
      block <- x[r1:r2, , drop = FALSE]
      idx <- which(abs(block) > tol, arr.ind = TRUE)

      if (!nrow(idx)) {
        next
      }

      rows <- idx[, 1L] + r1 - 1L
      cols <- idx[, 2L]
      nnz <- nnz + nrow(idx)

      z <- add_to_raster(z, rows, cols, nr, nc)
    }

    if (as.double(nnz) == total_cells) {
      message("No image plot: the covariance matrix is dense.")
      return(invisible(list(
        nrow = nr,
        ncol = nc,
        nonzero = nnz,
        mode = "dense"
      )))
    }

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)

    draw_raster(z)
    draw_limits(nr, nc)

    invisible(list(
      nrow = nr,
      ncol = nc,
      nonzero = nnz,
      mode = "raster",
      raster_dimension = dim(z)
    ))
  }

  ##########################################################################
  ## Dispatch
  ##########################################################################

  if (is_spam) {
    return(display_spam(covmat))
  }

  display_dense_matrix(covmat)
}