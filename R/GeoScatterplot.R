GeoScatterplot <- function(data, coordx, coordy = NULL, coordz = NULL, coordt = NULL,
                           coordx_dyn = NULL, distance = "Eucl", grid = FALSE,
                           maxdist = NULL, neighb = NULL,
                           times = NULL, numbins = 4, radius = 1,
                           bivariate = FALSE, ...) {

    call <- match.call()
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))

    type_dist <- CheckDistance(distance)
    spacetime <- FALSE
    maxtime <- NULL
    maxpairs <- Inf

    if (length(numbins) != 1 || is.na(numbins) || numbins < 1 ||
        numbins != as.integer(numbins)) {
        stop("numbins must be a positive integer")
    }
    numbins <- as.integer(numbins)


    if (bivariate) coordt <- c(0, 1)

    if (!is.null(coordt) && is.numeric(coordt) && !is.null(times) && is.numeric(times)) {
        if (length(coordt) > 1 && length(times) >= 1) spacetime <- TRUE
    }

    if (spacetime) {
        if (length(times) > 3) stop("Insert no more than three temporal instants")
        if (!is.null(maxtime) && maxtime < min(stats::dist(coordt))) {
            stop("Maximum temporal distance is too small")
        }
    }

    if (grid) {
        if (is.null(coordy)) stop("coordy must be provided when grid=TRUE")
        if (is.null(coordz)) {
            a <- as.matrix(expand.grid(coordx, coordy))
            coordx <- a[, 1]
            coordy <- a[, 2]
            coordz <- rep(0, length(coordx))
        } else {
            a <- as.matrix(expand.grid(coordx, coordy, coordz))
            coordx <- a[, 1]
            coordy <- a[, 2]
            coordz <- a[, 3]
        }
    }

    if (is.null(coordx_dyn)) {
        if (is.matrix(coordx) || is.data.frame(coordx)) {
            coordx <- as.matrix(coordx)
            if (ncol(coordx) == 1) {
                if (is.null(coordy)) stop("coordy must be provided when coordx has one column")
                if (is.null(coordz)) coordz <- rep(0, nrow(coordx))
                coords <- cbind(coordx[, 1], coordy, coordz)
            } else if (ncol(coordx) == 2) {
                coords <- cbind(coordx[, 1], coordx[, 2], 0)
            } else {
                coords <- coordx[, 1:3, drop = FALSE]
            }
        } else {
            if (is.null(coordy)) {
                stop("coordy must be provided when coordx is not a coordinate matrix")
            }
            if (is.null(coordz)) coordz <- rep(0, length(coordx))
            coords <- cbind(coordx, coordy, coordz)
        }

        coords <- as.matrix(coords)
        storage.mode(coords) <- "double"

        if (bivariate) coords <- rbind(coords, coords)

        data <- as.numeric(t(data))
        numcoord <- nrow(coords)

        if (bivariate) {
            ns <- rep(numcoord / 2, 2)
        } else if (!is.null(coordt)) {
            ns <- rep(numcoord, length(coordt))
        } else {
            ns <- numcoord
        }
    } else {
        coords <- do.call(rbind, coordx_dyn)
        coords <- as.matrix(coords)
        storage.mode(coords) <- "double"
        if (ncol(coords) == 2) coords <- cbind(coords, 0)

        data <- as.numeric(unlist(data))
        ns <- lengths(coordx_dyn) / 2
    }

    if (ncol(coords) == 2) coords <- cbind(coords, 0)
    if (ncol(coords) != 3) stop("coordinates must have two or three columns")

    if (length(data) != nrow(coords)) {
        stop("length(data) must match the number of coordinate rows after processing")
    }

    .panel_layout <- function(n) {
        nc <- ceiling(sqrt(n))
        nr <- ceiling(n / nc)
        c(nr, nc)
    }

    .check_pairs_memory <- function(numcoords) {
        n_pairs <- numcoords * (numcoords - 1) / 2
        if (is.finite(maxpairs) && n_pairs > maxpairs) {
            stop("Too many pairs requested. Reduce maxdist, use neighb, or increase maxpairs.")
        }
        as.integer(n_pairs)
    }

    .extract_pairs <- function(V) {
        v0 <- V[[8]]
        v1 <- V[[9]]
        v2 <- V[[10]]
        keep <- v0 != -1
        data.frame(v0 = v0[keep], v1 = v1[keep], v2 = v2[keep])
    }

    .call_pairs <- function(data, coords, bins, maxdist) {
        numcoords <- nrow(coords)
        n_pairs <- .check_pairs_memory(numcoords)

        v0 <- rep(-1, n_pairs)
        v1 <- rep(NA_real_, n_pairs)
        v2 <- rep(NA_real_, n_pairs)

        V <- .C("pairs",
                as.integer(numcoords),
                as.double(data),
                as.double(coords[, 1]),
                as.double(coords[, 2]),
                as.double(coords[, 3]),
                as.double(numbins),
                as.double(bins),
                as.double(v0),
                as.double(v1),
                as.double(v2),
                as.double(maxdist),
                as.integer(type_dist),
                as.double(radius),
                PACKAGE = "GeoModels",
                DUP = TRUE,
                NAOK = TRUE)

        .extract_pairs(V)
    }

    ##################################################
    ##### Scatterplot based on distances
    ##################################################

    if (is.null(neighb)) {
        if (is.null(maxdist)) stop("maxdist must be provided when neighb is NULL")

        if (!bivariate && !spacetime) {
            if (length(maxdist) != 1) stop("maxdist must be a scalar for the univariate spatial case")

            bins <- seq(0, maxdist, length.out = numbins + 1)
            vvv <- .call_pairs(data, coords, bins, maxdist)

            lay <- .panel_layout(numbins)
            par(mfrow = lay)

            for (i in seq_len(numbins)) {
                sel <- vvv$v0 == bins[i]
                v111 <- vvv$v1[sel]
                v222 <- vvv$v2[sel]

                main <- paste("(", signif(bins[i], 3), ",", signif(bins[i + 1], 3), "]")

                if (length(v111) == 0 || length(v222) == 0) {
                    plot.new()
                    title(main = main)
                    mtext("No pairs", side = 3, line = -2)
                    message("There are no points between ", signif(bins[i], 3),
                            " and ", signif(bins[i + 1], 3))
                } else {
                    plot(v111, v222, col = "#481567FF", xlab = "", ylab = "",
                         main = main, ...)
                    abline(0, 1)
                }
            }
        }

        if (bivariate) {
            if (length(maxdist) == 1) maxdist <- rep(maxdist, 2)
            if (length(maxdist) != 2) stop("maxdist must have length one or two for bivariate=TRUE")

            data1 <- data[seq_len(ns[1])]
            data2 <- data[(ns[1] + 1):(ns[1] + ns[2])]
            coords1 <- coords[seq_len(ns[1]), , drop = FALSE]
            coords2 <- coords[(ns[1] + 1):(ns[1] + ns[2]), , drop = FALSE]

            bins1 <- seq(0, maxdist[1], length.out = numbins + 1)
            bins2 <- seq(0, maxdist[2], length.out = numbins + 1)

            vvv1 <- .call_pairs(data1, coords1, bins1, maxdist[1])
            vvv2 <- .call_pairs(data2, coords2, bins2, maxdist[2])

            lay <- .panel_layout(2 * numbins)
            par(mfrow = lay)

            for (i in seq_len(numbins)) {
                sel <- vvv1$v0 == bins1[i]
                v111 <- vvv1$v1[sel]
                v222 <- vvv1$v2[sel]
                cc <- if (length(v111) > 1 && length(v222) > 1) stats::cor(v111, v222) else NA_real_
                main <- paste("cor =", signif(cc, 3),
                              paste("(", signif(bins1[i], 3), ",", signif(bins1[i + 1], 3), "]"))

                if (length(v111) == 0 || length(v222) == 0) {
                    plot.new(); title(main = main); mtext("No pairs", side = 3, line = -2)
                } else {
                    plot(v111, v222, col = "#481567FF", xlab = "", ylab = "",
                         main = main, ...)
                    abline(0, 1)
                }
            }

            for (i in seq_len(numbins)) {
                sel <- vvv2$v0 == bins2[i]
                v111 <- vvv2$v1[sel]
                v222 <- vvv2$v2[sel]
                cc <- if (length(v111) > 1 && length(v222) > 1) stats::cor(v111, v222) else NA_real_
                main <- paste("cor =", signif(cc, 3),
                              paste("(", signif(bins2[i], 3), ",", signif(bins2[i + 1], 3), "]"))

                if (length(v111) == 0 || length(v222) == 0) {
                    plot.new(); title(main = main); mtext("No pairs", side = 3, line = -2)
                } else {
                    plot(v111, v222, col = "#481567FF", xlab = "", ylab = "",
                         main = main, ...)
                    abline(0, 1)
                }
            }
        }
    }

    ##################################################
    ##### Scatterplot based on neighborhood
    ##################################################

    if (!is.null(neighb)) {
        ln <- length(neighb)
        lay <- .panel_layout(ln)

        if (!bivariate && !spacetime) {
            par(mfrow = lay)
            for (i in seq_len(ln)) {
                sel <- GeoNeighIndex(coordx = coords, neighb = neighb[i])
                data1 <- data[sel$colidx]
                data2 <- data[sel$rowidx]
                plot(data1, data2, col = "#481567FF", xlab = "", ylab = "",
                     main = paste("Neighb =", neighb[i]), ...)
                abline(0, 1)
            }
        }

        if (bivariate) {
            stop("Neighborhood scatterplot for bivariate=TRUE is not implemented yet")
        }

        if (spacetime) {
            stop("Neighborhood scatterplot for spacetime=TRUE is not implemented yet")
        }
    }

    invisible(call)
}