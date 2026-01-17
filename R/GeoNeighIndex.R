GeoNeighIndex <- function(coordx, coordy = NULL, coordz = NULL, coordt = NULL, coordx_dyn = NULL,
                         distance = "Eucl", neighb = 4, maxdist = NULL, maxtime = 1,
                         radius = 1, bivariate = FALSE, p_neighb = 1, thin_method = "iid")
{
  if (!is.numeric(neighb) || length(neighb) != 1L || !is.finite(neighb) || neighb <= 0 || neighb %% 1 != 0)
    stop("neighb must be a single positive integer.")
  K <- as.integer(neighb)

  thin_method <- match.arg(thin_method, c("iid", "match"))

  indices <- function(nn_idx, nn_dists, drop_first = FALSE) {
    if (drop_first) {
      nn_idx   <- nn_idx[, -1, drop = FALSE]
      nn_dists <- nn_dists[, -1, drop = FALSE]
    }
    k <- ncol(nn_idx)
    rowidx <- rep(seq_len(nrow(nn_idx)), times = k)  # aligns with as.vector() column-major
    colidx <- as.vector(nn_idx)
    lags   <- as.vector(nn_dists)
    list(rowidx = rowidx, colidx = colidx, lags = lags)
  }

  ##############################################################
  ##############################################################
  nn2Geo <- function(x, y, K = 1, distance = 0, maxdist = NULL, radius = 1,
                     drop_first = FALSE, x_equals_y = FALSE) {

    if (!is.matrix(x) || !is.matrix(y)) stop("x e y devono essere matrici")
    n_x <- nrow(x); n_y <- nrow(y)
    K <- min(K, n_x, n_y)

    project_coords <- function(coords) {
      prj <- mapproj::mapproject(coords[, 1], coords[, 2], projection = "sinusoidal")
      radius * cbind(prj$x, prj$y)
    }

    if (is.null(maxdist)) {
      nearest <- nabor::knn(x, y, k = K)
    } else {
      if (distance %in% c(1, 2)) {
        x_proj <- project_coords(x)
        y_proj <- project_coords(y)
        nearest <- nabor::knn(x_proj, y_proj, k = K, radius = maxdist)
      } else {
        nearest <- nabor::knn(x, y, k = K, radius = maxdist)
      }
    }

    ## Geodesic distances if requested (kept as in your logic)
    if (distance %in% c(1, 2)) {

      dists_mat <- matrix(NA_real_, nrow = n_y, ncol = K)

      for (i in seq_len(n_y)) {
        idx <- nearest$nn.idx[i, ]
        valid <- idx > 0
        if (!any(valid)) next

        di <- fields::rdist.earth.vec(
          x1 = matrix(y[i, ], nrow = 1),
          x2 = x[idx[valid], , drop = FALSE],
          miles = FALSE, R = 1
        )

        dists_mat[i, seq_len(length(di))] <- di
      }

      if (isTRUE(x_equals_y)) dists_mat[, 1] <- 0

      if (distance == 2) {
        nearest$nn.dists <- radius * dists_mat
      } else {
        nearest$nn.dists <- 2 * radius * sin(0.5 * dists_mat)
      }
    }

    sol <- indices(nearest$nn.idx, nearest$nn.dists, drop_first = drop_first)

    if (!is.null(maxdist)) {
      sel <- sol$colidx > 0
      return(list(lags = sol$lags[sel], rowidx = sol$rowidx[sel], colidx = sol$colidx[sel]))
    } else {
      return(list(lags = sol$lags, rowidx = sol$rowidx, colidx = sol$colidx))
    }
  }

  ##############################################################
  spacetime_index <- function(coords, coordx_dyn = NULL, N, K = 4, coordt = NULL,
                              numtime, maxtime = 1, maxdist = NULL,
                              distance = "Eucl", radius = 1) {

    m_s <- vector("list", numtime)

    # 1) Marginal spatial indexes (per time)
    if (is.null(coordx_dyn)) {
      inf <- nn2Geo(coords, coords, K + 1, distance, maxdist, radius,
                    drop_first = TRUE, x_equals_y = TRUE)
      aa  <- cbind(inf$rowidx, inf$colidx)
      lag <- inf$lags
      offset_seq <- N * (0:(numtime - 1))
      for (i in seq_len(numtime)) {
        m_s[[i]] <- cbind(aa + offset_seq[i], 0L, lag)
      }
    } else {
      ns <- vapply(coordx_dyn, nrow, FUN.VALUE = integer(1))
      for (i in seq_len(numtime)) {
        inf <- nn2Geo(coordx_dyn[[i]], coordx_dyn[[i]], K + 1, distance, maxdist, radius,
                      drop_first = TRUE, x_equals_y = TRUE)
        aa  <- cbind(inf$rowidx, inf$colidx)
        lag <- inf$lags
        m_s[[i]] <- cbind(aa + ns[i] * (i - 1L), 0L, lag)
      }
    }

    # 2) Temporal & spatio-temporal
    ct <- as.matrix(coordt)
    time_dists <- sort(unique(nabor::knn(ct, ct, k = round(maxtime) + 1)$nn.dists))
    nn <- time_dists[time_dists > 0]
    tnn <- length(nn)

    m_t  <- vector("list", tnn * (numtime - tnn))
    m_st <- vector("list", tnn * (numtime - tnn))

    counter <- 1L
    for (j in seq_along(nn)) {
      for (k in seq_len(numtime - j)) {
        n1 <- nrow(m_s[[k]])
        n2 <- nrow(m_s[[k + j]])
        bb <- min(n1, n2)

        m_t[[counter]] <- cbind(
          m_s[[k]][1:bb, 1],
          m_s[[k + j]][1:bb, 1],
          rep_len(nn[j], bb),
          rep_len(0L, bb)
        )

        m_st[[counter]] <- cbind(
          m_s[[k]][1:bb, 1],
          m_s[[k + j]][1:bb, 2],
          rep_len(nn[j], bb),
          m_s[[k]][1:bb, 4]
        )

        counter <- counter + 1L
      }
    }

    do.call(rbind, c(m_s, m_t, m_st))
  }

  ##############################################################
  bivariate_index <- function(coords, coordx_dyn = NULL, N, K = 4, maxdist, distance, radius) {
    if (length(K) == 3) { K1 <- K[1]; K2 <- K[2]; K3 <- K[3] } else { K1 <- K2 <- K3 <- K }
    if (length(maxdist) == 3) { maxdist1 <- maxdist[1]; maxdist2 <- maxdist[2]; maxdist3 <- maxdist[3] }
    else { maxdist1 <- maxdist2 <- maxdist3 <- maxdist }

    if (is.null(coordx_dyn)) {
      n_half <- as.integer(N / 2)
      cm <- coords[1:n_half, ]

      inf1 <- nn2Geo(cm, cm, K1 + 1, distance, maxdist1, radius, drop_first = TRUE, x_equals_y = TRUE)
      inf2 <- nn2Geo(cm, cm, K3 + 1, distance, maxdist3, radius, drop_first = TRUE, x_equals_y = TRUE)
      inf3 <- nn2Geo(cm, cm, K2 + 1, distance, maxdist2, radius, drop_first = TRUE, x_equals_y = TRUE)

      aa1 <- cbind(inf1$rowidx, inf1$colidx, 0L, 0L, inf1$lags)
      aa2 <- cbind(inf2$rowidx + n_half, inf2$colidx + n_half, 1L, 1L, inf2$lags)
      aa3 <- cbind(inf3$rowidx, inf3$colidx + n_half, 0L, 1L, inf3$lags)
      aa4 <- cbind(inf3$rowidx + n_half, inf3$colidx, 1L, 0L, inf3$lags)

      a5  <- nabor::knn(cm, k = K2)
      aa5 <- cbind(a5$nn.idx[, 1], a5$nn.idx[, 1] + n_half, 0L, 1L, 0)
      aa6 <- cbind(a5$nn.idx[, 1] + n_half, a5$nn.idx[, 1], 1L, 0L, 0)

    } else {
      ns <- vapply(coordx_dyn, nrow, integer(1))

      inf1 <- nn2Geo(coordx_dyn[[1]], coordx_dyn[[1]], K1 + 1, distance, maxdist1, radius,
                     drop_first = TRUE, x_equals_y = TRUE)
      inf2 <- nn2Geo(coordx_dyn[[2]], coordx_dyn[[2]], K3 + 1, distance, maxdist3, radius,
                     drop_first = TRUE, x_equals_y = TRUE)
      inf3 <- nn2Geo(coordx_dyn[[1]], coordx_dyn[[2]], K2 + 1, distance, maxdist2, radius,
                     drop_first = FALSE, x_equals_y = FALSE)

      aa1 <- cbind(inf1$rowidx, inf1$colidx, 0L, 0L, inf1$lags)
      aa2 <- cbind(inf2$rowidx + ns[1], inf2$colidx + ns[1], 1L, 1L, inf2$lags)
      aa3 <- cbind(inf3$rowidx, inf3$colidx + ns[1], 0L, 1L, inf3$lags)
      aa4 <- cbind(inf3$rowidx + ns[1], inf3$colidx, 1L, 0L, inf3$lags)

      a5  <- nabor::knn(coordx_dyn[[1]], coordx_dyn[[2]], k = K2 + 1)
      aa5 <- cbind(a5$nn.idx[, 1], a5$nn.idx[, 1] + ns[1], 0L, 1L, 0)
      aa6 <- cbind(a5$nn.idx[, 1] + ns[1], a5$nn.idx[, 1], 1L, 0L, 0)
    }

    rbind(aa1, aa2, aa3, aa4, aa5, aa6)
  }

  ##############################################################
  ## Hard-core greedy matching on candidate pairs
  ##############################################################
  hardcore_match_keep <- function(rowidx, colidx, n_nodes = NULL, priority = NULL) {
    stopifnot(length(rowidx) == length(colidx))
    d <- length(rowidx)
    if (d == 0L) return(integer(0))

    if (is.null(n_nodes)) n_nodes <- max(rowidx, colidx, na.rm = TRUE)

    if (is.null(priority)) priority <- stats::runif(d) # i.i.d. continuous marks
    ord <- order(priority)

    used <- rep(FALSE, n_nodes)
    keep <- logical(d)

    for (e in ord) {
      i <- rowidx[e]; j <- colidx[e]
      if (i <= 0L || j <= 0L) next
      if (!used[i] && !used[j]) {
        keep[e] <- TRUE
        used[i] <- TRUE
        used[j] <- TRUE
      }
    }

    which(keep)
  }

  ######################################
  ######### start ######################
  ######################################

  spatial   <- TRUE
  spacetime <- FALSE

  if (!is.null(coordx_dyn))
    if (!is.list(coordx_dyn)) stop(" coordx_dyn must be a list")

  if (!is.null(coordt))
    if (is.numeric(coordt) && is.numeric(maxtime))
      if (length(coordt) > 1 && length(maxtime) >= 1) spacetime <- TRUE

  distance <- CheckDistance(distance)
  spatial  <- !spacetime && !bivariate

  ## --- keep your coordinate parsing as-is to preserve behavior ---
  if (!is.null(coordy)) {
    if (is.null(coordz)) {
      coordy <- coordx[, 2]
      coordx <- coordx[, 1]
      coords <- cbind(coordx, coordy)
    } else {
      coordz <- coordx[, 3]
      coordy <- coordx[, 2]
      coordx <- coordx[, 1]
      coords <- cbind(coordx, coordy, coordz)
    }
    numcoord <- nrow(coords)
  } else {
    if (!bivariate) { coords <- coordx; numcoord <- nrow(coords) }
    if (bivariate) {
      if (is.null(coordx_dyn)) { coords <- coordx; numcoord <- nrow(coords) }
      else { coords <- 1; numcoord <- 1 }
    }
  }

  ##########################
  ## Build pair indices
  ##########################
  if (spatial) {
    sol <- nn2Geo(coords, coords, K + 1, distance, maxdist, radius,
                  drop_first = TRUE, x_equals_y = TRUE)
    gb <- list(
      colidx = sol$colidx,
      rowidx = sol$rowidx,
      lags   = sol$lags,
      maxdist = maxdist,
      neighb  = neighb
    )
  }

  if (spacetime) {
    numtime <- length(coordt)
    sol <- spacetime_index(coords, coordx_dyn, numcoord, K, coordt, numtime, maxtime, maxdist, distance, radius)
    gb <- list(
      colidx = sol[, 2],
      rowidx = sol[, 1],
      lags   = sol[, 4],
      lagt   = sol[, 3]
    )
  }

  if (bivariate) {
    sol <- bivariate_index(coords, coordx_dyn, numcoord, K, maxdist, distance, radius)
    gb <- list(
      colidx = sol[, 2],
      rowidx = sol[, 1],
      lags   = sol[, 5],
      first  = sol[, 3],
      second = sol[, 4],
      maxdist = maxdist,
      neighb  = neighb
    )
  }

  ## --- p_neighb checks (unchanged) ---
  if (!is.numeric(p_neighb) || length(p_neighb) != 1L || !is.finite(p_neighb))
    stop("p_neighb must be a finite numeric scalar.")
  if (p_neighb <= 0 || p_neighb > 1)
    stop("p_neighb must be in (0,1].")
  thinning_on <- (p_neighb < 1)

  ## Defensive: remove invalid indices (can happen with radius-based knn)
  if (length(gb$rowidx) > 0L) {
    valid <- (gb$rowidx > 0L) & (gb$colidx > 0L)
    if (!all(valid)) {
      gb$rowidx <- gb$rowidx[valid]
      gb$colidx <- gb$colidx[valid]
      gb$lags   <- gb$lags[valid]
      if (!is.null(gb$lagt))   gb$lagt   <- gb$lagt[valid]
      if (!is.null(gb$first))  gb$first  <- gb$first[valid]
      if (!is.null(gb$second)) gb$second <- gb$second[valid]
    }
  }

  #############################################################
  ################## thinning / matching #######################
  #############################################################

  ## --- hard-core matching: independent of p_neighb ---
  if (thin_method == "match") {

    if (p_neighb < 1) {
      warning("thin_method='match': p_neighb is ignored (matching selection is applied).", call. = FALSE)
    }

    n_all <- length(gb$rowidx)
    if (n_all == 0L) return(gb)

    keep_idx <- hardcore_match_keep(
      rowidx  = gb$rowidx,
      colidx  = gb$colidx,
      n_nodes = max(gb$rowidx, gb$colidx, na.rm = TRUE),
      priority = NULL
    )

    gb$colidx <- gb$colidx[keep_idx]
    gb$rowidx <- gb$rowidx[keep_idx]
    gb$lags   <- gb$lags[keep_idx]
    if (!is.null(gb$lagt))   gb$lagt   <- gb$lagt[keep_idx]
    if (!is.null(gb$first))  gb$first  <- gb$first[keep_idx]
    if (!is.null(gb$second)) gb$second <- gb$second[keep_idx]

    return(gb)
  }

  ## --- iid thinning: only if thin_method == "iid" AND p_neighb < 1 ---
  if (thin_method == "iid" && thinning_on) {

    if (!(spatial || spacetime || bivariate))
      stop("p_neighb thinning implemented only for spatial, spacetime, or bivariate.")

    ## fixed internal tuning (YOU choose)
    Kbins_s  <- 10L
    Kbins_t  <- 10L
    lambda_s <- 3
    lambda_t <- 2

    n_all <- length(gb$rowidx)
    if (n_all == 0L) return(gb)

    make_u <- function(x, Kbins) {
      Kbins <- as.integer(Kbins)
      Kbins <- max(1L, min(Kbins, length(x)))

      qbreaks <- stats::quantile(x, probs = seq(0, 1, length.out = Kbins + 1),
                                na.rm = TRUE, type = 7)
      if (any(!is.finite(qbreaks))) stop("Non-finite quantiles in thinning.")
      if (any(diff(qbreaks) == 0)) {
        x2 <- x + stats::runif(length(x), -1e-12, 1e-12)
        qbreaks <- stats::quantile(x2, probs = seq(0, 1, length.out = Kbins + 1),
                                  na.rm = TRUE, type = 7)
      }

      bin <- cut(x, breaks = qbreaks, include.lowest = TRUE, labels = FALSE)

      centers <- vapply(seq_len(Kbins), function(k) stats::median(x[bin == k]), numeric(1))
      ord  <- order(centers)
      rank <- match(seq_len(Kbins), ord)
      ubin <- (rank - 1) / max(1, (Kbins - 1))

      ubin[bin]
    }

    if (spatial) {

      h <- gb$lags
      Kall <- length(h)
      if (Kall == 0L) return(gb)

      Kbins <- as.integer(Kbins_s)
      Kbins <- max(1L, min(Kbins, Kall))

      qbreaks <- stats::quantile(
        h, probs = seq(0, 1, length.out = Kbins + 1),
        na.rm = TRUE, type = 7
      )
      if (any(!is.finite(qbreaks))) stop("Non-finite lag quantiles in thinning.")
      if (any(diff(qbreaks) == 0)) {
        h2 <- h + stats::runif(Kall, -1e-12, 1e-12)
        qbreaks <- stats::quantile(
          h2, probs = seq(0, 1, length.out = Kbins + 1),
          na.rm = TRUE, type = 7
        )
      }
      bin <- cut(h, breaks = qbreaks, include.lowest = TRUE, labels = FALSE)
      if (anyNA(bin)) bin[is.na(bin)] <- 1L
      centers <- vapply(seq_len(Kbins), function(k) stats::median(h[bin == k]), numeric(1))
      ordc <- order(centers)
      rankc <- match(seq_len(Kbins), ordc)
      u <- (rankc - 1) / max(1, (Kbins - 1))
      wbin <- exp(-lambda_s * u)
      nb <- tabulate(bin, nbins = Kbins)
      denom <- sum(wbin * nb)
      if (!is.finite(denom) || denom <= 0) denom <- sum(nb)
      pbin <- (p_neighb * Kall) * (wbin / denom)
      pbin <- pmin(1, pmax(0, pbin))

      p_ij <- pbin[bin]

      B <- stats::rbinom(Kall, 1, p_ij)
      keep_idx <- which(B == 1L)

      gb$colidx <- gb$colidx[keep_idx]
      gb$rowidx <- gb$rowidx[keep_idx]
      gb$lags   <- gb$lags[keep_idx]

    } else if (spacetime) {

      h  <- gb$lags
      tt <- gb$lagt

      isS  <- (tt == 0) & (h  > 0)
      isT  <- (h  == 0) & (tt > 0)
      isST <- (tt > 0) & (h  > 0)

      idx_h <- which(h  > 0)
      idx_t <- which(tt > 0)

      uh <- rep(NA_real_, n_all)
      ut <- rep(NA_real_, n_all)

      if (length(idx_h) > 0L) uh[idx_h] <- make_u(h[idx_h],  Kbins_s)
      if (length(idx_t) > 0L) ut[idx_t] <- make_u(tt[idx_t], Kbins_t)

      w <- rep(0, n_all)
      if (any(isS))  w[isS]  <- exp(-lambda_s * uh[isS])
      if (any(isT))  w[isT]  <- exp(-lambda_t * ut[isT])
      if (any(isST)) w[isST] <- exp(-lambda_s * uh[isST] - lambda_t * ut[isST])

      denom <- sum(w)
      if (!is.finite(denom) || denom <= 0) {
        w <- rep(1, n_all); denom <- n_all
      }

      p_ij <- (p_neighb * n_all) * (w / denom)
      p_ij <- pmin(1, pmax(0, p_ij))

      B <- stats::rbinom(n_all, 1, p_ij)
      keep_idx <- which(B == 1L)

      gb$colidx <- gb$colidx[keep_idx]
      gb$rowidx <- gb$rowidx[keep_idx]
      gb$lags   <- gb$lags[keep_idx]
      gb$lagt   <- gb$lagt[keep_idx]

    } else if (bivariate) {

      h <- gb$lags
      Kall <- length(h)
      if (Kall == 0L) return(gb)

      uh <- make_u(h, Kbins_s)
      w  <- exp(-lambda_s * uh)
      denom <- sum(w)
      if (!is.finite(denom) || denom <= 0) {
        w <- rep(1, Kall); denom <- Kall
      }

      p_ij <- (p_neighb * Kall) * (w / denom)
      p_ij <- pmin(1, pmax(0, p_ij))

      B <- stats::rbinom(Kall, 1, p_ij)
      keep_idx <- which(B == 1L)

      gb$colidx <- gb$colidx[keep_idx]
      gb$rowidx <- gb$rowidx[keep_idx]
      gb$lags   <- gb$lags[keep_idx]
      if (!is.null(gb$first))  gb$first  <- gb$first[keep_idx]
      if (!is.null(gb$second)) gb$second <- gb$second[keep_idx]
    }
  }

  return(gb)
}
