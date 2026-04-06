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

    if (!is.matrix(x) || !is.matrix(y))  stop("x and y must be matrices")
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

spacetime_index <- function(coords, coordx_dyn = NULL, N, K = 4, coordt = NULL,
                            numtime, maxtime = 1, maxdist = NULL,
                            distance = "Eucl", radius = 1) {

  if (is.null(coordt)) stop("coordt is required for spacetime_index().")
  if (numtime <= 1L) stop("numtime must be >= 2.")
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) || K < 1)
    stop("K must be a positive integer.")
  K <- as.integer(K)

  maxt_int <- as.integer(round(maxtime))
  if (!is.finite(maxt_int) || maxt_int < 1L) maxt_int <- 1L
  maxt_int <- min(maxt_int, numtime - 1L)
  if (maxt_int <= 0L) maxt_int <- 0L

  # per-time coords + offsets globali
  if (is.null(coordx_dyn)) {
    ns <- rep.int(as.integer(N), numtime)
    offset <- as.integer(N) * (0:(numtime - 1L))
    getX <- function(t) coords
  } else {
    if (!is.list(coordx_dyn) || length(coordx_dyn) != numtime)
      stop("coordx_dyn must be a list of length = numtime.")
    ns <- vapply(coordx_dyn, nrow, integer(1))
    if (any(ns <= 0L)) stop("Each coordx_dyn[[t]] must have at least 1 row.")
    offset <- c(0L, cumsum(ns))[seq_len(numtime)]
    offset <- as.integer(offset)
    getX <- function(t) coordx_dyn[[t]]
  }

  ct <- as.numeric(coordt)
  if (length(ct) != numtime) stop("length(coordt) must equal numtime.")
  dt <- function(t, j) abs(ct[t + j] - ct[t])

  glob <- function(idx_local, t) as.integer(idx_local + offset[t])

  # helper: nn2Geo gives (rowidx=query(target), colidx=neighbor in ref)
  # we return vectors target, neighbor, dist
  knn_dir <- function(ref, query, K, drop_self = FALSE, x_equals_y = FALSE) {
    inf <- nn2Geo(ref, query, if (drop_self) K + 1L else K,
                  distance, maxdist, radius,
                  drop_first = drop_self, x_equals_y = x_equals_y)

    if (length(inf$rowidx) == 0L)
      return(list(tar = integer(0), nei = integer(0), d = numeric(0)))

    tar <- inf$rowidx
    nei <- inf$colidx
    d   <- inf$lags

    ok <- is.finite(nei) & is.finite(tar) & nei > 0 & tar > 0 & is.finite(d)
    list(tar = tar[ok], nei = nei[ok], d = d[ok])
  }

  out <- list(); out_i <- 1L

  ##############################################################
  # 1) Spatial within-time
  #    OUTPUT columns are: [target, neighbor, lagt, lags]
  #    so that later gb$colidx=sol[,2], gb$rowidx=sol[,1] => (neighbor, target) in C
  ##############################################################
  for (t in seq_len(numtime)) {
    Xt <- getX(t)
    inf <- knn_dir(ref = Xt, query = Xt, K = K, drop_self = TRUE, x_equals_y = TRUE)

    if (length(inf$tar) > 0L) {
      out[[out_i]] <- cbind(
        glob(inf$tar, t),                   # target
        glob(inf$nei, t),                   # neighbor
        rep.int(0, length(inf$tar)),        # lagt
        as.numeric(inf$d)                   # lags
      )
      out_i <- out_i + 1L
    }
  }

  if (maxt_int == 0L) {
    sol <- do.call(rbind, out)
    ok <- is.finite(sol[,1]) & is.finite(sol[,2]) & sol[,1] > 0 & sol[,2] > 0
    return(sol[ok, , drop = FALSE])
  }

  ##############################################################
  # 2) Pure temporal same-site: neighbor @ t  -> target @ t+j
  #    lags=0, lagt>0
  ##############################################################
  for (j in seq_len(maxt_int)) {
    for (t in seq_len(numtime - j)) {
      bb <- min(ns[t], ns[t + j])
      if (bb <= 0L) next
      ii <- seq_len(bb)

      out[[out_i]] <- cbind(
        glob(ii, t + j),                    # target
        glob(ii, t),                        # neighbor
        rep.int(dt(t, j), bb),              # lagt
        rep.int(0, bb)                      # lags
      )
      out_i <- out_i + 1L
    }
  }

  ##############################################################
  # 3) Spatio-temporal cross-time: neighbor in t -> target in t+j
  #    keep ONLY lags>0 to avoid duplicating Step 2
  ##############################################################
  for (j in seq_len(maxt_int)) {
    for (t in seq_len(numtime - j)) {

      Xn <- getX(t)        # neighbor candidates at time t
      Xt <- getX(t + j)    # targets at time t+j
      lagt_val <- dt(t, j)

      inf <- knn_dir(ref = Xn, query = Xt, K = K, drop_self = FALSE, x_equals_y = FALSE)

      if (length(inf$tar) > 0L) {
        sel <- is.finite(inf$d) & (inf$d > 0)
        if (any(sel)) {
          out[[out_i]] <- cbind(
            glob(inf$tar[sel], t + j),      # target @ t+j
            glob(inf$nei[sel], t),          # neighbor @ t
            rep.int(lagt_val, sum(sel)),
            as.numeric(inf$d[sel])
          )
          out_i <- out_i + 1L
        }
      }
    }
  }

  sol <- do.call(rbind, out)
  ok <- is.finite(sol[, 1]) & is.finite(sol[, 2]) & sol[, 1] > 0 & sol[, 2] > 0
  sol[ok, , drop = FALSE]
}


  ##############################################################
 bivariate_index <- function(coords, coordx_dyn = NULL, N, K = 4,
                            maxdist, distance, radius) {

  # K può essere scalare o tripletta: (within1, cross, within2)
  if (length(K) == 3) { K1 <- K[1]; K2 <- K[2]; K3 <- K[3] } else { K1 <- K2 <- K3 <- K }

  # maxdist può essere scalare o tripletta
  if (length(maxdist) == 3) {
    maxdist1 <- maxdist[1]; maxdist2 <- maxdist[2]; maxdist3 <- maxdist[3]
  } else {
    maxdist1 <- maxdist2 <- maxdist3 <- maxdist
  }

  auto_tol <- function(A, B) {
    scaleA <- suppressWarnings(max(abs(range(A, finite = TRUE))))
    scaleB <- suppressWarnings(max(abs(range(B, finite = TRUE))))
    scale  <- max(scaleA, scaleB)
    if (!is.finite(scale) || scale <= 0) scale <- 1
    100 * .Machine$double.eps * scale
  }

  ##############################################################
  # CASE 1: coordx_dyn NULL (coords contiene var1 poi var2, co-located)
  ##############################################################
  if (is.null(coordx_dyn)) {

    n_half <- as.integer(N / 2)
    cm <- coords[1:n_half, , drop = FALSE]

    # within var1
    inf1 <- nn2Geo(cm, cm, K1 + 1, distance, maxdist1, radius,
                   drop_first = TRUE, x_equals_y = TRUE)
    aa1 <- cbind(inf1$rowidx, inf1$colidx, 0L, 0L, inf1$lags)

    # within var2
    inf2 <- nn2Geo(cm, cm, K3 + 1, distance, maxdist3, radius,
                   drop_first = TRUE, x_equals_y = TRUE)
    aa2 <- cbind(inf2$rowidx + n_half, inf2$colidx + n_half, 1L, 1L, inf2$lags)

    # cross-neighbors: ONLY 2 -> 1  (neighbors in var2 for each target in var1)
    inf21 <- nn2Geo(cm, cm, K2 + 1, distance, maxdist2, radius,
                    drop_first = TRUE, x_equals_y = TRUE)
    # qui inf21$lags > 0 per costruzione (drop_first=TRUE, x_equals_y=TRUE), ma filtriamo per sicurezza:
    sel <- is.finite(inf21$lags) & (inf21$lags > 0) & (inf21$colidx > 0) & (inf21$rowidx > 0)
    aa_cross <- cbind(
      inf21$rowidx[sel],                    # target var1
      inf21$colidx[sel] + n_half,           # neighbor var2 (shift)
      0L, 1L, inf21$lags[sel]
    )

    # collocated: ONLY 2 -> 1  (true collocated by construction)
    ii <- seq_len(n_half)
    aa_coll <- cbind(ii, ii + n_half, 0L, 1L, 0)

    return(rbind(aa1, aa2, aa_cross, aa_coll))
  }

  ##############################################################
  # CASE 2: coordx_dyn list length 2 (var1 grid, var2 grid) possibly different
  ##############################################################
  if (!is.list(coordx_dyn) || length(coordx_dyn) != 2)
    stop("For bivariate case with coordx_dyn, it must be a list of length 2: list(var1_coords, var2_coords).")

  X1 <- coordx_dyn[[1]]  # var1 targets
  X2 <- coordx_dyn[[2]]  # var2 neighbors
  n1 <- nrow(X1); n2 <- nrow(X2)

  # within var1
  inf1 <- nn2Geo(X1, X1, K1 + 1, distance, maxdist1, radius,
                 drop_first = TRUE, x_equals_y = TRUE)
  aa1 <- cbind(inf1$rowidx, inf1$colidx, 0L, 0L, inf1$lags)

  # within var2
  inf2 <- nn2Geo(X2, X2, K3 + 1, distance, maxdist3, radius,
                 drop_first = TRUE, x_equals_y = TRUE)
  aa2 <- cbind(inf2$rowidx + n1, inf2$colidx + n1, 1L, 1L, inf2$lags)

  # cross-neighbors: ONLY 2 -> 1
  # ref = X2, query = X1 => for each target in var1, pick K2 neighbors from var2
  inf_cross <- nn2Geo(X2, X1, K2, distance, maxdist2, radius,
                      drop_first = FALSE, x_equals_y = FALSE)

  sel <- is.finite(inf_cross$lags) & (inf_cross$lags > 0) &
         (inf_cross$colidx > 0) & (inf_cross$rowidx > 0)

  aa_cross <- cbind(
    inf_cross$rowidx[sel],          # target var1 (1..n1)
    inf_cross$colidx[sel] + n1,     # neighbor var2 (shift)
    0L, 1L, inf_cross$lags[sel]
  )

  # collocated: ONLY 2 -> 1, robust
  tol <- auto_tol(X1, X2)
  same_grid <- (n1 == n2) && isTRUE(all(abs(X1 - X2) <= tol))

  if (same_grid) {
    ii <- seq_len(n1)
    aa_coll <- cbind(ii, ii + n1, 0L, 1L, 0)
  } else {
    # nearest var2 for each var1 target, with real distance; set to 0 if within tol
    nn12 <- nn2Geo(X2, X1, K = 1, distance, maxdist2, radius,
                   drop_first = FALSE, x_equals_y = FALSE)

    l12 <- nn12$lags
    l12[is.finite(l12) & (l12 <= tol)] <- 0

    selc <- is.finite(nn12$rowidx) & is.finite(nn12$colidx) &
            (nn12$rowidx > 0) & (nn12$colidx > 0) & is.finite(l12)

    aa_coll <- cbind(nn12$rowidx[selc], nn12$colidx[selc] + n1, 0L, 1L, l12[selc])
  }

  rbind(aa1, aa2, aa_cross, aa_coll)
}
 
 ##############################################################
## Hard-core greedy matching - versione SEMPLICE
## Non cerca di matchare il numero di iid
##############################################################
hardcore_match_keep_target <- function(rowidx, colidx, m_target,
                                            n_nodes = NULL,
                                            batch_mult = 3L,
                                            max_rounds = 10L) {
  stopifnot(length(rowidx) == length(colidx))
  d <- length(rowidx)
  if (d == 0L || m_target <= 0L) return(integer(0))

  if (is.null(n_nodes)) n_nodes <- max(rowidx, colidx, na.rm = TRUE)

  m_cap <- as.integer(floor(n_nodes / 2))
  m_target <- as.integer(min(m_target, m_cap))
  if (m_target <= 0L) return(integer(0))

  used <- rep.int(FALSE, n_nodes)
  keep <- integer(m_target)
  cnt  <- 0L

  # per estrarre senza rimpiazzo in batch
  remaining <- d
  # lavoriamo su un vettore di indici da cui pescare
  pool <- seq_len(d)

  rounds <- 0L
  while (cnt < m_target && remaining > 0L && rounds < max_rounds) {
    rounds <- rounds + 1L

    # dimensione batch: proporzionale a quanto manca
    need <- m_target - cnt
    bsz  <- min(remaining, as.integer(batch_mult * need))
    if (bsz <= 0L) break

    # prendi un batch casuale dal pool (senza rimpiazzo)
    take_pos <- sample.int(remaining, size = bsz, replace = FALSE)
    batch <- pool[take_pos]

    # rimuovi batch dal pool (swap-delete veloce)
    pool[take_pos] <- pool[remaining - seq_len(bsz) + 1L]
    remaining <- remaining - bsz
    pool <- pool[seq_len(remaining)]

    # greedy sul batch
    for (e in batch) {
      i <- rowidx[e]; j <- colidx[e]
      if (i > 0L && j > 0L && !used[i] && !used[j]) {
        cnt <- cnt + 1L
        keep[cnt] <- e
        used[i] <- TRUE
        used[j] <- TRUE
        if (cnt >= m_target) break
      }
    }

    # se stai riempiendo lentamente, aumenta batch_mult dinamicamente
    # (opzionale, ma spesso aiuta)
    if (cnt < m_target && rounds >= 2L) batch_mult <- min(20L, batch_mult + 1L)
  }

  if (cnt > 0L) keep[seq_len(cnt)] else integer(0)
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
if (thin_method == "match") {

  n_all <- length(gb$rowidx)
  if (n_all == 0L) return(gb)

  n_nodes  <- max(gb$rowidx, gb$colidx, na.rm = TRUE)
  m_target <- as.integer(round(p_neighb * n_all))

  keep_idx <- hardcore_match_keep_target(
    rowidx = gb$rowidx,
    colidx = gb$colidx,
    m_target = m_target,
    n_nodes = n_nodes,
    batch_mult = 3L,
    max_rounds = 10L
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

      isS<- (tt == 0) & (h  > 0)
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
