GeoNeighIndex <- function(coordx, coordy = NULL, coordz = NULL, coordt = NULL,
                          coordx_dyn = NULL,
                          distance = "Eucl", neighb = 4, maxdist = NULL,
                          maxtime = 1, radius = 1, bivariate = FALSE,
                          p_neighb = 1, thin_method = "bernoulli")
{
  ##############################################################
  ## Checks
  ##############################################################

  if (!is.numeric(neighb) ||
      !(length(neighb) %in% c(1L, 3L)) ||
      any(!is.finite(neighb)) ||
      any(neighb <= 0) ||
      any(neighb %% 1 != 0)) {
    stop("neighb must be a positive integer, or a length-3 positive integer vector in the bivariate case.")
  }

  if (!bivariate && length(neighb) != 1L) {
    stop("neighb can have length 3 only when bivariate = TRUE.")
  }

  K <- as.integer(neighb)

  ## Canonical thinning-method names.
  ## Accepted inputs include "bernoulli", "Bernoulli",
  ## "fixedbudget", and "FixedBudget". The old aliases
  ## "targetbalanced" and "TargetBalanced" are also accepted
  ## for backward compatibility and mapped to "FixedBudget".
  ## Internally we use only "Bernoulli" and "FixedBudget".
  if (!is.character(thin_method) || length(thin_method) != 1L || is.na(thin_method)) {
    stop("thin_method must be a character scalar.")
  }

  thin_method <- tolower(thin_method)

  thin_method <- match.arg(
    thin_method,
    c("bernoulli", "fixedbudget", "targetbalanced")
  )

  thin_method <- switch(
    thin_method,
    bernoulli = "Bernoulli",
    fixedbudget = "FixedBudget",
    targetbalanced = "FixedBudget"
  )
 
  indices <- function(nn_idx, nn_dists, drop_first = FALSE) {

    if (drop_first) {
      nn_idx   <- nn_idx[, -1, drop = FALSE]
      nn_dists <- nn_dists[, -1, drop = FALSE]
    }

    k <- ncol(nn_idx)

    rowidx <- rep(seq_len(nrow(nn_idx)), times = k)
    colidx <- as.vector(nn_idx)
    lags   <- as.vector(nn_dists)

    list(
      rowidx = rowidx,
      colidx = colidx,
      lags = lags
    )
  }

  ##############################################################
  ## Spatial nearest-neighbor helper
  ##############################################################

  nn2Geo <- function(x, y, K = 1, distance = 0, maxdist = NULL, radius = 1,
                     drop_first = FALSE, x_equals_y = FALSE) {

    if (!is.matrix(x) || !is.matrix(y)) {
      stop("x and y must be matrices")
    }

    n_x <- nrow(x)
    n_y <- nrow(y)

    K <- as.integer(K)
    K <- min(K, n_x, n_y)

    project_coords <- function(coords) {

      prj <- mapproj::mapproject(
        coords[, 1],
        coords[, 2],
        projection = "sinusoidal"
      )

      radius * cbind(prj$x, prj$y)
    }

    if (is.null(maxdist)) {

      nearest <- nabor::knn(x, y, k = K)

    } else {

      if (distance %in% c(1, 2)) {

        x_proj <- project_coords(x)
        y_proj <- project_coords(y)

        nearest <- nabor::knn(
          x_proj,
          y_proj,
          k = K,
          radius = maxdist
        )

      } else {

        nearest <- nabor::knn(
          x,
          y,
          k = K,
          radius = maxdist
        )
      }
    }

    if (distance %in% c(1, 2)) {

      dists_mat <- matrix(NA_real_, nrow = n_y, ncol = K)

      for (i in seq_len(n_y)) {

        idx <- nearest$nn.idx[i, ]
        valid <- idx > 0

        if (!any(valid)) {
          next
        }

        di <- fields::rdist.earth.vec(
          x1 = matrix(y[i, ], nrow = 1),
          x2 = x[idx[valid], , drop = FALSE],
          miles = FALSE,
          R = 1
        )

        dists_mat[i, seq_len(length(di))] <- di
      }

      if (isTRUE(x_equals_y)) {
        dists_mat[, 1] <- 0
      }

      if (distance == 2) {
        nearest$nn.dists <- radius * dists_mat
      } else {
        nearest$nn.dists <- 2 * radius * sin(0.5 * dists_mat)
      }
    }

    sol <- indices(
      nearest$nn.idx,
      nearest$nn.dists,
      drop_first = drop_first
    )

    if (!is.null(maxdist)) {

      sel <- sol$colidx > 0

      return(list(
        lags = sol$lags[sel],
        rowidx = sol$rowidx[sel],
        colidx = sol$colidx[sel]
      ))

    } else {

      return(list(
        lags = sol$lags,
        rowidx = sol$rowidx,
        colidx = sol$colidx
      ))
    }
  }

  ##############################################################
  ## Space-time index
  ##############################################################

  spacetime_index <- function(coords, coordx_dyn = NULL, N, K = 4,
                              coordt = NULL, numtime,
                              maxtime = 1, maxdist = NULL,
                              distance = "Eucl", radius = 1) {

    if (is.null(coordt)) {
      stop("coordt is required for spacetime_index().")
    }

    if (numtime <= 1L) {
      stop("numtime must be >= 2.")
    }

    if (!is.numeric(K) || length(K) != 1L || !is.finite(K) || K < 1) {
      stop("K must be a positive integer.")
    }

    K <- as.integer(K)

    maxt_int <- as.integer(round(maxtime))

    if (!is.finite(maxt_int) || maxt_int < 1L) {
      maxt_int <- 1L
    }

    maxt_int <- min(maxt_int, numtime - 1L)

    if (maxt_int <= 0L) {
      maxt_int <- 0L
    }

    if (is.null(coordx_dyn)) {

      ns <- rep.int(as.integer(N), numtime)
      offset <- as.integer(N) * (0:(numtime - 1L))

      getX <- function(t) coords

    } else {

      if (!is.list(coordx_dyn) || length(coordx_dyn) != numtime) {
        stop("coordx_dyn must be a list of length = numtime.")
      }

      ns <- vapply(coordx_dyn, nrow, integer(1))

      if (any(ns <= 0L)) {
        stop("Each coordx_dyn[[t]] must have at least 1 row.")
      }

      offset <- c(0L, cumsum(ns))[seq_len(numtime)]
      offset <- as.integer(offset)

      getX <- function(t) coordx_dyn[[t]]
    }

    ct <- as.numeric(coordt)

    if (length(ct) != numtime) {
      stop("length(coordt) must equal numtime.")
    }

    dt <- function(t, j) {
      abs(ct[t + j] - ct[t])
    }

    glob <- function(idx_local, t) {
      as.integer(idx_local + offset[t])
    }

    knn_dir <- function(ref, query, K,
                        drop_self = FALSE,
                        x_equals_y = FALSE) {

      inf <- nn2Geo(
        ref,
        query,
        if (drop_self) K + 1L else K,
        distance,
        maxdist,
        radius,
        drop_first = drop_self,
        x_equals_y = x_equals_y
      )

      if (length(inf$rowidx) == 0L) {
        return(list(
          tar = integer(0),
          nei = integer(0),
          d = numeric(0)
        ))
      }

      tar <- inf$rowidx
      nei <- inf$colidx
      d   <- inf$lags

      ok <- is.finite(nei) &
        is.finite(tar) &
        nei > 0 &
        tar > 0 &
        is.finite(d)

      list(
        tar = tar[ok],
        nei = nei[ok],
        d = d[ok]
      )
    }

    out <- list()
    out_i <- 1L

    ##############################################################
    ## 1) Spatial within-time
    ##############################################################

    for (t in seq_len(numtime)) {

      Xt <- getX(t)

      inf <- knn_dir(
        ref = Xt,
        query = Xt,
        K = K,
        drop_self = TRUE,
        x_equals_y = TRUE
      )

      if (length(inf$tar) > 0L) {

        out[[out_i]] <- cbind(
          glob(inf$tar, t),
          glob(inf$nei, t),
          rep.int(0, length(inf$tar)),
          as.numeric(inf$d)
        )

        out_i <- out_i + 1L
      }
    }

    if (maxt_int == 0L) {

      sol <- do.call(rbind, out)

      ok <- is.finite(sol[, 1]) &
        is.finite(sol[, 2]) &
        sol[, 1] > 0 &
        sol[, 2] > 0

      return(sol[ok, , drop = FALSE])
    }

    ##############################################################
    ## 2) Pure temporal same-site
    ##############################################################

    for (j in seq_len(maxt_int)) {
      for (t in seq_len(numtime - j)) {

        bb <- min(ns[t], ns[t + j])

        if (bb <= 0L) {
          next
        }

        ii <- seq_len(bb)

        out[[out_i]] <- cbind(
          glob(ii, t + j),
          glob(ii, t),
          rep.int(dt(t, j), bb),
          rep.int(0, bb)
        )

        out_i <- out_i + 1L
      }
    }

    ##############################################################
    ## 3) Spatio-temporal cross-time
    ##############################################################

    for (j in seq_len(maxt_int)) {
      for (t in seq_len(numtime - j)) {

        Xn <- getX(t)
        Xt <- getX(t + j)

        lagt_val <- dt(t, j)

        inf <- knn_dir(
          ref = Xn,
          query = Xt,
          K = K,
          drop_self = FALSE,
          x_equals_y = FALSE
        )

        if (length(inf$tar) > 0L) {

          sel <- is.finite(inf$d) & (inf$d > 0)

          if (any(sel)) {

            out[[out_i]] <- cbind(
              glob(inf$tar[sel], t + j),
              glob(inf$nei[sel], t),
              rep.int(lagt_val, sum(sel)),
              as.numeric(inf$d[sel])
            )

            out_i <- out_i + 1L
          }
        }
      }
    }

    sol <- do.call(rbind, out)

    ok <- is.finite(sol[, 1]) &
      is.finite(sol[, 2]) &
      sol[, 1] > 0 &
      sol[, 2] > 0

    sol[ok, , drop = FALSE]
  }

  ##############################################################
  ## Bivariate index
  ##############################################################

  bivariate_index <- function(coords, coordx_dyn = NULL, N, K = 4,
                              maxdist, distance, radius) {

    if (length(K) == 3L) {
      K1 <- K[1]
      K2 <- K[2]
      K3 <- K[3]
    } else {
      K1 <- K2 <- K3 <- K
    }

    if (length(maxdist) == 3L) {
      maxdist1 <- maxdist[1]
      maxdist2 <- maxdist[2]
      maxdist3 <- maxdist[3]
    } else {
      maxdist1 <- maxdist2 <- maxdist3 <- maxdist
    }

    auto_tol <- function(A, B) {

      scaleA <- suppressWarnings(max(abs(range(A, finite = TRUE))))
      scaleB <- suppressWarnings(max(abs(range(B, finite = TRUE))))

      scale <- max(scaleA, scaleB)

      if (!is.finite(scale) || scale <= 0) {
        scale <- 1
      }

      100 * .Machine$double.eps * scale
    }

    ##############################################################
    ## CASE 1: coordx_dyn NULL
    ##############################################################

    if (is.null(coordx_dyn)) {

      n_half <- as.integer(N / 2)
      cm <- coords[1:n_half, , drop = FALSE]

      inf1 <- nn2Geo(
        cm,
        cm,
        K1 + 1,
        distance,
        maxdist1,
        radius,
        drop_first = TRUE,
        x_equals_y = TRUE
      )

      aa1 <- cbind(
        inf1$rowidx,
        inf1$colidx,
        0L,
        0L,
        inf1$lags
      )

      inf2 <- nn2Geo(
        cm,
        cm,
        K3 + 1,
        distance,
        maxdist3,
        radius,
        drop_first = TRUE,
        x_equals_y = TRUE
      )

      aa2 <- cbind(
        inf2$rowidx + n_half,
        inf2$colidx + n_half,
        1L,
        1L,
        inf2$lags
      )

      inf21 <- nn2Geo(
        cm,
        cm,
        K2 + 1,
        distance,
        maxdist2,
        radius,
        drop_first = TRUE,
        x_equals_y = TRUE
      )

      sel <- is.finite(inf21$lags) &
        (inf21$lags > 0) &
        (inf21$colidx > 0) &
        (inf21$rowidx > 0)

      aa_cross <- cbind(
        inf21$rowidx[sel],
        inf21$colidx[sel] + n_half,
        0L,
        1L,
        inf21$lags[sel]
      )

      ii <- seq_len(n_half)

      aa_coll <- cbind(
        ii,
        ii + n_half,
        0L,
        1L,
        0
      )

      return(rbind(
        aa1,
        aa2,
        aa_cross,
        aa_coll
      ))
    }

    ##############################################################
    ## CASE 2: coordx_dyn list length 2
    ##############################################################

    if (!is.list(coordx_dyn) || length(coordx_dyn) != 2L) {
      stop("For bivariate case with coordx_dyn, it must be a list of length 2: list(var1_coords, var2_coords).")
    }

    X1 <- coordx_dyn[[1]]
    X2 <- coordx_dyn[[2]]

    n1 <- nrow(X1)

    inf1 <- nn2Geo(
      X1,
      X1,
      K1 + 1,
      distance,
      maxdist1,
      radius,
      drop_first = TRUE,
      x_equals_y = TRUE
    )

    aa1 <- cbind(
      inf1$rowidx,
      inf1$colidx,
      0L,
      0L,
      inf1$lags
    )

    inf2 <- nn2Geo(
      X2,
      X2,
      K3 + 1,
      distance,
      maxdist3,
      radius,
      drop_first = TRUE,
      x_equals_y = TRUE
    )

    aa2 <- cbind(
      inf2$rowidx + n1,
      inf2$colidx + n1,
      1L,
      1L,
      inf2$lags
    )

    inf_cross <- nn2Geo(
      X2,
      X1,
      K2,
      distance,
      maxdist2,
      radius,
      drop_first = FALSE,
      x_equals_y = FALSE
    )

    sel <- is.finite(inf_cross$lags) &
      (inf_cross$lags > 0) &
      (inf_cross$colidx > 0) &
      (inf_cross$rowidx > 0)

    aa_cross <- cbind(
      inf_cross$rowidx[sel],
      inf_cross$colidx[sel] + n1,
      0L,
      1L,
      inf_cross$lags[sel]
    )

    tol <- auto_tol(X1, X2)

    same_grid <- (nrow(X1) == nrow(X2)) &&
      isTRUE(all(abs(X1 - X2) <= tol))

    if (same_grid) {

      ii <- seq_len(n1)

      aa_coll <- cbind(
        ii,
        ii + n1,
        0L,
        1L,
        0
      )

    } else {

      nn12 <- nn2Geo(
        X2,
        X1,
        K = 1,
        distance,
        maxdist2,
        radius,
        drop_first = FALSE,
        x_equals_y = FALSE
      )

      l12 <- nn12$lags

      l12[is.finite(l12) & (l12 <= tol)] <- 0

      selc <- is.finite(nn12$rowidx) &
        is.finite(nn12$colidx) &
        (nn12$rowidx > 0) &
        (nn12$colidx > 0) &
        is.finite(l12)

      aa_coll <- cbind(
        nn12$rowidx[selc],
        nn12$colidx[selc] + n1,
        0L,
        1L,
        l12[selc]
      )
    }

    rbind(
      aa1,
      aa2,
      aa_cross,
      aa_coll
    )
  }

  ##############################################################
  ## Calibrated Bernoulli probabilities
  ##############################################################

  calibrated_bernoulli_prob <- function(w, p) {

    n <- length(w)

    if (n == 0L) {
      return(numeric(0))
    }

    if (!is.numeric(p) ||
        length(p) != 1L ||
        !is.finite(p) ||
        p <= 0 ||
        p > 1) {
      stop("p_neighb must be in (0,1].")
    }

    if (p >= 1) {
      return(rep.int(1, n))
    }

    w <- as.numeric(w)

    good <- is.finite(w) & w > 0

    if (!any(good)) {
      w[] <- 1
    } else {
      w[!good] <- min(w[good])
    }

    target <- p * n

    f <- function(cst) {
      sum(pmin(1, cst * w)) - target
    }

    hi <- 1 / max(w)

    while (f(hi) < 0) {
      hi <- hi * 2
    }

    cstar <- stats::uniroot(
      f,
      interval = c(0, hi),
      tol = 1e-12
    )$root

    pi <- pmin(1, cstar * w)
    pi[!is.finite(pi)] <- 0

    pmin(1, pmax(0, pi))
  }

  ##############################################################
  ## Subset helper
  ##############################################################

  subset_gb_keep <- function(gb, keep_idx) {

    gb$colidx <- gb$colidx[keep_idx]
    gb$rowidx <- gb$rowidx[keep_idx]
    gb$lags   <- gb$lags[keep_idx]

    if (!is.null(gb$lagt)) {
      gb$lagt <- gb$lagt[keep_idx]
    }

    if (!is.null(gb$first)) {
      gb$first <- gb$first[keep_idx]
    }

    if (!is.null(gb$second)) {
      gb$second <- gb$second[keep_idx]
    }

    gb
  }

  ##############################################################
  ## Method II helper: fixed-budget thinning
  ##
  ## The option name remains "FixedBudget" for compatibility with
  ## the current simulation scripts, but the selection rule implemented
  ## here is fixed-budget thinning: exactly K_target candidate edges are
  ## sampled uniformly without replacement from the full NN candidate set.
  ##
  ## This is the fixed-size analogue of independent Bernoulli thinning:
  ## it preserves the same marginal sampling level approximately, but it
  ## removes the extra randomness in the total number of retained pairs.
  ##############################################################

  fixed_budget_keep <- function(n_edges, K_target) {

    if (n_edges <= 0L || K_target <= 0L) {
      return(integer(0))
    }

    K_target <- as.integer(round(K_target))
    K_target <- max(0L, min(K_target, as.integer(n_edges)))

    if (K_target <= 0L) {
      return(integer(0))
    }

    if (K_target >= n_edges) {
      return(seq_len(n_edges))
    }

    sample.int(n_edges, size = K_target, replace = FALSE)
  }

  ##############################################################
  ## Start
  ##############################################################

  spatial <- TRUE
  spacetime <- FALSE

  if (!is.null(coordx_dyn)) {
    if (!is.list(coordx_dyn)) {
      stop("coordx_dyn must be a list")
    }
  }

  if (!is.null(coordt)) {
    if (is.numeric(coordt) && is.numeric(maxtime)) {
      if (length(coordt) > 1L && length(maxtime) >= 1L) {
        spacetime <- TRUE
      }
    }
  }

  distance <- CheckDistance(distance)
  spatial <- !spacetime && !bivariate

  ##############################################################
  ## Coordinate parsing
  ##############################################################

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

    if (!bivariate) {
      coords <- coordx
      numcoord <- nrow(coords)
    }

    if (bivariate) {
      if (is.null(coordx_dyn)) {
        coords <- coordx
        numcoord <- nrow(coords)
      } else {
        coords <- 1
        numcoord <- 1
      }
    }
  }

  ##############################################################
  ## Build pair indices
  ##############################################################

  if (spatial) {

    K_sp <- as.integer(K[1])

    sol <- nn2Geo(
      coords,
      coords,
      K_sp + 1L,
      distance,
      maxdist,
      radius,
      drop_first = TRUE,
      x_equals_y = TRUE
    )

    gb <- list(
      colidx = sol$colidx,
      rowidx = sol$rowidx,
      lags = sol$lags,
      maxdist = maxdist,
      neighb = neighb
    )
  }

  if (spacetime) {

    K_st <- as.integer(K[1])
    numtime <- length(coordt)

    sol <- spacetime_index(
      coords,
      coordx_dyn,
      numcoord,
      K_st,
      coordt,
      numtime,
      maxtime,
      maxdist,
      distance,
      radius
    )

    gb <- list(
      colidx = sol[, 2],
      rowidx = sol[, 1],
      lags = sol[, 4],
      lagt = sol[, 3]
    )
  }

  if (bivariate) {

    sol <- bivariate_index(
      coords,
      coordx_dyn,
      numcoord,
      K,
      maxdist,
      distance,
      radius
    )

    gb <- list(
      colidx = sol[, 2],
      rowidx = sol[, 1],
      lags = sol[, 5],
      first = sol[, 3],
      second = sol[, 4],
      maxdist = maxdist,
      neighb = neighb
    )
  }

  ##############################################################
  ## p_neighb checks
  ##############################################################

  if (!is.numeric(p_neighb) ||
      length(p_neighb) != 1L ||
      !is.finite(p_neighb)) {
    stop("p_neighb must be a finite numeric scalar.")
  }

  if (p_neighb <= 0 || p_neighb > 1) {
    stop("p_neighb must be in (0,1].")
  }

  thinning_on <- p_neighb < 1

  ##############################################################
  ## Defensive cleaning
  ##############################################################

  if (length(gb$rowidx) > 0L) {

    valid <- (gb$rowidx > 0L) &
      (gb$colidx > 0L) &
      is.finite(gb$rowidx) &
      is.finite(gb$colidx)

    if (!all(valid)) {

      gb$rowidx <- gb$rowidx[valid]
      gb$colidx <- gb$colidx[valid]
      gb$lags   <- gb$lags[valid]

      if (!is.null(gb$lagt)) {
        gb$lagt <- gb$lagt[valid]
      }

      if (!is.null(gb$first)) {
        gb$first <- gb$first[valid]
      }

      if (!is.null(gb$second)) {
        gb$second <- gb$second[valid]
      }
    }
  }

  ##############################################################
  ## No thinning if p_neighb = 1
  ##############################################################

  if (!thinning_on) {
    return(gb)
  }

  ##############################################################
  ## Method I: calibrated Bernoulli thinning
  ## thin_method = "bernoulli"
  ##############################################################

  if (thin_method == "Bernoulli") {

    n_all <- length(gb$rowidx)

    if (n_all == 0L) {
      return(gb)
    }

    w <- rep.int(1, n_all)

    p_ij <- calibrated_bernoulli_prob(
      w = w,
      p = p_neighb
    )

    B <- stats::rbinom(
      n_all,
      size = 1L,
      prob = p_ij
    )

    keep_idx <- which(B == 1L)

    gb <- subset_gb_keep(gb, keep_idx)

    return(gb)
  }

  ##############################################################
  ## Method II: fixed-budget thinning
  ##
  ## The retained set has exact size
  ## K_target = round(p_neighb * number_of_candidate_edges).
  ##############################################################

  if (thin_method == "FixedBudget") {

    n_all <- length(gb$rowidx)

    if (n_all == 0L) {
      return(gb)
    }

    K_target <- as.integer(round(p_neighb * n_all))
    K_target <- max(1L, min(K_target, n_all))

    keep_idx <- fixed_budget_keep(
      n_edges = n_all,
      K_target = K_target
    )

    gb <- subset_gb_keep(gb, keep_idx)

    return(gb)
  }

  return(gb)
}