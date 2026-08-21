GeoScores <- function(data_to_pred,
                      probject = NULL,
                      pred = NULL,
                      mse = NULL,
                      score = c("pe", "crps", "intscore", "coverage"),
                      lower95 = NULL,
                      upper95 = NULL,
                      threshold = 0.5,
                      na.rm = TRUE) {

  score <- match.arg(score,
                     choices = c("brie", "crps", "lscore", "pit", "pe", "intscore", "coverage"),
                     several.ok = TRUE)

  if (!is.logical(na.rm) || length(na.rm) != 1L || is.na(na.rm)) {
    stop("na.rm must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold)) {
    stop("threshold must be a finite numeric scalar", call. = FALSE)
  }

  y <- as.numeric(c(data_to_pred))
  from_krig <- !is.null(probject) &&
    (inherits(probject, "GeoKrig") || inherits(probject, "GeoKrigloc"))
  from_simcond <- !is.null(probject) && inherits(probject, "GeoSimcond")

  if (!is.null(probject) && !from_krig && !from_simcond) {
    stop("probject must be an object of class GeoKrig, GeoKrigloc, or GeoSimcond", call. = FALSE)
  }

  if (xor(is.null(lower95), is.null(upper95))) {
    stop("lower95 and upper95 must be supplied together", call. = FALSE)
  }
  has_interval <- !is.null(lower95) && !is.null(upper95)

  check_mask <- function(mask, what) {
    if (!na.rm && !all(mask)) {
      stop("invalid or non-finite values found for ", what, call. = FALSE)
    }
    if (!any(mask)) {
      stop("no valid observations available for ", what, call. = FALSE)
    }
    mask
  }

  ## -----------------------------------------------------------------------
  ## Conditional-simulation objects: use the empirical predictive
  ## distribution represented by condsim. No Gaussian approximation is used.
  ## -----------------------------------------------------------------------
  if (from_simcond) {
    cs <- probject$condsim
    if (is.null(cs)) {
      stop("GeoSimcond object does not contain conditional simulations", call. = FALSE)
    }

    if (is.list(cs)) {
      if (length(cs) < 1L) {
        stop("GeoSimcond object does not contain conditional simulations", call. = FALSE)
      }
      lens <- vapply(cs, function(z) length(c(z)), integer(1L))
      if (any(lens != lens[1L])) {
        stop("conditional simulations have inconsistent lengths", call. = FALSE)
      }
      sim <- do.call(rbind, lapply(cs, function(z) as.numeric(c(z))))
    } else {
      sim <- as.matrix(cs)
      storage.mode(sim) <- "double"
    }

    if (length(dim(sim)) != 2L || nrow(sim) < 1L || ncol(sim) < 1L) {
      stop("GeoSimcond conditional simulations have an invalid shape", call. = FALSE)
    }

    ## GeoSimcond stores replications by rows and prediction locations by
    ## columns. Accept the transposed orientation defensively for external
    ## objects that preserve the same class contract.
    if (ncol(sim) != length(y) && nrow(sim) == length(y)) sim <- t(sim)
    if (ncol(sim) != length(y)) {
      stop("length of data_to_pred and GeoSimcond prediction locations does not match", call. = FALSE)
    }

    n <- ncol(sim)
    B <- nrow(sim)
    finite_draw <- is.finite(sim)
    if (!na.rm && !all(finite_draw)) {
      stop("invalid or non-finite conditional simulations found", call. = FALSE)
    }
    draw_count <- colSums(finite_draw)
    draw_ok <- draw_count > 0L

    if (!is.null(probject$cond_mean) && length(c(probject$cond_mean)) == n) {
      point_pred <- as.numeric(c(probject$cond_mean))
    } else {
      point_pred <- vapply(seq_len(n), function(j) {
        z <- sim[, j]
        z <- z[is.finite(z)]
        if (length(z)) mean(z) else NA_real_
      }, numeric(1L))
    }

    if (has_interval) {
      lower95 <- as.numeric(c(lower95))
      upper95 <- as.numeric(c(upper95))
      if (length(lower95) != n || length(upper95) != n) {
        stop("length of lower95/upper95 and prediction locations does not match", call. = FALSE)
      }
    }

    out <- list()

    if ("pe" %in% score) {
      idx <- check_mask(is.finite(y) & is.finite(point_pred), "point prediction scores")
      err <- y[idx] - point_pred[idx]
      out$MAE <- mean(abs(err))
      out$RMSPE <- sqrt(mean(err^2))
    }

    ## Empirical CRPS in O(B log B) per prediction location:
    ## E|X-y| - 1/(2)E|X-X'| for the empirical predictive distribution.
    if ("crps" %in% score) {
      crps_i <- rep(NA_real_, n)
      valid_y <- is.finite(y) & draw_ok
      if (!na.rm && !all(valid_y)) {
        stop("invalid values found for conditional-simulation CRPS", call. = FALSE)
      }
      for (j in which(valid_y)) {
        z <- sim[, j]
        z <- sort(z[is.finite(z)])
        m <- length(z)
        coeff <- 2 * seq_len(m) - m - 1
        crps_i[j] <- mean(abs(z - y[j])) - sum(coeff * z) / (m * m)
      }
      idx <- check_mask(is.finite(crps_i), "conditional-simulation CRPS")
      out$CRPS <- mean(crps_i[idx])
    }

    if ("pit" %in% score) {
      pit <- rep(NA_real_, n)
      valid_y <- is.finite(y) & draw_ok
      if (!na.rm && !all(valid_y)) {
        stop("invalid values found for conditional-simulation PIT", call. = FALSE)
      }
      for (j in which(valid_y)) {
        z <- sim[, j]
        z <- z[is.finite(z)]
        pit[j] <- mean(z <= y[j])
      }
      if (!any(is.finite(pit))) {
        stop("no valid observations available for conditional-simulation PIT", call. = FALSE)
      }
      out$PIT <- pit
    }

    if ("brie" %in% score) {
      brier_i <- rep(NA_real_, n)
      valid_y <- is.finite(y) & draw_ok
      if (!na.rm && !all(valid_y)) {
        stop("invalid values found for conditional-simulation Brier score", call. = FALSE)
      }
      for (j in which(valid_y)) {
        z <- sim[, j]
        z <- z[is.finite(z)]
        prob <- mean(z > threshold)
        obs <- as.numeric(y[j] > threshold)
        brier_i[j] <- (prob - obs)^2
      }
      idx <- check_mask(is.finite(brier_i), "conditional-simulation Brier score")
      out$Brier <- mean(brier_i[idx])
    }

    if ("intscore" %in% score || "coverage" %in% score) {
      if (has_interval) {
        l <- lower95
        u <- upper95
        interval_ok <- is.finite(y) & is.finite(l) & is.finite(u) & l <= u
      } else {
        l <- rep(NA_real_, n)
        u <- rep(NA_real_, n)
        for (j in which(draw_ok)) {
          z <- sim[, j]
          z <- z[is.finite(z)]
          qq <- as.numeric(stats::quantile(z, probs = c(0.025, 0.975),
                                           names = FALSE, type = 7))
          l[j] <- qq[1L]
          u[j] <- qq[2L]
        }
        interval_ok <- is.finite(y) & is.finite(l) & is.finite(u) & l <= u
      }

      idx <- check_mask(interval_ok, "conditional-simulation prediction-interval scores")
      yi <- y[idx]
      li <- l[idx]
      ui <- u[idx]
      alpha <- 0.05
      if ("intscore" %in% score) {
        out$IS95 <- mean((ui - li) +
                           (2 / alpha) * (li - yi) * (yi < li) +
                           (2 / alpha) * (yi - ui) * (yi > ui))
      }
      if ("coverage" %in% score) {
        out$Cvg95 <- mean(yi >= li & yi <= ui)
      }
    }

    if ("lscore" %in% score) {
      out$LogScore <- NA_real_
      warning(
        "LogScore is not computed automatically from a GeoSimcond object because conditional draws do not by themselves provide the predictive density at the observed value. Other requested simulation-based scores were computed.",
        call. = FALSE
      )
    }

    return(out)
  }

  ## -----------------------------------------------------------------------
  ## GeoKrig/GeoKrigloc or direct pred/mse input. Distributional scores are
  ## always computed from the Gaussian N(pred, mse) convention when requested.
  ## For a non-Gaussian GeoKrig object this is a Gaussian predictive
  ## approximation based on the OLP and its MSE; interpretation is left to the
  ## user. For a Gaussian model it is the natural conditional distribution.
  ## -----------------------------------------------------------------------
  if (from_krig) {
    pred <- as.numeric(c(probject$pred))
    mse <- if (is.null(probject$mse)) NULL else as.numeric(c(probject$mse))
  } else {
    if (is.null(pred)) {
      stop("GeoKrig/GeoKrigloc/GeoSimcond object or pred is mandatory", call. = FALSE)
    }
    pred <- as.numeric(c(pred))
    if (!is.null(mse)) mse <- as.numeric(c(mse))
  }

  if (length(y) != length(pred)) {
    stop("length of data_to_pred and pred does not match", call. = FALSE)
  }
  n <- length(pred)

  if (has_interval) {
    lower95 <- as.numeric(c(lower95))
    upper95 <- as.numeric(c(upper95))
    if (length(lower95) != n || length(upper95) != n) {
      stop("length of lower95/upper95 and pred does not match", call. = FALSE)
    }
  }

  point_requested <- "pe" %in% score
  distribution_scores <- intersect(score, c("brie", "crps", "lscore", "pit"))
  interval_scores <- intersect(score, c("intscore", "coverage"))
  needs_distribution <- length(distribution_scores) > 0L
  needs_interval <- length(interval_scores) > 0L
  needs_gaussian_distribution <- needs_distribution || (needs_interval && !has_interval)

  se <- NULL
  if (needs_gaussian_distribution) {
    if (!is.null(mse)) {
      if (length(mse) != n) stop("length of mse and pred does not match", call. = FALSE)
      se <- suppressWarnings(sqrt(mse))
    } else if (has_interval) {
      se <- (upper95 - lower95) / (2 * qnorm(0.975))
      mse <- se^2
    } else {
      stop("mse or both lower95 and upper95 are mandatory for Gaussian probabilistic scores", call. = FALSE)
    }
  }

  point_ok <- is.finite(y) & is.finite(pred)
  distribution_ok <- point_ok
  if (needs_gaussian_distribution) {
    distribution_ok <- distribution_ok & is.finite(se) & se > 0 &
      is.finite(mse) & mse > 0
  }
  interval_ok <- point_ok
  if (needs_interval && has_interval) {
    interval_ok <- interval_ok & is.finite(lower95) & is.finite(upper95) &
      lower95 < upper95
  } else if (needs_interval) {
    interval_ok <- distribution_ok
  }

  out <- list()

  if (point_requested) {
    idx <- check_mask(point_ok, "point prediction scores")
    err <- y[idx] - pred[idx]
    out$MAE <- mean(abs(err))
    out$RMSPE <- sqrt(mean(err^2))
  }

  if (needs_distribution) {
    idx <- check_mask(distribution_ok, "Gaussian predictive-distribution scores")
    yd <- y[idx]
    pd <- pred[idx]
    sed <- se[idx]
    msed <- mse[idx]
    errd <- yd - pd
    std <- errd / sed

    if ("brie" %in% score) {
      probs <- 1 - pnorm(threshold, mean = pd, sd = sed)
      obs <- as.numeric(yd > threshold)
      out$Brier <- mean((probs - obs)^2)
    }

    if ("lscore" %in% score) {
      out$LogScore <- mean(0.5 * (log(2 * pi * msed) + std^2))
    }

    if ("pit" %in% score) {
      out$PIT <- pnorm(yd, mean = pd, sd = sed)
    }

    if ("crps" %in% score) {
      out$CRPS <- mean(
        sed * (std * (2 * pnorm(std) - 1) +
                 2 * dnorm(std) - 1 / sqrt(pi))
      )
    }
  }

  if (needs_interval) {
    idx <- check_mask(interval_ok, "prediction-interval scores")
    yi <- y[idx]

    if (has_interval) {
      l <- lower95[idx]
      u <- upper95[idx]
    } else {
      sei <- se[idx]
      pi <- pred[idx]
      alpha <- 0.05
      l <- pi + sei * qnorm(alpha / 2)
      u <- pi + sei * qnorm(1 - alpha / 2)
    }

    alpha <- 0.05
    if ("intscore" %in% score) {
      out$IS95 <- mean((u - l) +
                         (2 / alpha) * (l - yi) * (yi < l) +
                         (2 / alpha) * (yi - u) * (yi > u))
    }
    if ("coverage" %in% score) {
      out$Cvg95 <- mean(yi >= l & yi <= u)
    }
  }

  out
}
