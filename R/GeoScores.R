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

  y <- as.numeric(c(data_to_pred))

  if (!is.null(probject) && (inherits(probject, "GeoKrig") || inherits(probject, "GeoKrigloc"))) {
    pred <- as.numeric(c(probject$pred))
    mse <- as.numeric(c(probject$mse))
  } else {
    if (is.null(pred)) {
      stop("GeoKrig object or pred is mandatory")
    }
    pred <- as.numeric(c(pred))
  }

  if (length(y) != length(pred)) {
    stop("length of data_to_pred and pred does not match")
  }

  has_interval <- !is.null(lower95) && !is.null(upper95)
  needs_distribution <- any(score %in% c("brie", "crps", "lscore", "pit")) ||
    (!has_interval && any(score %in% c("intscore", "coverage")))
  needs_interval <- any(score %in% c("intscore", "coverage"))

  if (has_interval) {
    lower95 <- as.numeric(c(lower95))
    upper95 <- as.numeric(c(upper95))
    if (length(lower95) != length(pred) || length(upper95) != length(pred)) {
      stop("length of lower95/upper95 and pred does not match")
    }
  } else {
    lower95 <- rep(NA_real_, length(pred))
    upper95 <- rep(NA_real_, length(pred))
  }

  if (!is.null(mse)) {
    mse <- as.numeric(c(mse))
    if (length(mse) != length(pred)) {
      stop("length of mse and pred does not match")
    }
    se <- sqrt(mse)
  } else if (has_interval) {
    # Approximation used only for normal-distribution scores when mse is absent.
    # Interval scores and coverage below use lower95/upper95 directly.
    se <- (upper95 - lower95) / (2 * qnorm(0.975))
    mse <- se^2
  } else if (needs_distribution || needs_interval) {
    stop("mse or both lower95 and upper95 are mandatory for probabilistic scores")
  } else {
    se <- rep(NA_real_, length(pred))
    mse <- rep(NA_real_, length(pred))
  }

  ok <- is.finite(y) & is.finite(pred)

  if (needs_distribution || (!has_interval && needs_interval)) {
    ok <- ok & is.finite(se) & se > 0 & is.finite(mse) & mse > 0
  }

  if (has_interval && needs_interval) {
    ok <- ok & is.finite(lower95) & is.finite(upper95) & lower95 < upper95
  }

  if (na.rm) {
    y <- y[ok]
    pred <- pred[ok]
    se <- se[ok]
    mse <- mse[ok]
    lower95 <- lower95[ok]
    upper95 <- upper95[ok]
  } else if (!all(ok)) {
    stop("non-finite values, invalid intervals, or non-positive predictive variances/standard errors found")
  }

  if (length(y) == 0) {
    stop("no valid observations available after removing missing/non-finite values")
  }

  err <- y - pred
  out <- list()

  if ("pe" %in% score) {
    out$MAE <- mean(abs(err))
    out$RMSPE <- sqrt(mean(err^2))
  }

  if (needs_distribution || needs_interval) {
    std <- err / se

    if ("brie" %in% score) {
      probs <- 1 - pnorm(threshold, mean = pred, sd = se)
      obs <- as.numeric(y > threshold)
      out$Brier <- mean((probs - obs)^2)
    }

    if ("lscore" %in% score) {
      out$LogScore <- mean(0.5 * (log(2 * pi * mse) + std^2))
    }

    if ("pit" %in% score) {
      if (!is.null(probject) &&
          (inherits(probject, "GeoKrig") || inherits(probject, "GeoKrigloc")) &&
          exists("GeoPit", mode = "function")) {
        probject$data_to_pred <- y
        probject$pred <- pred
        probject$mse <- mse
        out$PIT <- GeoPit(probject, type = "Uniform")$data
      } else {
        out$PIT <- pnorm(y, mean = pred, sd = se)
      }
    }

    if ("crps" %in% score) {
      out$CRPS <- mean(se * (std * (2 * pnorm(std) - 1) +
                               2 * dnorm(std) - 1 / sqrt(pi)))
    }

    if ("intscore" %in% score || "coverage" %in% score) {
      alpha <- 0.05

      if (has_interval) {
        l <- lower95
        u <- upper95
      } else {
        l <- pred + se * qnorm(alpha / 2)
        u <- pred + se * qnorm(1 - alpha / 2)
      }

      if ("intscore" %in% score) {
        out$IS95 <- mean((u - l) +
                           (2 / alpha) * (l - y) * (y < l) +
                           (2 / alpha) * (y - u) * (y > u))
      }

      if ("coverage" %in% score) {
        out$Cvg95 <- mean(y >= l & y <= u)
      }
    }
  }

  out
}
