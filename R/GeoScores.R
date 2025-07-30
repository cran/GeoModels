GeoScores <- function(data_to_pred, probject = NULL, pred = NULL, mse = NULL,
                      score = c("brie", "crps", "lscore", "pit", "pe", "intscore", "coverage")) {

  if (is.null(probject)) {
    if (is.null(pred) || is.null(mse)) stop("Geokrig object or (pred and mse) are mandatory")
  }

  data_to_pred <- c(data_to_pred)


  if (inherits(probject, "GeoKrig") || inherits(probject, "GeoKrigloc")) {
    pred <- c(probject$pred)
    mse  <- c(probject$mse)
  }

  stopifnot(is.numeric(pred))
  stopifnot(is.numeric(data_to_pred))
  if (!is.null(mse)) stopifnot(is.numeric(mse))

  N1 <- length(pred)
  N2 <- length(data_to_pred)
  if (N1 != N2) stop("length of data and predictions does not match")

  if (!is.null(mse)) sqrtvv <- sqrt(mse)

  # Inizializza
  rmse <- mae <- mad <- lscore <- crps <- pit <- brie <- intscore <- coverage <- NULL

  # PIT
  if (!is.null(mse) && "pit" %in% score) {
    if (!is.null(probject)) {
      probject$data_to_pred <- data_to_pred
      pit <- GeoPit(probject, type = "Uniform")$data
    } else {
      pit <- pnorm(data_to_pred, mean = pred, sd = sqrtvv)
    }
  }

  # Brier Score
  if (!is.null(mse) && "brie" %in% score) {
      threshold = 0.5
    probs <- 1 - pnorm(threshold, mean = pred, sd = sqrtvv)
    obs <- as.numeric(data_to_pred > threshold)
    brie <- mean((probs - obs)^2)
  }

  # Prediction Errors
  if ("pe" %in% score) {
    err <- data_to_pred - pred
    rmse <- sqrt(mean(err^2))
    mae  <- mean(abs(err))
    mad  <- median(abs(err))
  }

  # CRPS & LogScore
  if (!is.null(mse)) {
    err <- data_to_pred - pred
    std <- err / sqrtvv

    if ("lscore" %in% score) {

      lscore <- mean(0.5 * (log(2 * pi * mse) + std^2))
    }

    if ("crps" %in% score) {
      crps <- mean(sqrtvv * (std * (2 * pnorm(std) - 1) + 2 * dnorm(std) - 1 / sqrt(pi)))
    }
  }

  # Interval Score
  if (!is.null(mse) && "intscore" %in% score) {
    intervalScore <- function(Z, l, u, a = 0.1) {
      (u - l) + (2 / a) * (l - Z) * (Z < l) + (2 / a) * (Z - u) * (Z > u)
    }
    l <- pred + sqrtvv * qnorm(0.025)
    u <- pred + sqrtvv * qnorm(0.975)
    intscore <- mean(intervalScore(Z = data_to_pred, l = l, u = u))
  }

  # Coverage
  if (!is.null(mse) && "coverage" %in% score) {
    l <- pred + sqrtvv * qnorm(0.025)
    u <- pred + sqrtvv * qnorm(0.975)
    coverage <- mean(data_to_pred >= l & data_to_pred <= u)
  }

  # Output
  results <- list(
    brie = brie,
    rmse = rmse,
    mae = mae,
    mad = mad,
    lscore = lscore,
    crps = crps,
    pit = pit,
    intscore = intscore,
    coverage = coverage
  )

  return(results[!sapply(results, is.null)])
}
