####################################################
### File name: GeoTest.r
####################################################

### Procedures are in alphabetical order.

### Statistical hypothesis testing for nested models
GeoTests <- function(object1, object2, ..., statistic) {

  ## ------------------------------------------------------------------
  ## Internal utilities
  ## ------------------------------------------------------------------
  valid_statistics <- c("Wald", "Wilks", "WilksS")
  composite_likelihoods <- c("Marginal", "Conditional", "Difference")

  same_value <- function(x, y) {
    isTRUE(all.equal(x, y, check.attributes = FALSE))
  }

  as_named_numeric <- function(x) {
    if (is.null(x)) return(setNames(numeric(0), character(0)))
    unlist(x, use.names = TRUE)
  }

  ensure_named_matrix <- function(A, par_names, label) {
    if (is.null(A)) {
      stop(label, " is missing", call. = FALSE)
    }

    A <- as.matrix(A)
    storage.mode(A) <- "double"

    if (!all(dim(A) == c(length(par_names), length(par_names)))) {
      stop(label, " has incompatible dimensions", call. = FALSE)
    }

    if (is.null(rownames(A)) || is.null(colnames(A))) {
      rownames(A) <- colnames(A) <- par_names
    }

    if (!all(par_names %in% rownames(A)) ||
        !all(par_names %in% colnames(A))) {
      stop(label, " does not contain all estimated parameters", call. = FALSE)
    }

    A <- A[par_names, par_names, drop = FALSE]
    (A + t(A)) / 2
  }

  safe_solve <- function(A, label) {
    A <- as.matrix(A)
    storage.mode(A) <- "double"
    A <- (A + t(A)) / 2

    out <- tryCatch(
      solve(A),
      error = function(e) NULL
    )

    if (is.null(out) || any(!is.finite(out))) {
      stop(label, " is singular or numerically non-invertible", call. = FALSE)
    }

    (out + t(out)) / 2
  }

  finite_loglik <- function(model, label) {
    value <- if (!is.null(model$logCompLik)) {
      model$logCompLik
    } else {
      model$logLik
    }

    value <- as.numeric(value)

    if (length(value) != 1L || !is.finite(value)) {
      stop(label, " does not contain a finite log-likelihood value", call. = FALSE)
    }

    value
  }

  is_stochastic_thinning <- function(model) {
    p <- if (is.null(model$p_neighb)) 1 else as.numeric(model$p_neighb)
    method <- if (is.null(model$thin_method)) {
      "bernoulli"
    } else {
      tolower(as.character(model$thin_method))
    }

    method %in% c(
      "bernoulli", "fixedbudget", "fixed_budget",
      "targetbalanced", "match"
    ) && isTRUE(p < 1)
  }

  canonical_coordt <- function(model) {
    x <- model$coordt

    ## For purely spatial and bivariate fits, GeoModels may represent the
    ## absence of a time coordinate as NULL, numeric(0), or a scalar
    ## placeholder such as 0. These representations are equivalent.
    if (!isTRUE(model$spacetime) &&
        (is.null(x) || length(x) <= 1L)) {
      return(NULL)
    }

    unname(x)
  }

  check_common_fit_structure <- function(model1, model2) {
    fields <- c(
      "model", "corrmodel", "likelihood", "type", "copula",
      "distance", "grid", "n", "neighb", "maxdist", "maxtime",
      "p_neighb", "thin_method", "weighted"
    )

    bad <- fields[!vapply(
      fields,
      function(nm) same_value(model1[[nm]], model2[[nm]]),
      logical(1)
    )]

    if (length(bad) > 0L) {
      stop(
        "The fitted models are not directly comparable; different fields: ",
        paste(bad, collapse = ", "),
        call. = FALSE
      )
    }

    data_fields <- c(
      "data", "coordx", "coordy", "coordz",
      "coordx_dyn", "X"
    )

    bad_data <- data_fields[!vapply(
      data_fields,
      function(nm) same_value(model1[[nm]], model2[[nm]]),
      logical(1)
    )]

    if (!same_value(
      canonical_coordt(model1),
      canonical_coordt(model2)
    )) {
      bad_data <- c(bad_data, "coordt")
    }

    if (length(bad_data) > 0L) {
      stop(
        "The fitted models must use the same data, coordinates and design matrix; ",
        "different fields: ", paste(bad_data, collapse = ", "),
        call. = FALSE
      )
    }
  }

  nesting_information <- function(model1, model2) {
    par1 <- as_named_numeric(model1$param)
    par2 <- as_named_numeric(model2$param)
    fixed1 <- as_named_numeric(model1$fixed)
    fixed2 <- as_named_numeric(model2$fixed)

    if (is.null(names(par1)) || anyDuplicated(names(par1))) {
      stop("The unrestricted model has invalid parameter names", call. = FALSE)
    }

    if (is.null(names(par2)) || anyDuplicated(names(par2))) {
      stop("The restricted model has invalid parameter names", call. = FALSE)
    }

    if (!all(names(par2) %in% names(par1))) {
      stop(
        "Models are not nested: the restricted model estimates parameters ",
        "not estimated by the unrestricted model",
        call. = FALSE
      )
    }

    tested <- setdiff(names(par1), names(par2))

    if (length(tested) < 1L) {
      stop("Models are not nested or have the same number of parameters", call. = FALSE)
    }

    if (!all(tested %in% names(fixed2))) {
      missing_null <- setdiff(tested, names(fixed2))
      stop(
        "The restricted model must fix every tested parameter. Missing: ",
        paste(missing_null, collapse = ", "),
        call. = FALSE
      )
    }

    null_values <- fixed2[tested]

    if (any(!is.finite(null_values))) {
      stop("The null parameter values must be finite", call. = FALSE)
    }

    common_fixed <- setdiff(intersect(names(fixed1), names(fixed2)), tested)

    if (length(common_fixed) > 0L &&
        !same_value(fixed1[common_fixed], fixed2[common_fixed])) {
      stop(
        "Common fixed parameters must have identical values in the two models",
        call. = FALSE
      )
    }

    list(
      par1 = par1,
      par2 = par2,
      tested = tested,
      null_values = null_values,
      df = length(tested)
    )
  }

  check_same_composite_pairs <- function(model1, model2) {
    have_pairs1 <- !is.null(model1$rowidx) && !is.null(model1$colidx)
    have_pairs2 <- !is.null(model2$rowidx) && !is.null(model2$colidx)

    if (have_pairs1 && have_pairs2) {
      same_pairs <- identical(model1$rowidx, model2$rowidx) &&
        identical(model1$colidx, model2$colidx)

      if (!same_pairs) {
        stop(
          "Composite likelihood-ratio tests require the same retained pairs ",
          "in the unrestricted and restricted fits",
          call. = FALSE
        )
      }

      return(invisible(TRUE))
    }

    if (is_stochastic_thinning(model1) || is_stochastic_thinning(model2)) {
      stop(
        "The fits use stochastic thinning, but the retained pair indices are ",
        "not available to verify that the same pairs were used",
        call. = FALSE
      )
    }

    invisible(TRUE)
  }

  compute_statistic <- function(model1, model2, statistic, nesting) {
    likelihood <- as.character(model1$likelihood)
    is_composite <- likelihood %in% composite_likelihoods
    is_full <- identical(likelihood, "Full")

    if (!is_composite && !is_full) {
      stop("Unsupported likelihood type: ", likelihood, call. = FALSE)
    }

    if (statistic == "Wilks" && !is_full) {
      stop("'Wilks' is available only for full likelihood fits", call. = FALSE)
    }

    if (statistic == "WilksS" && !is_composite) {
      stop("'WilksS' is available only for composite likelihood fits", call. = FALSE)
    }

    tested <- nesting$tested
    theta <- nesting$par1[tested] - nesting$null_values
    df <- nesting$df

    if (statistic == "Wald") {
      V <- ensure_named_matrix(
        model1$varcov,
        names(nesting$par1),
        "The covariance matrix of the unrestricted model"
      )

      V_test <- V[tested, tested, drop = FALSE]
      V_test_inv <- safe_solve(V_test, "The covariance matrix for the tested parameters")

      W <- as.numeric(crossprod(theta, V_test_inv %*% theta))
      nu <- df
    }

    if (statistic == "Wilks") {
      W <- 2 * (
        finite_loglik(model1, "The unrestricted model") -
          finite_loglik(model2, "The restricted model")
      )
      nu <- df
    }

    if (statistic == "WilksS") {
      check_same_composite_pairs(model1, model2)

      W_raw <- 2 * (
        finite_loglik(model1, "The unrestricted model") -
          finite_loglik(model2, "The restricted model")
      )

      H <- ensure_named_matrix(
        model1$sensmat,
        names(nesting$par1),
        "The sensitivity matrix of the unrestricted model"
      )

      V <- ensure_named_matrix(
        model1$varcov,
        names(nesting$par1),
        "The Godambe covariance matrix of the unrestricted model"
      )

      ## Efficient sensitivity for the tested parameter block:
      ## H_eff = { (H^{-1})_{tested,tested} }^{-1}.
      H_inv <- safe_solve(H, "The sensitivity matrix")
      H_eff <- safe_solve(
        H_inv[tested, tested, drop = FALSE],
        "The efficient sensitivity matrix for the tested parameters"
      )

      V_test <- V[tested, tested, drop = FALSE]

      ## H_eff %*% V_test has positive real eigenvalues.  The symmetric
      ## similar matrix below is numerically more stable.
      R <- tryCatch(chol(H_eff), error = function(e) NULL)

      if (is.null(R)) {
        stop(
          "The efficient sensitivity matrix is not positive definite",
          call. = FALSE
        )
      }

      lambda_mat <- R %*% V_test %*% t(R)
      lambda_mat <- (lambda_mat + t(lambda_mat)) / 2
      lambda <- eigen(lambda_mat, symmetric = TRUE, only.values = TRUE)$values

      tol <- max(1, max(abs(lambda))) * sqrt(.Machine$double.eps)

      if (any(lambda < -tol) || any(!is.finite(lambda))) {
        stop("Invalid eigenvalues in the WilksS adjustment", call. = FALSE)
      }

      lambda[lambda < 0] <- 0

      if (!all(lambda > 0)) {
        stop(
          "Non-positive eigenvalues in the WilksS adjustment; ",
          "the tested parameters may be weakly identified",
          call. = FALSE
        )
      }

      sum_lambda <- sum(lambda)
      sum_lambda2 <- sum(lambda^2)

      nu <- sum_lambda^2 / sum_lambda2
      scale_factor <- sum_lambda2 / sum_lambda
      W <- W_raw / scale_factor
    }

    tol_W <- 100 * .Machine$double.eps *
      max(1, abs(finite_loglik(model1, "The unrestricted model")))

    if (W < -tol_W) {
      stop(
        "The test statistic is negative. Check optimizer convergence and ",
        "the ordering of unrestricted and restricted models",
        call. = FALSE
      )
    }

    W <- max(0, as.numeric(W))

    if (!is.finite(W) || !is.finite(nu) || nu <= 0) {
      stop("The test statistic or its degrees of freedom is invalid", call. = FALSE)
    }

    list(W = W, nu = as.numeric(nu))
  }

  ## ------------------------------------------------------------------
  ## Input checks and model collection
  ## ------------------------------------------------------------------
  if (missing(object1) || missing(object2)) {
    stop("Models one and two must be specified", call. = FALSE)
  }

  if (missing(statistic) || !is.character(statistic) || length(statistic) != 1L) {
    stop("statistic must be one of: Wald, Wilks, WilksS", call. = FALSE)
  }

  if (!statistic %in% valid_statistics) {
    stop(
      "Unknown statistic. Available tests are: ",
      paste(valid_statistics, collapse = ", "),
      call. = FALSE
    )
  }

  model_calls <- c(
    list(substitute(object1), substitute(object2)),
    as.list(substitute(list(...)))[-1L]
  )

  model_names <- vapply(
    model_calls,
    function(x) paste(deparse(x), collapse = ""),
    character(1)
  )

  models <- lapply(model_calls, eval, envir = parent.frame())

  if (!all(vapply(models, inherits, logical(1), what = "GeoFit"))) {
    stop("GeoTests can be used only with 'GeoFit' objects", call. = FALSE)
  }

  nummod <- length(models)
  numparam <- vapply(models, function(x) length(as_named_numeric(x$param)), integer(1))
  df <- nu <- W <- pvalue <- rep(NA_real_, nummod - 1L)

  ## ------------------------------------------------------------------
  ## Adjacent nested-model comparisons
  ## ------------------------------------------------------------------
  for (i in 2:nummod) {
    unrestricted <- models[[i - 1L]]
    restricted <- models[[i]]

    check_common_fit_structure(unrestricted, restricted)
    nesting <- nesting_information(unrestricted, restricted)

    stat <- compute_statistic(
      unrestricted,
      restricted,
      statistic,
      nesting
    )

    j <- i - 1L
    df[j] <- nesting$df
    nu[j] <- stat$nu
    W[j] <- stat$W
    pvalue[j] <- stats::pchisq(W[j], df = nu[j], lower.tail = FALSE)
  }

  ## ------------------------------------------------------------------
  ## Output table
  ## ------------------------------------------------------------------
  table <- data.frame(
    "Num.Par" = numparam,
    "Diff.Par" = c(NA_real_, df),
    "Df" = c(NA_real_, nu),
    "Chisq" = c(NA_real_, W),
    "Pr(>chisq)" = c(NA_real_, pvalue),
    check.names = FALSE
  )

  rownames(table) <- model_names

  structure(
    table,
    heading = "Statistical Hypothesis Test Table\n",
    class = c("data.frame")
  )
}