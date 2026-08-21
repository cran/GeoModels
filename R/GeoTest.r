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
  composite_likelihoods <- c("Marginal", "Conditional")

  same_value <- function(x, y) {
    isTRUE(all.equal(x, y, check.attributes = FALSE))
  }

  as_named_numeric <- function(x) {
    if (is.null(x)) {
      return(setNames(numeric(0), character(0)))
    }

    out <- unlist(x, use.names = TRUE)
    storage.mode(out) <- "double"
    out
  }

  canonical_model <- function(x) {
    x <- as.character(x)[1L]

    if (x %in% c("Gauss", "Gaussian")) {
      return("Gaussian")
    }
    if (x %in% c("SkewGauss", "SkewGaussian")) {
      return("SkewGaussian")
    }

    x
  }

  marginal_nesting_rule <- function(model1, model2) {
    unrestricted <- canonical_model(model1$model)
    restricted <- canonical_model(model2$model)

    if (identical(unrestricted, restricted)) {
      return(list(
        same_model = TRUE,
        null = setNames(numeric(0), character(0)),
        boundary = character(0)
      ))
    }

    rules <- list(
      "SkewGaussian->Gaussian" = list(
        null = c(skew = 0),
        boundary = character(0)
      ),
      "Tukeyh->Gaussian" = list(
        null = c(tail = 0),
        boundary = "tail = 0 in the Tukey-h model"
      ),
      "SinhAsinh->Gaussian" = list(
        null = c(skew = 0, tail = 1),
        boundary = character(0)
      )
    )

    key <- paste0(unrestricted, "->", restricted)
    rule <- rules[[key]]

    if (is.null(rule)) {
      stop(
        "The marginal models are not identical and no exact nesting relation ",
        "is implemented for ", unrestricted, " -> ", restricted,
        ". The first model must be the unrestricted model.",
        call. = FALSE
      )
    }

    list(
      same_model = FALSE,
      null = rule$null,
      boundary = rule$boundary
    )
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

  regression_parameter_names <- function(model) {
    if (isTRUE(model$bivariate)) {
      p <- as.integer(model$numbetas)

      if (length(p) != 2L || any(!is.finite(p)) || any(p < 1L)) {
        return(character(0))
      }

      block_names <- function(prefix, size) {
        if (size == 1L) {
          prefix
        } else {
          c(prefix, paste0(prefix, seq_len(size - 1L)))
        }
      }

      return(c(
        block_names("mean_1", p[1L]),
        block_names("mean_2", p[2L])
      ))
    }

    p <- if (!is.null(model$numbetas)) {
      as.integer(model$numbetas[1L])
    } else if (!is.null(model$X)) {
      ncol(as.matrix(model$X))
    } else {
      1L
    }

    if (length(p) != 1L || !is.finite(p) || p < 1L) {
      return(character(0))
    }

    if (p == 1L) {
      "mean"
    } else {
      c("mean", paste0("mean", seq_len(p - 1L)))
    }
  }

  regression_restrictions <- function(model1, model2, par1, par2, fixed1, fixed2) {
    X1 <- model1$X
    X2 <- model2$X

    if (same_value(X1, X2)) {
      return(list(
        different = FALSE,
        R = matrix(0, nrow = 0L, ncol = length(par1),
                   dimnames = list(NULL, names(par1))),
        r = numeric(0),
        reg1 = regression_parameter_names(model1),
        reg2 = regression_parameter_names(model2)
      ))
    }

    if (isTRUE(model1$bivariate) || isTRUE(model2$bivariate)) {
      stop(
        "Different design matrices are currently supported only for ",
        "univariate fits",
        call. = FALSE
      )
    }

    if (is.list(X1) || is.list(X2) || is.null(X1) || is.null(X2)) {
      stop(
        "The two fits have different design matrices that cannot be compared",
        call. = FALSE
      )
    }

    X1 <- as.matrix(X1)
    X2 <- as.matrix(X2)
    storage.mode(X1) <- storage.mode(X2) <- "double"

    if (nrow(X1) != nrow(X2)) {
      stop(
        "The unrestricted and restricted design matrices have different ",
        "numbers of rows",
        call. = FALSE
      )
    }

    p1 <- ncol(X1)
    p2 <- ncol(X2)

    if (p2 >= p1) {
      stop(
        "The restricted design matrix must contain fewer regression ",
        "parameters than the unrestricted design matrix",
        call. = FALSE
      )
    }

    reg1 <- regression_parameter_names(model1)
    reg2 <- regression_parameter_names(model2)

    if (length(reg1) != p1 || length(reg2) != p2) {
      stop(
        "The number of regression parameter names is inconsistent with X",
        call. = FALSE
      )
    }

    ## With different design matrices, the current implementation assumes
    ## that all regression coefficients in both models are estimated.
    if (!all(reg1 %in% names(par1)) || !all(reg2 %in% names(par2)) ||
        any(reg1 %in% names(fixed1)) || any(reg2 %in% names(fixed2))) {
      stop(
        "Regression nesting with different design matrices requires all ",
        "regression coefficients to be estimated in both models",
        call. = FALSE
      )
    }

    tol <- 1e-8
    rank1 <- qr(X1, tol = tol)$rank
    rank2 <- qr(X2, tol = tol)$rank

    if (rank1 != p1 || rank2 != p2) {
      stop(
        "The regression design matrices must have full column rank",
        call. = FALSE
      )
    }

    ## Find A such that X2 = X1 A. The restricted regression is nested if
    ## the columns of X2 belong to the column space of X1.
    A <- tryCatch(
      qr.solve(X1, X2, tol = tol),
      error = function(e) NULL
    )

    if (is.null(A)) {
      stop(
        "Unable to verify nesting of the regression design matrices",
        call. = FALSE
      )
    }

    A <- as.matrix(A)
    if (nrow(A) != p1) {
      A <- matrix(A, nrow = p1, ncol = p2)
    }

    residual <- X2 - X1 %*% A
    relative_error <- sqrt(sum(residual^2)) /
      max(1, sqrt(sum(X2^2)))

    if (!is.finite(relative_error) || relative_error > tol) {
      stop(
        "The restricted regression is not nested in the unrestricted ",
        "regression: col(X_restricted) is not contained in ",
        "col(X_unrestricted)",
        call. = FALSE
      )
    }

    if (qr(A, tol = tol)$rank != p2) {
      stop(
        "The restricted regression parametrization is rank deficient",
        call. = FALSE
      )
    }

    ## Under H0, beta_unrestricted belongs to col(A). An orthonormal basis
    ## of the complement of col(A) gives C beta_unrestricted = 0.
    QA <- qr.Q(qr(A, tol = tol), complete = TRUE)
    complement <- seq.int(p2 + 1L, p1)
    C <- t(QA[, complement, drop = FALSE])

    R <- matrix(
      0,
      nrow = nrow(C),
      ncol = length(par1),
      dimnames = list(
        paste0("regression_constraint_", seq_len(nrow(C))),
        names(par1)
      )
    )
    R[, reg1] <- C

    list(
      different = TRUE,
      R = R,
      r = setNames(rep(0, nrow(C)), rownames(R)),
      reg1 = reg1,
      reg2 = reg2
    )
  }

  check_common_fit_structure <- function(model1, model2) {
    ## The correlation family must always be identical. Individual
    ## correlation parameters may instead be fixed in the restricted fit.
    fields <- c(
      "corrmodel", "likelihood", "type", "copula",
      "distance", "grid", "n", "neighb", "maxdist", "maxtime",
      "p_neighb", "thin_method", "weighted", "bivariate",
      "spacetime", "radius", "local", "ns", "numpairs"
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
      "data", "coordx", "coordy", "coordz", "coordx_dyn"
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
        "The fitted models must use the same data and coordinates; ",
        "different fields: ", paste(bad_data, collapse = ", "),
        call. = FALSE
      )
    }
  }

  selector_restrictions <- function(parameter_names, null_values, par_names) {
    if (length(parameter_names) == 0L) {
      return(list(
        R = matrix(0, nrow = 0L, ncol = length(par_names),
                   dimnames = list(NULL, par_names)),
        r = numeric(0)
      ))
    }

    R <- matrix(
      0,
      nrow = length(parameter_names),
      ncol = length(par_names),
      dimnames = list(parameter_names, par_names)
    )

    R[cbind(seq_along(parameter_names), match(parameter_names, par_names))] <- 1

    list(
      R = R,
      r = unname(null_values[parameter_names])
    )
  }

  boundary_reasons <- function(parameter_names, null_values, marginal_rule) {
    reasons <- marginal_rule$boundary

    if (length(parameter_names) > 0L) {
      for (nm in parameter_names) {
        value <- as.numeric(null_values[nm])

        if (grepl("^nugget", nm) && isTRUE(all.equal(value, 0))) {
          reasons <- c(reasons, paste0(nm, " = 0"))
        }

        if (identical(nm, "sep") &&
            (isTRUE(all.equal(value, 0)) || isTRUE(all.equal(value, 1)))) {
          reasons <- c(reasons, paste0(nm, " = ", value))
        }

        if (grepl("^pcol", nm) && isTRUE(all.equal(abs(value), 1))) {
          reasons <- c(reasons, paste0(nm, " = ", value))
        }
      }
    }

    unique(reasons[nzchar(reasons)])
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

    if (length(par1) <= length(par2)) {
      stop(
        "The first model must be the unrestricted model and must have more ",
        "estimated parameters than the second model",
        call. = FALSE
      )
    }

    marginal_rule <- marginal_nesting_rule(model1, model2)

    regression <- regression_restrictions(
      model1, model2, par1, par2, fixed1, fixed2
    )

    if (regression$different) {
      nonreg1 <- setdiff(names(par1), regression$reg1)
      nonreg2 <- setdiff(names(par2), regression$reg2)

      if (!all(nonreg2 %in% nonreg1)) {
        stop(
          "Models are not nested: the restricted model estimates non-regression ",
          "parameters not estimated by the unrestricted model: ",
          paste(setdiff(nonreg2, nonreg1), collapse = ", "),
          call. = FALSE
        )
      }

      tested_parameters <- setdiff(nonreg1, nonreg2)
    } else {
      if (!all(names(par2) %in% names(par1))) {
        stop(
          "Models are not nested: the restricted model estimates parameters ",
          "not estimated by the unrestricted model: ",
          paste(setdiff(names(par2), names(par1)), collapse = ", "),
          call. = FALSE
        )
      }

      tested_parameters <- setdiff(names(par1), names(par2))
    }

    if (!marginal_rule$same_model) {
      required <- names(marginal_rule$null)
      missing_required <- setdiff(required, tested_parameters)

      if (length(missing_required) > 0L) {
        stop(
          "The unrestricted marginal model must estimate the parameter(s) ",
          "defining the reduction to the restricted model: ",
          paste(missing_required, collapse = ", "),
          call. = FALSE
        )
      }
    }

    null_values <- fixed2[tested_parameters]

    rule_names <- intersect(tested_parameters, names(marginal_rule$null))
    if (length(rule_names) > 0L) {
      both <- rule_names[rule_names %in% names(fixed2)]

      if (length(both) > 0L &&
          !same_value(fixed2[both], marginal_rule$null[both])) {
        stop(
          "The fixed values in the restricted fit are inconsistent with ",
          "the marginal nesting relation",
          call. = FALSE
        )
      }

      null_values[rule_names] <- marginal_rule$null[rule_names]
    }

    missing_null <- tested_parameters[
      !tested_parameters %in% names(null_values) |
        is.na(null_values[tested_parameters])
    ]

    if (length(missing_null) > 0L) {
      stop(
        "The restricted model must fix every additional non-regression ",
        "parameter of the unrestricted model. Missing null value(s): ",
        paste(missing_null, collapse = ", "),
        call. = FALSE
      )
    }

    null_values <- null_values[tested_parameters]

    if (any(!is.finite(null_values))) {
      stop(
        "The null parameter values must be finite. Limiting models such as ",
        "StudentT -> Gaussian are not supported by the standard tests",
        call. = FALSE
      )
    }

    common_fixed <- intersect(names(fixed1), names(fixed2))

    if (length(common_fixed) > 0L &&
        !same_value(fixed1[common_fixed], fixed2[common_fixed])) {
      stop(
        "Common fixed parameters must have identical values in the two models",
        call. = FALSE
      )
    }

    selector <- selector_restrictions(
      tested_parameters,
      null_values,
      names(par1)
    )

    R <- rbind(regression$R, selector$R)
    r <- c(regression$r, selector$r)

    if (nrow(R) < 1L) {
      stop("Models are not nested or have the same number of parameters", call. = FALSE)
    }

    df_expected <- length(par1) - length(par2)

    if (nrow(R) != df_expected) {
      stop(
        "The number of identified restrictions (", nrow(R),
        ") does not agree with the difference in estimated parameters (",
        df_expected, ")",
        call. = FALSE
      )
    }

    rank_R <- qr(t(R), tol = 1e-8)$rank
    if (rank_R != nrow(R)) {
      stop("The null restrictions are linearly dependent", call. = FALSE)
    }

    reasons <- boundary_reasons(
      tested_parameters,
      null_values,
      marginal_rule
    )

    list(
      par1 = par1,
      par2 = par2,
      R = R,
      r = as.numeric(r),
      df = nrow(R),
      tested = tested_parameters,
      null_values = null_values,
      boundary = reasons
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

    R_null <- nesting$R
    theta <- as.numeric(R_null %*% nesting$par1 - nesting$r)
    df <- nesting$df

    if (statistic == "Wald") {
      V <- ensure_named_matrix(
        model1$varcov,
        names(nesting$par1),
        "The covariance matrix of the unrestricted model"
      )

      V_null <- R_null %*% V %*% t(R_null)
      V_null_inv <- safe_solve(
        V_null,
        "The covariance matrix for the null restrictions"
      )

      W <- as.numeric(crossprod(theta, V_null_inv %*% theta))
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

      ## For general linear restrictions R theta = r, the efficient
      ## sensitivity and Godambe covariance of the restrictions are
      ## H_eff = (R H^{-1} R')^{-1} and V_eff = R V R'.
      H_inv <- safe_solve(H, "The sensitivity matrix")
      H_restriction_inv <- R_null %*% H_inv %*% t(R_null)
      H_eff <- safe_solve(
        H_restriction_inv,
        "The efficient sensitivity matrix for the null restrictions"
      )

      V_eff <- R_null %*% V %*% t(R_null)
      V_eff <- (V_eff + t(V_eff)) / 2

      ## H_eff %*% V_eff has positive real eigenvalues. The symmetric
      ## similar matrix below is numerically more stable.
      chol_H <- tryCatch(chol(H_eff), error = function(e) NULL)

      if (is.null(chol_H)) {
        stop(
          "The efficient sensitivity matrix is not positive definite",
          call. = FALSE
        )
      }

      lambda_mat <- chol_H %*% V_eff %*% t(chol_H)
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
          "the tested restrictions may be weakly identified",
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

    if (length(nesting$boundary) > 0L) {
      warning(
        "Some null values lie on the boundary of the parameter space (",
        paste(nesting$boundary, collapse = "; "),
        "). The usual chi-square reference distribution may not be valid.",
        call. = FALSE
      )
    }

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