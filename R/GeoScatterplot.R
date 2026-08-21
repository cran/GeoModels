####################################################
### File name: GeoScatterplot.R
####################################################

.GeoScatterplot_native_models <- c(
    "Gaussian", "SkewGaussian", "Tukeyh", "StudentT",
    "SinhAsinh", "LogGaussian", "Gamma", "Weibull"
)

## Continuous marginal models supported by the pairwise copula likelihood
## and by the shared bivariate-density engine.
.GeoScatterplot_copula_models <- c(
    "Gaussian", "StudentT", "LogGaussian", "Gamma", "Weibull",
    "Beta", "Beta2", "Kumaraswamy", "Kumaraswamy2",
    "Logistic", "SkewLaplace"
)

## Models for which the marginal CDF is continuous and can therefore be used
## on the Uniform/Gaussian diagnostic scales.
.GeoScatterplot_continuous_models <- c(
    "Gaussian", "SkewGaussian", "StudentT", "SkewStudentT",
    "Gaussian_misp_StudentT", "Gaussian_misp_SkewStudentT",
    "LogGaussian", "LogLogistic", "Gamma", "Weibull",
    "Beta", "Beta2", "Kumaraswamy", "Kumaraswamy2",
    "Logistic", "SkewLaplace", "SinhAsinh", "Tukeyh", "Tukeyh2",
    "TwoPieceGaussian", "TwoPieceStudentT",
    "TwoPieceTukeyh", "TwoPieceBimodal"
)

.GeoScatterplot_normalize_model <- function(model) {
    .GeoMarginalCanonicalModel(model)
}

.GeoScatterplot_three_columns <- function(z, validate = FALSE) {
    z <- as.matrix(z)
    storage.mode(z) <- "double"
    if (ncol(z) == 2L) z <- cbind(z, 0)
    if (isTRUE(validate) && ncol(z) != 3L) {
        stop("coordinates must have two or three columns")
    }
    z
}

.GeoScatterplot_representative_lag <- function(x, method, empty = 0) {
    x <- x[is.finite(x)]
    if (!length(x)) return(empty)
    if (identical(method, "median")) stats::median(x) else mean(x)
}

.GeoScatterplot_match_scale <- function(scale) {
    if (length(scale) > 1L) scale <- scale[1L]
    if (length(scale) != 1L || is.na(scale)) {
        stop("scale must be one of 'Original', 'Uniform', or 'Gaussian'")
    }

    key <- tolower(gsub("[[:blank:]_-]", "", as.character(scale)))
    aliases <- c(
        original = "Original",
        response = "Original",
        raw = "Original",
        uniform = "Uniform",
        pit = "Uniform",
        gaussian = "Gaussian",
        normal = "Gaussian"
    )

    if (!(key %in% names(aliases))) {
        stop("scale must be one of 'Original', 'Uniform', or 'Gaussian'")
    }
    unname(aliases[[key]])
}

.GeoScatterplot_parameters <- function(object) {
    fixed <- object$fixed

    ## An externally supplied mean can be a vector. It is not a model
    ## parameter and must not be flattened into the parameter vector.
    if (!is.null(.GeoMean_external(fixed))) fixed$mean <- NULL

    out <- c(
        unlist(object$param, use.names = TRUE),
        unlist(fixed, use.names = TRUE)
    )

    if (length(out)) {
        out <- out[!duplicated(names(out))]
    }
    out
}

.GeoScatterplot_fitted_location <- function(object, n) {
    external_mean <- .GeoMean_external(object$fixed, n = n)
    if (!is.null(external_mean)) return(external_mean)

    par <- .GeoScatterplot_parameters(object)
    if (is.null(object$X)) {
        X <- matrix(1, nrow = n, ncol = 1L)
    } else {
        X <- .GeoMean_design(object$X, n = n, name = "X")
    }
    beta <- .GeoMean_beta(par, ncol(X))
    as.numeric(X %*% beta)
}

.GeoScatterplot_coordinates <- function(object) {
    if (is.matrix(object$coordx) || is.data.frame(object$coordx)) {
        coords <- as.matrix(object$coordx)
    } else {
        if (is.null(object$coordy)) {
            stop("The GeoFit object does not contain valid spatial coordinates")
        }
        if (is.null(object$coordz)) {
            coords <- cbind(object$coordx, object$coordy)
        } else {
            coords <- cbind(object$coordx, object$coordy, object$coordz)
        }
    }

    storage.mode(coords) <- "double"
    if (ncol(coords) == 2L) coords <- cbind(coords, 0)
    if (ncol(coords) != 3L) {
        stop("The GeoFit object must contain two- or three-dimensional coordinates")
    }
    coords
}

.GeoScatterplot_spacetime_layout <- function(object) {
    y <- .GeoScatterplot_data_vector(object)

    if (!isTRUE(object$spacetime)) {
        coords <- .GeoScatterplot_coordinates(object)
        if (nrow(coords) != length(y)) {
            stop("The number of coordinates differs from the number of observations")
        }
        return(list(
            data = y, coords_obs = coords, coords_base = coords,
            coords_dyn = NULL, time_obs = rep(0, length(y)),
            time_values = 0, ns = nrow(coords), dynamic = FALSE
        ))
    }

    coordt <- as.numeric(object$coordt)
    if (length(coordt) < 2L || any(!is.finite(coordt))) {
        stop("The spatio-temporal GeoFit object does not contain valid temporal coordinates")
    }

    if (is.null(object$coordx_dyn)) {
        coords_base <- .GeoScatterplot_coordinates(object)
        nsite <- nrow(coords_base)
        nt <- length(coordt)
        expected <- nsite * nt
        if (length(y) != expected) {
            stop(
                "For fixed spatial locations, the number of observations must equal ",
                "the number of sites times the number of temporal instants"
            )
        }

        coords_obs <- coords_base[
            rep(seq_len(nsite), times = nt), , drop = FALSE
        ]
        time_obs <- rep(coordt, each = nsite)

        return(list(
            data = y, coords_obs = coords_obs, coords_base = coords_base,
            coords_dyn = NULL, time_obs = time_obs,
            time_values = coordt, ns = rep.int(nsite, nt), dynamic = FALSE
        ))
    }

    if (!is.list(object$coordx_dyn) ||
        length(object$coordx_dyn) != length(coordt)) {
        stop("coordx_dyn must be a list with one coordinate matrix per temporal instant")
    }

    coords_dyn <- lapply(object$coordx_dyn, function(z) {
        z <- as.matrix(z)
        storage.mode(z) <- "double"
        if (ncol(z) == 2L) z <- cbind(z, 0)
        if (ncol(z) != 3L) {
            stop("Each element of coordx_dyn must have two or three columns")
        }
        z
    })
    ns <- vapply(coords_dyn, nrow, integer(1))
    coords_obs <- do.call(rbind, coords_dyn)
    time_obs <- rep(coordt, times = ns)

    if (length(y) != nrow(coords_obs)) {
        stop("The dynamic coordinates and the observations have different lengths")
    }

    list(
        data = y, coords_obs = coords_obs, coords_base = NULL,
        coords_dyn = coords_dyn, time_obs = time_obs,
        time_values = coordt, ns = ns, dynamic = TRUE
    )
}

.GeoScatterplot_unique_tol <- function(x, tol = NULL) {
    x <- sort(as.numeric(x[is.finite(x)]))
    if (!length(x)) return(numeric(0))
    if (is.null(tol)) {
        tol <- sqrt(.Machine$double.eps) * max(1, abs(x))
    }
    out <- x[1L]
    if (length(x) == 1L) return(out)
    for (v in x[-1L]) {
        local_tol <- max(tol, sqrt(.Machine$double.eps) * max(1, abs(v), abs(out[length(out)])))
        if (abs(v - out[length(out)]) > local_tol) out <- c(out, v)
    }
    out
}

.GeoScatterplot_distinct_time_lags <- function(coordt, maxtime = NULL) {
    tt <- sort(as.numeric(coordt))
    if (!length(tt) || any(!is.finite(tt))) {
        stop("Temporal coordinates must be finite")
    }
    if (length(tt) == 1L) return(0)

    tol <- sqrt(.Machine$double.eps) * max(1, abs(tt))
    dd <- diff(tt)
    regular <- length(dd) <= 1L ||
        max(abs(dd - stats::median(dd))) <= max(tol, sqrt(.Machine$double.eps) * max(1, abs(dd)))

    if (regular) {
        ## Every lag of a regular grid is represented by a difference from
        ## the first time point; this avoids constructing a T by T matrix.
        possible <- c(0, tt[-1L] - tt[1L])
    } else {
        ## Irregular grids can genuinely have O(T^2) distinct lags.  Build
        ## them incrementally so memory is proportional to the distinct lag
        ## set rather than to a dense outer-difference matrix.
        possible <- 0
        nt <- length(tt)
        for (i in seq_len(nt - 1L)) {
            d <- tt[(i + 1L):nt] - tt[i]
            if (!is.null(maxtime) && length(maxtime) == 1L && is.finite(maxtime)) {
                d <- d[d <= maxtime + tol]
            }
            if (length(d)) {
                possible <- .GeoScatterplot_unique_tol(c(possible, d), tol = tol)
            }
        }
    }

    possible <- .GeoScatterplot_unique_tol(possible, tol = tol)
    if (!is.null(maxtime) && length(maxtime) == 1L && is.finite(maxtime) && maxtime >= 0) {
        possible <- possible[possible <= maxtime + tol]
    }
    possible
}

.GeoScatterplot_time_lags <- function(object, time.lag = NULL) {
    if (!isTRUE(object$spacetime)) return(0)

    coordt <- as.numeric(object$coordt)
    maxtime <- if (is.null(time.lag)) object$maxtime else NULL
    possible <- .GeoScatterplot_distinct_time_lags(coordt, maxtime = maxtime)
    tol <- sqrt(.Machine$double.eps) * max(1, abs(coordt), abs(possible))

    if (is.null(time.lag)) return(possible)

    if (!is.numeric(time.lag) || !length(time.lag) ||
        any(!is.finite(time.lag)) || any(time.lag < 0)) {
        stop("time.lag must contain non-negative finite temporal lags")
    }

    matched <- vapply(time.lag, function(u) {
        k <- which.min(abs(possible - u))
        local_tol <- max(tol, sqrt(.Machine$double.eps) * max(1, abs(u)))
        if (!length(k) || abs(possible[k] - u) > local_tol) {
            stop(
                "Temporal lag ", signif(u, 6),
                " is not available from the temporal coordinates in the GeoFit object"
            )
        }
        possible[k]
    }, numeric(1))

    .GeoScatterplot_unique_tol(matched, tol = tol)
}

.GeoScatterplot_neighbor_coordinates <- function(coords, distance = "Eucl") {
    distance_key <- tolower(gsub("[[:blank:]]", "", as.character(distance)[1L]))

    trim_one <- function(z) {
        z <- as.matrix(z)
        storage.mode(z) <- "double"
        if (distance_key != "eucl" && ncol(z) == 3L) {
            z <- z[, 1:2, drop = FALSE]
        }
        z
    }

    if (is.list(coords)) lapply(coords, trim_one) else trim_one(coords)
}

.GeoScatterplot_spatial_pairs_maxdist <- function(data, coords, maxdist,
                                                  bins, distance = "Eucl",
                                                  radius = 1) {
    coords <- as.matrix(coords)
    storage.mode(coords) <- "double"
    data <- as.numeric(data)
    n <- nrow(coords)

    if (length(data) != n) {
        stop("The data and coordinates must have the same length")
    }
    if (length(maxdist) != 1L || !is.finite(maxdist) || maxdist <= 0) {
        stop("maxdist must be a positive finite scalar")
    }
    if (!is.numeric(bins) || length(bins) < 2L ||
        any(!is.finite(bins)) || any(diff(bins) <= 0)) {
        stop("bins must be a strictly increasing finite numeric vector")
    }
    if (n < 2L) {
        return(list(
            bin = numeric(0), v1 = numeric(0), v2 = numeric(0),
            i = integer(0), j = integer(0), lags = numeric(0)
        ))
    }

    distance_key <- tolower(gsub("[[:blank:]]", "", as.character(distance)[1L]))
    if (!(distance_key %in% c("eucl", "geod", "chor"))) {
        stop("Unsupported distance in GeoScatterplot")
    }
    if (!is.numeric(radius) || length(radius) != 1L ||
        !is.finite(radius) || radius <= 0) {
        stop("radius must be a positive finite scalar")
    }

    max_cells <- getOption("GeoModels.scatterplot.maxdist_cells", 2000000L)
    if (!is.numeric(max_cells) || length(max_cells) != 1L ||
        !is.finite(max_cells) || max_cells < n) {
        max_cells <- max(n, 2000000L)
    }
    block_size <- max(1L, min(n, as.integer(floor(max_cells / n))))
    step <- bins[2L] - bins[1L]
    nbins <- length(bins) - 1L

    if (distance_key == "eucl") {
        coords_norm2 <- rowSums(coords * coords)
    } else {
        lon <- coords[, 1L] * pi / 180
        lat <- coords[, 2L] * pi / 180
        sin_lat <- sin(lat)
        cos_lat <- cos(lat)
    }

    out_bin <- out_v1 <- out_v2 <- out_lag <- list()
    out_i <- out_j <- list()
    out_k <- 0L

    starts <- seq.int(1L, n, by = block_size)
    for (start in starts) {
        rows <- start:min(n, start + block_size - 1L)

        if (distance_key == "eucl") {
            q <- coords[rows, , drop = FALSE]
            d2 <- outer(rowSums(q * q), coords_norm2, "+") -
                2 * tcrossprod(q, coords)
            d2[d2 < 0] <- 0
            dmat <- sqrt(d2)
        } else {
            cs <- outer(sin_lat[rows], sin_lat, "*") +
                outer(cos_lat[rows], cos_lat, "*") *
                cos(outer(lon[rows], lon, "-"))
            cs[cs > 1] <- 1
            cs[cs < -1] <- -1
            theta <- acos(cs)
            dmat <- if (distance_key == "geod") {
                radius * theta
            } else {
                2 * radius * sin(theta / 2)
            }
        }

        hit <- which(is.finite(dmat) & dmat <= maxdist, arr.ind = TRUE)
        if (!nrow(hit)) next
        i <- rows[hit[, 1L]]
        j <- hit[, 2L]
        keep <- i < j & is.finite(data[i]) & is.finite(data[j])
        if (!any(keep)) next

        i <- i[keep]
        j <- j[keep]
        lagv <- dmat[hit[keep, , drop = FALSE]]

        bin_idx <- floor((lagv - bins[1L]) / step) + 1L
        at_upper <- lagv == bins[length(bins)]
        bin_idx[at_upper] <- nbins
        safe_idx <- pmax(1L, pmin(nbins, bin_idx))
        valid_bin <- bin_idx >= 1L & bin_idx <= nbins &
            lagv >= bins[safe_idx] & lagv <= bins[safe_idx + 1L]
        if (!any(valid_bin)) next

        i <- i[valid_bin]
        j <- j[valid_bin]
        lagv <- lagv[valid_bin]
        bin_idx <- bin_idx[valid_bin]

        out_k <- out_k + 1L
        out_bin[[out_k]] <- bins[bin_idx]
        out_v1[[out_k]] <- data[i]
        out_v2[[out_k]] <- data[j]
        out_i[[out_k]] <- i
        out_j[[out_k]] <- j
        out_lag[[out_k]] <- lagv
    }

    if (!out_k) {
        return(list(
            bin = numeric(0), v1 = numeric(0), v2 = numeric(0),
            i = integer(0), j = integer(0), lags = numeric(0)
        ))
    }

    list(
        bin = unlist(out_bin, use.names = FALSE),
        v1 = unlist(out_v1, use.names = FALSE),
        v2 = unlist(out_v2, use.names = FALSE),
        i = as.integer(unlist(out_i, use.names = FALSE)),
        j = as.integer(unlist(out_j, use.names = FALSE)),
        lags = unlist(out_lag, use.names = FALSE)
    )
}

.GeoScatterplot_distance_label <- function(bins, k) {
    left <- "["
    right <- if (k == length(bins) - 1L) "]" else ")"
    paste0(left, signif(bins[k], 3), ", ", signif(bins[k + 1L], 3), right)
}

.GeoScatterplot_bin_membership <- function(lags, bins, k) {
    if (k < length(bins) - 1L) {
        lags >= bins[k] & lags < bins[k + 1L]
    } else {
        lags >= bins[k] & lags <= bins[k + 1L]
    }
}

.GeoScatterplot_nuisance <- function(object, model, par) {
    nb <- object$numbetas
    if (is.null(nb)) nb <- 1L

    ## The marginal nuisance vector is kept separate from the copula
    ## parameter. This avoids relying on model-specific copula indexing.
    nuisance_names <- NuisParam2(
        model = model,
        bivariate = FALSE,
        num_betas = nb,
        copula = NULL
    )
    nuisance_names <- nuisance_names[substr(nuisance_names, 1L, 4L) != "mean"]

    out <- numeric(length(nuisance_names))
    names(out) <- nuisance_names

    for (k in seq_along(nuisance_names)) {
        nm <- nuisance_names[k]
        if (nm %in% names(par)) {
            out[k] <- as.numeric(par[[nm]])
        } else if (nm == "sill" && model %in% c("Gamma", "Weibull", "Beta", "Beta2", "Kumaraswamy", "Kumaraswamy2")) {
            ## GeoFit removes structurally fixed sill values from its output.
            out[k] <- 1
        } else {
            stop("Parameter '", nm, "' is missing from the GeoFit object")
        }
    }

    out
}

.GeoScatterplot_copula_parameter <- function(object, model, par) {
    if (is.null(object$copula) || object$copula == "Gaussian") return(0)

    nb <- object$numbetas
    if (is.null(nb)) nb <- 1L

    base_names <- NuisParam2(
        model = model, bivariate = FALSE,
        num_betas = nb, copula = NULL
    )
    full_names <- NuisParam2(
        model = model, bivariate = FALSE,
        num_betas = nb, copula = object$copula
    )

    extra <- setdiff(full_names, base_names)
    extra <- extra[extra %in% names(par)]

    if (!length(extra)) {
        candidates <- if (object$copula == "Clayton") {
            c("nu", "copula_nu")
        } else {
            c("skew_cop", "copula_skew", "skewcop", "skew")
        }
        extra <- candidates[candidates %in% names(par)]
    }

    if (length(extra) != 1L) {
        stop(
            "Unable to identify the parameter of the ", object$copula,
            " copula in the GeoFit object"
        )
    }

    as.numeric(par[[extra]])
}

.GeoScatterplot_density_setup <- function(object, model) {
    par <- .GeoScatterplot_parameters(object)
    nuisance <- .GeoScatterplot_nuisance(object, model, par)

    ## A Marginal + Independence fit estimates only the univariate marginal
    ## model. Correlation and copula parameters, even when supplied in fixed,
    ## do not belong to the fitted independence likelihood and must therefore
    ## not enter the theoretical bivariate contour.
    independence <- identical(object$likelihood, "Marginal") &&
        identical(object$type, "Independence")

    if (independence) {
        return(list(
            model = as.integer(CkModel(model)),
            model_name = model,
            independence = TRUE,
            corrmodel_name = NULL,
            cormod = NA_integer_,
            type_cop = 0L,
            copula_name = "Independence",
            copula_par = 0,
            corr_names = character(0),
            corrpar = numeric(0),
            corrpar_named = numeric(0),
            nuisance = as.double(nuisance),
            nuisance_named = nuisance
        ))
    }

    corr_names <- CorrParam(object$corrmodel)

    if (!all(corr_names %in% names(par))) {
        missing <- corr_names[!(corr_names %in% names(par))]
        stop("Missing correlation parameter(s): ", paste(missing, collapse = ", "))
    }

    type_cop <- if (is.null(object$copula)) {
        0L
    } else {
        switch(
            object$copula,
            Gaussian = 1L,
            Clayton = 2L,
            SkewGaussian = 3L,
            stop("Unsupported copula in the GeoFit object")
        )
    }

    list(
        model = as.integer(CkModel(model)),
        model_name = model,
        independence = FALSE,
        corrmodel_name = object$corrmodel,
        cormod = as.integer(CkCorrModel(object$corrmodel)),
        type_cop = as.integer(type_cop),
        copula_name = object$copula,
        copula_par = .GeoScatterplot_copula_parameter(object, model, par),
        corr_names = corr_names,
        corrpar = as.double(par[corr_names]),
        corrpar_named = par[corr_names],
        nuisance = as.double(nuisance),
        nuisance_named = nuisance
    )
}

.GeoScatterplot_residual_n <- function(n, nobs, model) {
    if (is.null(n))
        stop("Residual-scale scatterplots for model '", model, "' require n.")
    z <- if (is.list(n)) unlist(n, use.names = FALSE) else as.numeric(n)
    if (length(z) == 1L) z <- rep(z, nobs)
    if (length(z) != nobs || any(!is.finite(z)) || any(z <= 0) ||
        any(abs(z - round(z)) > sqrt(.Machine$double.eps))) {
        stop("Residual-scale scatterplots for model '", model,
             "' require positive integer n with length 1 or the number of observations.")
    }
    as.numeric(z)
}

.GeoScatterplot_residual_scale <- function(data, mean, model, object = NULL) {
    if (model %in% c("LogGaussian", "Gamma", "Weibull")) {
        sc <- exp(mean)
        if (any(!is.finite(sc)) || any(sc <= 0)) {
            stop("Invalid fitted multiplicative mean")
        }
        return(data / sc)
    }

    if (model %in% .GeoScatterplot_native_models) {
        return(data - mean)
    }

    ## Count-model residual plots use the same deterministic Pearson scale as
    ## GeoResiduals().  This provides a meaningful common scale for pair plots
    ## without pretending that a continuous bivariate density is available.
    pearson <- function(obs, fitted_mean, fitted_var) {
        if (any(!is.finite(fitted_mean)) || any(!is.finite(fitted_var)) ||
            any(fitted_var < 0))
            stop("Invalid fitted count moments for residual-scale scatterplot")
        (obs - fitted_mean) / sqrt(pmax(fitted_var, .Machine$double.eps))
    }

    if (model == "Poisson") {
        fitted_mean <- exp(mean)
        return(pearson(data, fitted_mean, fitted_mean))
    }

    if (model %in% c("Binary", "Bernoulli")) {
        prob <- pnorm(mean)
        return(pearson(data, prob, prob * (1 - prob)))
    }

    if (model == "Binomial") {
        if (is.null(object)) stop("A GeoFit object is required for Binomial residuals")
        nn <- .GeoScatterplot_residual_n(object$n, length(data), model)
        prob <- pnorm(mean)
        fitted_mean <- nn * prob
        return(pearson(data, fitted_mean, nn * prob * (1 - prob)))
    }

    if (model %in% c("Geom", "Geometric", "BinomialNeg")) {
        if (is.null(object)) stop("A GeoFit object is required for Negative-Binomial residuals")
        nn <- if (model %in% c("Geom", "Geometric")) {
            rep(1, length(data))
        } else {
            .GeoScatterplot_residual_n(object$n, length(data), model)
        }
        prob <- pnorm(mean)
        fitted_mean <- nn * (1 - prob) / prob
        fitted_var <- nn * (1 - prob) / prob^2
        return(pearson(data, fitted_mean, fitted_var))
    }

    if (model == "PoissonGamma") {
        if (is.null(object)) stop("A GeoFit object is required for PoissonGamma residuals")
        par <- .GeoScatterplot_parameters(object)
        shape <- suppressWarnings(as.numeric(par["shape"]))
        if (length(shape) != 1L || !is.finite(shape) || shape <= 0)
            stop("PoissonGamma residual-scale scatterplots require a positive fitted shape parameter")
        fitted_mean <- exp(mean)
        fitted_var <- fitted_mean * (1 + fitted_mean / shape)
        return(pearson(data, fitted_mean, fitted_var))
    }

    stop("Residual-scale scatterplots are not implemented for model '", model, "'")
}

.GeoScatterplot_data_vector <- function(object) {
    if (is.list(object$data)) {
        out <- unlist(object$data, use.names = FALSE)
    } else {
        out <- c(t(object$data))
    }
    as.numeric(out)
}

.GeoScatterplot_pit_data <- function(object, scale) {
    if (!(scale %in% c("Uniform", "Gaussian"))) {
        stop("The PIT scale must be either 'Uniform' or 'Gaussian'")
    }

    pit_object <- GeoPit(object, type = scale)
    .GeoScatterplot_data_vector(pit_object)
}

.GeoScatterplot_nuisance_value <- function(nuisance, name,
                                                       default = NA_real_) {
    if (name %in% names(nuisance)) {
        value <- suppressWarnings(as.numeric(nuisance[[name]]))
        if (length(value) == 1L && is.finite(value)) return(value)
    }
    default
}

.GeoScatterplot_marginal_quantile <- function(u, mean, model, nuisance) {
    eps <- sqrt(.Machine$double.eps)
    u <- pmin(1 - eps, pmax(eps, u))
    sill <- .GeoScatterplot_nuisance_value(nuisance, "sill", 1)
    if (!is.finite(sill) || sill <= 0) {
        stop("The fitted marginal variance must be positive")
    }
    sd <- sqrt(sill)
    out <- NULL

    if (model == "Gaussian") {
        out <- mean + sd * stats::qnorm(u)
    }

    if (model == "SkewGaussian") {
        skew <- .GeoScatterplot_nuisance_value(nuisance, "skew")
        if (!is.finite(skew)) stop("The SkewGaussian skew parameter is missing")
        omega <- sqrt((skew^2 + sill) / sill)
        alpha <- skew / sd
        out <- mean + sd * sn::qsn(u, xi = 0, omega = omega, alpha = alpha)
    }

    if (model == "Tukeyh") {
        tail <- .GeoScatterplot_nuisance_value(nuisance, "tail")
        if (!is.finite(tail) || tail < 0) {
            stop("The Tukey-h tail parameter must be non-negative")
        }
        out <- mean + sd * .GeoMarginalQTukeyH(u, tail)
    }

    if (model == "StudentT") {
        inv_df <- .GeoScatterplot_nuisance_value(nuisance, "df")
        if (!is.finite(inv_df) || inv_df <= 0) {
            stop("The StudentT inverse degrees-of-freedom parameter must be positive")
        }
        out <- mean + sd * stats::qt(u, df = 1 / inv_df)
    }

    if (model == "SinhAsinh") {
        skew <- .GeoScatterplot_nuisance_value(nuisance, "skew")
        tail <- .GeoScatterplot_nuisance_value(nuisance, "tail")
        if (!is.finite(skew) || !is.finite(tail) || tail <= 0) {
            stop("The SinhAsinh skew and tail parameters are invalid")
        }
        z <- stats::qnorm(u)
        out <- mean + sd * sinh((asinh(z) + skew) / tail)
    }

    if (model == "LogGaussian") {
        out <- stats::qlnorm(u, meanlog = mean - sill / 2, sdlog = sd)
    }

    if (model == "Gamma") {
        shape <- .GeoScatterplot_nuisance_value(nuisance, "shape")
        if (!is.finite(shape) || shape <= 0) {
            stop("The Gamma shape parameter must be positive")
        }
        out <- stats::qgamma(
            u, shape = shape / 2,
            rate = shape / (2 * exp(mean))
        )
    }

    if (model == "Weibull") {
        shape <- .GeoScatterplot_nuisance_value(nuisance, "shape")
        if (!is.finite(shape) || shape <= 0) {
            stop("The Weibull shape parameter must be positive")
        }
        out <- stats::qweibull(
            u, shape = shape,
            scale = exp(mean) / gamma(1 + 1 / shape)
        )
    }

    if (model == "Beta") {
        shape1 <- .GeoScatterplot_nuisance_value(nuisance, "shape1")
        shape2 <- .GeoScatterplot_nuisance_value(nuisance, "shape2")
        lower <- .GeoScatterplot_nuisance_value(nuisance, "min")
        upper <- .GeoScatterplot_nuisance_value(nuisance, "max")
        if (!is.finite(shape1) || !is.finite(shape2) ||
            shape1 <= 0 || shape2 <= 0 || !is.finite(lower) ||
            !is.finite(upper) || upper <= lower) {
            stop("The fitted Beta parameters are invalid")
        }
        out <- lower + (upper - lower) * stats::qbeta(u, shape1, shape2)
    }

    if (model == "Beta2") {
        shape <- .GeoScatterplot_nuisance_value(nuisance, "shape")
        lower <- .GeoScatterplot_nuisance_value(nuisance, "min")
        upper <- .GeoScatterplot_nuisance_value(nuisance, "max")
        mu <- stats::plogis(mean)
        if (!is.finite(shape) || shape <= 0 || !is.finite(lower) ||
            !is.finite(upper) || upper <= lower) {
            stop("The fitted Beta2 parameters are invalid")
        }
        out <- lower + (upper - lower) * stats::qbeta(
            u, shape1 = mu * shape, shape2 = (1 - mu) * shape
        )
    }

    if (model == "Kumaraswamy") {
        shape1 <- .GeoScatterplot_nuisance_value(nuisance, "shape1")
        shape2 <- .GeoScatterplot_nuisance_value(nuisance, "shape2")
        lower <- .GeoScatterplot_nuisance_value(nuisance, "min")
        upper <- .GeoScatterplot_nuisance_value(nuisance, "max")
        if (!is.finite(shape1) || !is.finite(shape2) ||
            shape1 <= 0 || shape2 <= 0 || !is.finite(lower) ||
            !is.finite(upper) || upper <= lower) {
            stop("The fitted Kumaraswamy parameters are invalid")
        }
        qq <- (1 - (1 - u)^(1 / shape2))^(1 / shape1)
        out <- lower + (upper - lower) * qq
    }

    if (model == "Kumaraswamy2") {
        shape <- .GeoScatterplot_nuisance_value(nuisance, "shape")
        lower <- .GeoScatterplot_nuisance_value(nuisance, "min")
        upper <- .GeoScatterplot_nuisance_value(nuisance, "max")
        mu <- stats::plogis(mean)
        if (!is.finite(shape) || shape <= 0 || !is.finite(lower) ||
            !is.finite(upper) || upper <= lower || mu <= 0 || mu >= 1) {
            stop("The fitted Kumaraswamy2 parameters are invalid")
        }
        shape2 <- log(0.5) / log1p(-mu^shape)
        qq <- (1 - (1 - u)^(1 / shape2))^(1 / shape)
        out <- lower + (upper - lower) * qq
    }

    if (model == "Logistic") {
        out <- mean + sd * stats::qlogis(u)
    }

    if (model == "SkewLaplace") {
        skew <- .GeoScatterplot_nuisance_value(nuisance, "skew")
        if (!is.finite(skew) || skew <= 0 || skew >= 1) {
            stop("The SkewLaplace skew parameter must lie in (0,1)")
        }
        out <- mean + sd * .GeoMarginalQSkewLaplace(u, skew)
    }

    if (is.null(out)) {
        stop("Marginal quantiles are not implemented for model '", model, "'")
    }
    out
}

.GeoScatterplot_marginal_density <- function(x, mean, model, nuisance) {
    sill <- .GeoScatterplot_nuisance_value(nuisance, "sill", 1)
    if (!is.finite(sill) || sill <= 0) {
        stop("The fitted marginal variance must be positive")
    }
    sd <- sqrt(sill)
    out <- NULL

    if (model == "Gaussian") {
        out <- stats::dnorm(x, mean = mean, sd = sd)
    }

    if (model == "SkewGaussian") {
        skew <- .GeoScatterplot_nuisance_value(nuisance, "skew")
        if (!is.finite(skew)) stop("The SkewGaussian skew parameter is missing")
        omega <- sqrt((skew^2 + sill) / sill)
        alpha <- skew / sd
        out <- sn::dsn((x - mean) / sd, xi = 0, omega = omega, alpha = alpha) / sd
    }

    if (model == "Tukeyh") {
        tail <- .GeoScatterplot_nuisance_value(nuisance, "tail")
        if (!is.finite(tail) || tail < 0) {
            stop("The Tukey-h tail parameter must be non-negative and finite")
        }
        out <- .GeoMarginalDTukeyH((x - mean) / sd, tail) / sd
    }

    if (model == "StudentT") {
        inv_df <- .GeoScatterplot_nuisance_value(nuisance, "df")
        if (!is.finite(inv_df) || inv_df <= 0) {
            stop("The StudentT inverse degrees-of-freedom parameter must be positive")
        }
        out <- stats::dt((x - mean) / sd, df = 1 / inv_df) / sd
    }

    if (model == "SinhAsinh") {
        skew <- .GeoScatterplot_nuisance_value(nuisance, "skew")
        tail <- .GeoScatterplot_nuisance_value(nuisance, "tail")
        if (!is.finite(skew) || !is.finite(tail) || tail <= 0) {
            stop("The SinhAsinh skew and tail parameters are invalid")
        }
        xx <- (x - mean) / sd
        bb <- tail * asinh(xx) - skew
        z <- sinh(bb)
        jac <- tail * cosh(bb) / sqrt(1 + xx^2)
        out <- stats::dnorm(z) * jac / sd
    }

    if (model == "LogGaussian") {
        out <- stats::dlnorm(x, meanlog = mean - sill / 2, sdlog = sd)
    }

    if (model == "Gamma") {
        shape <- .GeoScatterplot_nuisance_value(nuisance, "shape")
        if (!is.finite(shape) || shape <= 0) {
            stop("The Gamma shape parameter must be positive")
        }
        out <- stats::dgamma(
            x, shape = shape / 2,
            rate = shape / (2 * exp(mean))
        )
    }

    if (model == "Weibull") {
        shape <- .GeoScatterplot_nuisance_value(nuisance, "shape")
        if (!is.finite(shape) || shape <= 0) {
            stop("The Weibull shape parameter must be positive")
        }
        out <- stats::dweibull(
            x, shape = shape,
            scale = exp(mean) / gamma(1 + 1 / shape)
        )
    }

    if (model == "Beta") {
        shape1 <- .GeoScatterplot_nuisance_value(nuisance, "shape1")
        shape2 <- .GeoScatterplot_nuisance_value(nuisance, "shape2")
        lower <- .GeoScatterplot_nuisance_value(nuisance, "min")
        upper <- .GeoScatterplot_nuisance_value(nuisance, "max")
        width <- upper - lower
        if (!is.finite(shape1) || !is.finite(shape2) ||
            shape1 <= 0 || shape2 <= 0 || !is.finite(width) || width <= 0) {
            stop("The fitted Beta parameters are invalid")
        }
        out <- stats::dbeta((x - lower) / width, shape1, shape2) / width
    }

    if (model == "Beta2") {
        shape <- .GeoScatterplot_nuisance_value(nuisance, "shape")
        lower <- .GeoScatterplot_nuisance_value(nuisance, "min")
        upper <- .GeoScatterplot_nuisance_value(nuisance, "max")
        width <- upper - lower
        mu <- stats::plogis(mean)
        if (!is.finite(shape) || shape <= 0 || !is.finite(width) || width <= 0) {
            stop("The fitted Beta2 parameters are invalid")
        }
        out <- stats::dbeta(
            (x - lower) / width,
            shape1 = mu * shape,
            shape2 = (1 - mu) * shape
        ) / width
    }

    if (model == "Kumaraswamy") {
        shape1 <- .GeoScatterplot_nuisance_value(nuisance, "shape1")
        shape2 <- .GeoScatterplot_nuisance_value(nuisance, "shape2")
        lower <- .GeoScatterplot_nuisance_value(nuisance, "min")
        upper <- .GeoScatterplot_nuisance_value(nuisance, "max")
        width <- upper - lower
        if (!is.finite(shape1) || !is.finite(shape2) ||
            shape1 <= 0 || shape2 <= 0 || !is.finite(width) || width <= 0) {
            stop("The fitted Kumaraswamy parameters are invalid")
        }
        z <- (x - lower) / width
        out <- ifelse(
            z > 0 & z < 1,
            shape1 * shape2 * z^(shape1 - 1) *
                (1 - z^shape1)^(shape2 - 1) / width,
            0
        )
    }

    if (model == "Kumaraswamy2") {
        shape <- .GeoScatterplot_nuisance_value(nuisance, "shape")
        lower <- .GeoScatterplot_nuisance_value(nuisance, "min")
        upper <- .GeoScatterplot_nuisance_value(nuisance, "max")
        width <- upper - lower
        mu <- stats::plogis(mean)
        if (!is.finite(shape) || shape <= 0 || !is.finite(width) ||
            width <= 0 || mu <= 0 || mu >= 1) {
            stop("The fitted Kumaraswamy2 parameters are invalid")
        }
        shape2 <- log(0.5) / log1p(-mu^shape)
        z <- (x - lower) / width
        out <- ifelse(
            z > 0 & z < 1,
            shape * shape2 * z^(shape - 1) *
                (1 - z^shape)^(shape2 - 1) / width,
            0
        )
    }

    if (model == "Logistic") {
        out <- stats::dlogis(x, location = mean, scale = sd)
    }

    if (model == "SkewLaplace") {
        skew <- .GeoScatterplot_nuisance_value(nuisance, "skew")
        if (!is.finite(skew) || skew <= 0 || skew >= 1) {
            stop("The SkewLaplace skew parameter must lie in (0,1)")
        }
        out <- .GeoMarginalDSkewLaplace((x - mean) / sd, skew) / sd
    }

    if (is.null(out)) {
        stop("Marginal densities are not implemented for model '", model, "'")
    }
    out
}

.GeoScatterplot_contour_levels <- function(z, nlevels,
                                                  scale = "Original") {
    zz <- z[is.finite(z) & z > 0]
    if (length(zz) < 2L || diff(range(zz)) <= 0) return(numeric(0))

    if (scale == "Gaussian") {
        ## Levels analogous to those used in the copula comparison figures:
        ## 0.02, 0.06, 0.10, ... . If many are available, retain an evenly
        ## spaced subset controlled by nlevels.
        upper <- max(zz)
        lev <- seq(0.02, upper, by = 0.04)
        lev <- lev[lev > min(zz) & lev < upper]
        if (length(lev) > nlevels) {
            idx <- unique(round(seq(1, length(lev), length.out = nlevels)))
            lev <- lev[idx]
        }
        if (length(lev)) return(lev)
    }

    ## Copula densities can increase sharply near the corners of the unit
    ## square. Quantile-based levels prevent a few extreme grid values from
    ## determining all automatically selected contours.
    if (scale == "Uniform") {
        probs <- seq(0.20, 0.90, length.out = nlevels)
        lev <- unique(as.numeric(stats::quantile(
            zz, probs = probs, names = FALSE, type = 8
        )))
        lev <- lev[lev > min(zz) & lev < max(zz)]
        if (length(lev)) return(lev)
    }

    lev <- pretty(range(zz), n = nlevels)
    lev <- lev[lev > min(zz) & lev < max(zz)]

    if (!length(lev)) {
        lev <- seq(min(zz), max(zz), length.out = nlevels + 2L)
        lev <- lev[-c(1L, length(lev))]
    }
    lev
}

.GeoScatterplot_response_density_grid <- function(setup, xgrid, ygrid,
                                                  lag, time_lag = 0,
                                                  mean1, mean2) {
    gx <- rep(xgrid, times = length(ygrid))
    gy <- rep(ygrid, each = length(xgrid))

    ans <- .C(
        "GeoBivDensity",
        model = as.integer(setup$model),
        cormod = as.integer(setup$cormod),
        type_cop = as.integer(setup$type_cop),
        copula_par = as.double(setup$copula_par),
        x = as.double(gx),
        y = as.double(gy),
        n = as.integer(length(gx)),
        lag = as.double(lag),
        lagt = as.double(time_lag),
        par = as.double(setup$corrpar),
        mean1 = as.double(mean1),
        mean2 = as.double(mean2),
        nuis = as.double(setup$nuisance),
        NN1 = as.integer(1L),
        NN2 = as.integer(1L),
        res = double(length(gx)),
        PACKAGE = "GeoModels",
        NAOK = TRUE
    )$res

    z <- matrix(ans, nrow = length(xgrid), ncol = length(ygrid))
    rho_info <- .GeoScatterplot_copula_rho(setup, lag, time_lag)
    attr(z, "rho") <- rho_info$rho
    attr(z, "rho_latent") <- rho_info$rho_latent
    attr(z, "nugget") <- rho_info$nugget
    z
}

.GeoScatterplot_copula_rho <- function(setup, lag, time_lag = 0) {
    if (length(lag) != 1L || !is.finite(lag) || lag < 0) {
        stop("The representative spatial lag must be a non-negative finite number")
    }
    if (length(time_lag) != 1L || !is.finite(time_lag) || time_lag < 0) {
        stop("The representative temporal lag must be a non-negative finite number")
    }

    corr_param <- as.list(setNames(
        as.numeric(setup$corrpar_named), names(setup$corrpar_named)
    ))
    param <- c(list(mean = 0, nugget = 0, sill = 1), corr_param)

    rho_latent <- as.numeric(GeoCorrFct(
        x = lag,
        t = time_lag,
        corrmodel = setup$corrmodel_name,
        model = "Gaussian",
        param = param,
        covariance = FALSE,
        variogram = FALSE
    )$corr)[1L]

    nugget <- .GeoScatterplot_nuisance_value(
        setup$nuisance_named, "nugget", 0
    )
    rho <- (1 - nugget) * rho_latent

    if (!is.finite(rho_latent) || !is.finite(rho)) {
        stop("Unable to evaluate the fitted correlation at the representative lag")
    }

    eps <- 1e-10
    rho <- max(-1 + eps, min(1 - eps, rho))
    rho_latent <- max(-1 + eps, min(1 - eps, rho_latent))

    list(rho = rho, rho_latent = rho_latent, nugget = nugget)
}

.GeoScatterplot_direct_copula_grid <- function(setup, xgrid, ygrid,
                                                     lag, time_lag = 0, scale) {
    if (setup$type_cop <= 0L) {
        stop("A fitted copula is required for direct copula-density evaluation")
    }

    if (scale == "Uniform") {
        eps <- sqrt(.Machine$double.eps)
        ugrid <- pmin(1 - eps, pmax(eps, xgrid))
        vgrid <- pmin(1 - eps, pmax(eps, ygrid))
    } else if (scale == "Gaussian") {
        ugrid <- stats::pnorm(xgrid)
        vgrid <- stats::pnorm(ygrid)
    } else {
        stop("Direct copula-density evaluation is available only on the Uniform or Gaussian scale")
    }

    gu <- rep(ugrid, times = length(vgrid))
    gv <- rep(vgrid, each = length(ugrid))
    rho_info <- .GeoScatterplot_copula_rho(setup, lag, time_lag)

    ans <- .C(
        "GeoCopulaDensity",
        type_cop = as.integer(setup$type_cop),
        copula_par = as.double(setup$copula_par),
        u = as.double(gu),
        v = as.double(gv),
        n = as.integer(length(gu)),
        rho = as.double(rho_info$rho),
        res = double(length(gu)),
        PACKAGE = "GeoModels",
        NAOK = TRUE
    )$res

    copula_density <- matrix(
        ans, nrow = length(ugrid), ncol = length(vgrid)
    )
    copula_density[
        !is.finite(copula_density) | copula_density < 0
    ] <- 0

    attr(copula_density, "rho") <- rho_info$rho
    attr(copula_density, "rho_latent") <- rho_info$rho_latent
    attr(copula_density, "nugget") <- rho_info$nugget
    attr(copula_density, "copula_parameter") <- setup$copula_par
    attr(copula_density, "copula") <- setup$copula_name

    if (scale == "Uniform") return(copula_density)

    gaussian_jacobian <- outer(
        stats::dnorm(xgrid), stats::dnorm(ygrid)
    )
    gaussian_density <- copula_density * gaussian_jacobian
    gaussian_density[
        !is.finite(gaussian_density) | gaussian_density < 0
    ] <- 0
    attr(gaussian_density, "rho") <- rho_info$rho
    attr(gaussian_density, "rho_latent") <- rho_info$rho_latent
    attr(gaussian_density, "nugget") <- rho_info$nugget
    attr(gaussian_density, "copula_parameter") <- setup$copula_par
    attr(gaussian_density, "copula") <- setup$copula_name
    gaussian_density
}

.GeoScatterplot_independence_density_grid <- function(setup, xgrid, ygrid,
                                                              mean1, mean2, scale) {
    if (scale == "Original") {
        fx <- .GeoScatterplot_marginal_density(
            x = xgrid, mean = mean1, model = setup$model_name,
            nuisance = setup$nuisance_named
        )
        fy <- .GeoScatterplot_marginal_density(
            x = ygrid, mean = mean2, model = setup$model_name,
            nuisance = setup$nuisance_named
        )
        z <- outer(fx, fy)
    } else if (scale == "Uniform") {
        ## Under independence, the copula density on (0,1)^2 is identically 1.
        z <- matrix(1, nrow = length(xgrid), ncol = length(ygrid))
    } else if (scale == "Gaussian") {
        ## After the marginal PIT and Gaussian-score transformation, an
        ## independence fit has joint density phi(z1) * phi(z2).
        z <- outer(stats::dnorm(xgrid), stats::dnorm(ygrid))
    } else {
        stop("Unsupported scale for an independence contour")
    }

    z[!is.finite(z) | z < 0] <- 0
    attr(z, "rho") <- 0
    attr(z, "rho_latent") <- 0
    attr(z, "nugget") <- NA_real_
    attr(z, "copula") <- "Independence"
    attr(z, "copula_parameter") <- NA_real_
    z
}

.GeoScatterplot_density_grid <- function(setup, xgrid, ygrid, lag,
                                         time_lag = 0, mean1, mean2, scale) {
    if (isTRUE(setup$independence)) {
        return(.GeoScatterplot_independence_density_grid(
            setup = setup, xgrid = xgrid, ygrid = ygrid,
            mean1 = mean1, mean2 = mean2, scale = scale
        ))
    }

    if (scale == "Original") {
        return(.GeoScatterplot_response_density_grid(
            setup = setup, xgrid = xgrid, ygrid = ygrid,
            lag = lag, time_lag = time_lag, mean1 = mean1, mean2 = mean2
        ))
    }

    ## For fitted copula models, evaluate the copula density directly.
    ## This is both simpler and more stable than transforming the complete
    ## response density back through the fitted marginal distribution.
    if (setup$type_cop > 0L) {
        return(.GeoScatterplot_direct_copula_grid(
            setup = setup, xgrid = xgrid, ygrid = ygrid,
            lag = lag, time_lag = time_lag, scale = scale
        ))
    }

    ## Native non-Gaussian random fields do not have a separate copula
    ## evaluator. For these models, recover the copula density from the
    ## fitted response density and fitted marginal densities.
    if (scale == "Uniform") {
        eps <- sqrt(.Machine$double.eps)
        ugrid <- pmin(1 - eps, pmax(eps, xgrid))
        vgrid <- pmin(1 - eps, pmax(eps, ygrid))
    } else {
        ugrid <- stats::pnorm(xgrid)
        vgrid <- stats::pnorm(ygrid)
    }

    x_response <- .GeoScatterplot_marginal_quantile(
        u = ugrid, mean = mean1, model = setup$model_name,
        nuisance = setup$nuisance_named
    )
    y_response <- .GeoScatterplot_marginal_quantile(
        u = vgrid, mean = mean2, model = setup$model_name,
        nuisance = setup$nuisance_named
    )

    response_density <- .GeoScatterplot_response_density_grid(
        setup = setup,
        xgrid = x_response,
        ygrid = y_response,
        lag = lag,
        time_lag = time_lag,
        mean1 = mean1,
        mean2 = mean2
    )

    fx <- .GeoScatterplot_marginal_density(
        x = x_response, mean = mean1, model = setup$model_name,
        nuisance = setup$nuisance_named
    )
    fy <- .GeoScatterplot_marginal_density(
        x = y_response, mean = mean2, model = setup$model_name,
        nuisance = setup$nuisance_named
    )

    denominator <- outer(fx, fy)
    copula_density <- response_density / denominator
    copula_density[!is.finite(copula_density) | copula_density < 0] <- 0
    attr(copula_density, "rho") <- attr(response_density, "rho")
    attr(copula_density, "rho_latent") <- attr(response_density, "rho_latent")
    attr(copula_density, "nugget") <- attr(response_density, "nugget")

    if (scale == "Uniform") return(copula_density)

    gaussian_jacobian <- outer(stats::dnorm(xgrid), stats::dnorm(ygrid))
    gaussian_density <- copula_density * gaussian_jacobian
    gaussian_density[!is.finite(gaussian_density) | gaussian_density < 0] <- 0
    attr(gaussian_density, "rho") <- attr(copula_density, "rho")
    attr(gaussian_density, "rho_latent") <- attr(copula_density, "rho_latent")
    attr(gaussian_density, "nugget") <- attr(copula_density, "nugget")
    gaussian_density
}

.GeoScatterplot_panel <- function(x, y, main, lag, time_lag = 0,
                                  setup, contour, mean1, mean2, scale,
                                  ngrid, nlevels, levels,
                                  contour.col, contour.lwd, contour.labels,
                                  point.col, dots) {
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]

    if (!length(x)) {
        plot.new()
        title(main = main)
        mtext("No pairs", side = 3, line = -2)
        return(list(
            x = x, y = y, lag = lag, time_lag = time_lag,
            density = NULL, xgrid = NULL, ygrid = NULL, levels = NULL
        ))
    }

    plot_args <- modifyList(
        list(
            x = x, y = y, col = point.col,
            xlab = "", ylab = "", main = main
        ),
        dots
    )
    do.call(graphics::plot, plot_args)
    graphics::abline(0, 1)

    if (!isTRUE(contour)) {
        return(list(
            x = x, y = y, lag = lag, time_lag = time_lag,
            density = NULL, xgrid = NULL, ygrid = NULL, levels = NULL
        ))
    }

    usr <- graphics::par("usr")
    if (scale == "Uniform") {
        eps <- 1 / (100 * ngrid)
        xlim <- c(max(0, usr[1L]), min(1, usr[2L]))
        ylim <- c(max(0, usr[3L]), min(1, usr[4L]))
        xgrid <- seq(max(eps, xlim[1L]), min(1 - eps, xlim[2L]),
                     length.out = ngrid)
        ygrid <- seq(max(eps, ylim[1L]), min(1 - eps, ylim[2L]),
                     length.out = ngrid)
    } else {
        xgrid <- seq(usr[1L], usr[2L], length.out = ngrid)
        ygrid <- seq(usr[3L], usr[4L], length.out = ngrid)
    }

    z <- .GeoScatterplot_density_grid(
        setup = setup,
        xgrid = xgrid,
        ygrid = ygrid,
        lag = lag,
        time_lag = time_lag,
        mean1 = mean1,
        mean2 = mean2,
        scale = scale
    )

    lev <- levels
    if (is.null(lev)) {
        lev <- .GeoScatterplot_contour_levels(z, nlevels, scale = scale)
    }

    if (length(lev)) {
        graphics::contour(
            xgrid, ygrid, z,
            levels = lev,
            add = TRUE,
            drawlabels = contour.labels,
            col = contour.col,
            lwd = contour.lwd
        )
    }

    rho <- attr(z, "rho")
    rho_latent <- attr(z, "rho_latent")
    nugget <- attr(z, "nugget")
    copula_parameter <- attr(z, "copula_parameter")
    copula_name <- attr(z, "copula")

    reflection <- NULL
    if (scale %in% c("Uniform", "Gaussian") &&
        nrow(z) == ncol(z) && length(xgrid) == length(ygrid)) {
        z_reflected <- z[nrow(z):1L, ncol(z):1L, drop = FALSE]
        delta <- abs(z - z_reflected)
        zmax <- max(z, na.rm = TRUE)
        zsum <- sum(z, na.rm = TRUE)
        reflection <- c(
            max_absolute = max(delta, na.rm = TRUE),
            max_relative = if (is.finite(zmax) && zmax > 0) {
                max(delta, na.rm = TRUE) / zmax
            } else NA_real_,
            L1_relative = if (is.finite(zsum) && zsum > 0) {
                sum(delta, na.rm = TRUE) / zsum
            } else NA_real_
        )
    }

    list(
        x = x, y = y, lag = lag, time_lag = time_lag,
        density = z, xgrid = xgrid, ygrid = ygrid, levels = lev,
        rho = rho, rho_latent = rho_latent, nugget = nugget,
        copula = copula_name, copula_parameter = copula_parameter,
        reflection = reflection
    )
}

.GeoScatterplot_GeoFit <- function(object, maxdist, neighb, time.lag,
                                   numbins, radius, contour, residuals,
                                   scale, gaussian.range, ngrid, nlevels,
                                   levels, contour.col, contour.lwd,
                                   contour.labels, point.col, lag.method,
                                   dots) {
    if (!inherits(object, "GeoFit")) stop("object must inherit from 'GeoFit'")
    if (isTRUE(object$bivariate)) stop("Bivariate GeoFit objects are not supported yet")

    spacetime <- isTRUE(object$spacetime)
    if (!spacetime && !is.null(time.lag)) {
        stop("time.lag can be used only with a spatio-temporal GeoFit object")
    }
    model <- .GeoScatterplot_normalize_model(object$model)
    scale <- .GeoScatterplot_match_scale(scale)
    independence <- identical(object$likelihood, "Marginal") &&
        identical(object$type, "Independence")

    if (scale != "Original" && identical(object$residual_type, "Pearson")) {
        stop(
            "Uniform and Gaussian scatterplot scales are not defined for Pearson-residual objects. Use the original discrete GeoFit object only for discrete PIT diagnostics, or use scale='Original' for Pearson residual scatterplots.",
            call. = FALSE
        )
    }
    if (scale != "Original" && !(model %in% .GeoScatterplot_continuous_models)) {
        stop(
            "Uniform and Gaussian scatterplot scales require a fitted continuous marginal model.",
            call. = FALSE
        )
    }

    graphics::par(pty = "s")

    if (length(gaussian.range) != 1L || !is.finite(gaussian.range) ||
        gaussian.range <= 0) {
        stop("gaussian.range must be a positive finite number")
    }

    if (scale != "Original" && isTRUE(residuals)) {
        stop("residuals can be used only with scale = 'Original'")
    }

    if (isTRUE(contour)) {
        if (isTRUE(object$missp)) {
            stop("Contours are not available for misspecified likelihood models")
        }

        if (isTRUE(independence)) {
            if (!(model %in% unique(c(
                .GeoScatterplot_native_models,
                .GeoScatterplot_copula_models
            )))) {
                stop(
                    "Independence contours are currently available only for ",
                    "continuous marginal models supported by GeoScatterplot"
                )
            }
        } else if (is.null(object$copula)) {
            if (!(model %in% .GeoScatterplot_native_models)) {
                stop(
                    "Native fitted contours are available only for: ",
                    paste(.GeoScatterplot_native_models, collapse = ", ")
                )
            }
        } else {
            if (!(object$copula %in% c("Gaussian", "SkewGaussian", "Clayton"))) {
                stop("Unsupported copula in the GeoFit object")
            }
            if (!(model %in% .GeoScatterplot_copula_models)) {
                stop(
                    "Copula contours are currently available for the continuous ",
                    "marginal models: ",
                    paste(.GeoScatterplot_copula_models, collapse = ", ")
                )
            }
        }
    }

    if (length(numbins) != 1L || is.na(numbins) || numbins < 1L ||
        numbins != as.integer(numbins)) {
        stop("numbins must be a positive integer")
    }
    numbins <- as.integer(numbins)

    if (isTRUE(contour)) {
        if (length(ngrid) != 1L || is.na(ngrid) || ngrid < 20L ||
            ngrid != as.integer(ngrid)) {
            stop("ngrid must be an integer not smaller than 20")
        }
        ngrid <- as.integer(ngrid)

        if (length(nlevels) != 1L || is.na(nlevels) || nlevels < 1L ||
            nlevels != as.integer(nlevels)) {
            stop("nlevels must be a positive integer")
        }
        nlevels <- as.integer(nlevels)
    }

    lag.method <- match.arg(lag.method, choices = c("median", "mean"))

    layout <- .GeoScatterplot_spacetime_layout(object)

    .GeoValidate3D(
        coordx = object$coordx, coordy = object$coordy,
        coordz = object$coordz, coordx_dyn = object$coordx_dyn,
        distance = object$distance, grid = isTRUE(object$grid),
        context = "GeoScatterplot"
    )
    yraw <- layout$data
    mu <- .GeoScatterplot_fitted_location(object, length(yraw))

    if (scale == "Original") {
        yplot <- if (isTRUE(residuals)) {
            .GeoScatterplot_residual_scale(yraw, mu, model, object = object)
        } else {
            yraw
        }
    } else {
        yplot <- .GeoScatterplot_pit_data(object, scale = scale)
        if (length(yplot) != length(yraw)) {
            stop("GeoPit() returned a transformed data vector with an invalid length")
        }
    }

    if (isTRUE(contour) && scale == "Original" && !isTRUE(residuals)) {
        mu_range <- range(mu, finite = TRUE)
        if (length(mu_range) != 2L || diff(mu_range) > sqrt(.Machine$double.eps)) {
            stop(
                "A single original-scale contour requires a constant fitted location; ",
                "use residuals = TRUE or select scale = 'Uniform' or 'Gaussian'"
            )
        }
    }

    contour_mean <- if (scale == "Original" && !isTRUE(residuals)) mu[1L] else 0

    if (is.null(neighb) && is.null(maxdist)) {
        if (!is.null(object$neighb) && length(object$neighb)) {
            neighb <- object$neighb
        } else if (!is.null(object$maxdist) && length(object$maxdist) == 1L &&
                   is.finite(object$maxdist)) {
            maxdist <- object$maxdist
        } else {
            stop("Specify either neighb or maxdist")
        }
    }

    setup <- if (isTRUE(contour)) {
        .GeoScatterplot_density_setup(object, model)
    } else {
        NULL
    }

    panel_dots <- dots
    if (scale == "Uniform") {
        if (is.null(panel_dots$xlim)) panel_dots$xlim <- c(0, 1)
        if (is.null(panel_dots$ylim)) panel_dots$ylim <- c(0, 1)
        if (is.null(panel_dots$xaxs)) panel_dots$xaxs <- "i"
        if (is.null(panel_dots$yaxs)) panel_dots$yaxs <- "i"
        if (is.null(panel_dots$asp)) panel_dots$asp <- 1
    }
    if (scale == "Gaussian") {
        if (is.null(panel_dots$xlim)) {
            panel_dots$xlim <- c(-gaussian.range, gaussian.range)
        }
        if (is.null(panel_dots$ylim)) {
            panel_dots$ylim <- c(-gaussian.range, gaussian.range)
        }
        if (is.null(panel_dots$xaxs)) panel_dots$xaxs <- "i"
        if (is.null(panel_dots$yaxs)) panel_dots$yaxs <- "i"
        if (is.null(panel_dots$asp)) panel_dots$asp <- 1
    }

    if (!is.null(object$anisopars)) {
        if (tolower(object$distance) != "eucl") {
            stop("Anisotropy with a non-Euclidean distance is not supported")
        }
        if (layout$dynamic) {
            layout$coords_dyn <- lapply(
                layout$coords_dyn,
                function(z) .GeoScatterplot_three_columns(GeoAniso(z, anisopars = object$anisopars))
            )
            layout$coords_obs <- do.call(rbind, layout$coords_dyn)
        } else {
            layout$coords_base <- .GeoScatterplot_three_columns(
                GeoAniso(layout$coords_base, anisopars = object$anisopars)
            )
            if (spacetime) {
                nsite <- nrow(layout$coords_base)
                layout$coords_obs <- layout$coords_base[
                    rep(seq_len(nsite), times = length(layout$time_values)),
                    , drop = FALSE
                ]
            } else {
                layout$coords_obs <- layout$coords_base
            }
        }
    }

    time_lags <- .GeoScatterplot_time_lags(object, time.lag)
    temporal_tol <- sqrt(.Machine$double.eps) *
        max(1, abs(layout$time_values), abs(time_lags))

    panel_output <- list()
    panel_index <- 1L

    if (!spacetime) {
        coords_work <- layout$coords_obs

        if (!is.null(neighb)) {
            if (any(!is.finite(neighb)) || any(neighb < 1L) ||
                any(neighb != as.integer(neighb))) {
                stop("neighb must contain positive integers")
            }
            neighb <- as.integer(neighb)

            nc <- ceiling(sqrt(length(neighb)))
            nr <- ceiling(length(neighb) / nc)
            graphics::par(mfrow = c(nr, nc))

            neigh_coords <- .GeoScatterplot_neighbor_coordinates(
                coords_work, object$distance
            )
            for (k in seq_along(neighb)) {
                sel <- GeoNeighIndex(
                    coordx = neigh_coords, distance = object$distance,
                    neighb = neighb[k], radius = radius
                )
                i <- as.integer(sel$colidx)
                j <- as.integer(sel$rowidx)
                lagv <- as.numeric(sel$lags)
                lag <- .GeoScatterplot_representative_lag(lagv, lag.method)

                main <- paste0("Neighb = ", neighb[k], "; h = ", signif(lag, 3))
                panel_output[[panel_index]] <- .GeoScatterplot_panel(
                    x = yplot[i], y = yplot[j], main = main,
                    lag = lag, time_lag = 0, setup = setup,
                    contour = contour, mean1 = contour_mean,
                    mean2 = contour_mean, scale = scale, ngrid = ngrid,
                    nlevels = nlevels, levels = levels,
                    contour.col = contour.col, contour.lwd = contour.lwd,
                    contour.labels = contour.labels, point.col = point.col,
                    dots = panel_dots
                )
                panel_output[[panel_index]]$indices <- cbind(colidx = i, rowidx = j)
                panel_output[[panel_index]]$pair_lags <- lagv
                panel_output[[panel_index]]$pair_time_lags <- rep(0, length(lagv))
                panel_index <- panel_index + 1L
            }
        } else {
            if (length(maxdist) != 1L || !is.finite(maxdist) || maxdist <= 0) {
                stop("maxdist must be a positive finite scalar")
            }

            bins <- seq(0, maxdist, length.out = numbins + 1L)
            pairs <- .GeoScatterplot_spatial_pairs_maxdist(
                data = yplot, coords = coords_work, maxdist = maxdist,
                bins = bins, distance = object$distance, radius = radius
            )

            nc <- ceiling(sqrt(numbins))
            nr <- ceiling(numbins / nc)
            graphics::par(mfrow = c(nr, nc))

            for (k in seq_len(numbins)) {
                use <- pairs$bin == bins[k]
                lag <- (bins[k] + bins[k + 1L]) / 2
                main <- paste0(
                    .GeoScatterplot_distance_label(bins, k),
                    "; h = ", signif(lag, 3)
                )

                panel_output[[panel_index]] <- .GeoScatterplot_panel(
                    x = pairs$v1[use], y = pairs$v2[use], main = main,
                    lag = lag, time_lag = 0, setup = setup,
                    contour = contour, mean1 = contour_mean,
                    mean2 = contour_mean, scale = scale, ngrid = ngrid,
                    nlevels = nlevels, levels = levels,
                    contour.col = contour.col, contour.lwd = contour.lwd,
                    contour.labels = contour.labels, point.col = point.col,
                    dots = panel_dots
                )
                panel_output[[panel_index]]$indices <- cbind(
                    colidx = pairs$i[use], rowidx = pairs$j[use]
                )
                panel_output[[panel_index]]$pair_lags <- pairs$lags[use]
                panel_output[[panel_index]]$pair_time_lags <- rep(0, sum(use))
                panel_index <- panel_index + 1L
            }
        }
    } else {
        if (!is.null(neighb)) {
            if (any(!is.finite(neighb)) || any(neighb < 1L) ||
                any(neighb != as.integer(neighb))) {
                stop(
                    "For spatio-temporal objects, ",
                    "'neighb' must contain positive integers"
                )
            }
            neighb <- as.integer(neighb)

            npanels <- length(neighb) * length(time_lags)
            nc <- ceiling(sqrt(npanels))
            nr <- ceiling(npanels / nc)
            graphics::par(mfrow = c(nr, nc))

            max_time <- if (length(time_lags)) max(time_lags) else 0
            selection_cache <- vector("list", length(neighb))

            for (kk in seq_along(neighb)) {
                K <- neighb[kk]
                spatial_input <- if (layout$dynamic) NULL else
                    .GeoScatterplot_neighbor_coordinates(
                        layout$coords_base, object$distance
                    )
                dynamic_input <- if (layout$dynamic)
                    .GeoScatterplot_neighbor_coordinates(
                        layout$coords_dyn, object$distance
                    ) else NULL
                sel <- GeoNeighIndex(
                    coordx = spatial_input,
                    coordt = layout$time_values,
                    coordx_dyn = dynamic_input,
                    distance = object$distance,
                    neighb = K,
                    maxdist = NULL,
                    maxtime = max_time,
                    radius = radius,
                    bivariate = FALSE
                )
                valid <- is.finite(sel$colidx) & is.finite(sel$rowidx) &
                    sel$colidx >= 1L & sel$rowidx >= 1L &
                    sel$colidx <= length(yplot) & sel$rowidx <= length(yplot) &
                    is.finite(sel$lags) & is.finite(sel$lagt)
                selection_cache[[kk]] <- list(
                    colidx = as.integer(sel$colidx[valid]),
                    rowidx = as.integer(sel$rowidx[valid]),
                    lags = as.numeric(sel$lags[valid]),
                    lagt = as.numeric(sel$lagt[valid])
                )
            }

            for (u in time_lags) {
                for (kk in seq_along(neighb)) {
                    sel <- selection_cache[[kk]]
                    i_all <- sel$colidx
                    j_all <- sel$rowidx
                    lag_all <- sel$lags
                    time_all <- sel$lagt
                    use <- abs(time_all - u) <= temporal_tol

                    i <- i_all[use]
                    j <- j_all[use]
                    lagv <- lag_all[use]
                    timev <- time_all[use]
                    lag <- .GeoScatterplot_representative_lag(lagv, lag.method, empty = 0)
                    u_rep <- .GeoScatterplot_representative_lag(timev, lag.method, empty = u)

                    main <- paste0(
                        "Neighb = ", neighb[kk], "; h = ", signif(lag, 3),
                        "; u = ", signif(u_rep, 3)
                    )

                    panel_output[[panel_index]] <- .GeoScatterplot_panel(
                        x = yplot[i], y = yplot[j], main = main,
                        lag = lag, time_lag = u_rep, setup = setup,
                        contour = contour, mean1 = contour_mean,
                        mean2 = contour_mean, scale = scale, ngrid = ngrid,
                        nlevels = nlevels, levels = levels,
                        contour.col = contour.col, contour.lwd = contour.lwd,
                        contour.labels = contour.labels, point.col = point.col,
                        dots = panel_dots
                    )
                    panel_output[[panel_index]]$indices <- cbind(colidx = i, rowidx = j)
                    panel_output[[panel_index]]$pair_lags <- lagv
                    panel_output[[panel_index]]$pair_time_lags <- timev
                    panel_output[[panel_index]]$neighb <- neighb[kk]
                    panel_index <- panel_index + 1L
                }
            }
        } else {
            if (length(maxdist) != 1L || !is.finite(maxdist) || maxdist <= 0) {
                stop("maxdist must be a positive finite scalar")
            }

            max_time <- if (length(time_lags)) max(time_lags) else 0
            K_all <- max(1L, max(layout$ns))
            sel <- GeoNeighIndex(
                coordx = if (layout$dynamic) NULL else
                    .GeoScatterplot_neighbor_coordinates(
                        layout$coords_base, object$distance
                    ),
                coordt = layout$time_values,
                coordx_dyn = if (layout$dynamic)
                    .GeoScatterplot_neighbor_coordinates(
                        layout$coords_dyn, object$distance
                    ) else NULL,
                distance = object$distance,
                neighb = K_all,
                maxdist = maxdist,
                maxtime = max_time,
                radius = radius,
                bivariate = FALSE
            )

            i_all <- as.integer(sel$colidx)
            j_all <- as.integer(sel$rowidx)
            lag_all <- as.numeric(sel$lags)
            time_all <- as.numeric(sel$lagt)
            valid <- is.finite(i_all) & is.finite(j_all) &
                i_all >= 1L & j_all >= 1L &
                i_all <= length(yplot) & j_all <= length(yplot) &
                is.finite(lag_all) & is.finite(time_all) & lag_all <= maxdist
            i_all <- i_all[valid]
            j_all <- j_all[valid]
            lag_all <- lag_all[valid]
            time_all <- time_all[valid]
            bins <- seq(0, maxdist, length.out = numbins + 1L)

            npanels <- numbins * length(time_lags)
            nc <- ceiling(sqrt(npanels))
            nr <- ceiling(npanels / nc)
            graphics::par(mfrow = c(nr, nc))

            for (u in time_lags) {
                time_use <- abs(time_all - u) <= temporal_tol
                for (k in seq_len(numbins)) {
                    spatial_use <- .GeoScatterplot_bin_membership(
                        lag_all, bins, k
                    )
                    use <- time_use & spatial_use
                    lag <- (bins[k] + bins[k + 1L]) / 2
                    main <- paste0(
                        .GeoScatterplot_distance_label(bins, k),
                        "; h = ", signif(lag, 3),
                        "; u = ", signif(u, 3)
                    )

                    panel_output[[panel_index]] <- .GeoScatterplot_panel(
                        x = yplot[i_all[use]], y = yplot[j_all[use]],
                        main = main, lag = lag, time_lag = u, setup = setup,
                        contour = contour, mean1 = contour_mean,
                        mean2 = contour_mean, scale = scale, ngrid = ngrid,
                        nlevels = nlevels, levels = levels,
                        contour.col = contour.col, contour.lwd = contour.lwd,
                        contour.labels = contour.labels, point.col = point.col,
                        dots = panel_dots
                    )
                    panel_output[[panel_index]]$indices <- cbind(
                        colidx = i_all[use], rowidx = j_all[use]
                    )
                    panel_output[[panel_index]]$pair_lags <- lag_all[use]
                    panel_output[[panel_index]]$pair_time_lags <- time_all[use]
                    panel_index <- panel_index + 1L
                }
            }
        }
    }

    invisible(list(
        call = match.call(), model = model,
        copula = if (isTRUE(independence)) "Independence" else object$copula,
        independence = isTRUE(independence),
        contour = isTRUE(contour), scale = scale,
        residuals = isTRUE(residuals), spacetime = spacetime,
        time_lags = time_lags, panels = panel_output
    ))
}

.GeoScatterplot_raw_time_indices <- function(coordt, times = NULL) {
    coordt <- as.numeric(coordt)
    if (length(coordt) < 2L || any(!is.finite(coordt))) {
        stop("coordt must contain at least two finite temporal coordinates")
    }
    if (anyDuplicated(coordt)) {
        stop("coordt must not contain duplicated temporal coordinates")
    }

    if (is.null(times)) return(seq_along(coordt))
    if (!is.numeric(times) || !length(times) || any(!is.finite(times))) {
        stop("times must contain finite temporal coordinates or valid time indices")
    }

    tol <- sqrt(.Machine$double.eps) * max(1, abs(coordt), abs(times))
    idx <- vapply(times, function(tt) {
        jj <- which.min(abs(coordt - tt))
        if (length(jj) && abs(coordt[jj] - tt) <= tol) return(as.integer(jj))

        if (tt == as.integer(tt) && tt >= 1L && tt <= length(coordt)) {
            return(as.integer(tt))
        }

        stop(
            "The temporal instant ", signif(tt, 8),
            " is neither present in coordt nor a valid time index"
        )
    }, integer(1L))

    idx <- sort(unique(idx))
    if (length(idx) < 2L) {
        stop("times must select at least two temporal instants")
    }
    idx
}

.GeoScatterplot_raw_time_lags <- function(coordt, time.lag = NULL) {
    coordt <- as.numeric(coordt)
    possible <- .GeoScatterplot_distinct_time_lags(coordt)
    tol <- sqrt(.Machine$double.eps) * max(1, abs(coordt), abs(possible))

    if (is.null(time.lag)) return(possible)
    if (!is.numeric(time.lag) || !length(time.lag) ||
        any(!is.finite(time.lag)) || any(time.lag < 0)) {
        stop("time.lag must contain non-negative finite temporal lags")
    }

    matched <- vapply(time.lag, function(u) {
        jj <- which.min(abs(possible - u))
        local_tol <- max(tol, sqrt(.Machine$double.eps) * max(1, abs(u)))
        if (!length(jj) || abs(possible[jj] - u) > local_tol) {
            stop(
                "Temporal lag ", signif(u, 8),
                " is not available from the selected temporal coordinates"
            )
        }
        possible[jj]
    }, numeric(1L))

    .GeoScatterplot_unique_tol(matched, tol = tol)
}

.GeoScatterplot_raw_spacetime_layout <- function(data, coords = NULL,
                                                  coordt, coordx_dyn = NULL,
                                                  times = NULL) {
    coordt <- as.numeric(coordt)
    time_idx <- .GeoScatterplot_raw_time_indices(coordt, times)
    coordt_use <- coordt[time_idx]
    nt <- length(coordt)

    if (is.null(coordx_dyn)) {
        coords <- .GeoScatterplot_three_columns(coords, validate = TRUE)
        nsite <- nrow(coords)

        if (is.matrix(data) || is.data.frame(data)) {
            data_mat <- as.matrix(data)
            if (!identical(dim(data_mat), c(nt, nsite))) {
                stop(
                    "For fixed spatio-temporal locations, data must be a ",
                    "length(coordt) by nrow(coordx) matrix"
                )
            }
        } else {
            data_vec <- as.numeric(data)
            if (length(data_vec) != nt * nsite) {
                stop(
                    "For fixed spatio-temporal locations, length(data) must equal ",
                    "length(coordt) * nrow(coordx)"
                )
            }
            data_mat <- matrix(data_vec, nrow = nt, ncol = nsite, byrow = TRUE)
        }

        data_mat <- data_mat[time_idx, , drop = FALSE]
        y <- as.numeric(t(data_mat))
        ns <- rep.int(nsite, length(time_idx))
        coords_obs <- coords[rep(seq_len(nsite), times = length(time_idx)),
                             , drop = FALSE]
        time_obs <- rep(coordt_use, each = nsite)

        return(list(
            data = y,
            coords_base = coords,
            coords_dyn = NULL,
            coords_obs = coords_obs,
            time_obs = time_obs,
            time_values = coordt_use,
            ns = ns,
            dynamic = FALSE
        ))
    }

    if (!is.list(coordx_dyn) || length(coordx_dyn) != nt) {
        stop("coordx_dyn must be a list with one coordinate matrix per temporal instant")
    }
    coords_dyn_all <- lapply(coordx_dyn, .GeoScatterplot_three_columns, validate = TRUE)
    ns_all <- vapply(coords_dyn_all, nrow, integer(1L))
    if (any(ns_all <= 0L)) {
        stop("Each element of coordx_dyn must contain at least one location")
    }

    if (is.list(data)) {
        if (length(data) != nt) {
            stop("For dynamic locations, data must have one component per temporal instant")
        }
        data_all <- lapply(seq_len(nt), function(k) {
            z <- as.numeric(data[[k]])
            if (length(z) != ns_all[k]) {
                stop("data[[", k, "]] and coordx_dyn[[", k, "]] have different lengths")
            }
            z
        })
    } else if ((is.matrix(data) || is.data.frame(data)) &&
               length(unique(ns_all)) == 1L &&
               identical(dim(as.matrix(data)), c(nt, ns_all[1L]))) {
        data_mat <- as.matrix(data)
        data_all <- lapply(seq_len(nt), function(k) as.numeric(data_mat[k, ]))
    } else {
        data_vec <- as.numeric(data)
        if (length(data_vec) != sum(ns_all)) {
            stop(
                "For dynamic locations, length(data) must equal the total ",
                "number of rows in coordx_dyn"
            )
        }
        ends <- cumsum(ns_all)
        starts <- c(1L, ends[-length(ends)] + 1L)
        data_all <- Map(function(a, b) data_vec[a:b], starts, ends)
    }

    coords_dyn <- coords_dyn_all[time_idx]
    data_use <- data_all[time_idx]
    ns <- ns_all[time_idx]
    y <- as.numeric(unlist(data_use, use.names = FALSE))
    coords_obs <- do.call(rbind, coords_dyn)
    time_obs <- rep(coordt_use, times = ns)

    list(
        data = y,
        coords_base = NULL,
        coords_dyn = coords_dyn,
        coords_obs = coords_obs,
        time_obs = time_obs,
        time_values = coordt_use,
        ns = ns,
        dynamic = TRUE
    )
}

.GeoScatterplot_raw_spacetime <- function(data, coords = NULL, coordt,
                                          coordx_dyn = NULL, times = NULL,
                                          time.lag = NULL, distance = "Eucl",
                                          maxdist = NULL, neighb = NULL,
                                          numbins = 4L, radius = 1,
                                          point.col, lag.method, dots,
                                          call) {
    layout <- .GeoScatterplot_raw_spacetime_layout(
        data = data, coords = coords, coordt = coordt,
        coordx_dyn = coordx_dyn, times = times
    )
    time_lags <- .GeoScatterplot_raw_time_lags(
        layout$time_values, time.lag = time.lag
    )
    temporal_tol <- sqrt(.Machine$double.eps) *
        max(1, abs(layout$time_values), abs(time_lags))

    lag.method <- match.arg(lag.method, choices = c("median", "mean"))
    if (is.null(neighb) && is.null(maxdist)) {
        stop("Specify either neighb or maxdist")
    }
    if (!is.null(neighb) && !is.null(maxdist)) {
        warning("Both neighb and maxdist were supplied; neighb is used")
    }

    max_time <- if (length(time_lags)) max(time_lags) else 0
    panel_index <- 1L
    panels <- list()

    if (!is.null(neighb)) {
        if (any(!is.finite(neighb)) || any(neighb < 1L) ||
            any(neighb != as.integer(neighb))) {
            stop("neighb must contain positive integers")
        }
        neighb <- as.integer(neighb)
        npanels <- length(neighb) * length(time_lags)
        nc <- ceiling(sqrt(npanels))
        nr <- ceiling(npanels / nc)
        graphics::par(mfrow = c(nr, nc))

        selection_cache <- lapply(neighb, function(K) {
            GeoNeighIndex(
                coordx = if (layout$dynamic) NULL else
                    .GeoScatterplot_neighbor_coordinates(layout$coords_base, distance),
                coordt = layout$time_values,
                coordx_dyn = if (layout$dynamic)
                    .GeoScatterplot_neighbor_coordinates(layout$coords_dyn, distance) else NULL,
                distance = distance,
                neighb = K,
                maxdist = NULL,
                maxtime = max_time,
                radius = radius,
                bivariate = FALSE
            )
        })

        for (u in time_lags) {
            for (kk in seq_along(neighb)) {
                sel <- selection_cache[[kk]]
                i_all <- as.integer(sel$colidx)
                j_all <- as.integer(sel$rowidx)
                lag_all <- as.numeric(sel$lags)
                time_all <- as.numeric(sel$lagt)

                valid <- is.finite(i_all) & is.finite(j_all) &
                    i_all >= 1L & j_all >= 1L &
                    i_all <= length(layout$data) & j_all <= length(layout$data) &
                    is.finite(lag_all) & is.finite(time_all)
                use <- valid & abs(time_all - u) <= temporal_tol

                i <- i_all[use]
                j <- j_all[use]
                lagv <- lag_all[use]
                timev <- time_all[use]
                lag <- .GeoScatterplot_representative_lag(lagv, lag.method, empty = 0)
                u_rep <- .GeoScatterplot_representative_lag(timev, lag.method, empty = u)
                main <- paste0(
                    "Neighb = ", neighb[kk], "; h = ", signif(lag, 3),
                    "; u = ", signif(u_rep, 3)
                )

                panels[[panel_index]] <- .GeoScatterplot_panel(
                    x = layout$data[i], y = layout$data[j], main = main,
                    lag = lag, time_lag = u_rep, setup = NULL,
                    contour = FALSE, mean1 = 0, mean2 = 0,
                    scale = "Original", ngrid = 20L, nlevels = 1L,
                    levels = NULL, contour.col = "red", contour.lwd = 1,
                    contour.labels = FALSE, point.col = point.col,
                    dots = dots
                )
                panels[[panel_index]]$indices <- cbind(colidx = i, rowidx = j)
                panels[[panel_index]]$pair_lags <- lagv
                panels[[panel_index]]$pair_time_lags <- timev
                panels[[panel_index]]$neighb <- neighb[kk]
                panel_index <- panel_index + 1L
            }
        }
    } else {
        if (length(maxdist) != 1L || !is.finite(maxdist) || maxdist <= 0) {
            stop("maxdist must be a positive finite scalar")
        }

        K_all <- max(1L, max(layout$ns))
        sel <- GeoNeighIndex(
            coordx = if (layout$dynamic) NULL else
                .GeoScatterplot_neighbor_coordinates(layout$coords_base, distance),
            coordt = layout$time_values,
            coordx_dyn = if (layout$dynamic)
                .GeoScatterplot_neighbor_coordinates(layout$coords_dyn, distance) else NULL,
            distance = distance,
            neighb = K_all,
            maxdist = maxdist,
            maxtime = max_time,
            radius = radius,
            bivariate = FALSE
        )

        i_all <- as.integer(sel$colidx)
        j_all <- as.integer(sel$rowidx)
        lag_all <- as.numeric(sel$lags)
        time_all <- as.numeric(sel$lagt)
        valid <- is.finite(i_all) & is.finite(j_all) &
            i_all >= 1L & j_all >= 1L &
            i_all <= length(layout$data) & j_all <= length(layout$data) &
            is.finite(lag_all) & is.finite(time_all) & lag_all <= maxdist
        i_all <- i_all[valid]
        j_all <- j_all[valid]
        lag_all <- lag_all[valid]
        time_all <- time_all[valid]

        bins <- seq(0, maxdist, length.out = numbins + 1L)
        npanels <- numbins * length(time_lags)
        nc <- ceiling(sqrt(npanels))
        nr <- ceiling(npanels / nc)
        graphics::par(mfrow = c(nr, nc))

        for (u in time_lags) {
            time_use <- abs(time_all - u) <= temporal_tol
            for (k in seq_len(numbins)) {
                spatial_use <- .GeoScatterplot_bin_membership(
                    lag_all, bins, k
                )
                use <- time_use & spatial_use
                lag <- (bins[k] + bins[k + 1L]) / 2
                main <- paste0(
                    .GeoScatterplot_distance_label(bins, k),
                    "; h = ", signif(lag, 3),
                    "; u = ", signif(u, 3)
                )

                panels[[panel_index]] <- .GeoScatterplot_panel(
                    x = layout$data[i_all[use]], y = layout$data[j_all[use]],
                    main = main, lag = lag, time_lag = u, setup = NULL,
                    contour = FALSE, mean1 = 0, mean2 = 0,
                    scale = "Original", ngrid = 20L, nlevels = 1L,
                    levels = NULL, contour.col = "red", contour.lwd = 1,
                    contour.labels = FALSE, point.col = point.col,
                    dots = dots
                )
                panels[[panel_index]]$indices <- cbind(
                    colidx = i_all[use], rowidx = j_all[use]
                )
                panels[[panel_index]]$pair_lags <- lag_all[use]
                panels[[panel_index]]$pair_time_lags <- time_all[use]
                panel_index <- panel_index + 1L
            }
        }
    }

    invisible(list(
        call = call, spacetime = TRUE, time_lags = time_lags,
        time_values = layout$time_values, panels = panels
    ))
}

GeoScatterplot <- function(data, coordx = NULL, coordy = NULL, coordz = NULL,
                           coordt = NULL, coordx_dyn = NULL,
                           distance = "Eucl", grid = FALSE,
                           maxdist = NULL, neighb = NULL,
                           times = NULL, time.lag = NULL, numbins = 4, radius = 1,
                           bivariate = FALSE,
                           contour = inherits(data, "GeoFit"), residuals = FALSE,
                           scale = c("Original", "Gaussian", "Uniform"),
                           gaussian.range = 3, ngrid = 80, nlevels = 6, levels = NULL,
                           contour.col = "red", contour.lwd = 1.5,
                           contour.labels = FALSE,
                           point.col = grDevices::adjustcolor("#481567FF", 0.45),
                           lag.method = c("median", "mean"), ...) {

    call <- match.call()
    opar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(opar))

    if (inherits(data, "GeoFit")) {
        fit_radius <- if (missing(radius)) data$radius else radius
        if (is.null(fit_radius) || !is.finite(fit_radius)) fit_radius <- 1

        return(.GeoScatterplot_GeoFit(
            object = data,
            maxdist = maxdist,
            neighb = neighb,
            time.lag = time.lag,
            numbins = numbins,
            radius = fit_radius,
            contour = contour,
            residuals = residuals,
            scale = scale,
            gaussian.range = gaussian.range,
            ngrid = ngrid,
            nlevels = nlevels,
            levels = levels,
            contour.col = contour.col,
            contour.lwd = contour.lwd,
            contour.labels = contour.labels,
            point.col = point.col,
            lag.method = lag.method,
            dots = list(...)
        ))
    }

    raw_coord_dim <- .GeoCoordinateDimension(
        coordx = coordx, coordy = coordy, coordz = coordz,
        coordx_dyn = coordx_dyn
    )

    scale <- .GeoScatterplot_match_scale(scale)
    if (scale != "Original") {
        stop("Uniform and Gaussian scales require a GeoFit object as the first argument")
    }

    if (isTRUE(contour)) {
        stop("Fitted contours require a GeoFit object as the first argument")
    }

    CheckDistance(distance)
    spacetime <- !isTRUE(bivariate) && !is.null(coordt) &&
        is.numeric(coordt) && length(coordt) > 1L

    if (length(numbins) != 1 || is.na(numbins) || numbins < 1 ||
        numbins != as.integer(numbins)) {
        stop("numbins must be a positive integer")
    }
    numbins <- as.integer(numbins)

    if (!spacetime && !is.null(time.lag)) {
        stop("time.lag can be used only for spatio-temporal data")
    }
    if (!spacetime && !is.null(times)) {
        stop("times can be used only for spatio-temporal data")
    }

    if (!is.null(coordx_dyn) && spacetime) {
        .GeoValidate3D(
            coordx_dyn = coordx_dyn, distance = distance, grid = grid,
            context = "GeoScatterplot"
        )
        return(.GeoScatterplot_raw_spacetime(
            data = data, coords = NULL, coordt = coordt,
            coordx_dyn = coordx_dyn, times = times,
            time.lag = time.lag, distance = distance,
            maxdist = maxdist, neighb = neighb, numbins = numbins,
            radius = radius, point.col = point.col,
            lag.method = lag.method, dots = list(...), call = call
        ))
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

        if (spacetime) {
            validate_coords <- if (identical(raw_coord_dim, 2L))
                coords[, 1:2, drop = FALSE] else coords
            .GeoValidate3D(
                coordx = validate_coords, distance = distance, grid = grid,
                context = "GeoScatterplot"
            )
            return(.GeoScatterplot_raw_spacetime(
                data = data, coords = coords, coordt = coordt,
                coordx_dyn = NULL, times = times,
                time.lag = time.lag, distance = distance,
                maxdist = maxdist, neighb = neighb, numbins = numbins,
                radius = radius, point.col = point.col,
                lag.method = lag.method, dots = list(...), call = call
            ))
        }

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
        ns <- vapply(coordx_dyn, nrow, integer(1L))
    }

    validate_coords <- if (identical(raw_coord_dim, 2L))
        coords[, 1:2, drop = FALSE] else coords
    .GeoValidate3D(
        coordx = validate_coords, distance = distance, grid = grid,
        context = "GeoScatterplot"
    )

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

    .call_pairs <- function(data, coords, bins, maxdist) {
        .GeoScatterplot_spatial_pairs_maxdist(
            data = data, coords = coords, maxdist = maxdist,
            bins = bins, distance = distance, radius = radius
        )
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
                sel <- vvv$bin == bins[i]
                v111 <- vvv$v1[sel]
                v222 <- vvv$v2[sel]

                main <- .GeoScatterplot_distance_label(bins, i)

                if (length(v111) == 0 || length(v222) == 0) {
                    plot.new()
                    title(main = main)
                    mtext("No pairs", side = 3, line = -2)
                    message("There are no points between ", signif(bins[i], 3),
                            " and ", signif(bins[i + 1], 3))
                } else {
                    plot_args <- modifyList(
                        list(x = v111, y = v222, col = point.col,
                             xlab = "", ylab = "", main = main),
                        list(...)
                    )
                    do.call(graphics::plot, plot_args)
                    graphics::abline(0, 1)
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
                sel <- vvv1$bin == bins1[i]
                v111 <- vvv1$v1[sel]
                v222 <- vvv1$v2[sel]
                cc <- if (length(v111) > 1 && length(v222) > 1) stats::cor(v111, v222) else NA_real_
                main <- paste("cor =", signif(cc, 3),
                              .GeoScatterplot_distance_label(bins1, i))

                if (length(v111) == 0 || length(v222) == 0) {
                    plot.new(); title(main = main); mtext("No pairs", side = 3, line = -2)
                } else {
                    plot_args <- modifyList(
                        list(x = v111, y = v222, col = point.col,
                             xlab = "", ylab = "", main = main),
                        list(...)
                    )
                    do.call(graphics::plot, plot_args)
                    graphics::abline(0, 1)
                }
            }

            for (i in seq_len(numbins)) {
                sel <- vvv2$bin == bins2[i]
                v111 <- vvv2$v1[sel]
                v222 <- vvv2$v2[sel]
                cc <- if (length(v111) > 1 && length(v222) > 1) stats::cor(v111, v222) else NA_real_
                main <- paste("cor =", signif(cc, 3),
                              .GeoScatterplot_distance_label(bins2, i))

                if (length(v111) == 0 || length(v222) == 0) {
                    plot.new(); title(main = main); mtext("No pairs", side = 3, line = -2)
                } else {
                    plot_args <- modifyList(
                        list(x = v111, y = v222, col = point.col,
                             xlab = "", ylab = "", main = main),
                        list(...)
                    )
                    do.call(graphics::plot, plot_args)
                    graphics::abline(0, 1)
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
            neigh_coords <- .GeoScatterplot_neighbor_coordinates(coords, distance)
            for (i in seq_len(ln)) {
                sel <- GeoNeighIndex(
                    coordx = neigh_coords, distance = distance,
                    neighb = neighb[i], radius = radius
                )
                data1 <- data[sel$colidx]
                data2 <- data[sel$rowidx]
                plot_args <- modifyList(
                    list(x = data1, y = data2, col = point.col,
                         xlab = "", ylab = "",
                         main = paste("Neighb =", neighb[i])),
                    list(...)
                )
                do.call(graphics::plot, plot_args)
                graphics::abline(0, 1)
            }
        }

        if (bivariate) {
            stop("Neighborhood scatterplot for bivariate=TRUE is not implemented yet")
        }

    }

    invisible(call)
}
