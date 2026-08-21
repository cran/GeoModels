####################################################
### File name: GeoFit2.R
####################################################

# GeoFit2 is a univariate starting-value wrapper around GeoFit.
#
# It first tries to estimate the marginal parameters under spatial
# independence and then calls the canonical GeoFit implementation with
# the updated starting values.  If an independence likelihood is not
# available for a model, the supplied starting values are used unchanged.
#

GeoFit2 <- function(data, coordx = NULL, coordy = NULL, coordz = NULL, coordt = NULL,
                    coordx_dyn = NULL, copula = NULL, corrmodel = NULL,
                    distance = "Eucl", fixed = NULL, anisopars = NULL,
                    est.aniso = c(FALSE, FALSE), grid = FALSE,
                    likelihood = "Marginal", lower = NULL, maxdist = Inf,
                    neighb = NULL, p_neighb = 1, maxtime = Inf,
                    memdist = TRUE, method = "cholesky", model = "Gaussian",
                    n = 1, onlyvar = FALSE, optimizer = "Nelder-Mead",
                    radius = 1, score = FALSE, sensitivity = FALSE,
                    sparse = FALSE, start = NULL,
                    thin_method = "bernoulli", type = "Pairwise",
                    upper = NULL, varest = FALSE, weighted = FALSE,
                    X = NULL, spobj = NULL, spdata = NULL,
                    independence_start = TRUE,
                    independence_optimizer = "Nelder-Mead",
                    warn_independence_failure = TRUE, check.duplicates = FALSE)
{
    call <- match.call()

    if (!isTRUE(memdist)) {
        warning(
            "'memdist = FALSE' is deprecated and ignored; GeoFit2 now always precomputes the selected pair structure.",
            call. = FALSE
        )
    }
    memdist <- TRUE

    start_was_null <- is.null(start)
    if (start_was_null) start <- list()

    # Collect only arguments accepted by GeoFit.  All values are forced
    # here, so do.call() does not depend on expressions in the caller.
    geofit_arg_names <- setdiff(
        names(formals(GeoFit2)),
        c(
            "independence_start",
            "independence_optimizer",
            "warn_independence_failure"
        )
    )

    fit_args <- mget(
        geofit_arg_names,
        envir = environment(),
        inherits = FALSE
    )

    start_original <- if (start_was_null) NULL else start
    start_used <- start
    independence_fit <- NULL
    independence_message <- NULL

    model_clean <- gsub("[[:blank:]]", "", model)
    type_clean <- gsub("[[:blank:]]", "", type)
    corrmodel_clean <- if (is.null(corrmodel)) NULL else
        gsub("[[:blank:]]", "", corrmodel)

    is_bivariate <- FALSE
    if (!is.null(corrmodel_clean)) {
        is_bivariate <- tryCatch(
            isTRUE(CheckBiv(CkCorrModel(corrmodel_clean))),
            error = function(e) FALSE
        )
    }

    if (is_bivariate) {
        stop(
            "GeoFit2 does not support bivariate correlation models. ",
            "Use GeoFit() directly for bivariate fitting.",
            call. = FALSE
        )
    }

    # Mean coefficients estimated by the preliminary independence fit are
    # useful also when the caller omitted them from start: GeoFit/WlsStart
    # will keep them in the estimated-parameter set.  The mapping follows the
    # columns of the user-supplied design matrix exactly; no intercept column
    # is added when X is supplied.
    mean_names <- if (is.null(X)) {
        "mean"
    } else if (is.list(X) && !is.data.frame(X)) {
        p <- unique(vapply(X, function(z) ncol(as.matrix(z)), integer(1L)))
        if (length(p) != 1L || is.na(p) || p < 1L) {
            stop(
                "All blocks of X must have the same positive number of columns.",
                call. = FALSE
            )
        }
        .GeoMean_names(p)
    } else {
        p <- ncol(as.matrix(X))
        if (is.null(p) || is.na(p) || p < 1L) {
            stop("X must contain at least one column.", call. = FALSE)
        }
        .GeoMean_names(p)
    }

    can_initialize <-
        isTRUE(independence_start) &&
        type_clean != "Independence"

    if (can_initialize) {

        ind_args <- fit_args

        # The preliminary fit estimates marginal parameters only.
        ind_args$likelihood <- "Marginal"
        ind_args$type <- "Independence"
        ind_args$optimizer <- independence_optimizer
        ind_args$onlyvar <- FALSE
        ind_args$score <- FALSE
        ind_args$sensitivity <- FALSE
        ind_args$varest <- FALSE
        ind_args$weighted <- FALSE
        ind_args$anisopars <- NULL
        ind_args$est.aniso <- c(FALSE, FALSE)
        ind_args$neighb <- NULL
        ind_args$maxdist <- Inf
        ind_args$maxtime <- Inf
        ind_args$p_neighb <- 1
        ind_args$thin_method <- "bernoulli"

        independence_fit <- tryCatch(
            do.call(GeoFit, ind_args),
            error = function(e) {
                independence_message <<- conditionMessage(e)
                NULL
            }
        )

        if (!is.null(independence_fit)) {
            estimated_marginal <- unlist(
                independence_fit$param,
                use.names = TRUE
            )

            start_names <- names(start_used)
            fixed_names <- names(unlist(fixed, use.names = TRUE))
            if (is.null(start_names)) start_names <- character(0L)
            if (is.null(fixed_names)) fixed_names <- character(0L)

            # Update marginal parameters already supplied in start and add
            # estimated regression coefficients that were intentionally
            # omitted from start.  Fixed parameters always have precedence.
            update_names <- if (start_was_null) {
                ## With start=NULL, use every eligible marginal estimate from
                ## the independence fit. Dependence parameters will keep the
                ## conservative automatic values prepared by GeoFit/WlsStart.
                names(estimated_marginal)
            } else {
                unique(c(
                    intersect(start_names, names(estimated_marginal)),
                    intersect(mean_names, names(estimated_marginal))
                ))
            }
            update_names <- setdiff(update_names, fixed_names)

            # In GeoModels, sill is the total marginal variance and nugget is
            # a correlation attenuation parameter.  Independence estimates
            # may initialize sill, but nugget must remain exactly as supplied
            # by the user and is never inferred from the marginal variance.
            update_names <- update_names[
                !grepl("^nugget($|_)", update_names)
            ]

            finite_update <- vapply(
                estimated_marginal[update_names],
                function(z) length(z) == 1L && is.finite(z),
                logical(1L)
            )
            update_names <- update_names[finite_update]

            if (length(update_names) > 0L) {
                start_used[update_names] <- as.list(
                    estimated_marginal[update_names]
                )
            } else {
                independence_message <- paste0(
                    "The independence fit succeeded for model '",
                    model_clean,
                    "', but it estimated no eligible parameter for the final start."
                )
            }
        }

        if (is.null(independence_fit) &&
            isTRUE(warn_independence_failure) &&
            !is.null(independence_message)) {
            warning(
                "Independence initialization was not available; ",
                "the supplied starting values are being used. ",
                "Reason: ", independence_message,
                call. = FALSE
            )
        }
    }

    fit_args$start <- start_used
    if (!is.null(independence_fit) && isTRUE(independence_fit$duplicate_locations_checked)) {
        fit_args$check.duplicates <- FALSE
    }
    fitted <- do.call(GeoFit, fit_args)

    # Keep the usual GeoFit class and add transparent initialization details.
    fitted$start_original <- start_original
    fitted$start_used <- start_used
    fitted$independence_fit <- independence_fit
    fitted$independence_start_message <- independence_message
    fitted$call <- call

    structure(fitted, class = "GeoFit")
}
