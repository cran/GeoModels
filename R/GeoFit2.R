####################################################
### File name: GeoFit2.R
####################################################

# GeoFit2 is a starting-value wrapper around GeoFit.
#
# It first tries to estimate the marginal parameters under spatial
# independence and then calls the canonical GeoFit implementation with
# the updated starting values.  If an independence likelihood is not
# available for a model, the supplied starting values are used unchanged.
#

GeoFit2 <- function(data, coordx, coordy = NULL, coordz = NULL, coordt = NULL,
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
                    warn_independence_failure = TRUE)
{
    call <- match.call()

    if (is.null(start))
        stop("Starting parameters are missing")

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

    start_original <- start
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

    can_initialize <-
        isTRUE(independence_start) &&
        type_clean != "Independence" &&
        !is_bivariate

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

            common <- intersect(
                names(start_used),
                names(estimated_marginal)
            )

            if (length(common) > 0L) {
                start_used[common] <- as.list(
                    estimated_marginal[common]
                )
            } else {
                independence_message <- paste0(
                    "The independence fit succeeded for model '",
                    model_clean,
                    "', but it estimated no parameter present in start."
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
    fitted <- do.call(GeoFit, fit_args)

    # Keep the usual GeoFit class and add transparent initialization details.
    fitted$start_original <- start_original
    fitted$start_used <- start_used
    fitted$independence_fit <- independence_fit
    fitted$independence_start_message <- independence_message
    fitted$call <- call

    structure(fitted, class = "GeoFit")
}
