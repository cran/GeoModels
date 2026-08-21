## Robust one-dimensional optimization used by GeoFit composite-likelihood
## paths.  The historical code called stats::optimize() over the complete
## automatic parameter interval.  For positive parameters that interval can
## be enormous (for example [.Machine$double.eps, 1e12]); even after a log
## transformation, optimize() assumes a unimodal objective and can converge
## to a remote numerical basin while ignoring the supplied starting value.
##
## The helper below first brackets a local minimum around the supplied start
## on the natural transformed scale (log scale for positive parameters), then
## refines only inside that bracket with stats::optimize().  This mirrors the
## start-dependent behaviour of the multidimensional optimizers used by
## GeoFit and avoids pathological jumps across many orders of magnitude.
.GeoModels_optimize_1d <- function(f, param, lower, upper, namesparam, ...) {
    if (length(param) != 1L)
        stop(".GeoModels_optimize_1d requires exactly one parameter", call. = FALSE)

    p0 <- as.numeric(param[1L])
    lo <- as.numeric(lower[1L])
    hi <- as.numeric(upper[1L])
    pname <- if (length(namesparam)) as.character(namesparam[1L]) else "parameter"

    if (!is.finite(p0))
        stop("Invalid starting value for one-parameter optimization: ", pname,
             call. = FALSE)

    ## optimize() needs finite bounds.  Normally GeoFit has already supplied
    ## them, but keep this helper safe for legacy/internal calls too.
    if (!is.finite(lo) || !is.finite(hi)) {
        if (p0 > 0 && (is.na(lo) || !is.finite(lo) || lo >= 0)) {
            if (!is.finite(lo)) lo <- max(.Machine$double.eps, p0 / 1e4)
            if (!is.finite(hi)) hi <- max(p0 * 1e4, p0 + 1)
        } else {
            span <- max(10, 20 * (1 + abs(p0)))
            if (!is.finite(lo)) lo <- p0 - span
            if (!is.finite(hi)) hi <- p0 + span
        }
    }

    if (!is.finite(lo) || !is.finite(hi) || lo >= hi)
        stop("Invalid interval for one-parameter optimization: ", pname,
             call. = FALSE)

    use_log <- is.finite(p0) && p0 > 0 && lo >= 0 && hi > 0

    if (use_log) {
        lo_pos <- max(lo, .Machine$double.eps)
        aa <- log(lo_pos)
        bb <- log(hi)
        z0 <- log(min(max(p0, lo_pos), hi))
        back <- exp
        ff <- function(z, ...) f(exp(z), namesparam = namesparam, ...)
        base_step <- log(2)
        parameter_scale <- "log"
    } else {
        aa <- lo
        bb <- hi
        z0 <- min(max(p0, lo), hi)
        back <- identity
        ff <- function(z, ...) f(z, namesparam = namesparam, ...)
        base_step <- max(0.25, 0.25 * (1 + abs(z0)))
        parameter_scale <- "linear"
    }

    if (!is.finite(aa) || !is.finite(bb) || aa >= bb)
        stop("Invalid transformed interval for one-parameter optimization: ",
             pname, call. = FALSE)

    ## Safe scalar objective evaluation for the bracketing stage.  Invalid
    ## trial points are treated as +Inf but genuine errors at the starting
    ## point are still exposed below if no finite bracket can be constructed.
    eval_one <- function(z) {
        val <- tryCatch(ff(z, ...), error = function(e) Inf)
        val <- suppressWarnings(as.numeric(val)[1L])
        if (!is.finite(val)) Inf else val
    }

    f0 <- eval_one(z0)
    if (!is.finite(f0))
        stop("Non-finite objective at the starting value for one-parameter optimization: ",
             pname, call. = FALSE)

    ## Start-centred adaptive bracketing.  The first probes are a factor of 2
    ## apart for positive parameters.  The radius is expanded only when the
    ## objective is still descending, so a nearby likelihood basin is not
    ## abandoned merely because the formal automatic bounds are extremely
    ## broad.
    step <- min(base_step, (bb - aa) / 4)
    if (!is.finite(step) || step <= 0) step <- (bb - aa) / 4

    z_left <- max(aa, z0 - step)
    z_right <- min(bb, z0 + step)
    f_left <- eval_one(z_left)
    f_right <- eval_one(z_right)

    bracket <- NULL

    if (f0 <= f_left && f0 <= f_right) {
        bracket <- c(z_left, z_right)
    } else {
        direction <- if (f_left < f_right) -1 else 1
        z_prev <- z0
        f_prev <- f0
        z_curr <- if (direction < 0) z_left else z_right
        f_curr <- if (direction < 0) f_left else f_right
        radius <- step

        for (ii in seq_len(40L)) {
            if (!is.finite(f_curr) || f_curr >= f_prev) {
                bracket <- sort(c(z_curr, z_prev))
                break
            }

            radius <- radius * 1.6
            z_next <- z0 + direction * radius
            z_next <- min(max(z_next, aa), bb)

            if (identical(z_next, z_curr) || abs(z_next - z_curr) <= .Machine$double.eps) {
                ## The objective kept improving up to an admissible boundary.
                edge2 <- if (direction < 0) min(z_prev, z_curr) else max(z_prev, z_curr)
                edge1 <- if (direction < 0) aa else bb
                bracket <- sort(c(edge1, edge2))
                break
            }

            f_next <- eval_one(z_next)

            if (!is.finite(f_next) || f_next >= f_curr) {
                bracket <- sort(c(z_prev, z_next))
                break
            }

            z_prev <- z_curr
            f_prev <- f_curr
            z_curr <- z_next
            f_curr <- f_next
        }
    }

    if (is.null(bracket) || length(bracket) != 2L ||
        !all(is.finite(bracket)) || bracket[1L] >= bracket[2L]) {
        ## Very defensive fallback: use a modest start-centred interval rather
        ## than the complete gigantic automatic interval.
        radius <- min(max(base_step * 4, (bb - aa) / 100), (bb - aa) / 2)
        bracket <- c(max(aa, z0 - radius), min(bb, z0 + radius))
    }

    ans <- stats::optimize(
        f = ff,
        interval = bracket,
        maximum = FALSE,
        ...
    )

    ans$minimum <- back(ans$minimum)
    ans$parameter_scale <- parameter_scale
    ans$search_interval <- back(bracket)
    ans$start <- p0
    ans
}
