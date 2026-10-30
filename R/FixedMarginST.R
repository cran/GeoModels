## Dimension and parameter handling for the fixed-margin models 85 and 89.
## No additional covariance parameter is introduced: sep is rho in [-1, 1].
.GeoSTCode <- function(corrmodel) {
  code <- if (is.character(corrmodel)) CkCorrModel(corrmodel) else corrmodel
  if (length(code) != 1L || is.na(code) || !(code %in% c(85L, 89L)))
    return(NA_integer_)
  as.integer(code)
}

.GeoSTCurrentDimension <- function() {
  .Call("GeoFixedMarginSTGetDimension", PACKAGE = "GeoModels")
}

.GeoSTDimension <- function(coordx = NULL, coordy = NULL, coordz = NULL,
                            coordx_dyn = NULL, distance = "Eucl",
                            context = "GeoModels") {
  metric <- if (is.numeric(distance)) distance else CheckDistance(distance)
  if (length(metric) != 1L || is.na(metric) || !(metric %in% 0:2))
    stop(context, ": invalid distance.", call. = FALSE)
  if (metric == 2L)
    stop(context, ": these fixed-margin models require Euclidean or chordal distances; ",
         "validity for geodesic distances is not established.", call. = FALSE)
  if (metric == 1L) return(3L)

  ## No public dimension argument or extra covariance parameter is needed.
  ## Coordinates determine the dimension. Lag-only calls inherit an enclosing
  ## coordinate-based evaluation; standalone Euclidean curves use the native
  ## default of two, restored after every scoped evaluation (also on errors).
  d <- .GeoCoordinateDimension(coordx, coordy, coordz, coordx_dyn)
  if (length(d) == 1L && is.na(d)) d <- .GeoSTCurrentDimension()
  if (!is.numeric(d) || length(d) != 1L || !is.finite(d) ||
      d != floor(d) || !(d %in% 1:3))
    stop(context, ": cannot determine a supported spatial dimension from the coordinates.",
         call. = FALSE)
  as.integer(d)
}

.GeoSTEnter <- function(corrmodel, coordx = NULL, coordy = NULL, coordz = NULL,
                        coordx_dyn = NULL, distance = "Eucl",
                        context = "GeoModels") {
  if (is.na(.GeoSTCode(corrmodel))) return(NULL)
  d <- .GeoSTDimension(coordx, coordy, coordz, coordx_dyn, distance,
                       context = context)
  .Call("GeoFixedMarginSTSetDimension", d, PACKAGE = "GeoModels")
}

.GeoSTRestore <- function(old) {
  if (!is.null(old))
    invisible(.Call("GeoFixedMarginSTSetDimension", as.integer(old),
                    PACKAGE = "GeoModels"))
  invisible(NULL)
}

.GeoSTValidate <- function(corrmodel, param, d = NULL, complete = TRUE,
                           context = "GeoModels") {
  code <- .GeoSTCode(corrmodel)
  if (is.na(code)) return(invisible(TRUE))
  if (is.null(d)) d <- .GeoSTCurrentDimension()
  wanted <- CorrelationPar(code)
  if (complete) {
    z <- .GeoValidateNamedScalars(param, wanted, context)
  } else {
    if (is.null(param) || !length(param)) return(invisible(TRUE))
    .GeoValidateParameterNames(param, context)
    z <- .GeoValidateNamedScalars(param, intersect(wanted, names(param)), context)
  }
  fail <- function(message) stop(context, ": ", message, call. = FALSE)
  has <- function(nm) nm %in% names(z)
  for (nm in intersect(c("scale_s", "scale_t"), names(z)))
    if (z[[nm]] <= 0) fail(paste0(nm, " must be positive."))
  if (has("sep") && abs(z[["sep"]]) > 1)
    fail("sep is the normalized interaction and must lie in [-1, 1].")
  for (axis in c("s", "t")) {
    sn <- paste0("smooth_", axis)
    pn <- paste0("power2_", axis)
    if (has(sn)) {
      bound <- if (code == 85L) -0.5 else 0
      if (z[[sn]] <= bound)
        fail(paste0(sn, " must be strictly greater than ", bound, "."))
    }
    if (code == 85L && has(pn)) {
      if (z[[pn]] <= 0) fail(paste0(pn, " must be positive."))
      if (has(sn)) {
        dj <- if (axis == "s") d else 1L
        minimum <- (dj + 1)/2 + z[[sn]] + 2
        if (z[[pn]] <= minimum)
          fail(paste0(pn, " must exceed (d_", axis, "+1)/2 + ", sn,
                      " + 2 = ", format(minimum), "; d_", axis, " = ", dj, "."))
      }
    }
  }
  invisible(TRUE)
}

.GeoSTRepairStarts <- function(initparam, corrmodel, start_names) {
  code <- .GeoSTCode(corrmodel)
  if (is.na(code)) return(initparam)
  par <- initparam$param
  fixed <- initparam$fixed
  auto <- setdiff(names(par), start_names)
  if ("sep" %in% auto) par[["sep"]] <- 0
  if (code == 85L) {
    d <- .GeoSTCurrentDimension()
    for (axis in c("s", "t")) {
      sn <- paste0("smooth_", axis)
      pn <- paste0("power2_", axis)
      base <- ((if (axis == "s") d else 1L) + 1)/2 + 2
      all <- c(par, fixed)
      ## An explicit mu constrains an automatic kappa; never modify user starts.
      if (sn %in% auto && !(pn %in% auto))
        par[[sn]] <- min(0.5, all[[pn]] - base - 0.25)
      all <- c(par, fixed)
      if (pn %in% auto) par[[pn]] <- base + all[[sn]] + 1
    }
  }
  initparam$param <- par
  .GeoSTValidate(code, c(par, fixed), complete = TRUE,
                  context = "GeoFit starting parameters")
  initparam
}
