## Validation at R/native boundaries.  Exact names and scalar lengths matter:
## unlist() silently drops absent list elements and expands vector elements.
.GeoValidateParameterNames <- function(param, context = "GeoModels") {
  if (!(is.list(param) || (is.numeric(param) && !is.complex(param))) ||
      is.data.frame(param))
    stop(context, ": parameters must be a named list or numeric vector.", call. = FALSE)
  nm <- names(param)
  if (length(param) && (is.null(nm) || anyNA(nm) || any(!nzchar(nm)) ||
                        anyDuplicated(nm)))
    stop(context, ": parameter names must be nonempty, unique and non-missing.",
         call. = FALSE)
  invisible(TRUE)
}

.GeoValidateNamedScalars <- function(param, required, context = "GeoModels") {
  .GeoValidateParameterNames(param, context)
  bad <- required[!vapply(required, function(nm) {
    if (!(nm %in% names(param))) return(FALSE)
    z <- param[[nm]]
    is.numeric(z) && !is.complex(z) && length(z) == 1L && is.finite(z)
  }, logical(1L))]
  if (length(bad))
    stop(context, ": missing or invalid finite scalar parameter(s): ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  setNames(vapply(required, function(nm) as.double(param[[nm]]), numeric(1L)),
           required)
}

.GeoValidateCorrelationParameters <- function(corrmodel, param, d = NULL,
                                                context = "GeoModels") {
  code <- if (is.character(corrmodel)) CkCorrModel(corrmodel) else corrmodel
  if (!is.numeric(code) || length(code) != 1L || !is.finite(code) ||
      code != floor(code) || is.null(CorrelationPar(code)))
    stop(context, ": invalid correlation model.", call. = FALSE)
  if (is.null(d)) d <- if (code %in% c(85L, 89L)) .GeoSTCurrentDimension() else 2L
  .GeoSTValidate(code, param, d = d, complete = TRUE, context = context)
  z <- .GeoValidateNamedScalars(param, CorrelationPar(code), context)
  fail <- function(msg) stop(context, ": ", msg, call. = FALSE)
  scales <- z[grepl("^scale($|_)", names(z))]
  ## Shkarofski admits a zero boundary in one of its scale parameters.
  if (code == 10L) {
    if (any(scales < 0) || all(scales == 0) ||
        (z[["scale_1"]] == 0 && z[["smooth"]] <= 0) ||
        (z[["scale_2"]] == 0 && z[["smooth"]] >= 0))
      fail("invalid Shkarofski scale/smoothness boundary.")
  } else if (any(scales <= 0)) fail("correlation scales must be positive.")
  if (any(z[grepl("^sill($|_)", names(z))] <= 0)) fail("sills must be positive.")
  ## Bivariate nugget parameters are variances, not proportions.
  if (any(z[grepl("^nugget", names(z))] < 0)) fail("nuggets must be nonnegative.")
  if ("pcol" %in% names(z) && abs(z[["pcol"]]) > 1)
    fail("pcol must be in [-1, 1].")
  if (code %in% c(14L, 20L) && z[["smooth"]] <= 0)
    fail("smoothness must be positive for this correlation model.")
  if (code %in% c(48L, 62L, 86L, 95L, 117L, 118L, 119L, 121L, 122L, 128L) &&
      any(z[grepl("^smooth", names(z))] <= 0))
    fail("Matern smoothness parameters must be positive.")
  if (code %in% c(1L, 5L, 8L, 24L, 25L) && z[["power2"]] <= 0)
    fail("power2 must be positive.")
  if (code %in% c(5L, 8L) && (z[["power1"]] <= 0 || z[["power1"]] > 2))
    fail("power1 must be in (0, 2].")
  if (code %in% c(12L, 17L) && (z[["power"]] <= 0 || z[["power"]] > 2))
    fail("power must be in (0, 2].")
  if (code %in% c(11L, 13L, 15L)) {
    minimum <- c(`11` = 1.5, `13` = 2.5, `15` = 3.5)[as.character(code)]
    if (z[["power2"]] < minimum) fail("invalid Wendland power2.")
  }
  if (code %in% c(26L, 29L)) {
    model_name <- if (code == 26L) "GenWend_Hole" else "GenWend_Matern_Hole"
    .GeoValidateHoleWendlandParameters(model_name, as.list(z), d, context)
  }
  unname(z)
}

.GeoValidateCurveNuisance <- function(param, model, bivariate,
                                     context = "GeoCorrFct") {
  required <- NuisParam(model, bivariate, c(1L, 1L))
  ## Bivariate Gaussian lag curves do not depend on the two means.
  if (bivariate) required <- required[!grepl("^mean", required)]
  z <- .GeoValidateNamedScalars(param, required, context)
  fail <- function(msg) stop(context, ": ", msg, call. = FALSE)
  if ("nugget" %in% names(z) && (z[["nugget"]] < 0 || z[["nugget"]] > 1))
    fail("nugget must be in [0, 1].")
  if ("sill" %in% names(z) && z[["sill"]] <= 0) fail("sill must be positive.")
  if (any(z[grepl("^shape", names(z))] <= 0)) fail("shape parameters must be positive.")
  if (model == "LogLogistic" && z[["shape"]] <= 2)
    fail("LogLogistic covariance requires shape > 2.")
  if (model %in% c("StudentT", "SkewStudentT", "TwoPieceStudentT") &&
      (z[["df"]] <= 0 || z[["df"]] >= 0.5))
    fail("finite-variance Student covariance requires inverse df in (0, 0.5).")
  if (model == "SinhAsinh" && z[["tail"]] <= 0)
    fail("SinhAsinh tail must be positive.")
  if (model %in% c("Tukeyh", "Tukeyh2", "Tukeygh", "TwoPieceTukeyh") &&
      any(z[grepl("^tail", names(z))] < 0 | z[grepl("^tail", names(z))] >= 0.5))
    fail("finite-variance Tukey covariance requires tail parameters in [0, 0.5).")
  if (model %in% c("SkewStudentT", "TwoPieceStudentT", "TwoPieceGaussian", "TwoPieceTukeyh") &&
      abs(z[["skew"]]) >= 1)
    fail("this model requires skew in (-1, 1).")
  invisible(TRUE)
}
