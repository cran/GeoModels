sp2Geo <- function(spobj, spdata = NULL)
{
  # spobj: sp or spacetime object
  # spdata: name of dependent variable;
  # one of the columns in the @data slot must be the dependent variable

  X <- NULL
  Y <- NULL
  coordt <- NULL

  cls <- class(spobj)[1]

  if (!(cls %in% c("SpatialPoints", "SpatialPointsDataFrame", "STFDF", "STIDF")))
    stop("It is not a sp or st object ")

  if (cls == "SpatialPoints") {
    coords <- sp::coordinates(spobj)
    pj <- sp::is.projected(spobj)
    if (is.na(pj)) warning("coordinate reference system is not defined")
  }

  if (cls == "SpatialPointsDataFrame") {
    if (is.null(spdata)) stop("Dependent variable name must be declared")
    if (!(spdata %in% names(spobj@data))) stop("Dependent variable name is not in the data slot")

    coords <- sp::coordinates(spobj)
    X <- spobj@data[, names(spobj@data) != spdata, drop = FALSE]
    Y <- spobj@data[, spdata]
    pj <- sp::is.projected(spobj)
    if (is.na(pj)) warning("coordinate reference system is not defined")
  }

  if (cls == "STFDF") {
    if (is.null(spdata)) stop("Dependent variable name must be declared")
    if (!(spdata %in% names(spobj@data))) stop("Dependent variable name is not in the data slot")

    coords <- sp::coordinates(spobj@sp)
    pj <- sp::is.projected(spobj@sp)
    coordt <- seq_along(spobj@time)
    Y <- spobj@data[, spdata]
    X <- spobj@data[, names(spobj@data) != spdata, drop = FALSE]
    if (is.na(pj)) warning("coordinate reference system is not defined")
  }

  if (cls == "STIDF") {
    if (is.null(spdata)) stop("Dependent variable name must be declared")

    spobj <- as(spobj, "STFDF")

    if (!(spdata %in% names(spobj@data))) stop("Dependent variable name is not in the data slot")

    coords <- sp::coordinates(spobj@sp)
    pj <- sp::is.projected(spobj@sp)
    coordt <- seq_along(spobj@time)
    Y <- spobj@data[, spdata]
    X <- spobj@data[, names(spobj@data) != spdata, drop = FALSE]
    if (is.na(pj)) warning("coordinate reference system is not defined")
  }

  return(list(coords = coords,
              coordt = coordt,
              X = X,
              Y = Y,
              pj = pj,
              projected = pj))
}