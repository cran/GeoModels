GeoVariogram <- function(data, coordx, coordy=NULL, coordz=NULL, coordt=NULL, coordx_dyn=NULL,
                         cloud=FALSE, distance="Eucl",
                         grid=FALSE, maxdist=NULL, neighb=NULL, maxtime=NULL, numbins=NULL,
                         radius=1, type='variogram', bivariate=FALSE, subsample=1, subsample_t=1)
{
  call <- match.call()
  corrmodel <- 'exponential'

    # -------------------------------------------------------------------
  # If coordx is omitted, derive a baseline coordx from coordx_dyn
  # (coordx_dyn is a list of coordinate matrices, one per time)
  # -------------------------------------------------------------------
  if (missing(coordx)) {
    if (is.null(coordx_dyn))
      stop("argument 'coordx' is missing: provide 'coordx' or 'coordx_dyn'.")
    if (!is.list(coordx_dyn) || length(coordx_dyn) < 1L)
      stop("'coordx_dyn' must be a non-empty list when 'coordx' is omitted.")
    base <- coordx_dyn[[1L]]

    if (is.data.frame(base)) base <- as.matrix(base)
    if (!is.matrix(base))
      stop("coordx_dyn[[1]] must be a matrix/data.frame (n x d).")

    # (optional but very useful) sanity check: list length vs coordt
    if (!is.null(coordt) && length(coordt) > 1L && length(coordx_dyn) != length(coordt)) {
      warning(sprintf("length(coordx_dyn)=%d differs from length(coordt)=%d. Check your dynamic coordinates.",
                      length(coordx_dyn), length(coordt)))
    }
    coordx <- base
  }

  if (is.data.frame(coordx)) coordx <- as.matrix(coordx)


  ### Check subsample parameters:
  if(!is.numeric(subsample) || length(subsample) != 1L || subsample <= 0 || subsample > 1)
    stop('subsample must be a numeric value between 0 and 1\n')
  if(!is.numeric(subsample_t) || length(subsample_t) != 1L || subsample_t <= 0 || subsample_t > 1)
    stop('subsample_t must be a numeric value between 0 and 1\n')

  if(is.null(type)) type <- 'variogram'

  # Set the type of model:
  if(type == 'variogram'){
    model <- 'Gaussian'
    fname <- 'Binned_Variogram'
  }

  # Checks if its a spatial or spatial-temporal random field:
  if(bivariate) coordt <- c(0,1)
  if(!is.null(coordt))
    if(is.numeric(coordt))
      if(length(coordt) > 1) corrmodel <- 'gneiting'
  if (is.data.frame(coordx)) coordx <- as.matrix(coordx)
  # -------------------------------------------------------------------
  # EARLY SUBSAMPLING (MUST BE BEFORE CkInput/StartParam WHICH SET C GLOBALS)
  # -------------------------------------------------------------------
  if (bivariate) subsample_t <- 1
  # EARLY temporal subsample (only if it's WELL-DEFINED: matrix cols or list length = length(coordt))
  if (!is.null(coordt) && length(coordt) > 1 && subsample_t < 1) {

    if (is.matrix(data) && ncol(data) == length(coordt)) {
      nt0 <- length(coordt)
      mt  <- max(2L, floor(subsample_t * nt0))
      mt  <- min(mt, nt0)
      idx_t <- sort(sample.int(nt0, mt, replace = FALSE))
      coordt <- coordt[idx_t]
      data   <- data[, idx_t, drop = FALSE]

    } else if (is.list(data) && length(data) == length(coordt)) {
      nt0 <- length(coordt)
      mt  <- max(2L, floor(subsample_t * nt0))
      mt  <- min(mt, nt0)
      idx_t <- sort(sample.int(nt0, mt, replace = FALSE))
      coordt <- coordt[idx_t]
      data   <- data[idx_t]

    } else {
      warning("subsample_t < 1 requested, but data/coordt format is not (ns x nt) matrix nor list-of-times; ignoring subsample_t.")
    }
  }

  # helper: subset data by spatial index idx0, robust to orientation (incl. bivariate 2 x n0)
  .subset_data_spatial <- function(data, idx0, n0) {
    if (is.matrix(data)) {
      if (nrow(data) == n0) return(data[idx0, , drop = FALSE])   # ns x p or ns x nt
      if (ncol(data) == n0) return(data[, idx0, drop = FALSE])   # p x ns (e.g. 2 x ns bivariate)
      stop("Early subsample: data is a matrix but neither nrow nor ncol matches number of locations (n0).")
    }
    if (is.list(data)) {
      if (!all(vapply(data, function(d) length(d) == n0, logical(1))))
        stop("Early subsample: data is a list but not all elements have length n0.")
      return(lapply(data, function(d) d[idx0]))
    }
    if (length(data) != n0)
      stop("Early subsample: data vector length does not match number of locations (n0).")
    data[idx0]
  }

  # EARLY spatial subsample (points/locations). For spacetime wide matrix, subsample rows (locations).
  if (subsample < 1) {

    # Case: coordx is matrix of spatial coordinates (n x d)
    if (is.matrix(coordx)) {
      n0 <- nrow(coordx)
      m  <- max(2L, floor(subsample * n0))
      m  <- min(m, n0)
      idx0 <- sort(sample.int(n0, m, replace = FALSE))

      coordx <- coordx[idx0, , drop = FALSE]

      # If coordy/coordz provided separately (rare when coordx is matrix), subset them too
      if (!is.null(coordy) && length(coordy) == n0) coordy <- coordy[idx0]
      if (!is.null(coordz) && length(coordz) == n0) coordz <- coordz[idx0]

      # data: robust to orientation (incl. bivariate 2 x n0 or n0 x 2)
      data <- .subset_data_spatial(data, idx0, n0)

      # If coordx_dyn is point-aligned (length n0), subset it
      if (!is.null(coordx_dyn) && length(coordx_dyn) == n0) coordx_dyn <- coordx_dyn[idx0]

    } else {
      # Case: coordx is vector, coordy optional vector, coordz optional vector
      n0 <- length(coordx)
      m  <- max(2L, floor(subsample * n0))
      m  <- min(m, n0)
      idx0 <- sort(sample.int(n0, m, replace = FALSE))

      coordx <- coordx[idx0]
      if (!is.null(coordy) && length(coordy) == n0) coordy <- coordy[idx0]
      if (!is.null(coordz) && length(coordz) == n0) coordz <- coordz[idx0]

      data <- .subset_data_spatial(data, idx0, n0)

      if (!is.null(coordx_dyn) && length(coordx_dyn) == n0) coordx_dyn <- coordx_dyn[idx0]
    }
  }

  # After early subsampling, do NOT subsample later (avoids mismatch with C globals)
  #subsample   <- 1
  #subsample_t <- 1

  # -------------------------------------------------------------------
  # Now run input checks + StartParam (sets C globals consistently)
  # -------------------------------------------------------------------
  checkinput <- CkInput(coordx, coordy, coordz, coordt, coordx_dyn, corrmodel, data, distance, "Fitting", NULL, grid,
                        'None', maxdist, maxtime, model, NULL, 'Nelder-Mead', NULL,
                        radius, NULL, NULL, NULL, 'GeoWLS', FALSE, FALSE, NULL, NULL)
  if(!is.null(checkinput$error))
    stop(checkinput$error)

  ### Specific checks of the Empirical Variogram:
  if (!is.logical(cloud) || length(cloud) != 1L)
    stop('insert a logical value (TRUE/FALSE) for the cloud parameter\n')

  if (is.null(numbins)) {
    numbins <- 13L
  } else {
    if (!is.numeric(numbins) || numbins < 2 || numbins %% 1 != 0)
      stop('insert an integer >= 2 for the number of bins\n')
    numbins <- as.integer(numbins)
  }

  if(bivariate) corrmodel <- 'Bi_matern_sep'

  n <- 1

  initparam <- StartParam(coordx, coordy, coordz, coordt, coordx_dyn, corrmodel, data, distance, "Fitting",
                          NULL, grid, 'None', maxdist, neighb, maxtime, model, n,
                          NULL, NULL, FALSE, radius, NULL, NULL, NULL,
                          'GeoWLS', 'GeoWLS', FALSE, NULL, NULL, FALSE, FALSE)

  spacetime_dyn <- NULL
  coordx <- initparam$coordx
  coordy <- initparam$coordy
  coordz <- initparam$coordz
  coordt <- initparam$coordt

  if(!is.null(initparam$error))
    stop(initparam$error)

  numvario <- numbins - 1L
  if(cloud){
    numbins <- numvario <- initparam$numpairs
    fname <- 'Cloud_Variogram'
  }

  # alloc
  bins <- double(numbins)
  moments <- double(numvario)
  lenbins <- integer(numvario)
  bint <- NULL
  lenbinst <- NULL
  lenbint <- NULL
  variogramst <- NULL
  variogramt <- NULL
  centert <- NULL
  regular <- NULL
  centers <- NULL
  variograms <- NULL
  memdist <- FALSE

  #***********************************************************************************************#
  #********************** bivariate ****************************************************************#
  #***********************************************************************************************#
  if(initparam$bivariate){
    memdist <- FALSE
    if(!is.null(neighb)) memdist <- TRUE

    n_var <- initparam$numtime
    spacetime_dyn <- FALSE
    if(!is.null(coordx_dyn)) spacetime_dyn <- TRUE
    ns <- initparam$ns
    NS <- cumsum(ns)

    if(!spacetime_dyn){
      data <- c(t(data))
      if(memdist){
        coordx <- rep(coordx, n_var)
        coordy <- rep(coordy, n_var)
        if(!is.null(coordz)) coordz <- rep(coordz, n_var)
      }
    }
    if(spacetime_dyn) data <- unlist(data)

    if(is.null(coordz)) coordz <- double(length(coordx))
    NS <- c(0,NS)[-(length(ns)+1)]

    if(!memdist){
      moments_marg  <- double(n_var * numvario)
      lenbins_marg  <- integer(n_var * numvario)
      moments_cross <- double(0.5 * n_var * (n_var - 1) * numvario)
      lenbins_cross <- integer(0.5 * n_var * (n_var - 1) * numvario)

      DEV <- dotCall64::.C64("Binned_Variogram_biv2",
                             bins=bins, coordx=coordx, coordy=coordy, coordz=coordz,
                             coordt=coordt, data=data,
                             lenbins_cross=lenbins_cross, moments_cross=moments_cross,
                             numbins=numbins,
                             lenbins_marg=lenbins_marg, moments_marg=moments_marg,
                             ns=ns, NS=NS,
                             SIGNATURE=c("double","double","double","double","double","double","integer","double",
                                         "integer","integer","double","integer","integer"),
                             INTENT=c("rw","r","r","r","r","r","rw","rw","r","rw","rw","r","r"),
                             NAOK=TRUE, PACKAGE="GeoModels", VERBOSE=0)

      bins <- DEV$bins
      lenbins_cross <- DEV$lenbins_cross; moments_cross <- DEV$moments_cross
      lenbins_marg  <- DEV$lenbins_marg;  moments_marg  <- DEV$moments_marg

      m_11 <- moments_marg[1:numvario]
      m_22 <- moments_marg[(numvario+1):(2*numvario)]
      m_12 <- moments_cross[1:numvario]
      l_11 <- lenbins_marg[1:numvario]
      l_22 <- lenbins_marg[(numvario+1):(2*numvario)]
      l_12 <- lenbins_cross[1:numvario]

      indbin_marg  <- l_11 > 0
      indbin_cross <- l_12 > 0

      centers_all <- bins[-length(bins)] + diff(bins)/2

      m_11 <- m_11[indbin_marg];  m_22 <- m_22[indbin_marg]
      l_11 <- l_11[indbin_marg];  l_22 <- l_22[indbin_marg]
      centers <- centers_all[indbin_marg]

      m_12 <- m_12[indbin_cross]; l_12 <- l_12[indbin_cross]

      variograms_11 <- m_11/l_11
      variograms_22 <- m_22/l_22
      variograms_12 <- m_12/l_12

      lenbins     <- rbind(l_11, l_22)
      lenbinst    <- l_12
      variograms  <- rbind(variograms_11, variograms_22)
      variogramst <- variograms_12

    } else {
      idx <- GeoNeighIndex(coordx=cbind(coordx,coordy,coordz), coordx_dyn=coordx_dyn,
                           distance=distance, neighb=neighb, maxdist=maxdist,
                           maxtime=1, radius=1, bivariate=TRUE)
      mm <- range(idx$lags)

      moments00 <- double(numvario); moments10 <- double(numvario); moments11 <- double(numvario)
      lenbins00 <- integer(numvario); lenbins10 <- integer(numvario); lenbins11 <- integer(numvario)

      DEV <- dotCall64::.C64("Binned_Variogram_biv2new",
                             SIGNATURE=c("double","integer","double","double",
                                         "double","double","double","double","double",
                                         "integer","integer","integer",
                                         "integer","integer","integer"),
                             bins=bins, length(idx$lags),
                             data[idx$colidx], data[idx$rowidx],
                             idx$lags, mm,
                             moments00=moments00, moments10=moments10, moments11=moments11,
                             lenbins00=lenbins00, lenbins10=lenbins10, lenbins11=lenbins11,
                             numbins, idx$first, idx$second,
                             INTENT=c("w","r","r","r",
                                      "r","r","w","w","w","w","w","w","r","r","r"),
                             NAOK=TRUE, PACKAGE="GeoModels", VERBOSE=0)

      bins <- DEV$bins
      centers_all <- bins[-length(bins)] + diff(bins)/2

      lenbins   <- rbind(DEV$lenbins00, DEV$lenbins11)
      lenbinst  <- DEV$lenbins10
      variograms <- rbind(DEV$moments00/DEV$lenbins00, DEV$moments11/DEV$lenbins11)
      variogramst <- DEV$moments10/DEV$lenbins10

      indbin <- DEV$lenbins00 > 0
      centers <- centers_all[indbin]
      lenbins <- lenbins[, indbin, drop=FALSE]
      variograms <- variograms[, indbin, drop=FALSE]
    }
  }

  #***********************************************************************************************#
  #********************** spacetime ****************************************************************#
  #***********************************************************************************************#
  if(initparam$spacetime){
    memdist <- FALSE
    if(!is.null(neighb)) memdist <- TRUE

    numtime <- initparam$numtime
    spacetime_dyn <- FALSE
    if(!is.null(coordx_dyn)) spacetime_dyn <- TRUE
    ns <- initparam$ns
    NS <- cumsum(ns)

    numbint <- numtime - 1L
    bint <- double(numbint)
    momentt <- double(numbint)
    lenbint <- integer(numbint)
    numbinst <- numvario * numbint
    momentst <- double(numbinst)
    lenbinst <- integer(numbinst)

    if(grid){
      if(is.null(coordz)) { a <- expand.grid(coordx, coordy); coordx <- a[,1]; coordy <- a[,2] }
      else { a <- expand.grid(coordx, coordy, coordz); coordx <- a[,1]; coordy <- a[,2]; coordz <- a[,3] }
    } else {
      if(!spacetime_dyn) data <- c(t(data))
      else data <- unlist(data)
    }

    NS <- c(0,NS)[-(length(ns)+1)]

    d <- diff(coordt)
    regular <- length(d) > 0 && length(unique(d)) == 1L
    if (regular) {
      u <- d[1]
      bint <- seq(0, max(coordt) - u, by = u)
    } else {
      regular <- FALSE
      step <- max(nabor::knn(matrix(coordt, ncol = 1), k = length(coordt))$nn.dists) / (numbint - 1)
      bint <- numeric(numbint); bint[1] <- 0
      for (h in 2:numbint) bint[h] <- bint[h - 1] + step
    }

    if(!memdist){
      if(!spacetime_dyn) {
        EV <- dotCall64::.C64(
          'Binned_Variogram_st2',
          bins=bins, bint=bint, coordx=coordx, coordy=coordy, coordz=coordz,
          coordt=coordt, data=data,
          lenbins=lenbins, lenbinst=lenbinst, lenbint=lenbint,
          moments=moments, momentst=momentst, momentt=momentt,
          numbins=numbins, numbint=numbint, ns=ns, NS=NS,
          SIGNATURE=c(rep("double", 7), rep("integer", 3), rep("double", 3), rep("integer", 4)),
          INTENT=c("rw","rw", rep("r", 5), rep("rw", 3), rep("rw", 3), rep("r", 4)),
          NAOK=TRUE, PACKAGE='GeoModels', VERBOSE=0
        )
      } else {
        EV <- dotCall64::.C64(
          'Binned_Variogram_st2_dyn',
          bins=bins, bint=bint, coordx=coordx, coordy=coordy, coordz=coordz,
          coordt=coordt, data=data,
          lenbins=lenbins, lenbinst=lenbinst, lenbint=lenbint,
          moments=moments, momentst=momentst, momentt=momentt,
          numbins=numbins, numbint=numbint, ns=ns, NS=NS,
          SIGNATURE=c(rep("double", 7), rep("integer", 3), rep("double", 3), rep("integer", 4)),
          INTENT=c("rw","rw", rep("r", 5), rep("rw", 3), rep("rw", 3), rep("r", 4)),
          NAOK=TRUE, PACKAGE='GeoModels', VERBOSE=0
        )
      }
    } else {
      stop("not implemented")
    }

    bins <- EV$bins; lenbins <- EV$lenbins
    bint <- EV$bint; lenbint <- EV$lenbint
    lenbinst <- EV$lenbinst
    moments <- EV$moments; momentt <- EV$momentt
    momentst <- EV$momentst

    centers_all <- bins[-length(bins)] + diff(bins)/2

    if (regular) centert <- bint
    else centert <- if (length(bint) > 1L) bint[-length(bint)] + diff(bint)/2 else NULL

    indbin   <- lenbins > 0
    indbint  <- lenbint > 0
    indbinst <- lenbinst > 0

    centers <- centers_all[indbin]
    moments <- moments[indbin]
    lenbins <- lenbins[indbin]

    if (!is.null(centert)) {
      bint    <- bint[indbint]
      centert <- centert[indbint]
      momentt <- momentt[indbint]
      lenbint <- lenbint[indbint]
    }

    momentst <- momentst[indbinst]
    lenbinst <- lenbinst[indbinst]

    variograms  <- moments/lenbins
    variogramt  <- if (!is.null(lenbint))  momentt/lenbint else NULL
    variogramst <- if (!is.null(lenbinst)) momentst/lenbinst else NULL
  }

  #***********************************************************************************************#
  #********************** spatial *****************************************************************#
  #***********************************************************************************************#
  if(!initparam$bivariate && !initparam$spacetime){
    memdist <- FALSE
    if(!is.null(neighb)) memdist <- TRUE

    if(grid){
      if(is.null(coordz)) { a <- expand.grid(coordx, coordy); coordx <- a[,1]; coordy <- a[,2] }
      else { a <- expand.grid(coordx, coordy, coordz); coordx <- a[,1]; coordy <- a[,2]; coordz <- a[,3] }
    }

    if(is.null(coordz)) coordz <- double(length(coordx))

    if(!memdist){
      fname <- paste(fname,"2",sep="")
      EV <- .C("Binned_Variogram2",
               bins=bins,
               as.double(coordx), as.double(coordy), as.double(coordz),
               as.double(coordt), as.double(data),
               lenbins=as.integer(lenbins),
               moments=as.double(moments),
               as.integer(numbins),
               PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)

      bins    <- EV$bins
      lenbins <- EV$lenbins
      moments <- EV$moments

      indbin <- lenbins > 0

      if (cloud) {
        bins    <- bins[indbin]
        centers <- bins
      } else {
        centers_all <- bins[-length(bins)] + diff(bins)/2
        centers <- centers_all[indbin]
      }

      moments <- moments[indbin]
      lenbins <- lenbins[indbin]
      variograms <- moments/lenbins

    } else {
      fname <- "Binned_Variogram2new"

      idx <- GeoNeighIndex(cbind(coordx, coordy, coordz),
                           distance=distance, neighb=neighb, maxdist=maxdist, radius=radius)

      mm <- range(idx$lags)

      bins    <- double(numbins)
      lenbins <- integer(numbins - 1L)
      moments <- double(numbins - 1L)

      data1 <- as.double(data[idx$rowidx])
      data2 <- as.double(data[idx$colidx])
      vdist <- as.double(idx$lags)
      np    <- as.integer(length(vdist))

      EV <- dotCall64::.C64(
        "Binned_Variogram2new",
        SIGNATURE=c("double","integer",rep("double",3),"integer","double","integer","double"),
        INTENT=c("rw","r","r","r","r","rw","rw","r","r"),
        bins=bins, np=np, data1=data1, data2=data2, vdist=vdist,
        lbins=lenbins, moms=moments,
        nbins=as.integer(numbins), mm=as.double(mm),
        NAOK=TRUE, PACKAGE="GeoModels", VERBOSE=0
      )

      bins    <- EV$bins
      lenbins <- EV$lbins
      moments <- EV$moms

      centers_all <- bins[-length(bins)] + diff(bins)/2
      indbin  <- lenbins > 0L
      centers <- centers_all[indbin]
      lenbins <- lenbins[indbin]
      moments <- moments[indbin]
      variograms <- moments / lenbins
    }
  }

  # Cleanup globale
  if(!memdist) .C('DeleteGlobalVar', PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
  else .C('DeleteGlobalVar2', PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)

  GeoVariogram <- list(bins=bins,
                       bint=bint,
                       bivariate=bivariate,
                       cloud=cloud,
                       centers=centers,
                       centert=centert,
                       lenbins=lenbins,
                       lenbinst=lenbinst,
                       lenbint=lenbint,
                       maxdist=maxdist,
                       maxtime=maxtime,
                       spacetime_dyn=spacetime_dyn,
                       subsample=subsample,
                       subsample_t=subsample_t,
                       variograms=variograms,
                       variogramst=variogramst,
                       variogramt=variogramt,
                       type=type)

  structure(c(GeoVariogram, call=call), class=c("GeoVariogram"))
}

####################################################
####################################################
####################################################
plot.GeoVariogram <- function(x, ...) {
  if (!inherits(x, "GeoVariogram"))
    stop("Enter an object obtained from the function GeoVariogram\n")

  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  dots <- list(...)

  xlab <- if (!is.null(dots$xlab)) dots$xlab else "Distance"
  ylab <- if (!is.null(dots$ylab)) dots$ylab else "Semi-Variogram"
  dots$xlab <- NULL
  dots$ylab <- NULL

  ispatim <- !is.null(x$bint)
  bivariate <- isTRUE(x$bivariate)

  lags <- c(0, x$centers)
  lagt <- if (ispatim) c(0, x$bint) else 0

  vario.main <- if (ispatim) "Space-time semi-variogram" else "Spatial semi-variogram"
  vario.zlab <- if (ispatim) ylab else NULL

  if (bivariate) {
    par(mfrow = c(2, 2))

    do.call(plot.default, c(list(x = x$centers, y = x$variograms[1,],
                                 main = "First semi-variogram",
                                 ylim = c(0, max(x$variograms[1,])),
                                 xlim = c(0, max(x$centers)),
                                 xlab = xlab, ylab = ylab), dots))

    ll1 <- if (min(x$variogramst) < 0) min(x$variogramst) else 0
    ll <- max(abs(x$variogramst))

    do.call(plot.default, c(list(x = x$centers, y = x$variogramst,
                                 main = "Cross semi-variogram",
                                 ylim = c(ll1, ll),
                                 xlim = c(0, max(x$centers)),
                                 xlab = xlab, ylab = ylab), dots))

    do.call(plot.default, c(list(x = x$centers, y = x$variogramst,
                                 main = "Cross semi-variogram",
                                 ylim = c(ll1, ll),
                                 xlim = c(0, max(x$centers)),
                                 xlab = xlab, ylab = ylab), dots))

    do.call(plot.default, c(list(x = x$centers, y = x$variograms[2,],
                                 main = "Second semi-variogram",
                                 ylim = c(0, max(x$variograms[2,])),
                                 xlim = c(0, max(x$centers)),
                                 xlab = xlab, ylab = ylab), dots))
  }

  if (ispatim) {
    par(mfrow = c(2, 2))

    if (!is.null(x$centert)) {
      evario <- matrix(x$variogramst, nrow = length(x$centers),
                       ncol = length(x$centert), byrow = TRUE)
      evario <- rbind(c(0, x$variogramt), cbind(x$variograms, evario))
      evario.grid <- expand.grid(c(0, x$centers), c(0, x$centert))
    } else {
      evario <- matrix(x$variogramst, nrow = length(x$centers),
                       ncol = length(x$bint), byrow = TRUE)
      evario <- rbind(c(0, x$variogramt), cbind(x$variograms, evario))
      evario.grid <- expand.grid(c(0, x$centers), c(0, x$bint))
    }

    scatterplot3d::scatterplot3d(evario.grid[, 1], evario.grid[, 2], c(evario),
                                 type = "h", highlight.3d = TRUE,
                                 cex.axis = .7, cex.lab = .7,
                                 main = paste("Empirical", vario.main),
                                 xlab = "Distance", ylab = "Time", zlab = vario.zlab,
                                 mar = c(2, 2, 2, 2), mgp = c(0, 0, 0))

    par(mai = c(.2, .2, .2, .2))
    persp(x = unique(evario.grid[, 1]), y = unique(evario.grid[, 2]), z = evario,
          xlab = "h", ylab = "u", zlab = expression(gamma(h, u)),
          ltheta = 90, shade = 0.75, ticktype = "detailed", phi = 30,
          theta = 30, main = "Smoothed space-time semi-variogram",
          cex.axis = .8, cex.lab = .8)

    par(mai = c(.5, .5, .5, .5), mgp = c(1.6, .6, 0))
    if (!is.null(x$centert)) {
      do.call(plot.default, c(list(x = x$centert, y = x$variogramt,
                                   xlab = "t", ylab = expression(gamma(t)),
                                   ylim = c(0, max(x$variogramt)),
                                   xlim = c(0, max(x$bint)),
                                   main = "Marginal temporal semi-variogram"), dots))
    } else {
      do.call(plot.default, c(list(x = x$bint, y = x$variogramt,
                                   xlab = "t", ylab = expression(gamma(t)),
                                   ylim = c(0, max(x$variogramt)),
                                   xlim = c(0, max(x$bint)),
                                   main = "Marginal temporal semi-variogram"), dots))
    }

    do.call(plot.default, c(list(x = x$centers, y = x$variograms,
                                 xlab = "h", ylab = expression(gamma(h)),
                                 ylim = c(0, max(x$variograms)),
                                 xlim = c(0, max(x$centers)),
                                 main = "Marginal spatial semi-variogram"), dots))
  }

  if (!ispatim && !bivariate) {
    do.call(plot.default, c(list(x = x$centers, y = x$variograms,
                                 xlab = xlab, ylab = ylab), dots))
  }

  invisible()
}
