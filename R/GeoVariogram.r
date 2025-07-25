####################################################
### File name: GeoVariogram.r
####################################################

### Procedures are in alphabetical order.

GeoVariogram <- function(data, coordx, coordy=NULL, coordz=NULL,coordt=NULL, coordx_dyn=NULL,cloud=FALSE, distance="Eucl",
                       grid=FALSE, maxdist=NULL,neighb=NULL, maxtime=NULL,numbins=NULL,
                       radius=1, type='variogram',bivariate=FALSE)
  {
    call <- match.call()
    corrmodel <- 'exponential'
    ### Check the parameters given in input:
    if(is.null(type))
      type <- 'variogram'
    # Set the type of model:
    if(type=='variogram'){
        model <- 'Gaussian'
        fname <- 'Binned_Variogram'}

    
    # Checks if its a spatial or spatial-temporal random field:
    if(bivariate) coordt=c(0,1)
    if(!is.null(coordt))
      if(is.numeric(coordt))
        if(length(coordt)>1) corrmodel <- 'gneiting'
    # Checks the input:

    checkinput <- CkInput(coordx, coordy, coordz,coordt, coordx_dyn, corrmodel, data, distance, "Fitting", NULL, grid,
                             'None', maxdist, maxtime, model,NULL, 'Nelder-Mead', NULL,
                             radius,  NULL, NULL,NULL, 'GeoWLS', FALSE, FALSE,NULL,NULL)
                             

    # Checks if there are errors in the input:
    if(!is.null(checkinput$error))
      stop(checkinput$error)
    ### START -- Specific checks of the Empirical Variogram:
    if(!is.null(cloud) & !is.logical(cloud))
      stop('insert a logical value (TRUE/FALSE) for the cloud parameter\n')
    
    if(!is.null(numbins) & !is.integer(numbins))
      if(numbins < 0)
        stop('insert a positive integer value for the number of bins\n')
    if(is.null(numbins))
      numbins <- 13
    ### END -- Specific checks of the Empirical Variogram
    if(bivariate)  corrmodel <- 'Bi_matern_sep'
  

    n=1

    initparam <- StartParam(coordx, coordy, coordz,coordt,coordx_dyn, corrmodel, data,distance, "Fitting",
                           NULL, grid, 'None', maxdist,neighb, maxtime, model, n, 
                           NULL, NULL, FALSE, radius, NULL, NULL, NULL,
                            'GeoWLS', 'GeoWLS', FALSE,NULL,NULL,FALSE,FALSE)



    spacetime_dyn=NULL
    coordx=initparam$coordx;coordy=initparam$coordy;coordz=initparam$coordz;
    coordt=initparam$coordt  
              
    # Checks if there are inconsistences:
    if(!is.null(initparam$error))
      stop(initparam$error)
    numvario <- numbins-1
    if(cloud){
        numbins <- numvario <- initparam$numpairs
        fname <- 'Cloud_Variogram'}
    ### Estimation of the empirical spatial or spatial-temporal variogram:
    bins <- double(numbins) # spatial bins
    moments <- double(numvario) # vector of spatial moments
    lenbins <- integer(numvario) # vector of spatial bin sizes
    bint <- NULL
    lenbinst <- NULL
    lenbint <- NULL
    variogramst <- NULL
    variogramt <- NULL  
    centert=NULL
    regular=NULL
  #***********************************************************************************************#
  #********************** bivariate *************************************************************************#
  #***********************************************************************************************#
    if(initparam$bivariate){
  memdist=FALSE
  if(!is.null(neighb)) memdist=TRUE

    n_var=initparam$numtime
    spacetime_dyn=FALSE
    if(!is.null(coordx_dyn)) spacetime_dyn=TRUE
    ns=initparam$ns
    NS=cumsum(ns)
    if(!spacetime_dyn){  data=c(t(data))
                         if(memdist){coordx=rep(coordx,n_var)
                                     coordy=rep(coordy,n_var)
                                     if(!is.null(coordz))  coordz=rep(coordz,n_var)
                                   }
                      }
    if(spacetime_dyn) {data=unlist(data)}
    

    if(is.null(coordz)) coordz=double(length(coordx))
    NS=c(0,NS)[-(length(ns)+1)]
    if(!memdist)  {  
               moments_marg<-double(n_var*numvario)   # vect of square differences for each component (n_var) 11  e 22
               lenbins_marg<-integer(n_var*numvario)  #
               moments_cross<-double(0.5*n_var*(n_var-1)*numvario)  # vect of square differences for cross components (12)
               lenbins_cross<-integer(0.5*n_var*(n_var-1)*numvario) #

         

  DEV=dotCall64::.C64("Binned_Variogram_biv2",bins=bins,coordx=coordx,coordy=coordy,coordz=coordz,
    coordt=coordt,data=data,lenbins_cross=lenbins_cross,moments_cross=moments_cross,
    numbins=numbins,lenbins_marg=lenbins_marg,moments_marg=moments_marg,ns=ns,NS=NS,
    SIGNATURE=c("double","double","double","double","double","double","integer","double",
      "integer","integer","double","integer","integer"),
    INTENT=c("rw","r","r","r","r","r","rw","rw","r","rw","rw","r","r"),
    NAOK=TRUE,PACKAGE="GeoModels",VERBOSE=0)




               bins=DEV$bins
               lenbins_cross=DEV$lenbins_cross; moments_cross=DEV$moments_cross
               lenbins_marg=DEV$lenbins_marg; moments_marg=DEV$moments_marg
               m_11=moments_marg[1:numvario];m_22=moments_marg[(numvario+1):(2*numvario)];m_12=moments_cross[1:numvario];
               l_11=lenbins_marg[1:numvario];l_22=lenbins_marg[(numvario+1):(2*numvario)];l_12=lenbins_cross[1:numvario];
               indbin_marg <- l_11>0;indbin_cross <- l_12>0  #
               bins<- bins[indbin_marg]; numbins <-length(bins)
               m_11 <- m_11[indbin_marg];m_22 <- m_22[indbin_marg];l_11 <- l_11[indbin_marg];l_22 <- l_22[indbin_marg]
               m_12 <- m_12[indbin_cross];l_12 <- l_12[indbin_cross]
               variograms_11 <- m_11/l_11;variograms_22 <- m_22/l_22
               variograms_12 <- m_12/l_12   
               centers <-   bins[1:(numbins[1]-1)]+diff(bins)/2
               lenbins=rbind(l_11,l_22);lenbinst=l_12
               variograms=rbind(variograms_11,variograms_22) 
               variogramst=variograms_12
            }   
             else {
             idx=GeoNeighIndex(coordx=cbind(coordx,coordy,coordz),coordx_dyn=coordx_dyn,
                  distance = distance, neighb = neighb, maxdist = maxdist,maxtime=1,radius=1,bivariate=TRUE)
             mm=range(idx$lags)
             moments00=double(numvario);moments10=double(numvario);moments11=double(numvario)
             lenbins00=double(numvario);lenbins10=double(numvario);lenbins11=double(numvario);

DEV=dotCall64::.C64("Binned_Variogram_biv2new", 
                    SIGNATURE = c("double","integer","double","double", 
                      "double","double","double","double","double",
                                        "integer","integer","integer",
                                        "integer","integer","integer"),
                    bins=bins, length(idx$lags),data[idx$colidx],data[idx$rowidx],
                    idx$lags, mm,  moments00=moments00, moments10=moments10, moments11=moments11, 
                                   lenbins00=lenbins00,lenbins10=lenbins10,lenbins11=lenbins11,
                    numbins,idx$first,idx$second,
                    INTENT = c("w","r","r","r",
                        "r","r","w","w","w","w","w","w","r","r","r"), 
                    NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0)
     #DEV$lenbins10=DEV$lenbins10*2
     bins=DEV$bins
     centers = bins[1:(numbins[1]-1)]+diff(bins)/2
     lenbins=rbind(DEV$lenbins00,DEV$lenbins11)
     lenbinst=DEV$lenbins10
     variograms=rbind(DEV$moments00/DEV$lenbins00,DEV$moments11/DEV$lenbins11) 
     variogramst=DEV$moments10/DEV$lenbins10
                  }  
    }
  #***********************************************************************************************#
  #********************** spacetime *************************************************************************#
  #***********************************************************************************************#
if(initparam$spacetime){
     memdist=FALSE
   if(!is.null(neighb)) memdist=TRUE
      numtime=initparam$numtime
      spacetime_dyn=FALSE
      if(!is.null(coordx_dyn)) spacetime_dyn=TRUE
      ns=initparam$ns
      NS=cumsum(ns)
      numbint <- initparam$numtime-1 # number of temporal bins
      bint <- double(numbint)        # temporal bins
      momentt <- double(numbint)     # vector of temporal moments
      lenbint <- integer(numbint)    # vector of temporal bin sizes
      numbinst <- numvario*numbint   # number of spatial-temporal bins
      binst <- double(numbinst)      # spatial-temporal bins
      momentst <- double(numbinst)   # vector of spatial-temporal moments
      lenbinst <- integer(numbinst)  # vector of spatial-temporal bin sizes
      #if(cloud) fname <- 'Cloud_Variogram_st' else 
    

      if(grid)     {
          if(is.null(coordz)) {a=expand.grid(coordx,coordy);coordx=a[,1];coordy=a[,2]}
          else                {a=expand.grid(coordx,coordy,coordz);coordx=a[,1];coordy=a[,2];coordz=a[,3]}
           }
      else{
              if(!spacetime_dyn) data=c(t(data))
              else               data=unlist(data)
      }
   

       NS=c(0,NS)[-(length(ns)+1)]
    
##### temporal bint #############
    if(var(diff(coordt))==0)           ### regular times   
         { regular=TRUE
           #bint=c(0,coordt[-numtime]);bint=unique(bint)
           u=unique(diff(coordt)); bint=seq(0,max(coordt)-u,u)
           
         }
     else  {     ### not regular times
           regular=FALSE
           step=max(nabor::knn(coordt,  k=length(coordt))$nn.dists)/(numbint-1);
           bint[1]=0;
           for(h in 2:(numbint)) bint[h]=bint[h-1]+step
         }
#########################
 if(is.null(coordz))  coordz=double(length(coordx)*(length(coordx)-1))
  if(!memdist)  { 


    
if(!spacetime_dyn) 
  EV = dotCall64::.C64(
    'Binned_Variogram_st2',
    bins    = bins,bint    = bint,coordx  = coordx,
    coordy  = coordy,coordz  = coordz,
    coordt  = coordt,data    = data,
    lenbins = lenbins,lenbinst= lenbinst,
    lenbint = lenbint,moments = moments,
    momentst= momentst,momentt = momentt,
    numbins = numbins,numbint = numbint,ns= ns, NS= NS,
    SIGNATURE = c(rep("double", 7),rep("integer", 3),rep("double", 3),rep("integer", 4)),
    INTENT = c("rw","rw",rep("r", 5), rep("rw", 3),rep("rw", 3), rep("r", 4)   ),
    NAOK = TRUE, PACKAGE = 'GeoModels', VERBOSE = 0)
if(spacetime_dyn) 
  EV = dotCall64::.C64(
    'Binned_Variogram_st2_dyn',
    bins    = bins,bint    = bint,coordx  = coordx,
    coordy  = coordy,coordz  = coordz,coordt  = coordt,data    = data,lenbins = lenbins,
    lenbinst= lenbinst,lenbint = lenbint,moments = moments,momentst= momentst,
    momentt = momentt,numbins = numbins,numbint = numbint,ns      = ns,NS      = NS,
    SIGNATURE = c(rep("double", 7),rep("integer", 3),rep("double", 3),rep("integer", 4)),
    INTENT = c("rw", "rw", rep("r", 5),rep("rw", 3), rep("rw", 3),rep("r", 4)),
    NAOK = TRUE, PACKAGE = 'GeoModels', VERBOSE = 0
  )

                     
    }
    else {
             stop("not implemented")
             idx=GeoNeighIndex(cbind(coordx,coordy,coordz),coordt=coordt,coordx_dyn=coordx_dyn,distance = distance, neighb = neighb, maxdist = maxdist,maxtime=maxtime,radius=radius)
        }
 

       bins=EV$bins;   lenbins=EV$lenbins
       bint=EV$bint;  lenbint=EV$lenbint
 
       #print(bint)
       lenbinst=EV$lenbinst
       moments=EV$moments; momentt=EV$momentt
       #print(momentt)
       momentst=EV$momentst
       centers <- bins[1:numvario]+diff(bins)/2
       if(!regular) centert <- bint[1:(numbint-1)]+diff(bint)/2
       ###  make sure that lengths in bins/t >0
      indbin <- lenbins>0; indbint <- lenbint>0; indbinst <- lenbinst>0
      bins <- bins[indbin]; bint <- bint[indbint]
      centers <- centers[indbin];    centert <- centert[indbint]; 
      moments <- moments[indbin]; lenbins <- lenbins[indbin]
      momentt <- momentt[indbint]; lenbint <- lenbint[indbint]; 
      momentst <- momentst[indbinst]; lenbinst <- lenbinst[indbinst]
      ####
      variograms <- moments/lenbins; variogramt <- momentt/lenbint; variogramst <- momentst/lenbinst
    }
#***********************************************************************************************#
#********************** spatial *************************************************************************#
#***********************************************************************************************#
 
    if(!initparam$bivariate&&!initparam$spacetime){  
      memdist=FALSE
  if(!is.null(neighb)) memdist=TRUE
     
    if(grid)     {
                      if(is.null(coordz)) {a=expand.grid(coordx,coordy);coordx=a[,1];coordy=a[,2];}
                      else {a=expand.grid(coordx,coordy,coordz);coordx=a[,1];coordy=a[,2];coordz=a[,3];}

                  }

 if(is.null(coordz)) coordz=double(length(coordx))


    #print(head(coordx))
    #  print(head(coordy))
    #   print(head(coordz))
    #   print(memdist)

     if(!memdist)  { 
     fname <- paste(fname,"2",sep="") 
     # Computes the spatial moments
  
      EV=.C("Binned_Variogram2", bins=bins,  as.double(coordx),as.double(coordy),as.double(coordz),as.double(coordt),as.double(data), 
        lenbins=as.integer(lenbins), moments=as.double(moments), as.integer(numbins),PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
       }
      else {
       fname="Binned_Variogram2new"
         idx=GeoNeighIndex(cbind(coordx,coordy,coordz),distance = distance, neighb = neighb, maxdist = maxdist,radius=radius)
         mm=range(idx$lags)
         EV=dotCall64::.C64("Binned_Variogram2",bins=bins,coordx=coordx,coordy=coordy,coordz=coordz,
                   coordt=coordt,data=data,lenbins=lenbins,moments=moments,numbins=numbins,
                   SIGNATURE=c("double","double","double","double","double","double","integer","double","integer"),
                   INTENT=c("rw","r","r","r","r","r","rw","rw","r"),NAOK=TRUE,PACKAGE="GeoModels",VERBOSE=0)


       }
       bins=EV$bins
       lenbins=EV$lenbins
       moments=EV$moments
      # Computes the spatial variogram:
      indbin <- lenbins>0
      bins <- bins[indbin]
      numbins <- length(bins)
      # check if cloud or binned variogram:
      if(cloud) centers <- bins else centers <- bins[1:(numbins-1)]+diff(bins)/2
      moments <- moments[indbin]
      lenbins <- lenbins[indbin]
      variograms <- moments/lenbins}
    # Start --- compute the extremal coefficient
    if(!memdist).C('DeleteGlobalVar', PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
    else .C('DeleteGlobalVar2', PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
 
    GeoVariogram <- list(bins=bins,
                       bint=bint,
                       bivariate=bivariate,
                       cloud=cloud,
                       centers=centers,
                       centert= centert,
                       lenbins=lenbins,
                       lenbinst=lenbinst,
                       lenbint=lenbint,
                       maxdist =maxdist,
                       maxtime = maxtime,
                       spacetime_dyn=spacetime_dyn,
                       variograms=variograms,
                       variogramst=variogramst,
                       variogramt=variogramt,
                       type=type)

    structure(c(GeoVariogram, call = call), class = c("GeoVariogram"))

  }










plot.GeoVariogram <- function(x, ...) {
  if (!inherits(x, "GeoVariogram"))
    stop("Enter an object obtained from the function GeoVariogram\n")

  # Salva impostazioni grafiche
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

  ########## bivariate
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

  ########## spacetime
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

  ########## spatial
  if (!ispatim && !bivariate) {
    do.call(plot.default, c(list(x = x$centers, y = x$variograms,
                                 xlab = xlab, ylab = ylab), dots))
  }

  invisible()
}
