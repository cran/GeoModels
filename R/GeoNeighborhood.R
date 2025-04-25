GeoNeighborhood = function(data=NULL, coordx=NULL, coordy=NULL,coordz=NULL, coordt=NULL, coordx_dyn=NULL,bivariate=FALSE, 
            distance="Eucl", grid=FALSE, loc, neighb=NULL,maxdist=NULL,maxtime=NULL,
                 radius=6371, time=NULL, X=NULL,M=NULL,spobj=NULL,spdata=NULL,parallel=FALSE,ncores=NULL)
{
  ################################################################################## 
  ############## internal function   ##############################################
  ################################################################################## 
  nbor=function(coords,loc,distance,maxdist,neighb){
  if(distance=="Geod"||distance=="Chor")
  {
    coords_p=coords; loc_p=loc
    prj=mapproj::mapproject(coords_p[,1], coords_p[,2], projection="sinusoidal") 
    coords_p=radius*cbind(prj$x,prj$y)
    prjloc=mapproj::mapproject(loc_p[,1], loc_p[,2], projection="sinusoidal")
    loc_p=radius*cbind(prjloc$x,prjloc$y)
  }
  ##################################################################################
  searchtype = "standard"
  if(!is.null(maxdist)){searchtype = "radius"} else {maxdist=0}
  if(is.null(neighb)) neighb=min(100, NN)
  
  ## computing neigh indexes
  if(distance=="Geod"||distance=="Chor")  {#out<- RANN::nn2(coords_p,loc_p, k=neighb,searchtype=searchtype,radius=maxdist)
                                           out<- nabor::knn(coords_p,loc_p, k=neighb,radius=maxdist)
                                           }
  if(distance=="Eucl")                    {# out<- RANN::nn2(coords,loc,     k=neighb,searchtype=searchtype,radius=maxdist)
                                           out<- nabor::knn(coords,loc, k=neighb,radius=maxdist)
                                          }
    return(out)                                       
  }
  ################################################################################## 
  ################################################################################## 
  ##################################################################################    
  XX=NULL
  MM=NULL
  numtime=1
  sel_ss=1
  sel_tt=1
  coords=NULL

  if(is.vector(loc))    loc=t(as.matrix(loc)) ## case of 1 location sites given as vector
  if(!is.matrix(loc))   stop("loc parameter must be a matrix")
  if(!is.logical(bivariate))   stop("bivariate must be logical")
  #if(!(ncol(loc)==2))   stop("loc parameter must be a matrix  N X 2")
  if(is.null(neighb)&&is.null(maxdist))     stop("maxdist (maxtime) and/or neighb must  be specified")
  if(!is.null(neighb)) neighb=round(neighb)

  spacetime=FALSE

  if(!is.null(coordt)) {spacetime=TRUE}
  if(spacetime) if(!is.vector(time))  stop("time parameter is missing")

  space=!spacetime&&!bivariate 
##################
  if(!is.null(spobj)) {
   if(space||bivariate){
        a=sp2Geo(spobj,spdata); coordx=a$coords 
       if(!a$pj) {if(distance!="Chor") distance="Geod"}
    }
   if(spacetime){
        a=sp2Geo(spobj,spdata); coordx=a$coords ; coordt=a$coordt 
        if(!a$pj) {if(distance!="Chor") distance="Geod"}
     }
   if(!is.null(a$Y)&&!is.null(a$X)) {data=a$Y ; X=a$X }
}
################ 
  dyn=FALSE
  if(!is.null(coordx_dyn))  dyn=TRUE  



  ## handling spatial coordinates

  if(!is.null(coordx)){
  if(is.null(coordy)&&is.null(coordz)) {coords=as.matrix(coordx)}
     else{
          if(grid) {coords=as.matrix(expand.grid(coordx,coordy,coordz))}
          else    { coords=as.matrix(cbind(coordx,coordy,coordz))  }
  }}
  



  Nloc=nrow(loc) # number of location sites
  NN=nrow(coords)
  #####################################
  sel_tt=NULL
  colnames(loc)=NULL;colnames(coords)=NULL;
  ##################################################################
  ##################################################################

  if(space){

    sel_ss=data_sel=numpoints=XX=MM=list()
    out= nbor(coords,loc,distance,maxdist,neighb)

   
    for(i in 1:Nloc)
    {
      ss=out$nn.idx[i,];
      sel_ss[[i]]=coords[ss,]
      numpoints[[i]]=nrow(sel_ss[[i]])
      if(!is.null(data)) data_sel[[i]]=data[ss]
      if(!is.null(M)) MM[[i]]=M[ss]
      if(!is.null(X)) XX[[i]]=X[ss,]
    }
  }
#####################################################################################
if (spacetime) {

  if (!dyn) {  # caso non dinamico

    Tloc <- length(time)
    Nloc <- nrow(loc)
    TT <- length(coordt)
    NN <- nrow(coords)

    if (!is.null(M)) {
      Mma <- matrix(M, nrow = TT, ncol = NN, byrow = TRUE)
    }
    if (!is.null(X)) {
      Xxa <- vector("list", ncol(X))
      for (i in seq_len(ncol(X))) {
        Xxa[[i]] <- matrix(X[, i], nrow = TT, ncol = NN, byrow = TRUE)
      }
    }

    out <- nbor(coords, loc, distance, maxdist, neighb)

    out_t <- nabor::knn(coordt, time, k = length(coordt), radius = maxtime)
    qq <- out_t$nn.idx[1, ]
    ss <- which(qq != 0)
    out_t <- list(ind = cbind(as.vector(out_t$nn.idx[, ss]), time))

    total_points <- Nloc * Tloc
    sel_ss <- vector("list", total_points)
    numpoints <- integer(total_points)
    data_sel <- vector("list", total_points)
    sel_tt <- vector("list", total_points)
    XX <- vector("list", total_points)
    MM <- vector("list", total_points)

    k <- 1
    for (i in seq_len(Nloc)) {
      for (j in seq_len(Tloc)) {

        sel_s <- out$nn.idx[i, ]
        sel_ss[[k]] <- matrix(coords[sel_s, ], ncol = 2)
        numpoints[k] <- nrow(sel_ss[[k]])

        sel_t <- out_t$ind[, 1][out_t$ind[, 2] == time[j]]
        if (length(sel_t) == 0) {
          sel_tt[[k]] <- numeric(0)
          data_sel[[k]] <- NULL
          MM[[k]] <- NULL
          XX[[k]] <- NULL
          k <- k + 1
          next
        }

        pp <- order(coordt[sel_t])
        sel_tt[[k]] <- coordt[sel_t][pp]

        if (!is.null(data)) {
          data_sel[[k]] <- (data[sel_t, sel_s])[pp, , drop = FALSE]
        }
        if (!is.null(M)) {
          QQ <- Mma[sel_t, sel_s, drop = FALSE]
          MM[[k]] <- as.numeric(t(QQ[pp, , drop = FALSE]))
        }
        if (!is.null(X)) {
          KK <- NULL
          for (l in seq_len(ncol(X))) {
            AA <- Xxa[[l]]
            QQ <- AA[sel_t, sel_s, drop = FALSE]
            KK <- cbind(KK, as.numeric(t(QQ[pp, , drop = FALSE])))
          }
          XX[[k]] <- KK
        }

        k <- k + 1
      }
    }
  }  # fine non dinamico

  if (dyn) {  # caso dinamico

    if (length(coordx_dyn) != length(data))
      stop("Length of coordx_dyn and data must match in dynamic spacetime case.")

    Tloc <- length(time)
    Nloc <- nrow(loc)

    total_points <- Nloc * Tloc
    sel_ss <- vector("list", total_points)
    sel_tt <- vector("list", total_points)
    data_sel <- vector("list", total_points)
    XX <- vector("list", total_points)
    MM <- vector("list", total_points)
    numpoints <- integer(total_points)

    out_temp <- nabor::knn(coordt, time, k = length(coordt), radius = maxtime)
    sel_t_all <- out_temp$nn.idx[1, ]
    sel_t_all <- sel_t_all[sel_t_all > 0]

    coordt_sel <- coordt[sel_t_all]
    ord_t <- order(coordt_sel)
    sel_t_ord <- sel_t_all[ord_t]
    coordt_sel_ord <- coordt_sel[ord_t]

    k <- 1
    for (i in seq_len(Nloc)) {
      for (j in seq_len(Tloc)) {

        coords_list <- vector("list", length(sel_t_ord))
        data_list <- vector("list", length(sel_t_ord))
        MM_list <- vector("list", length(sel_t_ord))
        XX_list <- vector("list", length(sel_t_ord))

        for (idx in seq_along(sel_t_ord)) {
          tt <- sel_t_ord[idx]
          coords_t <- coordx_dyn[[tt]]
          data_t <- data[[tt]]
          M_t <- if (!is.null(M)) M[[tt]] else NULL
          X_t <- if (!is.null(X)) X[[tt]] else NULL

          out_spat <- nbor(coords_t, matrix(loc[i, ], ncol = ncol(coords_t)), distance, maxdist, neighb)
          ss <- out_spat$nn.idx[1, ]

          coords_list[[idx]] <- coords_t[ss, , drop = FALSE]
          data_list[[idx]] <- data_t[ss]
          if (!is.null(M)) MM_list[[idx]] <- M_t[ss]
          if (!is.null(X)) XX_list[[idx]] <- X_t[ss, , drop = FALSE]
        }

        sel_ss[[k]] <- coords_list[[1]]
        sel_tt[[k]] <- coordt_sel_ord
        numpoints[k] <- nrow(coords_list[[1]])

        if (!is.null(data)) data_sel[[k]] <- do.call(rbind, data_list)
        if (!is.null(M)) MM[[k]] <- unlist(MM_list)
        if (!is.null(X)) XX[[k]] <- do.call(rbind, XX_list)

        k <- k + 1
      }
    }
  }  # fine dinamico

}  # fine spacetime


#####################################################################################
  if(bivariate)   
  {
    if(dyn){
    Nloc=nrow(loc)
    if(dyn) coords=rbind(coords[[1]],coords[[2]])
    sel_ss=numpoints=data_sel=sel_tt=XX=MM=list()
    if(dyn)      {out1= nbor(coords[[1]],loc,distance,maxdist,neighb)
                 out2= nbor(coords[[2]],loc,distance,maxdist,neighb)}
    else  {out= nbor(coords,loc,distance,maxdist,neighb)}
     
    for(i in 1:Nloc){
      sel=out$nn.idx[i,]
      sel_ss[[i]]=coords[sel,]
      numpoints[[i]]=ncol(sel_ss[[i]])
      if(!is.null(M)) {MM[[i]]=cbind(M[ss,1],M[ss,2])}
      if(!is.null(data))data_sel[[i]]=matrix(data[,sel],nrow=2)
      if(!is.null(X))   XX[[i]]=X[c(sel,2*sel),]
    }
  }
  #####################################################################################
  if(!dyn)   #### to fix case dyn
  {
    Nloc=nrow(loc)
    if(dyn) coords=do.call(rbind,coordx_dyn)
    sel_ss=numpoints=data_sel=sel_tt=XX=MM=list()
    if(dyn)      {out1= nbor(coords[[1]],loc,distance,maxdist,neighb)
                 out2= nbor(coords[[2]],loc,distance,maxdist,neighb)}
    else  {out= nbor(coords,loc,distance,maxdist,neighb)}
     
    for(i in 1:Nloc){
      sel=out$nn.idx[i,]
      sel_ss[[i]]=coords[sel,]
      numpoints[[i]]=ncol(sel_ss[[i]])
      if(!is.null(M)) {MM[[i]]=cbind(M[ss,1],M[ss,2])}
      if(!is.null(data))data_sel[[i]]=matrix(data[,sel],nrow=2)
      if(!is.null(X))   XX[[i]]=X[c(sel,2*sel),]
    }
  }
}
  ##################################################################################
  ##################################################################################
  ##################################################################################
  if(length(MM)==0)  MM=NULL
  if(length(XX)==0) XX=NULL
  if(length(data_sel)==0) data_sel=NULL
 
  return(list(data=data_sel,coordx=sel_ss,coordt=sel_tt,distance=distance, 
              numpoints=numpoints,numtime=numtime,radius=radius,spacetime=spacetime,X=XX,M=MM))
} 