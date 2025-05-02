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

  ##if(space){
   # sel_ss=data_sel=numpoints=XX=MM=list()
   # out= nbor(coords,loc,distance,maxdist,neighb)
    #for(i in 1:Nloc)
   # {
   #   ss=out$nn.idx[i,];
    #  sel_ss[[i]]=coords[ss,]
   #   numpoints[[i]]=nrow(sel_ss[[i]])
   #   if(!is.null(data)) data_sel[[i]]=data[ss]
    #  if(!is.null(M)) MM[[i]]=M[ss]
   #   if(!is.null(X)) XX[[i]]=X[ss,]
   # }
  #}

  if(space) {
  sel_ss <- data_sel <- numpoints <- XX <- MM <- vector("list", Nloc)
  out <- nbor(coords, loc, distance, maxdist, neighb)
  nn_idx <- out$nn.idx
  has_data <- !is.null(data)
  has_M <- !is.null(M)
  has_X <- !is.null(X)
  for(i in seq_len(Nloc)) {
    ss <- nn_idx[i, ]
    sel_ss[[i]] <- coords[ss, , drop = FALSE]
    numpoints[[i]] <- nrow(sel_ss[[i]])
    if(has_data) data_sel[[i]] <- data[ss]
    if(has_M) MM[[i]] <- M[ss]
    if(has_X) XX[[i]] <- X[ss, , drop = FALSE]
  }
}
#####################################################################################
if (spacetime) {
  if (!dyn) {  # caso non dinamico ottimizzato
    Tloc <- length(time)
    Nloc <- nrow(loc)
    TT <- length(coordt)
    NN <- nrow(coords)
    
    # Pre-allocazione e preprocessing
    Mma <- if (!is.null(M)) matrix(M, nrow = TT, ncol = NN, byrow = TRUE) else NULL
    Xxa <- if (!is.null(X)) lapply(seq_len(ncol(X)), function(i) matrix(X[, i], nrow = TT, ncol = NN, byrow = TRUE)) else NULL
    
    out <- nbor(coords, loc, distance, maxdist, neighb)
    out_t <- nabor::knn(coordt, time, k = length(coordt), radius = maxtime)
    
    # Preprocessamento degli indici temporali
    qq <- out_t$nn.idx[1, ]
    ss <- which(qq != 0)
    out_t_ind <- cbind(as.vector(out_t$nn.idx[, ss]), time)
    
    total_points <- Nloc * Tloc
    sel_ss <- vector("list", total_points)
    numpoints <- integer(total_points)
    data_sel <- vector("list", total_points)
    sel_tt <- vector("list", total_points)
    XX <- vector("list", total_points)
    MM <- vector("list", total_points)
    
    # Funzioni helper per evitare ripetizioni
    get_time_indices <- function(j) {
      sel_t <- out_t_ind[out_t_ind[, 2] == time[j], 1]
      if (length(sel_t) == 0) return(NULL)
      pp <- order(coordt[sel_t])
      list(indices = sel_t[pp], times = coordt[sel_t][pp])
    }
    
    # Loop ottimizzato
    k <- 1
    for (i in seq_len(Nloc)) {
      sel_s <- out$nn.idx[i, ]
      coords_sel <- matrix(coords[sel_s, ], ncol = 2)
      
      for (j in seq_len(Tloc)) {
        time_info <- get_time_indices(j)
        
        if (is.null(time_info)) {
          sel_ss[[k]] <- coords_sel
          numpoints[k] <- nrow(coords_sel)
          k <- k + 1
          next
        }
        
        sel_ss[[k]] <- coords_sel
        numpoints[k] <- nrow(coords_sel)
        sel_tt[[k]] <- time_info$times
        
        if (!is.null(data)) {
          data_sel[[k]] <- data[time_info$indices, sel_s, drop = FALSE][order(coordt[time_info$indices]), , drop = FALSE]
        }
        
        if (!is.null(M)) {
          MM[[k]] <- as.numeric(t(Mma[time_info$indices, sel_s, drop = FALSE][order(coordt[time_info$indices]), , drop = FALSE]))
        }
        
        if (!is.null(X)) {
          XX[[k]] <- do.call(cbind, lapply(Xxa, function(x) 
            as.numeric(t(x[time_info$indices, sel_s, drop = FALSE][order(coordt[time_info$indices]), , drop = FALSE]))))
        }
        
        k <- k + 1
      }
    }
  }  # fine non dinamico
  
  if (dyn) {  # caso dinamico ottimizzato
    if (length(coordx_dyn) != length(data))
      stop("Length of coordx_dyn and data must match in dynamic spacetime case.")
    
    Tloc <- length(time)
    Nloc <- nrow(loc)
    total_points <- Nloc * Tloc
    
    # Preprocessamento temporale
    out_temp <- nabor::knn(coordt, time, k = length(coordt), radius = maxtime)
    sel_t_all <- out_temp$nn.idx[1, ]
    sel_t_all <- sel_t_all[sel_t_all > 0]
    ord_t <- order(coordt[sel_t_all])
    sel_t_ord <- sel_t_all[ord_t]
    coordt_sel_ord <- coordt[sel_t_ord]
    
    # Preallocazione
    result <- vector("list", total_points)
    loc_matrix <- as.matrix(loc)
    
    # Funzione per processare un singolo punto spazio-temporale
    process_point <- function(i, j) {
      coords_list <- lapply(sel_t_ord, function(tt) {
        coords_t <- coordx_dyn[[tt]]
        out_spat <- nbor(coords_t, matrix(loc_matrix[i, ], ncol = ncol(coords_t)), 
                         distance, maxdist, neighb)
        ss <- out_spat$nn.idx[1, ]
        
        list(
          coords = coords_t[ss, , drop = FALSE],
          data = if (!is.null(data)) data[[tt]][ss] else NULL,
          M = if (!is.null(M)) M[[tt]][ss] else NULL,
          X = if (!is.null(X)) X[[tt]][ss, , drop = FALSE] else NULL
        )
      })
      
      list(
        sel_ss = coords_list[[1]]$coords,
        sel_tt = coordt_sel_ord,
        numpoints = nrow(coords_list[[1]]$coords),
        data_sel = if (!is.null(data)) do.call(rbind, lapply(coords_list, `[[`, "data")) else NULL,
        MM = if (!is.null(M)) unlist(lapply(coords_list, `[[`, "M")) else NULL,
        XX = if (!is.null(X)) do.call(rbind, lapply(coords_list, `[[`, "X")) else NULL
      )
    }
    
    # Applicazione parallela (se possibile) o con lapply
    result <- lapply(seq_len(total_points), function(k) {
      i <- ((k - 1) %/% Tloc) + 1
      j <- ((k - 1) %% Tloc) + 1
      process_point(i, j)
    })
    
    # Riorganizzazione dei risultati
    sel_ss <- lapply(result, `[[`, "sel_ss")
    sel_tt <- lapply(result, `[[`, "sel_tt")
    numpoints <- vapply(result, `[[`, integer(1), "numpoints")
    data_sel <- if (!is.null(data)) lapply(result, `[[`, "data_sel") else vector("list", total_points)
    MM <- if (!is.null(M)) lapply(result, `[[`, "MM") else vector("list", total_points)
    XX <- if (!is.null(X)) lapply(result, `[[`, "XX") else vector("list", total_points)
  }  # fine dinamico
}  # fine spacetime


#####################################################################################
if(bivariate) {
  Nloc <- nrow(loc)
  # Preallocazione delle liste base
  sel_ss <- numpoints <- data_sel <- sel_tt <- vector("list", Nloc)
  if(!is.null(M)) {MM <- vector("list", Nloc)} 
  else {MM <- if(dyn) list(NULL, NULL) else NULL}
  if(!is.null(X)) {XX <- vector("list", Nloc)} 
  else {XX <- if(dyn) list(NULL, NULL) else NULL}
  if(dyn) {
    # CASO DINAMICO ----------------------------------------------------------
    coords1 <- coords[[1]]
    coords2 <- coords[[2]]
    n1 <- nrow(coords1)
    out1 <- nbor(coords1, loc, distance, maxdist, neighb)
    out2 <- nbor(coords2, loc, distance, maxdist, neighb)
    has_data <- !is.null(data)
    for(i in seq_len(Nloc)) 
    {
      sel1 <- out1$nn.idx[i, ]
      sel2 <- out2$nn.idx[i, ]
      sel_ss[[i]] <- rbind(coords1[sel1, ], coords2[sel2, ])
      numpoints[[i]] <- length(sel1) + length(sel2)
      if(has_data) { data_sel[[i]] <- matrix(c(data[1, sel1], data[2, sel2]), nrow = 2)}
      if(!is.null(M)) {MM[[i]] <- list(M[[1]][sel1], M[[2]][sel2])}
      if(!is.null(X)) {XX[[i]] <- list(X[[1]][sel1, , drop = FALSE], X[[2]][sel2, , drop = FALSE])}
    }
  }
  else {
    # CASO NON DINAMICO ------------------------------------------------------
    out <- nbor(coords, loc, distance, maxdist, neighb)
    n <- nrow(coords)
    has_data <- !is.null(data)
    for(i in seq_len(Nloc)) {
      sel <- out$nn.idx[i, ]
      sel_ss[[i]] <- coords[sel, ]
      numpoints[[i]] <- length(sel)
      if(has_data) {
        data_sel[[i]] <- matrix(data[, sel], nrow = 2)
      }
      if(!is.null(M)) { MM[[i]] <- c(M[sel], M[n + sel]) }  # Combina i due vettori
      if(!is.null(X)) { XX[[i]] <- rbind(X[sel, , drop = FALSE], X[n + sel, , drop = FALSE])}
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