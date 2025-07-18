####################################################
GeoNeighIndex<-function(coordx,coordy=NULL,coordz=NULL,coordt=NULL,coordx_dyn=NULL,
                              distance="Eucl",neighb=4,maxdist=NULL,maxtime=1,radius=1,bivariate=FALSE)
{

indices <- function(X, Y) {
  # Versione vettorizzata senza loop
  n <- ncol(X) - 1
  xy <- cbind(rep(X[, 1], times = n),as.vector(X[, -1]))
  d <- as.vector(Y[, -1])
  list(xy = xy, d = d)
}
##############################################################
nn2Geo <- function(x, y, K = 1, distance = 0, maxdist = NULL, radius=1) {
  # Controlli input base
  if (!is.matrix(x) || !is.matrix(y)) stop("x e y devono essere matrici")
  n_x <- nrow(x)
  n_y <- nrow(y)
  
  # Imposta K massimo per evitare errori
  K <- min(K, n_x, n_y)
  
  # Funzione helper per proiezione sinusoidale
  project_coords <- function(coords) {
    prj <- mapproj::mapproject(coords[,1], coords[,2], projection = "sinusoidal")
    radius * cbind(prj$x, prj$y)
  }

  # Calcolo nearest neighbors con o senza maxdist
  if (is.null(maxdist)) {
    nearest <- nabor::knn(x, y, k = K)
  } else {
    if (distance %in% c(1, 2)) {
      x_proj <- project_coords(x)
      y_proj <- project_coords(y)
      nearest <- nabor::knn(x_proj, y_proj, k = K, radius = maxdist)
    } else {
      nearest <- nabor::knn(x, y, k = K, radius = maxdist)
    }
  }
  # Calcolo distanze geografiche se richiesto
  if (distance %in% c(1, 2)) {
    # Funzione vettoriale per calcolare distanze geografiche
    calc_dist_row <- function(i) {
      idx <- nearest$nn.idx[i, ]
      valid_idx <- idx > 0
      if (!any(valid_idx)) return(numeric(0))
      fields::rdist.earth.vec(
        x1 = matrix(x[i, ], nrow = 1),
        x2 = x[idx[valid_idx], , drop = FALSE],
        miles = FALSE, R = 1
      )
    }
    # Applico la funzione a tutte le righe (vettorizzazione parziale)
    dists_list <- lapply(seq_len(n_x), calc_dist_row)
    # Preallocazione matrice distanze con NA
    dists_mat <- matrix(NA_real_, nrow = n_x, ncol = K)
    for (i in seq_len(n_x)) {
      valid_len <- length(dists_list[[i]])
      if (valid_len > 0) dists_mat[i, seq_len(valid_len)] <- dists_list[[i]]
    }
    # Imposta distanza a zero per primo vicino (se è se stesso)
    dists_mat[, 1] <- 0
    # Calcola distanza finale in metri o approssimata
    if (distance == 2) {
      nearest$nn.dists <- radius * dists_mat
    } else {
      nearest$nn.dists <- 2 * radius * sin(0.5 * dists_mat)
    }
  }
  # Estrai indici e distanze tramite funzione indices (presupposta definita)
  sol <- indices(nearest$nn.idx, nearest$nn.dists)
  # Filtra in base a maxdist se specificato
  if (!is.null(maxdist)) {
    sel <- sol$xy[, 2] > 0
    return(list(lags = sol$d[sel], rowidx = sol$xy[sel, 1], colidx = sol$xy[sel, 2]))
  } else {
    return(list(lags = sol$d, rowidx = sol$xy[, 1], colidx = sol$xy[, 2]))
  }
}


##############################################################
spacetime_index <- function(coords, coordx_dyn = NULL, N, K = 4, coordt = NULL,
                            numtime, maxtime = 1, maxdist = NULL, 
                            distance = "Eucl", radius=1) {
  
  m_s <- vector("list", numtime)
  
  # -------------------------------
  # 1. Marginal Spatial Indexes
  # -------------------------------
  if (is.null(coordx_dyn)) {
    inf <- nn2Geo(coords, coords, K + 1, distance, maxdist, radius)
    aa <- cbind(inf$rowidx, inf$colidx)
    lag <- inf$lags
    offset_seq <- N * (0:(numtime - 1))
    for (i in seq_len(numtime)) {
      temp <- cbind(aa + offset_seq[i], 0L, lag)
      m_s[[i]] <- temp
    }
  } else {
    ns <- vapply(coordx_dyn, nrow, FUN.VALUE = integer(1))
    for (i in seq_len(numtime)) {
      inf <- nn2Geo(coordx_dyn[[i]], coordx_dyn[[i]], K + 1, distance, maxdist, radius)
      aa <- cbind(inf$rowidx, inf$colidx)
      lag <- inf$lags
      temp <- cbind(aa + ns[i] * (i - 1L), 0L, lag)
      m_s[[i]] <- temp
    }
  }

  # -------------------------------
  # 2. Temporal & Spatio-Temporal
  # -------------------------------
  time_dists <- sort(unique(nabor::knn(coordt, coordt, k = round(maxtime) + 1)$nn.dists))
  nn <- time_dists[time_dists > 0]
  tnn <- length(nn)
  
  m_t <- vector("list", tnn * (numtime - tnn))
  m_st <- vector("list", tnn * (numtime - tnn))
  
  counter <- 1L
  for (j in seq_len(tnn)) {
    for (k in seq_len(numtime - j)) {
      n1 <- nrow(m_s[[k]])
      n2 <- nrow(m_s[[k + j]])
      bb <- min(n1, n2)

      m_t[[counter]] <- cbind(
        m_s[[k]][1:bb, 1], 
        m_s[[k + j]][1:bb, 1], 
        rep_len(nn[j], bb), 
        rep_len(0L, bb)
      )
      m_st[[counter]] <- cbind(
        m_s[[k]][1:bb, 1], 
        m_s[[k + j]][1:bb, 2], 
        rep_len(nn[j], bb), 
        m_s[[k]][1:bb, 4]
      )
      counter <- counter + 1L
    }
  }

  # -------------------------------
  # 3. Final Result
  # -------------------------------
  final <- do.call(rbind, c(m_s, m_t, m_st))
  return(final)
}


bivariate_index <- function(coords, coordx_dyn = NULL, N, K = 4, maxdist, distance, radius) {
  # Gestione parametri K e maxdist
  if (length(K) == 3) {
    K1 <- K[1]; K2 <- K[2]; K3 <- K[3]
  } else {
    K1 <- K2 <- K3 <- K
  }
  
  if (length(maxdist) == 3) {
    maxdist1 <- maxdist[1]; maxdist2 <- maxdist[2]; maxdist3 <- maxdist[3]
  } else {
    maxdist1 <- maxdist2 <- maxdist3 <- maxdist
  }
  
  if (is.null(coordx_dyn)) {
    # Caso statico
    n_half <- as.integer(N / 2)
    cm <- coords[1:n_half, ]
    # Calcolo una volta sola i risultati che vengono riutilizzati
    inf1 <- nn2Geo(cm, cm, K1 + 1, distance, maxdist1, radius)
    inf2 <- nn2Geo(cm, cm, K3 + 1, distance, maxdist3, radius)
    inf3 <- nn2Geo(cm, cm, K2 + 1, distance, maxdist2, radius)
    # Preallocazione e costruzione matrici
    aa1 <- cbind(inf1$rowidx, inf1$colidx, 0L, 0L, inf1$lags)
    aa2 <- cbind(inf2$rowidx + n_half, inf2$colidx + n_half, 1L, 1L, inf2$lags)
    aa3 <- cbind(inf3$rowidx, inf3$colidx + n_half, 0L, 1L, inf3$lags)
    aa4 <- cbind(inf3$rowidx + n_half, inf3$colidx, 1L, 0L, inf3$lags)
    
    # Ottimizzazione: calcolo knn una sola volta
    a5 <- nabor::knn(cm, k = K2)
    aa5 <- cbind(a5$nn.idx[, 1], a5$nn.idx[, 1] + n_half, 0L, 1L, 0)
    aa6 <- cbind(a5$nn.idx[, 1] + n_half, a5$nn.idx[, 1], 1L, 0L, 0)
  } else {
    # Caso dinamico
    ns <- vapply(coordx_dyn, nrow, integer(1))
    # Calcolo una volta sola i risultati che vengono riutilizzati
    inf1 <- nn2Geo(coordx_dyn[[1]], coordx_dyn[[1]], K1 + 1, distance, maxdist1, radius)
    inf2 <- nn2Geo(coordx_dyn[[2]], coordx_dyn[[2]], K3 + 1, distance, maxdist3, radius)
    inf3 <- nn2Geo(coordx_dyn[[1]], coordx_dyn[[2]], K2 + 1, distance, maxdist2, radius)
    # Preallocazione e costruzione matrici
    aa1 <- cbind(inf1$rowidx, inf1$colidx, 0L, 0L, inf1$lags)
    aa2 <- cbind(inf2$rowidx + ns[1], inf2$colidx + ns[1], 1L, 1L, inf2$lags)
    aa3 <- cbind(inf3$rowidx, inf3$colidx + ns[1], 0L, 1L, inf3$lags)
    aa4 <- cbind(inf3$rowidx + ns[1], inf3$colidx, 1L, 0L, inf3$lags)
    # Ottimizzazione: calcolo knn una sola volta
    a5 <- nabor::knn(coordx_dyn[[1]], coordx_dyn[[2]], k = K2 + 1)
    aa5 <- cbind(a5$nn.idx[, 1], a5$nn.idx[, 1] + ns[1], 0L, 1L, 0)
    aa6 <- cbind(a5$nn.idx[, 1] + ns[1], a5$nn.idx[, 1], 1L, 0L, 0)
  }
  
  # Combinazione finale
  rbind(aa1, aa2, aa3, aa4, aa5, aa6)
}



######################################
######### start ######################
######################################

    spatial=TRUE
    spacetime=FALSE
    
    if(!is.null(coordx_dyn))  
                if(!is.list(coordx_dyn)) stop(" coordx_dyn must be a list")
    ### Check the parameters given in input
    # Checks if its a spatial or spatial-temporal random field:
    if(!is.null(coordt))
    if(is.numeric(coordt)&&is.numeric(maxtime)) 
    if(length(coordt)>1&&length(maxtime)>=1)  spacetime=TRUE
    distance=CheckDistance(distance)
    spatial= !spacetime&&!bivariate

K=neighb
## for spacetime or bivariate
#if(!is.null(coordx_dyn)){
      if(!is.null(coordy)){
               if(is.null(coordz)){
                          coordy <- coordx[,2]
                          coordx <- coordx[,1]
                          coords=cbind(coordx,coordy)
                                  }
                else {    coordz <- coordx[,3]
                          coordy <- coordx[,2]
                          coordx <- coordx[,1]
                              coords=cbind(coordx,coordy,coordz)
                     }        
                numcoord=nrow(coords)
           }
      else {
                       if(!bivariate)
                       { coords=coordx;  numcoord=nrow(coords) }
                       if(bivariate) {
                        if(is.null(coordx_dyn)) {coords=coordx; 
                                                 numcoord=nrow(coords) }
                        else {coords=1;numcoord=1}
                      }
           }  
#}             
##########################
########################## 
if(spatial)   #  spatial case
{
##########################################

  sol = nn2Geo(coords,coords,K+1 ,distance,maxdist,radius) ##### 
  gb=list(); gb$colidx=sol$colidx;
             gb$rowidx=sol$rowidx ;
             gb$lags=sol$lags
            # gb$lagt=NULL
             gb$maxdist=maxdist
             gb$neighb=neighb
} #### end spatial case 
##############################################   
if(spacetime)   #  space time  case
{ 
  numtime=length(coordt)
  sol=spacetime_index(coords,coordx_dyn, numcoord,K,coordt,numtime,maxtime,maxdist,distance,radius)
  gb=list(); gb$colidx=sol[,2];
             gb$rowidx=sol[,1] ;
             gb$lags=sol[,4]
             gb$lagt=sol[,3]
} 
if(bivariate)  { #space bivariate  case
   sol=bivariate_index(coords,coordx_dyn,numcoord,K,maxdist,distance,radius)

   gb=list(); gb$colidx=sol[,2];
             gb$rowidx=sol[,1] ;
             gb$lags=sol[,5]
             gb$first=sol[,3]
             gb$second=sol[,4]
             gb$maxdist=maxdist
             gb$neighb=neighb
} 
return(gb)
}

