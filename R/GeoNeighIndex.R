####################################################
GeoNeighIndex<-function(coordx,coordy=NULL,coordx_dyn=NULL,coordt=NULL,
                              distance="Eucl",neighb=4,maxdist=NULL,maxtime=1,radius=6371,bivariate=FALSE)
{

fxy <- function(x,y, tol = 15){
 
 
  xx = atan(x)/(pi/2)
  yy = atan(y)/(pi/2)

 
  xdig = (as.numeric(strsplit(as.character(xx), "")[[1]][-(1:2)]))
  ydig = (as.numeric(strsplit(as.character(yy), "")[[1]][-(1:2)]))
  # length(xdig);length(ydig);
  if (length(xdig) < tol){
    xdig = (as.numeric(strsplit(as.character(xx-10^(-tol)), "")[[1]][-(1:2)]))
  }
  if (length(ydig) < tol){
    ydig = (as.numeric(strsplit(as.character(yy-10^-(tol)), "")[[1]][-(1:2)]))
  }
 
  xdig = xdig[1:tol]
  ydig = ydig[1:tol]
 
  if (length(xdig) > tol){print(c(x,y))}
  if (length(ydig) > tol){print(c(x,y))}
 
 
  if (y>=x){
   
    z = paste0(c('0.',as.vector(rbind(xdig,ydig))), collapse = "")
  }else{
    z = paste0(c('0.',as.vector(rbind(ydig,xdig))), collapse = "")
   
  }  
  bol = (z)
  if(is.na(bol))
  {
    cat("\n input: \n",c(x,y),"\n")
    cat("xdig: \n",c(xdig),"\n")
    cat("ydig: \n",c(ydig),"\n\n")
  }
  return(bol )
}
fxy <- Vectorize(fxy)
#########################################
indices <- function(X,Y)
 {
             res = NULL;res_d = NULL
             for(i in 2:ncol(X))
             {
                sol = cbind(X[,1],X[,i])
                res = rbind(res,sol)
                sol_d = cbind(Y[,1],Y[,i])
                res_d = rbind(res_d,sol_d)
             }
         
            return(list(xy = res,d = res_d[,2]))
 }

##########################################
nn2Geo <- function(x,y, K = 1,distance=0,maxdist=NULL,radius=6371)  
  {
            if(is.null(maxdist)) 
               {
               #nearest = RANN::nn2(x,y,k = K,treetype = c("kd"))} ### case neighboord
               nearest = nabor::knn(x,y,k = K)} ### case neighboord
            else     {
                    
                     K=min(K-1,nrow(x)) # case of  maxdist 
                    # nearest = RANN::nn2(x,y,searchtype = c("radius"),
                     #          treetype = c("kd"),radius = maxdist,k=K  )
                       nearest = nabor::knn(x,y,radius = maxdist,k=K  )

                     }
            #########  cases geod (2) or chordal (1) distances :  to improve this  code!!
            if(distance==2||distance==1){
                  nn=nrow(x); 
                  nnd=ncol(nearest$nn.dists)
                  mm=matrix(0,nrow=nn,ncol=nnd)
                  for(i in 1:nn){   ## can we improve that?
                  si=nearest$nn.idx[i,];sel1=si[si>0]
                  a=fields::rdist.earth.vec(x1=matrix(x[i,],ncol=2),
                                            x2=matrix(x[sel1,],ncol=2), miles = FALSE, R = 1)
                  mm[i,][1:length(a)]=a
                
                 }
              mm[,1]=0 # just to be sure
             if(distance==2)  mm=radius*mm   # geodesic
             if(distance==1)  mm=2*radius*sin(0.5*mm)   # chordal  
             nearest$nn.dists=mm
             }
            ########################################### 
            sol = indices(nearest$nn.idx,nearest$nn.dists)
            if(is.null(maxdist)) lags <- sol$d;rowidx <- sol$xy[,1];colidx <- sol$xy[,2]
            if(!is.null(maxdist)){
                                    sel = sol$xy[,2]>0
                                    lags=sol$d[sel];rowidx <- sol$xy[,1][sel];colidx <- sol$xy[,2][sel]
                                 }
         return(list (lags=lags, rowidx = rowidx, colidx = colidx))
   }
#########################

spacetime_index=function(coords,coordx_dyn=NULL,N,K=4,coordt=NULL
                         ,numtime,maxtime=1,maxdist=NULL,distance="Eucl",radius=6371)
{
  
  ##############
  m_s=list();m_t=m_st=NULL;
  ##############         
  ## building marginal spatial indexes
 # tt0 <- proc.time()
  if(is.null(coordx_dyn)) 
  {
    inf=nn2Geo(coords,coords,K+1,distance,maxdist,radius)
    aa=cbind(inf$rowidx,inf$colidx)   ## spatial index (fixed coordinates) 
    for(i in 1:numtime) {
      # i = 1
      # repito las coordenadas numtime veces en elementos de una lista
      m_s[[i]]=data.frame(cbind(aa+N*(i-1),0,inf$lags))
      
    }
  }

  if(!is.null(coordx_dyn))
  {        ns=lengths(coordx_dyn)/2 
  for(i in 1:numtime){
    inf=nn2Geo(coordx_dyn[[i]],coordx_dyn[[i]],K+1,distance,maxdist,radius)
    aa=cbind(inf$rowidx,inf$colidx)
    m_s[[i]]=cbind( aa+ns[i]*(i-1),0,inf$lags)    ## spatial index (dynamic coordinates)
  }
  }
  ## building  temporal  and spatiotemporal indexes
  
  ## temporal distances (not zero distance)
  #nn=sort(unique(c(RANN::nn2(coordt,coordt,k=round(maxtime)+1,treetype = c("kd"))$nn.dists)))[-1]  
  #nn=sort(unique(c(nabor::knn(coordt,coordt,k=round(maxtime)+1))$nn.dists))[-1]  

   ## first way
    #qq=c(nabor::knn(coordt,coordt,k=length(coordt),radius=maxtime)$nn.dists)
    #qq=qq[is.finite(qq)];a=sort(unique(qq)); 
    #nn=a[a>0]

    ## second way
      a=sort(unique(c(nabor::knn(coordt,coordt,k=round(maxtime)+1)$nn.dists)))
   nn=a[a>0]


  tnn=length(nn)   
  # sol <- NULL
  m_t <- list()
  m_st <- list()
  contador <- 1
  
  for(j in 1:tnn){
    for(k in 1:(numtime-tnn)){
      # j = 1;k = 1
      bb=nrow(m_s[[k]])
      m_t[[contador]] =data.frame(cbind( m_s[[k]][,1], m_s[[k+j]][,1], rep(nn[j],bb),rep(0,bb)) )
      m_st[[contador]]=data.frame(cbind( m_s[[k]][,1], m_s[[k+j]][,2], rep(nn[j],bb), m_s[[k]][,4]) )
      contador  <- contador +1
     
    }
  }
  ######
  SS = data.table::rbindlist(m_s)
  TT = data.table::rbindlist(m_t)
  ST = data.table::rbindlist(m_st)
  ##final space-time indexes and distances
  final=data.table::rbindlist(list(SS,TT,ST))
  return(as.matrix(final))
}

#spacetime_index=function(coords,coordx_dyn,N,K,coordt,numtime,maxtime,maxdist,distance,radius)
#  {
##############
#m_s=list();m_t=m_st=NULL;
##############         
## building marginal spatial indexes
#if(is.null(coordx_dyn)) 
#   {
#        
#         inf=nn2Geo(coords,coords,K+1,distance,maxdist,radius)
#         aa=cbind(inf$rowidx,inf$colidx)   ## spatial index (fixed coordinates)
#         for(i in 1:numtime) {
#                  m_s[[i]]=cbind(aa+N*(i-1),0,inf$lags)
#                  }
#   }
#if(!is.null(coordx_dyn))
#  {        ns=lengths(coordx_dyn)/2 
#           for(i in 1:numtime){
#                  inf=nn2Geo(coordx_dyn[[i]],coordx_dyn[[i]],K+1,distance,maxdist,radius)
#                  aa=cbind(inf$rowidx,inf$colidx)
#                  m_s[[i]]=cbind( aa+ns[i]*(i-1),0,inf$lags)    ## spatial index (dynamic coordinates)
#                  }
 # }
         ## building  temporal  and spatiotemporal indexes
         
         ## temporal distances (not zero distance)
  #       nn=sort(unique(c(RANN::nn2(coordt,coordt,k=maxtime+1,treetype = c("kd"))$nn.dists)))[-1]  
  #       tnn=length(nn)   
  #       for(j in 1:tnn){
  #        for(k in 1:(numtime-tnn)){
  #          bb=nrow(m_s[[k]])
  #         m_t =rbind(m_t, cbind( m_s[[k]][,1], m_s[[k+j]][,1], rep(nn[j],bb)) )
  #         m_st=rbind(m_st,cbind( m_s[[k]][,1], m_s[[k+j]][,2], rep(nn[j],bb), m_s[[k]][,4]) )
  #       }}
  ##      ######
    #    TT=cbind(m_t,rep(0,nrow(m_t)))  
    #    SS=do.call(rbind,args=c(m_s));
     #   ST=m_st
      #  ##final space-time indexes and distances
       # final=rbind(SS,TT,ST)
       # return(final)
  #}

bivariate_index=function(coords,coordx_dyn,N,K,maxdist,distance,radius)
  {

 if(length(K)==1)  K1=K2=K3=K
 if(length(K)==3) {K1=K[1];K2=K[2];K3=K[3]}
if(is.null(coordx_dyn)) 
   {  
         inf1=nn2Geo(coords[1:(N/2),],  coords[1:(N/2),],  K1+1,distance,maxdist,radius)
         inf2=nn2Geo(coords[(N/2+1):N,],coords[(N/2+1):N,],K2+1,distance,maxdist,radius)
         inf3=nn2Geo(coords[1:(N/2),],  coords[(N/2+1):N,],K3+1,distance,maxdist,radius)
          aa1=cbind(inf1$rowidx,  inf1$colidx,0,0,inf1$lags)
          aa2=cbind(inf2$rowidx+N/2,inf2$colidx+N/2,1,1,inf2$lags)
          aa3=cbind(inf3$rowidx  ,inf3$colidx+N/2,0,1,inf3$lags)
          SS= cbind(rbind(aa1,aa2,aa3))            
   }
if(!is.null(coordx_dyn))
  {       
     ns=lengths(coordx_dyn)/2
    inf1=nn2Geo(coordx_dyn[[1]],  coordx_dyn[[1]],  K1+1,distance,maxdist,radius)
    inf2=nn2Geo(coordx_dyn[[2]],  coordx_dyn[[2]],  K2+1,distance,maxdist,radius)
    inf3=nn2Geo(coordx_dyn[[1]],  coordx_dyn[[2]],  K3+1,distance,maxdist,radius)
          aa1=cbind(inf1$rowidx,  inf1$colidx,0,0,inf1$lags)
          aa2=cbind(inf2$rowidx+ns[1],inf2$colidx+ns[1],1,1,inf2$lags)
          aa3=cbind(inf3$rowidx  ,inf3$colidx+ns[1],0,1,inf3$lags)
          SS= cbind(rbind(aa1,aa2,aa3))   
  }
##final bivariate  indexes and distances
return(SS)
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
      if(!is.null(coordy)){coordy <- coordx[,2]
                          coordx <- coordx[,1]
                          coords=cbind(coordx,coordy)
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
             gb$lagt=NULL
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
             gb$lagt=NULL
             gb$first=sol[,3]
             gb$second=sol[,4]
} 
return(gb)
}

