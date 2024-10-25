####################################################
### File name: GeoNeighbSelect.r
####################################################

GeoNeighbSelect <- function(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL,
  copula=NULL,corrmodel=NULL, distance="Eucl",fixed=NULL,anisopars=NULL,
  est.aniso=c(FALSE,FALSE), grid=FALSE, likelihood='Marginal', 
  lower=NULL,maxdist_v=Inf,neighb=c(1,2,3,4,5),
  maxtime=Inf, maxtime_v=Inf, memdist=TRUE,model='Gaussian',n=1,numbins=20, ncores=NULL,
  optimizer='Nelder-Mead', parallel=FALSE, bivariate=FALSE,
  radius=6371, start=NULL,  type='Pairwise', upper=NULL,  weighted=FALSE,
  X=NULL,nosym=FALSE,spobj=NULL,spdata=NULL)
{

if(!is.numeric(neighb))  stop("neighb must be a numeric vector")
if(sum(neighb-floor(neighb))) stop("neighb must be a positive integer numeric vector")

estimates=best_T=best_K= NULL

############## computing semivariogram ###############
 semiv=GeoVariogram(data=data, coordx=coordx, coordy=coordy, coordt=coordt, 
coordx_dyn=coordx_dyn, distance=distance,
              grid=grid, maxdist=maxdist_v,neighb=NULL,
              maxtime=maxtime_v, numbins=numbins, 
              radius=radius, type='variogram',bivariate=bivariate)
######################################################
if(!is.null(try(CkCorrModel(corrmodel), TRUE)))
{
    bivariate<-CheckBiv(CkCorrModel(corrmodel))
    spacetime<-CheckST(CkCorrModel(corrmodel))
    space=!spacetime&&!bivariate
}else {stop("correlation model is not valid\n")}


####################
coremax=parallel::detectCores()
if(is.na(coremax)||coremax==1) parallel=FALSE
####################


K=length(neighb); res=double(K)

P=NULL
if(spacetime){
   if(!is.numeric(maxtime))  stop("neighb must be a numeric vector")
   P=length(maxtime)
   res=double(K*P)
}
#######################################################################################
##################################### SPATIAL #########################################
#######################################################################################
if(space) {

M=1
##################################### spatial not parallel ############################
if(!parallel)
{
for(M in 1:K) {
  aa=GeoFit(data=data, coordx=coordx, coordy=coordy, coordt=coordt, coordx_dyn=coordx_dyn,copula=copula,corrmodel=corrmodel, distance=distance,
                         fixed=fixed,anisopars=anisopars,est.aniso=est.aniso, grid=grid, likelihood=likelihood, 
                         lower=lower,neighb=neighb[M],
                          maxtime=maxtime, memdist=memdist,model=model,n=n, 
                          optimizer=optimizer, 
                         radius=radius, start=start,  
                         type=type, upper=upper,  weighted=weighted,X=X,nosym=nosym,spobj=spobj,spdata=spdata)
 clest=  aa$param
 estimates=rbind(estimates,unlist(clest))

 

 if(aa$convergence=="Successful")
 {
  ## first method slightly faster
 #cc=GeoCorrFct(semiv$centers,t=semiv$centert,corrmodel=corrmodel, model=model,distance=distance, param=c(aa$param,aa$fixed),radius=radius,n=n,covariance=TRUE,variogram=TRUE)$corr
 # res[M]=sum(cc - semiv$variograms)^2

 ## second method slighty slower but with graphics..
 cc=GeoCovariogram(fitted=aa,distance=distance,show.vario=TRUE, vario=semiv,pch=20,invisible=TRUE)
 res[M]=cc 
 }
 else { res[M]=Inf}

 }
}
##################################### end spatial not parallel ############################
##################################### spatial parallel ############################
if(parallel)
{
if(is.null(ncores)){ n.cores <- coremax - 1 }
else
{  if(!is.numeric(ncores)) stop("number of cores not valid\n")
   if(ncores>coremax||ncores<1) stop("number of cores not valid\n")
   n.cores=ncores
}
cat("Performing",K,"estimations using",n.cores,"cores...\n")

future::plan(multisession, workers = n.cores)

progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")
pb <- progressr::progressor(along = 1:K)
k=1
xx=foreach::foreach(k = 1:K,.combine = rbind,
                           .options.future = list(seed = TRUE)) %dofuture% 
    {  

aa=GeoFit(data=data, coordx=coordx, coordy=coordy, coordt=coordt, coordx_dyn=coordx_dyn,copula=copula,corrmodel=corrmodel, distance=distance,
                         fixed=fixed,anisopars=anisopars,est.aniso=est.aniso, grid=grid, likelihood=likelihood, 
                         lower=lower,neighb=neighb[k],
                          maxtime=maxtime, memdist=memdist,model=model,n=n, 
                          optimizer=optimizer,  
                         radius=radius, start=start,  
                         type=type, upper=upper,  weighted=weighted,X=X,nosym=nosym,spobj=spobj,spdata=spdata)
  pb(sprintf("k=%g", k)) 
     xx=GeoCovariogram(fitted=aa,distance=distance,show.vario=TRUE, vario=semiv,pch=20,invisible=TRUE)

     }   
res=as.numeric(xx)
future::plan(sequential)
}

}
#######################################################################################
##################################### SPATIO TEMPORAL ################################
#######################################################################################


if(spacetime){
######################################## Space time  NO parallel ############################
if(!parallel)
{
  F=1
for(L in 1:P) { 
for(M in 1:K) {
  aa=GeoFit(data=data, coordx=coordx, coordy=coordy, coordt=coordt, coordx_dyn=coordx_dyn,copula=copula,corrmodel=corrmodel, distance=distance,
                         fixed=fixed,anisopars=anisopars,est.aniso=est.aniso, grid=grid, likelihood=likelihood, 
                         lower=lower,neighb=neighb[M],
                          maxtime=maxtime[L], memdist=memdist,model=model,n=n, 
                          optimizer=optimizer, 
                         radius=radius, start=start,  
                         type=type, upper=upper,  weighted=weighted,X=X,nosym=nosym,spobj=spobj,spdata=spdata)
 clest=  aa$param
 estimates=rbind(estimates,unlist(clest))
 
 if(aa$convergence=="Successful")
 {
  ## first method slightly faster
 #cc=GeoCorrFct(semiv$centers,t=semiv$centert,corrmodel=corrmodel, model=model,distance=distance, param=c(aa$param,aa$fixed),radius=radius,n=n,covariance=TRUE,variogram=TRUE)$corr
 # res[M]=sum(cc - semiv$variograms)^2

 ## second method slighty slower but with graphics..
 cc=GeoCovariogram(fitted=aa,distance=distance,show.vario=TRUE, vario=semiv,pch=20,fix.lagt=1,fix.lags=1,invisible=TRUE)
 res[F]=cc; F=F+1 
 }
 else { res[F]=Inf; F=F+1 }
}}
}
######################################## end Space time NO parallel ############################
######################################## Space time   parallel ############################
if(!parallel)
{

}
######################################## end Space time   parallel ############################
}
################# end SPATIOTEMPORAL ######################################################

if(space) {
indexmin=which.min(res)
bestK=neighb[indexmin]
}
if(spacetime)
{
  indexmin=which.min(res)
  bestKT=c(as.matrix(expand.grid(neighb,maxtime))[indexmin,])
  bestK=as.numeric(bestKT[1]);best_T=as.numeric(bestKT[2]);
}

a=list(best_neighb=bestK,best_maxtime=best_T,res=res,estimates=estimates,best_est=estimates[indexmin,])

return(a)
}