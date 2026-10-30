#include "header.h"


double clamp_corr(double x)
{
    const double eps = 1e-10;

    if (!R_FINITE(x)) return R_NaN;

    if (x >  1.0 - eps) return  1.0 - eps;
    if (x < -1.0 + eps) return -1.0 + eps;

    return x;
}
 
/***********************Legendre polinomial****************************************/

void lgnd (int *lmax,double *x, double *p){

    int l;
    p[0] = 1;
    p[1] = x[0];
    for (l = 1; l < *lmax; ++l)
    {
        p[l+1] = ((2*l+1)*x[0]*p[l]-l*p[l-1])/(l+1);
    }
  }
/*************************************************************************************/
/**************** matrix functions ***************************************************/
/*************************************************************************************/


void ludcmp(double **a, int n, int *indx, double *dd)
{
    int i,imax=0,j,k;
    double   big,dum,sum,temp,*vv;
    vv=(double *) R_Calloc(n,double);
    *dd = 1.0;
    for (i=0;i<n;i++)
    {
        big = 0.0;
        for (j=0;j<n;j++)
        {if ((temp=fabs(a[i][j])) > big) big = temp;}
        vv[i] = 1.0/big;}
    for (j=0;j<n;j++){
        for (i=0;i<j;i++){
            sum = a[i][j];
            for (k=0;k<i;k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;}
        big = 0.0;
        for (i=j;i<n;i++){
            sum = a[i][j];
            for (k=0;k<j;k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            if ((dum=vv[i]*fabs(sum)) >= big){
                big = dum;
                imax = i;}}
        if (j != imax){
            for (k=0;k<n;k++){
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;}
            *dd = -(*dd);
            vv[imax] = vv[j];}
        indx[j] = imax;
        if (a[j][j] == 0.0) a[j][j] = MAXERR;
        if (j != n-1){
            dum = 1.0 / a[j][j];
            for (i=j+1;i<n;i++) a[i][j] *= dum;}
    }
    R_Free(vv);
}


void lubksb(double **a, int n, int *indx, double *b)
{
    int i,ip,j,ii=-1;
    double   sum;
    for (i=0;i<n;i++){
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii>=0)
            for (j=ii;j<i;j++) sum -= a[i][j] * b[j];
        else if (sum) ii = i;
        b[i] = sum;}
    for (i=n-1;i>=0;i--){
        sum = b[i];
        for (j=i+1;j<n;j++) sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];}
}



/*compute the quadratic form given a simmetric square matrix*/
double QFORM(double **A,double *x,double *y,int n)
{
    int i,j;double qf = 0;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            qf=qf + x[i]*A[i][j]*y[j];}}
    return(qf);
}


/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

double hypot3d(double x, double y, double z) {
    double max_val = fmax(fabs(x), fmax(fabs(y), fabs(z)));
    if (max_val == 0) return 0;  // Se tutti i valori sono zero, la distanza è zero.
    
    // Scala i valori rispetto al massimo per evitare overflow o underflow
    x /= max_val;
    y /= max_val;
    z /= max_val;

    return max_val * sqrt(x * x + y * y + z * z);
}

double Dist_geodesic(double loni, double lati, double lonj, double latj, double radius) {
    if (loni == lonj && lati == latj) return 0.0;

    // Convert coordinates to radians once
    double ai = lati * M_PI / 180.0;
    double bi = loni * M_PI / 180.0;
    double aj = latj * M_PI / 180.0;
    double bj = lonj * M_PI / 180.0;

    // Calculate the geodesic distance using the spherical law of cosines
    double val = sin(ai) * sin(aj) + cos(ai) * cos(aj) * cos(bi - bj);
    val = fmin(fmax(val, -1.0), 1.0);  // Ensure the value is within the domain of acos
    return acos(val) * radius;
}




double Dist_chordal(double loni, double lati, double lonj, double latj, double radius) {
    // Compute the geodesic distance first and then use it to compute the chordal distance
    double geodesic_dist = Dist_geodesic(loni, lati, lonj, latj, radius);
    return radius * 2.0 * sin(geodesic_dist / (2.0 * radius));
}



double dist(int type_dist, double coordx, double locx, double coordy, double locy, double coordz, double locz, double radius) {
    double lags = 0.0;
    switch(type_dist) {

          case 0: // Euclidean distance
    if (fabs(coordz) < EPS1 && fabs(locz) < EPS1) {
        // 2D case
        lags = hypot(coordx - locx, coordy - locy);
    } else {
        // 3D case
        lags = hypot3d(coordx - locx, coordy - locy, coordz - locz);
    }
    break;

     /*   case 0: 
            lags = hypot(coordx - locx, coordy - locy);
            if (fabs(coordz - locz) < EPS1) {
                lags = hypot3d(coordx - locx, coordy - locy, coordz - locz);
            }
            break;*/
        case 1: // Chordal distance
            lags = Dist_chordal(coordx, coordy, locx, locy, radius);
            break;
        case 2: // Geodesic distance
            lags = Dist_geodesic(coordx, coordy, locx, locy, radius);
            break;
    }
    return lags;
}



// ---- wrapper che costruisce matrice distanze ----
void geo_distances(double *coords, int *ncoords, int *p, int *type_dist, double *radius,
                     double *out) {
    int n = *ncoords;
    int dim = *p;   // 2 o 3
    int td = *type_dist;
    double R = *radius;

    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            double d;
            if (dim == 2) {
                d = dist(td,
                         coords[i], coords[j],
                         coords[i + n], coords[j + n],
                         0.0, 0.0,
                         R);
            } else {
                d = dist(td,
                         coords[i], coords[j],
                         coords[i + n], coords[j + n],
                         coords[i + 2*n], coords[j + 2*n],
                         R);
            }
            out[i + j*n] = d;
            out[j + i*n] = d; // simmetria
        }
    }
}


void Maxima_Minima_dist(double *res,double *coordx,double *coordy,double *coordz,int *nsize,int *type_dist,double *radius)
{
  double res1=0.0,res2=-LOW,lags=0.0;
  int i=0,j=0;
    for(i=0; i<(*nsize-1);i++){
        for(j=(i+1); j<*nsize;j++){
          lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],coordz[i],coordz[j],*radius);
                   // lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],coordz[i],coordz[j],*radius);
          res1 = fmax(res1, lags);
          res2 = fmin(res2, lags);
}}
res[0]=res2;res[1]=res1;
return;
}





/*void RangeDist(double *max, double *min)
{
  *max=*maximdista;
  *min=*minimdista;
  return;
}*/


/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
// Computes the spatial distances 
//void Space_Dist(double *coordx,double *coordy,int *ia,int *idx,
//		int *ismal,int *ja,int *colidx,int *rowidx ,double thres) 

void Space_Dist(double *coordx,double *coordy,double *coordz,int *ia,int *idx,
        int *ismal,int *ja,int *colidx,int *rowidx ,double thres) 

{
  int i=0,h=0,j=0;
  double dij=0.0;
  /******************************************************************************/
  if(*istap){   // tapering case
      ia[0]=1;
	for(i=0;i<ncoord[0];i++)
	  for(j=0;j<ncoord[0];j++){
	    dij=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],coordz[i],coordz[j],*REARTH);
          //  dij=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],coordz[i],coordz[j],*REARTH);
	    if(dij<= thres){
	      tlags[h]=dij;
	      ja[h]=j+1;
	      idx[h]=i*(ncoord[0])+j+1;
	      ia[i+1]=ia[i+1]+1;
	      h++;}}
	for(i=0;i<ncoord[0];i++)
	  ia[i+1]=ia[i+1]+ia[i];
}  //end tapering
  /******************************************************************************/
else{  //no tapering
     h=0;
	  for(i=0;i<(ncoord[0]-1);i++){
	    for(j=(i+1);j<ncoord[0];j++){
      dij=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],coordz[i],coordz[j],*REARTH);
     if(dij<= thres){
        tlags[h]=dij;
        colidx[h]=i;  rowidx[h]=j; 
        h++;}
	      }}
   }
 // saving  spatial distances   
    *npairs=h;
    lags= (double *) R_Calloc(*npairs,double);
    for(i=0;i<*npairs;i++)  lags[i]=tlags[i];  
     //R_Free(tlags); tlags = NULL;
	return;
}



// Computes the spatial-temporal distances on regular and irregular grid:
void SpaceTime_Dist(double *coordx,double *coordy,double *coordz,double *coordt,int *ia,int *idx,int *ismal,int *ja,
                    int *tapmodel,int *ns, int  *NS,int *colidx,int *rowidx ,double *thres,double *thret)
{
  int i=0,cc=0,j=0,h,k=0,t=0,v=0;
  double dij=0.0,dtv=0.0;


  if (*istap) {// start tapering case

  double thre[2], c_supp[2];
  thre[0]=thres[1];thre[1]=thret[1];
  
  h=0;

 if(isst[0]){  // space time case
        ia[0] = 1;
        for(t=0;t<*ntime;t++){
        for(i=0;i<ncoord[0];i++){
          cc=0;
           for(v=0;v<*ntime;v++){
               dtv=fabs(coordt[t]-coordt[v]);
          for(j=0;j<ncoord[0];j++){
         // dij=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            dij=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],coordz[i],coordz[j],*REARTH);
                // Rprintf("%f %f %f %f %f %f -  %f   %d %d %d %d \n",coordx[i],coordx[j],coordy[i],coordy[j],coordz[i],coordz[j],dij,
                  //  t,i,v,j); 
                Comp_supp(c_supp,tapmodel, dij, dtv,thre);
                  if((dij<c_supp[0]||is_equal(dij,c_supp[0]))&&(dtv<c_supp[1] ||  is_equal(dtv,c_supp[1]))){    
                               tlags[h]=dij;
                               tlagt[h]=dtv;
                               idx[h] =(t * (ncoord[0]) * (ncoord[0]) * ntime[0]) +  (i*  ntime[0] *  *ncoord) +  (1+j+ *ncoord * v);
                               ja[h]=1+j+(*ncoord) * v;
                               cc=cc+1;
                               h = h +1 ;
                }}}
                ia[k+1]=ia[k]+cc;
                k=k+1;}}
   }   // end space time case
}    // end tapering case
  else {   // no tapering
  h=0;
  for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[v];j++){
           dij=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],
                       coordz[(i+NS[t])],coordz[(j+NS[v])],*REARTH);
     
                      if(dij<=thres[1]){
                          tlags[h]=dij;tlagt[h]=0;
                           colidx[h]=i+NS[t];  rowidx[h]=j+NS[v];
                          h++;
                                    }}}
               else {
          dtv=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
           dij=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],
                   coordz[(i+NS[t])],coordz[(j+NS[v])],*REARTH);
                          if(dij<=thres[1] && dtv<=thret[1]){
                   
                            tlags[h]=dij;tlagt[h]=dtv;
                           colidx[h]=i+NS[t];  rowidx[h]=j+NS[v];  
                           h++; 
                             }}}}}}    
}  // end no tapering case

    *npairs=h;
    lags= (double *) R_Calloc(*npairs,double);
    lagt= (double *) R_Calloc(*npairs,double);
    for(i=0;i<*npairs;i++)  
        {lags[i]=tlags[i];lagt[i]=tlagt[i];  }
    // R_Free(tlags); tlags = NULL; R_Free(tlagt); tlagt = NULL;   

  return;
}


/*********************************************************************************************/
/*********************************************************************************************/

// Computes the spatial  distances bivariate case
void SpaceBiv_Dist(double *coordx,double *coordy,double *coordz,double *coordt,int *ia,int *idx,int *ismal,int *ja,
                    int *tapmodel,int *ns, int  *NS,int *colidx,int *rowidx ,double *thres)
{
  int i=0,cc=0,j=0,k=0,t=0,v=0;
  double dij=0.0;
  
      int count=0;

  if (*istap) {// start tapering case


       ia[0] = 1;
   for(t=0;t<*ntime;t++){
   for(i=0;i<ncoord[0];i++){  //for(i=0;i<ns[t];i++){
          cc=0;
           for(v=0;v<*ntime;v++){  
             for(j=0;j<ncoord[0];j++){   //  for(j=0;j<ns[v];j++){
                    //dij=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                      dij=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],coordz[i],coordz[j],*REARTH);
                    if(dij<=dista[t][v]){
                  if(j<i)  {tfirst[count]=v;tsecond[count]=t;}
                  else     {tfirst[count]=t;tsecond[count]=v;}
                               tlags[count]=dij;
                               idx[count] =(t * ncoord[0] * ncoord[0] * ntime[0]) +  (i*  ntime[0] *  ncoord[0]) +  (1+j+ ncoord[0] * v);
                               ja[count]=1+j+ncoord[0] * v;
                               cc=cc+1;
                               count = count +1 ;
                }}}
                ia[k+1]=ia[k]+cc;
                k=k+1;}}

}    // end tapering case

  else {   // no tapering
        count=0;

  for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){    
      if(t==v){
         for(j=i+1;j<ns[v];j++){

          dij=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],
            coordz[(i+NS[t])],coordz[(j+NS[v])],*REARTH);        
                        if(dij<=dista[t][v]){

                           tlags[count]=dij;
                           colidx[count]=i+NS[t];  rowidx[count]=j+NS[v];
                           tfirst[count]=t;tsecond[count]=v;

                           count++;
                                }
        }}
             else {  
         for(j=0;j<ns[v];j++){
         dij=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],
            coordz[(i+NS[t])],coordz[(j+NS[v])],*REARTH);
                        if(dij<=dista[t][v]){

                            tlags[count]=dij;
                           colidx[count]=i+NS[t];  rowidx[count]=j+NS[v];
                           tfirst[count]=t;tsecond[count]=v;
                           count++;
                                 }}}}}}
}  // end no tapering case


  *npairs=count;
  first_1 =(int *)  R_Calloc(count,int);
  second_1=(int *)  R_Calloc(count,int);
  lags_1=(double *) R_Calloc(count,double);
  for(i=0;i<count;i++) { lags_1[i]=tlags[i];
                         first_1[i]=tfirst[i];
                         second_1[i]=tsecond[i];}
  //R_Free(tlags); tlags = NULL;R_Free(tfirst); tfirst = NULL;R_Free(tsecond); tsecond = NULL;
  return;
}






/*********************************************************************************************/
/*********************************************************************************************/
/*********************************************************************************************/

// comp support function:
void Comp_supp(double *c_supp,int *cormod, double h,double u, double *par)
{
    switch(*cormod) // computing comp supports
    {
        
        case 30:// Wendland1
        case 32:// Wendland2
        case 34:// Wendland3
        c_supp[0]=par[0];
        c_supp[1]=-LOW;
        break;
        case 200: // separable spacetime wendland
        case 202:
        case 204:
        case 206:
        case 208:
        case 210:
        case 212:
        case 214:
        case 216:
        c_supp[0]=par[0];
        c_supp[1]=par[1];
        break;
        case  218:   // quasi taper in time
        case  220:
        case  222:
        case 64:
        case 66:
        case 68:
        c_supp[0]=-LOW;
        c_supp[1]=par[1]*pow(1+pow(h/par[0],tapsep[1]/2),-tapsep[0]/(tapsep[1]/2));// arg=pow(1+pow(h/scale_s,power_s/2),-sep/(power_s/2));
        break;
        case  224: // quasi taper in space
        case  226:
        case  228:
        case 63:
        case 65:
        case 67:

        c_supp[0]=tapsep[2]*pow(1+pow(u/tapsep[1],tapsep[1]/2),-tapsep[4]/(tapsep[1]/2));
        c_supp[1]=-LOW;
        break;
        //case  226: // quasi taper in space x quasi taper in time;
        //c_supp[0]=par[0]*pow(1+pow(u/par[1],1),-tapsep[0]);
        //c_supp[1]=par[1]*pow(1+pow(h/par[0],1),-tapsep[1]);
        break;
        case 230:
        case 69:
        case 70:
        case 71:
        case 72:
        case 73:
        case 74:
        case 75:
        case 76:
        case 77:
        case 78:
        case 85: /* fixed rectangular support, independent of sep */
            c_supp[0]=par[0];
            c_supp[1]=par[1];
            break;
            
    }
}




// check if 'val1' is equal to 'val2'  when both are double
int is_equal(double val1, double val2)
{
  return fabs(val1-val2)<MAXERR;
}





//return x*(x-1)*...j  j<x
int fmax_int(int u,int v)
{
  if(u>=v) return u;
  else return(v);
}
/**********************************************/
int fmin_int(int u,int v)
{
  if(u<=v) return u;
  else return(v);
}

// sterling approx



/**********************************************/
/**********************************************/








// Set the global variables for the spatial and spatial-temporal fitting:



void SetGlobalVar(int *biv,double *coordx,double *coordy,double *coordz,double *coordt,int *grid,int *ia,
		  int *idx,int *ismal,int *ja,int *mem, int *nsite,int *nsitex,int *nsitey,int *nsitez,
		  int *npair,double *radius,double *srange, double *sep,int *st, int *times,double *trange,
		  int *tap,int *tapmodel,int *tp,int *weighted, int *colidx,int *rowidx, int *ns, int *NS, int *dyn)
{
  DeleteGlobalVar();



  //Spatial settings: //Spatial settings:
  maxdist=(double *) R_Calloc(3,double);//spatial threshould
  if(maxdist==NULL) {*ismal=0; return;}

  maxtime=(double *) R_Calloc(1,double);//temporal threshold
    if(maxtime==NULL) {*ismal=0; return;}

  ntime=(int *) R_Calloc(1,int);//number of times
    if(ntime==NULL) {*ismal=0; return;}
    *ntime=*times;

/**********************************/
  ncoord=(int *) R_Calloc(1,int);//number of total spatial coordinates
  if(ncoord==NULL) {*ismal=0; return;}
  ncoord[0]=*nsite;

  ncoordx=(int *) R_Calloc(1,int);//number of the first spatial coordinates
  if(ncoordx==NULL) {*ismal=0; return;}
  *ncoordx=*nsitex;

  ncoordy=(int *) R_Calloc(1,int);//number of the second spatial coordinates
  if(ncoordy==NULL) {*ismal=0; return;}
  *ncoordy=*nsitey;

   ncoordz=(int *) R_Calloc(1,int);//number of the third spatial coordinates
  if(ncoordz==NULL) {*ismal=0; return;}
  *ncoordz=*nsitez;


/************************/
  npairs=(int *) R_Calloc(1,int);//effective number of pairs
  //npairs=(int *) R_alloc(1, sizeof(int));
  if(npairs==NULL) {*ismal=0; return;}

  isbiv=(int *) R_Calloc(1,int);//is a bivariate random field?
  if(isbiv==NULL) {*ismal=0; return;}
  isbiv[0]=biv[0];  

  ismem=(int *) R_Calloc(1,int);//is distances computed using memory allocation
  if(ismem==NULL) {*ismal=0; return;}
  *ismem=mem[0];

  isst=(int *) R_Calloc(1,int);//is a spatio-temporal random field?
  if(isst==NULL) {*ismal=0; return;}
  isst[0]=st[0]; 
    
  cdyn=(int *) R_Calloc(1,int);//dynamic coords
  if(cdyn==NULL || dyn==NULL) {*ismal=0; return;}
  cdyn[0]=dyn[0]; 

  istap=(int *) R_Calloc(1,int);//is tapering?
  if(istap==NULL) {*ismal=0; return;}
  istap[0]=tap[0];

  type=(int *) R_Calloc(1,int);//type of distance
  if(type==NULL) {*ismal=0; return;}

  REARTH=(double *) R_Calloc(1,double);//radius of hearth
   if(REARTH==NULL) {*ismal=0; return;}
    REARTH[0]=radius[0];

  //double *REARTH=*radius;//definition
    

  tapsep=(double *) R_Calloc(5,double);
           if(tapsep==NULL){*ismal=0; return;}

  
  *type=*tp;
  


 // case of not saving distances
  if(!ismem[0]) {

                   if(srange[1]) maxdist[0]=srange[1];
                   else maxdist[0]=-LOW;

                   if(isst[0]){ if(trange[1]) maxtime[0]=trange[1];
                                else maxtime[0]=-LOW;
                       }
                   if(isbiv[0]){  
                                   int i=0;
                                   dista= (double **) R_Calloc(ntime[0],double *);
                                   if(dista==NULL) {*ismal=0; return;}
                                   for(i=0;i<ntime[0];i++){
                                   dista[i]=(double *) R_Calloc(ntime[0],double);
                                   if(dista[i]==NULL) {*ismal=0; return;}}
                                  
                                  if(srange[1]) maxdist[0]=srange[1];
                                  else maxdist[0]=-LOW;
                                  if(srange[2]) maxdist[1]=srange[2];
                                  else maxdist[1]=-LOW;
                                  if(srange[3]) maxdist[2]=srange[3];
                                  else maxdist[2]=-LOW;
                                  dista[0][0]=maxdist[0];dista[0][1]=maxdist[1];dista[1][0]=dista[0][1];dista[1][1]=maxdist[2];  
                    }
                  return;    // return here for the not  saving distance  
                } // end not saving distances
  else 
  {     // start  saving distances
  /***********************************************************/  
if(!isst[0]&&!isbiv[0]) {// spatial case

           // settting compact support
            if(srange[1]) maxdist[0]=srange[1];
            else maxdist[0]=-LOW;  

          if(istap[0])  // tapering case
            {
              *npairs=(int)( (*ncoord)  * (*ncoord));
              tlags=(double *) R_Calloc(*npairs,double);
              if(tlags==NULL){*ismal=0; return;}
            }  // end tapering case
          else 
            { // distances composite likelihood  
            *npairs= (int)( 0.5 * (*ncoord) * (*ncoord-1));
            tlags= (double *) R_Calloc(*npairs,double);
            if(tlags==NULL) {*ismal=0; return;}
            } // end  no tapering case


 // computing spatial distances and indexes      
 Space_Dist(coordx,coordy,coordz,ia,idx,ismal,ja,colidx,rowidx,srange[1]);
  //Space_Dist(coordx,coordy,coordz,ia,idx,ismal,ja,colidx,rowidx,srange[1]);
 R_Free(tlags); tlags = NULL;
      if(!ismal[0]) return;
  /***********************************************************/  
}  // end spatial case

else {  //spatio temporal case or bivariate case

  
        int qq;

        if (istap[0])  qq=(*ncoord) * (*ntime);
        else        qq=(*ncoord);
   
    // setting compact supports for space-time and bivariate case
       if(isst[0]){ 
                   if(srange[1]) maxdist[0]=srange[1];
                   else maxdist[0]=-LOW;
                   if(trange[1]) maxtime[0]=trange[1];
                   else maxtime[0]=-LOW; 
        }
       if(isbiv[0])              {
                                   int i=0;
                                   dista= (double **) R_Calloc(ntime[0],double *);
           
                                   if(dista==NULL) {*ismal=0; return;}
                                   for(i=0;i<ntime[0];i++){
                                   dista[i]=(double *) R_Calloc(ntime[0],double);
                                   if(dista[i]==NULL) {*ismal=0; return;}}
           
                                  if(srange[1]) maxdist[0]=srange[1];
                                  else maxdist[0]=-LOW;
                                  if(srange[2]) maxdist[1]=srange[2];
                                  else maxdist[1]=-LOW;
                                  if(srange[3]) maxdist[2]=srange[3];
                                  else maxdist[2]=-LOW;
           
                                  dista[0][0]=maxdist[0];dista[0][1]=maxdist[1];dista[1][0]=dista[0][1];dista[1][1]=maxdist[2]; 
                                  }
     if(istap[0])  // tapering case
        {
          // allocating vectors
           *npairs=(int)( qq * qq );
           tlags=(double *) R_Calloc(*npairs,double);
           if(tlags==NULL){*ismal=0; return;}
         
         if(isst[0])    { //spatio temporal case
                           tlagt=(double *) R_Calloc(*npairs,double);
                           if(tlagt==NULL){*ismal=0; return;}
                           }
         if(isbiv[0])         { //bivariate case
                           tfirst=(int *) R_Calloc(*npairs,int);
                           if(tfirst==NULL){*ismal=0; return;}
                           tsecond=(int *) R_Calloc(*npairs,int);
                           if(tsecond==NULL){*ismal=0; return;}
                              }
           tapsep[0]=sep[0];tapsep[1]=sep[1];tapsep[2]=sep[2];tapsep[3]=sep[3];tapsep[4]=sep[4];
           if(tapsep[0]==1) tapsep[0]=0.99999999;
       }  // end tapering
else {  // distance for composite likelihood
            
    
               if(isst[0])  npairs[0]=(int)(qq * (qq-1) * 0.5);
               if(isbiv[0]) npairs[0]=(int)(qq * (qq-1) * 0.5);
    
               tlags= (double *) R_Calloc(*npairs,double);
    
              if(tlags==NULL) {*ismal=0; return;}
          // allocates the matrix of temporal distances:
          if(isst[0]) {
             //memory allocation of matrix temporal distances
              
             tlagt= (double *) R_Calloc(*npairs,double);
              
             if(tlagt==NULL) {*ismal=0; return;}
                     }
          if(isbiv[0]) {
              
            tfirst=(int *) R_Calloc(*npairs,int);
              
            if(tfirst==NULL){*ismal=0; return;}
              
            tsecond=(int *) R_Calloc(*npairs,int);
              
            if(tsecond==NULL){*ismal=0; return;}
                      }
       }
if(isst[0])  {SpaceTime_Dist(coordx,coordy,coordz,coordt,ia,idx,ismal,ja,tapmodel,
                                        ns,NS,colidx,rowidx,srange,trange);
               R_Free(tlags); tlags = NULL; R_Free(tlagt); tlagt = NULL;}
if(isbiv[0]) {
    
    SpaceBiv_Dist(coordx,coordy,coordz,coordt,ia,idx,ismal,ja,tapmodel,
                                     ns,NS,colidx,rowidx,srange);
    
               R_Free(tlags); tlags = NULL;R_Free(tfirst); tfirst = NULL;R_Free(tsecond); tsecond = NULL;
    
}
  if(!ismal[0]) return;
 } //end spatio temporal case or bivariate case

  /* Bohman weights are fixed by pair geometry and support cutoffs.
     Cache them once rather than recomputing them at every optimizer step. */
  if (pairweights) { R_Free(pairweights); pairweights = NULL; }
  if (weighted && weighted[0] && npairs && npairs[0] > 0) {
      pairweights = (double *)R_Calloc(npairs[0], double);
      if (!pairweights) { *ismal = 0; return; }
      for (int ii = 0; ii < npairs[0]; ++ii) {
          const double hh = isbiv[0] ? lags_1[ii] : lags[ii];
          double ww = CorFunBohman(hh, maxdist[0]);
          if (isst[0]) ww *= CorFunBohman(lagt[ii], maxtime[0]);
          pairweights[ii] = (R_FINITE(ww) && ww > 0.0) ? ww : 0.0;
      }
  }

  npair[0]=npairs[0];
      
  return;
     }
}
/****************************************************************/
/* Both setup paths own the same process-local context.  Release every
 * buffer independently, including partially initialized and bivariate state.
 * The two public cleanup entry points are intentionally idempotent. */
void DeleteGlobalVar(void)
{
#define GEO_FREE_CONTEXT(p) do { if ((p) != NULL) { R_Free(p); (p) = NULL; } } while (0)
    if (dista != NULL) {
        const int rows = ntime != NULL ? ntime[0] : 0;
        for (int i = 0; i < rows; ++i) GEO_FREE_CONTEXT(dista[i]);
        GEO_FREE_CONTEXT(dista);
    }
    GEO_FREE_CONTEXT(lags); GEO_FREE_CONTEXT(lagt);
    GEO_FREE_CONTEXT(lags_1); GEO_FREE_CONTEXT(first_1); GEO_FREE_CONTEXT(second_1);
    GEO_FREE_CONTEXT(first); GEO_FREE_CONTEXT(second);
    GEO_FREE_CONTEXT(tlags); GEO_FREE_CONTEXT(tlagt);
    GEO_FREE_CONTEXT(tfirst); GEO_FREE_CONTEXT(tsecond);
    GEO_FREE_CONTEXT(pairweights);
    GEO_FREE_CONTEXT(pg_count1); GEO_FREE_CONTEXT(pg_count2); GEO_FREE_CONTEXT(pg_count_n);
    GEO_FREE_CONTEXT(maxdist); GEO_FREE_CONTEXT(maxtime);
    GEO_FREE_CONTEXT(ncoord); GEO_FREE_CONTEXT(ncoordx);
    GEO_FREE_CONTEXT(ncoordy); GEO_FREE_CONTEXT(ncoordz);
    GEO_FREE_CONTEXT(npairs); GEO_FREE_CONTEXT(ntime);
    GEO_FREE_CONTEXT(type); GEO_FREE_CONTEXT(REARTH); GEO_FREE_CONTEXT(tapsep);
    GEO_FREE_CONTEXT(isbiv); GEO_FREE_CONTEXT(isst); GEO_FREE_CONTEXT(ismem);
    GEO_FREE_CONTEXT(istap); GEO_FREE_CONTEXT(cdyn);
#undef GEO_FREE_CONTEXT
}
/*#######################################################################*/

void SetGlobalVar2(int *nsite, int *times,          // 2
                   double *h, int *nn, double *maxh,    // 5  
                   double *u, int *tt, double *maxu,    // 8
                   int *st, int *biv, int *one, int *two) // 12
{
    int i = 0;
    
    // Input validation
    if (!nsite || !times || !h || !nn || !maxh || !u || !maxu || 
        !st || !biv || !one || !two) {
        error("Invalid NULL pointer passed to SetGlobalVar2");
        return;
    }
    
    if (*nsite <= 0 || *times <= 0 || *nn <= 0) {
        error("Invalid negative or zero values for nsite, times, or nn");
        return;
    }
    
    /* Reinitialization is allowed; release a previous pair cache first. */
    if (ncoord || ntime || maxdist || maxtime || npairs || isbiv || isst ||
        lags || lags_1 || lagt || pairweights || pg_count1 || pg_count2 ||
        pg_count_n)
        DeleteGlobalVar2();

    // Allocate basic global variables
    ncoord = (int *)R_Calloc(1, int);
    if (!ncoord) {
        error("Memory allocation failed for ncoord");
        return;
    }
    ncoord[0] = *nsite;
    
    ntime = (int *)R_Calloc(1, int);
    if (!ntime) {
        R_Free(ncoord);
        error("Memory allocation failed for ntime");
        return;
    }
    ntime[0] = *times;
    
    maxdist = (double *)R_Calloc(1, double);
    if (!maxdist) {
        R_Free(ncoord);
        R_Free(ntime);
        error("Memory allocation failed for maxdist");
        return;
    }
    maxdist[0] = *maxh;
    
    maxtime = (double *)R_Calloc(1, double);
    if (!maxtime) {
        R_Free(ncoord);
        R_Free(ntime);
        R_Free(maxdist);
        error("Memory allocation failed for maxtime");
        return;
    }
    maxtime[0] = *maxu;
    
    npairs = (int *)R_Calloc(1, int);
    if (!npairs) {
        R_Free(ncoord);
        R_Free(ntime);
        R_Free(maxdist);
        R_Free(maxtime);
        error("Memory allocation failed for npairs");
        return;
    }
    npairs[0] = *nn;
    
    isbiv = (int *)R_Calloc(1, int);
    if (!isbiv) {
        R_Free(ncoord);
        R_Free(ntime);
        R_Free(maxdist);
        R_Free(maxtime);
        R_Free(npairs);
        error("Memory allocation failed for isbiv");
        return;
    }
    isbiv[0] = *biv;
    
    isst = (int *)R_Calloc(1, int);
    if (!isst) {
        R_Free(ncoord);
        R_Free(ntime);
        R_Free(maxdist);
        R_Free(maxtime);
        R_Free(npairs);
        R_Free(isbiv);
        error("Memory allocation failed for isst");
        return;
    }
    isst[0] = *st;

    // Process arrays based on model type
    // Note: bivariate and spatio-temporal are mutually exclusive
    if (!isst[0] && !isbiv[0]) {  
        // Pure spatial case
        lags = (double *)R_Calloc(npairs[0], double);
        if (!lags) {
            // Clean up basic globals on error
            R_Free(ncoord); R_Free(ntime); R_Free(maxdist);
            R_Free(maxtime); R_Free(npairs); R_Free(isbiv); R_Free(isst);
            error("Memory allocation failed for spatial lags");
            return;
        }
        
        for (i = 0; i < npairs[0]; i++) {
            lags[i] = h[i];
        }
    }
    else if (isst[0]) {  
        // Spatio-temporal case (isbiv[0] must be 0)
        lags = (double *)R_Calloc(npairs[0], double);
        lagt = (double *)R_Calloc(npairs[0], double);
        
        if (!lags || !lagt) {
            if (lags) R_Free(lags);
            if (lagt) R_Free(lagt);
            // Clean up basic globals
            R_Free(ncoord); R_Free(ntime); R_Free(maxdist);
            R_Free(maxtime); R_Free(npairs); R_Free(isbiv); R_Free(isst);
            error("Memory allocation failed for spatio-temporal lags");
            return;
        }
        
        for (i = 0; i < npairs[0]; i++) {
            lags[i] = h[i];
            lagt[i] = u[i];
        }
    }
    else if (isbiv[0]) {  
        // Spatial bivariate case (isst[0] must be 0)
        lags_1 = (double *)R_Calloc(npairs[0], double);
        second_1 = (int *)R_Calloc(npairs[0], int);
        first_1 = (int *)R_Calloc(npairs[0], int);
        
        if (!lags_1 || !second_1 || !first_1) {
            if (lags_1) R_Free(lags_1);
            if (second_1) R_Free(second_1);
            if (first_1) R_Free(first_1);
            // Clean up basic globals
            R_Free(ncoord); R_Free(ntime); R_Free(maxdist);
            R_Free(maxtime); R_Free(npairs); R_Free(isbiv); R_Free(isst);
            error("Memory allocation failed for bivariate arrays");
            return;
        }
        
        for (i = 0; i < npairs[0]; i++) {
            lags_1[i] = h[i];
            first_1[i] = one[i];
            second_1[i] = two[i];
        }
    }

    /* Pair weights are initialized lazily by GeoCachedPairWeight().
       This avoids allocating an additional O(npairs) vector for unweighted
       fits while still computing every fixed Bohman weight only once. */
    pairweights = NULL;

    return;
}

/*#######################################################################*/
/* Release the pair cache.  Free every pointer independently so this is */
/* also safe after a partial allocation failure.                         */
/*#######################################################################*/
void DeleteGlobalVar2(void)
{
    DeleteGlobalVar();
}

void SetPGPairCounts(double *data1, double *data2, int *n, int *ok)
{
    if (ok) ok[0] = 0;
    if (!data1 || !data2 || !n || n[0] < 0) return;

    if (pg_count1) { R_Free(pg_count1); pg_count1 = NULL; }
    if (pg_count2) { R_Free(pg_count2); pg_count2 = NULL; }
    if (pg_count_n) { R_Free(pg_count_n); pg_count_n = NULL; }

    pg_count_n = (int *)R_Calloc(1, int);
    if (!pg_count_n) return;
    pg_count_n[0] = n[0];
    if (n[0] == 0) {
        if (ok) ok[0] = 1;
        return;
    }

    pg_count1 = (int *)R_Calloc(n[0], int);
    pg_count2 = (int *)R_Calloc(n[0], int);
    if (!pg_count1 || !pg_count2) goto fail;

    for (int i = 0; i < n[0]; ++i) {
        if (ISNAN(data1[i]) || ISNAN(data2[i])) {
            pg_count1[i] = INT_MIN;
            pg_count2[i] = INT_MIN;
            continue;
        }
        if (!R_FINITE(data1[i]) || !R_FINITE(data2[i])) goto fail;
        const double r1 = nearbyint(data1[i]);
        const double r2 = nearbyint(data2[i]);
        if (r1 < 0.0 || r2 < 0.0 || r1 > INT_MAX || r2 > INT_MAX ||
            fabs(data1[i] - r1) > 1e-10 ||
            fabs(data2[i] - r2) > 1e-10)
            goto fail;
        pg_count1[i] = (int)r1;
        pg_count2[i] = (int)r2;
    }
    if (ok) ok[0] = 1;
    return;

fail:
    if (pg_count1) { R_Free(pg_count1); pg_count1 = NULL; }
    if (pg_count2) { R_Free(pg_count2); pg_count2 = NULL; }
    if (pg_count_n) { R_Free(pg_count_n); pg_count_n = NULL; }
}

double GeoCachedPairWeight(int index, int weighted, int spacetime)
{
    if (!weighted) return 1.0;
    if (index < 0 || !npairs || index >= npairs[0] || npairs[0] <= 0)
        return R_NaN;

    if (!pairweights) {
        const int use_st = spacetime || (isst && isst[0]);
        const double *hvec = (isbiv && isbiv[0]) ? lags_1 : lags;
        if (!hvec || !maxdist || maxdist[0] <= 0.0 ||
            (use_st && (!lagt || !maxtime || maxtime[0] <= 0.0)))
            return R_NaN;

        pairweights = (double *)R_Calloc(npairs[0], double);
        if (!pairweights) return R_NaN;
        for (int i = 0; i < npairs[0]; ++i) {
            double out = CorFunBohman(hvec[i], maxdist[0]);
            if (use_st) out *= CorFunBohman(lagt[i], maxtime[0]);
            pairweights[i] = (R_FINITE(out) && out > 0.0) ? out : 0.0;
        }
    }
    return pairweights[index];
}

/*#######################################################################*/
/* Funzione helper per controllo stato delle variabili globali */
/*#######################################################################*/






// ============= TEST QQNORM

#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)		\
if (log_p) {					\
if(p > 0)					\
return NAN;				\
if(p == 0) /* upper bound*/			\
return lower_tail ? _RIGHT_ : _LEFT_;	\
if(p == -INFINITY)				\
return lower_tail ? _LEFT_ : _RIGHT_;	\
}							\
else { /* !log_p */					\
if(p < 0 || p > 1)				\
return NAN;				\
if(p == 0)					\
return lower_tail ? _LEFT_ : _RIGHT_;	\
if(p == 1)					\
return lower_tail ? _RIGHT_ : _LEFT_;	\
}

#define R_D_Cval(p)	(lower_tail ? (0.5 - (p) + 0.5) : (p))	/*  1 - p */

#define R_DT_CIv(p)	(log_p ? (lower_tail ? -expm1(p) : exp(p)) \
: R_D_Cval(p))

#define R_D_Lval(p)	(lower_tail ? (p) : (0.5 - (p) + 0.5))	/*  p  */


#define R_DT_qIv(p)	(log_p ? (lower_tail ? exp(p) : - expm1(p)) \
: R_D_Lval(p))


double qnorm55(double p, double mu, double sigma, int lower_tail, int log_p)
{
    double p_, q, r, val;
    
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma))
        return p + mu + sigma;
#endif
    R_Q_P01_boundaries(p, -INFINITY, INFINITY);
    
    if(sigma  < 0)	return NAN;
    if(sigma == 0)	return mu;
    
    p_ = R_DT_qIv(p);/* real lower_tail prob. p */
    q = p_ - 0.5;
    
#ifdef DEBUG_qnorm
    REprintf("qnorm(p=%10.7g, m=%g, s=%g, l.t.= %d, log= %d): q = %g\n",
             p,mu,sigma, lower_tail, log_p, q);
#endif
    
    
    /*-- use AS 241 --- */
    /* double ppnd16_(double *p, long *ifault)*/
    /*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
     
     Produces the normal deviate Z corresponding to a given lower
     tail area of P; Z is accurate to about 1 part in 10**16.
     
     (original fortran code used PARAMETER(..) for the coefficients
     and provided hash codes for checking them...)
     */
    if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
        r = .180625 - q * q;
        val =
        q * (((((((r * 2509.0809287301226727 +
                   33430.575583588128105) * r + 67265.770927008700853) * r +
                 45921.953931549871457) * r + 13731.693765509461125) * r +
               1971.5909503065514427) * r + 133.14166789178437745) * r +
             3.387132872796366608)
        / (((((((r * 5226.495278852854561 +
                 28729.085735721942674) * r + 39307.89580009271061) * r +
               21213.794301586595867) * r + 5394.1960214247511077) * r +
             687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    }
    else { /* closer than 0.075 from {0,1} boundary */
        
        /* r = min(p, 1-p) < 0.075 */
        if (q > 0)
            r = R_DT_CIv(p);/* 1-p */
        else
            r = p_;/* = R_DT_Iv(p) ^=  p */
        
        r = sqrt(- ((log_p &&
                     ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ?
                    p : /* else */ log(r)));
        /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */
#ifdef DEBUG_qnorm
        REprintf("\t close to 0 or 1: r = %7g\n", r);
#endif
        
        if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                         .0227238449892691845833) * r + .24178072517745061177) *
                       r + 1.27045825245236838258) * r +
                      3.64784832476320460504) * r + 5.7694972214606914055) *
                    r + 4.6303378461565452959) * r +
                   1.42343711074968357734)
            / (((((((r *
                     1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                    r + .0151986665636164571966) * r +
                   .14810397642748007459) * r + .68976733498510000455) *
                 r + 1.6763848301838038494) * r +
                2.05319162663775882187) * r + 1.);
        }
        else { /* very close to  0 or 1 */
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                         2.71155556874348757815e-5) * r +
                        .0012426609473880784386) * r + .026532189526576123093) *
                      r + .29656057182850489123) * r +
                     1.7848265399172913358) * r + 5.4637849111641143699) *
                   r + 6.6579046435011037772)
            / (((((((r *
                     2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                    r + 1.8463183175100546818e-5) * r +
                   7.868691311456132591e-4) * r + .0148753612908506148525)
                 * r + .13692988092273580531) * r +
                .59983220655588793769) * r + 1.);
        }
        
        if(q < 0.0)
            val = -val;
        /* return (q >= 0.)? r : -r ;*/
    }
    return mu + sigma * val;
}






void qnorm55_call(double *p, double *mu, double *sigma, int *lower_tail, int *log_p, double *res)
{
    *res=   qnorm55( *p,  *mu,  *sigma,  *lower_tail,  *log_p);
    
}
