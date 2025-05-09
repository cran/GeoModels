#include "header.h"


static int find_bin(double* bins, int nbins, double lags) {
    int left = 0, right = nbins - 2; //
    while (left <= right) {
        int mid = (left + right) / 2;
        if (lags >= bins[mid] && lags < bins[mid + 1]) {
            return mid;
        } else if (lags < bins[mid]) {
            right = mid - 1;
        } else {
            left = mid + 1;
        }
    }
    return -1; // not found
}



// Calculates the pairs of coordinates and stores the distances within the given bins:
void pairs(int *ncoords, double *data, double *coordx, double *coordy, double *coordz,
           double *numbins, double *bins, double *v0, double *v1, double *v2,
           double *maxdist, int *typ, double *radius) 
{
    int ncrd = *ncoords;
    int numbin = (int)(*numbins);
    double max_dist = *maxdist;
    int k = 0;

    for (int i = 0; i < ncrd - 1; i++) {
        double xi = coordx[i];
        double yi = coordy[i];
        double zi = coordz[i];
        double di = data[i];

        for (int j = i + 1; j < ncrd; j++) {
            double distance = dist(typ[0], xi, coordx[j], yi, coordy[j], zi, coordz[j], *radius);

            if (distance <= max_dist) {
                int bin_idx = find_bin(bins, numbin - 1, distance);
                if (bin_idx != -1) {
                    v0[k] = bins[bin_idx];
                    v1[k] = di;
                    v2[k] = data[j];
                    k++;
    }}}}}

// binned spatial variogram:
void Binned_Variogram2(double *bins, double *coordx, double *coordy, double *coordz, double *coordt,
                       double *data, int *lbins, double *moms, int *nbins)
{
    int npoints = ncoord[0]; // Assumo ncoord è globale o definito altrove
    double mm[2];
    Maxima_Minima_dist(mm, coordx, coordy, coordz, ncoord, type, REARTH);
    if (maxdist[0] < mm[1]) mm[1] = maxdist[0];

    double step = (mm[1] - mm[0]) / (*nbins - 1);

    bins[0] = mm[0];
    for (int h = 1; h < *nbins; h++) {
        bins[h] = bins[h - 1] + step;
    }
    for (int h = 0; h < *nbins - 1; h++) {
        moms[h] = 0.0;
        lbins[h] = 0;
    }
    for (int i = 0; i < npoints; i++) {
        double xi = coordx[i], yi = coordy[i], zi = coordz[i];
        double xi_data = data[i]; 
        if (ISNAN(xi_data)) continue;
        for (int j = i + 1; j < npoints; j++) {
            double lags = dist(type[0], xi, coordx[j], yi, coordy[j], zi, coordz[j], *REARTH);
            if (lags > *maxdist) continue;
            int bin_idx = find_bin(bins, *nbins, lags);
            if (bin_idx == -1) continue;
            double xj_data = data[j];
            if (ISNAN(xj_data)) continue;

            double diff = xi_data - xj_data;
            moms[bin_idx] += 0.5 * diff * diff;
            lbins[bin_idx]++;
}}}
/***********************************************************************************************************************************/
// binned spatial variogram:
void Binned_Variogram2new(double *bins, int *np, double *data1, double *data2, 
                          double *vdist, int *lbins, double *moms, int *nbins, double *mm)
{
    int n = *np;
    int nbin = *nbins;

    // Calcola step e bins
    double step = (mm[1] - mm[0]) / (nbin - 1);
    bins[0] = mm[0];
    for (int h = 1; h < nbin; h++) {
        bins[h] = bins[h - 1] + step;
    }

    // Azzeramento momenti e conteggi (opzionale, se non già fatto)
    for (int h = 0; h < nbin - 1; h++) {
        moms[h] = 0.0;
        lbins[h] = 0;
    }

    // Loop principale con ricerca binaria
    for (int k = 0; k < n; k++) {
        double dist_val = vdist[k];
        int bin_idx = find_bin(bins, nbin, dist_val);
        if (bin_idx == -1) continue;  // fuori range

        double x = data1[k];
        double y = data2[k];
        if (ISNAN(x) || ISNAN(y)) continue;

        double diff = x - y;
        moms[bin_idx] += 0.5 * diff * diff;
        lbins[bin_idx]++;
    }
}
/***********************************************************************************************************************************/
void Binned_Variogram_22(double *bins, double *coordx, double *coordy,double *coordz, double *coordt,double *data, int *lbins, double *moms, int *nbins)
{
  int h=0, i=0, j=0, n=0,p=0;
  double x,y,step=0.0,*mm;
  //Set the binnes step:
  mm=(double *) R_alloc(2, sizeof(double));
  Maxima_Minima_dist(mm, coordx, coordy,coordz, ncoord,type,REARTH);
  if(maxdist[0]<mm[1]) mm[1]=maxdist[0];
  step=(mm[1]-mm[0])/(*nbins-1);
  bins[0]= mm[0];
  //define bins:
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //Computes the binned moments:
  for(i=0;i<(ncoord[0]-1);i++){
    for(j=(i+1);j<ncoord[0];j++){
      if(lags[p]<=*maxdist){
  for(h=0;h<(*nbins-1);h++)
    if((bins[h]<=lags[p]) && (lags[p]<bins[h+1])){
            x=data[n+i ]; y=data[n+j ];
            if(!(ISNAN(x)||ISNAN(y))){
        moms[h]+=0.5*pow(x-y,2);
        lbins[h]+=1;}}
        p++;}}}
  return;
}




// binned spatial-temporal variogram:
void Binned_Variogram_st2(double *bins, double *bint, double *coordx, double *coordy, double *coordz, double *coordt,
                         double *data, int *lbins, int *lbinst, int *lbint,
                         double *moms, double *momst, double *momt,
                         int *nbins, int *nbint, int *ns, int *NS)
{
  int h, i, j, q, t, u, v;
  double lags = 0.0, lagt = 0.0, step = 0.0, x, y, diff2;
  int nbins1 = *nbins - 1;
  int nbint1 = *nbint - 1;

  double *mm;
  mm = (double *) R_alloc(2, sizeof(double));
  Maxima_Minima_dist(mm, coordx, coordy, coordz, ncoord, type, REARTH);

  if (maxdist[0] < mm[1]) mm[1] = maxdist[0];

  step = mm[1] / nbins1;
  bins[0] = mm[0];
  for (h = 1; h < *nbins; h++)
    bins[h] = bins[h - 1] + step;

  for (t = 0; t < ntime[0]; t++) {
    for (i = 0; i < ns[t]; i++) {
      for (v = t; v < ntime[0]; v++) {
        if (t == v) { // marginal spatial variogram
          for (j = i + 1; j < ns[v]; j++) {
            lags = dist(type[0], coordx[i + NS[t]], coordx[j + NS[v]],
                        coordy[i + NS[t]], coordy[j + NS[v]],
                        coordz[i + NS[t]], coordz[j + NS[v]], *REARTH);
            if (lags <= *maxdist) {
              x = data[i + NS[t]];
              y = data[j + NS[v]];
              if (!(ISNAN(x) || ISNAN(y))) {
                diff2 = 0.5 * (x - y) * (x - y);
                for (h = 0; h < nbins1; h++) {
                  if (bins[h] <= lags && lags < bins[h + 1]) {
                    moms[h] += diff2;
                    lbins[h]++;
                    break; // bin found, exit loop
                  }
                }
              }
            }
          }
        } else {
          lagt = fabs(coordt[t] - coordt[v]);
          for (j = 0; j < ns[v]; j++) {
            x = data[i + NS[t]];
            y = data[j + NS[v]];
            if (i == j) { // marginal temporal variogram
              if (lagt <= *maxtime && !(ISNAN(x) || ISNAN(y))) {
                diff2 = 0.5 * (x - y) * (x - y);
                for (u = 0; u < nbint1; u++) {
                  if (bint[u] <= lagt && lagt < bint[u + 1]) {
                    momt[u] += diff2;
                    lbint[u]++;
                    break;
                  }
                }
              }
            } else { // spatial-temporal variogram
              lags = dist(type[0], coordx[i + NS[t]], coordx[j + NS[v]],
                          coordy[i + NS[t]], coordy[j + NS[v]],
                          coordz[i + NS[t]], coordz[j + NS[v]], *REARTH);
              if (lags <= *maxdist && lagt <= *maxtime && !(ISNAN(x) || ISNAN(y))) {
                diff2 = 0.5 * (x - y) * (x - y);
                q = 0;
                for (h = 0; h < nbins1; h++) {
                  for (u = 0; u < nbint1; u++) {
                    if (bins[h] <= lags && lags < bins[h + 1] &&
                        bint[u] <= lagt && lagt < bint[u + 1]) {
                      momst[q] += diff2;
                      lbinst[q]++;
                      break;
                    }
                    q++;
                  }}}}}}}}}}


// binned spatial-temporal variogram:
void Binned_Variogram_st2_dyn(
    double *bins, double *bint, double *coordx, double *coordy, double *coordz, double *coordt, double *data,
    int *lbins, int *lbinst, int *lbint, double *moms, double *momst, double *momt,
    int *nbins, int *nbint, int *ns, int *NS
) {
    int h, i, j, q, t, u, v;
    double x, y, lags = 0.0, lagt = 0.0, step = 0.0, diff2;
    int nbins1 = *nbins - 1;
    int nbint1 = *nbint - 1;

    // defines the spatial bins:
    double *mm = (double *) R_alloc(2, sizeof(double));
    Maxima_Minima_dist(mm, coordx, coordy, coordz, ncoord, type, REARTH); // computing max and min distances
    if (maxdist[0] < mm[1]) mm[1] = maxdist[0];

    // Set the bin step:
    step = mm[1] / nbins1;
    bins[0] = mm[0];
    for (h = 1; h < *nbins; h++)
        bins[h] = bins[h - 1] + step;

    for (t = 0; t < ntime[0]; t++) {
        for (i = 0; i < ns[t]; i++) {
            for (v = t; v < ntime[0]; v++) {
                if (t == v) { // marginal spatial variogram
                    for (j = i + 1; j < ns[v]; j++) {
                        lags = dist(type[0], coordx[i + NS[t]], coordx[j + NS[v]],
                                    coordy[i + NS[t]], coordy[j + NS[v]],
                                    coordz[i + NS[t]], coordz[j + NS[v]], *REARTH);
                        if (lags <= *maxdist) {
                            x = data[i + NS[t]];
                            y = data[j + NS[v]];
                            if (!(ISNAN(x) || ISNAN(y))) {
                                diff2 = 0.5 * (x - y) * (x - y);
                                for (h = 0; h < nbins1; h++) {
                                    if (bins[h] <= lags && lags < bins[h + 1]) {
                                        moms[h] += diff2;
                                        lbins[h]++;
                                        break; // bin found, exit loop
                }}}}}
                } else {
                    lagt = fabs(coordt[t] - coordt[v]);
                    for (j = 0; j < ns[v]; j++) {
                        lags = dist(type[0], coordx[i + NS[t]], coordx[j + NS[v]],
                                    coordy[i + NS[t]], coordy[j + NS[v]],
                                    coordz[i + NS[t]], coordz[j + NS[v]], *REARTH);

                        // marginal temporal variogram (special bin)
                        if ((bins[0] / 2 <= lags) && (lags < bins[1] / 2) && lagt <= *maxtime) {
                            x = data[i + NS[t]];
                            y = data[j + NS[v]];
                            if (!(ISNAN(x) || ISNAN(y))) {
                                diff2 = 0.5 * (x - y) * (x - y);
                                for (u = 0; u < nbint1; u++) {
                                    if (bint[u] <= lagt && lagt < bint[u + 1]) {
                                        momt[u] += diff2;
                                        lbint[u]++;
                                        break;
                                    }
                                }
                            }
                        }

                        // spatial-temporal variogram
                        if (lags <= *maxdist && lagt <= *maxtime) {
                            x = data[i + NS[t]];
                            y = data[j + NS[v]];
                            if (!(ISNAN(x) || ISNAN(y))) {
                                diff2 = 0.5 * (x - y) * (x - y);
                                q = 0;
                                for (h = 0; h < nbins1; h++) {
                                    for (u = 0; u < nbint1; u++) {
                                        if (bins[h] <= lags && lags < bins[h + 1] &&
                                            bint[u] <= lagt && lagt < bint[u + 1]) {
                                            momst[q] += diff2;
                                            lbinst[q]++;
                                            break;
                                        }
                                        q++;
}}}}}}}}}}





void Binned_Variogram_biv2new(double *bins, int *np,double *data1, double *data2,  double *vdist, double *mm,
     double *moms00,double *moms10,double *moms11,
       int *lbins00,int *lbins10,int *lbins11,
     int *nbins,int *first, int *second)
{
  int h=0,  k=0;
  double step=0.0; 
  double a=0.0,b=0.0,c=0.0,d=0.0;
  step=(mm[1]-mm[0])/(*nbins-1);
  bins[0]= mm[0];
  //define bins:
  for(h=1;h<*nbins;h++) bins[h]=bins[h-1]+step;

            for(k=0;k<*np;k++){  
                  for(h=0;h<(*nbins-1);h++){
                 if((bins[h]<=vdist[k]) && (vdist[k]<bins[h+1])){ 
    
if(!(ISNAN(data1[k])||ISNAN(data2[k]))) {     

if(!first[k]) a=  data1[k]-data2[k];
if(first[k])  b=  data1[k]-data2[k];
if(!second[k])c=  data1[k]-data2[k];
if(second[k]) d=  data1[k]-data2[k];

                  
                if(!first[k]&&!second[k]) { moms00[h]+=0.5*a*c; lbins00[h]+=1;}

                
               if(first[k]&&!second[k])  { moms10[h]+=0.5*b*c;lbins10[h]+=1;}
                if(!first[k]&&second[k])  { moms10[h]+=0.5*a*d;lbins10[h]+=1;}

             
                 if(first[k]&&second[k])      {moms11[h]+=0.5*b*d; lbins11[h]+=1;}
     }}}}
  return;

}


void Binned_Variogram_biv2(
    double *bins, double *coordx, double *coordy, double *coordz, double *coordt, double *data,
    int *cross_lbins, double *cross_moms, int *nbins,
    int *marg_lbins, double *marg_moms, int *ns, int *NS)
{
    int h, i, j, t, v;
    double x, y, a, b, lags = 0.0, step = 0.0, *mm, md;

    mm = (double *) R_alloc(2, sizeof(double));
    Maxima_Minima_dist(mm, coordx, coordy, coordz, ncoord, type, REARTH);
    md = fmax(dista[0][1], fmax(dista[1][1], dista[0][0]));
    if (md < mm[1]) mm[1] = md;

    if (*nbins < 2) return;

    step = mm[1] / (*nbins - 1);
    bins[0] = 0.0;
    for (h = 1; h < *nbins; h++) {
        bins[h] = bins[h - 1] + step;
    }

    for (t = 0; t < ntime[0]; t++) {
        for (i = 0; i < ns[t]; i++) {
            int idx_i_t = i + NS[t];

            for (v = t; v < ntime[0]; v++) {
                if (t == v) {
                    for (j = i + 1; j < ns[t]; j++) {
                        int idx_j_t = j + NS[t];
                        lags = dist(type[0],
                                    coordx[idx_i_t], coordx[idx_j_t],
                                    coordy[idx_i_t], coordy[idx_j_t],
                                    coordz[idx_i_t], coordz[idx_j_t],
                                    *REARTH);

                        if (lags <= dista[t][v]) {
                            for (h = 0; h < *nbins - 1; h++) {
                                if (lags >= bins[h] && lags < bins[h + 1]) {
                                    x = data[idx_i_t];
                                    y = data[idx_j_t];
                                    if (!(ISNAN(x) || ISNAN(y))) {
                                        double diff = x - y;
                                        marg_moms[h + t * (*nbins - 1)] += 0.5 * diff * diff;
                                        marg_lbins[h + t * (*nbins - 1)] += 1;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                } else {
                    for (j = 0; j < ns[v]; j++) {
                        int idx_j_v = j + NS[v];
                        lags = dist(type[0],
                                    coordx[idx_i_t], coordx[idx_j_v],
                                    coordy[idx_i_t], coordy[idx_j_v],
                                    coordz[idx_i_t], coordz[idx_j_v],
                                    *REARTH);

                        if (lags <= dista[t][v]) {
                            for (h = 0; h < *nbins - 1; h++) {
                                if (lags >= bins[h] && lags < bins[h + 1]) {
                                    // corretti: due punti da campo t, due da campo v
                                    x = data[idx_i_t];       // campo t
                                    y = data[j + NS[t]];     // campo t
                                    a = data[idx_i_t + (NS[v] - NS[t])]; // campo v
                                    b = data[idx_j_v];       // campo v

                                    if (!(ISNAN(x) || ISNAN(y) || ISNAN(a) || ISNAN(b))) {
                                        cross_moms[h + (v - t - 1) * (*nbins - 1)] += 0.5 * (x - y) * (a - b);
                                        cross_lbins[h + (v - t - 1) * (*nbins - 1)] += 1;
                                    }
                                    break;
                                }}}}}}}}}


/***********************************************************************************************************************************/
// variogram cloud:
void Cloud_Variogram2(double *bins, double *coordx, double *coordy,double *coordz, double *coordt,double *data, int *lbins, double *moms, int *nbins)
{
  int  h=0,i=0, j=0, n=0;double lags=0.0,x,y;
 //Computes the cloud moments:
  for(i=0;i<(ncoord[0]-1);i++)
    for(j=(i+1);j<ncoord[0];j++){
          dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],coordz[i],coordz[j],*REARTH);
      bins[h]=lags;
        x=data[n+i ];  y=data[n+j ];
        if(!(ISNAN(x)||ISNAN(y))){
	        moms[h]+=0.5*pow(x-y,2);
        lbins[h]=1;
        h++;}}
  return;
}
/***********************************************************************************************************************************/
/***********************************************************************************************************************************/

// Least square method for Gaussian spatial-temporal random field:
void LeastSquare_G(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		      int *nbins, int *nbint, double *nuis, double *par, double *res)
{
  int h=0, i=0, u=0;
  double vario=0.0, varhat=0.0;
  //Checks the nuisance parameters (nugget, sill and corr):
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2) {
    *res=LOW; return;}
  // Computes the least squares:
  for(u=0;u<*nbint;u++)
    for(h=0;h<(*nbins-1);h++){
      vario=moms[i]/lbins[i];// Computes the empirical variogram
      varhat=Variogram(cormod,0.5*(bins[h]+bins[h+1]),bint[u],nuis[1],nuis[2],par); 
      *res=*res-pow(varhat-vario,2);// Computes the least squares
      i++;}
  return;
}
// Weighted least square method for Gaussian spatial-temporal random field:
void WLeastSquare_G(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		       int *nbins, int *nbint, double *nuis, double *par, double *res)
{
  int h=0,i=0,u=0;
  double vario=0.0,varhat=0.0;
  //Checks the nuisance parameters (nugget, sill and corr):
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2) {
    *res=LOW; return;}
  // Computes the least squares:
  for(u=0;u<*nbint;u++)
    for(h=0;h<(*nbins-1);h++){
      vario=moms[i]/lbins[i];// Computes the empirical variogram
      varhat=Variogram(cormod,0.5*(bins[h]+bins[h+1]),bint[u],nuis[1],nuis[2],par);
      if(vario) *res=*res-pow(varhat-vario,2)*(lbins[i]/pow(vario,2));// Computes the weighted least squares
      i++;}
  return;
}
