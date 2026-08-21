#include "header.h"


static int find_bin(double *bins, int nbins, double lags)
{
    int left, right;

    if (nbins < 2 || !R_FINITE(lags)) return -1;
    if (lags < bins[0] || lags > bins[nbins - 1]) return -1;

    /* The final bin is closed on the right so that the maximum lag is kept. */
    if (lags == bins[nbins - 1]) return nbins - 2;

    left = 0;
    right = nbins - 2;
    while (left <= right) {
        int mid = (left + right) / 2;
        if (lags >= bins[mid] && lags < bins[mid + 1]) {
            return mid;
        }
        if (lags < bins[mid]) right = mid - 1;
        else left = mid + 1;
    }
    return -1;
}



void pairs(int *ncoords, double *data, double *coordx, double *coordy, double *coordz,
           double *numbins, double *bins, double *v0, double *v1, double *v2,
           double *maxdist, int *typ, double *radius)
{
    int ncrd = *ncoords;
    int numbin = (int)(*numbins);
    double max_dist = *maxdist;
    double b0, step;
    int k = 0;
    if (ncrd < 2 || numbin < 1 || max_dist < 0.0) return;

    b0 = bins[0];
    step = bins[1] - bins[0];
    if (step <= 0.0) return;

    for (int i = 0; i < ncrd - 1; i++) {
        double xi = coordx[i];
        double yi = coordy[i];
        double zi = coordz[i];
        double di = data[i];
        if (!R_FINITE(di)) continue;
        for (int j = i + 1; j < ncrd; j++) {
            double dj = data[j];
            double distance;
            int bin_idx;
            if (!R_FINITE(dj)) continue;
            distance = dist(typ[0], xi, coordx[j],
                            yi, coordy[j],
                            zi, coordz[j],
                            *radius);
            if (distance > max_dist) continue;
            bin_idx = (int)((distance - b0) / step);
            /* Include exactly max_dist in the last bin. */
            if (bin_idx == numbin && distance == max_dist) bin_idx = numbin - 1;
            if (bin_idx >= 0 && bin_idx < numbin &&
                distance >= bins[bin_idx] && distance <= bins[bin_idx + 1]) {
                v0[k] = bins[bin_idx];
                v1[k] = di;
                v2[k] = dj;
                k++;
            }
        }
    }
}

// binned spatial variogram:
void Binned_Variogram2(double *bins, double *coordx, double *coordy, double *coordz, double *coordt,
                       double *data, int *lbins, double *moms, int *nbins)
{
    int npoints = ncoord[0];
    double mm[2], upper, step;

    if (*nbins < 2) return;
    for (int h = 0; h < *nbins - 1; h++) {
        moms[h] = 0.0;
        lbins[h] = 0;
    }

    Maxima_Minima_dist(mm, coordx, coordy, coordz, ncoord, type, REARTH);
    upper = mm[1];
    if (maxdist[0] < upper) upper = maxdist[0];

    /* If the cutoff lies below the smallest observed positive distance,
       return a valid empty binning rather than a decreasing grid. */
    if (!R_FINITE(mm[0]) || !R_FINITE(upper) || upper < mm[0]) {
        double hi = (R_FINITE(maxdist[0]) && maxdist[0] > 0.0) ? maxdist[0] : 1.0;
        step = hi / (*nbins - 1);
        bins[0] = 0.0;
        for (int h = 1; h < *nbins; h++) bins[h] = bins[h - 1] + step;
        bins[*nbins - 1] = hi;
        return;
    }

    if (upper <= mm[0]) {
        double delta = sqrt(DBL_EPSILON) * fmax(1.0, fabs(mm[0]));
        bins[0] = fmax(0.0, mm[0] - delta);
        upper = mm[0] + delta;
    } else {
        bins[0] = mm[0];
    }

    step = (upper - bins[0]) / (*nbins - 1);
    for (int h = 1; h < *nbins; h++) bins[h] = bins[h - 1] + step;
    bins[*nbins - 1] = upper;

    for (int i = 0; i < npoints; i++) {
        double xi = coordx[i], yi = coordy[i], zi = coordz[i];
        double xi_data = data[i];
        if (!R_FINITE(xi_data)) continue;
        for (int j = i + 1; j < npoints; j++) {
            double lags = dist(type[0], xi, coordx[j], yi, coordy[j], zi, coordz[j], *REARTH);
            int bin_idx;
            double xj_data, diff;
            if (!R_FINITE(lags) || lags > *maxdist) continue;
            bin_idx = find_bin(bins, *nbins, lags);
            if (bin_idx < 0) continue;
            xj_data = data[j];
            if (!R_FINITE(xj_data)) continue;
            diff = xi_data - xj_data;
            moms[bin_idx] += 0.5 * diff * diff;
            lbins[bin_idx]++;
        }
    }
}
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
    bins[nbin - 1] = mm[1];

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
        if (!R_FINITE(x) || !R_FINITE(y)) continue;

        double diff = x - y;
        moms[bin_idx] += 0.5 * diff * diff;
        lbins[bin_idx]++;
    }
}


// binned spatial-temporal variogram:
void Binned_Variogram_st2(double *bins, double *bint, double *coordx, double *coordy, double *coordz, double *coordt,
                         double *data, int *lbins, int *lbinst, int *lbint,
                         double *moms, double *momst, double *momt,
                         int *nbins, int *nbint, int *ns, int *NS)
{
    int i, j, t, v;
    int nbins1 = *nbins - 1;
    int nbint1 = *nbint - 1;
    double mm[2];
    double step;

    if (nbins1 <= 0 || nbint1 <= 0) return;

    for (int h = 0; h < nbins1; h++) { moms[h] = 0.0; lbins[h] = 0; }
    for (int u = 0; u < nbint1; u++) { momt[u] = 0.0; lbint[u] = 0; }
    for (int q = 0; q < nbins1 * nbint1; q++) { momst[q] = 0.0; lbinst[q] = 0; }

    Maxima_Minima_dist(mm, coordx, coordy, coordz, ncoord, type, REARTH);
    {
        double upper = mm[1];
        if (maxdist[0] < upper) upper = maxdist[0];
        if (!R_FINITE(mm[0]) || !R_FINITE(upper) || upper < mm[0]) {
            double hi = (R_FINITE(maxdist[0]) && maxdist[0] > 0.0) ? maxdist[0] : 1.0;
            step = hi / nbins1;
            bins[0] = 0.0;
            for (int h = 1; h < *nbins; h++) bins[h] = bins[h - 1] + step;
            bins[*nbins - 1] = hi;
        } else {
            if (upper <= mm[0]) {
                double delta = sqrt(DBL_EPSILON) * fmax(1.0, fabs(mm[0]));
                bins[0] = fmax(0.0, mm[0] - delta);
                upper = mm[0] + delta;
            } else bins[0] = mm[0];
            step = (upper - bins[0]) / nbins1;
            for (int h = 1; h < *nbins; h++) bins[h] = bins[h - 1] + step;
            bins[*nbins - 1] = upper;
        }
    }

    for (t = 0; t < ntime[0]; t++) {
        for (i = 0; i < ns[t]; i++) {
            for (v = t; v < ntime[0]; v++) {
                if (t == v) {
                    for (j = i + 1; j < ns[v]; j++) {
                        double x = data[i + NS[t]];
                        double y = data[j + NS[v]];
                        double lags;
                        int hb;

                        if (!R_FINITE(x) || !R_FINITE(y)) continue;
                        lags = dist(type[0], coordx[i + NS[t]], coordx[j + NS[v]],
                                    coordy[i + NS[t]], coordy[j + NS[v]],
                                    coordz[i + NS[t]], coordz[j + NS[v]], *REARTH);
                        if (lags > *maxdist) continue;
                        hb = find_bin(bins, *nbins, lags);
                        if (hb < 0) continue;
                        {
                            double d = x - y;
                            moms[hb] += 0.5 * d * d;
                            lbins[hb]++;
                        }
                    }
                } else {
                    double lagt = fabs(coordt[t] - coordt[v]);
                    double ttol = 100.0 * DBL_EPSILON *
                                  fmax(1.0, fmax(fabs(coordt[t]), fabs(coordt[v])));
                    int ub = find_bin(bint, *nbint, lagt);
                    if (lagt > *maxtime + ttol || ub < 0) continue;

                    for (j = 0; j < ns[v]; j++) {
                        double x = data[i + NS[t]];
                        double y = data[j + NS[v]];
                        double lags;
                        double diff2;

                        if (!R_FINITE(x) || !R_FINITE(y)) continue;
                        diff2 = 0.5 * (x - y) * (x - y);

                        if (i == j && ns[t] == ns[v]) {
                            momt[ub] += diff2;
                            lbint[ub]++;
                            continue;
                        }

                        lags = dist(type[0], coordx[i + NS[t]], coordx[j + NS[v]],
                                    coordy[i + NS[t]], coordy[j + NS[v]],
                                    coordz[i + NS[t]], coordz[j + NS[v]], *REARTH);
                        if (lags <= *maxdist) {
                            int hb = find_bin(bins, *nbins, lags);
                            if (hb >= 0) {
                                int q = hb * nbint1 + ub;
                                momst[q] += diff2;
                                lbinst[q]++;
                            }
                        }
                    }
                }
            }
        }
    }
}


// binned spatial-temporal variogram:
void Binned_Variogram_st2_dyn(
    double *bins, double *bint, double *coordx, double *coordy, double *coordz, double *coordt, double *data,
    int *lbins, int *lbinst, int *lbint, double *moms, double *momst, double *momt,
    int *nbins, int *nbint, int *ns, int *NS)
{
    int i, j, t, v;
    int nbins1 = *nbins - 1;
    int nbint1 = *nbint - 1;
    double mm[2];
    double step;
    double temporal_tol;

    if (nbins1 <= 0 || nbint1 <= 0) return;

    for (int h = 0; h < nbins1; h++) { moms[h] = 0.0; lbins[h] = 0; }
    for (int u = 0; u < nbint1; u++) { momt[u] = 0.0; lbint[u] = 0; }
    for (int q = 0; q < nbins1 * nbint1; q++) { momst[q] = 0.0; lbinst[q] = 0; }

    Maxima_Minima_dist(mm, coordx, coordy, coordz, ncoord, type, REARTH);
    {
        double upper = mm[1];
        if (maxdist[0] < upper) upper = maxdist[0];
        if (!R_FINITE(mm[0]) || !R_FINITE(upper) || upper < mm[0]) {
            double hi = (R_FINITE(maxdist[0]) && maxdist[0] > 0.0) ? maxdist[0] : 1.0;
            step = hi / nbins1;
            bins[0] = 0.0;
            for (int h = 1; h < *nbins; h++) bins[h] = bins[h - 1] + step;
            bins[*nbins - 1] = hi;
        } else {
            if (upper <= mm[0]) {
                double delta = sqrt(DBL_EPSILON) * fmax(1.0, fabs(mm[0]));
                bins[0] = fmax(0.0, mm[0] - delta);
                upper = mm[0] + delta;
            } else bins[0] = mm[0];
            step = (upper - bins[0]) / nbins1;
            for (int h = 1; h < *nbins; h++) bins[h] = bins[h - 1] + step;
            bins[*nbins - 1] = upper;
        }
    }

    /* The temporal margin is gamma(0,u): only collocated points, up to
       numerical distance tolerance, are used.  It no longer depends on the
       width of the first spatial bin. */
    temporal_tol = 100.0 * DBL_EPSILON * fmax(1.0, fabs(mm[1]));

    for (t = 0; t < ntime[0]; t++) {
        for (i = 0; i < ns[t]; i++) {
            for (v = t; v < ntime[0]; v++) {
                if (t == v) {
                    for (j = i + 1; j < ns[v]; j++) {
                        double x = data[i + NS[t]];
                        double y = data[j + NS[v]];
                        double lags;
                        int hb;

                        if (!R_FINITE(x) || !R_FINITE(y)) continue;
                        lags = dist(type[0], coordx[i + NS[t]], coordx[j + NS[v]],
                                    coordy[i + NS[t]], coordy[j + NS[v]],
                                    coordz[i + NS[t]], coordz[j + NS[v]], *REARTH);
                        if (lags > *maxdist) continue;
                        hb = find_bin(bins, *nbins, lags);
                        if (hb < 0) continue;
                        {
                            double d = x - y;
                            moms[hb] += 0.5 * d * d;
                            lbins[hb]++;
                        }
                    }
                } else {
                    double lagt = fabs(coordt[t] - coordt[v]);
                    double ttol = 100.0 * DBL_EPSILON *
                                  fmax(1.0, fmax(fabs(coordt[t]), fabs(coordt[v])));
                    int ub = find_bin(bint, *nbint, lagt);
                    if (lagt > *maxtime + ttol || ub < 0) continue;

                    for (j = 0; j < ns[v]; j++) {
                        double x = data[i + NS[t]];
                        double y = data[j + NS[v]];
                        double lags;
                        double diff2;

                        if (!R_FINITE(x) || !R_FINITE(y)) continue;
                        lags = dist(type[0], coordx[i + NS[t]], coordx[j + NS[v]],
                                    coordy[i + NS[t]], coordy[j + NS[v]],
                                    coordz[i + NS[t]], coordz[j + NS[v]], *REARTH);
                        diff2 = 0.5 * (x - y) * (x - y);

                        if (R_FINITE(lags) && lags <= temporal_tol) {
                            momt[ub] += diff2;
                            lbint[ub]++;
                            continue;
                        }

                        if (R_FINITE(lags) && lags <= *maxdist) {
                            int hb = find_bin(bins, *nbins, lags);
                            if (hb >= 0) {
                                int q = hb * nbint1 + ub;
                                momst[q] += diff2;
                                lbinst[q]++;
                            }
                        }
                    }
                }
            }
        }
    }
}





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
    double mm[2], md, step;

    Maxima_Minima_dist(mm, coordx, coordy, coordz, ncoord, type, REARTH);
    md = fmax(dista[0][1], fmax(dista[1][1], dista[0][0]));
    if (md < mm[1]) mm[1] = md;
    if (*nbins < 2) return;

    step = mm[1] / (*nbins - 1);
    bins[0] = 0.0;
    for (h = 1; h < *nbins; h++) bins[h] = bins[h - 1] + step;
    bins[*nbins - 1] = mm[1];

    for (t = 0; t < ntime[0]; t++) {
        for (i = 0; i < ns[t]; i++) {
            int idx_i_t = i + NS[t];
            for (v = t; v < ntime[0]; v++) {
                if (t == v) {
                    for (j = i + 1; j < ns[t]; j++) {
                        int idx_j_t = j + NS[t];
                        double lags = dist(type[0],
                                           coordx[idx_i_t], coordx[idx_j_t],
                                           coordy[idx_i_t], coordy[idx_j_t],
                                           coordz[idx_i_t], coordz[idx_j_t], *REARTH);
                        if (lags > dista[t][v]) continue;
                        h = find_bin(bins, *nbins, lags);
                        if (h >= 0) {
                            double x = data[idx_i_t];
                            double y = data[idx_j_t];
                            if (R_FINITE(x) && R_FINITE(y)) {
                                double d = x - y;
                                marg_moms[h + t * (*nbins - 1)] += 0.5 * d * d;
                                marg_lbins[h + t * (*nbins - 1)] += 1;
                            }
                        }
                    }
                } else {
                    /* Classical cross-semivariogram on aligned supports:
                       use each unordered positive-lag spatial pair once.
                       Collocated i==j pairs contribute the trivial zero at h=0
                       and are deliberately kept out of the first positive bin. */
                    int ncommon = ns[t] < ns[v] ? ns[t] : ns[v];
                    for (j = i + 1; j < ncommon; j++) {
                        int idx_j_v = j + NS[v];
                        double lags = dist(type[0],
                                           coordx[idx_i_t], coordx[idx_j_v],
                                           coordy[idx_i_t], coordy[idx_j_v],
                                           coordz[idx_i_t], coordz[idx_j_v], *REARTH);
                        if (!R_FINITE(lags) || lags <= 0.0 || lags > dista[t][v]) continue;
                        h = find_bin(bins, *nbins, lags);
                        if (h >= 0) {
                            double x = data[idx_i_t];
                            double y = data[j + NS[t]];
                            double a = data[i + NS[v]];
                            double b = data[idx_j_v];
                            if (R_FINITE(x) && R_FINITE(y) && R_FINITE(a) && R_FINITE(b)) {
                                cross_moms[h + (v - t - 1) * (*nbins - 1)] += 0.5 * (x - y) * (a - b);
                                cross_lbins[h + (v - t - 1) * (*nbins - 1)] += 1;
                            }
                        }
                    }
                }
            }
        }
    }
}


/***********************************************************************************************************************************/
// variogram cloud:
void Cloud_Variogram2(double *bins, double *coordx, double *coordy, double *coordz, double *coordt,
                      double *data, int *lbins, double *moms, int *nbins)
{
    int h = 0;
    int npoints = ncoord[0];

    for (int i = 0; i < npoints - 1; i++) {
        double x = data[i];
        if (!R_FINITE(x)) continue;
        for (int j = i + 1; j < npoints; j++) {
            double y = data[j];
            double lags;
            double d;

            if (!R_FINITE(y)) continue;
            lags = dist(type[0], coordx[i], coordx[j],
                        coordy[i], coordy[j], coordz[i], coordz[j], *REARTH);
            if (lags > *maxdist) continue;
            if (h >= *nbins) return;

            d = x - y;
            bins[h] = lags;
            moms[h] = 0.5 * d * d;
            lbins[h] = 1;
            h++;
        }
    }
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
      if (!(lbins[i] > 0.0) || !R_FINITE(moms[i])) { i++; continue; }
      vario=moms[i]/lbins[i];// Computes the empirical variogram
      if (!R_FINITE(vario)) { i++; continue; }
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
      if (!(lbins[i] > 0.0) || !R_FINITE(moms[i])) { i++; continue; }
      vario=moms[i]/lbins[i];// Computes the empirical variogram
      if (!R_FINITE(vario)) { i++; continue; }
      varhat=Variogram(cormod,0.5*(bins[h]+bins[h+1]),bint[u],nuis[1],nuis[2],par);
      if(vario) *res=*res-pow(varhat-vario,2)*(lbins[i]/pow(vario,2));// Computes the weighted least squares
      i++;}
  return;
}
