#include "header.h"

/* Robust finite-difference helpers.
   Keep hard domain failures as LOW, but avoid artificial objective jumps
   when a single pair has an underflowed or non-positive probability. */
#ifndef LOG_FLOOR
#define LOG_FLOOR -700.0
#endif
#ifndef PROB_FLOOR
#define PROB_FLOOR 1e-300
#endif

static R_INLINE double Geo_safe_log_prob(double p)
{
    return (!R_FINITE(p) || p <= 0.0) ? LOG_FLOOR : log(p);
}

/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
/******************************************************************************************/

// Composite marginal (pariwise) log-likelihood for the spatial Gaussian model:
void Comp_Pair_Gauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                        double *par, int *weigthed, double *res, double *mean1, double *mean2,
                        double *nuis,  int *type_cop, int *cond)
{
    // Controllo rapido dei parametri
    const double sill = nuis[1];
    const double nugget = nuis[0];
    if(sill < 0 || nugget < 0 || nugget > 1) {*res = LOW; return;}

    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double nn = 1.0 - nugget;
    
    double total = 0.0;

    /* Opt-in OpenMP path for Gaussian spatial Matérn pairwise likelihood. */
    if (GeoModels_gauss_pair_omp(0, *cormod, data1, data2, par, weighted,
                                 mean1, mean2, n_pairs, max_dist,
                                 sill, nugget, &total)) {
        *res = R_FINITE(total) ? total : LOW;
        return;
    }

    // Loop principale separato per i due casi (con e senza pesi)
    if(weighted) {
        // Caso con pesi
        for(int i = 0; i < n_pairs; i++) {
            const double d1 = data1[i];
            const double d2 = data2[i];
            
            if(!ISNAN(d1) && !ISNAN(d2)) {
                const double lag = lags[i];
                const double corr = CorFct(cormod, lag, 0, par, 0, 0);
                const double weights = GeoCachedPairWeight(i, 1, 0);
                
                // Calcolo della verosimiglianza
                total += log_biv_Norm(nn * clamp_corr(corr), d1, d2, mean1[i], mean2[i], sill, 0) * weights;
            }
        }
    } else {
        // Caso senza pesi
        for(int i = 0; i < n_pairs; i++) {
            const double d1 = data1[i];
            const double d2 = data2[i];
            
            if(!ISNAN(d1) && !ISNAN(d2)) {
                const double lag = lags[i];
                const double corr = CorFct(cormod, lag, 0, par, 0, 0);
                
                // Calcolo della verosimiglianza
                total += log_biv_Norm(nn * clamp_corr(corr), d1, d2, mean1[i], mean2[i], sill, 0);
            }
        }
    }

    // Assegnazione finale con controllo
    *res = R_FINITE(total) ? total : LOW;
}
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
void Comp_Pair_WrapGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                           double *par, int *weigthed, double *res, double *mean1, double *mean2,
                           double *nuis,  int *type_cop, int *cond)
{
    // Controllo rapido dei parametri
    const double sill = nuis[1];
    const double nugget = nuis[0];
       if(sill <= 0 || nugget < 0 || nugget >= 1) {*res = LOW; return;}
    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];

    const double nn = 1.0 - nugget;  // Rinomina da "scale" a "nn"
    const double alfa = 3.0;
    double total = 0.0;

    // Loop principale separato per i due casi (con e senza pesi)
    if(weighted) {
        // Caso con pesi
        for(int i = 0; i < n_pairs; i++) {
            const double d1 = data1[i];
            const double d2 = data2[i];
            
            if(!ISNAN(d1) && !ISNAN(d2)) {
                const double lag = lags[i];
                const double corr = CorFct(cormod, lag, 0, par, 0, 0);
                const double wrap_gauss = biv_wrapped_log(alfa, d1, d2, mean1[i], mean2[i], 
                                                   nugget, sill, nn * clamp_corr(corr));  // Utilizzo di nn
                const double weights = GeoCachedPairWeight(i, 1, 0);
                total += wrap_gauss * weights;
            }
        }
    } else {
        // Caso senza pesi
        for(int i = 0; i < n_pairs; i++) {
            const double d1 = data1[i];
            const double d2 = data2[i];
            
            if(!ISNAN(d1) && !ISNAN(d2)) {
                const double lag = lags[i];
                const double corr = CorFct(cormod, lag, 0, par, 0, 0);
                const double wrap_gauss = biv_wrapped_log(alfa, d1, d2, mean1[i], mean2[i], 
                                                   nugget, sill, nn * clamp_corr(corr));  // Utilizzo di nn
                total += wrap_gauss;
            }
        }
    }
    // Assegnazione finale con controllo
    *res = R_FINITE(total) ? total : LOW;
}
/******************************************************************************************/
void Comp_Pair_SinhGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                              double *par, int *weigthed, double *res, double *mean1, double *mean2,
                              double *nuis,  int *type_cop, int *cond)
{
    const double nuis0 = nuis[0];
    const double nuis1 = nuis[1];
    const double nuis2 = nuis[2];
    const double nuis3 = nuis[3];

    if (nuis3 < 0 || nuis1 < 0 || nuis0 < 0 || nuis0 >= 1) {*res = LOW;return;}

    const int weighted = *weigthed;
    const int n_pairs = npairs[0];

    const double scale = 1.0 - nuis0;
    double total = 0.0;

    if (weighted) {
        for (int i = 0; i < n_pairs; i++) {
            const double d1 = data1[i];
            const double d2 = data2[i];

            if (!ISNAN(d1) && !ISNAN(d2)) {
                const double lag = lags[i];
                const double corr = CorFct(cormod, lag, 0, par, 0, 0);
                const double sinh_val = biv_sinh(scale * clamp_corr(corr), d1, d2, mean1[i], mean2[i], nuis2, nuis3, nuis1);
                const double weights = GeoCachedPairWeight(i, 1, 0);
                total += weights * Geo_safe_log_prob(sinh_val);
            }
        }
    } else {
        for (int i = 0; i < n_pairs; i++) {
            const double d1 = data1[i];
            const double d2 = data2[i];

            if (!ISNAN(d1) && !ISNAN(d2)) {
                const double lag = lags[i];
                const double corr = CorFct(cormod, lag, 0, par, 0, 0);
                const double sinh_val = biv_sinh(scale * clamp_corr(corr), d1, d2, mean1[i], mean2[i], nuis2, nuis3, nuis1);
                total += Geo_safe_log_prob(sinh_val);
            }
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}


void Comp_Pair_SkewLaplace2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                              double *par, int *weigthed, double *res, double *mean1, double *mean2,
                              double *nuis,  int *type_cop, int *cond)
{
    // Controllo precoce dei parametri
    const double nugget = nuis[0];
      const double skew=nuis[2];
      const double sill=nuis[1];

    
    if(nugget < 0 || nugget >= 1 || skew < 0 || skew >= 1) {
        *res = LOW;
        return;
    }

    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];

    const double scale = 1.0 - nugget;
  
    double total = 0.0; double val =0.0;
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
       // Rprintf("%f %f -%f %f %f \n",mean1[1],mean1[2],par[1],sill,skew);
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Calcolo correlazione
            const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
            
            const double weights = weighted ? GeoCachedPairWeight(i, 1, 0) : 1.0;

         val=log_biv_skewlaplace(scale*clamp_corr(corr),d1,d2,mean1[i],mean2[i],sill, skew);
   
            total +=  val*weights;
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}

/******************************************************************************************/
void Comp_Pair_SkewGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                               double *par, int *weigthed, double *res, double *mean1, double *mean2,
                               double *nuis,  int *type_cop, int *cond)
{
    const double sill = nuis[1];
    const double nugget = nuis[0];
      if(sill < 0 || nugget < 0 || nugget > 1) {*res = LOW; return;}
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];

    const double skew_param = nuis[2];
    double total = 0.0;

    if (weighted) {
        for (int i = 0; i < n_pairs; i++) {
            const double d1 = data1[i];
            const double d2 = data2[i];

            if (!ISNAN(d1) && !ISNAN(d2)) {
                const double lag = lags[i];
                const double corr = CorFct(cormod, lag, 0, par, 0, 0);

                const double skew_val = biv_skew(clamp_corr(corr), d1, d2, mean1[i], mean2[i], sill, skew_param, nugget);
                const double weights = GeoCachedPairWeight(i, 1, 0);
                total += weights * Geo_safe_log_prob(skew_val);
            }
        }
    } else {
        for (int i = 0; i < n_pairs; i++) {
            const double d1 = data1[i];
            const double d2 = data2[i];

            if (!ISNAN(d1) && !ISNAN(d2)) {
                const double lag = lags[i];
                const double corr = CorFct(cormod, lag, 0, par, 0, 0);

                const double skew_val = biv_skew(clamp_corr(corr), d1, d2, mean1[i], mean2[i], sill, skew_param, nugget);
                total += Geo_safe_log_prob(skew_val);
            }
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}

/*********************************************************/
void Comp_Pair_Gamma2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                         double *par, int *weigthed, double *res, double *mean1, double *mean2,
                         double *nuis,  int *type_cop, int *cond)
{
    const double nugget      = nuis[0];
    const double shape_param = nuis[2];

    // guardie più corrette
    if (!R_FINITE(nugget) || nugget < 0.0 || nugget >= 1.0 ||
        !R_FINITE(shape_param) || shape_param <= 0.0) {
        *res = LOW;
        return;
    }

    const int weighted   = *weigthed;
    const int n_pairs    = npairs[0];

    const double scale   = 1.0 - nugget;

    double total = 0.0;

    for (int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];

        if (ISNAN(d1) || ISNAN(d2)) continue;

        const double lag  = lags[i];
        double corr0 = CorFct(cormod, lag, 0, par, 0, 0);
        if (!R_FINITE(corr0)) continue;
        // corr effettiva con nugget e clamp numerico
        const double corr = clamp_corr(scale * corr0);
        // peso
        const double w = weighted ? GeoCachedPairWeight(i, 1, 0) : 1.0;

        // LOG-densità (no log(biv_gamma()))
        double logg = biv_gamma(corr, d1, d2, mean1[i], mean2[i], shape_param);

        // se ancora non finita, metti un floor "morbido" invece di far collassare tutto
        if (!R_FINITE(logg)) logg = LOG_FLOOR;

        total += w * logg;
      
      
    }

    // normalizzazione (consigliata): stessa massimizzazione, migliore numerica
    *res = total; /// sumw;
}
/*********************************************************/
void Comp_Pair_Weibull2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                            double *par, int *weigthed, double *res, double *mean1, double *mean2,
                            double *nuis,  int *type_cop, int *cond)
{
    // Controllo precoce dei parametri
    const double nugget = nuis[0];
    const double shape_param = nuis[2];
    
    if (nugget < 0 || nugget >= 1 || shape_param < 0) {*res = LOW;return; }

    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];

    const double scale = 1.0 - nugget;
    double total = 0.0;
    double weights = 1.0;  // Valore di default

    if (weighted) {
        for (int i = 0; i < n_pairs; i++) {
            const double d1 = data1[i];
            const double d2 = data2[i];
            
            if (!ISNAN(d1) && !ISNAN(d2)) {
                const double lag = lags[i];
                const double corr = CorFct(cormod, lag, 0, par, 0, 0);

                weights = GeoCachedPairWeight(i, 1, 0); // Calcolo peso

                const double wbl = biv_Weibull(scale * corr, d1, d2, mean1[i], mean2[i], shape_param);
                total += weights * wbl;
            }
        }
    } else {
        for (int i = 0; i < n_pairs; i++) {
            const double d1 = data1[i];
            const double d2 = data2[i];
            
            if (!ISNAN(d1) && !ISNAN(d2)) {
                const double lag = lags[i];
                const double corr = CorFct(cormod, lag, 0, par, 0, 0);

                const double wbl = biv_Weibull(scale * corr, d1, d2, mean1[i], mean2[i], shape_param);
                total += wbl;
            }
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}

/*********************************************************/
void Comp_Pair_Kumaraswamy2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i;double corr,zi,zj,weights=1.0,bl;
    double nugget=nuis[0];
         double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                     if(*weigthed) weights=GeoCachedPairWeight(i, 1, 0);
                  bl=biv_Kumara((1-nugget)*clamp_corr(corr),zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);

        *res += weights * Geo_safe_log_prob(bl);
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Kumaraswamy22mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i;double corr,zi,zj,weights=1.0,bl;
    double nugget=nuis[0];
         double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                     if(*weigthed) weights=GeoCachedPairWeight(i, 1, 0);
                  bl=biv_Kumara2((1-nugget)*clamp_corr(corr),zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);

        *res += weights * Geo_safe_log_prob(bl);
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Beta2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i;double corr,zi,zj,weights=1.0,bl;
    double nugget=nuis[0];
    double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max)  {*res=LOW;  return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                     if(*weigthed) weights=GeoCachedPairWeight(i, 1, 0);
                  bl=biv_beta((1-nugget)*clamp_corr(corr),zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);
        *res += weights * Geo_safe_log_prob(bl);
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_LogGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                             double *par, int *weigthed, double *res, double *mean1, double *mean2,
                             double *nuis,  int *type_cop, int *cond)
{
    // Controllo precoce dei parametri
    const double sill = nuis[1];
    const double nugget = nuis[0];
    if (sill < 0 || nugget < 0 || nugget > 1) {*res = LOW;return;}

    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];

    const double scale = 1.0 - nugget;
    double total = 0.0;
    double weights = 1.0; // Valore di default per i pesi

    if (weighted) {
        for (int i = 0; i < n_pairs; i++) {
            const double d1 = data1[i];
            const double d2 = data2[i];
            
            if (!ISNAN(d1) && !ISNAN(d2)) {
                const double lag = lags[i];
                const double corr = CorFct(cormod, lag, 0, par, 0, 0);
                
                weights = GeoCachedPairWeight(i, 1, 0); // Calcolo peso

                const double lognorm_val = d2lognorm(d1, d2, sill, nugget, mean1[i], mean2[i], scale * clamp_corr(corr));
                total += weights * Geo_safe_log_prob(lognorm_val);
            }
        }
    } else {
        for (int i = 0; i < n_pairs; i++) {
            const double d1 = data1[i];
            const double d2 = data2[i];
            
            if (!ISNAN(d1) && !ISNAN(d2)) {
                const double lag = lags[i];
                const double corr = CorFct(cormod, lag, 0, par, 0, 0);

                const double lognorm_val = d2lognorm(d1, d2, sill, nugget, mean1[i], mean2[i], scale * clamp_corr(corr));
                total += Geo_safe_log_prob(lognorm_val);
            }
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}

/*********************************************************/
void Comp_Pair_PoisbinnegGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                 double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                 double *nuis,  int *type_cop, int *cond)
{
    // Controllo precoce dei parametri
    const double nugget = nuis[0];
    if(nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }

    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];

    const int N = N1[0];  // Parametro N della Poisson-Binomiale Negativa
    const double scale = 1.0 - nugget;
    
    double total = 0.0;
    double weights = 1.0;

    // Loop principale ottimizzato
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Calcolo correlazione e probabilità
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            const double ai = mean1[i];
            const double aj = mean2[i];
            // Calcolo probabilità congiunta
            const double p11 = pbnorm22(ai, aj, scale * clamp_corr(corr));
            // Calcolo probabilità marginali
            const double p1 = pnorm(ai, 0, 1, 1, 0);
            const double p2 = pnorm(aj, 0, 1, 1, 0);
            // Conversione a interi
            const int uu = (int)d1;
            const int vv = (int)d2;
            // Calcolo pesi
            weights = weighted ? GeoCachedPairWeight(i, 1, 0) : 1.0; 
            // Calcolo verosimiglianza e accumulo
            const double binneg_val = biv_poisbinneg(N, uu, vv, p1, p2, p11);
            total += weights * log(binneg_val);
        }
    }
    // Assegnazione finale con controllo
    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Pair_PoisbinGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if(nugget<0||nugget>=1){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                   ai=mean1[i];aj=mean2[i];
                   corr=CorFct(cormod,lags[i],0,par,0,0);
                    p11=pbnorm22(ai,aj,(1-nugget)*clamp_corr(corr));
                p1=pnorm(ai,0,1,1,0);p2=pnorm(aj,0,1,1,0);
                    u=data1[i];v=data2[i];
                       if(*weigthed) weights=GeoCachedPairWeight(i, 1, 0);
                          uu=(int) u; vv=(int) v;
                        bl=biv_poisbin(N1[0],uu,vv,p1,p2,p11);

                    *res += weights * Geo_safe_log_prob(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}


/*********************************************************/



void Comp_Pair_BinomnegGauss2mem(
    int *cormod,
    double *data1,
    double *data2,
    int *N1,
    int *N2,
    double *par,
    int *weigthed,
    double *res,
    double *mean1,
    double *mean2,
    double *nuis,
    int *type_cop,
    int *cond)
{
    int i;
    int used_pairs = 0;
    int failed = 0;

    int N;
    int cache_size;
    double *lgamma_cache;

    double nugget;
    double one_minus_nugget;

    const double P_CLAMP = 1e-10;

    (void) type_cop;
    (void) cond;

    *res = 0.0;

    nugget = nuis[0];

    if (!R_FINITE(nugget) ||
        nugget < 0.0 ||
        nugget >= 1.0) {
        *res = LOW;
        return;
    }

    one_minus_nugget = 1.0 - nugget;

    N = N1[0];

    if (N <= 0 || N == INT_MAX) {
        *res = LOW;
        return;
    }

    /*
     * This model assumes one common integer r = N at all locations.
     * Do not silently ignore site-specific values.
     */
    for (i = 0; i < npairs[0]; i++) {
        if (N1[i] != N || N2[i] != N) {
            *res = LOW;
            return;
        }
    }

    /*
     * All cached factorial arguments used by the log-pmf are between
     * 0 and N-1.  N+1 entries are sufficient and avoid scanning the data.
     */
    cache_size = N + 1;

    lgamma_cache =
        (double *) R_alloc((size_t) cache_size, sizeof(double));

    for (i = 0; i < cache_size; i++)
        lgamma_cache[i] = lgammafn(i + 1.0);

    for (i = 0; i < npairs[0]; i++) {

        double weights = 1.0;
        double u_raw;
        double v_raw;
        int u;
        int v;

        double ai;
        double aj;
        double corr;
        double p1;
        double p2;
        double p11;
        double lower_p11;
        double upper_p11;
        double log_bl;

        if (ISNAN(data1[i]) || ISNAN(data2[i]))
            continue;

        u_raw = data1[i];
        v_raw = data2[i];

        /*
         * Negative-binomial observations must be non-negative integers.
         * Avoid silently truncating non-integer values.
         */
        if (!R_FINITE(u_raw) ||
            !R_FINITE(v_raw) ||
            u_raw < 0.0 ||
            v_raw < 0.0 ||
            u_raw > (double) INT_MAX ||
            v_raw > (double) INT_MAX ||
            fabs(u_raw - nearbyint(u_raw)) > 1e-8 ||
            fabs(v_raw - nearbyint(v_raw)) > 1e-8) {
            failed = 1;
            break;
        }

        u = (int) nearbyint(u_raw);
        v = (int) nearbyint(v_raw);

        if (*weigthed) {
            weights = GeoCachedPairWeight(i, 1, 0);

            if (!R_FINITE(weights)) {
                failed = 1;
                break;
            }

            if (weights <= 0.0)
                continue;
        }

        ai = mean1[i];
        aj = mean2[i];

        if (!R_FINITE(ai) || !R_FINITE(aj)) {
            failed = 1;
            break;
        }

        corr = clamp_corr(
            CorFct(cormod, lags[i], 0, par, 0, 0)
        );

        if (!R_FINITE(corr)) {
            failed = 1;
            break;
        }

        p1 = pnorm(ai, 0.0, 1.0, 1, 0);
        p2 = pnorm(aj, 0.0, 1.0, 1, 0);

        p1 = fmax(P_CLAMP, fmin(1.0 - P_CLAMP, p1));
        p2 = fmax(P_CLAMP, fmin(1.0 - P_CLAMP, p2));

        p11 = pbnorm22(
            ai,
            aj,
            one_minus_nugget * corr
        );

        if (!R_FINITE(p11)) {
            failed = 1;
            break;
        }

        /*
         * Enforce only the exact Frechet bounds.  Do not subtract
         * P_CLAMP from min(p1,p2), because equality is admissible.
         */
        lower_p11 = fmax(0.0, p1 + p2 - 1.0);
        upper_p11 = fmin(p1, p2);

        p11 = fmin(
            upper_p11,
            fmax(lower_p11, p11)
        );

        log_bl = log_biv_binomneg(
            N,
            u,
            v,
            p1,
            p2,
            p11,
            lgamma_cache,
            cache_size
        );

        /*
         * With the log-sum-exp implementation, underflow is no longer a
         * valid reason for a zero probability.  A non-finite log-pmf now
         * indicates an invalid parameter/data combination.
         */
        if (!R_FINITE(log_bl)) {
            failed = 1;
            break;
        }

        *res += weights * log_bl;
        used_pairs++;
    }


    if (failed ||
        used_pairs == 0 ||
        !R_FINITE(*res)) {
        *res = LOW;
    }
}

/******************************************************/
void Comp_Pair_BinomnegBinary2mem(int *cormod, double *data1, double *data2,
    int *N1, int *N2, double *par, int *weigthed, double *res,
    double *mean1, double *mean2, double *nuis, int *type_cop, int *cond)
{
    *res = 0.0;
 
    double nugget = nuis[0];
    if (nugget < 0.0 || nugget >= 1.0) { *res = LOW; return; }
 
    int N = N1[0];
    if (N <= 0) { *res = LOW; return; }
 
#define P_CLAMP 1e-10
 
    int    i, uu, vv;
    double u, v, bl, weights, ai, aj, corr, p1, p2, p11;
 
    for (i = 0; i < npairs[0]; i++) {
 
        if (ISNAN(data1[i]) || ISNAN(data2[i])) continue;
 
        weights = 1.0;
        if (*weigthed) {
            weights = GeoCachedPairWeight(i, 1, 0);
            if (weights <= 0.0) continue;
        }
 
        ai   = mean1[i];
        aj   = mean2[i];
        corr = clamp_corr(CorFct(cormod, lags[i], 0, par, 0, 0));
        p11  = pbnorm22(ai, aj, (1 - nugget) * corr);
        p1   = fmax(P_CLAMP, fmin(1.0 - P_CLAMP, pnorm(ai, 0, 1, 1, 0)));
        p2   = fmax(P_CLAMP, fmin(1.0 - P_CLAMP, pnorm(aj, 0, 1, 1, 0)));
 
        /* FIX 4: evita p11 == min(p1,p2) esattamente */
        p11  = fmax(fmax(0.0, p1 + p2 - 1.0),
                    fmin(fmin(p1, p2) - P_CLAMP, p11));
 
        u  = data1[i];
        v  = data2[i];
        uu = (int)u;
        vv = (int)v;
 
        /* NOTA: quando biv_binegbinary sara' aggiornata con la cache,
         * aggiungere qui la stessa logica di malloc/free delle altre funzioni */
        bl = biv_binegbinary(N, uu, vv, p1, p2, p11);
 
        if (!R_FINITE(bl) || bl <= 0.0) {
            *res += weights * LOG_FLOOR;
            continue;
        }
    
        *res += weights * log(bl);
    }
 
#undef P_CLAMP
 
    /* FIX 2: soglia robusta per npairs piccolo */

    if (!R_FINITE(*res)) *res = LOW;
}
/*********************************************************/
/*********************************************************/
void Comp_Pair_BinomnegGaussZINB2mem(int *cormod, double *data1, double *data2,
    int *N1, int *N2, double *par, int *weigthed, double *res,
    double *mean1, double *mean2, double *nuis, int *type_cop, int *cond)
{
    *res = 0.0;
 
    const double nugget1 = nuis[0];
    const double nugget2 = nuis[1];
    const double mup     = nuis[2];
 
    if (nugget1 < 0.0 || nugget1 >= 1.0 ||
        nugget2 < 0.0 || nugget2 >= 1.0) { *res = LOW; return; }
 
    const int    n_pairs  = npairs[0];

    const int    N        = N1[0];
    if (N <= 0) { *res = LOW; return; }
 
    /* FIX 1: scansione preliminare per dimensionare la cache */
    int MAX_UV = 0;
    for (int k = 0; k < n_pairs; k++) {
        if (!ISNAN(data1[k]) && data1[k] > MAX_UV) MAX_UV = (int)data1[k];
        if (!ISNAN(data2[k]) && data2[k] > MAX_UV) MAX_UV = (int)data2[k];
    }
    /* FIX 3: dimensione sicura */
    int cache_size = N + 2 * MAX_UV + 10;
    double *lgamma_cache = (double *) R_alloc((size_t) cache_size, sizeof(double));
    for (int k = 0; k < cache_size; k++)
        lgamma_cache[k] = lgammafn(k + 1);
 
    double total       = 0.0;
;
 
    for (int i = 0; i < n_pairs; i++) {
 
        if (ISNAN(data1[i]) || ISNAN(data2[i])) continue;
 
        const double lag     = lags[i];
        const double corr    = clamp_corr(CorFct(cormod, lag, 0, par, 0, 0));
        const double weights = *weigthed ? GeoCachedPairWeight(i, 1, 0) : 1.0;
        if (weights <= 0.0) continue;
 
        const int uu = (int)data1[i];
        const int vv = (int)data2[i];
 
        /* FIX 1: passa la cache allocata una sola volta */
        const double zinb_val = biv_binomnegZINB(
            N, corr, uu, vv,
            mean1[i], mean2[i],
            nugget1, nugget2, mup,
            lgamma_cache, cache_size
        );
 
        if (!R_FINITE(zinb_val) || zinb_val <= 0.0) {
            total += weights * LOG_FLOOR;
            continue;
        }
    
        total += weights * log(zinb_val);
    }
 
 


    *res = R_FINITE(total) ? total : LOW;
}

void Comp_Pair_BinomGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                               double *par, int *weigthed, double *res, double *mean1, double *mean2,
                               double *nuis, int *type_cop, int *cond)
{
    const double nugget = nuis[0];
    if (nugget >= 1 || nugget < 0) { *res = LOW; return; }

    const int    weighted = *weigthed;
    const int    n_pairs  = npairs[0];

    const int    N        = N1[0];
    const double scale    = 1.0 - nugget;
    double       total    = 0.0;

    for (int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i], d2 = data2[i];
        if (ISNAN(d1) || ISNAN(d2)) continue;

        const double lag = lags[i];

        /* correlazione con clipping */
        const double corr        = CorFct(cormod, lag, 0, par, 0, 0);
        const double scaled_corr = scale * clamp_corr(corr);

        /* parametri lineari con clipping */
        const double ai = fmax(-5.0, fmin(5.0, mean1[i]));
        const double aj = fmax(-5.0, fmin(5.0, mean2[i]));

        /* probabilità con clipping */
        const double p1  = fmax(1e-15, fmin(1-1e-15, pnorm(ai, 0, 1, 1, 0)));
        const double p2  = fmax(1e-15, fmin(1-1e-15, pnorm(aj, 0, 1, 1, 0)));
        const double p11 = fmax(1e-15, fmin(fmin(p1, p2), pbnorm22(ai, aj, scaled_corr)));

        const double binom_val = biv_binom(N, (int)d1, (int)d2, p1, p2, p11);
        const double log_val   = (binom_val > 1e-300 && R_FINITE(binom_val))
                                     ? log(binom_val)
                                     : -700.0;

        /* unica differenza tra i due rami originali */
        const double w = weighted ? GeoCachedPairWeight(i, 1, 0) : 1.0;
        total += w * log_val;
    }

    *res = R_FINITE(total) ? total : LOW;
}

/******************************************************************************************/
void Comp_Pair_BinomLogi2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 u=data1[i];v=data2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);
                p11=pblogi22(ai,aj,(1-nugget)*clamp_corr(corr));
                 p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                        if(*weigthed) weights=GeoCachedPairWeight(i, 1, 0);
                          uu=(int) u; vv=(int) v;
                        bl=biv_binom (N1[0],uu,vv,p1,p2,p11);
                    *res += weights * Geo_safe_log_prob(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Pair_BinomNNGauss2mem(int *cormod, double *data1, double *data2,
        int *N1, int *N2, double *par, int *weighted, double *res,
        double *mean1, double *mean2, double *nuis, int *type_cop, int *cond)
{
    *res = 0.0;  /* FIX 5: reset esplicito */

    double nugget = nuis[0];
    if (nugget < 0.0 || nugget >= 1.0) { *res = LOW; return; }
    double one_minus_nugget = 1.0 - nugget;

#define P_CLAMP 1e-10

    int    i, uu, vv;
    double corr, p1, p2, p11, bl, weights, ai, aj;

    for (i = 0; i < npairs[0]; i++) {

        if (ISNAN(data1[i]) || ISNAN(data2[i])) continue;

        weights = 1.0;
        if (*weighted) {
            weights = GeoCachedPairWeight(i, 1, 0);
            if (weights <= 0.0) continue;
        }

        if (data1[i] < 0.0 || data2[i] < 0.0) continue;

        /* FIX 4: N check esplicito */
        if (N1[i] <= 0 || N2[i] <= 0) continue;

        uu = (int) data1[i];
        vv = (int) data2[i];

        ai = mean1[i];
        aj = mean2[i];

        corr = clamp_corr(CorFct(cormod, lags[i], 0, par, 0, 0));
        p1 = fmax(P_CLAMP, fmin(1.0 - P_CLAMP, pnorm(ai, 0.0, 1.0, 1, 0)));
        p2 = fmax(P_CLAMP, fmin(1.0 - P_CLAMP, pnorm(aj, 0.0, 1.0, 1, 0)));

        p11 = pbnorm22(ai, aj, one_minus_nugget * corr);
        p11 = fmax(fmax(0.0, p1 + p2 - 1.0), fmin(fmin(p1, p2), p11));

        bl = biv_binom222(N1[i], N2[i], uu, vv, p1, p2, p11);

        if (!R_FINITE(bl) || bl <= 0.0) {
            *res += weights * LOG_FLOOR;
            continue;
        }
        *res += weights * log(bl);
    }



    if (!R_FINITE(*res)) *res = LOW;
}

void Comp_Pair_BinomNNGauss_misp2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, N=2,n1,n2;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success

    double M0[2] = {0.0, 0.0};
    double M1[2] = {0.0, 0.0};
    double *M[2] = {M0, M1};
    double dat[2] = {0.0, 0.0};

    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);
                 p11=pbnorm22(ai,aj,(1-nugget)*clamp_corr(corr));
                 p1=pnorm(ai,0,1,1,0);
                 p2=pnorm(aj,0,1,1,0);
                 u=data1[i];v=data2[i];
                 n1=N1[i];n2=N2[i];
                 if(*weigthed) weights=GeoCachedPairWeight(i, 1, 0);
                 M[0][0]=n1*p1*(1-p1);   M[1][1]=n2*p2*(1-p2);  // var1 var2
                 M[0][1]= fmin_int(n1,n2)*(p11-p1*p2) ;       // covariance
                 M[1][0]= M[0][1];
                 dat[0]=u-n1*p1;dat[1]=v-n2*p2; 
                   //#####
                 bl=log_dNnorm(N,M,dat);
                 *res += weights * bl;       
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}




void Comp_Pair_BinomNNLogi2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu=0,vv=0,n1,n2;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);

                 p11=pblogi22(ai,aj,(1-nugget)*clamp_corr(corr));
                 p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                 u=data1[i];v=data2[i];
                 n1=N1[i];n2=N2[i];
                 if(*weigthed) weights=GeoCachedPairWeight(i, 1, 0);
                 uu=(int) u; vv=(int) v;
                 bl=biv_binom222(n1,n2,uu,vv,p1,p2,p11);
                 *res += weights * Geo_safe_log_prob(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_LogLogistic2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                              double *par, int *weigthed, double *res, double *mean1, double *mean2,
                              double *nuis,  int *type_cop, int *cond)
{
    // Controllo precoce dei parametri
    const double nugget = nuis[0];
    const double shape_param = nuis[2]; // Parametro shape della Log-Logistica
    
    if(nugget < 0 || nugget >= 1 || shape_param <= 2) {
        *res = LOW;
        return;
    }

    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];

    const double scale = 1.0 - nugget;
    
    double total = 0.0;  // Variabile accumulatore

    // Loop principale ottimizzato
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Calcolo correlazione
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            
            // Calcolo pesi se necessario
            const double weights = weighted ? GeoCachedPairWeight(i, 1, 0) : 1.0;
            
            // Calcolo della Log-Logistica bivariata e accumulo risultato
            const double loglogistic_val = biv_LogLogistic(
                scale * clamp_corr(corr), d1, d2, 
                mean1[i], mean2[i], 
                shape_param
            );
            total += weights * Geo_safe_log_prob(loglogistic_val);
        }
    }

    // Assegnazione finale con controllo
    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Pair_Logistic2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,nugget=0.0;
        nugget=nuis[0];
    if(nugget>=1||nugget<0.0 ) {*res=LOW;  return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=(data1[i]);zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                if(*weigthed) weights=GeoCachedPairWeight(i, 1, 0);
                    bl= biv_Logistic((1-nugget)*clamp_corr(corr),zi,zj,mean1[i],mean2[i],nuis[1]);
                    *res += weights * Geo_safe_log_prob(bl);
                }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Pois2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                       double *par, int *weigthed, double *res, double *mean1, double *mean2,
                       double *nuis,  int *type_cop, int *cond)
{
    // Early parameter validation
    const double nugget = nuis[0];
    if(nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }

    // Precompute frequently used values
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];

    const double scale = 1.0 - nugget;
    
    double total = 0.0;  // Accumulator variable

    // Optimized main loop
    if (weighted) {
        // Con pesi
        for(int i = 0; i < n_pairs; i++) {
            const double d1 = data1[i];
            const double d2 = data2[i];
            
            if(!ISNAN(d1) && !ISNAN(d2)) {
                // Compute correlation
                const double lag = lags[i];
                const double corr = CorFct(cormod, lag, 0, par, 0, 0);
                
                // Compute means (exp transformed)
                const double mui = exp(mean1[i]);
                const double muj = exp(mean2[i]);
                
                // Convert to integers
                const int uu = (int)d1;
                const int ww = (int)d2;
                
                // Compute the weight
                const double weights = GeoCachedPairWeight(i, 1, 0);
                
                // Compute bivariate Poisson and accumulate result
                const double poisson_val = biv_Poisson(scale * clamp_corr(corr), uu, ww, mui, muj);
                total += poisson_val * weights;
            }
        }
    } else {
        // Senza pesi
        for(int i = 0; i < n_pairs; i++) {
            const double d1 = data1[i];
            const double d2 = data2[i];
            
            if(!ISNAN(d1) && !ISNAN(d2)) {
                // Compute correlation
                const double lag = lags[i];
                const double corr = CorFct(cormod, lag, 0, par, 0, 0);
                
                // Compute means (exp transformed)
                const double mui = exp(mean1[i]);
                const double muj = exp(mean2[i]);
                
                // Convert to integers
                const int uu = (int)d1;
                const int ww = (int)d2;
                
                // Compute bivariate Poisson and accumulate result
                const double poisson_val = biv_Poisson(scale * clamp_corr(corr), uu, ww, mui, muj);
                total += poisson_val;  // No weight applied here
            }
        }
    }

    // Final assignment with check
    *res = R_FINITE(total) ? total : LOW;
}

/*********************************************************/

void Comp_Pair_PoisGamma2mem(int *cormod, double *data1, double *data2,
                             int *N1, int *N2, double *par, int *weigthed,
                             double *res, double *mean1, double *mean2,
                             double *nuis, int *type_cop, int *cond)
{
    (void)N1; (void)N2; (void)type_cop; (void)cond;
    const double nugget = nuis[0], shape = nuis[2];
    if (!R_FINITE(nugget) || nugget < 0.0 || nugget >= 1.0 ||
        !R_FINITE(shape) || !(shape > 0.0) || !npairs || !lags) {
        *res = LOW; return;
    }
    GeoKahanSum total; GeoKahanInit(&total);
    int used = 0;
    for (int i = 0; i < npairs[0]; ++i) {
        int r, t; double mi, mj;
        const int prep = GeoPGPreparePairAt(i, data1[i], data2[i],
                                             mean1[i], mean2[i],
                                             &r, &t, &mi, &mj);
        if (prep == 0) continue;
        if (prep < 0) { *res = LOW; return; }
        const double weight = GeoCachedPairWeight(i, weigthed[0], 0);
        if (!R_FINITE(weight)) { *res = LOW; return; }
        if (weight <= 0.0) continue;
        double base = CorFct(cormod, lags[i], 0.0, par, 0, 0);
        if (!R_FINITE(base) || fabs(base) > 1.0 + 1e-8) { *res = LOW; return; }
        double rho = lags[i] <= 0.0 ? fmax2(-1.0, fmin2(1.0, base))
                                    : (1.0 - nugget) * base;
        if (rho >= 1.0 && lags[i] > 0.0) rho = 1.0 - 1e-15;
        if (rho <= -1.0 && lags[i] > 0.0) rho = -1.0 + 1e-15;
        const double ll = GeoPGLogJoint(rho, r, t, mi, mj, shape);
        if (!R_FINITE(ll)) { *res = LOW; return; }
        GeoKahanAdd(&total, weight * ll); ++used;
    }
    *res = (used > 0 && R_FINITE(total.sum)) ? total.sum : LOW;
}



void Comp_Pair_PoisGammaZIP2mem(int *cormod, double *data1, double *data2,
                                int *N1, int *N2, double *par, int *weigthed,
                                double *res, double *mean1, double *mean2,
                                double *nuis, int *type_cop, int *cond)
{
    (void)N1; (void)N2; (void)type_cop; (void)cond;
    const double nugget1 = nuis[0], nugget2 = nuis[1];
    const double mup = nuis[2], shape = nuis[4];
    if (!R_FINITE(nugget1) || !R_FINITE(nugget2) ||
        nugget1 < 0.0 || nugget1 >= 1.0 || nugget2 < 0.0 || nugget2 >= 1.0 ||
        !R_FINITE(mup) || !R_FINITE(shape) || !(shape > 0.0) ||
        !npairs || !lags) { *res = LOW; return; }
    GeoKahanSum total; GeoKahanInit(&total); int used = 0;
    for (int i = 0; i < npairs[0]; ++i) {
        int r, t; double mi, mj;
        const int prep = GeoPGPreparePairAt(i, data1[i], data2[i],
                                             mean1[i], mean2[i],
                                             &r, &t, &mi, &mj);
        if (prep == 0) continue;
        if (prep < 0) { *res = LOW; return; }
        const double weight = GeoCachedPairWeight(i, weigthed[0], 0);
        if (!R_FINITE(weight)) { *res = LOW; return; }
        if (weight <= 0.0) continue;
        double rho = CorFct(cormod, lags[i], 0.0, par, 0, 0);
        if (!R_FINITE(rho) || fabs(rho) > 1.0 + 1e-8) { *res = LOW; return; }
        rho = fmax2(-1.0, fmin2(1.0, rho));
        const double ll = GeoPGZIPLogJoint(rho, r, t, mi, mj, mup,
                                            nugget1, nugget2, shape);
        if (!R_FINITE(ll)) { *res = LOW; return; }
        GeoKahanAdd(&total, weight * ll); ++used;
    }
    *res = (used > 0 && R_FINITE(total.sum)) ? total.sum : LOW;
}



void Comp_Pair_Gauss_misp_PoisGamma2mem(int *cormod, double *data1, double *data2,
                                        int *N1, int *N2, double *par,
                                        int *weigthed, double *res,
                                        double *mean1, double *mean2,
                                        double *nuis, int *type_cop, int *cond)
{
    (void)N1; (void)N2; (void)type_cop; (void)cond;
    const double nugget = nuis[0], shape = nuis[2];
    if (!R_FINITE(nugget) || nugget < 0.0 || nugget >= 1.0 ||
        !R_FINITE(shape) || !(shape > 0.0) || !npairs || !lags) {
        *res = LOW; return;
    }
    GeoKahanSum total; GeoKahanInit(&total); int used = 0;
    for (int i = 0; i < npairs[0]; ++i) {
        if (ISNAN(data1[i]) || ISNAN(data2[i])) continue;
        if (!R_FINITE(data1[i]) || !R_FINITE(data2[i]) ||
            !R_FINITE(mean1[i]) || !R_FINITE(mean2[i]) ||
            mean1[i] > 700.0 || mean2[i] > 700.0) { *res = LOW; return; }
        const double mi = exp(mean1[i]), mj = exp(mean2[i]);
        const double weight = GeoCachedPairWeight(i, weigthed[0], 0);
        if (!R_FINITE(weight)) { *res = LOW; return; }
        if (weight <= 0.0) continue;
        double base = CorFct(cormod, lags[i], 0.0, par, 0, 0);
        if (!R_FINITE(base) || fabs(base) > 1.0 + 1e-8) { *res = LOW; return; }
        double rho = lags[i] <= 0.0 ? fmax2(-1.0, fmin2(1.0, base))
                                    : (1.0 - nugget) * base;
        if (rho >= 1.0) rho = 1.0 - 1e-15;
        if (rho <= -1.0) rho = -1.0 + 1e-15;
        const double ll = GeoPGGaussLogJoint(rho, data1[i], data2[i],
                                              mi, mj, shape);
        if (!R_FINITE(ll)) { *res = LOW; return; }
        GeoKahanAdd(&total, weight * ll); ++used;
    }
    *res = (used > 0 && R_FINITE(total.sum)) ? total.sum : LOW;
}



/*********************************************************/
void Comp_Pair_PoisZIP2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu,ww;
    double weights=1.0,corr,mui,muj,bl;
   double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];


      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}

      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],0,par,0,0);
                    // if(fabs(corr)>1|| !R_FINITE(corr)) {*res=LOW; return;}
                        if(*weigthed) weights=GeoCachedPairWeight(i, 1, 0);
                      uu=(int) data1[i];  ww=(int) data2[i];
                      bl=biv_PoissonZIP(clamp_corr(corr),uu,ww,mui, muj,mup,nugget1,nugget2);
                      *res += weights * Geo_safe_log_prob(bl);
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}



/*********************************************************/
void Comp_Pair_Gauss_misp_Pois2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                 double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                 double *nuis,  int *type_cop, int *cond)
{
    // Early parameter validation
    const double nugget = nuis[0];
    if(nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }

    // Constants and precomputations
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];

    const int N = 2; // Fixed matrix dimension
    const double scale = 1.0 - nugget;
    
    double total = 0.0;  // Accumulator variable

    // Fixed-size covariance workspace (no heap allocation).
    double M0[2], M1[2], dat[2];
    double *M[2] = {M0, M1};
// Main processing loop
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Compute means
            const double mui = exp(mean1[i]);
            const double muj = exp(mean2[i]);
            
            // Compute correlations
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0) * scale;
            const double corr1 = corr_pois(clamp_corr(corr), mui, muj);
            
            // Compute weights if needed
            const double weights = weighted ? GeoCachedPairWeight(i, 1, 0) : 1.0;
            
            // Build covariance matrix efficiently
            const double sqrt_mui_muj = sqrt(mui * muj);
            M[0][0] = mui;
            M[1][1] = muj;
            M[0][1] = M[1][0] = sqrt_mui_muj * corr1;
            
            // Center data
            dat[0] = d1 - mui;
            dat[1] = d2 - muj;
            
            // Compute density and accumulate result
            total += log_dNnorm(N, M, dat) * weights;
        }
    }

    // No workspace cleanup is required.


    // Final assignment with check
    *res = R_FINITE(total) ? total : LOW;
}

void Comp_Pair_Gauss_misp_PoisZIP2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;
    double weights=1.0,corr,mui,muj,bl;
    double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];

      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],0,par,0,0);
                      if(*weigthed) weights=GeoCachedPairWeight(i, 1, 0);
                      bl=biv_Mis_PoissonZIP(clamp_corr(corr),data1[i],data2[i],mui, muj,mup,nugget1,nugget2);
                      *res += weights * Geo_safe_log_prob(bl);
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_Gauss_misp_SkewT2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;
    double weights=1.0,sill,nugget,skew,corr,corr2,df,bl;


    df=1/nuis[0];
    nugget=nuis[1];
    sill=nuis[2];
    skew=nuis[3];

    if(ISNAN(df)||df<=2||!R_FINITE(skew)||fabs(skew)>1||
       !R_FINITE(sill)||sill<=0||!R_FINITE(nugget)||nugget<0||nugget>=1) {
        *res=LOW; return;
    }
    const double log_factor=GeoStudentCorrLogFactor(df);
    const double factor=exp(log_factor);
    const double student_variance=1.0/(1.0-2.0/df);
    const double MM=skew*sqrt(student_variance*(2.0/M_PI)*factor);
    const double FF=student_variance*(1.0-skew*skew*(2.0/M_PI)*factor);

     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                     corr=clamp_corr(CorFct(cormod,lags[i],0,par,0,0));
                     corr2=GeoSkewStudentCorrFromLogFactor(corr,(1-nugget)*corr,df,skew,log_factor);
                     if(*weigthed) weights=GeoCachedPairWeight(i, 1, 0);
                          bl=log_biv_Norm(corr2,data1[i],data2[i],mean1[i]+sqrt(sill)*MM,
                                                                 mean2[i]+sqrt(sill)*MM,sill*FF,0);
                        *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Gauss_misp_T2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                               double *par, int *weigthed, double *res, double *mean1, double *mean2,
                               double *nuis,  int *type_cop, int *cond)
{
    const double df_param = nuis[0], nugget = nuis[1], sill = nuis[2];
    if (!R_FINITE(df_param) || !R_FINITE(nugget) || !R_FINITE(sill) ||
        df_param < 0.0 || df_param >= 0.5 || nugget < 0.0 ||
        nugget >= 1.0 || sill <= 0.0) {
        *res = LOW;
        return;
    }
    const double nu = df_param == 0.0 ? R_PosInf : 1.0 / df_param;
    const double variance = sill / (1.0 - 2.0 * df_param);
    const double log_factor = GeoStudentCorrLogFactor(nu);
    const double scale = 1.0 - nugget;
    double total = 0.0;
    for (int i = 0; i < npairs[0]; ++i) {
        if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
            const double lag = lags[i];
            const double weights = *weigthed ? GeoCachedPairWeight(i, 1, 0) : 1.0;
            if (weights == 0.0) continue;
            const double rho = CorFct(cormod, lag, 0.0, par, 0, 0);
            const double corr = GeoStudentCorrFromLogFactor(
                rho, scale * rho, nu, log_factor);
            double log_pair = log_biv_Norm(clamp_corr(corr), data1[i], data2[i],
                                           mean1[i], mean2[i], variance, 0);
            total += weights * log_pair;
        }
    }
    *res = R_FINITE(total) ? total : LOW;
}


/*********************************************************/
 void Comp_Pair_Tukeyhh2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                          double *par, int *weigthed, double *res, double *mean1, double *mean2,
                          double *nuis,  int *type_cop, int *cond)
{
    // Early parameter validation
    const double sill = nuis[1];
    const double nugget = nuis[0];
    const double h1 = nuis[3];
    const double h2 = nuis[2];
    
    if(sill <= 0 || h1 < 0 || h1 >= 0.5 || h2 < 0 || h2 >= 0.5 || nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }

    // Precompute frequently used values
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];

    const double scale = 1.0 - nugget;
    
    double total = 0.0;  // Accumulator variable

    // Main processing loop
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Compute correlation
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            
            // Compute weights if needed
            const double weights = weighted ? GeoCachedPairWeight(i, 1, 0) : 1.0;
            
            // Compute bivariate Tukey h-h likelihood
            const double tukey_val = biv_tukey_hh(
                scale * clamp_corr(corr), d1, d2, 
                mean1[i], mean2[i], 
                sill, h1, h2
            );
            
            // Accumulate result
            total += weights * Geo_safe_log_prob(tukey_val);
        }
    }

    // Final assignment with check
    *res = R_FINITE(total) ? total : LOW;
}




/*********************************************************/
void Comp_Pair_Tukeyh2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                         double *par, int *weigthed, double *res, double *mean1, double *mean2,
                         double *nuis,  int *type_cop, int *cond)
{
    // Early parameter validation
    const double sill = nuis[1];
    const double nugget = nuis[0];
    const double tail = nuis[2];
    
    if(sill < 0 || tail < 0 || tail > 0.5 || nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }

    // Precompute frequently used values
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];

    const double scale = 1.0 - nugget;
    
    double total = 0.0;  // Accumulator variable

    // Main processing loop
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Compute correlation
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0) * scale;
            
            // Compute weights if needed
            const double weights = weighted ? GeoCachedPairWeight(i, 1, 0) : 1.0;
            
            // Compute bivariate Tukey h likelihood and accumulate
            const double tukey_val = biv_tukey_h(
                clamp_corr(corr), d1, d2, 
                mean1[i], mean2[i], 
                tail, sill
            );
            total += weights * Geo_safe_log_prob(tukey_val);
        }
    }

    // Final assignment with check
    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Pair_Gauss_misp_Tukeygh2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i;double bl,corr,corr2,zi,zj,weights=1.0,eta,tail,sill,nugget,u,eta2,mu,vv;
    eta  = nuis[2];  //skewness parameter
    tail = nuis[3];  //tail parameter
    sill =nuis[1];
    nugget=nuis[0];

    eta2=eta*eta;
    u=1-tail;
    mu=(exp(eta2/(2*u))-1)/(eta*sqrt(u));
    vv=((exp(2*eta2/(1-2*tail))-2*exp(eta2/(2*(1-2*tail)))+1)/(eta2*
                           sqrt(1-2*tail))-mu*mu);
        if(fabs(eta)<1e-5)
           {
           mu=0.0;
           vv=R_pow(1-2*tail,-3/2);
           }
         if(sill<0||nugget<0||nugget>=1||tail<0||tail>0.5) {*res=LOW;  return;}
   for(i=0;i<npairs[0];i++){
          zi=data1[i];zj=data2[i];
if(!ISNAN(zi)&&!ISNAN(zj) ){

                    corr=(1-nugget)*CorFct(cormod,lags[i],0,par,0,0);
                    corr2=corr_tukeygh(clamp_corr(corr),eta,tail);
                    if(*weigthed) weights=GeoCachedPairWeight(i, 1, 0);
                bl=log_biv_Norm(corr2,zi,zj,mean1[i]+sqrt(sill)*mu,
                                            mean2[i]+sqrt(sill)*mu, sill*vv,0);
                    *res+= weights*bl;
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_T2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                    double *par, int *weigthed, double *res, double *mean1, double *mean2,
                    double *nuis,  int *type_cop, int *cond) 
{
    // Controllo precoce dei parametri
    const double sill = nuis[2];
    const double nugget = nuis[1];
    const double df = nuis[0];
    if(sill <= 0 || nugget < 0 || nugget >= 1 || df <= 0 || df > 0.5) {
        *res = LOW; 
        return;
    }
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];

    const double inv_sqrt_sill = 1.0 / sqrt(sill);
    const double inv_sill = 1.0 / sill;
    double total = 0.0;  // Variabile accumulatore
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Standardizzazione dei dati
            const double std_d1 = (d1 - mean1[i]) * inv_sqrt_sill;
            const double std_d2 = (d2 - mean2[i]) * inv_sqrt_sill;
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            const double weights = weighted ? GeoCachedPairWeight(i, 1, 0) : 1.0;
            const double t_val = biv_T(clamp_corr(corr), std_d1, std_d2, df, nugget) * inv_sill;
            total += weights * Geo_safe_log_prob(t_val);
        }
    }

    // Assegnazione finale con controllo
    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Pair_TWOPIECETukeyh2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,eta,tail,qq,sill,nugget;
    eta  = nuis[2];  //skewness parameter
    tail = nuis[3];  //tail parameter
    sill =nuis[1];
    nugget=nuis[0];

     if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1||tail<0||tail>0.5) {*res=LOW;  return;}

       qq=qnorm((1-eta)/2,0,1,1,0);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);
                      p11=pbnorm22(qq,qq,clamp_corr(corr));
                    if(*weigthed) weights=GeoCachedPairWeight(i, 1, 0);
                    bl=biv_two_pieceTukeyh((1-nugget)*clamp_corr(corr),zi,zj,sill,eta,tail,p11,mean1[i],mean2[i]);
                    *res += weights * Geo_safe_log_prob(bl);
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_TWOPIECET2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,qq;

    double eta=nuis[3];  //skewness parameter
    double df=nuis[0];
    double sill=nuis[2];
    double nugget=nuis[1];
    if(sill<0||nugget<0||nugget>=1 ||fabs(eta)>1|| df >0.5||df<0) {*res=LOW;  return;}

       qq=qnorm((1-eta)/2,0,1,1,0);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                      p11=pbnorm22(qq,qq,clamp_corr(corr));
                    if(*weigthed) weights=GeoCachedPairWeight(i, 1, 0);
                    /********************************************************/
    bl=biv_two_pieceT(clamp_corr(corr),zi,zj,sill,df,eta,p11,mean1[i],mean2[i],nugget);
                    /********************************************************/
                           *res += weights * Geo_safe_log_prob(bl);
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_TWOPIECEGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                double *nuis,  int *type_cop, int *cond)
{
    // Early parameter validation
    const double eta = nuis[2];  // skewness parameter
    const double sill = nuis[1];
    const double nugget = nuis[0];
    
    const double qq = qnorm((1 - eta) / 2, 0, 1, 1, 0);
    if(fabs(eta) >= 1 || sill <= 0 || nugget < 0 || nugget >= 1) {  // Changed > to >= for eta
        *res = LOW;
        return;
    }
    // Precompute frequently used values
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];

    const double scale = 1.0 - nugget;
    double total = 0.0;  // Accumulator variable
    // Main processing loop
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Compute correlation
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            const double p11 = pbnorm22(qq, qq, clamp_corr(corr));
            const double weights = weighted ? GeoCachedPairWeight(i, 1, 0) : 1.0;
            const double density = biv_two_pieceGaussian(
                scale * clamp_corr(corr), d1, d2, 
                sill, eta, p11, 
                mean1[i], mean2[i]
            );
            total += weights * log(density);
        }
    }

    // Final assignment with check
    *res = R_FINITE(total) ? total : LOW;
}

void Comp_Pair_TWOPIECEBIMODAL2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,eta,qq,sill,df,nugget,delta;

    eta=nuis[4];  //skewness parameter
    delta=nuis[3];
    sill=nuis[2];
    nugget=nuis[1];
    df=nuis[0];
 if( fabs(eta)>1||df<0||nugget<0||nugget>=1||delta<0||sill<0) {*res=LOW;  return;}

    qq=qnorm((1-eta)/2,0,1,1,0);

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                        p11=pbnorm22(qq,qq,clamp_corr(corr));
                    if(*weigthed) weights=GeoCachedPairWeight(i, 1, 0);
                    /********************************************************/
                   bl=biv_two_piece_bimodal((1-nugget)*clamp_corr(corr),zi,zj,sill,df,delta,eta,p11,mean1[i],mean2[i]);
                    /********************************************************/
                           *res += weights * Geo_safe_log_prob(bl);
                }}}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}

/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/********************* SPATI0-TEMPORAL CASE ***********************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/

void Comp_Pair_Gauss_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                           double *par, int *weigthed, double *res, double *mean1, double *mean2,
                           double *nuis,  int *type_cop, int *cond)
{
    // Estrai valori una sola volta
    const int n_pairs = npairs[0];

    const double sill = nuis[1];
    const double nugget = nuis[0];
    
    // Controllo preliminare dei parametri
    if(sill < 0 || nugget < 0 || nugget > 1) {
        *res = LOW; 
        return;
    }

    double total = 0.0;
    const int weighted = *weigthed;
    const double one_minus_nugget = 1.0 - nugget;

    // Ottimizzazione: puntatori per accesso continuo alla memoria
    const double *lags_ptr = lags;
    const double *lagt_ptr = lagt;
    const double *data1_ptr = data1;
    const double *data2_ptr = data2;
    const double *mean1_ptr = mean1;
    const double *mean2_ptr = mean2;

    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1_ptr[i];
        const double d2 = data2_ptr[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double current_lag = lags_ptr[i];
            const double current_lagt = lagt_ptr[i];
            
            // Calcola la correlazione
            const double corr = CorFct(cormod, current_lag, current_lagt, par, 0, 0);
            const double adjusted_corr = one_minus_nugget * clamp_corr(corr);
            
            // Calcola la log verosimiglianza bivariata
            const double bl = log_biv_Norm(adjusted_corr, d1, d2, mean1_ptr[i], mean2_ptr[i], sill, 0);
            
            // Calcola i pesi se necessario
            double weights = 1.0;
            if(weighted) {
                weights = GeoCachedPairWeight(i, 1, 1);
            }
            
            total += bl * weights;
        }
    }

    *res = R_FINITE(total) ? total : LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_WrapGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;
    double  u=0.0, w=0.0,weights=1.0;
    double bl=0.0,corr=0.0;
    double alfa=3.0;  double nugget=nuis[0];  double sill =nuis[1];
   if(sill<=0 || nugget<0||nugget>=1){*res=LOW; return;}
   for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                u=data1[i];w=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                    bl=biv_wrapped_log(alfa,u,w,mean1[i],mean2[i],nugget,sill,(1-nugget)*clamp_corr(corr));
                                if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                             *res += weights * bl;
                                }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_T_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                       double *par, int *weigthed, double *res, double *mean1, double *mean2,
                       double *nuis,  int *type_cop, int *cond)
{
    // Estrai valori una sola volta
    const int n_pairs = npairs[0];

    const double sill = nuis[2];
    const double nugget = nuis[1];
    const double df = nuis[0];
    
    // Controllo preliminare dei parametri
    if(sill < 0 || nugget < 0 || nugget >= 1 || df < 0 || df > 0.5) {
        *res = LOW; 
        return;
    }

    double total = 0.0;
    const int weighted = *weigthed;
    const double inv_sqrt_sill = 1.0 / sqrt(sill);
    const double inv_sill = 1.0 / sill;

    // Ottimizzazione: puntatori per accesso continuo alla memoria
    const double *lags_ptr = lags;
    const double *lagt_ptr = lagt;
    const double *data1_ptr = data1;
    const double *data2_ptr = data2;
    const double *mean1_ptr = mean1;
    const double *mean2_ptr = mean2;

    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1_ptr[i];
        const double d2 = data2_ptr[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Calcola le differenze standardizzate
            const double z1 = (d1 - mean1_ptr[i]) * inv_sqrt_sill;
            const double z2 = (d2 - mean2_ptr[i]) * inv_sqrt_sill;
            
            // Calcola la correlazione
            const double corr = CorFct(cormod, lags_ptr[i], lagt_ptr[i], par, 0, 0);
            
            // Calcola i pesi se necessario
            double weights = 1.0;
            if(weighted) {
                weights = GeoCachedPairWeight(i, 1, 1);
            }
            
            // Calcola la densità bivariata t
            const double bt = biv_T(corr, z1, z2, df, nugget);
            total += weights * log(bt * inv_sill);
        }
    }

    *res = R_FINITE(total) ? total : LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Gauss_misp_T_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                  double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                  double *nuis,  int *type_cop, int *cond)
{
    const double df_param = nuis[0], nugget = nuis[1], sill = nuis[2];
    if (!R_FINITE(df_param) || !R_FINITE(nugget) || !R_FINITE(sill) ||
        df_param < 0.0 || df_param >= 0.5 || nugget < 0.0 ||
        nugget >= 1.0 || sill <= 0.0) {
        *res = LOW;
        return;
    }
    const double nu = df_param == 0.0 ? R_PosInf : 1.0 / df_param;
    const double variance = sill / (1.0 - 2.0 * df_param);
    const double log_factor = GeoStudentCorrLogFactor(nu);
    const double scale = 1.0 - nugget;
    double total = 0.0;
    for (int i = 0; i < npairs[0]; ++i) {
        if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
            const double lag = lags[i];
            const double weights = *weigthed ? GeoCachedPairWeight(i, 1, 1) : 1.0;
            if (weights == 0.0) continue;
            const double rho = CorFct(cormod, lag, lagt[i], par, 0, 0);
            const double corr = GeoStudentCorrFromLogFactor(
                rho, scale * rho, nu, log_factor);
            double log_pair = log_biv_Norm(clamp_corr(corr), data1[i], data2[i],
                                           mean1[i], mean2[i], variance, 0);
            total += weights * log_pair;
        }
    }
    *res = R_FINITE(total) ? total : LOW;
}

/******************************************************************************************/
void Comp_Pair_Pois_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                          double *par, int *weigthed, double *res, double *mean1, double *mean2,
                          double *nuis,  int *type_cop, int *cond)
{
    // Pre-carica le variabili globali in variabili locali
    const int npairs_val = npairs[0];
    const int weighted_flag = *weigthed;

    const double nugget = nuis[0];
    
    // Controllo iniziale sul nugget
    if(nugget < 0 || nugget >= 1) {
        *res = LOW; 
        return;
    }
    double total = 0.0;
    double weights = 1.0;
    
    // Loop principale ottimizzato
    for(int i = 0; i < npairs_val; i++) {
        const double u = data1[i];
        const double w = data2[i];
        if(!ISNAN(u) && !ISNAN(w)) {
            // Calcola la correlazione
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0);
            const double mui = exp(mean1[i]);
            const double muj = exp(mean2[i]);
            const int uu = (int)u;
            const int ww = (int)w;
            const double bl = biv_Poisson((1.0 - nugget) * clamp_corr(corr), uu, ww, mui, muj);
            if(weighted_flag) {
                weights = GeoCachedPairWeight(i, 1, 1);
            }
            // Accumula il risultato
            total += bl * weights;
        }
    }
    
    // Assegna il risultato finale
    *res = R_FINITE(total) ? total : LOW;
}

/******************************************************************************************/
void Comp_Pair_PoisGamma_st2mem(int *cormod, double *data1, double *data2,
                                int *N1, int *N2, double *par, int *weigthed,
                                double *res, double *mean1, double *mean2,
                                double *nuis, int *type_cop, int *cond)
{
    (void)N1; (void)N2; (void)type_cop; (void)cond;
    const double nugget = nuis[0], shape = nuis[2];
    if (!R_FINITE(nugget) || nugget < 0.0 || nugget >= 1.0 ||
        !R_FINITE(shape) || !(shape > 0.0) || !npairs || !lags || !lagt) {
        *res = LOW; return;
    }
    GeoKahanSum total; GeoKahanInit(&total); int used = 0;
    for (int i = 0; i < npairs[0]; ++i) {
        int r, t; double mi, mj;
        const int prep = GeoPGPreparePairAt(i, data1[i], data2[i],
                                             mean1[i], mean2[i],
                                             &r, &t, &mi, &mj);
        if (prep == 0) continue;
        if (prep < 0) { *res = LOW; return; }
        const double weight = GeoCachedPairWeight(i, weigthed[0], 1);
        if (!R_FINITE(weight)) { *res = LOW; return; }
        if (weight <= 0.0) continue;
        double base = CorFct(cormod, lags[i], lagt[i], par, 0, 0);
        if (!R_FINITE(base) || fabs(base) > 1.0 + 1e-8) { *res = LOW; return; }
        const int zero = lags[i] <= 0.0 && lagt[i] <= 0.0;
        double rho = zero ? fmax2(-1.0, fmin2(1.0, base))
                          : (1.0 - nugget) * base;
        if (!zero && rho >= 1.0) rho = 1.0 - 1e-15;
        if (!zero && rho <= -1.0) rho = -1.0 + 1e-15;
        const double ll = GeoPGLogJoint(rho, r, t, mi, mj, shape);
        if (!R_FINITE(ll)) { *res = LOW; return; }
        GeoKahanAdd(&total, weight * ll); ++used;
    }
    *res = (used > 0 && R_FINITE(total.sum)) ? total.sum : LOW;
}




/******************************************************************************************/
void Comp_Pair_Gauss_misp_PoisGamma_st2mem(int *cormod, double *data1,
                                           double *data2, int *N1, int *N2,
                                           double *par, int *weigthed,
                                           double *res, double *mean1,
                                           double *mean2, double *nuis,
                                           int *type_cop, int *cond)
{
    (void)N1; (void)N2; (void)type_cop; (void)cond;
    const double nugget = nuis[0], shape = nuis[2];
    if (!R_FINITE(nugget) || nugget < 0.0 || nugget >= 1.0 ||
        !R_FINITE(shape) || !(shape > 0.0) || !npairs || !lags || !lagt) {
        *res = LOW; return;
    }
    GeoKahanSum total; GeoKahanInit(&total); int used = 0;
    for (int i = 0; i < npairs[0]; ++i) {
        if (ISNAN(data1[i]) || ISNAN(data2[i])) continue;
        if (!R_FINITE(data1[i]) || !R_FINITE(data2[i]) ||
            !R_FINITE(mean1[i]) || !R_FINITE(mean2[i]) ||
            mean1[i] > 700.0 || mean2[i] > 700.0) { *res = LOW; return; }
        const double mi = exp(mean1[i]), mj = exp(mean2[i]);
        const double weight = GeoCachedPairWeight(i, weigthed[0], 1);
        if (!R_FINITE(weight)) { *res = LOW; return; }
        if (weight <= 0.0) continue;
        double base = CorFct(cormod, lags[i], lagt[i], par, 0, 0);
        if (!R_FINITE(base) || fabs(base) > 1.0 + 1e-8) { *res = LOW; return; }
        const int zero = lags[i] <= 0.0 && lagt[i] <= 0.0;
        double rho = zero ? fmax2(-1.0, fmin2(1.0, base))
                          : (1.0 - nugget) * base;
        if (rho >= 1.0) rho = 1.0 - 1e-15;
        if (rho <= -1.0) rho = -1.0 + 1e-15;
        const double ll = GeoPGGaussLogJoint(rho, data1[i], data2[i],
                                              mi, mj, shape);
        if (!R_FINITE(ll)) { *res = LOW; return; }
        GeoKahanAdd(&total, weight * ll); ++used;
    }
    *res = (used > 0 && R_FINITE(total.sum)) ? total.sum : LOW;
}


/******************************************************************************************/
void Comp_Pair_PoisZIP_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)

{
    int i=0,uu,ww;
     double weights=1.0,corr,mui,muj,bl,u=0.0, w=0.0;


       double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];


      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
    // Computes the log-likelihood:
  for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                          u=data1[i];      w=data2[i];
                     corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                     mui=exp(mean1[i]);
                     muj=exp(mean2[i]);
                          uu=(int) u;  ww=(int) w;
                      bl=biv_PoissonZIP(clamp_corr(corr),uu,ww,mui, muj,mup,nugget1,nugget2);
               if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                       *res += weights * Geo_safe_log_prob(bl);

                                    }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Gauss_misp_Pois_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0,N=2;
     double weights=1.0,corr,corr2,mui,muj,bl,u=0.0, w=0.0;
   double nugget=nuis[0];
     if(nugget<0||nugget>=1){*res=LOW; return;}
// Fixed-size workspace kept on the stack.
    double M0[2] = {0.0, 0.0};
    double M1[2] = {0.0, 0.0};
    double *M[2] = {M0, M1};
    double dat[2] = {0.0, 0.0};
    for(i=0;i<npairs[0];i++){
       u=data1[i];w=data2[i];
if(!ISNAN(u)&&!ISNAN(w) ){
                            corr=CorFct(cormod,lags[i],lagt[i],par,0,0)*(1-nugget);
                            mui=exp(mean1[i]);muj=exp(mean2[i]);
                            corr2=corr_pois(clamp_corr(corr),mui, muj);

                            M[0][0]=mui; M[1][1]=muj;M[0][1]=sqrt(mui*muj)*corr2;M[1][0]= M[0][1];
                           dat[0]=u-mui;dat[1]=w-muj;
                          if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                              bl=log_dNnorm(N,M,dat);
                              *res += weights * bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Gauss_misp_PoisZIP_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;
     double weights=1.0,corr,mui,muj,bl;
     double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];


      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
    // Computes the log-likelihood:
  for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                     corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                     mui=exp(mean1[i]);
                     muj=exp(mean2[i]);

               bl=biv_Mis_PoissonZIP(clamp_corr(corr),data1[i],data2[i],mui, muj,mup,nugget1,nugget2);
               if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                       *res += weights * Geo_safe_log_prob(bl);

                                    }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Tukeyh_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{


    int i=0;
    double corr,zi,zj,weights=1.0,bl;
    double sill=nuis[1];
    double nugget=nuis[0];
    double tail=nuis[2];
      if( sill<0||tail<0||tail>0.5||nugget<0||nugget>=1){*res=LOW; return;}

for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);

                                   if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);

 bl=biv_tukey_h((1-nugget)*clamp_corr(corr),zi,zj,mean1[i],mean2[i],tail,sill);
                             *res += weights * Geo_safe_log_prob(bl);

                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Tukeyhh_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{


    int i=0;
    double corr,zi,zj,weights=1.0,bl;
double sill=nuis[1];
    double nugget=nuis[0];
    double h1=nuis[3];
    double h2=nuis[2];
      if( sill<=0||h1<0||h1>=0.5||h2<0||h2>=0.5||nugget<0||nugget>=1){*res=LOW; return;}

for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
            if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
bl=biv_tukey_hh((1-nugget)*clamp_corr(corr),zi,zj,mean1[i],mean2[i],sill,h1,h2);
                             *res += weights * Geo_safe_log_prob(bl);
                         }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void  Comp_Pair_TWOPIECEGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
       int i=0;
    double corr,zi,zj,weights=1.0,p11,eta,qq,sill,bl,nugget;
   eta=nuis[2];  //skewness parameter
    sill=nuis[1];
    nugget=nuis[0];

     if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1) {*res=LOW;  return;}
  qq=qnorm((1-eta)/2,0,1,1,0);

for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                   if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                        p11=pbnorm22(qq,qq,clamp_corr(corr));
  bl=biv_two_pieceGaussian((1-nugget)*clamp_corr(corr),zi,zj,sill,eta,p11,mean1[i],mean2[i]);

                           *res += weights * Geo_safe_log_prob(bl);
                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}

/******************************************************************************************/
void  Comp_Pair_TWOPIECETukeyh_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
       int i=0;
    double corr,zi,zj,weights=1.0,p11,eta,qq,sill,bl,nugget,tail;
     eta  = nuis[2];  //skewness parameter
    tail = nuis[3];  //tail parameter
    sill =nuis[1];
    nugget=nuis[0];

    if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1||tail<0||tail>0.5) {*res=LOW;  return;}


     qq=qnorm((1-eta)/2,0,1,1,0);
      //   if( fabs(eta)>1  || tail<=0) {*res=LOW;  return;}

for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                  if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                      
                        p11=pbnorm22(qq,qq,clamp_corr(corr));
                         bl=biv_two_pieceTukeyh((1-nugget)*clamp_corr(corr),zi,zj,sill,eta,tail,p11,mean1[i],mean2[i]);

                           *res += weights * Geo_safe_log_prob(bl);

                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}

/******************************************************************************************/
void  Comp_Pair_TWOPIECET_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
       int i=0;
    double corr,zi,zj,weights=1.0,p11,qq,bl;
   double eta=nuis[3];  //skewness parameter
    double df=nuis[0];
    double sill=nuis[2];
    double nugget=nuis[1];
    if(sill<0||nugget<0||nugget>=1 ||fabs(eta)>1|| df >0.5||df<0) {*res=LOW;  return;}
       qq=qnorm((1-eta)/2,0,1,1,0);
       //  if( fabs(eta)>1|| sill<0||df >0.5||df<0) {*res=LOW;  return;}


for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                         if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                        p11=pbnorm22(qq,qq,clamp_corr(corr));
                         bl=biv_two_pieceT(clamp_corr(corr),zi,zj,sill,df,eta,p11,mean1[i],mean2[i],nugget);

                           *res += weights * Geo_safe_log_prob(bl);

                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}

/******************************************************************************************/
void Comp_Pair_PoisbinGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu=0,ww=0;
    double dens=0.0,weights=1.0,u,w,corr=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
     double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}

     for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
              corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                p11=pbnorm22(mean1[i],mean2[i],(1-nugget)*clamp_corr(corr));
                p1=pnorm(mean1[i],0,1,1,0);
                p2=pnorm(mean2[i],0,1,1,0);
                u=data1[i];w=data2[i];
                                     uu=(int) u;
                                     ww=(int) w;
                                    if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                                    dens=biv_poisbin (N1[0],uu,ww,p1,p2,p11);
                                 *res+=log(dens)*weights;
                               }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_PoisbinnegGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0,uu=0,ww=0;
    double dens=0.0,weights=1.0,u,w,corr;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
 double nugget=nuis[0];
      if(nugget<0||nugget>=1){*res=LOW; return;}

           for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                 p11=pbnorm22(mean1[i],mean2[i],(1-nugget)*clamp_corr(corr));
                 p1=pnorm((mean1[i]),0,1,1,0);
                 p2=pnorm((mean2[i]),0,1,1,0);
                                u=data1[i];w=data2[i];

                                     uu=(int) u; ww=(int) w;
                                    if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                                    dens=biv_poisbinneg (N1[0],uu,ww,p1,p2,p11);
                                 *res+=log(dens)*weights;}}

    if(!R_FINITE(*res))*res = LOW;
    return;
}


void Comp_Pair_BinomGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu=0,ww=0;
    double bl=0.0,weights=1.0,u=0.0,w=0.0,a=0.0,b=0.0,corr=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success

      double nugget=nuis[0];
      if(nugget<0||nugget>=1){*res=LOW; return;}
           for(i=0;i<npairs[0];i++){
                   u=data1[i];w=data2[i];
             if(!ISNAN(u)&&!ISNAN(w) ){
                            corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                            a=mean1[i];b=mean2[i];
                            psj=pbnorm22(a,b,(1-nugget)*clamp_corr(corr));
                            p1=pnorm(a,0,1,1,0);
                            p2=pnorm(b,0,1,1,0);
                            uu=(int) u;  ww=(int) w;
                          if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                                    bl=biv_binom (N1[0],uu,ww,p1,p2,psj);
                                 *res += weights * Geo_safe_log_prob(bl);
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}
void Comp_Pair_BinomLogi_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu=0,ww=0;
    double bl=0.0,weights=1.0,u=0.0,w=0.0,a=0.0,b=0.0,corr=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success

      double nugget=nuis[0];
      if(nugget<0||nugget>=1){*res=LOW; return;}
           for(i=0;i<npairs[0];i++){
                   u=data1[i];w=data2[i];
             if(!ISNAN(u)&&!ISNAN(w) ){
                            corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                            a=mean1[i];b=mean2[i];
                            psj=pblogi22(a,b,(1-nugget)*clamp_corr(corr));
                            p1=1/(1+exp(-a));p2=1/(1+exp(-b));
                            uu=(int) u;  ww=(int) w;
                          if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                                    bl=biv_binom (N1[0],uu,ww,p1,p2,psj);
                                 *res += weights * Geo_safe_log_prob(bl);
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}
void Comp_Pair_BinomNNGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu=0,ww=0,n1,n2;
    double bl=0.0,weights=1.0,u=0.0,w=0.0,a=0.0,b=0.0,corr=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success

      double nugget=nuis[0];
      if(nugget<0||nugget>=1){*res=LOW; return;}

           for(i=0;i<npairs[0];i++){
                   u=data1[i];w=data2[i];
             if(!ISNAN(u)&&!ISNAN(w) ){
                            corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                            a=mean1[i];b=mean2[i];
                            psj=pbnorm22(a,b,(1-nugget)*clamp_corr(corr));
                            p1=pnorm(a,0,1,1,0);
                            p2=pnorm(b,0,1,1,0);
                            uu=(int) u;  ww=(int) w;
                            n1=N1[i];n2=N2[i];
                          if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                                    bl=biv_binom222(n1,n2,uu,ww,p1,p2,psj);
                                 *res += weights * Geo_safe_log_prob(bl);
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}

void Comp_Pair_BinomNNLogi_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu=0,ww=0,n1,n2;
    double bl=0.0,weights=1.0,u=0.0,w=0.0,a=0.0,b=0.0,corr=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success

      double nugget=nuis[0];
      if(nugget<0||nugget>=1){*res=LOW; return;}

           for(i=0;i<npairs[0];i++){
                   u=data1[i];w=data2[i];
             if(!ISNAN(u)&&!ISNAN(w) ){
                            corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                            a=mean1[i];b=mean2[i];
                            psj=pblogi22(a,b,(1-nugget)*clamp_corr(corr));
                            p1=1/(1+exp(-a));p2=1/(1+exp(-b));
                            uu=(int) u;  ww=(int) w;
                            n1=N1[i];n2=N2[i];
                          if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                                    bl=biv_binom222(n1,n2,uu,ww,p1,p2,psj);
                                 *res += weights * Geo_safe_log_prob(bl);
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}

void Comp_Pair_BinomnegGauss_st2mem(int *cormod, double *data1, double *data2,
    int *N1, int *N2, double *par, int *weigthed, double *res,
    double *mean1, double *mean2, double *nuis, int *type_cop, int *cond)
{
    *res = 0.0;
 
    double nugget = nuis[0];
    if (nugget < 0.0 || nugget >= 1.0) { *res = LOW; return; }
 
    int N = N1[0];
    if (N <= 0) { *res = LOW; return; }
 
    /* FIX 1: scansione preliminare per dimensionare la cache */
    int MAX_UV = 0;
    for (int k = 0; k < npairs[0]; k++) {
        if (!ISNAN(data1[k]) && data1[k] > MAX_UV) MAX_UV = (int)data1[k];
        if (!ISNAN(data2[k]) && data2[k] > MAX_UV) MAX_UV = (int)data2[k];
    }
    /* FIX 3: dimensione sicura */
    int cache_size = N + 2 * MAX_UV + 10;
    double *lgamma_cache = (double *) R_alloc((size_t) cache_size, sizeof(double));
    for (int k = 0; k < cache_size; k++)
        lgamma_cache[k] = lgammafn(k + 1);
 
    int    i, uu, ww;
    double bl, weights, u, w, a, b, corr, p1, p2, psj;
 
#define P_CLAMP 1e-10
 
    for (i = 0; i < npairs[0]; i++) {
 
        u = data1[i];
        w = data2[i];
        if (ISNAN(u) || ISNAN(w)) continue;
 
        weights = 1.0;
        if (*weigthed)
            weights = GeoCachedPairWeight(i, 1, 1);
        if (weights <= 0.0) continue;
 
        corr = clamp_corr(CorFct(cormod, lags[i], lagt[i], par, 0, 0));
        a    = mean1[i];
        b    = mean2[i];
        psj  = pbnorm22(a, b, (1 - nugget) * corr);
        p1   = fmax(P_CLAMP, fmin(1.0 - P_CLAMP, pnorm(a, 0, 1, 1, 0)));
        p2   = fmax(P_CLAMP, fmin(1.0 - P_CLAMP, pnorm(b, 0, 1, 1, 0)));
 
        /* FIX 4: evita psj == min(p1,p2) esattamente */
        psj  = fmax(fmax(0.0, p1 + p2 - 1.0),
                    fmin(fmin(p1, p2) - P_CLAMP, psj));
 
        uu = (int)u;
        ww = (int)w;
 
        /* FIX 1: passa la cache allocata una sola volta */
        bl = biv_binomneg(N, uu, ww, p1, p2, psj, lgamma_cache, cache_size);
 
        if (!R_FINITE(bl) || bl <= 0.0) {
            *res += weights * LOG_FLOOR;
            continue;
        }
        *res += weights * log(bl);
    }
 
 
#undef P_CLAMP
 
    /* FIX 2: soglia robusta per npairs piccolo */

    if (!R_FINITE(*res)) *res = LOW;
}
 
void Comp_Pair_BinomnegLogi_st2mem(int *cormod, double *data1, double *data2,
    int *N1, int *N2, double *par, int *weigthed, double *res,
    double *mean1, double *mean2, double *nuis, int *type_cop, int *cond)
{
    *res = 0.0;
 
    double nugget = nuis[0];
    if (nugget < 0.0 || nugget >= 1.0) { *res = LOW; return; }
 
    int N = N1[0];
    if (N <= 0) { *res = LOW; return; }
 
    /* FIX 1: scansione preliminare per dimensionare la cache */
    int MAX_UV = 0;
    for (int k = 0; k < npairs[0]; k++) {
        if (!ISNAN(data1[k]) && data1[k] > MAX_UV) MAX_UV = (int)data1[k];
        if (!ISNAN(data2[k]) && data2[k] > MAX_UV) MAX_UV = (int)data2[k];
    }
    /* FIX 3: dimensione sicura */
    int cache_size = N + 2 * MAX_UV + 10;
    double *lgamma_cache = (double *) R_alloc((size_t) cache_size, sizeof(double));
    for (int k = 0; k < cache_size; k++)
        lgamma_cache[k] = lgammafn(k + 1);
 
    int    i, uu, ww;
    double bl, weights, u, w, a, b, corr, p1, p2, psj;
 
    for (i = 0; i < npairs[0]; i++) {
 
        u = data1[i];
        w = data2[i];
        if (ISNAN(u) || ISNAN(w)) continue;
 
        weights = 1.0;
        if (*weigthed)
            weights = GeoCachedPairWeight(i, 1, 1);
        if (weights <= 0.0) continue;
 
        corr = clamp_corr(CorFct(cormod, lags[i], lagt[i], par, 0, 0));
        a    = mean1[i];
        b    = mean2[i];
        psj  = pblogi22(a, b, (1 - nugget) * corr);
        p1   = 1.0 / (1.0 + exp(-a));
        p2   = 1.0 / (1.0 + exp(-b));
        uu   = (int)u;
        ww   = (int)w;
 
        /* FIX 1: passa la cache allocata una sola volta */
        bl = biv_binomneg(N, uu, ww, p1, p2, psj, lgamma_cache, cache_size);
 
        if (!R_FINITE(bl) || bl <= 0.0) {
            *res += weights * LOG_FLOOR;
            continue;
        }
        *res += weights * log(bl);
    }
 
 
    if (!R_FINITE(*res)) *res = LOW;
}

void Comp_Pair_BinomnegGaussZINB_st2mem(int *cormod, double *data1, double *data2,
    int *N1, int *N2, double *par, int *weigthed, double *res,
    double *mean1, double *mean2, double *nuis, int *type_cop, int *cond)
{
    *res = 0.0;
 
    double nugget1 = nuis[0];
    double nugget2 = nuis[1];
    double mup     = nuis[2];
 
    if (nugget1 < 0.0 || nugget1 >= 1.0 ||
        nugget2 < 0.0 || nugget2 >= 1.0) { *res = LOW; return; }
 
    int N = N1[0];
    if (N <= 0) { *res = LOW; return; }
 
    /* FIX 1: scansione preliminare per dimensionare la cache */
    int MAX_UV = 0;
    for (int k = 0; k < npairs[0]; k++) {
        if (!ISNAN(data1[k]) && data1[k] > MAX_UV) MAX_UV = (int)data1[k];
        if (!ISNAN(data2[k]) && data2[k] > MAX_UV) MAX_UV = (int)data2[k];
    }
    /* FIX 3: dimensione sicura */
    int cache_size = N + 2 * MAX_UV + 10;
    double *lgamma_cache = (double *) R_alloc((size_t) cache_size, sizeof(double));
    for (int k = 0; k < cache_size; k++)
        lgamma_cache[k] = lgammafn(k + 1);
 
    int    i, uu, ww;
    double bl, weights, u, w, a, b, corr;
 
    for (i = 0; i < npairs[0]; i++) {
 
        u = data1[i];
        w = data2[i];
        if (ISNAN(u) || ISNAN(w)) continue;
 
        weights = 1.0;
        if (*weigthed)
            weights = GeoCachedPairWeight(i, 1, 1);
        if (weights <= 0.0) continue;
 
        corr = clamp_corr(CorFct(cormod, lags[i], lagt[i], par, 0, 0));
        a    = mean1[i];
        b    = mean2[i];
        uu   = (int)u;
        ww   = (int)w;
 
        /* FIX 1: passa la cache allocata una sola volta */
        bl = biv_binomnegZINB(N, corr, uu, ww, a, b,
                              nugget1, nugget2, mup,
                              lgamma_cache, cache_size);
 
        if (!R_FINITE(bl) || bl <= 0.0) {
            *res += weights * LOG_FLOOR;
            continue;
        }
        *res += weights * log(bl);
    }
 
 
  
    if (!R_FINITE(*res)) *res = LOW;
}

/*********************************************************************************/

// Composite marginal (difference) log-likelihood for the spatial-temporal Gaussian model:
/******************************************************************************************/
void Comp_Pair_SkewGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i=0;
    double bl,corr,zi,zj,weights=1.0;
    double nugget=nuis[0];
    double sill=nuis[1];
     if(nugget<0|| nugget>=1||sill<0){*res=LOW;  return;}

       for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                                    bl=biv_skew(clamp_corr(corr),zi,zj,mean1[i],mean2[i],nuis[1],nuis[2],nugget);
                             *res += weights * Geo_safe_log_prob(bl);
                         }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}


void Comp_Pair_SkewLaplace_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                              double *par, int *weigthed, double *res, double *mean1, double *mean2,
                              double *nuis,  int *type_cop, int *cond)
{
    // Controllo precoce dei parametri
    const double nugget = nuis[0];
      const double skew=nuis[2];

    
    if(nugget < 0 || nugget >= 1 || skew < 0 || skew >= 1) {
        *res = LOW;
        return;
    }

    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];

    const double scale = 1.0 - nugget;

  
    const double sill = nuis[1];
    double total = 0.0; 

    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Calcolo correlazione
            const double corr = CorFct(cormod,lags[i], lagt[i], par, 0, 0);
            // Calcolo pesi se necessario
            const double weights = weighted ? (GeoCachedPairWeight(i, 1, 1)) : 1.0;
 
           const double val=log_biv_skewlaplace(scale*clamp_corr(corr),d1,d2,mean1[i],mean2[i],sill, skew);
            total += weights * val;
        }
    }
    *res = R_FINITE(total) ? total : LOW;
}



/******************************************************************************************/
void Comp_Pair_SinhGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i=0;
    double bl,corr,zi,zj,weights=1.0;
       if(nuis[3]<0||nuis[1]<0||nuis[0]<0||nuis[0]>=1) {*res=LOW;  return;}

       for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){

                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0)*(1-nuis[0]);
                                    if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                                    bl=biv_sinh(clamp_corr(corr),zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],nuis[1]);

                             *res += weights * Geo_safe_log_prob(bl);
                         }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Gamma_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
{

    int i=0;
    double corr,zi,zj,weights=1.0,bl=1.0;
        double nugget=nuis[0];
     if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
    for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                            if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                                  bl=biv_gamma((1-nugget)*clamp_corr(corr),zi,zj,mean1[i],mean2[i],nuis[2]);
                                   if (!R_FINITE(bl)) bl = LOG_FLOOR;
                                    *res+= weights*bl;
                         }}}}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Kumaraswamy_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{


    int i=0;
    double corr,zi,zj,weights=1.0,bl;
    //double sill=nuis[1];
    double nugget=nuis[0];
    double min=nuis[4];
     double max=nuis[5];
  if(nuis[2]<0||nuis[3]<0||min>max)  {*res=LOW;  return;}
   for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                    if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                                   bl= biv_Kumara((1-nugget)*clamp_corr(corr),zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);
                             *res += weights * Geo_safe_log_prob(bl);
                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}

/******************************************************************************************/
void Comp_Pair_Kumaraswamy2_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{


    int i=0;
    double corr,zi,zj,weights=1.0,bl;
    //double sill=nuis[1];
    double nugget=nuis[0];
    double min=nuis[4];
     double max=nuis[5];
  if(nuis[2]<0||nuis[3]<0||min>max)  {*res=LOW;  return;}
   for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                    if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                                   bl= biv_Kumara2((1-nugget)*clamp_corr(corr),zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);
                             *res += weights * Geo_safe_log_prob(bl);
                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Beta_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;
    double corr,zi,zj,weights=1.0,bl;
    //double sill=nuis[1];
    double nugget=nuis[0];
    double min=nuis[4];
     double max=nuis[5];
  if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}
   for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                    if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                                   bl= biv_beta((1-nugget)*clamp_corr(corr),zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);
                             *res += weights * Geo_safe_log_prob(bl);
                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Weibull_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i=0;
    double corr,zi,zj,weights=1.0,bl=0.0;
  double nugget=nuis[0];
     if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
       for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                   if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                                      bl=biv_Weibull((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]);
                                      *res+= weights*bl;
                         }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_LogGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i=0;
    double corr,zi,zj,weights=1.0,bl=0.0;
     double sill=nuis[1];double nugget=nuis[0];

    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}
     for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                              if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                  bl=d2lognorm(zi,zj,sill,nugget, mean1[i], mean2[i],clamp_corr(corr));
                             *res += weights * Geo_safe_log_prob(bl);
                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Pair_LogLogistic_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i=0;
    double bl,corr,zi,zj,weights=1.0,nugget=0.0;
      nugget=nuis[0];
     if(nugget<0||nugget>=1||nuis[2]<=2) {*res=LOW;  return;}

       for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                            zi=data1[i];zj=data2[i];

                                corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                                bl=biv_LogLogistic((1-nugget)*clamp_corr(corr),zi,zj,mean1[i],mean2[i],nuis[2]);
                                *res += weights * Geo_safe_log_prob(bl);
                            }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/

void Comp_Pair_Logistic_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i=0;
    double bl,corr=0.0,zi=0.0,zj=0.0,weights=1.0;
    if( nuis[1]<=0)  {*res=LOW;  return;}

    double sill=1-nuis[0];

          for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                            zi=data1[i];zj=data2[i];

                                corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                                bl=biv_Logistic(sill*clamp_corr(corr),zi,zj,mean1[i],mean2[i],nuis[1]);
                                *res += weights * Geo_safe_log_prob(bl);
                            }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
void Comp_Pair_TWOPIECEBIMODAL_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i=0;
    double p11,qq,bl,corr,zi,zj,weights=1.0;
   double eta=nuis[4];  //skewness parameter
   double delta=nuis[3];
   double sill=nuis[2];
   double nugget=nuis[1];
   double df=nuis[0];
 if( fabs(eta)>1||df<0||nugget<0||nugget>=1||delta<0||sill<0) {*res=LOW;  return;}

           qq=qnorm((1-eta)/2,0,1,1,0);

          for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                            zi=data1[i];zj=data2[i];
                                corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                p11=pbnorm22(qq,qq,clamp_corr(corr));
                                if(*weigthed) weights=GeoCachedPairWeight(i, 1, 1);
                                bl=biv_two_piece_bimodal((1-nugget)*clamp_corr(corr),zi,zj,sill,df,delta,eta,p11,mean1[i],mean2[i]);
                             *res += weights * Geo_safe_log_prob(bl);
                         }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/********************* BIVARIATE CASE *****************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/

/* pairwise for bivariate GRF*/
void Comp_Pair_Gauss_biv2mem(int *cormod, double *data1,double *data2,int *NN,
    double *par, int *weigthed,double *res,double *mean1,double *mean2,
    double *nuis)
{
    int i=0;

    double  dens=0.0,weights=1.0;
    if(  par[0]<0|| par[1]<0|| par[2]<0|| par[3]<0) {*res=LOW;  return;}
 for(i=0;i<npairs[0];i++){
  if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
          dens=log_biv2gauss(cormod,lags_1[i],par, data1[i]-mean1[i], data2[i]-mean2[i],first_1[i], second_1[i]);
          *res+= dens*weights;
                                }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}


/* pairwise  skew for bivariate GRF*/
void Comp_Pair_SkewGauss_biv2mem(int *cormod, double *data1,double *data2,int *NN,
    double *par, int *weigthed,double *res,double *mean1,double *mean2,
    double *nuis)

{
    int i=0;
    double u=0.0, w=0.0, rhotv=0.0,weights=1.0;
      double vari[2] = {par[0], par[1]};  /// variances of the skew gaussian
    par[0]=1;par[1]=1;
    if(vari[0]<0||vari[1]<0){*res=LOW; return;}
      weights=1;

        for(i=0;i<npairs[0];i++){
                             u=data1[i];
                             w=data2[i];
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                    rhotv=CorFct(cormod,lags_1[i],0,par,first_1[i],second_1[i]);
                     *res+= log(biv_skew2(rhotv,u,w,vari[first_1[i]],vari[second_1[i]],1,nuis[first_1[i]],nuis[second_1[i]]))*weights;
                                }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
/******************************************************************************************/
/********************* SPATIAL Gaussian and Clayton and skegaussian COPULA *****************************************************/
/******************************************************************************************/
/******************************************************************************************/

// Composite marginal (pariwise) log-likelihood for the spatial Gaussian model:
void Comp_Pair_GaussCop2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                            double *par, int *weigthed, double *res, double *mean1, double *mean2,
                            double *nuis,  int *type_cop, int *cond)
{
    // === 1. Initialization and Validation ===
    const double sill = nuis[1];
    const double nugget = nuis[0];
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;

    const int model = 1;
    const int type_cop_val = type_cop[0];
    const int cond_val = cond[0];

    // Validation of input parameters
    if (!R_FINITE(sill) || sill <= 0.0 || nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }

    double total = 0.0;

    // === 2. Main Computation ===
    if (weigthed_val) {
        // Loop with weights calculation
        for (int i = 0; i < npairs_val; i++) {
            const double d1 = data1[i];
            const double d2 = data2[i];
            if (!ISNAN(d1) && !ISNAN(d2)) {
                const double lag_i = lags[i];
                const double corr = CorFct(cormod, lag_i, 0, par, 0, 0); // Pre-calculate correlation once
                const double weights = GeoCachedPairWeight(i, 1, 0);  // Calculate weights

                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          d1, d2, mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl * weights;  // Update total
            }
        }
    } else {
        // Loop without weights calculation
        for (int i = 0; i < npairs_val; i++) {
            const double d1 = data1[i];
            const double d2 = data2[i];

            if (!ISNAN(d1) && !ISNAN(d2)) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0); // Pre-calculate correlation once
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          d1, d2, mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl;  // Update total
            }
        }
    }

    // === 3. Finalization ===
    *res = R_FINITE(total) ? total : LOW;
}

// Composite marginal (pariwise) log-likelihood for the spatial Gaussian model:
void Comp_Pair_TCop2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                        double *par, int *weigthed, double *res, double *mean1, double *mean2,
                        double *nuis,  int *type_cop, int *cond)
{
    int model = 12;
    const double df = nuis[0];
    const double nugget = nuis[1];
    const double sill = nuis[2];
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;

    const int type_cop_val = type_cop[0];
    const int cond_val = cond[0];

    if (!R_FINITE(df) || df <= 0.0 || !R_FINITE(sill) || sill <= 0.0 ||
        nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }

    double total = 0.0;

    if (weigthed_val) {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double lag = lags[i];
                const double corr = CorFct(cormod, lag, 0, par, 0, 0);
                const double weights = GeoCachedPairWeight(i, 1, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val, data1[i], data2[i],
                                          mean1[i], mean2[i], nuis, model, N1[i], N2[i]);
                total += bl * weights;
            }
        }
    } else {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val, data1[i], data2[i],
                                          mean1[i], mean2[i], nuis, model, N1[i], N2[i]);
                total += bl;
            }
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}

void Comp_Pair_BetaCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
     int model=28; 
    /*############*/
    int i=0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill <= 0 || nugget < 0 || nugget >= 1 ||
       !R_FINITE(nuis[2]) || !R_FINITE(nuis[3]) || nuis[2] <= 0 || nuis[3] <= 0 ||
       !R_FINITE(nuis[4]) || !R_FINITE(nuis[5]) || nuis[4] >= nuis[5]) { *res=LOW; return; }

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=GeoCachedPairWeight(i, 1, 0);
           bl=biv_cop(clamp_corr(corr),type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]) ;
                        *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Pair_Beta2Cop2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                            double *par, int *weigthed, double *res, double *mean1, double *mean2,
                            double *nuis,  int *type_cop, int *cond)
{
    int model = 50;
    const double sill = nuis[1];
    const double nugget = nuis[0];
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;

    const int type_cop_val = type_cop[0];
    const int cond_val = cond[0];

    if (sill <= 0 || nugget < 0 || nugget >= 1 || !R_FINITE(nuis[2]) || nuis[2] <= 0 ||
        !R_FINITE(nuis[3]) || !R_FINITE(nuis[4]) || nuis[3] >= nuis[4]) {
        *res = LOW;
        return;
    }

    double total = 0.0;

    if (weigthed_val) {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double lag = lags[i];
                const double corr = CorFct(cormod, lag, 0, par, 0, 0);
                const double weights = GeoCachedPairWeight(i, 1, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl * weights;
            }
        }
    } else {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl;
            }
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}
void Comp_Pair_KumaraswamyCop2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                  double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                  double *nuis,  int *type_cop, int *cond)
{
    int model = 33;
    const double sill = nuis[1];
    const double nugget = nuis[0];
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;

    const int type_cop_val = type_cop[0];
    const int cond_val = cond[0];

    if (sill <= 0 || nugget < 0 || nugget >= 1 ||
        !R_FINITE(nuis[2]) || !R_FINITE(nuis[3]) || nuis[2] <= 0 || nuis[3] <= 0 ||
        !R_FINITE(nuis[4]) || !R_FINITE(nuis[5]) || nuis[4] >= nuis[5]) {
        *res = LOW;
        return;
    }

    double total = 0.0;

    if (weigthed_val) {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double lag = lags[i];
                const double corr = CorFct(cormod, lag, 0, par, 0, 0);
                const double weights = GeoCachedPairWeight(i, 1, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl * weights;
            }
        }
    } else {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl;
            }
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}

void Comp_Pair_Kumaraswamy2Cop2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                   double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                   double *nuis, int *type_cop, int *cond)
{
    const int model = 42;
    const double nugget = nuis[0];
    const double shape = nuis[2];
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;

    const int type_cop_val = type_cop[0];
    const int cond_val = cond[0];

    if (nugget < 0.0 || nugget >= 1.0 || !R_FINITE(shape) || shape <= 0.0 ||
        !R_FINITE(nuis[3]) || !R_FINITE(nuis[4]) || nuis[3] >= nuis[4]) {
        *res = LOW;
        return;
    }

    double total = 0.0;
    for (int i = 0; i < npairs_val; ++i) {
        if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
            const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
            const double weights = weigthed_val ? GeoCachedPairWeight(i, 1, 0) : 1.0;
            const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                      data1[i], data2[i], mean1[i], mean2[i],
                                      nuis, model, N1[i], N2[i]);
            total += bl * weights;
        }
    }
    *res = R_FINITE(total) ? total : LOW;
}

void Comp_Pair_WeibullCop2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                              double *par, int *weigthed, double *res, double *mean1, double *mean2,
                              double *nuis,  int *type_cop, int *cond)
{
    const int model = 26;
    const double nugget = nuis[0];
    const double shape = nuis[2];
    const int npairs_val = npairs[0];

    const int weigthed_val = *weigthed;
    const int type_cop_val = type_cop[0];
    const int cond_val = cond[0];

    if (nugget < 0 || nugget >= 1 || !R_FINITE(shape) || shape <= 0.0) {
        *res = LOW;
        return;
    }

    double total = 0.0;
    double weights = 1.0;

    if (weigthed_val) {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                weights = GeoCachedPairWeight(i, 1, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl * weights;
            }
        }
    } else {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl;
            }
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}


void Comp_Pair_GammaCop2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                            double *par, int *weigthed, double *res, double *mean1, double *mean2,
                            double *nuis,  int *type_cop, int *cond)
{
    const int model = 21;
    const double nugget = nuis[0];
    const int npairs_val = npairs[0];

    const int weigthed_val = *weigthed;
    const int type_cop_val = type_cop[0];
    const int cond_val = cond[0];

    if (nugget < 0 || nugget >= 1 || !R_FINITE(nuis[2]) || nuis[2] <= 0.0) {
        *res = LOW;
        return;
    }

    double total = 0.0;
    double weights = 1.0;

    if (weigthed_val) {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                weights = GeoCachedPairWeight(i, 1, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl * weights;
            }
        }
    } else {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl;
            }
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}

void Comp_Pair_LogGaussCop2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                               double *par, int *weigthed, double *res, double *mean1, double *mean2,
                               double *nuis,  int *type_cop, int *cond)
{
    const int model = 22;
    const double sill = nuis[1];
    const double nugget = nuis[0];
    const int npairs_val = npairs[0];

    const int weigthed_val = *weigthed;
    const int type_cop_val = type_cop[0];
    const int cond_val = cond[0];

    if (sill <= 0 || nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }

    double total = 0.0;
    double weights = 1.0;

    if (weigthed_val) {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                weights = GeoCachedPairWeight(i, 1, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl * weights;
            }
        }
    } else {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl;
            }
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}

void Comp_Pair_SkewLaplaceCop2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                               double *par, int *weigthed, double *res, double *mean1, double *mean2,
                               double *nuis,  int *type_cop, int *cond)
{
    const int model = 59;
    const double nugget = nuis[0];
    const int npairs_val = npairs[0];

    const int weigthed_val = *weigthed;
    const int type_cop_val = type_cop[0];
    const int cond_val = cond[0];

    if (nuis[1] <= 0 || nugget < 0 || nugget >= 1 ||
        !R_FINITE(nuis[2]) || nuis[2] <= 0.0 || nuis[2] >= 1.0) {
        *res = LOW;
        return;
    }

    double total = 0.0;
    double weights = 1.0;

    if (weigthed_val) {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                weights = GeoCachedPairWeight(i, 1, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl * weights;
            }
        }
    } else {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl;
            }
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}

void Comp_Pair_PoisCop2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                           double *par, int *weigthed, double *res, double *mean1, double *mean2,
                           double *nuis,  int *type_cop, int *cond)
{
    const int model = 30;
    const double nugget = nuis[0];
    const int npairs_val = npairs[0];

    const int weigthed_val = *weigthed;
    const int type_cop_val = type_cop[0];
    const int cond_val = cond[0];

    if (nugget < 0 || nugget > 1) {
        *res = LOW;
        return;
    }

    double total = 0.0;
    double weights = 1.0;

    if (weigthed_val) {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                weights = GeoCachedPairWeight(i, 1, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl * weights;
            }
        }
    } else {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl;
            }
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}

void Comp_Pair_BinomNNGaussCop2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                   double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                   double *nuis,  int *type_cop, int *cond)
{
    const int model = 11;
    const double nugget = nuis[0];
    const int npairs_val = npairs[0];

    const int weigthed_val = *weigthed;
    const int type_cop_val = type_cop[0];
    const int cond_val = cond[0];

    if (nugget < 0 || nugget > 1) {
        *res = LOW;
        return;
    }

    double total = 0.0;
    double weights = 1.0;

    if (weigthed_val) {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                weights = GeoCachedPairWeight(i, 1, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl * weights;
            }
        }
    } else {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl;
            }
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}

void Comp_Pair_BinomnegGaussCop2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                    double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                    double *nuis,  int *type_cop, int *cond)
{
    const int model = 16;
    const double nugget = nuis[0];
    const int npairs_val = npairs[0];

    const int weigthed_val = *weigthed;
    const int type_cop_val = type_cop[0];
    const int cond_val = cond[0];

    if (nugget < 0 || nugget > 1) {
        *res = LOW;
        return;
    }

    double total = 0.0;
    double weights = 1.0;

    if (weigthed_val) {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                weights = GeoCachedPairWeight(i, 1, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl * weights;
            }
        }
    } else {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl;
            }
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}

void Comp_Pair_LogisticCop2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                double *nuis,  int *type_cop, int *cond)
{
    const int model = 25;
    const double nugget = nuis[0];
    const int npairs_val = npairs[0];

    const int weigthed_val = *weigthed;
    const int type_cop_val = type_cop[0];
    const int cond_val = cond[0];

    if (nugget >= 1 || nugget < 0.0 || !R_FINITE(nuis[1]) || nuis[1] <= 0.0) {
        *res = LOW;
        return;
    }

    double total = 0.0;
    double weights = 1.0;

    if (weigthed_val) {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                weights = GeoCachedPairWeight(i, 1, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl * weights;
            }
        }
    } else {
        for (int i = 0; i < npairs_val; i++) {
            if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
                const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
                const double bl = biv_cop(clamp_corr(corr), type_cop_val, cond_val,
                                          data1[i], data2[i], mean1[i], mean2[i],
                                          nuis, model, N1[i], N2[i]);
                total += bl;
            }
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}
void Comp_Pair_PoisGammaZIP_st2mem(int *cormod, double *data1, double *data2,
                                   int *N1, int *N2, double *par, int *weigthed,
                                   double *res, double *mean1, double *mean2,
                                   double *nuis, int *type_cop, int *cond)
{
    (void)N1; (void)N2; (void)type_cop; (void)cond;
    const double nugget1 = nuis[0], nugget2 = nuis[1];
    const double mup = nuis[2], shape = nuis[4];
    if (!R_FINITE(nugget1) || !R_FINITE(nugget2) ||
        nugget1 < 0.0 || nugget1 >= 1.0 || nugget2 < 0.0 || nugget2 >= 1.0 ||
        !R_FINITE(mup) || !R_FINITE(shape) || !(shape > 0.0) ||
        !npairs || !lags || !lagt) { *res = LOW; return; }
    GeoKahanSum total; GeoKahanInit(&total); int used = 0;
    for (int i = 0; i < npairs[0]; ++i) {
        int r, t; double mi, mj;
        const int prep = GeoPGPreparePairAt(i, data1[i], data2[i],
                                             mean1[i], mean2[i],
                                             &r, &t, &mi, &mj);
        if (prep == 0) continue;
        if (prep < 0) { *res = LOW; return; }
        const double weight = GeoCachedPairWeight(i, weigthed[0], 1);
        if (!R_FINITE(weight)) { *res = LOW; return; }
        if (weight <= 0.0) continue;
        double rho = CorFct(cormod, lags[i], lagt[i], par, 0, 0);
        if (!R_FINITE(rho) || fabs(rho) > 1.0 + 1e-8) { *res = LOW; return; }
        rho = fmax2(-1.0, fmin2(1.0, rho));
        const double ll = GeoPGZIPLogJoint(rho, r, t, mi, mj, mup,
                                            nugget1, nugget2, shape);
        if (!R_FINITE(ll)) { *res = LOW; return; }
        GeoKahanAdd(&total, weight * ll); ++used;
    }
    *res = (used > 0 && R_FINITE(total.sum)) ? total.sum : LOW;
}
