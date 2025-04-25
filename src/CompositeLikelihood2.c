#include "header.h"

/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/

// Composite marginal (pariwise) log-likelihood for the spatial Gaussian model:
void Comp_Pair_Gauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                        double *par, int *weigthed, double *res, double *mean1, double *mean2,
                        double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Controllo rapido dei parametri
    const double sill = nuis[1];
    const double nugget = nuis[0];
    if(sill < 0 || nugget < 0 || nugget > 1) {
        *res = LOW; 
        return;
    }
    
    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[ 0];
    const double max_dist = maxdist[0];
    const double scale = 1.0 - nugget;
    
    double total = 0.0;
    
    // Loop principale ottimizzato
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Calcolo correlazione
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            
            // Pesi se necessario
            double weights = 1.0;
            if(weighted) {
                weights = CorFunBohman(lag, max_dist);
            }
            
            // Calcolo della verosimiglianza
            total += log_biv_Norm(scale * corr, d1, d2, mean1[i], mean2[i], sill, 0) * weights;
        }
    }
    
    // Assegnazione finale
    *res = R_FINITE(total) ? total : LOW;
}

/******************************************************************************************/
void Comp_Diff_Gauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
  int i=0;
  double vario=0.0,u,v,weights=1.0;

 double nugget=nuis[0];
    double sill=nuis[1];

      if(sill<0 || nugget<0||nugget>=1){*res=LOW; return;}

 for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
  vario=Variogram(cormod,lags[i],0,nuis[0],nuis[1],par);
         u=data1[i];v=data2[i];
            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
    *res+= -0.5*(log(2*M_PI)+log(vario)+
                   R_pow(u-v,2)/(2*vario))*weights;}}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/******************************************************************************************/
/******************************************************************************************/
void Comp_Pair_WrapGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                           double *par, int *weigthed, double *res, double *mean1, double *mean2,
                           double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Controllo rapido dei parametri
    const double sill = nuis[1];
    const double nugget = nuis[0];
    if(sill < 0 || nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }

    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double scale = 1.0 - nugget;
    const double alfa = 2.0;
    double total = 0.0;

    // Loop principale ottimizzato
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Calcolo correlazione
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            const double wrap_gauss = biv_wrapped(alfa, d1, d2, mean1[i], mean2[i], 
                                               nugget, sill, scale * corr);
            double weights = 1.0;
            if(weighted) {
                weights = CorFunBohman(lag, max_dist);
            }
            total += log(wrap_gauss) * weights;
        }
    }

    // Assegnazione finale con controllo
    *res = R_FINITE(total) ? total : LOW;
}
/******************************************************************************************/
void Comp_Pair_SinhGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                            double *par, int *weigthed, double *res, double *mean1, double *mean2,
                            double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Controllo rapido dei parametri
    const double nuis0 = nuis[0];
    const double nuis1 = nuis[1];
    const double nuis2 = nuis[2];
    const double nuis3 = nuis[3];
    
    if(nuis3 < 0 || nuis1 < 0 || nuis0 < 0 || nuis0 >= 1) {
        *res = LOW;
        return;
    }

    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double scale = 1.0 - nuis0;
    double total = 0.0;
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Calcolo correlazione
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            double weights = 1.0;
            if(weighted) {
                weights = CorFunBohman(lag, max_dist);
            }
            const double sinh_val = biv_sinh(scale * corr, d1, d2, mean1[i], mean2[i], 
                                          nuis2, nuis3, nuis1);
            total += weights * log(sinh_val);
        }
    }

    // Assegnazione finale con controllo
    *res = R_FINITE(total) ? total : LOW;
}
/******************************************************************************************/
void Comp_Pair_SkewGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                            double *par, int *weigthed, double *res, double *mean1, double *mean2,
                            double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Controllo rapido dei parametri
    const double sill = nuis[1];
    const double nugget = nuis[0];
    if(nugget < 0 || nugget >= 1 || sill < 0) {
        *res = LOW;
        return;
    }

    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double skew_param = nuis[2];  // Parametro di skewness
    
    double total = 0.0;

    // Loop principale ottimizzato
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Calcolo correlazione
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            
            // Applicazione pesi se necessario
            double weights = 1.0;
            if(weighted) {
                weights = CorFunBohman(lag, max_dist);
            }
            
            // Calcolo biv_skew e accumulo risultato
            const double skew_val = biv_skew(corr, d1, d2, mean1[i], mean2[i], 
                                          sill, skew_param, nugget);
            total += weights * log(skew_val);
        }
    }

    // Assegnazione finale con controllo
    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Pair_Gamma2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                        double *par, int *weigthed, double *res, double *mean1, double *mean2,
                        double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Controllo precoce dei parametri
    const double nugget = nuis[0];
    const double shape_param = nuis[2]; // Parametro shape della gamma
    
    if(nugget < 0 || nugget >= 1 || shape_param < 0) {
        *res = LOW;
        return;
    }

    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double scale = 1.0 - nugget; // (1-nugget) precalcolato
    
    double total = 0.0;
                double weights = 1.0;

    // Loop principale ottimizzato
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Calcolo correlazione
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            
            // Calcolo della bivariata gamma
            const double gamma_val = biv_gamma(scale * corr, d1, d2, 
                                             mean1[i], mean2[i], shape_param);
            
            // Applicazione pesi se necessario

            if(weighted) {
                weights = CorFunBohman(lag, max_dist);
            }
            
            // Accumulo risultato
            total += weights * log(gamma_val);
        }
    }

    // Assegnazione finale con controllo
    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Pair_Weibull2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                          double *par, int *weigthed, double *res, double *mean1, double *mean2,
                          double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Controllo precoce dei parametri
    const double nugget = nuis[0];
    const double shape_param = nuis[2]; // Parametro shape della Weibull
    
    if(nugget < 0 || nugget >= 1 || shape_param < 0) {
        *res = LOW;
        return;
    }

    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double scale = 1.0 - nugget;
    
    double total = 0.0;
    double weights = 1.0;  // Valore di default

    // Loop principale ottimizzato
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Calcolo correlazione
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            
            // Reset del peso a default
            weights = 1.0;
            
            // Calcolo peso se necessario
            if(weighted) {
                weights = CorFunBohman(lag, max_dist);
            }
            
            // Calcolo della Weibull bivariata e accumulo risultato
            const double wbl = biv_Weibull(scale * corr, d1, d2, 
                                         mean1[i], mean2[i], shape_param);
            total += weights * log(wbl);
        }
    }

    // Assegnazione finale con controllo
    *res = R_FINITE(total) ? total : LOW;
}

/*********************************************************/
void Comp_Pair_Kumaraswamy2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                     if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                  bl=biv_Kumara((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);

        *res+= weights*log(bl);
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}


/*********************************************************/
void Comp_Pair_Kumaraswamy22mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                     if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                  bl=biv_Kumara2((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);

        *res+= weights*log(bl);
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Beta2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                     if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                  bl=biv_beta((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);

        *res+= weights*log(bl);
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_LogGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                          double *par, int *weigthed, double *res, double *mean1, double *mean2,
                          double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Controllo precoce dei parametri
    const double sill = nuis[1];
    const double nugget = nuis[0];
    
    if(sill < 0 || nugget < 0 || nugget > 1) {
        *res = LOW;
        return;
    }

    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double scale = 1.0 - nugget;
    
    double total = 0.0;  // Variabile accumulatore
    double weights = 1.0; // Valore di default per i pesi

    // Loop principale ottimizzato
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Calcolo correlazione
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            
            // Reset del peso a default
            weights = 1.0;
            
            // Calcolo peso se necessario
            if(weighted) {
                weights = CorFunBohman(lag, max_dist);
            }
            
            // Calcolo della log-normale bivariata e accumulo risultato
            const double lognorm_val = d2lognorm(d1, d2, sill, nugget, 
                                               mean1[i], mean2[i], scale * corr);
            total += weights * log(lognorm_val);
        }
    }

    // Assegnazione finale con controllo
    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Pair_PoisbinnegGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                 double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                 double *nuis, int *local, int *GPU, int *type_cop, int *cond)
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
    const double max_dist = maxdist[0];
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
            const double p11 = pbnorm22(ai, aj, scale * corr);
            
            // Calcolo probabilità marginali
            const double p1 = pnorm(ai, 0, 1, 1, 0);
            const double p2 = pnorm(aj, 0, 1, 1, 0);
            
            // Conversione a interi
            const int uu = (int)d1;
            const int vv = (int)d2;
            
            // Calcolo pesi
            weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            
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
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                    p11=pbnorm22(ai,aj,(1-nugget)*corr);
                p1=pnorm(ai,0,1,1,0);p2=pnorm(aj,0,1,1,0);
                    u=data1[i];v=data2[i];
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u; vv=(int) v;
                        bl=biv_poisbin(N1[0],uu,vv,p1,p2,p11);

                    *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}


/*********************************************************/
void Comp_Pair_BinomnegGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Controllo precoce dei parametri
    const double nugget = nuis[0];
    if(nugget >= 1 || nugget < 0) {
        *res = LOW;
        return;
    }

    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const int N = N1[0];  // Parametro N della Binomiale Negativa
    const double scale = 1.0 - nugget;
    
    double total = 0.0;  // Variabile accumulatore

    // Loop principale ottimizzato
    for(int i = 0; i < n_pairs; i++) {
        if(!ISNAN(data1[i]) && !ISNAN(data2[i])) {
            // Calcolo correlazione e probabilità
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            const double ai = mean1[i];
            const double aj = mean2[i];
            
            // Calcolo probabilità congiunta e marginali
            const double p11 = pbnorm22(ai, aj, scale * corr);
            const double p1 = pnorm(ai, 0, 1, 1, 0);
            const double p2 = pnorm(aj, 0, 1, 1, 0);
            
            // Conversione a interi e calcolo pesi
            const int uu = (int)data1[i];
            const int vv = (int)data2[i];
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            
            // Calcolo verosimiglianza e accumulo
            const double binomneg_val = biv_binomneg(N, uu, vv, p1, p2, p11);
            total += weights * log(binomneg_val);
        }
    }

    // Assegnazione finale con controllo
    *res = R_FINITE(total) ? total : LOW;
}
/******************************************************/
void Comp_Pair_BinomnegBinary2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
       if( nugget>=1||nugget<0){*res=LOW; return;}
    

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                  corr=CorFct(cormod,lags[i],0,par,0,0);
                p11=pbnorm22(ai,aj,(1-nugget)*corr);
                p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);
                u=data1[i];v=data2[i];
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                uu=(int) u;vv=(int) v;
                bl=log(biv_binegbinary(N1[0],uu,vv,p1,p2,p11));

            *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_BinomnegLogi2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
       if( nugget>=1 || nugget<0){*res=LOW; return;}
    //compute the composite log-likelihood:

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);

                    p11=pblogi22(ai,aj,(1-nugget)*corr);
                    p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));

                    u=data1[i];v=data2[i];
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u;
                         vv=(int) v;
                        bl=biv_binomneg (N1[0],uu,vv,p1,p2,p11);
                    *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_BinomnegGaussZINB2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                    double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                    double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Controllo precoce dei parametri
    const double nugget1 = nuis[0];
    const double nugget2 = nuis[1];
    const double mup = nuis[2];
    
    if(nugget1 < 0 || nugget1 >= 1 || nugget2 < 0 || nugget2 >= 1) {
        *res = LOW;
        return;
    }

    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const int N = N1[0];  // Parametro N della Binomiale Negativa Zero-Inflated
    
    double total = 0.0;  // Variabile accumulatore

    // Loop principale ottimizzato
    for(int i = 0; i < n_pairs; i++) {
        if(!ISNAN(data1[i]) && !ISNAN(data2[i])) {
            // Calcolo correlazione
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            
            // Calcolo pesi
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            
            // Conversione a interi
            const int uu = (int)data1[i];
            const int vv = (int)data2[i];
            
            // Calcolo verosimiglianza ZINB
            const double zinb_val = biv_binomnegZINB(
                N, corr, uu, vv, 
                mean1[i], mean2[i], 
                nugget1, nugget2, mup
            );
            
            // Accumulo risultato
            total += weights * log(zinb_val);
        }
    }

    // Assegnazione finale con controllo
    *res = R_FINITE(total) ? total : LOW;
}
/******************************************************************************************/
void Comp_Pair_BinomGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                             double *par, int *weigthed, double *res, double *mean1, double *mean2,
                             double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Controllo precoce dei parametri
    const double nugget = nuis[0];
    if(nugget >= 1 || nugget < 0) {
        *res = LOW;
        return;
    }

    // Variabili precalcolate
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const int N = N1[0];  // Parametro N della Binomiale
    const double scale = 1.0 - nugget;
    
    double total = 0.0;  // Variabile accumulatore

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
            
            // Calcolo probabilità congiunta e marginali
            const double p11 = pbnorm22(ai, aj, scale * corr);
            const double p1 = pnorm(ai, 0, 1, 1, 0);
            const double p2 = pnorm(aj, 0, 1, 1, 0);
            
            // Conversione a interi e calcolo pesi
            const int uu = (int)d1;
            const int vv = (int)d2;
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            
            // Calcolo verosimiglianza e accumulo
            const double binom_val = biv_binom(N, uu, vv, p1, p2, p11);
            total += weights * log(binom_val);
        }
    }

    // Assegnazione finale con controllo
    *res = R_FINITE(total) ? total : LOW;
}
/******************************************************************************************/
void Comp_Pair_BinomLogi2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                p11=pblogi22(ai,aj,(1-nugget)*corr);
                 p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u; vv=(int) v;
                        bl=biv_binom (N1[0],uu,vv,p1,p2,p11);
                    *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Pair_BinomNNGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                 corr=CorFct(cormod,lags[i],0,par,0,0);
                 p11=pbnorm22(ai,aj,(1-nugget)*corr);
                 p1=pnorm(ai,0,1,1,0);
                 p2=pnorm(aj,0,1,1,0);
                 u=data1[i];v=data2[i];
                 if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                 uu=(int) u; vv=(int) v;
                 bl=biv_binom222(N1[i],N2[i],uu,vv,p1,p2,p11);
                 *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}



void Comp_Pair_BinomNNGauss_misp2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0, N=2,n1,n2;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success

    double **M;
    M= (double **) R_Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) R_Calloc(N,double);}
    double *dat;
    dat=(double *) R_Calloc(N,double);

    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);
                 p11=pbnorm22(ai,aj,(1-nugget)*corr);
                 p1=pnorm(ai,0,1,1,0);
                 p2=pnorm(aj,0,1,1,0);
                 u=data1[i];v=data2[i];
                 n1=N1[i];n2=N2[i];
                 if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                 M[0][0]=n1*p1*(1-p1);   M[1][1]=n2*p2*(1-p2);  // var1 var2
                 M[0][1]= fmin_int(n1,n2)*(p11-p1*p2) ;       // covariance
                 M[1][0]= M[0][1];
                 dat[0]=u-n1*p1;dat[1]=v-n2*p2; 
                 //Rprintf("%d %f %f %f \n",fmin_int(n1,n2),p1,p2,p11 );
                   //#####
                 bl=dNnorm(N,M,dat);
                 *res+= log(bl)*weights;       
                }}
           for(i=0;i<N;i++)  {R_Free(M[i]);}
    R_Free(M);
    if(!R_FINITE(*res))*res = LOW;
    return;
}




void Comp_Pair_BinomNNLogi2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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

                 p11=pblogi22(ai,aj,(1-nugget)*corr);
                 p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                 u=data1[i];v=data2[i];
                 n1=N1[i];n2=N2[i];
                 if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                 uu=(int) u; vv=(int) v;
                 bl=biv_binom222(n1,n2,uu,vv,p1,p2,p11);
                 *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_LogLogistic2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                              double *par, int *weigthed, double *res, double *mean1, double *mean2,
                              double *nuis, int *local, int *GPU, int *type_cop, int *cond)
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
    const double max_dist = maxdist[0];
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
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            
            // Calcolo della Log-Logistica bivariata e accumulo risultato
            const double loglogistic_val = biv_LogLogistic(
                scale * corr, d1, d2, 
                mean1[i], mean2[i], 
                shape_param
            );
            total += weights * log(loglogistic_val);
        }
    }

    // Assegnazione finale con controllo
    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Pair_Logistic2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,nugget=0.0;
        nugget=nuis[0];
    if(nugget>=1||nugget<0.0 ) {*res=LOW;  return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=(data1[i]);zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    bl= biv_Logistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[1]);
                    *res+= weights*log(bl);
                }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Pois2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                       double *par, int *weigthed, double *res, double *mean1, double *mean2,
                       double *nuis, int *local, int *GPU, int *type_cop, int *cond)
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
    const double max_dist = maxdist[0];
    const double scale = 1.0 - nugget;
    
    double total = 0.0;  // Accumulator variable

    // Optimized main loop
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
            
            // Compute weights if needed
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            
            // Compute bivariate Poisson and accumulate result
            const double poisson_val = biv_Poisson(scale * corr, uu, ww, mui, muj);
            total += log(poisson_val) * weights;
        }
    }

    // Final assignment with check
    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Pair_PoisGamma2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                            double *par, int *weigthed, double *res, double *mean1, double *mean2,
                            double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Early parameter validation
    const double nugget = nuis[0];
    const double gamma_param = nuis[2]; // Gamma distribution parameter
    if(nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }

    // Precompute frequently used values
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double scale = 1.0 - nugget;
    
    double total = 0.0;  // Accumulator variable

    // Optimized main loop
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
            
            // Compute weights if needed
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            
            // Compute Poisson-Gamma and accumulate result
            const double poisgamma_val = biv_PoissonGamma(
                scale * corr, uu, ww, mui, muj, gamma_param);
            total += log(poisgamma_val) * weights;
        }
    }

    // Final assignment with check
    *res = R_FINITE(total) ? total : LOW;
}

void Comp_Pair_PoisGammaZIP2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{

    int i=0, uu,vv;
    double weights=1.0,corr,mui,muj,bl,u,v;
   double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2]; double shape=nuis[3];


      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],0,par,0,0);
                        u=data1[i];v=data2[i];
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u;vv=(int) v;
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                      bl=log(biv_PoissonGammaZIP(corr,uu,vv,mui, muj,mup,nugget1,nugget2,shape));
                      *res+= bl*weights;
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}


void Comp_Pair_Gauss_misp_PoisGamma2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                      double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                      double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Early parameter validation
    const double nugget = nuis[0];
    const double gamma_param = nuis[2]; // Gamma distribution parameter
    if(nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }

    // Constants and precomputations
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const int N = 2; // Fixed matrix dimension
    const double scale = 1.0 - nugget;
    
    double total = 0.0;  // Accumulator variable

    // Allocate memory for covariance matrix and data vector
    double **M = (double **) R_Calloc(N, double *);
    double *dat = (double *) R_Calloc(N, double);
    for(int i = 0; i < N; i++) {
        M[i] = (double *) R_Calloc(N, double);
    }

    // Main processing loop
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Compute means and variances
            const double mui = exp(mean1[i]);
            const double muj = exp(mean2[i]);
            const double vvi = mui * (1.0 + mui/gamma_param);
            const double vvj = muj * (1.0 + muj/gamma_param);
            
            // Compute correlations
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0) * scale;
            const double corr1 = corr_pois_gen(corr, mui, muj, gamma_param);
            
            // Compute weights if needed
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            
            // Build covariance matrix
            M[0][0] = vvi;
            M[1][1] = vvj;
            M[0][1] = M[1][0] = sqrt(vvi * vvj) * corr1;
            
            // Center data
            dat[0] = d1 - mui;
            dat[1] = d2 - muj;
            
            // Compute density and accumulate result
            const double density = dNnorm(N, M, dat);
            total += log(density) * weights;
        }
    }

    // Clean up allocated memory
    for(int i = 0; i < N; i++) {
        R_Free(M[i]);
    }
    R_Free(M);
    R_Free(dat);

    // Final assignment with check
    *res = R_FINITE(total) ? total : LOW;
}


/*********************************************************/
void Comp_Pair_PoisZIP2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0, uu,ww;
    double weights=1.0,corr,mui,muj,bl;
   double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];


      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
  // Rprintf("%d   \n",npairs[0]);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],0,par,0,0);
                    // if(fabs(corr)>1|| !R_FINITE(corr)) {*res=LOW; return;}
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                      uu=(int) data1[i];  ww=(int) data2[i];
                      //Rprintf("%f %f \n",mui,muj);
                      bl=biv_PoissonZIP(corr,uu,ww,mui, muj,mup,nugget1,nugget2);
                      *res+= log(bl)*weights;
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Gauss_misp_Pois2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                 double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                 double *nuis, int *local, int *GPU, int *type_cop, int *cond)
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
    const double max_dist = maxdist[0];
    const int N = 2; // Fixed matrix dimension
    const double scale = 1.0 - nugget;
    
    double total = 0.0;  // Accumulator variable

    // Allocate memory for covariance matrix and data vector
    double **M = (double **) R_Calloc(N, double *);
    double *dat = (double *) R_Calloc(N, double);
    for(int i = 0; i < N; i++) {
        M[i] = (double *) R_Calloc(N, double);
    }

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
            const double corr1 = corr_pois(corr, mui, muj);
            
            // Compute weights if needed
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            
            // Build covariance matrix efficiently
            const double sqrt_mui_muj = sqrt(mui * muj);
            M[0][0] = mui;
            M[1][1] = muj;
            M[0][1] = M[1][0] = sqrt_mui_muj * corr1;
            
            // Center data
            dat[0] = d1 - mui;
            dat[1] = d2 - muj;
            
            // Compute density and accumulate result
            total += log(dNnorm(N, M, dat)) * weights;
        }
    }

    // Clean up allocated memory
    for(int i = 0; i < N; i++) {
        R_Free(M[i]);
    }
    R_Free(M);
    R_Free(dat);

    // Final assignment with check
    *res = R_FINITE(total) ? total : LOW;
}

void Comp_Pair_Gauss_misp_PoisZIP2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0;
    double weights=1.0,corr,mui,muj,bl;
    double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];

      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}

  // Rprintf("%d   \n",npairs[0]);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],0,par,0,0);
                      if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                      bl=biv_Mis_PoissonZIP(corr,data1[i],data2[i],mui, muj,mup,nugget1,nugget2);
                  //    Rprintf("%f %f\n",bl,mup);
                      *res+= log(bl)*weights;
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_Gauss_misp_SkewT2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0;
    double weights=1.0,sill,nugget,skew,corr,corr2,df,bl;


    df=1/nuis[0];
    nugget=nuis[1];
    sill=nuis[2];
    skew=nuis[3];

    if(df<2||fabs(skew)>1||sill<0||nugget<0||nugget>=1){*res=LOW; return;}
    //auxuliary variables
    double D1=(df-1)/2;
    double D2=df/2;
    //double delta=skew/sqrt(1-skew*skew);
    double MM=sqrt(df)*gammafn(D1)*skew/(sqrt(M_PI)*gammafn(D2));
    double FF=(df/(df-2)-MM*MM);

     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                     corr=CorFct(cormod,lags[i],0,par,0,0)*(1-nugget);
                     corr2=corr_skewt(corr,df,skew);
                     if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
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
                               double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Early parameter validation
    const double sill = nuis[2];
    const double nugget = nuis[1];
    const double df_param = nuis[0];
    
    if(sill < 0 || nugget < 0 || nugget >= 1 || df_param <= 0 || df_param > 0.5) {
        *res = LOW;
        return;
    }

    // Precompute frequently used values
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double df = 1.0 / df_param;
    const double sill_scale = sill * df / (df - 2);
    const double scale = 1.0 - nugget;
    
    double total = 0.0;  // Accumulator variable

    // Precompute log terms for correlation calculation
    const double log_df_minus_2 = log(df - 2);
    const double two_lgamma_df_minus_1_half = 2 * lgammafn(0.5 * (df - 1));
    const double two_lgamma_df_half = 2 * lgammafn(0.5 * df);
    const double log_constant = log_df_minus_2 + two_lgamma_df_minus_1_half - (log(2) + two_lgamma_df_half);

    // Main processing loop
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Compute correlation
            const double lag = lags[i];
            double corr = CorFct(cormod, lag, 0, par, 0, 0);
            
            // Transform correlation for T distribution
            if(fabs(corr) > 0) {
                const double corr_sq = corr * corr;
                corr = exp(log_constant + 
                          log(hypergeo(0.5, 0.5, 0.5 * df, corr_sq)) + 
                          log(corr * scale));
            }
            
            // Compute weights if needed
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            
            // Compute bivariate normal log-likelihood and accumulate
            total += log_biv_Norm(corr, d1, d2, mean1[i], mean2[i], sill_scale, 0) * weights;
        }
    }

    // Final assignment with check
    *res = R_FINITE(total) ? total : LOW;
}


/*********************************************************/
 void Comp_Pair_Tukeyhh2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                          double *par, int *weigthed, double *res, double *mean1, double *mean2,
                          double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Early parameter validation
    const double sill = nuis[1];
    const double nugget = nuis[0];
    const double h1 = nuis[3];
    const double h2 = nuis[2];
    
    if(sill < 0 || h1 < 0 || h1 > 0.5 || h2 < 0 || h2 > 0.5 || nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }

    // Precompute frequently used values
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
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
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            
            // Compute bivariate Tukey h-h likelihood
            const double tukey_val = biv_tukey_hh(
                scale * corr, d1, d2, 
                mean1[i], mean2[i], 
                sill, h1, h2
            );
            
            // Accumulate result
            total += weights * log(tukey_val);
        }
    }

    // Final assignment with check
    *res = R_FINITE(total) ? total : LOW;
}




/*********************************************************/
void Comp_Pair_Tukeyh2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                         double *par, int *weigthed, double *res, double *mean1, double *mean2,
                         double *nuis, int *local, int *GPU, int *type_cop, int *cond)
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
    const double max_dist = maxdist[0];
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
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            
            // Compute bivariate Tukey h likelihood and accumulate
            const double tukey_val = biv_tukey_h(
                corr, d1, d2, 
                mean1[i], mean2[i], 
                tail, sill
            );
            total += weights * log(tukey_val);
        }
    }

    // Final assignment with check
    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Pair_Gauss_misp_Tukeygh2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                    corr2=corr_tukeygh(corr,eta,tail);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
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
                    double *nuis, int *local, int *GPU, int *type_cop, int *cond) 
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
    const double max_dist = maxdist[0];
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
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            const double t_val = biv_T(corr, std_d1, std_d2, df, nugget) * inv_sill;
            total += weights * log(t_val);
        }
    }

    // Assegnazione finale con controllo
    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Pair_TWOPIECETukeyh2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                      p11=pbnorm22(qq,qq,corr);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    bl=biv_two_pieceTukeyh((1-nugget)*corr,zi,zj,sill,eta,tail,p11,mean1[i],mean2[i]);
                    *res+= weights*log(bl);
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_TWOPIECET2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                      p11=pbnorm22(qq,qq,corr);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    /********************************************************/
    bl=biv_two_pieceT(corr,zi,zj,sill,df,eta,p11,mean1[i],mean2[i],nugget);
  //Rprintf("%f- %f- %f %f %f %f %f %f %f %f %f  \n",lags[i],corr,eta,sill,zi,zj,df,eta,p11,mean1[i],mean2[i]);
                    /********************************************************/
                           *res+= weights*log(bl);
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_TWOPIECEGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                double *nuis, int *local, int *GPU, int *type_cop, int *cond)
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
    const double max_dist = maxdist[0];
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
            const double p11 = pbnorm22(qq, qq, corr);
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            const double density = biv_two_pieceGaussian(
                scale * corr, d1, d2, 
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
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                        p11=pbnorm22(qq,qq,corr);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    /********************************************************/
                   bl=biv_two_piece_bimodal((1-nugget)*corr,zi,zj,sill,df,delta,eta,p11,mean1[i],mean2[i]);
          // Rprintf("%f %f  --%f  %f %f %f %f  -%f %f \n",bl,lags[i],df,delta,eta,sill,corr,par[0],par[1]);
                    /********************************************************/
                           *res+= weights*log(bl);
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
                           double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Estrai valori una sola volta
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double max_time = maxtime[0];
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
            const double adjusted_corr = one_minus_nugget * corr;
            
            // Calcola la log verosimiglianza bivariata
            const double bl = log_biv_Norm(adjusted_corr, d1, d2, mean1_ptr[i], mean2_ptr[i], sill, 0);
            
            // Calcola i pesi se necessario
            double weights = 1.0;
            if(weighted) {
                weights = CorFunBohman(current_lag, max_dist) * CorFunBohman(current_lagt, max_time);
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
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0;
    double  u=0.0, w=0.0,weights=1.0;
    double bl=0.0,corr=0.0;
    double alfa=2.0;  double nugget=nuis[0];  double sill =nuis[1];
   if(sill<0 || nugget<0||nugget>=1){*res=LOW; return;}
   for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                u=data1[i];w=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                    bl=biv_wrapped(alfa,u,w,mean1[i],mean2[i],nugget,sill,(1-nugget)*corr);
                                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                             *res+= weights*log(bl);
                                }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_T_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                       double *par, int *weigthed, double *res, double *mean1, double *mean2,
                       double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Estrai valori una sola volta
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double max_time = maxtime[0];
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
                weights = CorFunBohman(lags_ptr[i], max_dist) * 
                         CorFunBohman(lagt_ptr[i], max_time);
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
                                  double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Estrai valori una sola volta
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double max_time = maxtime[0];
    const double sill = nuis[2];
    const double nugget = nuis[1];
    const double nu = nuis[0];  // Parametro di degrees of freedom
    
    // Controllo preliminare dei parametri
    if(sill < 0 || nugget < 0 || nugget >= 1 || nu < 0 || nu > 0.5) {
        *res = LOW;
        return;
    }

    double total = 0.0;
    const int weighted = *weigthed;
    const double df = 1.0 / nu;  // df = 1/nu
    const double sill_adjusted = sill * df / (df - 2.0);
    const double one_minus_nugget = 1.0 - nugget;

    // Ottimizzazione: puntatori per accesso continuo alla memoria
    const double *lags_ptr = lags;
    const double *lagt_ptr = lagt;
    const double *data1_ptr = data1;
    const double *data2_ptr = data2;
    const double *mean1_ptr = mean1;
    const double *mean2_ptr = mean2;

    // Pre-calcola costanti per la trasformazione della correlazione
    const double log_df_minus_2 = log(df - 2.0);
    const double two_lgamma_half_df_minus_1 = 2.0 * lgammafn(0.5 * (df - 1.0));
    const double log_2_plus_two_lgamma_half_df = log(2.0) + 2.0 * lgammafn(0.5 * df);

    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1_ptr[i];
        const double d2 = data2_ptr[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Calcola la correlazione base
            double corr = CorFct(cormod, lags_ptr[i], lagt_ptr[i], par, 0, 0);
            
            // Trasforma la correlazione (parte più costosa)
            const double corr_sq = corr * corr;
            const double log_hypergeo = log(hypergeo(0.5, 0.5, 0.5 * df, corr_sq));
            const double log_corr_part = log(corr * one_minus_nugget);
            
            corr = exp(log_df_minus_2 + two_lgamma_half_df_minus_1 - 
                      log_2_plus_two_lgamma_half_df + log_hypergeo + log_corr_part);
            
            // Calcola i pesi se necessario
            double weights = 1.0;
            if(weighted) {
                weights = CorFunBohman(lags_ptr[i], max_dist) * 
                         CorFunBohman(lagt_ptr[i], max_time);
            }
            
            // Calcola la log-verosimiglianza
            const double bl = log_biv_Norm(corr, d1, d2, mean1_ptr[i], mean2_ptr[i], sill_adjusted, 0);
            total += bl * weights;
        }
    }

    *res = R_FINITE(total) ? total : LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Gauss_misp_Pois_stmem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                   double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                   double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Estrai valori una sola volta
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double max_time = maxtime[0];
    const double nugget = nuis[0];
    
    // Controllo preliminare dei parametri
    if(nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }

    double total = 0.0;
    const int weighted = *weigthed;
    const double one_minus_nugget = 1.0 - nugget;
    const int N = 2;  // Dimensione fissa della matrice

    // Allocazione memoria una sola volta
    double **M = (double **) R_Calloc(N, double *);
    double *dat = (double *) R_Calloc(N, double);
    for(int i = 0; i < N; i++) {
        M[i] = (double *) R_Calloc(N, double);
    }

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
            // Calcola la correlazione
            const double corr = one_minus_nugget * CorFct(cormod, lags_ptr[i], lagt_ptr[i], par, 0, 0);
            
            // Calcola le medie
            const double mui = exp(mean1_ptr[i]);
            const double muj = exp(mean2_ptr[i]);
            const double sqrt_mui_muj = sqrt(mui * muj);
            
            // Calcola la correlazione di Poisson
            const double corr2 = corr_pois(corr, mui, muj);
            
            // Riempimento matrice M
            M[0][0] = mui;
            M[1][1] = muj;
            M[0][1] = sqrt_mui_muj * corr2;
            M[1][0] = M[0][1];
            
            // Calcola le differenze
            dat[0] = d1 - mui;
            dat[1] = d2 - muj;
            
            // Calcola i pesi se necessario
            double weights = 1.0;
            if(weighted) {
                weights = CorFunBohman(lags_ptr[i], max_dist) * 
                         CorFunBohman(lagt_ptr[i], max_time);
            }
            
            // Calcola la densità normale multivariata
            const double bl = dNnorm(N, M, dat);
            total += log(bl) * weights;
        }
    }

    // Libera memoria allocata
    for(int i = 0; i < N; i++) {
        R_Free(M[i]);
    }
    R_Free(M);
    R_Free(dat);

    *res = R_FINITE(total) ? total : LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Pois_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                          double *par, int *weigthed, double *res, double *mean1, double *mean2,
                          double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Pre-carica le variabili globali in variabili locali
    const int npairs_val = npairs[0];
    const int weighted_flag = *weigthed;
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];
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
            const double bl = biv_Poisson((1.0 - nugget) * corr, uu, ww, mui, muj);
            if(weighted_flag) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                          CorFunBohman(lagt[i], maxtime_val);
            }
            // Accumula il risultato
            total += log(bl) * weights;
        }
    }
    
    // Assegna il risultato finale
    *res = R_FINITE(total) ? total : LOW;
}

/******************************************************************************************/
void Comp_Pair_PoisGamma_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                               double *par, int *weigthed, double *res, double *mean1, double *mean2,
                               double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Pre-carica variabili globali e parametri
    const int npairs_val = npairs[0];
    const double nugget = nuis[0];
    const double nuis2 = nuis[2];  // Parametro aggiuntivo per PoissonGamma
    //const int weighted_flag = *weigthed;
    
    // Controllo iniziale sul nugget
    if(nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }
    
    double total = 0.0;
    const double weights = 1.0;  // Costante poiché il peso è commentato nel codice originale
    
    // Loop principale ottimizzato
    for(int i = 0; i < npairs_val; i++) {
        const double u = data1[i];
        const double w = data2[i];
        
        if(!ISNAN(u) && !ISNAN(w)) {
            // Calcola la correlazione
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0);
            
            // Calcola le medie
            const double mui = exp(mean1[i]);
            const double muj = exp(mean2[i]);
            
            // Converti i valori in interi
            const int uu = (int)u;
            const int ww = (int)w;
            
            // Calcola la probabilità bivariata Poisson-Gamma
            const double bl = biv_PoissonGamma((1.0 - nugget) * corr, uu, ww, mui, muj, nuis2);
            
            // Nota: la parte dei pesi è commentata nel codice originale
            // if(weighted_flag) {
            //     const double weights = CorFunBohman(lags[i], maxdist[0]) * 
            //                          CorFunBohman(lagt[i], maxtime[0]);
            // }
            
            // Accumula il risultato
            total += log(bl) * weights;
        }
    }
    
    // Assegna il risultato finale
    *res = R_FINITE(total) ? total : LOW;
}



/******************************************************************************************/
void Comp_Pair_Gauss_misp_PoisGamma_st2mem(int *cormod, double *data1, double *data2, 
                                          int *N1, int *N2, double *par, int *weigthed, 
                                          double *res, double *mean1, double *mean2,
                                          double *nuis, int *local, int *GPU, 
                                          int *type_cop, int *cond)
{
    // Pre-carica variabili globali e parametri
    const int npairs_val = npairs[0];
    const double nugget = nuis[0];
    const double nuis2 = nuis[2];  // Parametro gamma per la varianza
    //const int weighted_flag = *weigthed;
    
    // Controllo iniziale sul nugget
    if(nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }

    // Allocazione matrici e vettori (ottimizzata)
    const int N = 2;
    double **M = (double **) R_Calloc(N, double *);
    for(int i = 0; i < N; i++) {
        M[i] = (double *) R_Calloc(N, double);
    }
    double *dat = (double *) R_Calloc(N, double);
    
    double total = 0.0;
    const double weights = 1.0;  // Costante poiché il peso è commentato
    
    // Variabili riutilizzabili
    double corr, corr2, mui, muj, vvi, vvj, u, w;
    
    // Loop principale ottimizzato
    for(int i = 0; i < npairs_val; i++) {
        u = data1[i];
        w = data2[i];
        
        if(!ISNAN(u) && !ISNAN(w)) {
            // Calcola correlazione e medie
            corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0) * (1.0 - nugget);
            mui = exp(mean1[i]);
            muj = exp(mean2[i]);
            
            // Calcola varianze
            vvi = mui * (1.0 + mui / nuis2);
            vvj = muj * (1.0 + muj / nuis2);
            
            // Calcola correlazione per Poisson-Gamma
            corr2 = corr_pois_gen(corr, mui, muj, nuis2);
            
            // Riempimento matrice di covarianza simmetrica
            M[0][0] = vvi;
            M[1][1] = vvj;
            M[0][1] = M[1][0] = sqrt(vvi * vvj) * corr2;  // Assegnazione simultanea
            
            // Calcolo scostamenti
            dat[0] = u - mui;
            dat[1] = w - muj;
            
            // Calcolo densità normale multivariata
            const double bl = dNnorm(N, M, dat);
            total += log(bl) * weights;
        }
    }
    
    // Assegna il risultato finale
    *res = R_FINITE(total) ? total : LOW;
    
    // Libera memoria allocata
    for(int i = 0; i < N; i++) {
        R_Free(M[i]);
    }
    R_Free(M);
    R_Free(dat);
}

/******************************************************************************************/
void Comp_Pair_PoisZIP_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)

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
                      bl=biv_PoissonZIP(corr,uu,ww,mui, muj,mup,nugget1,nugget2);
                //   Rprintf("%d %d--%f %f %f  %f \n",uu,ww,lags[i],lagt[i],corr,bl);
               // if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                       *res+= log(bl)*weights;

                                    }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Gauss_misp_Pois_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0,N=2;
     double weights=1.0,corr,corr2,mui,muj,bl,u=0.0, w=0.0;
   double nugget=nuis[0];
     if(nugget<0||nugget>=1){*res=LOW; return;}
// ###
    double **M;
        M= (double **) R_Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) R_Calloc(N,double);}
    double *dat;
    dat=(double *) R_Calloc(N,double);
// ###
    for(i=0;i<npairs[0];i++){
       u=data1[i];w=data2[i];
if(!ISNAN(u)&&!ISNAN(w) ){
                            corr=CorFct(cormod,lags[i],lagt[i],par,0,0)*(1-nugget);
                            mui=exp(mean1[i]);muj=exp(mean2[i]);
                            corr2=corr_pois(corr,mui, muj);

                            M[0][0]=mui; M[1][1]=muj;M[0][1]=sqrt(mui*muj)*corr2;M[1][0]= M[0][1];
                           dat[0]=u-mui;dat[1]=w-muj;
                          //if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                              bl=dNnorm(N,M,dat);
                              *res+= log(bl)*weights;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Gauss_misp_PoisZIP_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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

               bl=biv_Mis_PoissonZIP(corr,data1[i],data2[i],mui, muj,mup,nugget1,nugget2);
                //   Rprintf("%d %d--%f %f %f  %f \n",uu,ww,lags[i],lagt[i],corr,bl);
               // if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                       *res+= log(bl)*weights;

                                    }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Tukeyh_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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

                                   if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);

 bl=biv_tukey_h((1-nugget)*corr,zi,zj,mean1[i],mean2[i],tail,sill);
                             *res+= weights*log(bl);

                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Tukeyhh_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{


    int i=0;
    double corr,zi,zj,weights=1.0,bl;
double sill=nuis[1];
    double nugget=nuis[0];
    double h1=nuis[3];
    double h2=nuis[2];
      if( sill<0||h1<0||h1>0.5||h2<0||h2>0.5||nugget<0||nugget>=1){*res=LOW; return;}

for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
bl=biv_tukey_hh((1-nugget)*corr,zi,zj,mean1[i],mean2[i],sill,h1,h2);
                             *res+= weights*log(bl);
                         }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void  Comp_Pair_TWOPIECEGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                                   if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                        p11=pbnorm22(qq,qq,corr);
  bl=biv_two_pieceGaussian((1-nugget)*corr,zi,zj,sill,eta,p11,mean1[i],mean2[i]);

                           *res+= weights*log(bl);
                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}

/******************************************************************************************/
void  Comp_Pair_TWOPIECETukeyh_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                                  if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                      // p11=pbnorm(cormod,lags[i],lagt[i],qq,qq,nugget,1,par,0);
                        p11=pbnorm22(qq,qq,corr);
                         bl=biv_two_pieceTukeyh((1-nugget)*corr,zi,zj,sill,eta,tail,p11,mean1[i],mean2[i]);

                           *res+= weights*log(bl);

                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}

/******************************************************************************************/
void  Comp_Pair_TWOPIECET_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                        p11=pbnorm22(qq,qq,corr);
                         bl=biv_two_pieceT(corr,zi,zj,sill,df,eta,p11,mean1[i],mean2[i],nugget);

                           *res+= weights*log(bl);

                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}

/******************************************************************************************/
void Comp_Pair_PoisbinGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                //psj=pbnorm(cormod,lags[i],lagt[i],mean1[i],mean2[i],nuis[0],1,par,0);
                p11=pbnorm22(mean1[i],mean2[i],(1-nugget)*corr);
                p1=pnorm(mean1[i],0,1,1,0);
                p2=pnorm(mean2[i],0,1,1,0);
                u=data1[i];w=data2[i];
                                     uu=(int) u;
                                     ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                    dens=biv_poisbin (N1[0],uu,ww,p1,p2,p11);
                                 *res+=log(dens)*weights;
                               }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_PoisbinnegGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                 p11=pbnorm22(mean1[i],mean2[i],(1-nugget)*corr);
                 p1=pnorm((mean1[i]),0,1,1,0);
                 p2=pnorm((mean2[i]),0,1,1,0);
                                u=data1[i];w=data2[i];

                                     uu=(int) u; ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                    dens=biv_poisbinneg (N1[0],uu,ww,p1,p2,p11);
                                 *res+=log(dens)*weights;}}

    if(!R_FINITE(*res))*res = LOW;
    return;
}


void Comp_Pair_BinomGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                            psj=pbnorm22(a,b,(1-nugget)*corr);
                            p1=pnorm(a,0,1,1,0);
                            p2=pnorm(b,0,1,1,0);
                            uu=(int) u;  ww=(int) w;
                          if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                    bl=biv_binom (N1[0],uu,ww,p1,p2,psj);
                                 *res+=log(bl)*weights;
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}
void Comp_Pair_BinomLogi_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                            psj=pblogi22(a,b,(1-nugget)*corr);
                            p1=1/(1+exp(-a));p2=1/(1+exp(-b));
                            uu=(int) u;  ww=(int) w;
                          if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                    bl=biv_binom (N1[0],uu,ww,p1,p2,psj);
                                 *res+=log(bl)*weights;
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}
void Comp_Pair_BinomNNGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                            psj=pbnorm22(a,b,(1-nugget)*corr);
                            p1=pnorm(a,0,1,1,0);
                            p2=pnorm(b,0,1,1,0);
                            uu=(int) u;  ww=(int) w;
                            n1=N1[i];n2=N2[i];
                          if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                    bl=biv_binom222(n1,n2,uu,ww,p1,p2,psj);
                                 *res+=log(bl)*weights;
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}

void Comp_Pair_BinomNNLogi_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                            psj=pblogi22(a,b,(1-nugget)*corr);
                            p1=1/(1+exp(-a));p2=1/(1+exp(-b));
                            uu=(int) u;  ww=(int) w;
                            n1=N1[i];n2=N2[i];
                          if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                    bl=biv_binom222(n1,n2,uu,ww,p1,p2,psj);
                                 *res+=log(bl)*weights;
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}

void Comp_Pair_BinomnegGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0, uu=0,ww=0;
    double bl=0.0,weights=1.0,u=0.0,w=0.0,a=0.0,b=0.0,corr=0.0,p1=0.0,p2=0.0,psj=0.0;

      double nugget=nuis[0];
      if(nugget<0||nugget>=1){*res=LOW; return;}

           for(i=0;i<npairs[0];i++){
                   u=data1[i];w=data2[i];
             if(!ISNAN(u)&&!ISNAN(w) ){

                             corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                            a=mean1[i];b=mean2[i];
                            psj=pbnorm22(a,b,(1-nugget)*corr);
                            p1=pnorm(a,0,1,1,0);
                            p2=pnorm(b,0,1,1,0);
                            uu=(int) u;  ww=(int) w;
                          if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
             bl=biv_binomneg (N1[0],uu,ww,p1,p2,psj);
                                 *res+=log(bl)*weights;
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}


void Comp_Pair_BinomnegLogi_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0, uu=0,ww=0;
    double bl=0.0,weights=1.0,u=0.0,w=0.0,a=0.0,b=0.0,corr=0.0,p1=0.0,p2=0.0,psj=0.0;

      double nugget=nuis[0];
      if(nugget<0||nugget>=1){*res=LOW; return;}

           for(i=0;i<npairs[0];i++){
                   u=data1[i];w=data2[i];
             if(!ISNAN(u)&&!ISNAN(w) ){

                             corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                            a=mean1[i];b=mean2[i];

                            psj=pblogi22(a,b,(1-nugget)*corr);
                            p1=1/(1+exp(-a));p2=1/(1+exp(-b));

                            uu=(int) u;  ww=(int) w;
                          if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
             bl=biv_binomneg (N1[0],uu,ww,p1,p2,psj);
                                 *res+=log(bl)*weights;
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}

void Comp_Pair_BinomnegGaussZINB_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0, uu=0,ww=0;
    double bl=0.0,weights=1.0,u=0.0,w=0.0,a=0.0,b=0.0,corr=0.0;

     double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];
      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
           for(i=0;i<npairs[0];i++){
                   u=data1[i];w=data2[i];
             if(!ISNAN(u)&&!ISNAN(w) ){
                            corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                            a=mean1[i];b=mean2[i];
                          uu=(int) u;  ww=(int) w;
                          if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                              bl=biv_binomnegZINB(N1[0],corr,uu,ww,a,b,nugget1,nugget2,mup);

                                 *res+=log(bl)*weights;
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}




/*********************************************************************************/

// Composite marginal (difference) log-likelihood for the spatial-temporal Gaussian model:
void Comp_Diff_Gauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=00;
    double u,w,vario=0.0,weights=1.0;
    double nugget=nuis[0];
    double sill=nuis[1];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}

      for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                vario=Variogram(cormod,lags[i],lagt[i],nuis[0],nuis[1],par);
                                      u=data1[i];w=data2[i];
                                              if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                        *res+= (-0.5*(log(2*M_PI)+log(vario)+
                                                     R_pow(u-w,2)/(2*vario)))*weights;}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/******************************************************************************************/
void Comp_Pair_SkewGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                    bl=biv_skew(corr,zi,zj,mean1[i],mean2[i],nuis[1],nuis[2],nugget);
                             *res+= weights*log(bl);
                         }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_SinhGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{

    int i=0;
    double bl,corr,zi,zj,weights=1.0;
       if(nuis[3]<0||nuis[1]<0||nuis[0]<0||nuis[0]>=1) {*res=LOW;  return;}

       for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){

                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0)*(1-nuis[0]);
                                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                    bl=biv_sinh(corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],nuis[1]);

                             *res+= weights*log(bl);
                         }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Gamma_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                  bl=biv_gamma((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]);
                                    *res+= weights*log(bl);
                         }}}}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Kumaraswamy_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                   bl= biv_Kumara((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);
                             *res+= weights*log(bl);
                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}

/******************************************************************************************/
void Comp_Pair_Kumaraswamy2_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                   bl= biv_Kumara2((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);
                             *res+= weights*log(bl);
                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Beta_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                   bl= biv_beta((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);
                             *res+= weights*log(bl);
                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Weibull_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{

    int i=0;
    double corr,zi,zj,weights=1.0,bl=0.0;
  double nugget=nuis[0];
     if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
    //Rprintf("%d\n",npairs[0]);
       for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                   // if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                      bl=biv_Weibull((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]);
                                      *res+= weights*log(bl);
                         }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_LogGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                              if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                  bl=d2lognorm(zi,zj,sill,nugget, mean1[i], mean2[i],corr);
                             *res+= weights*log(bl);
                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Pair_LogLogistic_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{

    int i=0;
    double bl,corr,zi,zj,weights=1.0,nugget=0.0;
      nugget=nuis[0];
     if(nugget<0||nugget>=1||nuis[2]<=2) {*res=LOW;  return;}

       for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                            zi=data1[i];zj=data2[i];

                                corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                bl=biv_LogLogistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]);
                                *res+= weights*log(bl);
                            }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/

void Comp_Pair_Logistic_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{

    int i=0;
    double bl,corr=0.0,zi=0.0,zj=0.0,weights=1.0;
    if( nuis[1]<=0)  {*res=LOW;  return;}

    double sill=1-nuis[0];

          for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                            zi=data1[i];zj=data2[i];

                                corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                bl=biv_Logistic(sill*corr,zi,zj,mean1[i],mean2[i],nuis[1]);
                                *res+= weights*log(bl);
                            }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}



void Comp_Pair_TWOPIECEBIMODAL_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
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
                                //p11=pbnorm(cormod,lags[i],lagt[i],qq,qq,nugget,1,par,0);
                                p11=pbnorm22(qq,qq,corr);
                                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                bl=biv_two_piece_bimodal((1-nugget)*corr,zi,zj,sill,df,delta,eta,p11,mean1[i],mean2[i]);
                             *res+= weights*log(bl);
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
    double *nuis, int *local,int *GPU)
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
    double *nuis, int *local,int *GPU)

{
    int i=0;
    double u=0.0, w=0.0, rhotv=0.0,weights=1.0;
    int N=2;
      double *vari;vari=(double *) R_Calloc(N,double);vari[0]=par[0];vari[1]=par[1];  /// variances of the skew gaussian
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

        R_Free(vari);
    if(!R_FINITE(*res))*res = LOW;
    return;
}














































/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/********************* SPATIAL Gaussian and Clayton COPULA *****************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/

// Composite marginal (pariwise) log-likelihood for the spatial Gaussian model:
void Comp_Pair_GaussCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
     int model=1; 
    /*############*/
    int i=0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]) ;
                      *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

// Composite marginal (pariwise) log-likelihood for the spatial Gaussian model:
void Comp_Pair_TCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)

{
     int model=12;
    /*############*/
    int i=0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]) ;
                      *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}




void Comp_Pair_BetaCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
     int model=28; 
    /*############*/
    int i=0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]) ;
                        *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Pair_Beta2Cop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop,int *cond)
{
    int model=50; 
    /*############*/
    int i=0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]) ;
                       *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Pair_KumaraswamyCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop,int *cond)
{
     int model=33; 
    /*############*/
    int i=0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]) ;
                        *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}





void Comp_Pair_WeibullCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop,int *cond)
{
    int model=26; 
    /*############*/
   
    double  weights=1.0,nugget,corr,bl; int i;
    nugget=nuis[0];
    //if(nugget<0||nugget>1){*res=LOW; return;}
     if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
//Rprintf("%d -----\n", npairs[0]);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       //if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]) ;
                       *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}


void Comp_Pair_GammaCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop,int *cond)
{
    int model=21; 
    /*############*/
  
    double  weights=1.0,nugget,corr,bl;  int i;
    nugget=nuis[0];
    if(nugget<0||nugget>1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                      // if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]) ;
                       *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Pair_LogGaussCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop,int *cond)
{
  int model=24; 
    /*############*/
    int i=0;
    double corr,weights=1.0,bl;
    double sill=nuis[1];double nugget=nuis[0];

    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]) ;
                       *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

     
void Comp_Pair_PoisCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop,int *cond)

{
    int model=30; 
    /*############*/
    int i=0;
    double  weights=1.0,nugget,corr,bl;
    nugget=nuis[0];
    if(nugget<0||nugget>1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]) ;
                       *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}



void Comp_Pair_BinomNNGaussCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop,int *cond)

{
    int model=11; 
    /*############*/
    int i=0;
    double  weights=1.0,nugget,corr,bl;
    nugget=nuis[0];
    if(nugget<0||nugget>1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]) ;
                       *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}



void Comp_Pair_BinomnegGaussCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop,int *cond)

{
    int model=16; 
    /*############*/
    int i=0;
    double  weights=1.0,nugget,corr,bl;
    nugget=nuis[0];
    if(nugget<0||nugget>1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]) ;
                       *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}


void Comp_Pair_LogisticCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop,int *cond)

{
     int model=25; 
    /*############*/
    int i;
    double bl,corr,weights=1.0,nugget;
        nugget=nuis[0];
    if(nugget>=1||nugget<0.0 ) {*res=LOW;  return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
            bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]) ;
                       *res+= bl*weights;
                }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
