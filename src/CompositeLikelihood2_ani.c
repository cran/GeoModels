#include "header.h"
/******************************************************************************************/
/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
/******************************************************************************************/


// Composite marginal (pariwise) log-likelihood for the spatial Gaussian model:
void Comp_Pair_Gauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(!R_FINITE(sill) || sill<=0 || nugget<0 || nugget>=1){*res=LOW; return;}
    for(i=0;i<(npairs[0]);i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
  lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                      corr=CorFct(cormod,lag,0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                      bl=log_biv_Norm((1-nugget)*corr,data1[i],data2[i],mean1[i],mean2[i],sill,0);
                      *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/******************************************************************************************/
void Comp_Diff_Gauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
  int i=0;double lag=0.0;
  double vario=0.0,u,v,weights=1.0;

 double nugget=nuis[0];
    double sill=nuis[1];

      if(sill<0 || nugget<0||nugget>=1){*res=LOW; return;}

 for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
  vario=Variogram(cormod,lag,0,nuis[0],nuis[1],par);
         u=data1[i];v=data2[i];
            if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
    *res+= -0.5*(log(2*M_PI)+log(vario)+
                   R_pow(u-v,2)/(2*vario))*weights;}}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_WrapGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;
    double  u=0.0,v=0.0,weights=1.0,corr=0.0;
    double wrap_gauss;
    double alfa=3.0;
     double nugget=nuis[0];
    double sill=nuis[1];

      if(sill<=0 || nugget<0||nugget>=1){*res=LOW; return;}
 for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                u=data1[i];
                v=data2[i];
                corr=(1.0-nugget)*clamp_corr(CorFct(cormod,lag,0,par,0,0));
                wrap_gauss=biv_wrapped_log(alfa,u,v,mean1[i],mean2[i],nuis[0],nuis[1],corr);
                    if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                    *res+=wrap_gauss*weights ;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_SinhGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i=0;double lag=0.0;double corr,zi,zj,bb=0.0,weights=1.0;
           if(nuis[3]<0||nuis[1]<0||nuis[0]<0||nuis[0]>=1) {*res=LOW;  return;}

   for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=data1[i];zj=data2[i];
                      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    corr=CorFct(cormod,lag,0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                    bb=log(biv_sinh((1-nuis[0])*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],nuis[1]));
                    *res+= weights*bb;
                 }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_SkewGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

      double sill=nuis[1];
      double nugget=nuis[0];
     if(nugget<0|| nugget>=1||sill<0){*res=LOW;  return;}

    int i=0;double lag=0.0;double corr,zi,zj,weights=1.0;
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=data1[i];zj=data2[i];
                    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    corr=CorFct(cormod,lag,0,par,0,0);
                             if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                     *res+= weights*log(biv_skew(corr,zi,zj,mean1[i],mean2[i],sill,nuis[2],nuis[0]));
                 }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Pair_SkewLaplace2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1, double *data2, int *N1, int *N2,
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
    const double max_dist = maxdist[0];
    const double scale = 1.0 - nugget;


    const double sill = nuis[1];
    double total = 0.0; 

    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Calcolo correlazione
             const double lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            
            // Calcolo pesi se necessario
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            
            // log bivariate skew laplace mean1[i], mean2[i]
           const double val=log_biv_skewlaplace(scale*corr,d1,d2,mean1[i],mean2[i],sill, skew);
            total += weights * val;
        }
    }
    *res = R_FINITE(total) ? total : LOW;
}

/*********************************************************/

void Comp_Pair_Gamma2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i;double corr,zi,zj,weights=1.0,bl=0.0,lag=0.0;
      double nugget=nuis[0];
     if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}

     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=(data1[i]);zj=(data2[i]);
                lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    corr=CorFct(cormod,lag,0,par,0,0);
                        if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);

                    bl=biv_gamma((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]);
                    //Rprintf("%f\n",bl);
                     *res+= weights*log(bl);

                }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_Weibull2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i;double corr,zi,zj,weights=1.0,bl=0.0,lag=0.0;
      double nugget=nuis[0];
     if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}

     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=(data1[i]);zj=(data2[i]);
                lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    corr=CorFct(cormod,lag,0,par,0,0);
                        if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);

                    bl=biv_Weibull((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]);
                    //Rprintf("%f\n",bl);
                     *res+= weights*log(bl);

                }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_Kumaraswamy2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i;double corr,zi,zj,weights=1.0,bl,lag=0.0;
    double nugget=nuis[0];
         double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    corr=CorFct(cormod,lag,0,par,0,0);
                     if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                  bl=biv_Kumara((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);

        *res+= weights*log(bl);
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}


/*********************************************************/
void Comp_Pair_Kumaraswamy22mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i;double corr,zi,zj,weights=1.0,bl,lag=0.0;
    double nugget=nuis[0];
         double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    corr=CorFct(cormod,lag,0,par,0,0);
                     if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                  bl=biv_Kumara2((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);

        *res+= weights*log(bl);
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Beta2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i;double corr,zi,zj,weights=1.0,bl,lag=0.0;
    double nugget=nuis[0];
    double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max)  {*res=LOW;  return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    corr=CorFct(cormod,lag,0,par,0,0);
                     if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                  bl=biv_beta((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);

        *res+= weights*log(bl);
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_LogGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i=0;double lag=0.0;double corr,zi,zj,weights=1.0,bl=0.0;
    double sill=nuis[1];double nugget=nuis[0];

    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=data1[i];zj=data2[i];
                    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    corr=CorFct(cormod,lag,0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                    bl=d2lognorm(zi,zj,sill,nugget, mean1[i], mean2[i],(1-nugget)*corr);
                    *res+= weights*log(bl);
                    }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_PoisbinnegGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0,  uu=0,vv=0;
    double bl,u,v,weights=1.0,ai=0.0,aj=0.0,corr=0.0,lag=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
     double nugget=nuis[0];
     if(nugget<0||nugget>=1){*res=LOW; return;}

        for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                   ai=mean1[i];aj=mean2[i];
                   lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    corr=CorFct(cormod,lag,0,par,0,0);
                    p11=pbnorm22(ai,aj,(1-nugget)*corr);

                p1=pnorm(ai,0,1,1,0);p2=pnorm(aj,0,1,1,0);
                    u=data1[i];v=data2[i];
    if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                          uu=(int) u;  vv=(int) v;
                        bl=biv_poisbinneg(N1[0],uu,vv,p1,p2,p11);
                    *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_PoisbinGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,lag=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
 double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                   ai=mean1[i];aj=mean2[i];
                   corr=CorFct(cormod,lag,0,par,0,0);
                   lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    p11=pbnorm22(ai,aj,(1-nugget)*corr);
                p1=pnorm(ai,0,1,1,0);p2=pnorm(aj,0,1,1,0);
                    u=data1[i];v=data2[i];
                       if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                          uu=(int) u; vv=(int) v;
                        bl=biv_poisbin(N1[0],uu,vv,p1,p2,p11);

                    *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}


/*********************************************************/
void Comp_Pair_BinomnegGauss2mem_aniso(int *cormod, double *coord1, double *coord2,
    double *data1, double *data2, int *N1, int *N2,
    double *par, int *weigthed, double *res, double *mean1, double *mean2,
    double *nuis, int *type_cop, int *cond)
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
 
#define P_CLAMP 1e-10
 
    int    i, uu, vv, valid_pairs = 0;
    double u, v, bl, weights, ai, aj, corr, lag, p1, p2, p11;
 
    for (i = 0; i < npairs[0]; i++) {
 
        if (ISNAN(data1[i]) || ISNAN(data2[i])) continue;
 
        /* lag calcolato dalle coordinate anisotrope (invariato rispetto originale) */
        lag  = hypot(coord1[2*i]   - coord2[2*i],
                     coord1[2*i+1] - coord2[2*i+1]);
 
        weights = 1.0;
        if (*weigthed) {
            weights = CorFunBohman(lag, maxdist[0]);
            if (weights <= 0.0) continue;
        }
 
        corr = clamp_corr(CorFct(cormod, lag, 0, par, 0, 0));
        ai   = mean1[i];
        aj   = mean2[i];
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
 
        /* FIX 1: passa la cache allocata una sola volta */
        bl = biv_binomneg(N, uu, vv, p1, p2, p11, lgamma_cache, cache_size);
 
        if (!R_FINITE(bl) || bl <= 0.0) continue;
        valid_pairs++;
        *res += weights * log(bl);
    }
 
 
#undef P_CLAMP
 
    /* FIX 2: soglia robusta per npairs piccolo */
    if (valid_pairs < (npairs[0] + 1) / 2) { *res = LOW; return; }
    if (!R_FINITE(*res)) *res = LOW;
}
/*********************************************************
void Comp_Pair_BinomnegLogi2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,lag=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
       if( nugget>=1 || nugget<0){*res=LOW; return;}
    //compute the composite log-likelihood:

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                  lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                 corr=CorFct(cormod,lag,0,par,0,0);

                    p11=pblogi22(ai,aj,(1-nugget)*corr);
                    p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));

                    u=data1[i];v=data2[i];
                         if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                          uu=(int) u;
                         vv=(int) v;
                        bl=biv_binomneg (N1[0],uu,vv,p1,p2,p11);
                    *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}*/
/********************************************************
void Comp_Pair_BinomnegGaussZINB2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,lag=0.0;
    double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];


      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                  lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                 corr=CorFct(cormod,lag,0,par,0,0);
                    u=data1[i];v=data2[i];
                         if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                          uu=(int) u;
                         vv=(int) v;
                        bl=biv_binomnegZINB(N1[0],corr,uu,vv,ai,aj,nugget1,nugget2,mup);
                    *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}*/
/******************************************************************************************/
void Comp_Pair_BinomGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,lag=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                 corr=CorFct(cormod,lag,0,par,0,0);
                    p11=pbnorm22(ai,aj,(1-nugget)*corr);
                p1=pnorm(ai,0,1,1,0);
                p2=pnorm(aj,0,1,1,0);
                u=data1[i];v=data2[i];
                        if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                          uu=(int) u; vv=(int) v;
                        bl=biv_binom (N1[0],uu,vv,p1,p2,p11);
                    *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_BinomLogi2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis ,int *type_cop, int *cond)
{
    int i=0, uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,lag=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 u=data1[i];v=data2[i];
                 lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                 corr=CorFct(cormod,lag,0,par,0,0);
                p11=pblogi22(ai,aj,(1-nugget)*corr);
                 p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                        if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                          uu=(int) u; vv=(int) v;
                        bl=biv_binom (N1[0],uu,vv,p1,p2,p11);
                    *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Pair_BinomNNGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu=0,vv=0,n1,n2;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,lag=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                 corr=CorFct(cormod,lag,0,par,0,0);
                 p11=pbnorm22(ai,aj,(1-nugget)*corr);
                 p1=pnorm(ai,0,1,1,0);
                 p2=pnorm(aj,0,1,1,0);
                 u=data1[i];v=data2[i];
                 n1=N1[i];n2=N2[i];
                 if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                 uu=(int) u; vv=(int) v;
                 bl=biv_binom222(n1,n2,uu,vv,p1,p2,p11);
                 *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}



void Comp_Pair_BinomNNGauss_misp2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, N=2,n1,n2;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,lag=0.0;
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
                 lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                 corr=CorFct(cormod,lag,0,par,0,0);
                 p11=pbnorm22(ai,aj,(1-nugget)*corr);
                 p1=pnorm(ai,0,1,1,0);
                 p2=pnorm(aj,0,1,1,0);
                 u=data1[i];v=data2[i];
                 n1=N1[i];n2=N2[i];
                 if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                 M[0][0]=n1*p1*(1-p1);   M[1][1]=n2*p2*(1-p2);  // var1 var2
                 M[0][1]= fmin_int(n1,n2)*(p11-p1*p2) ;       // covariance
                 M[1][0]= M[0][1];
                 dat[0]=u-n1*p1;dat[1]=v-n2*p2; 
                 //Rprintf("%d %f %f %f \n",fmin_int(n1,n2),p1,p2,p11 );
                   //#####
                 bl=log_dNnorm(N,M,dat);
                 *res += bl*weights;       
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}




void Comp_Pair_BinomNNLogi2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu=0,vv=0,n1,n2;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,lag=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                 corr=CorFct(cormod,lag,0,par,0,0);

                 p11=pblogi22(ai,aj,(1-nugget)*corr);
                 p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                 u=data1[i];v=data2[i];
                 n1=N1[i];n2=N2[i];
                 if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                 uu=(int) u; vv=(int) v;
                 bl=biv_binom222(n1,n2,uu,vv,p1,p2,p11);
                 *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_LogLogistic2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i;double bl,corr,zi,zj,weights=1.0,nugget=0.0,lag=0.0;
     nugget=nuis[0];
     if(nugget<0||nugget>=1||nuis[2]<=2) {*res=LOW;  return;}
   for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    corr=CorFct(cormod,lag,0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                    bl=biv_LogLogistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]);
                    *res+= weights*log(bl);
                }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Logistic2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,nugget=0.0,lag=0.0;
        nugget=nuis[0];
    if(nugget>=1||nugget<0.0 ) {*res=LOW;  return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=(data1[i]);zj=(data2[i]);
                lag=hypot(coord1[2*i]-coord2[2*i],
                          coord1[2*i+1]-coord2[2*i+1]);
                corr=CorFct(cormod,lag,0,par,0,0);
                if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                    bl= biv_Logistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[1]);
                    *res+= weights*log(bl);
                }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Pois2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu,ww;
    double weights=1.0,corr,mui,muj,bl,lag=0.0;
    double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}
  // Rprintf("%d   \n",npairs[0]);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                     corr=CorFct(cormod,lag,0,par,0,0);
                    // if(fabs(corr)>1|| !R_FINITE(corr)) {*res=LOW; return;}
                        if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                      uu=(int) data1[i];  ww=(int) data2[i];
                      bl=biv_Poisson((1-nugget)*corr,uu,ww,mui, muj);
                      *res+= bl*weights;
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_PoisGamma2mem_aniso(int *cormod, double *coord1, double *coord2,
                                   double *data1, double *data2, int *N1, int *N2,
                                   double *par, int *weigthed, double *res,
                                   double *mean1, double *mean2, double *nuis,
                                   int *type_cop, int *cond)
{
    (void)N1; (void)N2; (void)type_cop; (void)cond;
    const double nugget = nuis[0], shape = nuis[2];
    if (!R_FINITE(nugget) || nugget < 0.0 || nugget >= 1.0 ||
        !R_FINITE(shape) || !(shape > 0.0) || !npairs) { *res = LOW; return; }
    GeoKahanSum total; GeoKahanInit(&total); int used = 0;
    for (int i = 0; i < npairs[0]; ++i) {
        int r, t; double mi, mj;
        const int prep = GeoPGPreparePairAt(i, data1[i], data2[i],
                                             mean1[i], mean2[i],
                                             &r, &t, &mi, &mj);
        if (prep == 0) continue;
        if (prep < 0) { *res = LOW; return; }
        const double lag = hypot(coord1[2*i] - coord2[2*i],
                                 coord1[2*i+1] - coord2[2*i+1]);
        double weight = 1.0;
        if (weigthed[0]) weight = CorFunBohman(lag, maxdist[0]);
        if (!R_FINITE(weight)) { *res = LOW; return; }
        if (weight <= 0.0) continue;
        double base = CorFct(cormod, lag, 0.0, par, 0, 0);
        if (!R_FINITE(base) || fabs(base) > 1.0 + 1e-8) { *res = LOW; return; }
        double rho = lag <= 0.0 ? fmax2(-1.0, fmin2(1.0, base))
                                : (1.0 - nugget) * base;
        if (lag > 0.0 && rho >= 1.0) rho = 1.0 - 1e-15;
        if (lag > 0.0 && rho <= -1.0) rho = -1.0 + 1e-15;
        const double ll = GeoPGLogJoint(rho, r, t, mi, mj, shape);
        if (!R_FINITE(ll)) { *res = LOW; return; }
        GeoKahanAdd(&total, weight * ll); ++used;
    }
    *res = (used > 0 && R_FINITE(total.sum)) ? total.sum : LOW;
}



void Comp_Pair_Gauss_misp_PoisGamma2mem_aniso(int *cormod, double *coord1,
                                              double *coord2, double *data1,
                                              double *data2, int *N1, int *N2,
                                              double *par, int *weigthed,
                                              double *res, double *mean1,
                                              double *mean2, double *nuis,
                                              int *type_cop, int *cond)
{
    (void)N1; (void)N2; (void)type_cop; (void)cond;
    const double nugget = nuis[0], shape = nuis[2];
    if (!R_FINITE(nugget) || nugget < 0.0 || nugget >= 1.0 ||
        !R_FINITE(shape) || !(shape > 0.0) || !npairs) { *res = LOW; return; }
    GeoKahanSum total; GeoKahanInit(&total); int used = 0;
    for (int i = 0; i < npairs[0]; ++i) {
        if (ISNAN(data1[i]) || ISNAN(data2[i])) continue;
        if (!R_FINITE(data1[i]) || !R_FINITE(data2[i]) ||
            !R_FINITE(mean1[i]) || !R_FINITE(mean2[i]) ||
            mean1[i] > 700.0 || mean2[i] > 700.0) { *res = LOW; return; }
        const double mi = exp(mean1[i]), mj = exp(mean2[i]);
        const double lag = hypot(coord1[2*i] - coord2[2*i],
                                 coord1[2*i+1] - coord2[2*i+1]);
        double weight = weigthed[0] ? CorFunBohman(lag, maxdist[0]) : 1.0;
        if (!R_FINITE(weight)) { *res = LOW; return; }
        if (weight <= 0.0) continue;
        double base = CorFct(cormod, lag, 0.0, par, 0, 0);
        if (!R_FINITE(base) || fabs(base) > 1.0 + 1e-8) { *res = LOW; return; }
        double rho = lag <= 0.0 ? fmax2(-1.0, fmin2(1.0, base))
                                : (1.0 - nugget) * base;
        if (rho >= 1.0) rho = 1.0 - 1e-15;
        if (rho <= -1.0) rho = -1.0 + 1e-15;
        const double ll = GeoPGGaussLogJoint(rho, data1[i], data2[i], mi, mj, shape);
        if (!R_FINITE(ll)) { *res = LOW; return; }
        GeoKahanAdd(&total, weight * ll); ++used;
    }
    *res = (used > 0 && R_FINITE(total.sum)) ? total.sum : LOW;
}



/*********************************************************/
void Comp_Pair_PoisZIP2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu,ww;
    double weights=1.0,corr,mui,muj,bl,lag=0.0;
   double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];


      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
  // Rprintf("%d   \n",npairs[0]);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                     corr=CorFct(cormod,lag,0,par,0,0);
                    // if(fabs(corr)>1|| !R_FINITE(corr)) {*res=LOW; return;}
                        if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                      uu=(int) data1[i];  ww=(int) data2[i];
                      bl=biv_PoissonZIP(corr,uu,ww,mui, muj,mup,nugget1,nugget2);
                      *res+= log(bl)*weights;
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Gauss_misp_Pois2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0,N=2;
    double  weights=1.0,corr,corr1,mui,muj,bl,lag=0.0;
    double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}
    double M0[2] = {0.0, 0.0};
    double M1[2] = {0.0, 0.0};
    double *M[2] = {M0, M1};
    double dat[2] = {0.0, 0.0};
    for(i=0;i<npairs[0];i++){

                  if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
             //***********/
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                     corr=CorFct(cormod,lag,0,par,0,0)*(1-nugget);
                      corr1=corr_pois(corr,mui, muj);
                      if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                        M[0][0]=mui; M[1][1]=muj;M[0][1]=sqrt(mui*muj)*corr1;M[1][0]= M[0][1];
                        dat[0]=data1[i]-mui;dat[1]=data2[i]-muj;

                      bl=log_dNnorm(N,M,dat);
                      *res += bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Pair_Gauss_misp_PoisZIP2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;
    double weights=1.0,corr,mui,muj,bl;
    double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];

      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}

  // Rprintf("%d   \n",npairs[0]);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                     corr=CorFct(cormod,lag,0,par,0,0);
                      if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                      bl=biv_Mis_PoissonZIP(corr,data1[i],data2[i],mui, muj,mup,nugget1,nugget2);
                  //    Rprintf("%f %f\n",bl,mup);
                      *res+= log(bl)*weights;
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_Gauss_misp_SkewT2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;
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
    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                     corr=clamp_corr(CorFct(cormod,lag,0,par,0,0));
                     corr2=GeoSkewStudentCorrFromLogFactor(corr,(1-nugget)*corr,df,skew,log_factor);
                     if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                          bl=log_biv_Norm(corr2,data1[i],data2[i],mean1[i]+sqrt(sill)*MM,
                                                                 mean2[i]+sqrt(sill)*MM,sill*FF,0);
                        *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Gauss_misp_T2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
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
            const double lag = hypot(coord1[2*i]-coord2[2*i], coord1[2*i+1]-coord2[2*i+1]);
            const double weights = *weigthed ? CorFunBohman(lag, maxdist[0]) : 1.0;
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
void Comp_Pair_Tukeyhh2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,lag=0.0;
    double sill=nuis[1];
    double nugget=nuis[0];
    double h1=nuis[3];
    double h2=nuis[2];
      if( sill<=0||h1<0||h1>=0.5||h2<0||h2>=0.5||nugget<0||nugget>=1){*res=LOW; return;}
        for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    corr=CorFct(cormod,lag,0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                   bl=biv_tukey_hh((1-nugget)*corr,zi,zj,mean1[i],mean2[i],sill,h1,h2);
                             *res+= weights*log(bl);
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}





/*********************************************************/
void Comp_Pair_Tukeyh2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,lag=0.0;
    double sill=nuis[1];
    double nugget=nuis[0];
    double tail=nuis[2];
      if( sill<0||tail<0||tail>0.5||nugget<0||nugget>=1){*res=LOW; return;}
        for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    corr=CorFct(cormod,lag,0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                   bl=biv_tukey_h((1-nugget)*corr,zi,zj,mean1[i],mean2[i],tail,sill);
                             *res+= weights*log(bl);
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Gauss_misp_Tukeygh2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i;double bl,corr,corr2,zi,zj,weights=1.0,eta,tail,sill,nugget,u,eta2,mu,vv,lag=0.0;
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
lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    corr=(1-nugget)*CorFct(cormod,lag,0,par,0,0);
                    corr2=corr_tukeygh(corr,eta,tail);
                 //   if(corr2<0) Rprintf("%f %f %f \n",corr2,par[0],par[1]);
                    if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                bl=log_biv_Norm(corr2,zi,zj,mean1[i]+sqrt(sill)*mu,
                                            mean2[i]+sqrt(sill)*mu, sill*vv,0);
            //  Rprintf("%f %f-- %f %f \n",mean1[i],mean2[i],zi,zj);
                    *res+= weights*bl;
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_T2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,lag=0.0;
    double sill=nuis[2];
    double nugget=nuis[1];
    double df=nuis[0];
      if( sill<0||nugget<0||nugget>=1||df<0||df>0.5){*res=LOW; return;}

   for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    corr=CorFct(cormod,lag,0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                   bl=biv_T(corr,(zi-mean1[i])/sqrt(sill),
                                            (zj-mean2[i])/sqrt(sill),df,nugget)/sill;
                             *res+= weights*log(bl);
                }}


    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_TWOPIECETukeyh2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,eta,tail,qq,sill,nugget,lag=0.0;
    eta  = nuis[2];  //skewness parameter
    tail = nuis[3];  //tail parameter
    sill =nuis[1];
    nugget=nuis[0];

     if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1||tail<0||tail>0.5) {*res=LOW;  return;}

       qq=qnorm((1-eta)/2,0,1,1,0);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                 corr=CorFct(cormod,lag,0,par,0,0);
                      p11=pbnorm22(qq,qq,corr);
                    if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                    bl=biv_two_pieceTukeyh((1-nugget)*corr,zi,zj,sill,eta,tail,p11,mean1[i],mean2[i]);
                    *res+= weights*log(bl);
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_TWOPIECET2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,qq,lag=0.0;

    double eta=nuis[3];  //skewness parameter
    double df=nuis[0];
    double sill=nuis[2];
    double nugget=nuis[1];
    if(sill<0||nugget<0||nugget>=1 ||fabs(eta)>1|| df >0.5||df<0) {*res=LOW;  return;}

       qq=qnorm((1-eta)/2,0,1,1,0);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    corr=CorFct(cormod,lag,0,par,0,0);
                      p11=pbnorm22(qq,qq,corr);
                    if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                    /********************************************************/
    bl=biv_two_pieceT(corr,zi,zj,sill,df,eta,p11,mean1[i],mean2[i],nugget);
  //Rprintf("%f- %f- %f %f %f %f %f %f %f %f %f  \n",lag,corr,eta,sill,zi,zj,df,eta,p11,mean1[i],mean2[i]);
                    /********************************************************/
                           *res+= weights*log(bl);
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_TWOPIECEGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i;double bl,corr,zi,zj,weights=1.0,p11,eta,qq,sill,nugget,lag=0.0;
    eta=nuis[2];  //skewness parameter
    sill=nuis[1];
    nugget=nuis[0];
       qq=qnorm((1-eta)/2,0,1,1,0);
         if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1) {*res=LOW;  return;}
   for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    corr=CorFct(cormod,lag,0,par,0,0);
                      p11=pbnorm22(qq,qq,corr);
                       if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                    bl=biv_two_pieceGaussian((1-nugget)*corr,zi,zj,sill,eta,p11,mean1[i],mean2[i]);
                    *res+= weights*log(bl);
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}

void Comp_Pair_TWOPIECEBIMODAL2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,eta,qq,sill,df,nugget,delta,lag=0.0;

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
                    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                    corr=CorFct(cormod,lag,0,par,0,0);
                        p11=pbnorm22(qq,qq,corr);
                    if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                    /********************************************************/
                   bl=biv_two_piece_bimodal((1-nugget)*corr,zi,zj,sill,df,delta,eta,p11,mean1[i],mean2[i]);
          // Rprintf("%f %f  --%f  %f %f %f %f  -%f %f \n",bl,lag,df,delta,eta,sill,corr,par[0],par[1]);
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
































/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/********************* SPATIAL Gaussian and clayton COPULA *****************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/

// Composite marginal (pariwise) log-likelihood for the spatial Gaussian model:
void Comp_Pair_GaussCop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
  int model=1; 
    /*############*/
    int i=0;double lag=0.0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}


    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                      corr=clamp_corr(CorFct(cormod,lag,0,par,0,0));
                       if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
           bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]);
                      *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Pair_TCop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
  int model=12;
    /* StudentT nuisance order: df, nugget, sill, [nu]. */
    int i=0;double lag=0.0;
    double  weights=1.0,sill,nugget,corr,bl,df;
    df=nuis[0];nugget=nuis[1];sill=nuis[2];
    if(!R_FINITE(df)||df<=0||!R_FINITE(sill)||sill<=0||nugget<0||nugget>=1){*res=LOW; return;}


    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                      corr=clamp_corr(CorFct(cormod,lag,0,par,0,0));
                       if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
           bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]);
                      *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Pair_BetaCop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
   int model=28; 
    /*############*/
    int i=0;double lag=0.0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(!R_FINITE(sill)||sill<=0 || nugget<0||nugget>=1 ||
       !R_FINITE(nuis[2])||!R_FINITE(nuis[3])||nuis[2]<=0||nuis[3]<=0 ||
       !R_FINITE(nuis[4])||!R_FINITE(nuis[5])||nuis[4]>=nuis[5]) {*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                      corr=clamp_corr(CorFct(cormod,lag,0,par,0,0));
                       if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
           bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]);
                        *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Pair_Beta2Cop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
     int model=50; 
    /*############*/
    int i=0;double lag=0.0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(!R_FINITE(sill)||sill<=0 || nugget<0||nugget>=1 ||
       !R_FINITE(nuis[2])||nuis[2]<=0 ||
       !R_FINITE(nuis[3])||!R_FINITE(nuis[4])||nuis[3]>=nuis[4]) {*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                      corr=clamp_corr(CorFct(cormod,lag,0,par,0,0));
                       if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
           bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]);
                       *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Pair_KumaraswamyCop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
   int model=33; 
    /*############*/
    int i=0;double lag=0.0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(!R_FINITE(sill)||sill<=0 || nugget<0||nugget>=1 ||
       !R_FINITE(nuis[2])||!R_FINITE(nuis[3])||nuis[2]<=0||nuis[3]<=0 ||
       !R_FINITE(nuis[4])||!R_FINITE(nuis[5])||nuis[4]>=nuis[5]) {*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                      corr=clamp_corr(CorFct(cormod,lag,0,par,0,0));
                       if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
           bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]);
                        *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
void Comp_Pair_Kumaraswamy2Cop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int model=42; 
    /*############*/
    int i=0;double lag=0.0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(!R_FINITE(sill)||sill<=0 || nugget<0||nugget>=1 ||
       !R_FINITE(nuis[2])||nuis[2]<=0 ||
       !R_FINITE(nuis[3])||!R_FINITE(nuis[4])||nuis[3]>=nuis[4]) {*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                      corr=clamp_corr(CorFct(cormod,lag,0,par,0,0));
                       if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
           bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]);
                      *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Pair_WeibullCop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    const int model=26;
    const double nugget=nuis[0], shape=nuis[2];
    if(nugget<0||nugget>=1||!R_FINITE(shape)||shape<=0){*res=LOW;return;}
    double total=0.0;
    for(int i=0;i<npairs[0];++i){
        if(!ISNAN(data1[i])&&!ISNAN(data2[i])){
            double lag=hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
            double corr=clamp_corr(CorFct(cormod,lag,0,par,0,0));
            double w=*weigthed?CorFunBohman(lag,maxdist[0]):1.0;
            total += w*biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]);
        }
    }
    *res=R_FINITE(total)?total:LOW;
}

void Comp_Pair_GammaCop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    const int model=21;
    const double nugget=nuis[0], shape=nuis[2];
    if(nugget<0||nugget>=1||!R_FINITE(shape)||shape<=0){*res=LOW;return;}
    double total=0.0;
    for(int i=0;i<npairs[0];++i){
        if(!ISNAN(data1[i])&&!ISNAN(data2[i])){
            double lag=hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
            double corr=clamp_corr(CorFct(cormod,lag,0,par,0,0));
            double w=*weigthed?CorFunBohman(lag,maxdist[0]):1.0;
            total += w*biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]);
        }
    }
    *res=R_FINITE(total)?total:LOW;
}

void Comp_Pair_LogGaussCop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    const int model=22;
    const double nugget=nuis[0], sill=nuis[1];
    if(nugget<0||nugget>=1||!R_FINITE(sill)||sill<=0){*res=LOW;return;}
    double total=0.0;
    for(int i=0;i<npairs[0];++i){
        if(!ISNAN(data1[i])&&!ISNAN(data2[i])){
            double lag=hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
            double corr=clamp_corr(CorFct(cormod,lag,0,par,0,0));
            double w=*weigthed?CorFunBohman(lag,maxdist[0]):1.0;
            total += w*biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]);
        }
    }
    *res=R_FINITE(total)?total:LOW;
}

void Comp_Pair_LogisticCop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    const int model=25;
    const double nugget=nuis[0], sill=nuis[1];
    if(nugget<0||nugget>=1||!R_FINITE(sill)||sill<=0){*res=LOW;return;}
    double total=0.0;
    for(int i=0;i<npairs[0];++i){
        if(!ISNAN(data1[i])&&!ISNAN(data2[i])){
            double lag=hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
            double corr=clamp_corr(CorFct(cormod,lag,0,par,0,0));
            double w=*weigthed?CorFunBohman(lag,maxdist[0]):1.0;
            total += w*biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]);
        }
    }
    *res=R_FINITE(total)?total:LOW;
}

void Comp_Pair_SkewLaplaceCop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    const int model=59;
    const double nugget=nuis[0], sill=nuis[1], skew=nuis[2];
    if(nugget<0||nugget>=1||!R_FINITE(sill)||sill<=0||!R_FINITE(skew)||skew<=0||skew>=1){*res=LOW;return;}
    double total=0.0;
    for(int i=0;i<npairs[0];++i){
        if(!ISNAN(data1[i])&&!ISNAN(data2[i])){
            double lag=hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
            double corr=clamp_corr(CorFct(cormod,lag,0,par,0,0));
            double w=*weigthed?CorFunBohman(lag,maxdist[0]):1.0;
            total += w*biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]);
        }
    }
    *res=R_FINITE(total)?total:LOW;
}

void Comp_Pair_BinomNNGaussCop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    const int model=11;
    const double nugget=nuis[0];
    if(nugget<0||nugget>=1){*res=LOW;return;}
    double total=0.0;
    for(int i=0;i<npairs[0];++i){
        if(!ISNAN(data1[i])&&!ISNAN(data2[i])){
            double lag=hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
            double corr=clamp_corr(CorFct(cormod,lag,0,par,0,0));
            double w=*weigthed?CorFunBohman(lag,maxdist[0]):1.0;
            total += w*biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]);
        }
    }
    *res=R_FINITE(total)?total:LOW;
}

void Comp_Pair_BinomnegGaussCop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    const int model=16;
    const double nugget=nuis[0];
    if(nugget<0||nugget>=1){*res=LOW;return;}
    double total=0.0;
    for(int i=0;i<npairs[0];++i){
        if(!ISNAN(data1[i])&&!ISNAN(data2[i])){
            double lag=hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
            double corr=clamp_corr(CorFct(cormod,lag,0,par,0,0));
            double w=*weigthed?CorFunBohman(lag,maxdist[0]):1.0;
            total += w*biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]);
        }
    }
    *res=R_FINITE(total)?total:LOW;
}

void Comp_Pair_PoisCop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
   int model=30; 
    /*############*/
    int i=0;double lag=0.0;
    double  weights=1.0,nugget,corr,bl;
    nugget=nuis[0];
    if(nugget<0||nugget>=1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]);
                      corr=clamp_corr(CorFct(cormod,lag,0,par,0,0));
                       if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
           bl=biv_cop(corr,type_cop[0],cond[0],data1[i],data2[i],mean1[i],mean2[i],nuis,model,N1[i],N2[i]);
                        *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
void Comp_Pair_PoisGammaZIP2mem_aniso(int *cormod, double *coord1, double *coord2,
                                      double *data1, double *data2, int *N1, int *N2,
                                      double *par, int *weigthed, double *res,
                                      double *mean1, double *mean2, double *nuis,
                                      int *type_cop, int *cond)
{
    (void)N1; (void)N2; (void)type_cop; (void)cond;
    const double nugget1 = nuis[0], nugget2 = nuis[1];
    const double mup = nuis[2], shape = nuis[4];
    if (!R_FINITE(nugget1) || !R_FINITE(nugget2) ||
        nugget1 < 0.0 || nugget1 >= 1.0 || nugget2 < 0.0 || nugget2 >= 1.0 ||
        !R_FINITE(mup) || !R_FINITE(shape) || !(shape > 0.0) || !npairs) {
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
        const double lag = hypot(coord1[2*i] - coord2[2*i],
                                 coord1[2*i+1] - coord2[2*i+1]);
        double weight = weigthed[0] ? CorFunBohman(lag, maxdist[0]) : 1.0;
        if (!R_FINITE(weight)) { *res = LOW; return; }
        if (weight <= 0.0) continue;
        double rho = CorFct(cormod, lag, 0.0, par, 0, 0);
        if (!R_FINITE(rho) || fabs(rho) > 1.0 + 1e-8) { *res = LOW; return; }
        rho = fmax2(-1.0, fmin2(1.0, rho));
        const double ll = GeoPGZIPLogJoint(rho, r, t, mi, mj, mup,
                                            nugget1, nugget2, shape);
        if (!R_FINITE(ll)) { *res = LOW; return; }
        GeoKahanAdd(&total, weight * ll); ++used;
    }
    *res = (used > 0 && R_FINITE(total.sum)) ? total.sum : LOW;
}
