#include "header.h"
/******************************************************************************************/
/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
/******************************************************************************************/
void Comp_Cond_Gauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                        double *par, int *weigthed, double *res, double *mean1, double *mean2,
                        double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    const double sill = nuis[1];
    const double nugget = nuis[0];
    if(sill <= 0 || nugget < 0 || nugget > 1) {  // Changed sill check to <= 0
        *res = LOW;
        return;
    }
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double scale = 1.0 - nugget;
    const double sqrt_sill = sqrt(sill);
    double total = 0.0;  // Accumulator variable

    // Main processing loop
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Compute correlation
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            const double log_biv = log_biv_Norm(scale * corr, d1, d2, mean1[i], mean2[i], sill, 0);
            const double log_cond = dnorm(d2, mean2[i], sqrt_sill, 1);
            total += (log_biv - log_cond) * weights;
        }
    }
    // Final assignment with check
    *res = R_FINITE(total) ? total : LOW;
}

/******************************************************************************************/
void Comp_Cond_WrapGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                           double *par, int *weigthed, double *res, double *mean1, double *mean2,
                           double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    const double alfa = 2.0;
    const double nugget = nuis[0];
    const double sill = nuis[1];
    if(sill <= 0 || nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double scale = 1.0 - nugget;
    double total = 0.0;
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            const double biv = biv_wrapped(alfa, d1, d2, mean1[i], mean2[i], nugget, sill, scale * corr);
            const double cond = one_log_wrapped(alfa, d2, mean2[i], sill);
            
            total += (log(biv) - cond) * weights;
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}

/*********************************************************/
void Comp_Cond_Tukeyh2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                         double *par, int *weigthed, double *res, double *mean1, double *mean2,
                         double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    const double sill = nuis[1];
    const double nugget = nuis[0];
    const double tail = nuis[2];
    if(sill <= 0 || tail <= 0 || tail >= 0.5 || nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double scale = 1.0 - nugget;
    double total = 0.0;
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            const double biv = biv_tukey_h(scale * corr, d1, d2, mean1[i], mean2[i], tail, sill);
            const double cond = one_log_tukeyh(d2, mean2[i], sill, tail);
            total += (log(biv) - cond) * weights;
        }
    }
    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Cond_Tukeyhh2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                          double *par, int *weigthed, double *res, double *mean1, double *mean2,
                          double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    const double sill = nuis[1];
    const double nugget = nuis[0];
    const double h1 = nuis[3];
    const double h2 = nuis[2];
    if(sill <= 0 || h1 <= 0 || h1 >= 0.5 || h2 <= 0 || h2 >= 0.5 || nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double scale = 1.0 - nugget;
    double total = 0.0;
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            const double biv = biv_tukey_hh(scale * corr, d1, d2, mean1[i], mean2[i], sill, h1, h2);
            const double cond = one_log_tukeyhh(d2, mean2[i], sill, h1, h2);
            total += (log(biv) - cond) * weights;
        }
    }
    *res = R_FINITE(total) ? total : LOW;
}

/******************************************************************************************/
void Comp_Cond_SkewGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                            double *par, int *weigthed, double *res, double *mean1, double *mean2,
                            double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Parameter validation
    const double nugget = nuis[0];
    const double sill = nuis[1];
    const double skew = nuis[2];
    if(nugget < 0 || nugget >= 1 || sill <= 0) {
        *res = LOW;
        return;
    }
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    double total = 0.0;
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Compute correlation
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            const double biv_dens = biv_skew(corr, d1, d2, mean1[i], mean2[i], sill, skew, nugget);
            const double cond_dens = one_log_SkewGauss(d2, mean2[i], sill, skew);
            total += (log(biv_dens) - cond_dens) * weights;
        }
    }
    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Cond_T2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                    double *par, int *weigthed, double *res, double *mean1, double *mean2,
                    double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    const double sill = nuis[2];
    const double nugget = nuis[1];
    const double df = nuis[0];
    const double df1 = 1.0/nuis[0];
    if(sill <= 0 || nugget < 0 || nugget >= 1 || df <= 0 || df > 0.5) {
        *res = LOW;
        return;
    }

    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double inv_sqrt_sill = 1.0/sqrt(sill);
    const double inv_sill = 1.0/sill;
    double total = 0.0;

    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double qi = (d1 - mean1[i]) * inv_sqrt_sill;
            const double qj = (d2 - mean2[i]) * inv_sqrt_sill;
            const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
            const double weights = weighted ? CorFunBohman(lags[i], max_dist) : 1.0;
            const double l2 = one_log_T(d2, mean2[i], sill, df1);
            const double biv = biv_T(corr, qi, qj, df, nugget) * inv_sill;
            total += (log(biv) - l2) * weights;
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Cond_Gauss_misp_T2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                               double *par, int *weigthed, double *res, double *mean1, double *mean2,
                               double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    const double sill = nuis[2];
    const double nugget = nuis[1];
    const double df_param = nuis[0];
    
    if(sill <= 0 || nugget < 0 || nugget >= 1 || df_param <= 0 || df_param > 0.5) {
        *res = LOW;
        return;
    }

    const double df = 1.0 / df_param;
    const double df_ratio = df / (df - 2);
    const double sill_scaled = sill * df_ratio;
    const double sqrt_sill_scaled = sqrt(sill_scaled);
    const double scale = 1.0 - nugget;
    
    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    
    const double log_df_minus_2 = log(df - 2);
    const double two_lgamma_df_minus_1_half = 2 * lgammafn(0.5 * (df - 1));
    const double two_lgamma_df_half = 2 * lgammafn(0.5 * df);
    const double log_const_part = log_df_minus_2 + two_lgamma_df_minus_1_half - (log(2) + two_lgamma_df_half);
    double total = 0.0;
    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double lag = lags[i];
            double corr = CorFct(cormod, lag, 0, par, 0, 0);
            if(fabs(corr) > 0) {
                const double corr_sq = corr * corr;
                corr = exp(log_const_part + 
                          log(hypergeo(0.5, 0.5, 0.5 * df, corr_sq)) + 
                          log(corr * scale));
            }
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            const double l2 = dnorm(d2, mean2[i], sqrt_sill_scaled, 1);
            const double biv = log_biv_Norm(corr, d1, d2, mean1[i], mean2[i], sill_scaled, 0);
            total += (biv - l2) * weights;
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}

/*********************************************************/
void Comp_Cond_Gauss_misp_SkewT2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0;
    double weights=1.0,sill,nugget,skew,corr,corr2,df,bl,l2;


    df=1/nuis[0];
    nugget=nuis[1];
    sill=nuis[2];
    skew=nuis[3];

    if(df<2||fabs(skew)>1||sill<0||nugget<0||nugget>=1){*res=LOW; return;}
    //auxuliary variables
    double D1=(df-1)/2;
    double D2=df/2;
    //double delta=skew/sqrt(1-skew*skew);
    double MM=(sqrt(df)*skew)/(sqrt(M_PI))*exp(lgammafn(D1)-lgammafn(D2));
    double FF=(df/(df-2)-MM*MM);

     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                     corr=CorFct(cormod,lags[i],0,par,0,0)*(1-nugget);
                     corr2=corr_skewt(corr,df,skew);
                     if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          bl=log_biv_Norm(corr2,data1[i],data2[i],mean1[i]+sqrt(sill)*MM,
                                                                  mean2[i]+sqrt(sill)*MM,
                                                                  sill*FF,0);
                         // l1=dnorm(data1[i],mean1[i]+sqrt(sill)*MM,sqrt(sill*FF),1);
                          l2=dnorm(data2[i],mean2[i]+sqrt(sill)*MM,sqrt(sill*FF),1);
                        //*res+= (2*bl-l1-l2)*weights;
                        *res+= (bl-l2)*weights;


                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
/*********************************************************/
void Comp_Cond_Gauss_misp_Tukeygh2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i;double bl,corr,corr2,zi,zj,weights=1.0,eta,tail,sill,nugget,u,eta2,mu,vv,l2;
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
                 //   if(corr2<0) Rprintf("%f %f %f \n",corr2,par[0],par[1]);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                bl=log_biv_Norm(corr2,zi,zj,mean1[i]+sqrt(sill)*mu,
                                            mean2[i]+sqrt(sill)*mu, sill*vv,0);

                     // l1= dnorm(zi, mean1[i]+sqrt(sill)*mu,sqrt(sill*vv),1);
                      l2= dnorm(zj, mean2[i]+sqrt(sill)*mu,sqrt(sill*vv),1);
                     // *res+= (2*bl-l1-l2)*weights;
                      *res+= (bl-l2)*weights;

                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Cond_SinhGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                            double *par, int *weigthed, double *res, double *mean1, double *mean2,
                            double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    const double nugget = nuis[0];
    const double sill = nuis[1];
    const double param2 = nuis[2];
    const double param3 = nuis[3];
    if(nugget < 0 || nugget >= 1 || sill <= 0 || param3 < 0) {
        *res = LOW;
        return;
    }

    const int weighted = *weigthed;
    const int n_pairs = npairs[0];
    const double max_dist = maxdist[0];
    const double scale = 1.0 - nugget;
    const double inv_sqrt_sill = 1.0/sqrt(sill);
    const double inv_sill = 1.0/sill;
    double total = 0.0;

    for(int i = 0; i < n_pairs; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            const double weights = weighted ? CorFunBohman(lag, max_dist) : 1.0;
            const double std_d1 = (d1 - mean1[i]) * inv_sqrt_sill;
            const double std_d2 = (d2 - mean2[i]) * inv_sqrt_sill;
            const double biv = biv_sinh(scale * corr, std_d1, std_d2, 0, 0, param2, param3, 1) * inv_sill;
            const double cond = one_log_sas(d2, mean2[i], param2, param3, sill);
            total += weights * (log(biv) - cond);
        }
    }

    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Cond_Gamma2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                        double *par, int *weigthed, double *res, double *mean1, double *mean2,
                        double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    const double nugget = nuis[0];
    const double nuis2 = nuis[2];
    const int weighted_flag = *weigthed;
    const int npairs_val = npairs[0];
    const double maxdist_val = maxdist[0];

    if(nugget < 0 || nugget >= 1 || nuis2 < 0) {
        *res = LOW;
        return;
    }
    double total = 0.0;  // Accumulatore locale
    // Pre-calcola (1-nugget) per riutilizzarlo
    const double one_minus_nugget = 1.0 - nugget;
    for(int i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double zi = d1;
            const double zj = d2;
            const double m1 = mean1[i];
            const double m2 = mean2[i];
            const double corr = CorFct(cormod, lags[i], 0, par, 0, 0);
            const double scaled_corr = one_minus_nugget * corr;
            const double l2 = one_log_gamma(zj, m2, nuis2);
            const double bl = log(biv_gamma(scaled_corr, zi, zj, m1, m2, nuis2)) - l2;
            double weights = 1.0;
            if(weighted_flag) {
                weights = CorFunBohman(lags[i], maxdist_val);
            }
            total += weights * bl;
        }
    }
    // Assegna il risultato finale
    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Cond_Weibull2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                          double *par, int *weigthed, double *res, double *mean1, double *mean2,
                          double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    const double nugget = nuis[0], nuis2 = nuis[2], maxdist_val = maxdist[0];
    const int weighted_flag = *weigthed, npairs_val = npairs[0];
    if(nugget < 0 || nugget >= 1 || nuis2 < 0) {
        *res = LOW;
        return;
    }
    double total = 0.0, one_minus_nugget = 1.0 - nugget;
    for(int i = 0; i < npairs_val; i++) {
        const double d1 = data1[i], d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double zi = d1, zj = d2, m1 = mean1[i], m2 = mean2[i], lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            const double scaled_corr = one_minus_nugget * corr;
            const double l2 = one_log_weibull(zj, m2, nuis2);
            const double bl = log(biv_Weibull(scaled_corr, zi, zj, m1, m2, nuis2)) - l2;
            double weights = weighted_flag ? CorFunBohman(lag, maxdist_val) : 1.0;
            total += weights * bl;
        }
    }
    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Cond_LogGauss2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                           double *par, int *weigthed, double *res, double *mean1, double *mean2,
                           double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    const double nugget = nuis[0], sill = nuis[1], maxdist_val = maxdist[0];
    const int weighted_flag = *weigthed, npairs_val = npairs[0];
    if(sill < 0 || nugget < 0 || nugget > 1) {
        *res = LOW;
        return;
    }
    double total = 0.0, one_minus_nugget = 1.0 - nugget;
    for(int i = 0; i < npairs_val; i++) {
        const double d1 = data1[i], d2 = data2[i];
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double zi = d1, zj = d2, m1 = mean1[i], m2 = mean2[i], lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0);
            const double scaled_corr = one_minus_nugget * corr;
            const double l2 = one_log_loggaussian(zj, m2, sill);
            const double bl = log(d2lognorm(zi, zj, sill, nugget, m1, m2, scaled_corr)) - l2;
            double weights = weighted_flag ? CorFunBohman(lag, maxdist_val) : 1.0;
            total += weights * bl;
        }
    }
    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Cond_Beta2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i;double corr,zi,zj,weights=1.0,bl,l2=0.0;
    double nugget=nuis[0];
    double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max)  {*res=LOW;  return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=data1[i]; zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                     if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                  //  l1=one_log_beta(zi,nuis[2],nuis[3],min,max);
                    l2=one_log_beta(zj,nuis[2],nuis[3],min,max);
                //  bl=2*log(biv_beta((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-(l1+l2);
                bl=log(biv_beta((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-l2;
        *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Kumaraswamy2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i;double corr,zi,zj,weights=1.0,bl,l2=0.0;
    double nugget=nuis[0];
         double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
               zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);

                    //l1=one_log_kumma(zi,mean1[i],nuis[2],nuis[3],min,max);
                    l2=one_log_kumma(zj,mean2[i],nuis[2],nuis[3],min,max);

                     if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                  //bl=2*log(biv_Kumara((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-(l1+l2);
                  bl=log(biv_Kumara((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-l2;
        *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Kumaraswamy22mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i;double corr,zi,zj,weights=1.0,bl,l2=0.0;
    double nugget=nuis[0];
         double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                   // l1=one_log_kumma2(zi,mean1[i],nuis[2],nuis[3],min,max);
                    l2=one_log_kumma2(zj,mean2[i],nuis[2],nuis[3],min,max);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
    //bl=2*log(biv_Kumara2((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-(l1+l2);
    bl=log(biv_Kumara2((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-l2;
        *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Gauss_misp_Pois2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                  double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                  double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    const int N = 2;
    const double nugget = nuis[0], maxdist_val = maxdist[0];
    const int weighted_flag = *weigthed, npairs_val = npairs[0];
    
    if(nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }
    // Allocate memory for covariance matrix and data vector
    double **M = (double **) R_Calloc(N, double *);
    for(int i = 0; i < N; i++) {
        M[i] = (double *) R_Calloc(N, double);
    }
    double *dat = (double *) R_Calloc(N, double);
    double total = 0.0;
    for(int i = 0; i < npairs_val; i++) {
        if(!ISNAN(data1[i]) && !ISNAN(data2[i])) {
            const double mui = exp(mean1[i]);
            const double muj = exp(mean2[i]);
            const double lag = lags[i];
            const double corr = CorFct(cormod, lag, 0, par, 0, 0) * (1 - nugget);
            const double corr1 = corr_pois(corr, mui, muj);
            M[0][0] = mui;
            M[1][1] = muj;
            M[0][1] = M[1][0] = sqrt(mui * muj) * corr1;
            
            // Set up data vector
            dat[0] = data1[i] - mui;
            dat[1] = data2[i] - muj;
            const double l2 = dnorm(data2[i], muj, sqrt(muj), 1);
            const double bl = log(dNnorm(N, M, dat)) - l2;
            const double weights = weighted_flag ? CorFunBohman(lag, maxdist_val) : 1.0;
            total += bl * weights;
        }
    }
    for(int i = 0; i < N; i++) {
        R_Free(M[i]);
    }
    R_Free(M);
    R_Free(dat);
    
    *res = R_FINITE(total) ? total : LOW;
}
/*********************************************************/
void Comp_Cond_BinomNNGauss_misp2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0, N=2,n1,n2;
    double u,v,m1,m2,l2,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0;
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
                 m1=n1*p1;m2=n2*p2;
                 if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                 M[0][0]=m1*(1-p1);   M[1][1]=m2*(1-p2);  // var1 var2
                 M[0][1]= fmin_int(n1,n2)*(p11-p1*p2) ;       // covariance
                 M[1][0]= M[0][1];
                 dat[0]=u-m1;dat[1]=v-m2; 
                 //Rprintf("%d %f %f %f \n",fmin_int(n1,n2),p1,p2,p11 );
                 //l1=dnorm(u,m1,sqrt(m1*(1-p1)),1);
                 l2=dnorm(v,m2,sqrt(m2*(1-p2)),1);;
                // bl= 2*log(dNnorm(N,M,dat)) -(l1+l2); 
                 bl= log(dNnorm(N,M,dat)) -l2; 
                 *res+= bl*weights;       
                }}
    for(i=0;i<N;i++)  {R_Free(M[i]);}
    R_Free(M);
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Cond_Gauss_misp_PoisGamma2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0,N=2;
    double  weights=1.0,corr,corr1,mui,muj,bl,l2,bi,bj,vvi,vvj;
    double nugget=nuis[0];
      if(nugget<0||nugget>=1){*res=LOW; return;}
double **M;
        M= (double **) R_Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) R_Calloc(N,double);}
    double *dat;
    dat=(double *) R_Calloc(N,double);
    for(i=0;i<npairs[0];i++){

                  if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
             //***********/
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                    bi= nuis[2]/mui; bj= nuis[2]/muj;
                    vvi= mui*(1+1/bi); vvj= muj*(1+1/bj);
                     corr=CorFct(cormod,lags[i],0,par,0,0)*(1-nugget);
                      corr1=corr_pois_gen(corr,mui, muj, nuis[2]);
                      if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                        M[0][0]=vvi; M[1][1]=vvj;M[0][1]=sqrt(vvi*vvj)*corr1;M[1][0]= M[0][1];
                        dat[0]=data1[i]-mui;dat[1]=data2[i]-muj;
                     // l1=dnorm(data1[i],mui,sqrt(vvi),1);
                      l2=dnorm(data2[i],muj,sqrt(vvj),1);;
                     // bl=2*log(dNnorm(N,M,dat))-(l1+l2);
                      bl=log(dNnorm(N,M,dat))-l2;

                      *res+= bl*weights;
                    }}
   for(i=0;i<N;i++)  {R_Free(M[i]);}
    R_Free(M);

    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
void Comp_Cond_PoisGamma2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0, uu,ww;
    double weights=1.0,corr,mui,muj,bl,l2;
    double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}
  // Rprintf("%d   \n",npairs[0]);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],0,par,0,0);
                    // if(fabs(corr)>1|| !R_FINITE(corr)) {*res=LOW; return;}
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                      uu=(int) data1[i];  ww=(int) data2[i];
                    //  l1=one_log_dpoisgamma(uu,mui,nuis[2]);
                      l2=one_log_dpoisgamma(ww,muj,nuis[2]);
         //bl=2*log(biv_PoissonGamma((1-nugget)*corr,uu,ww,mui, muj,nuis[2]))  - (l1+l2);
         bl=log(biv_PoissonGamma((1-nugget)*corr,uu,ww,mui, muj,nuis[2]))  - l2;
                      *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Cond_PoisGammaZIP2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{

    int i=0, uu,vv;
    double weights=1.0,corr,mui,muj,bl,l2=0.0,u,v;
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
                    l2=one_log_PoisgammaZIP(vv,muj,mup,shape);
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                      bl=log(biv_PoissonGammaZIP(corr,uu,vv,mui, muj,mup,nugget1,nugget2,shape))-l2;
                      *res+= bl*weights;
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
void Comp_Cond_Pois2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0, uu,ww;
    double weights=1.0,corr,mui,muj,bl,l2;
    double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}
  // Rprintf("%d   \n",npairs[0]);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],0,par,0,0);
                    // if(fabs(corr)>1|| !R_FINITE(corr)) {*res=LOW; return;}
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                      uu=(int) data1[i];  ww=(int) data2[i];
                     // l1=dpois(uu,mui,1);
                      l2=dpois(ww,muj,1);
                     // bl=2*log(biv_Poisson((1-nugget)*corr,uu,ww,mui, muj)) - (l1+l2);
                      bl=log(biv_Poisson((1-nugget)*corr,uu,ww,mui, muj)) - l2;
                      *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Cond_BinomGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0, uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0;
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
                          l2=dbinom(vv,N1[0],p2,1);
                        bl=log(biv_binom (N1[0],uu,vv,p1,p2,p11))-l2;
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Cond_BinomLogi2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0, uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);
    p11=pblogi22(log(exp(ai)+1),log(exp(aj)+1),(1-nugget)*corr);
    //Rprintf("p11: %f\n",p11);
                p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                u=data1[i];v=data2[i];
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u; vv=(int) v;
                          l2=dbinom(vv,N1[0],p2,1);
                        bl=log(biv_binom (N1[0],uu,vv,p1,p2,p11))-l2;
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
void Comp_Cond_BinomNNGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0, uu=0,vv=0,n1,n2;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0;
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
                 n1=N1[i];n2=N2[i];
                 if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                 uu=(int) u; vv=(int) v;
                 l2=dbinom(vv,n2,p2,1);
                 bl=log(biv_binom222(n1,n2,uu,vv,p1,p2,p11))-l2;
                 *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Cond_BinomNNLogi2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0, uu=0,vv=0,n1,n2;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 n1=N1[i];n2=N2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);
                 p11=pblogi22(log(exp(ai)+1),log(exp(aj)+1),(1-nugget)*corr);
                 p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                 u=data1[i];v=data2[i];
                 if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                 uu=(int) u; vv=(int) v;
                 l2=dbinom(vv,n2,p2,1);
                 bl=log(biv_binom222(n1,n2,uu,vv,p1,p2,p11))-l2;
                 *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_BinomnegGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
       if( nugget>=1 || nugget<0){*res=LOW; return;}
    //compute the composite log-likelihood:

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);
                    p11=pbnorm22(ai,aj,(1-nugget)*corr);
                p1=pnorm(ai,0,1,1,0);p2=pnorm(aj,0,1,1,0);
                    u=data1[i];v=data2[i];
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u;vv=(int) v;
                         l2=one_log_negbinom_marg(vv,N1[0],p2);
                         bl=log(biv_binomneg(N1[0],uu,vv,p1,p2,p11))-l2;
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/*********************************************************/
void Comp_Cond_BinomnegLogi2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
       if( nugget>=1 || nugget<0){*res=LOW; return;}
    //compute the composite log-likelihood:

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);
                    p11=pblogi22(log(exp(ai)+1),log(exp(aj)+1),(1-nugget)*corr);
                p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                    u=data1[i];v=data2[i];
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u;vv=(int) v;
                         l2=one_log_negbinom_marg(vv,N1[0],p2);
                        bl=log(biv_binomneg(N1[0],uu,vv,p1,p2,p11))-l2;
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}


/******************************************************/
void Comp_Cond_BinomnegBinary2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0;
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
                l2=dbinom(vv,1,1-pow(p2,N1[0]),1);
                bl=log(biv_binegbinary(N1[0],uu,vv,p1,p2,p11))-l2;

            *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}


void Comp_Cond_TWOPIECETukeyh2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,eta,tail,qq,sill,nugget,l2=0.0;
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

            //l1=one_log_two_pieceTukey(zi,mean1[i],sill,tail,eta);
            l2=one_log_two_pieceTukey(zj,mean2[i],sill,tail,eta);

           p11=pbnorm22(qq,qq,corr);
           if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
          // bl=2*log(biv_two_pieceTukeyh((1-nugget)*corr,zi,zj,sill,eta,tail,p11,mean1[i],mean2[i]))-(l1+l2);
           bl=log(biv_two_pieceTukeyh((1-nugget)*corr,zi,zj,sill,eta,tail,p11,mean1[i],mean2[i]))-l2;
               *res+= weights*bl;
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_TWOPIECET2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,qq,l2=0.0;
    double eta=nuis[3];  //skewness parameter
    double sill=nuis[2];
    double nugget=nuis[1];
    double df=nuis[0];
    if(sill<0||nugget<0||nugget>=1 ||fabs(eta)>1|| df >0.5||df<0) {*res=LOW;  return;}
       qq=qnorm((1-eta)/2,0,1,1,0);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                corr=CorFct(cormod,lags[i],0,par,0,0);
                //l1=one_log_two_pieceT(zi,mean1[i],sill,df,eta);
                l2=one_log_two_pieceT(zj,mean2[i],sill,df,eta);
                 p11=pbnorm22(qq,qq,corr);
                 if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                 /********************************************************/
                // bl=2*log(biv_two_pieceT(corr,zi,zj,sill,df,eta,p11,mean1[i],mean2[i],nugget)) -(l1+l2);
                 bl=log(biv_two_pieceT(corr,zi,zj,sill,df,eta,p11,mean1[i],mean2[i],nugget)) -l2;
                 /********************************************************/
                         *res+= weights*bl;
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_TWOPIECEGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,eta,qq,sill,nugget,l2=0.0;
    eta=nuis[2];  //skewness parameter
    sill=nuis[1];
    nugget=nuis[0];
       qq=qnorm((1-eta)/2,0,1,1,0);
         if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1) {*res=LOW;  return;}
   for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);
               // l1=one_log_two_pieceGauss(zi,mean1[i],sill,eta);
                l2=one_log_two_pieceGauss(zj,mean2[i],sill,eta);

                p11=pbnorm22(qq,qq,corr);
                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                   // bl=2*log(biv_two_pieceGaussian((1-nugget)*corr,zi,zj,sill,eta,p11,mean1[i],mean2[i]))-l1-l2;
                    bl=log(biv_two_pieceGaussian((1-nugget)*corr,zi,zj,sill,eta,p11,mean1[i],mean2[i]))-l2;
                 
                    *res+= weights*bl;
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/********************************************************/
void Comp_Cond_TWOPIECEBIMODAL2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,eta,qq,sill,df,nugget,delta ,l2=0.0;

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
             //   l1=one_log_bomidal(zi,mean1[i],sill,df,delta,eta);
                l2=one_log_bomidal(zj,mean2[i],sill,df,delta,eta);
                p11=pbnorm22(qq,qq,corr);
                 if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    /********************************************************/
                 // bl=2*log(biv_two_piece_bimodal((1-nugget)*corr,zi,zj,sill,df,delta,eta,p11,mean1[i],mean2[i]))-(l1+l2);
                  bl=log(biv_two_piece_bimodal((1-nugget)*corr,zi,zj,sill,df,delta,eta,p11,mean1[i],mean2[i]))-l2;
                           *res+= weights*bl;
                }}}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/********************************************************/
void Comp_Cond_BinomnegGaussZINB2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0;
    double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];
      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);
                    u=data1[i];v=data2[i];
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u;
                         vv=(int) v;
                    l2=one_log_BinomnegZIP(vv,N1[0],aj,mup);
                    bl=log(biv_binomnegZINB(N1[0],corr,uu,vv,ai,aj,nugget1,nugget2,mup))-l2;
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;

}
/************************************************/
void Comp_Cond_PoisZIP2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{

    int i=0, uu,vv;
    double weights=1.0,corr,mui,muj,bl,l2=0.0,u,v;
   double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];


      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],0,par,0,0);
                        u=data1[i];v=data2[i];
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u;vv=(int) v;

                 //   l1=one_log_PoisZIP(uu,mui,mup);
                    l2=one_log_PoisZIP(vv,muj,mup);
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                   //  bl=2*log(biv_PoissonZIP(corr,uu,vv,mui, muj,mup,nugget1,nugget2))-(l1+l2);
                      bl=log(biv_PoissonZIP(corr,uu,vv,mui, muj,mup,nugget1,nugget2))-l2;
                      *res+= bl*weights;
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Cond_Gauss_misp_PoisZIP2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0;
    double weights=1.0,corr,mui,muj,bl ,l2=0.0;
    double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];
    double p=pnorm(mup,0,1,1,0);

      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}

  // Rprintf("%d   \n",npairs[0]);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],0,par,0,0);


                      //l1=dnorm(data1[i],(1-p)*mui,sqrt(mui*(1-p)*(1+p*mui)),1);
                      l2=dnorm(data2[i],(1-p)*muj,sqrt(muj*(1-p)*(1+p*muj)),1);


                      if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                      //bl=2*log(biv_Mis_PoissonZIP(corr,data1[i],data2[i],mui, muj,mup,nugget1,nugget2))-(l1+l2);
                      bl=log(biv_Mis_PoissonZIP(corr,data1[i],data2[i],mui, muj,mup,nugget1,nugget2))-l2;

                      *res+= bl*weights;
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Cond_LogLogistic2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{

    int i;double corr,zi,zj,weights=1.0,bl=1.0,l2=0.0;
    double nugget=nuis[0];
    if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);

                   // l1=one_log_loglogistic(zi,exp(mean1[i]),nuis[2]);
                    l2=one_log_loglogistic(zj,exp(mean2[i]),nuis[2]);

                  //  bl=2*log(biv_LogLogistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]))-(l1+l2);
                  bl=log(biv_LogLogistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]))-l2;
                           
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);

  *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}



void Comp_Cond_Logistic2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{

    int i;double corr,zi,zj,weights=1.0,bl=1.0,l2=0.0;
    double nugget=nuis[0];
    if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);

                    //l1=one_log_logistic(zi,mean1[i],nuis[1]) ;
                    l2=one_log_logistic(zj,mean2[i],nuis[1])  ;

                   // bl=2*log(biv_Logistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[1])) -(l1+l2);
                    bl=log(biv_Logistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[1])) -l2;
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);

  *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/********************* SPACE TIME CASE *****************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
void Comp_Cond_Gauss_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                           double *par, int *weigthed, double *res, double *mean1, double *mean2,
                           double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    const double nugget = nuis[0];
    const double sill = nuis[1];
    const int npairs_val = npairs[0];
    
    if(sill < 0 || nugget < 0 || nugget >= 1) {
        *res = LOW; 
        return;
    }
    const double sill_sqrt = sqrt(sill);
    const double sill_factor = 1 - nugget;
    
    double total = 0.0;
    int i;
    for(i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double m1 = mean1[i];
            const double m2 = mean2[i];
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0);
            const double bl = log_biv_Norm(sill_factor * corr, d1, d2, m1, m2, sill, 0);
            const double l2 = dnorm(d2, m2, sill_sqrt, 1);
            
            total += (bl - l2);  // weights=1.0 quindi non serve moltiplicare
        }
    }
    *res = total;
    if(!R_FINITE(*res)) {
        *res = LOW;
    }
}
/******************************************************************************************/
void Comp_Cond_T_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                       double *par, int *weigthed, double *res, double *mean1, double *mean2,
                       double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Estrai i parametri una sola volta
    const double df = nuis[0];
    const double nugget = nuis[1];
    const double sill = nuis[2];
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];
    
    // Controllo rapido dei parametri
    if(sill < 0 || nugget < 0 || nugget >= 1 || df < 0 || df > 0.5) {
        *res = LOW;
        return;
    }
    
    // Pre-calcola valori utili
    const double sill_sqrt = sqrt(sill);
    const double sill_inv = 1.0 / sill;
    
    double total = 0.0;
    int i;
    
    for(i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double m1 = mean1[i];
            const double m2 = mean2[i];
            const double z1 = (d1 - m1) / sill_sqrt;
            const double z2 = (d2 - m2) / sill_sqrt;
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0);
            double weights = 1.0;
            if(weigthed_val) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                          CorFunBohman(lagt[i], maxtime_val);
            }
            const double biv_t = biv_T(corr, z1, z2, df, nugget);
            const double bl = log(biv_t * sill_inv);
            const double l2 = one_log_T(d2, m2, sill, 1/df);
            
            total += (bl - l2) * weights;
        }
    }
    *res = total;
    if(!R_FINITE(*res)) {
        *res = LOW;
    }
}
/******************************************************************************************/
void Comp_Cond_Gauss_misp_T_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                  double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                  double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Estrai i parametri una sola volta
    const double df_param = nuis[0];
    const double nugget = nuis[1];
    const double sill = nuis[2];
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];
    
    // Controllo rapido dei parametri
    if(sill < 0 || nugget < 0 || nugget >= 1 || df_param < 0 || df_param > 0.5) {
        *res = LOW;
        return;
    }
    const double df = 1.0 / df_param;
    const double var = sill * df / (df - 2.0);
    const double var_sqrt = sqrt(var);
    const double log_df_minus_2 = log(df - 2.0);
    const double two_lgamma_df_minus_1_half = 2.0 * lgammafn(0.5 * (df - 1.0));
    const double two_lgamma_df_half = 2.0 * lgammafn(0.5 * df);
    const double log_2 = log(2.0);
    double total = 0.0;
    int i;
    for(i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double m1 = mean1[i];
            const double m2 = mean2[i];
            double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0);
            const double corr_sq = corr * corr;

            corr = exp(log_df_minus_2 + two_lgamma_df_minus_1_half - (log_2 + two_lgamma_df_half) + 
                      log(hypergeo(0.5, 0.5, 0.5 * df, corr_sq)) + log(corr * (1.0 - nugget)));

            double weights = 1.0;
            if(weigthed_val) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                          CorFunBohman(lagt[i], maxtime_val);
            }
            const double bl = log_biv_Norm(corr, d1, d2, m1, m2, var, 0);
            const double l2 = dnorm(d2, m2, var_sqrt, 1);
            
            total += (bl - l2) * weights;
        }
    }
    
    *res = total;
    if(!R_FINITE(*res)) {
        *res = LOW;
    }
}
/*********************************************************/
void Comp_Cond_Tukeyh_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                            double *par, int *weigthed, double *res, double *mean1, double *mean2,
                            double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Estrai i parametri una sola volta
    const double nugget = nuis[0];
    const double sill = nuis[1];
    const double tail = nuis[2];
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];
    const double nugget_factor = 1.0 - nugget;

    if(sill < 0 || tail < 0 || tail > 0.5 || nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }
    double total = 0.0;
    int i;
    for(i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double m1 = mean1[i];
            const double m2 = mean2[i];
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0) * nugget_factor;
            double weights = 1.0;
            if(weigthed_val) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                          CorFunBohman(lagt[i], maxtime_val);
            }
            const double l2 = one_log_tukeyh(d2, m2, sill, tail);
            const double biv = biv_tukey_h(corr, d1, d2, m1, m2, tail, sill);
            const double bl = log(biv) - l2;
            
            total += weights * bl;
        }
    }
    *res = total;
    if(!R_FINITE(*res)) {
        *res = LOW;
    }
}
/*********************************************************/
void Comp_Cond_Tukeyhh_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                             double *par, int *weigthed, double *res, double *mean1, double *mean2,
                             double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Estrai i parametri una sola volta
    const double nugget = nuis[0];
    const double sill = nuis[1];
    const double h2 = nuis[2];
    const double h1 = nuis[3];
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];
    const double nugget_factor = 1.0 - nugget;
    if(sill < 0 || h1 < 0 || h1 > 0.5 || h2 < 0 || h2 > 0.5 || nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }
    double total = 0.0;
    int i;
    for(i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double m1 = mean1[i];
            const double m2 = mean2[i];
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0) * nugget_factor;
            double weights = 1.0;
            if(weigthed_val) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                          CorFunBohman(lagt[i], maxtime_val);
            }
            const double l2 = one_log_tukeyhh(d2, m2, sill, h1, h2);
            const double biv = biv_tukey_hh(corr, d1, d2, m1, m2, sill, h1, h2);
            const double bl = log(biv) - l2;
            
            total += weights * bl;
        }
    }
    *res = total;
    if(!R_FINITE(*res)) {
        *res = LOW;
    }
}
/******************************************************************************************/
void Comp_Cond_SkewGauss_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                               double *par, int *weigthed, double *res, double *mean1, double *mean2,
                               double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Estrai parametri una volta sola
    const double nugget = nuis[0];
    const double sill = nuis[1];
    const double skew = nuis[2];
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];

    if(nugget < 0 || nugget >= 1 || sill < 0) {
        *res = LOW;
        return;
    }
    
    double total = 0.0;
    int i;
    for(i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double m1 = mean1[i];
            const double m2 = mean2[i];
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0);
            double weights = 1.0;
            if(weigthed_val) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                          CorFunBohman(lagt[i], maxtime_val);
            }
            const double l2 = one_log_SkewGauss(d2, m2, sill, skew);
            const double biv = biv_skew(corr, d1, d2, m1, m2, sill, skew, nugget);
            const double bb = log(biv) - l2;
            
            total += weights * bb;
        }
    }
    *res = total;
    if(!R_FINITE(*res)) {
        *res = LOW;
    }
}
/******************************************************************************************/
void Comp_Cond_SinhGauss_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                               double *par, int *weigthed, double *res, double *mean1, double *mean2,
                               double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // Estrai parametri una volta sola
    const double nugget = nuis[0];
    const double sill = nuis[1];
    const double param1 = nuis[2];  // Parametro 1 della distribuzione Sinh
    const double param2 = nuis[3];  // Parametro 2 della distribuzione Sinh
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];
    const double nugget_factor = 1.0 - nugget;
    
    // Controllo parametri
    if(param2 < 0 || sill < 0 || nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }
    double total = 0.0;
    int i;
    for(i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double m1 = mean1[i];
            const double m2 = mean2[i];
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0) * nugget_factor;
            double weights = 1.0;
            if(weigthed_val) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                          CorFunBohman(lagt[i], maxtime_val);
            }
            const double l2 = one_log_sas(d2, m2, param1, param2, sill);
            const double biv = biv_sinh(corr, d1, d2, m1, m2, param1, param2, sill);
            const double bb = log(biv) - l2;
            
            total += weights * bb;
        }
    }
    
    *res = total;
    if(!R_FINITE(*res)) {
        *res = LOW;
    }
}
/*********************************************************/
void Comp_Cond_Gamma_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                           double *par, int *weigthed, double *res, double *mean1, double *mean2,
                           double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // ### 1. Estrazione e pre-calcolo parametri ###
    const double nugget = nuis[0];
    const double shape = nuis[2];  // Parametro shape della Gamma
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];
    const double nugget_factor = 1.0 - nugget;
    if (nugget < 0 || nugget >= 1 || shape < 0) {
        *res = LOW;
        return;
    }
    double total = 0.0;
    for (int i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];

        if (!ISNAN(d1) && !ISNAN(d2)) {
            const double m1 = mean1[i];
            const double m2 = mean2[i];
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0) * nugget_factor;
            double weights = 1.0;
            if (weigthed_val) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                          CorFunBohman(lagt[i], maxtime_val);
            }
            const double l2 = one_log_gamma(d2, m2, shape);
            const double biv = biv_gamma(corr, d1, d2, m1, m2, shape);
            const double bl = log(biv) - l2;

            total += weights * bl;
        }
    }
    *res = total;
    if (!R_FINITE(*res)) {
        *res = LOW;
    }
}
/*********************************************************/
void Comp_Cond_Weibull_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                             double *par, int *weigthed, double *res, double *mean1, double *mean2,
                             double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // ### 1. Estrazione parametri e pre-calcolo ###
    const double nugget = nuis[0];
    const double shape = nuis[2];  // Parametro shape della Weibull
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];
    const double nugget_factor = 1.0 - nugget;
    if (nugget < 0 || nugget >= 1 || shape < 0) {
        *res = LOW;
        return;
    }
    double total = 0.0;
    for (int i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];

        if (!ISNAN(d1) && !ISNAN(d2)) {
            const double m1 = mean1[i];
            const double m2 = mean2[i];
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0) * nugget_factor;
            double weights = 1.0;
            if (weigthed_val) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                          CorFunBohman(lagt[i], maxtime_val);
            }
            const double l2 = one_log_weibull(d2, m2, shape);
            const double biv = biv_Weibull(corr, d1, d2, m1, m2, shape);
            const double bl = log(biv) - l2;

            total += weights * bl;
        }
    }
    *res = total;
    if (!R_FINITE(*res)) {
        *res = LOW;
    }
}
/*********************************************************/
void Comp_Cond_LogGauss_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                              double *par, int *weigthed, double *res, double *mean1, double *mean2,
                              double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    const double nugget = nuis[0];
    const double sill = nuis[1];
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];
    const double nugget_factor = 1.0 - nugget;
    if(sill < 0 || nugget < 0 || nugget > 1) {
        *res = LOW; 
        return;
    }
    double total = 0.0;
    
    for(int i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double m1 = mean1[i];
            const double m2 = mean2[i];
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0) * nugget_factor;
            double weights = 1.0;
            if(weigthed_val) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                          CorFunBohman(lagt[i], maxtime_val);
            }
            const double l2 = one_log_loggaussian(d2, m2, sill);
            const double biv = d2lognorm(d1, d2, sill, nugget, m1, m2, corr);
            const double bl = log(biv) - l2;
            
            total += weights * bl;
        }
    }
    *res = total;
    if(!R_FINITE(*res)) {
        *res = LOW;
    }
}

/******************************************************************************************/
void Comp_Cond_WrapGauss_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                               double *par, int *weigthed, double *res, double *mean1, double *mean2,
                               double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    const double alfa = 2.0;  // Valore costante
    const double nugget = nuis[0];
    const double sill = nuis[1];
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];
    const double nugget_factor = 1.0 - nugget;
    if(sill < 0 || nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }
    double total = 0.0;
    for(int i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double m1 = mean1[i];
            const double m2 = mean2[i];
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0) * nugget_factor;
            double weights = 1.0;
            if(weigthed_val) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                          CorFunBohman(lagt[i], maxtime_val);
            }
            const double biv = biv_wrapped(alfa, d1, d2, m1, m2, nugget, sill, corr);
            const double l2 = one_log_wrapped(alfa, d2, m2, sill);
            const double bl = log(biv) - l2;

            total += weights * bl;
        }
    }
    *res = total;
    if(!R_FINITE(*res)) {
        *res = LOW;
    }
}

/*********************************************************/
void Comp_Cond_Beta_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{

    int i;double corr,zi,zj,weights=1.0,bl,l2=0.0;
    double nugget=nuis[0];
    double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max)  {*res=LOW;  return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=data1[i]; zj=data2[i];
                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                              if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                  //  l1=one_log_beta(zi,nuis[2],nuis[3],min,max);
                    l2=one_log_beta(zj,nuis[2],nuis[3],min,max);
                 // bl=2*log(biv_beta((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-(l1+l2);
                 bl=log(biv_beta((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-l2;

        *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/

void Comp_Cond_Kumaraswamy_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{

    int i;double corr,zi,zj,weights=1.0,bl,l2=0.0;
    double nugget=nuis[0];
         double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
               zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);

                    //l1=one_log_kumma(zi,mean1[i],nuis[2],nuis[3],min,max);
                    l2=one_log_kumma(zj,mean2[i],nuis[2],nuis[3],min,max);

                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                  //bl=2*log(biv_Kumara((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-(l1+l2);
                  bl=log(biv_Kumara((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-l2;

        *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Kumaraswamy2_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i;double corr,zi,zj,weights=1.0,bl,l2=0.0;
    double nugget=nuis[0];
         double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                  //  l1=one_log_kumma2(zi,mean1[i],nuis[2],nuis[3],min,max);
                    l2=one_log_kumma2(zj,mean2[i],nuis[2],nuis[3],min,max);
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                  //bl=2*log(biv_Kumara2((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-(l1+l2);
                  bl=log(biv_Kumara2((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-l2;
        *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}



/*********************************************************/
void Comp_Cond_Gauss_misp_Pois_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                    double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                    double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    const int N = 2;
    const double nugget = nuis[0];
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];
    const double nugget_factor = 1.0 - nugget;
    // Check parameters
    if(nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }
    // Allocate memory once
    double **M = (double **) R_Calloc(N, double *);
    double *dat = (double *) R_Calloc(N, double);
    for(int i = 0; i < N; i++) {
        M[i] = (double *) R_Calloc(N, double);
    }

    // =============== 2. Main Computation Loop ===============
    double total = 0.0;

    for(int i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];

        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Precompute means
            const double mui = exp(mean1[i]);
            const double muj = exp(mean2[i]);
            const double sqrt_muj = sqrt(muj);
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0) * nugget_factor;
            const double corr1 = corr_pois(corr, mui, muj);
            double weights = 1.0;
            if(weigthed_val) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                         CorFunBohman(lagt[i], maxtime_val);
            }
            M[0][0] = mui;
            M[1][1] = muj;
            M[0][1] = M[1][0] = sqrt(mui * muj) * corr1;
            dat[0] = d1 - mui;
            dat[1] = d2 - muj;
            const double l2 = dnorm(d2, muj, sqrt_muj, 1);
            const double bl = log(dNnorm(N, M, dat)) - l2;

            total += bl * weights;
        }
    }

    // =============== 3. Cleanup and Finalization ===============
    for(int i = 0; i < N; i++) {
        R_Free(M[i]);
    }
    R_Free(M);
    R_Free(dat);
    *res = total;
    if(!R_FINITE(*res)) {
        *res = LOW;
    }
}

void Comp_Cond_Gauss_misp_PoisGamma_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                         double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                         double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // ============= 1. INITIALIZATION AND PARAMETER SETUP =============
    const int N = 2;
    const double nugget = nuis[0];
    const double gamma_param = nuis[2];  // Gamma distribution parameter
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];
    const double nugget_factor = 1.0 - nugget;
    if(nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }
    double **M = (double **) R_Calloc(N, double *);
    double *dat = (double *) R_Calloc(N, double);
    for(int i = 0; i < N; i++) {
        M[i] = (double *) R_Calloc(N, double);
    }
    double total = 0.0;

    for(int i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];

        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Precompute means and variances
            const double mui = exp(mean1[i]);
            const double muj = exp(mean2[i]);
            const double bi = gamma_param / mui;
            const double bj = gamma_param / muj;
            const double vvi = mui * (1.0 + 1.0/bi);
            const double vvj = muj * (1.0 + 1.0/bj);
            const double sqrt_vvj = sqrt(vvj);
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0) * nugget_factor;
            const double corr1 = corr_pois_gen(corr, mui, muj, gamma_param);
            double weights = 1.0;
            if(weigthed_val) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                         CorFunBohman(lagt[i], maxtime_val);
            }
            M[0][0] = vvi;
            M[1][1] = vvj;
            M[0][1] = M[1][0] = sqrt(vvi * vvj) * corr1;
            dat[0] = d1 - mui;
            dat[1] = d2 - muj;
            const double l2 = dnorm(d2, muj, sqrt_vvj, 1);
            const double bl = log(dNnorm(N, M, dat)) - l2;

            total += bl * weights;
        }
    }
    // ============= 4. CLEANUP AND FINALIZATION =============
    for(int i = 0; i < N; i++) {
        R_Free(M[i]);
    }
    R_Free(M);
    R_Free(dat);
    *res = total;
    if(!R_FINITE(*res)) {
        *res = LOW;
    }
}

/*********************************************************/
void Comp_Cond_PoisGamma_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                               double *par, int *weigthed, double *res, double *mean1, double *mean2,
                               double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // ============= 1. INITIALIZATION AND PARAMETER SETUP =============
    const double nugget = nuis[0];
    const double gamma_param = nuis[2];  // Gamma distribution parameter
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];
    const double nugget_factor = 1.0 - nugget;

    // Parameter validation
    if(nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }
    double total = 0.0;

    for(int i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Precompute means and convert data to integers
            const double mui = exp(mean1[i]);
            const double muj = exp(mean2[i]);
            const int uu = (int)d1;
            const int ww = (int)d2;
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0) * nugget_factor;
            double weights = 1.0;
            if(weigthed_val) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                         CorFunBohman(lagt[i], maxtime_val);
            }
            const double l2 = one_log_dpoisgamma(ww, muj, gamma_param);
            const double biv = biv_PoissonGamma(corr, uu, ww, mui, muj, gamma_param);
            const double bl = log(biv) - l2;

            total += bl * weights;
        }
    }
    *res = total;
    if(!R_FINITE(*res)) {
        *res = LOW;
    }
}
/*********************************************************/
void Comp_Cond_Pois_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                          double *par, int *weigthed, double *res, double *mean1, double *mean2,
                          double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // ==================== 1. INITIALIZATION ====================
    const double nugget = nuis[0];
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];
    const double nugget_factor = 1.0 - nugget;
    if(nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }
    double total = 0.0;
    for(int i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double mui = exp(mean1[i]);
            const double muj = exp(mean2[i]);
            const int uu = (int)d1;
            const int ww = (int)d2;
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0) * nugget_factor;
            double weights = 1.0;
            if(weigthed_val) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                         CorFunBohman(lagt[i], maxtime_val);
            }

            const double l2 = dpois(ww, muj, 1);
            const double biv = biv_Poisson(corr, uu, ww, mui, muj);
            const double bl = log(biv) - l2;

            total += weights * bl;
        }
    }
    *res = total;
    if(!R_FINITE(*res)) {
        *res = LOW;
    }
}
/******************************************************************************************/
void Comp_Cond_BinomGauss_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    const double nugget = nuis[0];
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;
    const int trials = N1[0];  // Number of binomial trials
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];
    const double nugget_factor = 1.0 - nugget;
    if(nugget >= 1 || nugget < 0) {
        *res = LOW;
        return;
    }
    double total = 0.0;
    for(int i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Precompute probabilities
            const double ai = mean1[i];
            const double aj = mean2[i];
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0) * nugget_factor;
            const double p11 = pbnorm22(ai, aj, corr);
            const double p1 = pnorm(ai, 0, 1, 1, 0);
            const double p2 = pnorm(aj, 0, 1, 1, 0);
            const int uu = (int)d1;
            const int vv = (int)d2;
            double weights = 1.0;
            if(weigthed_val) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                         CorFunBohman(lagt[i], maxtime_val);
            }
            const double l2 = dbinom(vv, trials, p2, 1);
            const double biv = biv_binom(trials, uu, vv, p1, p2, p11);
            const double bl = log(biv) - l2;

            total += weights * bl;
        }
    }
    *res = total;
    if(!R_FINITE(*res)) {
        *res = LOW;
    }
}
/******************************************************************************************/
void Comp_Cond_BinomLogi_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0, uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    p11=pblogi22(log(exp(ai)+1),log(exp(aj)+1),(1-nugget)*corr);
                p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                u=data1[i];v=data2[i];
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                          uu=(int) u; vv=(int) v;
                           l2=dbinom(vv,N1[0],p2,1);
                         bl=log(biv_binom (N1[0],uu,vv,p1,p2,p11))-l2;
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Cond_BinomNNGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0, uu=0,vv=0,n1,n2;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                  n1=N1[i];n2=N2[i];
                 corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    p11=pbnorm22(ai,aj,(1-nugget)*corr);
                p1=pnorm(ai,0,1,1,0);
                p2=pnorm(aj,0,1,1,0);
                u=data1[i];v=data2[i];
                                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                            uu=(int) u; vv=(int) v;
                           l2=dbinom(vv,n2,p2,1);
                           bl=log(biv_binom222(n1,n2,uu,vv,p1,p2,p11))-l2;
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/******************************************************************************************/
void Comp_Cond_BinomNNLogi_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0, uu=0,vv=0,n1,n2;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                      n1=N1[i];n2=N2[i];
                 corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    p11=pblogi22(log(exp(ai)+1),log(exp(aj)+1),(1-nugget)*corr);
                p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                u=data1[i];v=data2[i];
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                            uu=(int) u; vv=(int) v;
                           l2=dbinom(vv,n2,p2,1);
                        bl=log(biv_binom222(n1,n2,uu,vv,p1,p2,p11))-l2;
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_BinomnegGauss_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                   double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                   double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{

    const double nugget = nuis[0];
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;
    const int trials = N1[0];  // Number of trials for negative binomial
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];
    const double nugget_factor = 1.0 - nugget;
    if(nugget >= 1 || nugget < 0) {
        *res = LOW;
        return;
    }
    double total = 0.0;
    for(int i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        if(!ISNAN(d1) && !ISNAN(d2)) {
            // Precompute probabilities
            const double ai = mean1[i];
            const double aj = mean2[i];
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0) * nugget_factor;
            const double p11 = pbnorm22(ai, aj, corr);
            const double p1 = pnorm(ai, 0, 1, 1, 0);
            const double p2 = pnorm(aj, 0, 1, 1, 0);
            const int uu = (int)d1;
            const int vv = (int)d2;
            double weights = 1.0;
            if(weigthed_val) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                         CorFunBohman(lagt[i], maxtime_val);
            }
            const double l2 = one_log_negbinom_marg(vv, trials, p2);
            const double biv = biv_binomneg(trials, uu, vv, p1, p2, p11);
            const double bl = log(biv) - l2;

            total += weights * bl;
        }
    }
    *res = total;
    if(!R_FINITE(*res)) {
        *res = LOW;
    }
}



/*********************************************************/
void Comp_Cond_BinomnegLogi_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
       if( nugget>=1 || nugget<0){*res=LOW; return;}
    //compute the composite log-likelihood:

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    p11=pblogi22(log(exp(ai)+1),log(exp(aj)+1),(1-nugget)*corr);
                p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                    u=data1[i];v=data2[i];
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                          uu=(int) u;  vv=(int) v;
                         l2=one_log_negbinom_marg(vv,N1[0],p2);
                          bl=log(biv_binomneg(N1[0],uu,vv,p1,p2,p11))-l2;
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}



void Comp_Cond_TWOPIECETukeyh_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,eta,tail,qq,sill,nugget,l2=0.0;
    eta  = nuis[2];  //skewness parameter
    tail = nuis[3];  //tail parameter
    sill =nuis[1];
    nugget=nuis[0];

     if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1||tail<0||tail>0.5) {*res=LOW;  return;}

       qq=qnorm((1-eta)/2,0,1,1,0);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
           zi=data1[i];zj=data2[i];


           corr=CorFct(cormod,lags[i],lagt[i],par,0,0);

           // l1=one_log_two_pieceTukey(zi,mean1[i],sill,tail,eta);
            l2=one_log_two_pieceTukey(zj,mean2[i],sill,tail,eta);

           p11=pbnorm22(qq,qq,corr);
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
          // bl=2*log(biv_two_pieceTukeyh((1-nugget)*corr,zi,zj,sill,eta,tail,p11,mean1[i],mean2[i]))-(l1+l2);
           bl=log(biv_two_pieceTukeyh((1-nugget)*corr,zi,zj,sill,eta,tail,p11,mean1[i],mean2[i]))-l2;
               *res+= weights*bl;
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}void Comp_Cond_TWOPIECET_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    int i;
    double bl, corr, zi, zj, weights = 1.0, p11, qq, l2 = 0.0;
    double eta = nuis[3];  // skewness parameter
    double sill = nuis[2];
    double nugget = nuis[1];
    double df = nuis[0];

    // Pre-check to avoid unnecessary computations
    if (sill < 0 || nugget < 0 || nugget >= 1 || fabs(eta) > 1 || df > 0.5 || df < 0) {
        *res = LOW;
        return;
    }

    qq = qnorm((1 - eta) / 2, 0, 1, 1, 0);

    // Pre-calculate npairs[0] once if it is used multiple times
    int npairs_val = npairs[0];
    
    for (i = 0; i < npairs_val; i++) {
        // Check for NaN in data1 and data2 before processing
        if (!ISNAN(data1[i]) && !ISNAN(data2[i])) {
            zi = data1[i];
            zj = data2[i];
            corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0);
            l2 = one_log_two_pieceT(zj, mean2[i], sill, df, eta);
            p11 = pbnorm22(qq, qq, corr);
            if (*weigthed) {
                weights = CorFunBohman(lags[i], maxdist[0]) * CorFunBohman(lagt[i], maxtime[0]);
            }
            bl = log(biv_two_pieceT(corr, zi, zj, sill, df, eta, p11, mean1[i], mean2[i], nugget)) - l2;
            *res += weights * bl;
        }
    }
    if (!R_FINITE(*res)) {
        *res = LOW;
    }
    return;
}

void Comp_Cond_TWOPIECEGauss_st2mem(int *cormod, double *data1, double *data2, int *N1, int *N2,
                                   double *par, int *weigthed, double *res, double *mean1, double *mean2,
                                   double *nuis, int *local, int *GPU, int *type_cop, int *cond)
{
    // ============= 1. INITIALIZATION AND PARAMETER SETUP =============
    const double eta = nuis[2];    // Skewness parameter
    const double sill = nuis[1];
    const double nugget = nuis[0];
    const int npairs_val = npairs[0];
    const int weigthed_val = *weigthed;
    const double maxdist_val = maxdist[0];
    const double maxtime_val = maxtime[0];
    const double nugget_factor = 1.0 - nugget;
    const double qq = qnorm((1.0 - eta)/2.0, 0.0, 1.0, 1, 0);
    if(fabs(eta) > 1 || sill < 0 || nugget < 0 || nugget >= 1) {
        *res = LOW;
        return;
    }
    double total = 0.0;

    for(int i = 0; i < npairs_val; i++) {
        const double d1 = data1[i];
        const double d2 = data2[i];
        if(!ISNAN(d1) && !ISNAN(d2)) {
            const double m1 = mean1[i];
            const double m2 = mean2[i];
            const double corr = CorFct(cormod, lags[i], lagt[i], par, 0, 0) * nugget_factor;
            const double p11 = pbnorm22(qq, qq, corr);
            double weights = 1.0;
            if(weigthed_val) {
                weights = CorFunBohman(lags[i], maxdist_val) * 
                         CorFunBohman(lagt[i], maxtime_val);
            }
            const double l2 = one_log_two_pieceGauss(d2, m2, sill, eta);
            const double biv = biv_two_pieceGaussian(corr, d1, d2, sill, eta, p11, m1, m2);
            const double bl = log(biv) - l2;

            total += weights * bl;
        }
    }
    *res = total;
    if(!R_FINITE(*res)) {
        *res = LOW;
    }
}


void Comp_Cond_TWOPIECEBIMODAL_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,eta,qq,sill,df,nugget,delta ,l2=0.0;

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
                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
              //  l1=one_log_bomidal(zi,mean1[i],sill,df,delta,eta);
                l2=one_log_bomidal(zj,mean2[i],sill,df,delta,eta);
                p11=pbnorm22(qq,qq,corr);
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                    /********************************************************/
                 // bl=2*log(biv_two_piece_bimodal((1-nugget)*corr,zi,zj,sill,df,delta,eta,p11,mean1[i],mean2[i]))-(l1+l2);
                  bl=log(biv_two_piece_bimodal((1-nugget)*corr,zi,zj,sill,df,delta,eta,p11,mean1[i],mean2[i]))-l2;

                    /********************************************************/
                           *res+= weights*bl;
                }}}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}

/********************************************************/
void Comp_Cond_BinomnegGaussZINB_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0;
    double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];
      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    u=data1[i];v=data2[i];
                                  if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                          uu=(int) u;vv=(int) v;
                    l2=one_log_BinomnegZIP(vv,N1[0],aj,mup);

                    bl=log(biv_binomnegZINB(N1[0],corr,uu,vv,ai,aj,nugget1,nugget2,mup))-l2;
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/************************************************/
void Comp_Cond_PoisZIP_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0, uu,vv;
    double weights=1.0,corr,mui,muj,bl,l2=0.0,u,v;
   double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];


      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                        u=data1[i];v=data2[i];
                                      if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                          uu=(int) u;
                         vv=(int) v;

                   // l1=one_log_PoisZIP(uu,mui,mup);
                    l2=one_log_PoisZIP(vv,muj,mup);
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                  //   bl=2*log(biv_PoissonZIP(corr,uu,vv,mui, muj,mup,nugget1,nugget2))-(l1+l2);
                       bl=log(biv_PoissonZIP(corr,uu,vv,mui, muj,mup,nugget1,nugget2))-l2;
                      *res+= bl*weights;
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Cond_Gauss_misp_PoisZIP_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{
    int i=0;
    double weights=1.0,corr,mui,muj,bl ,l2=0.0;
    double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];
    double p=pnorm(mup,0,1,1,0);

      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],lagt[i],par,0,0);


                   //   l1=dnorm(data1[i],(1-p)*mui,sqrt(mui*(1-p)*(1+p*mui)),1);
                      l2=dnorm(data2[i],(1-p)*muj,sqrt(muj*(1-p)*(1+p*muj)),1);
if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
           // bl=2*log(biv_Mis_PoissonZIP(corr,data1[i],data2[i],mui, muj,mup,nugget1,nugget2))-(l1+l2);
            bl=log(biv_Mis_PoissonZIP(corr,data1[i],data2[i],mui, muj,mup,nugget1,nugget2))-l2;

                      *res+= bl*weights;
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}


void Comp_Cond_LogLogistic_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{

    int i;double corr,zi,zj,weights=1.0,bl=1.0,l2=0.0;
    double nugget=nuis[0];
    if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);

                   // l1=one_log_loglogistic(zi,exp(mean1[i]),nuis[2]);
                    l2=one_log_loglogistic(zj,exp(mean2[i]),nuis[2]);

        //bl=2*log(biv_LogLogistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2])) -(l1+l2);
        bl=log(biv_LogLogistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2])) -l2;
             if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);

  *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}





void Comp_Cond_Logistic_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond)
{

    int i;double corr,zi,zj,weights=1.0,bl=1.0,l2=0.0;
    double nugget=nuis[0];
    if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    //l1=one_log_logistic(zi,mean1[i],nuis[1]) ;
                    l2=one_log_logistic(zj,mean2[i],nuis[1])  ;
                    //bl=2*log(biv_Logistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[1])) -(l1+l2);
                 bl=log(biv_Logistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[1])) -l2;
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);

  *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}






