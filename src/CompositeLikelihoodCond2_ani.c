#include "header.h"

//******************************************************************************************/
/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
/******************************************************************************************/
void Comp_Cond_Gauss2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
      int i=0;double lag=0.0;
    double  weights=1.0,sill,nugget,corr,bl,l2;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}

    for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                      corr=CorFct(cormod,lag,0,par,0,0);
                      //Rprintf(" %d- %f %f \n",i,lag,corr);
                       if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                             //weights=CorFunW_gen(lag,6, 3, maxdist[0]);
                      bl=log_biv_Norm((1-nugget)*corr,data1[i],data2[i],mean1[i],mean2[i],sill,0);
                      l2= dnorm(data2[i], mean2[i],sqrt(sill),1);
                       *res+= (bl-l2)*weights;
                    }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/*********************************************************/
void Comp_Cond_Tukeyh2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;double bl,corr,weights=1.0,l2=0.0;
    double sill=nuis[1];
    double nugget=nuis[0];
    double tail=nuis[2];
      if( sill<0||tail<0||tail>0.5||nugget<0||nugget>=1){*res=LOW; return;}
    for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                corr=CorFct(cormod,lag,0,par,0,0);
                if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                 bl=log(biv_tukey_h((1-nugget)*corr,data1[i],data2[i],mean1[i],mean2[i],tail,sill));
                 l2=one_log_tukeyh(data2[i],mean2[i],sill,tail);
                           *res+= (bl-l2)*weights;
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Tukeyhh2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;double bl,corr,weights=1.0 ,l2=0.0;
    double sill=nuis[1];
    double nugget=nuis[0];
    double h1=nuis[3];
    double h2=nuis[2];
      if( sill<0||h1<0||h1>0.5||h2<0||h2>0.5||nugget<0||nugget>=1){*res=LOW; return;}
        for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
             
  lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                corr=CorFct(cormod,lag,0,par,0,0);
                l2=one_log_tukeyhh(data2[i],mean2[i],sill,h1,h2);
               if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
               bl=log(biv_tukey_hh((1-nugget)*corr,data1[i],data2[i],mean1[i],mean2[i],sill,h1,h2))-l2;
                             *res+= weights*bl;
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}

/******************************************************************************************/
void Comp_Cond_SkewGauss2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
      double sill=nuis[1],  skew=nuis[2], nugget=nuis[0],l2=0.0,bb=0.0;
     if(nugget<0|| nugget>=1||sill<0){*res=LOW;  return;}
    int i=0;double lag=0.0;double corr,zi,zj,weights=1.0;
      for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=data1[i];zj=data2[i];
                      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                    corr=CorFct(cormod,lag,0,par,0,0);
    l2=one_log_SkewGauss(zj,mean2[i],sill,skew);
                    if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
        bb=log(biv_skew(corr,zi,zj,mean1[i],mean2[i],sill,skew,nugget))-l2;
                  *res+= weights*bb;
                 }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_T2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;double bl,corr,zi,zj,qi,qj,weights=1.0,l2=0.0;
    double sill=nuis[2];
    double nugget=nuis[1];
    double df=nuis[0];
    double df1=1/nuis[0];
      if( sill<0||nugget<0||nugget>=1||df<0||df>0.5){*res=LOW; return;}

   for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                qi=(zi-mean1[i])/sqrt(sill);qj=(zj-mean2[i])/sqrt(sill);
                  lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                corr=CorFct(cormod,lag,0,par,0,0);
                  l2=one_log_T(zj,mean2[i],sill,df1);
              if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                    bl=log(biv_T(corr,qi,qj,df,nugget)/sill)-l2;
                *res+= weights*bl;
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Gauss_misp_T2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis ,int *type_cop, int *cond)
{
    int i=0;double lag=0.0;
    double weights=1.0,corr,df=0.0,bl,l2=0.0,var=0.0;

     double sill=nuis[2];
    double nugget=nuis[1];


    if( sill<0||nugget<0||nugget>=1||nuis[0]<0||nuis[0]>0.5){*res=LOW; return;}
    df=1/nuis[0];
    var=sill*df/(df-2);
     for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
  lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
           corr=(1-nugget)*CorFct(cormod,lag,0,par,0,0);
           corr=exp(log(df-2)+2*lgammafn(0.5*(df-1))-(log(2)+2*lgammafn(df/2))+log(hypergeo(0.5,0.5, df/2,corr*corr))+log(corr*(1-nugget)));
         l2=dnorm(data2[i],data2[i],sqrt(sill*df/(df-2)),1);
         if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                    bl=log_biv_Norm(corr,data1[i],data2[i],mean1[i],mean2[i],var,0)-l2;
                       *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
void Comp_Cond_Gauss_misp_SkewT2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;
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

     for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                     corr=CorFct(cormod,lag,0,par,0,0)*(1-nugget);
                     corr2=corr_skewt(corr,df,skew);
                     if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                          bl=log_biv_Norm(corr2,data1[i],data2[i],mean1[i]+sqrt(sill)*MM,
                                                                  mean2[i]+sqrt(sill)*MM,
                                                                  sill*FF,0);
                          l2=dnorm(data2[i],mean2[i]+sqrt(sill)*MM,sqrt(sill*FF),1);
                        *res+= (bl-l2)*weights;


                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
/*********************************************************/
void Comp_Cond_Gauss_misp_Tukeygh2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;double bl,corr,corr2,zi,zj,weights=1.0,eta,tail,sill,nugget,u,eta2,mu,vv,l2;
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
   for(i=0;i<*npairs;i++){
          zi=data1[i];zj=data2[i];
if(!ISNAN(zi)&&!ISNAN(zj) ){
  lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                    corr=(1-nugget)*CorFct(cormod,lag,0,par,0,0);
                    corr2=corr_tukeygh(corr,eta,tail);
                    if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                bl=log_biv_Norm(corr2,zi,zj,mean1[i]+sqrt(sill)*mu,
                                            mean2[i]+sqrt(sill)*mu, sill*vv,0);
                      l2= dnorm(zj, mean2[i]+sqrt(sill)*mu,sqrt(sill*vv),1);
                      *res+= (bl-l2)*weights;

                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Cond_SinhGauss2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

        int i=0;double lag=0.0;double corr,zi,zj,bb=0.0,l2=0.0,weights=1.0;
           if(nuis[3]<0||nuis[1]<0||nuis[0]<0||nuis[0]>=1) {*res=LOW;  return;}

   for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=data1[i];zj=data2[i];
                      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                    corr=CorFct(cormod,lag,0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                    bb=log(biv_sinh((1-nuis[0])*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],nuis[1]));
                    l2=one_log_sas(zj,mean2[i],nuis[2],nuis[3],nuis[1]);
                    *res+= (weights*bb-l2);
                 }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Gamma2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;double corr,zi,zj,weights=1.0,bl=0.0,l2=0.0;
      double nugget=nuis[0];
     if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}

     for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=(data1[i]);zj=(data2[i]);
                  lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                    corr=CorFct(cormod,lag,0,par,0,0);
                        if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                    l2=one_log_gamma(zj,mean2[i],nuis[2]);
             bl=log(biv_gamma((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]))- l2;
                     *res+= weights*bl;
                }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
void Comp_Cond_Weibull2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;double corr,zi,zj,weights=1.0,bl=0.0,l2=0.0;
      double nugget=nuis[0];
     if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}

     for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=(data1[i]);zj=(data2[i]);
                  lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                    corr=CorFct(cormod,lag,0,par,0,0);
                        if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                    l2=one_log_weibull(zj,mean2[i],nuis[2]);
             bl=log(biv_Weibull((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]))- l2;
                     *res+= weights*bl;
                }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_LogGauss2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;double corr,zi,zj,weights=1.0,bl=0.0,l2=0.0;
    double sill=nuis[1];double nugget=nuis[0];

    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}
      for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                    corr=CorFct(cormod,lag,0,par,0,0);
                     zi=(data1[i]);zj=(data2[i]);
                    l2=one_log_loggaussian(zj,mean2[i],sill);
                    if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                 bl=log(d2lognorm(zi,zj,sill,nugget, mean1[i], mean2[i],(1-nugget)*corr))-l2;
                    *res+= weights*bl;
                    }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Beta2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;double corr,zi,zj,weights=1.0,bl,l2=0.0;
    double nugget=nuis[0];
    double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max)  {*res=LOW;  return;}
    for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=data1[i]; zj=data2[i];
                      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                    corr=CorFct(cormod,lag,0,par,0,0);
                     if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                    l2=one_log_beta(zj,nuis[2],nuis[3],min,max);
                bl=log(biv_beta((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-l2;
        *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Kumaraswamy2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;double corr,zi,zj,weights=1.0,bl,l2=0.0;
    double nugget=nuis[0];
         double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}

    for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
               zi=data1[i];zj=data2[i];
                 lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                    corr=CorFct(cormod,lag,0,par,0,0);
                    l2=one_log_kumma(zj,mean2[i],nuis[2],nuis[3],min,max);
                     if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                  bl=log(biv_Kumara((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-l2;
        *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Kumaraswamy22mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;double corr,zi,zj,weights=1.0,bl,l2=0.0;
    double nugget=nuis[0];
         double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}
    for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                    corr=CorFct(cormod,lag,0,par,0,0);
                    l2=one_log_kumma2(zj,mean2[i],nuis[2],nuis[3],min,max);
                    if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
    bl=log(biv_Kumara2((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-l2;
        *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Gauss_misp_Pois2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0,N=2;
    double  weights=1.0,corr,corr1,mui,muj,bl,l2,lag=0.0;
    double nugget=nuis[0];
      if(nugget<0||nugget>=1){*res=LOW; return;}
double **M;
        M= (double **) R_Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) R_Calloc(N,double);}
    double *dat;
    dat=(double *) R_Calloc(N,double);
    for(i=0;i<*npairs;i++){

                  if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
             //***********/
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                     corr=CorFct(cormod,lag,0,par,0,0)*(1-nugget);
                      corr1=corr_pois(corr,mui, muj);
                      if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                        M[0][0]=mui; M[1][1]=muj;M[0][1]=sqrt(mui*muj)*corr1;M[1][0]= M[0][1];
                        dat[0]=data1[i]-mui;dat[1]=data2[i]-muj;
                      l2=dnorm(data2[i],muj,sqrt(muj),1);;
                       bl=log(dNnorm(N,M,dat))-l2;

                      *res+= bl*weights;
                    }}
   for(i=0;i<N;i++)  {R_Free(M[i]);}
    R_Free(M);

    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_BinomNNGauss_misp2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, N=2,n1,n2;
    double u,v,m1,m2,l2,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,lag=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success

    double **M;
    M= (double **) R_Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) R_Calloc(N,double);}
    double *dat;
    dat=(double *) R_Calloc(N,double);

    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                   lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                 corr=CorFct(cormod,lag,0,par,0,0);
                 p11=pbnorm22(ai,aj,(1-nugget)*corr);
                 p1=pnorm(ai,0,1,1,0);
                 p2=pnorm(aj,0,1,1,0);
                 u=data1[i];v=data2[i];
                 n1=N1[i];n2=N2[i];
                 m1=n1*p1;m2=n2*p2;
                 if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                 M[0][0]=m1*(1-p1);   M[1][1]=m2*(1-p2);  // var1 var2
                 M[0][1]= fmin_int(n1,n2)*(p11-p1*p2) ;       // covariance
                 M[1][0]= M[0][1];
                 dat[0]=u-m1;dat[1]=v-m2; 
                 l2=dnorm(v,m2,sqrt(m2*(1-p2)),1);; 
                 bl= log(dNnorm(N,M,dat)) -l2; 
                 *res+= bl*weights;       
                }}
    for(i=0;i<N;i++)  {R_Free(M[i]);}
    R_Free(M);
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Cond_Gauss_misp_PoisGamma2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0,N=2;
    double  weights=1.0,corr,corr1,mui,muj,bl,l2,bi,bj,vvi,vvj,lag=0.0;
    double nugget=nuis[0];
      if(nugget<0||nugget>=1){*res=LOW; return;}
double **M;
        M= (double **) R_Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) R_Calloc(N,double);}
    double *dat;
    dat=(double *) R_Calloc(N,double);
    for(i=0;i<*npairs;i++){

                  if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
             //***********/
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                    bi= nuis[2]/mui; bj= nuis[2]/muj;
                    vvi= mui*(1+1/bi); vvj= muj*(1+1/bj);
                      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                     corr=CorFct(cormod,lag,0,par,0,0)*(1-nugget);
                      corr1=corr_pois_gen(corr,mui, muj, nuis[2]);
                      if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                        M[0][0]=vvi; M[1][1]=vvj;M[0][1]=sqrt(vvi*vvj)*corr1;M[1][0]= M[0][1];
                        dat[0]=data1[i]-mui;dat[1]=data2[i]-muj;
                      l2=dnorm(data2[i],muj,sqrt(vvj),1);;
                      bl=log(dNnorm(N,M,dat))-l2;
                      *res+= bl*weights;
                    }}
   for(i=0;i<N;i++)  {R_Free(M[i]);}
    R_Free(M);

    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
void Comp_Cond_PoisGamma2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu,ww;
    double weights=1.0,corr,mui,muj,bl,l2,lag=0.0;
    double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}
  // Rprintf("%d   \n",*npairs);
      for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                     corr=CorFct(cormod,lag,0,par,0,0);
                    // if(fabs(corr)>1|| !R_FINITE(corr)) {*res=LOW; return;}
                        if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
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
/*********************************************************/
void Comp_Cond_Pois2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu,ww;
    double weights=1.0,corr,mui,muj,bl,l2,lag=0.0;
    double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}
  // Rprintf("%d   \n",*npairs);
      for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                     corr=CorFct(cormod,lag,0,par,0,0);
                    // if(fabs(corr)>1|| !R_FINITE(corr)) {*res=LOW; return;}
                        if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
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
void Comp_Cond_BinomGauss2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0,lag=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}

    for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                   lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                 corr=CorFct(cormod,lag,0,par,0,0);
    p11=pbnorm22(ai,aj,(1-nugget)*corr);
    //Rprintf("p11: %f\n",p11);
                p1=pnorm(ai,0,1,1,0);
                p2=pnorm(aj,0,1,1,0);
                u=data1[i];v=data2[i];
                        if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                          uu=(int) u; vv=(int) v;
                          l2=dbinom(vv,N1[0],p2,1);
                        bl=log(biv_binom (N1[0],uu,vv,p1,p2,p11))-l2;
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Cond_BinomLogi2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0,lag=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}

    for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                   lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                 corr=CorFct(cormod,lag,0,par,0,0);
    p11=pblogi22(log(exp(ai)+1),log(exp(aj)+1),(1-nugget)*corr);
                p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                u=data1[i];v=data2[i];
                        if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                          uu=(int) u; vv=(int) v;
                          l2=dbinom(vv,N1[0],p2,1);
                        bl=log(biv_binom (N1[0],uu,vv,p1,p2,p11))-l2;
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
void Comp_Cond_BinomNNGauss2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu=0,vv=0,n1,n2;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0,lag=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    for(i=0;i<*npairs;i++){
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
                 l2=dbinom(vv,n2,p2,1);
                 bl=log(biv_binom222(n1,n2,uu,vv,p1,p2,p11))-l2;
                 *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Cond_BinomNNLogi2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0, uu=0,vv=0,n1,n2;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0,lag=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                   lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                 corr=CorFct(cormod,lag,0,par,0,0);
                 p11=pblogi22(log(exp(ai)+1),log(exp(aj)+1),(1-nugget)*corr);
                 p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                 u=data1[i];v=data2[i];
                 n1=N1[i];n2=N2[i];
                 if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                 uu=(int) u; vv=(int) v;
                 //l1=dbinom(uu,n1,p1,1);
                 l2=dbinom(vv,n2,p2,1);
                 //bl=2*log(biv_binom222(n1,n2,uu,vv,p1,p2,p11))-(l1+l2);
                 bl=log(biv_binom222(n1,n2,uu,vv,p1,p2,p11))-l2;
                 *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_BinomnegGauss2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0,lag=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
       if( nugget>=1 || nugget<0){*res=LOW; return;}
    //compute the composite log-likelihood:

    for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                 corr=CorFct(cormod,lag,0,par,0,0);
                    p11=pbnorm22(ai,aj,(1-nugget)*corr);
                p1=pnorm(ai,0,1,1,0);p2=pnorm(aj,0,1,1,0);
                    u=data1[i];v=data2[i];
                         if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                          uu=(int) u;vv=(int) v;
                         l2=one_log_negbinom_marg(vv,N1[0],p2);
                         bl=log(biv_binomneg(N1[0],uu,vv,p1,p2,p11))-l2;
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/*********************************************************/
void Comp_Cond_BinomnegLogi2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0,lag=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
       if( nugget>=1 || nugget<0){*res=LOW; return;}
    //compute the composite log-likelihood:

    for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                 corr=CorFct(cormod,lag,0,par,0,0);
                    p11=pblogi22(log(exp(ai)+1),log(exp(aj)+1),(1-nugget)*corr);
                p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                    u=data1[i];v=data2[i];
                         if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                          uu=(int) u;vv=(int) v;
                         l2=one_log_negbinom_marg(vv,N1[0],p2);
                        bl=log(biv_binomneg(N1[0],uu,vv,p1,p2,p11))-l2;
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Cond_TWOPIECETukeyh2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;double bl,corr,zi,zj,weights=1.0,p11,eta,tail,qq,sill,nugget,l2=0.0;
    eta  = nuis[2];  //skewness parameter
    tail = nuis[3];  //tail parameter
    sill =nuis[1];
    nugget=nuis[0];

     if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1||tail<0||tail>0.5) {*res=LOW;  return;}

       qq=qnorm((1-eta)/2,0,1,1,0);
      for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
           zi=data1[i];zj=data2[i];

       lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
           corr=CorFct(cormod,lag,0,par,0,0);

            //l1=one_log_two_pieceTukey(zi,mean1[i],sill,tail,eta);
            l2=one_log_two_pieceTukey(zj,mean2[i],sill,tail,eta);

           p11=pbnorm22(qq,qq,corr);
           if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
          // bl=2*log(biv_two_pieceTukeyh((1-nugget)*corr,zi,zj,sill,eta,tail,p11,mean1[i],mean2[i]))-(l1+l2);
           bl=log(biv_two_pieceTukeyh((1-nugget)*corr,zi,zj,sill,eta,tail,p11,mean1[i],mean2[i]))-l2;
               *res+= weights*bl;
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_TWOPIECET2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;double bl,corr,zi,zj,weights=1.0,p11,qq,l2=0.0;
    double eta=nuis[3];  //skewness parameter
    double sill=nuis[2];
    double nugget=nuis[1];
    double df=nuis[0];
    if(sill<0||nugget<0||nugget>=1 ||fabs(eta)>1|| df >0.5||df<0) {*res=LOW;  return;}
       qq=qnorm((1-eta)/2,0,1,1,0);
      for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                  lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                corr=CorFct(cormod,lag,0,par,0,0);
                //l1=one_log_two_pieceT(zi,mean1[i],sill,df,eta);
                l2=one_log_two_pieceT(zj,mean2[i],sill,df,eta);
                 p11=pbnorm22(qq,qq,corr);
                 if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
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
void Comp_Cond_TWOPIECEGauss2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;double bl,corr,zi,zj,weights=1.0,p11,eta,qq,sill,nugget,l2=0.0;
    eta=nuis[2];  //skewness parameter
    sill=nuis[1];
    nugget=nuis[0];
       qq=qnorm((1-eta)/2,0,1,1,0);
         if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1) {*res=LOW;  return;}
   for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                  lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                    corr=CorFct(cormod,lag,0,par,0,0);
               // l1=one_log_two_pieceGauss(zi,mean1[i],sill,eta);
                l2=one_log_two_pieceGauss(zj,mean2[i],sill,eta);

                p11=pbnorm22(qq,qq,corr);
                if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                   // bl=2*log(biv_two_pieceGaussian((1-nugget)*corr,zi,zj,sill,eta,p11,mean1[i],mean2[i]))-l1-l2;
                    bl=log(biv_two_pieceGaussian((1-nugget)*corr,zi,zj,sill,eta,p11,mean1[i],mean2[i]))-l2;
                 
                    *res+= weights*bl;
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/********************************************************/
void Comp_Cond_TWOPIECEBIMODAL2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;double bl,corr,zi,zj,weights=1.0,p11,eta,qq,sill,df,nugget,delta ,l2=0.0;

    eta=nuis[4];  //skewness parameter
    delta=nuis[3];
    sill=nuis[2];
    nugget=nuis[1];
    df=nuis[0];
    //alpha=2*(delta+1)/df;
 if( fabs(eta)>1||df<0||nugget<0||nugget>=1||delta<0||sill<0) {*res=LOW;  return;}
    qq=qnorm((1-eta)/2,0,1,1,0);
    for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                if(!ISNAN(zi)&&!ISNAN(zj) ){
                      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                    corr=CorFct(cormod,lag,0,par,0,0);
             //   l1=one_log_bomidal(zi,mean1[i],sill,df,delta,eta);
                l2=one_log_bomidal(zj,mean2[i],sill,df,delta,eta);
                p11=pbnorm22(qq,qq,corr);
                 if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                    /********************************************************/
                 // bl=2*log(biv_two_piece_bimodal((1-nugget)*corr,zi,zj,sill,df,delta,eta,p11,mean1[i],mean2[i]))-(l1+l2);
                  bl=log(biv_two_piece_bimodal((1-nugget)*corr,zi,zj,sill,df,delta,eta,p11,mean1[i],mean2[i]))-l2;
                           *res+= weights*bl;
                }}}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/********************************************************/
void Comp_Cond_BinomnegGaussZINB2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l2=0.0,lag=0.0;
    double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];
      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
    for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                    lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                 corr=CorFct(cormod,lag,0,par,0,0);
                    u=data1[i];v=data2[i];
                         if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                          uu=(int) u;vv=(int) v;
                    l2=one_log_BinomnegZIP(vv,N1[0],aj,mup);
                    bl=log(biv_binomnegZINB(N1[0],corr,uu,vv,ai,aj,nugget1,nugget2,mup))-l2;
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;

}
/************************************************/
void Comp_Cond_PoisZIP2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i=0, uu,vv;
    double weights=1.0,corr,mui,muj,bl,l2=0.0,u,v,lag=0.0;
   double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];


      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
      for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                     corr=CorFct(cormod,lag,0,par,0,0);
                        u=data1[i];v=data2[i];
                         if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                          uu=(int) u;vv=(int) v;

                 //   l1=one_log_PoisZIP(uu,mui,mup);
                    l2=one_log_PoisZIP(vv,muj,mup);
                        if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                   //  bl=2*log(biv_PoissonZIP(corr,uu,vv,mui, muj,mup,nugget1,nugget2))-(l1+l2);
                      bl=log(biv_PoissonZIP(corr,uu,vv,mui, muj,mup,nugget1,nugget2))-l2;
                      *res+= bl*weights;
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Cond_Gauss_misp_PoisZIP2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{
    int i=0;double lag=0.0;
    double weights=1.0,corr,mui,muj,bl ,l2=0.0;
    double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];
    double p=pnorm(mup,0,1,1,0);

      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}

  // Rprintf("%d   \n",*npairs);
      for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                     corr=CorFct(cormod,lag,0,par,0,0);


                      //l1=dnorm(data1[i],(1-p)*mui,sqrt(mui*(1-p)*(1+p*mui)),1);
                      l2=dnorm(data2[i],(1-p)*muj,sqrt(muj*(1-p)*(1+p*muj)),1);


                      if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
                      //bl=2*log(biv_Mis_PoissonZIP(corr,data1[i],data2[i],mui, muj,mup,nugget1,nugget2))-(l1+l2);
                      bl=log(biv_Mis_PoissonZIP(corr,data1[i],data2[i],mui, muj,mup,nugget1,nugget2))-l2;

                      *res+= bl*weights;
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Cond_LogLogistic2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i=0;double lag=0.0;double corr,zi,zj,weights=1.0,bl=1.0,l2=0.0;
    double nugget=nuis[0];
    if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
     for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                    corr=CorFct(cormod,lag,0,par,0,0);
                    l2=one_log_loglogistic(zj,exp(mean2[i]),nuis[2]);
                  bl=log(biv_LogLogistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]))-l2;
                           
                         if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);

  *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}



void Comp_Cond_Logistic2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *type_cop, int *cond)
{

    int i=0;double lag=0.0;double corr,zi,zj,weights=1.0,bl=1.0,l2=0.0;
    double nugget=nuis[0];
    if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
     for(i=0;i<*npairs;i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                      lag= hypot(coord1[2*i]-coord2[2*i],coord1[2*i+1]-coord2[2*i+1]); 
                    corr=CorFct(cormod,lag,0,par,0,0);
                    l2=one_log_logistic(zj,mean2[i],nuis[1])  ;
                    bl=log(biv_Logistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[1])) -l2;
                         if(*weigthed) weights=CorFunBohman(lag,maxdist[0]);
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


