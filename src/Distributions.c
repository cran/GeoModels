#include "header.h"
//******************************************************************************
double polevl(double x, const double coef[], int N)

{
    double ans;
    int i;
    const double *p;
    p = coef;
    ans = *p++;
    i = N;
    do
        ans = ans * x + *p++;
    while (--i);   
    return (ans);
}
double p1evl(double x, const double coef[], int N)
{
    double ans;
    const double *p;
    int i;
    p = coef;
    ans = x + *p++;
    i = N - 1;
    do
        ans = ans * x + *p++;
    while (--i);
    return (ans);
}
/************************************************************************/
 int is_nonpos_int(double x)
{
    return x <= 0 && x == ceil(x) && fabs(x) < 1e13;
}
 
double poch(double a, double m)
{
    double r;
    r = 1.0;
    /* Recurse down */
    while (m >= 1.0) {
        if (a + m == 1) {
            break;
        }
        m -= 1.0;
        r *= (a + m);
        if (!R_FINITE(r) || r == 0) {
            break;
        }
    }

    /* Recurse up */
    while (m <= -1.0) {
        if (a + m == 0) {
            break;
        }
        r /= (a + m);
        m += 1.0;
        if (!R_FINITE(r) || r == 0) {
            break;
        }
    }

    if (m == 0) {
        /* Easy case */
        return r;
    }
    else if (a > 1e4 && fabs(m) <= 1) {
        /* Avoid loss of precision */
        return r * R_pow(a, m) * (
            1
            + m*(m-1)/(2*a)
            + m*(m-1)*(m-2)*(3*m-1)/(24*a*a)
            + m*m*(m-1)*(m-1)*(m-2)*(m-3)/(48*a*a*a)
            );
    }

    /* Check for infinity */
    if (is_nonpos_int(a + m) && !is_nonpos_int(a) && a + m != m) {
        return INFINITY;
    }

    /* Check for zero */
    if (!is_nonpos_int(a + m) && is_nonpos_int(a)) {
        return 0;
    }

    return(r * exp(lgammafn(a + m) - lgammafn(a)) * sign(gammafn(a + m)) * sign(gammafn(a)));
}
/************************************************************************/






// ===================================== START: Bivariate Normal  =====================================//

#define HSQRT 1.414213562373095048801688724209698078569671
//#define HSQRT 1.4142

double Phi(double x)
{
    double val =(1+     (1-erfc(x/HSQRT) )    )/2;
    return ( val );
}
double Phi2diag( double x, double a, double px, double pxs )
{
    double sol=NAN;
    if( a <= 0.0 ) sol = px;
    if( a >= 1.0 ) sol =  px * px;
    double b = 2.0 - a, sqrt_ab = sqrt( a * b );
    double c1 = 6.36619772367581343e-001;
    double c2 = 1.25331413731550025;
    double c3 = 1.57079632679489662;
    double c4 = 1.591549430918953358e-001;
    //double asr = ( a > 0.1 ? asin( 1.0 - a ) : acos( sqrt_ab ) );
    double asr;
    if(a > 0.1)
    {
        asr = asin( 1.0 - a );
    }else
    {
        asr = acos( sqrt_ab );
    }
    
    double comp = px * pxs;
    if( comp * ( 1.0 - a - c1 * asr ) < 5e-17 ) sol =  b * comp;
    double tmp = c2 * x;
    double alpha = a * x * x / b;
    double a_even = -tmp * a;
    double a_odd = -sqrt_ab * alpha;
    double beta = x * x;
    double b_even = tmp * sqrt_ab;
    double b_odd = sqrt_ab * beta;
    double delta = 2.0 * x * x / b;
    double d_even = ( 1.0 - a ) * c3 - asr;
    double d_odd = tmp * ( sqrt_ab - a );
    double res = 0.0, res_new = d_even + d_odd;
    int k = 2;
    double cond = fabs(res-res_new);
    while( cond>DEPSILON )
    {
        d_even = ( a_odd + b_odd + delta * d_even ) / k;
        a_even *= alpha / k;
        b_even *= beta / k;
        k++;
        a_odd *= alpha / k;
        b_odd *= beta / k;
        d_odd = ( a_even + b_even + delta * d_odd ) / k;
        k++;
        res = res_new;
        res_new += d_even + d_odd;
        cond = fabs(res-res_new);
    }
    double sol1;
    if(isnan(sol))
    {
        res *= exp( -x * x / b ) * c4;
        sol1 =  fmax( ( 1.0 + c1 * asr ) * comp, b * comp - fmax( 0.0, res ) );
    }else{
        sol1 = sol;
    }
    return sol1;
}


double Phi2help( double x, double y, double rho )
{
    double s = sqrt( ( 1.0 - rho ) * ( 1.0 + rho ) );
    double a = 0.0, b1 = -fabs( x ), b2 = 0.0;
    if( rho > 0.99 )
    {
        double tmp = sqrt( ( 1.0 - rho ) / ( 1.0 + rho ) );
        b2 = -fabs( ( x - y ) / s - x * tmp );
        a = R_pow( ( x - y ) / x / s - tmp,2 );
    }
    else if( rho < -0.99 )
    {
        double tmp = sqrt( ( 1.0 + rho ) / ( 1.0 - rho ) );
        b2 = -fabs( ( x + y ) / s - x * tmp );
        a = R_pow( ( x + y ) / x / s - tmp,2 );
    }
    else
    {
        b2 = -fabs( rho * x - y ) / s;
        a = R_pow( b2 / x ,2);
    }
    
    double p1 = Phi( b1 ), p2 = Phi( b2 ), q = 0.0;
    if( a <= 1.0 )
        q = 0.5 * Phi2diag( b1, 2.0 * a / ( 1.0 + a ), p1, p2 );
    else
        q = p1 * p2 - 0.5 * Phi2diag( b2, 2.0 / ( 1.0 + a ), p2, p1 );
    int c1 = ( y / x >= rho ), c2 = ( x < 0.0 ), c3 = c2 && ( y >= 0.0 );
    
    bool c13 = (c1 && c3);
    bool c12 = (c1 && c2);
    double sol;
    if(c13) {sol = (q - 0.5);}
    else if(c12) {sol = (q);}
    else if(c1 ) {sol = (0.5 - p1 + q);}
    else if(c3 ) {sol = (p1 - q - 0.5);}
    else if(c2 ) {sol = (p1 - q);}
    else {sol = (0.5 - q);}
    
    //return (0.5 - p1 + q );
    return (sol );
    //return ( c1 && c3 ? q - 0.5
    //        : c1 && c2 ? q
    //        : c1 ? 0.5 - p1 + q
    //        : c3 ? p1 - q - 0.5
    //        : c2 ? p1 - q
    //        : 0.5 - q );
}

double Phi2( double x, double y, double rho )
{
    double sol = NAN;
    if( ( 1.0 - rho ) * ( 1.0 + rho ) <= 0.0 )
    {
        if( rho > 0.0 )
        {
            //return (Phi( min( x, y ) ));
            sol = (Phi( fmin( x, y ) ));
        }
        else
        {
            //return (max( 0.0, min( 1.0, Phi( x ) + Phi( y ) - 1.0 ) ));
            sol = (fmax( 0.0, fmin( 1.0, Phi( x ) + Phi( y ) - 1.0 ) ));
        }
    }
    if( x == 0.0 && y == 0.0 )
    {
        if( rho > 0.0 )
        {
            //return (Phi2diag( 0.0, 1.0 - rho, 0.5, 0.5 ));
            sol = (Phi2diag( 0.0, 1.0 - rho, 0.5, 0.5 ));
        }
        else
        {
            //return (0.5 - Phi2diag( 0.0, 1.0 + rho, 0.5, 0.5 ));
            sol = (0.5 - Phi2diag( 0.0, 1.0 + rho, 0.5, 0.5 ));
        }
    }
    else
    {
        sol = (fmax( 0.0,
                   fmin( 1.0,
                       Phi2help( x, y, rho ) + Phi2help( y, x, rho ) ) ));
    }
    
    return (sol);
}



/*for bivariate t distributions*/
double A[] = {
    8.11614167470508450300E-4,
    -5.95061904284301438324E-4,
    7.93650340457716943945E-4,
    -2.77777777730099687205E-3,
    8.33333333333331927722E-2
};

double B[] = {
    -1.37825152569120859100E3,
    -3.88016315134637840924E4,
    -3.31612992738871184744E5,
    -1.16237097492762307383E6,
    -1.72173700820839662146E6,
    -8.53555664245765465627E5
};

double C[] = {
    /* 1.00000000000000000000E0, */
    -3.51815701436523470549E2,
    -1.70642106651881159223E4,
    -2.20528590553854454839E5,
    -1.13933444367982507207E6,
    -2.53252307177582951285E6,
    -2.01889141433532773231E6
};





 /****** integrand  in  generalized wendland function ************************/
double int_onef2(double x,double a, double b1,double b2,double y)
{
    double res=0.0,of1;
    of1=gammafn(b1)*R_pow(x*y,1-b1/2)*bessel_i(2*sqrt(x*y),b1-1,1);
    res=R_pow(1-x,b2-a-1)*R_pow(x,a-1)*of1;
    return (res);
}
void integr_onef2(double *x, int n, void *ex){
    int i;double a,b1,b2,y;
    a =    ((double*)ex)[0];  //smooth
    b1 = ((double*)ex)[1];  //alpha
    b2 =     ((double*)ex)[2];  //h
    y =     ((double*)ex)[3];  //h
    for (i=0;i<n;i++) {x[i]=int_onef2(x[i],a,b1,b2,y);}
    return;
}
// function computing generalized wendland
double onef2integral(double x, double *param) {
    
    double ex[4], lower, upper, epsabs, epsrel, result, abserr, *work;
    int neval, ier, subdiv, lenw, last, *iwork;
    subdiv = 100;
    epsabs = R_pow(DBL_EPSILON, 0.25);
    epsrel = epsabs;
    lenw = 4 * subdiv;           /* as instructed in WRE */
    iwork =   (int *) R_Calloc(subdiv, int);  /* idem */
    work = (double *) R_Calloc(lenw, double); /* idem */
    ex[0] = param[0]; ex[1] = param[1];ex[2]= param[2];ex[3]=x;
    lower=0;
    upper=1;
    // Compute the integral
    Rdqags(integr_onef2, (void *) &ex,
               &lower, &upper, &epsabs, &epsrel, &result,
               &abserr, &neval, &ier, &subdiv, &lenw, &last, iwork, work);
    R_Free(iwork);R_Free(work);
    return(result);
}

/******************************************************************************/
/******************************************************************************/
/* Wendland covariance */

/* integrand  in  generalized wendland function*/
double int_gen(double x,double mu, double alpha,double lag,double supp)
{
    double res=0.0,y;
    y=lag/supp;
    res=R_pow(1-x,mu-1)*R_pow(x*x-y*y,alpha)/beta(2*alpha+1,mu);
    return (res);///(R_pow(2,alpha-1)*gamma(alpha)*R_pow(supp,2*alpha)));
}
void integr_gen(double *x, int n, void *ex){
    int i;double mu,alpha,beta,y;
    mu =    ((double*)ex)[0];  //mu
    alpha = ((double*)ex)[1];  //alpha
    beta =     ((double*)ex)[2];  //csupp
    y =     ((double*)ex)[3];  //h
    for (i=0;i<n;i++) {x[i]=int_gen(x[i],mu,alpha,y,beta);}
    return;
}
// function computing generalized wendland
double wendintegral(double x, double *param) {
    double ex[4], lower, upper, epsabs, epsrel, result, abserr, *work;
    int neval, ier, subdiv, lenw, last, *iwork;
    subdiv = 100;
    epsabs = R_pow(DBL_EPSILON, 0.25);
    epsrel = epsabs;
    lenw = 4 * subdiv;           /* as instructed in WRE */
    iwork =   (int *) R_Calloc(subdiv, int);  /* idem */
    work = (double *) R_Calloc(lenw, double); /* idem */
    ex[0] = param[0]; ex[1] = param[1]; ex[2] = param[2];ex[3]=x;
    lower=x/param[2];
    upper=1;
    // Compute the integral
    if(x<=param[2]) {
    Rdqags(integr_gen, (void *) &ex,
               &lower, &upper, &epsabs, &epsrel, &result,
               &abserr, &neval, &ier, &subdiv, &lenw, &last, iwork, work);

    }else   {result=0;}
    R_Free(iwork);R_Free(work);
    return(result);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
// Integrand function (derivatives of the Student-t cdf):
double int_pt(double x, double df)
  {
    double res=0.0, x2=0.0, y=0.0;
    x2=R_pow(x,2);
    y=1+x2/df;
    res=0.5*dt(x,df,0)*((df+1)*x2/R_pow(df,2)/y-log(y));
    return res;
  }
// Vectorised integrand function:
void integr_pt(double *x, int n, void *ex)
{
  int i=0;
  double d=0.0;
  d=*((double *)ex);
  for(i=0;i<n;i++)
    x[i]=int_pt(x[i],d);
  return;
}


// bivariate of Gaussian bivariate rf
double log_biv2gauss(int *cormod, double dij,double *par, double data1, double data2, int first,int second)
{
double rhott,rhovv,rhotv,det,dens=0.0;
rhott=CorFct(cormod,0,0,par,first,first);
rhovv=CorFct(cormod,0,0,par,second,second);
rhotv=CorFct(cormod,dij,0,par,first,second);
det=rhott*rhovv-R_pow(rhotv,2);
dens=-0.5*(2*log(2*M_PI)+log(det)+(rhovv*R_pow(data1,2)+rhott*R_pow(data2,2)-2*(data1*data2)*rhotv)/det);
return dens;
}


// compute  bivariate normal standard pdf:
double d2norm(double x, double y, double rho)
{
  double res=0.0, omr=1-R_pow(rho,2);
  res=(1/(2*M_PI))*exp(-0.5*(R_pow(x,2)-2*rho*x*y+R_pow(y,2))/omr)/sqrt(omr);
  return(res);
}
// compute  bivariate normal standard pdf:
double d22norm(double x, double y,double v11,double v22,double v12)
{
  double res=0.0;
  double cc=sqrt(v11*v22);
  double rho=v12/cc;
  double omr=1-R_pow(rho,2);
  double aa= 2*M_PI*cc*sqrt(omr);
  double zz=R_pow(x,2)/v11  + R_pow(y,2)/v22-2*rho*x*y/cc;
  res=exp(-0.5*zz/omr)/aa;
  return(res);
}

/*
double d2lognorm(double x, double y, double sill,double nugget, double mux,double muy,double rho)
{
  double KK=exp(sill/2);
  x=x*KK; y=y*KK;
  double res=0.0, q=0.0, omr=R_pow(sill,2)-R_pow(rho*sill,2);
  q=(sill*R_pow((log(x)-mux),2) + sill*R_pow((log(y)-muy),2)
    -2*rho*sill*(log(x)-mux)*(log(y)-muy))/omr;
  res=exp(-q/2)/(2*x*y*M_PI*sqrt(omr));
  return(res*R_pow(KK,2));
}
*/



double d2lognorm(double x, double y, double sill, double nugget, double mux, double muy, double rho)
{

    
    double sigma = sqrt(sill);
    double sigma2 = sill;
    double mu_i = exp(mux);double mu_j = exp(muy);
    double v_i = x / mu_i;double v_j = y / mu_j;
    double rho_eff = rho;  // Correlazione effettiva
    double one_minus_rho2 = 1.0 - rho_eff * rho_eff;
    double term_i = (log(v_i) + sigma2/2.0) / sigma;
    double term_j = (log(v_j) + sigma2/2.0) / sigma;
    double quad = (term_i*term_i + term_j*term_j - 2.0*rho_eff*term_i*term_j) / (2.0 * one_minus_rho2);
    double prefactor = 1.0 / (2.0 * M_PI * sigma2 * v_i * v_j * sqrt(one_minus_rho2));
    double density_V = prefactor * exp(-quad);
    double result = density_V / (mu_i * mu_j);
    return result;
}


double biv_sinh(double corr,double zi,double zj,double mi,double mj,double skew,double tail,double vari){

double xi=(zi-mi)/sqrt(vari);
double xj=(zj-mj)/sqrt(vari);
double asinh_xi=asinh(xi);
double asinh_xj=asinh(xj);
double b1=tail*asinh_xi-skew;
double b2=tail*asinh_xj-skew;
double sinh_b1=sinh(b1);
double sinh_b2=sinh(b2);
double cosh_b1=cosh(b1);
double cosh_b2=cosh(b2);
double Z1=sinh_b1;
double Z2=sinh_b2;
double one_minus_corr2=1.0-corr*corr;
double xi2p1=xi*xi+1.0;
double xj2p1=xj*xj+1.0;
double sqrt_k=sqrt(one_minus_corr2);
double denom=2.0*M_PI*sqrt_k;
double A=(cosh_b1*cosh_b2*tail*tail)/(denom*sqrt(xi2p1*xj2p1));
double exponent=-(Z1*Z1+Z2*Z2-2.0*corr*Z1*Z2)/(2.0*one_minus_corr2);
double B=exp(exponent);
return(A*B)/vari;
}


// cdf of  a bivariate Gausssian distribution
void cdf_norm_call(double *lim1,double *lim2,double *a11,double *a12, double *res)
{
    *res=   cdf_norm(lim1[0],lim2[0],a11[0],a12[0]);
    //*res=   cdf_norm(0,0,.5,1);
}

void CBessel(double *xxx, double *nuu, int *expscale,double *res, int *tipo)
{
    switch (*tipo) {
        case 1:
            *res = bessel_k(*xxx,*nuu,*expscale);
            break;
        case 2:
            *res = bessel_i(*xxx,*nuu,*expscale);
        default:
            break;
}}


double biv_chisqu2(double corr,double zi,double zj, double shape)
{
double KK1,KK2,KK3,rr;
   rr=1-corr*corr;  
   KK1=R_pow(gammafn(shape/2),2)*R_pow(rr,shape/2);
   KK2=R_pow(2,-shape)*R_pow(zi*zj,shape/2-1)*exp(-0.5*(zi+zj)/rr)/KK1;
   KK3= KK2* gammafn(shape/2)*R_pow(0.5*fabs(corr)*sqrt(zi*zj)/rr,1-shape/2)*bessel_i(fabs(corr)*sqrt(zi*zj)/rr ,shape/2-1,1);
   return(KK3);
}


double biv_two_piece_bimodal(double rho,double zi,double zj,double sill,double nu,double delta,double eta,
             double p11,double mui,double muj)
{
double res=0.0;
double alpha=2*(delta+1)/nu;
double etamas=1+eta;
double etamos=1-eta;
double zistd=(zi-mui)/sqrt(sill);
double zjstd=(zj-muj)/sqrt(sill);

double nn=R_pow(2,1-alpha/2);

if(zi>=mui&&zj>=muj)
{res=          (R_pow(zistd*zjstd,alpha-1)*p11/R_pow(etamos,2*alpha))*biv_gamma_gen(rho,R_pow(zistd/etamos,alpha),R_pow(zjstd/etamos,alpha),0,0,nu,nn);  }
if(zi>=mui&&zj<muj)
{res=R_pow(-zistd*zjstd,alpha-1)*(1-eta-2*p11)/(R_pow(etamos,alpha)*R_pow(etamas,alpha))*biv_gamma_gen(rho,R_pow(zistd/etamos,alpha),R_pow(-zjstd/etamas,alpha),0,0,nu,nn);}
if(zi<mui&&zj>=muj)
{res=R_pow(-zistd*zjstd,alpha-1)*(1-eta-2*p11)/(R_pow(etamos,alpha)*R_pow(etamas,alpha))*biv_gamma_gen(rho,R_pow(-zistd/etamas,alpha),R_pow(zjstd/etamos,alpha),0,0,nu,nn);}
if(zi<mui&&zj<muj)
{res=     (R_pow(zistd*zjstd,alpha-1)*(p11+eta)/R_pow(etamas,2*alpha))*biv_gamma_gen(rho,R_pow(-zistd/etamas,alpha),R_pow(-zjstd/etamas,alpha),0,0,nu,nn);}
return(R_pow(alpha,2)*res/sill);
}



double biv_Weibull2(double rho12,double zi,double zj,double mi,double mj,double shape1,double shape2){

double a1=gammafn(1.0+1.0/shape1);
double a2=gammafn(1.0+1.0/shape2);
double mui=exp(mi);
double muj=exp(mj);
double k=1.0-rho12*rho12;
double zi_mui=zi/mui;
double zj_muj=zj/muj;
double a1zi=zi_mui*a1;
double a2zj=zj_muj*a2;
double pow1=R_pow(a1zi,shape1);
double pow2=R_pow(a2zj,shape2);
double pow1_sqrt=R_pow(a1zi,shape1/2.0);
double pow2_sqrt=R_pow(a2zj,shape2/2.0);
double a=shape1*shape2*R_pow(a1,shape1)*R_pow(a2,shape2)*R_pow(zi,shape1-1.0)*R_pow(zj,shape2-1.0)/(R_pow(mui,shape1)*R_pow(muj,shape2)*k);
double b=exp(-(pow1+pow2)/k);
double c=bessel_i(2.0*fabs(rho12)*pow1_sqrt*pow2_sqrt/k,0,1);
return a*b*c;
}
double asy_log_besselI(double z,double nu)
{
     double val;
     double K=4*R_pow(nu,2);
     val =  (z-0.5*log(2*M_PI*z))+
           log((1-(K-1)/(8*z)*(1-(K-9)/(2*8*z)*(1-(K-25)/(3*8*z)))));
     return(val);
}



double biv_Weibull(double corr,double zi,double zj,double mui,double muj,double shape){

double ci=exp(mui);
double cj=exp(muj);
double k=1.0/gammafn(1.0+1.0/shape);
double ui=zi/ci;
double uj=zj/cj;
double a=1.0-corr*corr;
double uipow=R_pow(ui,shape);
double ujpow=R_pow(uj,shape);
double uipow_half=R_pow(ui,shape/2.0);
double ujpow_half=R_pow(uj,shape/2.0);
double kpow=R_pow(k,-shape);
double z=2.0*fabs(corr)*uipow_half*ujpow_half*kpow/a;
double A=2.0*log(shape)-2.0*shape*log(k)+(shape-1.0)*log(ui*uj)-log(a);
double B=-kpow*(uipow+ujpow)/a;
double res=A+B+log(bessel_i(z,0,2))+z-(mui+muj);
return exp(res);
}
/*******************************************/

/***

double psi(int i, int j, double p1, double p2, double p12){
    double aux= p1 + p2 - p12;
    if(i==0 || j==0){return 0.0;}
    if(i==1 && j==1){return (p1+p2-p1*p2)/(p2*p1*aux);}
    double aux1= (psi( i-1, j-1, p1, p2, p12)*p12+psi( i, j-1, p1, p2, p12)*(p1-p12)+psi( i-1, j, p1, p2, p12)*(p2-p12))/aux;
    double aux2= ( (i-p2/aux)/p2 + (j-p1/aux)/p1 )/aux;
    double aux3= (2-aux)/R_pow(aux,2);
    return(aux1+aux2+aux3);
}
***/

double psi(int i,int j,double p1,double p2,double p12,double**cache){
    if(i==0||j==0)return 0.0;
    if(i==1&&j==1)return(p1+p2-p1*p2)/(p2*p1*(p1+p2-p12));
    if(cache[i][j]!=-1.0)return cache[i][j];

    double aux=p1+p2-p12;
    double psi11=psi(i-1,j-1,p1,p2,p12,cache);
    double psi01=psi(i,j-1,p1,p2,p12,cache);
    double psi10=psi(i-1,j,p1,p2,p12,cache);
    double aux1=(psi11*p12+psi01*(p1-p12)+psi10*(p2-p12))/aux;
    double aux2=((i-p2/aux)/p2+(j-p1/aux)/p1)/aux;
    double aux3=(2.0-aux)/(aux*aux);
    cache[i][j]=aux1+aux2+aux3;
    return cache[i][j];
}




double** create_cache(int max_i,int max_j){
    double** cache=(double**)malloc((max_i+1)*sizeof(double*));
    for(int i=0;i<=max_i;i++){
        cache[i]=(double*)malloc((max_j+1)*sizeof(double));
        for(int j=0;j<=max_j;j++)cache[i][j]=-1.0;
    }
    return cache;
}



double cov_binom_neg(int m,double p11,double p1, double p2) {

double** cache=create_cache(m,m);
double a;
if(fabs(p11-p1*p2)<1e-32) return(0.0); //indipendence case
else    
 a=psi( m, m, p1, p2, p11,cache)- m*m/(p1*p2);
return(a);
}

/*****************************************/
/********* poisson gamma correlation******/
/*****************************************/
double corrPGs(double corr,double mean,double a){

double corr2=corr*corr;
double nu=a/mean;
double KK=nu*(1-corr2);
double cc=4.0/R_pow(2.0+KK,2);
double dd=exp(log(nu)+0.5*(log(nu)+log1p(-corr2))+a*log(2.0+KK)-log1p(nu)-(a+0.5)*log(4.0+KK));
double aa=hypergeo((1.0-a)/2.0,-a/2.0,1,cc);
double bb=((a+1.0)/(2.0+KK))*hypergeo((2.0-a)/2.0,(1.0-a)/2.0,2,cc);
return corr2*(1.0-dd*(aa+bb));
}
   

double CorrPGns(double corr,double mean_i,double mean_j,double a){

double rho2=corr*corr;
double nu_i=a/mean_i;
double nu_j=a/mean_j;
double lnuij=log(nu_i)+log(nu_j);
double lpnuij=log1p(nu_i)+log1p(nu_j);
double mnui=-1.0/nu_i;
double mnuj=-1.0/nu_j;
double sum=0.0;
int iter1=1000,iter2=1000;
for(int r=0;r<iter1;r++){
int r2=r+2;
for(int m=0;m<iter2;m++){
double a1m=1.0-a-m;
double h1=hypergeo(1.0,a1m,r2,mnui);
double h2=hypergeo(1.0,a1m,r2,mnuj);
if(h1<=0.0||h2<=0.0)break;
double H=log(h1)+log(h2)-2.0*lgammafn(r2);
double aux1=2.0*m*log(corr)+m*lnuij+2.0*lgammafn(r+m+1+a);
double aux2=(r+m)*lpnuij+lgammafn(m+a)+lgammafn(m);
double term=H+aux1-aux2;
if(!R_finite(term)||fabs(term)<1e-100)break;
sum+=term;
}
}
double aux3=exp(log(rho2)+(a+1.0)*log1p(-rho2)+(a-0.5)*lnuij-lgammafn(a)-(a+0.5)*lpnuij-0.5*log(mean_i*mean_j));
double aux4=exp(log(rho2)+0.5*(lnuij-lpnuij));
return aux4+aux3*sum;
}

/******************************************/
double corr_pois_gen(double corr,double mean_i, double mean_j, double a){ 

      //if( (corr>(1-1e-56)) &&  corr<=1 &&  fabs(mean_i-mean_j)<1e-320){return(1.0);}
if(fabs(corr)<1e-200){return(0.0);}
    else{
    double res=0.0;
    if (fabs(mean_i-mean_j)<MACHEP) res=corrPGs(corr,mean_i,a);
    else                            res=CorrPGns(corr,mean_i,mean_j,a);
    return(res);
}
}

/******************************************/
/******** Poisson correlation *************/
/******************************************/
double corr_pois_s(double rho,double mi)
{
double res=0.0;
double rho2=R_pow(rho,2);
double gg=2*mi/(1-rho2);
res=rho2*(1-exp(-gg)*(bessel_i(gg,0,1)+bessel_i(gg,1,1)));
return(res);
}   
/******************************************/
double corr_pois_ns(double rho,double mi,double mj){
int r=0;
double res0=0.0,sum=0.0;
double rho2=rho*rho;
double ki=mi/(1.0-rho2);
double kj=mj/(1.0-rho2);
double K=rho2*(1.0-rho2)/sqrt(mi*mj);
while(r<8000){
    double aa=r+1;
    sum+=exp(log(igam(aa,ki))+log(igam(aa,kj)));
    if(fabs(sum-res0)<MACHEP) break;
    res0=sum;
    r++;
}
return sum*K;
}
/******************************************/
double corr_pois(double rho,double mi,double mj){
    if( (rho>(1-1e-6)) &&  rho<=1 &&  fabs(mi-mj)<1e-320){return(1.0);}
if(fabs(rho)<1e-10){return(0.0);}
    else{
    double res=0.0;
    if (fabs(mi-mj)<MACHEP)         res=corr_pois_s(rho,mi);
    else                            res=corr_pois_ns(rho,mi,mj);
    return(res);
}
}
/*****************************************/
/*****************************************/

double corr_tukeygh(double rho,double eta,double tail){

    if(fabs(rho)<1e-16) return 0.0;
    double rho2 = rho * rho;
    double tail2 = tail * tail;
    double eta2 = eta * eta;
    double u = 1.0 - tail;
    double a = 1.0 + rho;
    double mu, cova, vari, rho1 = 0.0;

    if(fabs(eta) > MACHEP){
        double A1 = exp(a * eta2 / (1.0 - tail * a));
        double A2 = 2.0 * exp(0.5 * eta2 * (1.0 - tail * (1.0 - rho2)) / (u * u - tail2 * rho2));
        double A3 = eta2 * sqrt(u * u - rho2 * tail2);
        mu = (exp(eta2 / (2.0 * u)) - 1.0) / (eta * sqrt(u));
        cova = ((A1 - A2 + 1.0) / A3 - mu * mu);
        vari = ((exp(2.0 * eta2 / (1.0 - 2.0 * tail)) - 2.0 * exp(eta2 / (2.0 * (1.0 - 2.0 * tail))) + 1.0) / 
               (eta2 * sqrt(1.0 - 2.0 * tail)) - mu * mu);
        rho1 = cova / vari;
    } else {
        double inv_tail = R_pow(1.0 - 2.0 * tail, -1.5);
        rho1 = (-rho / ((1.0 + tail * (rho - 1.0)) * (-1.0 + tail + tail * rho) * 
                        sqrt(1.0 + tail * (-2.0 + tail - tail * rho2)))) / inv_tail;
    }

    return rho1;
}




double corr_skewt(double corr,double df,double skew)
{
if(fabs(corr)<MACHEP){return(0.0);}
    else{
    double d,w,corr1,skew2,CorSkew,nu,l,y,KK;
    skew2=skew*skew; nu=df; l=df/2;  d=(df-1)/2; w=sqrt(1-skew2);
    y=corr;
    CorSkew=(2*skew2/(M_PI*w*w+skew2*(M_PI-2)))*(sqrt(1-y*y)+y*asin(y)-1)+w*w*y/(w*w+skew2*(1-2/M_PI)) ;
   
    if(df<99){
    KK=exp(log(M_PI)+log(nu-2)+2*lgammafn(d)
    - log(2)-log(M_PI*R_pow(gammafn(l),2)*(1+skew2) - skew2*(nu-2)*R_pow(gammafn(d),2)));
    corr1=KK * (hypergeo(0.5,0.5,l,y*y)*((1+skew2*(1-2/M_PI))*CorSkew+2*skew2/M_PI) - 2*skew2/M_PI);}
    else {corr1=CorSkew;}
return(corr1);
}
}

double biv_gamma(double corr,double zi,double zj,double mui,double muj,double shape) {
    double a = 1.0 - corr * corr;
    double ci = zi / exp(mui);
    double cj = zj / exp(muj);
    double gam = gammafn(shape / 2.0);
    double res = 0.0;

    if (fabs(corr) > MACHEP) {
        double z = shape * fabs(corr) * sqrt(ci * cj) / a;
        double A = (shape / 2.0 - 1.0) * log(ci * cj) - shape * (ci + cj) / (2.0 * a);
        double C = (1.0 - shape / 2.0) * log(z / 2.0);
        double B = log(gam) + (shape / 2.0) * log(a) + shape * log(2.0) - shape * log(shape);
        double D = log(bessel_i(z, shape / 2.0 - 1.0, 2)) + z;
        res = (A + C + D) - (mui + muj + B);
        return exp(res);
    } else {
        double B = R_pow(shape / (2.0 * exp(mui)), shape / 2.0) * R_pow(zi, shape / 2.0 - 1.0) * exp(-shape * ci / 2.0) / gam;
        double C = R_pow(shape / (2.0 * exp(muj)), shape / 2.0) * R_pow(zj, shape / 2.0 - 1.0) * exp(-shape * cj / 2.0) / gam;
        res = B * C;
    }
    return res;
}


double biv_gamma_gen(double corr,double zi,double zj,double mui, double muj, double shape,double n)
{
    double a=0.0,A=0.0,D=0.0,z=0.0,res=0.0,B=0.0,C=0.0;
    double ci=zi/exp(mui);double cj=zj/exp(muj);
    double gam = gammafn(shape/2);
    a=1-R_pow(corr,2);  
  if (fabs(corr) > MACHEP) {
        z=n*fabs(corr)*sqrt(ci*cj)/a;
        A=(shape/2-1)*log(ci*cj) -n*(ci+cj)/(2*a); ///ok
        C=(1-shape/2)*log(z/2); //ok
        B=log(gam)+(shape/2)*log(a)+shape*log(2)-shape*log(n);  
        D=log(bessel_i(z,shape/2-1,2))+z; //ok
        res=(A+C+D)-(mui+muj+B);
        return(exp(res));
    }
    else
    {
        B=(R_pow((n/(2*exp(mui))),shape/2)*R_pow(zi,shape/2-1)*exp(-(n*ci/(2))))/gam;
        C=(R_pow((n/(2*exp(muj))),shape/2)*R_pow(zj,shape/2-1)*exp(-(n*cj/(2))))/gam;
        res=B*C;    
    }
return(res);
}



void biv_gamma_call(double *corr,double *zi,double *zj,double *mui, double *muj, double *shape, double *res)
{
    *res = biv_gamma(*corr,*zi,*zj,*mui,*muj,*shape);
}

/***********/
double Hyp_conf_laplace_approx(double a, double b, double z)
{
double y=0.0,r=0.0,A=0.0;
  y=2*a/(b-z+sqrt(R_pow(z-b,2)+ 4*a*z ));

  r=(R_pow(y,2)/a)+(R_pow(1-y,2)/(b-a));

  A=R_pow(b,b-0.5)*R_pow(r,-0.5)*R_pow((y/a),a)*R_pow((1-y)/(b-a),b-a)*exp(z*y);
  return(A);
}
/*****************************************************/
double asym_aprox_1F1(double a, double b, double z)
{
double k1,k2,k3,res=0.0;
double alpha=(b-a-1)/z;
if(b<z+a+1)
{
  k1=(gammafn(b)*exp(z)*R_pow(1-alpha,a-1))/(gammafn(a)*gammafn(b-a));
  k2=sqrt(2*M_PI*alpha/z)*R_pow(alpha/exp(1),alpha*z);
  k3=(2-a*alpha)*(1-a)/(2*R_pow(1-alpha,2));
  res=k1*k2*(1+ (k3+ 1/(12*alpha))/z);  
}
if(b>z+a+1)
{
k1=gammafn(b)/(gammafn(b-a)*R_pow(b-a-z-1,a));
k2=a*(a+1)*(a+1-b)/(2*R_pow(b-a-z-1,2));
res=k1*(1-k2);
}
if(b==z+a+1)
{k1=gammafn(b)/(2*gammafn(a)*gammafn(b-a));
 k2=R_pow(2/(b-a-1),a/2);
 k3=2*gammafn(0.5*(a+3))*sqrt(2/(b-a-1))/3;
res=k1*k2*(gammafn(a/2)-k3);}
return(res);
}

          
double biv_skew(double corr, double z1, double z2, double mi, double mj, double vari, double skew, double nugget) {
    double zi = z1 - mi;
    double zj = z2 - mj;

    double skew2 = skew * skew;
    double skew3 = skew2 * skew;
    double skew4 = skew2 * skew2;
    double vari2 = vari;
    double vari4 = vari2 * vari2;
    double corr1 = (1 - nugget) * corr;
    double corr2 = corr * corr;
    double corr12 = corr1 * corr1;

    double s2v2 = skew2 * vari2;
    double pp1 = skew * vari2 * corr1;
    double pp = skew * vari2 * corr;
    double ff = pp1 + pp;
    double mm = pp1 - pp;
    double aux1 = vari2 + skew2;

    double aux2 = corr1 * vari2 + corr * skew2;
    double pdf1 = d22norm(zi, zj, aux1, aux1, aux2);

    double bb = vari4 * corr12 + 2 * s2v2 * corr * corr1 + skew4 * corr2 - vari4 - 2 * s2v2 - skew4;
    double lim1 = ((mm) * zj + (pp * corr1 + skew3 * corr2 - skew * vari2 - skew3) * zi) / bb;
    double lim2 = ((mm) * zi + (pp * corr1 + skew3 * corr2 - skew * vari2 - skew3) * zj) / bb;
    double a11 = (vari4 * corr12 + s2v2 * corr2 - vari4 - s2v2) / bb;
    double a12 = (vari4 * corr * corr12 + (s2v2 * corr2 - s2v2) * corr1 - vari4 * corr) / bb;
    double cdf1 = cdf_norm(lim1, lim2, a11, a12);

    // PDF 2
    aux2 = vari2 * corr1 - skew2 * corr;
    double pdf2 = d22norm(zi, zj, aux1, aux1, aux2);

    bb = vari4 * corr12 - 2 * s2v2 * corr * corr1 + skew4 * corr2 - vari4 - 2 * s2v2 - skew4;
    lim1 = ((ff) * zj + (-pp * corr1 + skew3 * corr2 - skew * vari2 - skew3) * zi) / bb;
    lim2 = -(-(ff) * zi + (pp * corr1 - skew3 * corr2 + skew * vari2 + skew3) * zj) / bb;
    a11 = (vari4 * corr12 + s2v2 * corr2 - vari4 - s2v2) / bb;
    a12 = -(vari4 * corr * corr12 + (s2v2 - s2v2 * corr2) * corr1 - vari4 * corr) / bb;
    double cdf2 = cdf_norm(lim1, lim2, a11, a12);

    double dens = 2 * (pdf1 * cdf1 + pdf2 * cdf2);
    return dens;
}





double biv_wrapped(double alfa, double u, double v, double mi, double mj, double nugget, double sill, double corr) {
    double x = u - 2 * atan(mi) - M_PI;
    double y = v - 2 * atan(mj) - M_PI;

    double s1 = sill;
    double s12 = sill * corr;
    double det = R_pow(s1, 2) - R_pow(s12, 2);
    double wrap_gauss = 0.0;

    double k1 = -alfa, k2 = -alfa;

    double common_factor = 1 / (2 * M_PI * sqrt(det)); // Precompute common factor

    while (k1 <= alfa) {
        while (k2 <= alfa) {
            double term1 = x + 2 * k2 * M_PI;
            double term2 = y + 2 * k1 * M_PI;
            double quadr = -0.5 * (1.0 / det) * (s1 * R_pow(term2, 2) + s1 * R_pow(term1, 2) - 2.0 * s12 * term1 * term2);

            wrap_gauss += common_factor * exp(quadr);
            k2++;
        }
        k1++;
        k2 = -alfa;
    }
    return wrap_gauss;
}



double biv_skew2(double corr, double zi, double zj, double vari1, double vari2, double rho, double skew1, double skew2) {
    // Calcoli preliminari
    double taul = R_pow(vari1, 0.5);
    double taur = R_pow(vari2, 0.5);
    double c1 = 1.0 / (1 - R_pow(corr, 2));
    
    // Variabili di calcolo del determinante
    double aux1 = R_pow(taul, 2) + R_pow(skew1, 2);
    double aux2 = taul * taur * corr + skew1 * skew2 * corr;
    double aux3 = R_pow(taur, 2) + R_pow(skew2, 2);
    double det = aux1 * aux3 - R_pow(aux2, 2);

    // Calcolo quadratico e PDF 1
    double quadr = (1.0 / det) * (aux1 * R_pow(zj, 2) + aux3 * R_pow(zi, 2) - 2.0 * aux2 * zi * zj);
    double pdf1 = (0.5 / M_PI) * (1 / sqrt(det)) * exp(-0.5 * quadr);

    // Calcoli intermedi per la CDF 1
    aux1 = (R_pow(skew1 / taul, 2) + 1) * c1;
    aux2 = (skew1 * skew2) / (taul * taur) + 1 * c1 * corr;
    aux3 = (R_pow(skew2 / taur, 2) + 1) * c1;
    det = aux1 * aux3 - R_pow(aux2, 2);
    double a11 = (1 / det) * aux3;
    double a12 = (1 / det) * aux2;
    double a22 = (1 / det) * aux1;
    double factor = c1 / det;
    double fact1 = factor * (aux3 * (skew1 / R_pow(taul, 2)) - aux2 * (skew2 / (taul * taur)) * corr);
    double fact2 = factor * (-aux3 * (skew1 / (taul * taur)) * corr + aux2 * (skew2 / R_pow(taur, 2)));
    double fact3 = factor * (aux2 * (skew1 / R_pow(taul, 2)) - aux1 * (skew2 / (taul * taur)) * corr);
    double fact4 = factor * (-aux2 * (skew1 / (taul * taur)) * corr + aux1 * (skew2 / R_pow(taur, 2)));

    double lim1 = fact1 * zi + fact2 * zj;
    double lim2 = fact3 * zi + fact4 * zj;
    double cdf1 = cdf_norm2(lim1, lim2, a11, a12, a22);

    // Calcoli per il secondo PDF
    aux1 = R_pow(taul, 2) + R_pow(skew1, 2);
    aux2 = taul * taur * corr - skew1 * skew2 * corr;
    aux3 = R_pow(taur, 2) + R_pow(skew2, 2);
    det = aux1 * aux3 - R_pow(aux2, 2);
    quadr = (1.0 / det) * (aux1 * R_pow(zj, 2) + aux3 * R_pow(zi, 2) - 2.0 * aux2 * zi * zj);
    double pdf2 = (0.5 / M_PI) * (1 / sqrt(det)) * exp(-0.5 * quadr);

    // Calcoli intermedi per la CDF 2
    aux1 = R_pow(skew1 / taul, 2) * c1 + c1;
    aux2 = (skew1 * skew2) / (taul * taur) * c1 * corr - c1 * corr;
    aux3 = R_pow(skew2 / taur, 2) * c1 + c1;
    det = aux1 * aux3 - R_pow(aux2, 2);
    a11 = (1 / det) * aux3;
    a12 = (1 / det) * aux2;
    a22 = (1 / det) * aux1;
    factor = c1 / det;
    fact1 = factor * (aux3 * (skew1 / R_pow(taul, 2)) - aux2 * (skew2 / (taul * taur)) * corr);
    fact2 = factor * (-aux3 * (skew1 / (taul * taur)) * corr + aux2 * (skew2 / R_pow(taur, 2)));
    fact3 = factor * (aux2 * (skew1 / R_pow(taul, 2)) - aux1 * (skew2 / (taul * taur)) * corr);
    fact4 = factor * (-aux2 * (skew1 / (taul * taur)) * corr + aux1 * (skew2 / R_pow(taur, 2)));

    lim1 = fact1 * zi + fact2 * zj;
    lim2 = fact3 * zi + fact4 * zj;
    double cdf2 = cdf_norm2(lim1, lim2, a11, a12, a22);

    // Calcolo finale della densità
    double dens = 2 * (pdf1 * cdf1 + pdf2 * cdf2);
    return dens;
}



double log_biv_Norm(double corr, double zi, double zj, double mi, double mj, double vari, double nugget) {
    double u = zi - mi;
    double v = zj - mj;
    double u2 = u * u;
    double v2 = v * v;
    double s1 = vari;
    double s12 = s1 * corr * (1 - nugget);
    double det = s1 * s1 - s12 * s12;
    double dens = -0.5 * (2 * log(2 * M_PI) + log(det) + (s1 * (u2 + v2) - 2 * s12 * u * v) / det);
    
    return dens;
}

double biv_Norm(double corr, double zi, double zj, double mi, double mj, double vari1, double vari2, double nugget) {
    double u = zi - mi;
    double v = zj - mj;
    double u2 = u * u; // u^2
    double v2 = v * v; // v^2
    double s1 = sqrt(vari1 * vari2);
    double s12 = s1 * corr * (1 - nugget);
    double det = s1 * s1 - s12 * s12;
    double dens = -0.5 * (2 * log(2 * M_PI) + log(det) + (s1 * (u2 + v2) - 2 * s12 * u * v) / det);
    return exp(dens); // ritorna l'esponenziale della densità
}



/*multivariate gaussian PDF*/
double dNnorm(int N, double **M, double *dat) {
    double dd, pdf = 0.0;
    int m, n, i;
    int *indx = (int *) R_Calloc(N, int);
    double *col = (double *) R_Calloc(N, double);
    double **inverse = (double **) R_Calloc(N, double *);

    for (i = 0; i < N; i++) {
        inverse[i] = (double *) R_Calloc(N, double);
    }
    ludcmp(M, N, indx, &dd);
    dd = 1.0;
    for (m = 0; m < N; m++) {
        dd *= M[m][m];
    }
    for (n = 0; n < N; n++) {
        for (m = 0; m < N; m++) {
            col[m] = 0.0;
        }
        col[n] = 1.0;
        lubksb(M, N, indx, col);
        for (m = 0; m < N; m++) {
            inverse[m][n] = col[m];
        }
    }
    pdf = R_pow(2 * M_PI, -N / 2) * R_pow(dd, -0.5) * exp(-0.5 * QFORM(inverse, dat, dat, N));
    R_Free(indx);
    R_Free(col);
    for (i = 0; i < N; i++) {
        R_Free(inverse[i]);
    }
    R_Free(inverse);
    return pdf;
}


// cdf of  a bivariate Gausssian distribution
double cdf_norm(double lim1,double lim2,double a11,double a12)
{
    double res=0.0;
    double  lowe[2]={0,0}, uppe[2]={lim1/sqrt(a11),lim2/sqrt(a11)}, corre[1] ={a12/a11};
    int infin[2]={0,0};
    double auxil=1-R_pow(corre[0],2);
    double det=R_pow(a11,2)-R_pow(a12,2) ;
    res= a11* sqrt(auxil/det) *  F77_CALL(bvnmvn)(lowe,uppe,infin,corre);
    return(res);
    //return(0.5);
}


// compute the inverse lambert w transformation
double inverse_lamb(double x,double tail)
{
  double sign,value=0.0;
  value = sqrt(LambertW(tail*x*x)/tail);
  if (x >= 0) sign= 1;
  //if(x==0.0) sign=1;
  if (x < 0) sign= -1;
   return(sign*value);
}

double biv_tukey_hh(double corr, double data_i, double data_j, double mui, double muj,
                    double sill, double hl, double hr) {

    const double sqrt_sill = sqrt(sill);
    const double zi = (data_i - mui) / sqrt_sill;
    const double zj = (data_j - muj) / sqrt_sill;
    const double zi2 = zi * zi;
    const double zj2 = zj * zj;

    // Calcoli per hl e hr
    const double hl_zi = inverse_lamb(zi, hl);
    const double hl_zj = inverse_lamb(zj, hl);
    const double hr_zi = inverse_lamb(zi, hr);
    const double hr_zj = inverse_lamb(zj, hr);

    const double Lhl_zi = 1.0 + LambertW(hl * zi2);
    const double Lhl_zj = 1.0 + LambertW(hl * zj2);
    const double Lhr_zi = 1.0 + LambertW(hr * zi2);
    const double Lhr_zj = 1.0 + LambertW(hr * zj2);

    double res;

    if (fabs(corr) > 1e-30) {
        if (zi >= 0 && zj >= 0) {
            res = dbnorm(hr_zi, hr_zj, 0, 0, 1, corr) *
                  hr_zi * hr_zj / (zi * zj * Lhr_zi * Lhr_zj);
        } else if ((zi >= 0 && zj < 0) || (zi < 0 && zj >= 0)) {
            res = dbnorm(hr_zi, hl_zj, 0, 0, 1, corr) *
                  hr_zi * hl_zj / (zi * zj * Lhr_zi * Lhl_zj);
        } else {
            res = dbnorm(hl_zi, hl_zj, 0, 0, 1, corr) *
                  hl_zi * hl_zj / (zi * zj * Lhl_zi * Lhl_zj);
        }
    } else {
        double A = (zi >= 0) ? dnorm(hr_zi, 0, 1, 0) * hr_zi / (zi * Lhr_zi)
                             : dnorm(hl_zi, 0, 1, 0) * hl_zi / (zi * Lhl_zi);

        double B = (zj >= 0) ? dnorm(hr_zj, 0, 1, 0) * hr_zj / (zj * Lhr_zj)
                             : dnorm(hl_zj, 0, 1, 0) * hl_zj / (zj * Lhl_zj);

        res = A * B;
    }

    return res / sill;
}



// cdf of  a bivariate Gausssian distribution
double cdf_norm2(double lim1,double lim2,double a11,double a12, double a22)
{
  double res=0.0;

  double  lowe[2]={0,0}, uppe[2]={lim1/sqrt(a11),lim2/sqrt(a22)}, corre[1] ={a12/sqrt(a11*a22)};
  int infin[2]={0,0};
  double auxil=1-R_pow(corre[0],2);
  double det=a11*a22-R_pow(a12,2) ;
  res= sqrt(a11*a22)* sqrt(auxil/det) *  F77_CALL(bvnmvn)(lowe,uppe,infin,corre);
  return(res);
}



// compute the bivariate normal cdf for the bernoulli RF:
double pbnorm(int *cormod, double h, double u, double mean1, double mean2, 
  double nugget, double var,double *par, double thr)
{
  double res=0;
  double lim_inf[2]={0,0};//lower bound for the integration
  double lim_sup[2]={mean1,mean2};
  int infin[2]={0,0};//set the bounds for the integration
  double corr[1]={(1-nugget)*CorFct(cormod,h,u,par,0,0)};
    res=F77_CALL(bvnmvn)(lim_inf,lim_sup,infin,corr);
    //res = Phi2(lim_sup[0],lim_sup[1],corr[0]);
  return(res);
}



// compute a sequence (depending on space and time) of bivariate normal cdf:
void vpbnorm(int *cormod, double *h, double *u, int *nlags, int *nlagt,
         double *nuis, double *par, double *rho, double *thr)
{
  int i,j,t=0;
  for(j=0;j<*nlagt;j++)
    for(i=0;i<*nlags;i++){
      rho[t]=pbnorm(cormod,h[i],u[j],nuis[0],nuis[0],nuis[1],nuis[2],par,thr[0]);
      t++;}
  return;
}

/**********************************************************************/
double aux_biv_binomneg_simple(int NN, int u, double p01, double p10, double p11) {
    int i = 0;
    double kk1 = 0.0, dens = 0.0;
    double lgamma_NN_u_1 = lgammafn(NN - 1 + u + 1);
    double lgamma_u_NN_2 = lgammafn(u - NN + 2);
    for (i = fmax_int(0, NN - u - 1); i <= NN - 1; i++) {
        double lgamma_i1 = lgammafn(i + 1);
        double lgamma_NN_i = lgammafn(NN - i);
        kk1 = exp(lgamma_NN_u_1 - (lgamma_i1 + lgamma_NN_i + lgamma_NN_i + lgamma_u_NN_2 + i));
        double term1 = R_pow(p11, i + 1);
        double term2 = R_pow(1 + p11 - (p01 + p10), u - NN + 1 + i);
        double term3 = R_pow(p01 - p11, NN - 1 - i);
        double term4 = R_pow(p10 - p11, NN - 1 - i);
        dens += kk1 * term1 * term2 * term3 * term4;
    }
    return dens;
}

double aux_biv_binomneg(int NN, int u, int v, double x, double y, double p11) {
    int a = 0, i = 0;
    double kk1 = 0.0, kk2 = 0.0, dens1 = 0.0, dens2 = 0.0;
    
    // Pre-calcoliamo i termini di lgammafn per evitare chiamate ripetute
    double lgamma_NN_u = lgammafn(NN + u);
    double lgamma_v_u = lgammafn(v - u);
    for (a = fmax_int(0, u - v + NN - 1); a <= NN - 2; a++) {
        double lgamma_i1 = lgammafn(1); // lgammafn(i+1) quando i=0
        double lgamma_NN_i = lgammafn(NN);
        for (i = fmax_int(0, a - u); i <= fmin_int(a, NN - 1); i++) {
            lgamma_i1 = lgammafn(i + 1);  // lgammafn(i+1)
            lgamma_NN_i = lgammafn(NN - i); // lgammafn(NN-i)
            double lgamma_a_i_1 = lgammafn(a - i + 1); // lgammafn(a-i+1)
            double lgamma_u_a_i_1 = lgammafn(u - a + i + 1); // lgammafn(u-a+i+1)
            kk1 = exp(lgamma_NN_u - (lgamma_i1 + lgamma_NN_i + lgamma_a_i_1 + lgamma_u_a_i_1));
            double lgamma_v_a_NN_u = lgammafn(v + a - NN - u + 2);
            double lgamma_NN_a_1 = lgammafn(NN - a - 1);
            kk2 = exp(lgamma_v_u - (lgamma_v_a_NN_u + lgamma_NN_a_1));
            
            dens1 += kk1 * kk2 * R_pow(p11, i + 1) * R_pow(1 + p11 - (x + y), u - a + i) *
                     R_pow(x - p11, NN - i - 1) * R_pow(y - p11, a - i) * R_pow(1 - y, v - u - NN + a + 1) * R_pow(y, NN - a - 1);
        }
    }

    for (a = fmax_int(0, u - v + NN); a <= NN - 1; a++) {
        double lgamma_i1 = lgammafn(1); // lgammafn(i+1) quando i=0
        double lgamma_NN_1_i = lgammafn(NN - 1); // lgammafn(NN-1-i) quando i=0
        for (i = fmax_int(0, a - u); i <= fmin_int(a, NN - 1); i++) {
            // Pre-calcoliamo i valori di lgammafn all'interno dei cicli
            lgamma_i1 = lgammafn(i + 1);  // lgammafn(i+1)
            lgamma_NN_1_i = lgammafn(NN - 1 - i + 1); // lgammafn(NN-1-i)
            double lgamma_a_i_1 = lgammafn(a - i + 1); // lgammafn(a-i+1)
            double lgamma_u_a_i_1 = lgammafn(u - a + i + 1); // lgammafn(u-a+i+1)
            kk1 = exp(lgamma_NN_u - (lgamma_i1 + lgamma_NN_1_i + lgamma_a_i_1 + lgamma_u_a_i_1));
            double lgamma_v_a_NN_u = lgammafn(v + a - NN - u + 1);
            double lgamma_NN_a = lgammafn(NN - a);
            kk2 = exp(lgamma_v_u - (lgamma_v_a_NN_u + lgamma_NN_a));
            
            dens2 += kk1 * kk2 * R_pow(p11, i) * R_pow(1 + p11 - (x + y), u - a + i) *
                     R_pow(x - p11, NN - i) * R_pow(y - p11, a - i) * R_pow(1 - y, v - u - NN + a) * R_pow(y, NN - a);
        }
    }
    return dens1 + dens2;
}


void biv_binomneg_call(int *NN,int *u, int *v, double *p01, double *p10,double *p11,double *res)
{
    *res = biv_binomneg (*NN,*u,*v,*p01, *p10,*p11);
}




double biv_binomneg (int NN, int u, int v, double p01,double p10,double p11)
{
double dens=0.0;
if(u<v)    dens=aux_biv_binomneg(NN,u,v,p01,p10,p11);

if(u==v)      dens=aux_biv_binomneg_simple(NN,v,p10,p01,p11);

if(u>v)    dens=aux_biv_binomneg(NN,v,u,p10,p01,p11);
return(dens);
}
/**********************************************************************/
double bin_aux(int a,int NN,int u,int v,double p1, double p2,double p11)
{
  double kk,dens;
kk=exp(lgammafn(NN+1)-(lgammafn(a+1)+lgammafn(u-a+1)+lgammafn(v-a+1)+lgammafn(NN-u-v+a+1)));
dens=kk*(R_pow(p11,a)*R_pow(p1-p11,u-a)*R_pow(p2-p11,v-a)*R_pow(1+p11-(p1+p2),NN-u-v+a));
return(dens);
}


void biv_binom_call(int *NN,int *u, int *v, double *p01, double *p10,double *p11,double *res)
{
    *res = biv_binom (*NN,*u,*v,*p01, *p10,*p11);
}


double biv_binom(int NN, int u, int v, double p01, double p10, double p11) {   
    int a;
    double kk = 0.0, dens = 0.0;
    double lgamma_NN_1 = lgammafn(NN + 1);
    for (a = fmax_int(0, u + v - NN); a <= fmin_int(u, v); a++) {
        double lgamma_a1 = lgammafn(a + 1);
        double lgamma_u_a1 = lgammafn(u - a + 1);
        double lgamma_v_a1 = lgammafn(v - a + 1);
        double lgamma_NN_uv_a1 = lgammafn(NN - u - v + a + 1);
        kk = exp(lgamma_NN_1 - (lgamma_a1 + lgamma_u_a1 + lgamma_v_a1 + lgamma_NN_uv_a1));
        double term1 = R_pow(p11, a);
        double term2 = R_pow(p01 - p11, u - a);
        double term3 = R_pow(p10 - p11, v - a);
        double term4 = R_pow(1 + p11 - (p01 + p10), NN - u - v + a);
        dens += kk * term1 * term2 * term3 * term4;
    }

    return dens;
}



double aux_biv_binom(int n1,int n2, int u, int v,double p01,double p10,double p11)
{

int k;
double dens=0.0;
int N=n1-n2;
for(k=0;k<=N;k++)
dens+= exp( lgammafn(N+1)-lgammafn(k+1)-lgammafn(N-k+1)+k*log(p01)+(N-k)*log1p(-p01)+log(biv_binom(n2, u-k,  v, p01,p10,p11)));
 return(dens);   
}

double biv_binom222(int n1,int n2, int u, int v, double p01,double p10,double p11)
{
  
 double res=0.0;
 if(n1> n2)  res=aux_biv_binom(n1,n2,u,v,p01,p10,p11);
 if(n2> n1)  res=aux_biv_binom(n2,n1,v,u,p10,p01,p11);  
 if(n1==n2)  res=biv_binom (n1,u,v,p01,p10,p11);
return(res);
}

/// biv binomial type II


double  biv_binom2 (int NN_i,int NN_j, int k, int u, int v, double p01,double p10,double p11)
{

int a,i,j;
double const1=0.0,const2=0.0,const3=0.0,dens=0.0,dens1;
double P10=p01-p11;
double P01=p10-p11;
double P00=1+p11-(p01+p10);
double P11=p11;

for(i=0;i<=fmin_int(NN_i-k,u);i++){
   for(j=0;j<=fmin_int(NN_j-k,v);j++){
     
for(a=fmax_int(0,u+v-k-i-j);a<=fmin_int(u-i,v-j);a++){         
       const1=exp(lgammafn(k+1)-(lgammafn(a+1)+lgammafn(u-i-a+1)+lgammafn(v-j-a+1)+lgammafn(k-u-v+i+j+a+1)));
        dens1=const1*R_pow(P11,a)*R_pow(P00,k-u-v+i+j+a)*
             R_pow(P10,u-i-a)*R_pow(P01,v-j-a);
   
       const2=exp(lgammafn(NN_i-k+1)-(lgammafn(NN_i-k-i+1)+lgammafn(i+1)));
       const3=exp(lgammafn(NN_j-k+1)-(lgammafn(NN_j-k-j+1)+lgammafn(j+1)));

       dens+=dens1*const2*const3*
             R_pow(P11+P10,i) * R_pow(P11+P01,j) *
             R_pow(P00+P01,NN_i-k-i) * R_pow(P00+P10,NN_j-k-j);
  }}}
    return(dens);
}



// bivariate pois-binomial 
double biv_poisbin (int NN, int u, int  v, double p01,double p10,double p11)
{
    int a=0,a_u=0;
    double dens=0.0, kk=0.0;
    a_u=fmin(u,v);
    dens=0;
    for(a=0;a<=a_u;a++){

    kk=exp(-(lgammafn(a+1)+lgammafn(u-a+1)+lgammafn(v-a+1)));
        dens+=kk*(R_pow(NN*p11,a)*R_pow(NN*(p01-p11),u-a)*R_pow(NN*(p10-p11),v-a));
    }
    return(exp(-NN*(p01+p10-p11))*dens);
}



// bivariate pois-bineg 
double biv_poisbinneg (int NN, int u, int v, double p01,double p10,double p11)
{
    int a=0,a_u=0;
    double dens=0.0,pp=0.0,bb=0.0, kk=0.0;
    a_u=fmin_int(u,v);
    dens=0;
    pp=p11*p01*p10;
    bb=(1-p11)/p11;
    for(a=0;a<=a_u;a++){

       kk=exp(-(lgammafn(a+1)+lgammafn(u-a+1)+lgammafn(v-a+1)));
       dens+=kk*(R_pow(NN*((p11*(p01+p10)-p01*p10*(p11+1))/pp),a)*R_pow(NN*((p10-p11)/(p10*p11)),u-a)*R_pow(NN*bb,v-a));
    }
    return(exp(-NN*bb)*dens);
}



// marginal binomial
double marg_binom(int n,double x,double p)
{
 double pr=0.0;double kk;
 kk=fac(n,n-x+1)/fac(x,1);
 pr=kk*R_pow(p,x)*R_pow(1-p,n-x);
return(pr);
}
//marginal geom
double marg_geom(int x,double p)
{
 double pr=0.0;
 pr=p*R_pow(1-p,x);
    return(pr);
}
//marginal binom neg
double marg_binomneg(int n,int x,double p)
{
 double pr=0.0;
 pr=(fac(x+n-1,1)/(fac(n-1,1)*fac(x,1)))*R_pow(p,n)*R_pow(1-p,x);
    return(pr);
}
//marginal poisson
double marg_pois(int n,double x,double p)
{
 double pr=0.0,lambda=0.0;
 lambda=n*p;
 pr=R_pow(lambda,x)*exp(-lambda)/fac(x,1);
    return(pr);
}
       
/* cdf inflated poisson*/
double ppoisinflated(double z, double m, double p)
{
 double res;
 double ans = ppois(z, m,1,0);
 if(z < 0) res=0;
 else     res=p + (1 - p) * ans;
 return(res);
}
/* cdf inflated binom neg*/
double pbneginflated(double z, int N,double m, double p)
{
 double res;   
 double ans = pnbinom(z,N,m,1,0);
 if(z < 0) res=0;
 else     res=p + (1 - p) * ans;
 return(res);
}


double marg_p(double categ_0,double psm,int *model,int n)
{
    double res=0.0;
    if(model[0]==2||model[0]==11) res=marg_binom (n-1,categ_0,psm);
    if(model[0]==14)           res=marg_geom (categ_0,psm);
    if(model[0]==15)           res=marg_pois (n-1,categ_0,psm);
    return(res);
}


/*********** ratio of gamma function type 1********************/
double aprox(int k,double a, double b) 
{
return(R_pow(k,a-b)*(1+ (a-b)*(a+b-1)/(2*k)));
}
/************ ratio of gamma function type 2 *******************/
double aprox1(double a) {
  double res=0.0;
    res=R_pow(R_pow(a,5)+(5*R_pow(a,4)/4)+(25*R_pow(a,3)/32)+(35*R_pow(a,2)/128)+(75*a/2048)+3/8192,1/10);
    return(res);
 }





/**********************************************************/

/* Logarithm of Gamma function 
double lgam(double x)
{
    int sign;
    return lgam_sgn(x, &sign);
}
*/
double lgam_sgn(double x, int *sign)
{
    double p, q, u, w, z;
    int i;
    
    *sign = 1;
    
    if (!R_FINITE(x))
        return x;
    
    if (x < -34.0) {
        q = -x;
        w = lgam_sgn(q, sign);
        p = floor(q);
        if (p == q) {
        lgsing:
            //mtherr("lgam", SING);
           // printf("SIGN\n");
            return (NPY_INFINITY);
        }
        i = p;
        if ((i & 1) == 0)
            *sign = -1;
        else
            *sign = 1;
        z = q - p;
        if (z > 0.5) {
            p += 1.0;
            z = p - q;
        }
        z = q * sin(NPY_PI * z);
        if (z == 0.0)
            goto lgsing;

        z = LOGPI - log(z) - w;
        return (z);
    }
    
    if (x < 13.0) {
        z = 1.0;
        p = 0.0;
        u = x;
        while (u >= 3.0) {
            p -= 1.0;
            u = x + p;
            z *= u;
        }
        while (u < 2.0) {
            if (u == 0.0)
                goto lgsing;
            z /= u;
            p += 1.0;
            u = x + p;
        }
        if (z < 0.0) {
            *sign = -1;
            z = -z;
        }
        else
            *sign = 1;
        if (u == 2.0)
            return (log(z));
        p -= 2.0;
        x = x + p;
        p = x * polevl(x, B, 5) / p1evl(x, C, 6);
        return (log(z) + p);
    }
    
    if (x > MAXLGM) {
        return (*sign * NPY_INFINITY);
    }
    
    q = (x - 0.5) * log(x) - x + LS2PI;
    if (x > 1.0e8)
        return (q);
    
    p = 1.0 / (x * x);
    if (x >= 1000.0)
        q += ((7.9365079365079365079365e-4 * p
               - 2.7777777777777777777778e-3) * p
              + 0.0833333333333333333333) / x;
    else
        q += polevl(p, A, 4) / x;
    return (q);
}












/***********************************************************/
/*********** bivariate T distribution*********************/
double biv_T(double rho, double zi, double zj, double nuu, double nugget)
{
    int k = 0;
    double nu = 1 / nuu;
    double res0 = 0.0, RR = 0.0, pp1 = 0.0, pp2 = 0.0;
    double r2 = rho * rho;
    double bb1, bb2;
    double x = zi, y = zj;
    double cc = (nu + 1) / 2;
    double nu2 = nu / 2;
    double rho1 = rho * (1 - nugget);
    double rn = rho1 * rho1;
    double rho2 = pow1p(-r2, -nu2 - 1); // pow1p(x, y) computes (1 + x)^y
    double rho12 = pow1p(-rho1 * rho1, -nu - 0.5);
    double x1 = (x * x * (1 - r2) + nu * (1 - rn));
    double y1 = (y * y * (1 - r2) + nu * (1 - rn));
    double C, B;
    double b1 = exp(nu * log(nu) - cc * log(x1 * y1) + 2 * lgammafn(cc));
    double c1 = exp(log(M_PI) + 2 * lgammafn(nu / 2) + log(rho12) + log(rho2));
    double b2 = rho1 * x * y * R_pow(nu, nu + 2) * R_pow(x1 * y1, -nu2 - 1);
    double c2 = 2 * M_PI * pow1p(-rn, -nu - 0.5) * pow1p(-r2, -nu2 - 2);
    double a1 = 0, a2 = 0;
    double aux = R_pow(rho1 * x * y, 2) * pow1p(-r2, 2) / (x1 * y1);
    double aux1 = R_pow(rho * nu, 2) * pow1p(-rn, 2) / (x1 * y1);
    if (rho > DEPSILON)
    {
        while (k <= 3000)
        {
            // pp1 = hypergeo(cc + k, cc + k, 0.5, aux);
            pp1 = (0.5 - 2 * (cc + k)) * log1p(-aux) + log(hypergeo(0.5 - (cc + k), 0.5 - (cc + k), 0.5, aux)); // Euler
            bb1 = pp1 + k * log(aux1) + 2 * (lgammafn(cc + k) - lgammafn(cc)) - lgammafn(k + 1) - lgammafn(nu2 + k) + lgammafn(nu2);
            a1 += exp(bb1);
            // pp2 = hypergeo(nu2 + 1 + k, nu2 + 1 + k, 1.5, aux);
            pp2 = (1.5 - 2 * (nu2 + 1 + k)) * log1p(-aux) + log(hypergeo(1.5 - (nu2 + 1 + k), 1.5 - (nu2 + 1 + k), 1.5, aux)); // Euler
            bb2 = pp2 + k * log(aux1) + 2 * log((1 + k / nu2)) + lgammafn(nu2 + k) - lgammafn(k + 1) - lgammafn(nu2);
            a2 += exp(bb2);
            RR = (b1 / c1) * a1 + (b2 / c2) * a2;
            if (fabs(RR - res0) < 1e-10 || !R_FINITE(RR)) { break; }
            else { res0 = RR; }
            k++;
        }
        if (!R_finite(RR)) RR = 1e-320;
    }

    if (rho <= DEPSILON)
    {
        C = lgammafn(cc) + log(R_pow((1 + x * x / nu), -cc)) - log(sqrt(M_PI * nu)) - lgammafn(nu / 2);
        B = lgammafn(cc) + log(R_pow((1 + y * y / nu), -cc)) - log(sqrt(M_PI * nu)) - lgammafn(nu / 2);
        RR = exp(B) * exp(C);
    }

    return RR;
}

/*********** Appell F4 function ********/



void appellF4_call(double *a,double *b,double *c,double *d,double *x,double *y, double *res)
{
    *res = appellF4(*a,*b,*c,*d,*x,*y);
}



double appellF4(double a,double b,double c,double d,double x,double y)
{
double RR=0.0,bb=0.0;int k=0;
  while( k<=900 )
    {
    bb=exp(k*log(y)+(lgammafn(a+k)+lgammafn(b+k)+lgammafn(d))
               -(lgammafn(a)+lgammafn(b)+lgammafn(d+k)+lgammafn(k+1))
               +(c-(a+k)-(b+k))*log1p(-x)+log(hypergeo(c-a-k,c-b-k,c,x))); //euler
              // +log(hypergeo(a+k,b+k,c,x));
    if((fabs(bb)<1e-7||!R_FINITE(bb))  ) {break;}
        RR=RR+bb;
        k++;
    }
    if(!R_finite(RR)) RR=1e-320;
return(RR);
}




double appellF42211(double x,double y)
{
double RR=0.0,bb=0.0;int k=0;
  while( k<=5000 )
    {
    bb=exp(k*log(y)+2*(lgammafn(2+k)-lgammafn(1+k))
       +(1-2*(2+k))*log1p(-x)+log(hypergeo(-1-k,-1-k,1,x))); 
    if((fabs(bb)<1e-10||!R_FINITE(bb))  ) {break;}
        RR=RR+bb;
        k++;
    }
    if(!R_finite(RR)) RR=1e-320;
return(RR);
}





double appellF4_mod(double nu,double rho,double x,double y,double nugget)
{
  double x2,y2,arg,arg1,pp1,pp2,app;
  double rho1=rho*(1-nugget);
  double rho22=R_pow(1-rho*rho,-nu/2-1);
  double rho12=R_pow(1-rho1*rho1,-nu-0.5);

x2=(x*x*(1-rho*rho)+nu*(1-rho1*rho1));
y2=(y*y*(1-rho*rho)+nu*(1-rho1*rho1));

arg=(nu+1)/2;  arg1=nu/2;
pp1=exp(nu*log(nu)-arg*log(x2*y2)+2*lgammafn(arg));
pp2=exp(log(M_PI)+2*lgammafn(arg1)+log(rho12)+log(rho22));
app=appellF4(arg,arg,0.5,arg1,R_pow(rho1*x*y*(1-rho*rho),2)/(x2*y2), R_pow(rho*nu*(1-rho1*rho1),2)/(x2*y2));
return(4*pp1*app/pp2);
    //return(50);


}

/*********** bivariate two piece-T distribution********************/ 

double biv_two_pieceT(double rho,double zi,double zj,double sill,double nuu,double eta,
             double p11,double mui,double muj,double nugget)
{
double res=0.0;
double nu=1/nuu;
double etamas=1+eta;
double etamos=1-eta;
double zistd=(zi-mui)/sqrt(sill);
double zjstd=(zj-muj)/sqrt(sill);

    if(rho>DBL_EPSILON){
    if(zistd>=0&&zjstd>=0)
{res=          (p11/R_pow(etamos,2))*appellF4_mod(nu,rho,zistd/etamos,zjstd/etamos,nugget);}
    if(zistd>=0&&zjstd<0)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*appellF4_mod(nu,rho,zistd/etamos,-zjstd/etamas,nugget);}
      if(zistd<0&&zjstd>=0)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*appellF4_mod(nu,rho,-zistd/etamas,zjstd/etamos,nugget);}
    if(zistd<0&&zjstd<0)
{res=    ((p11+eta)/R_pow(etamas,2))*appellF4_mod(nu,rho,-zistd/etamas,-zjstd/etamas,nugget);}

}
    if(rho<DBL_EPSILON){
    
    if(zi>=mui)
         {res=0.5*2*exp((nu/2)*log(nu)+lgammafn((nu+1)/2)-((nu+1)/2)*log(R_pow(zistd/etamos,2)+nu)-0.5*log(M_PI)-lgammafn(nu/2));}
         if(zj<muj)
         {res=0.5*2*exp((nu/2)*log(nu)+lgammafn((nu+1)/2)-((nu+1)/2)*log(R_pow(zjstd/etamas,2)+nu)-0.5*log(M_PI)-lgammafn(nu/2));}
      }
return(res/sill);
    
}





/***** bivariate half gaussian ****/     
double biv_half_Gauss(double rho,double zi,double zj)
{
double kk=0, dens=0,a=0,b=0,rho2=1-rho*rho;
  kk=(M_PI)*sqrt(rho2);
  a=exp(- (1/(2*(rho2)))*(R_pow(zi,2)+R_pow(zj,2)-2*rho*zi*zj));
  b=exp(- (1/(2*(rho2)))*(R_pow(zi,2)+R_pow(zj,2)+2*rho*zi*zj));
  dens=(a+b)/kk;
  return(dens);
}
/***** bivariate two piece gaussian *****/
double biv_two_pieceGaussian(double rho, double zi, double zj, double sill, double eta,
                              double p11, double mui, double muj)
{
    double res = 0.0;
    double etamas = 1 + eta;
    double etamos = 1 - eta;
    double zistd = (zi - mui) / sqrt(sill);
    double zjstd = (zj - muj) / sqrt(sill);
    
    double z1 = (zistd >= 0) ? zistd / etamos : -zistd / etamas;
    double z2 = (zjstd >= 0) ? zjstd / etamos : -zjstd / etamas;
    
    if (zistd >= 0 && zjstd >= 0) {
        res = (p11 / R_pow(etamos, 2)) * biv_half_Gauss(rho, z1, z2);
    } else if (zistd >= 0 && zjstd < 0) {
        res = ((1 - eta - 2 * p11) / (2 * (1 - eta * eta))) * biv_half_Gauss(rho, z1, -z2);
    } else if (zistd < 0 && zjstd >= 0) {
        res = ((1 - eta - 2 * p11) / (2 * (1 - eta * eta))) * biv_half_Gauss(rho, -z1, z2);
    } else if (zistd < 0 && zjstd < 0) {
        res = ((p11 + eta) / R_pow(etamas, 2)) * biv_half_Gauss(rho, -z1, -z2);
    }
    
    return res / sill;
}





double biv_beta(double rho,double zi,double zj,double ai,double aj,double shape1,double shape2,double min,double max)
{
  double ki=0.0,kj=0.0,p1=0.0,p2=0.0,rho2=0.0,res=0.0,p3;
  double dd=max-min;
  zi=(zi-min)/dd;zj=(zj-min)/dd;
 ki=1-zi; kj=1-zj;
   double aa=0.5*(shape1+shape2);
if(rho) {
  rho2=rho*rho;

    p1=R_pow(zi*zj,shape1/2-1)*R_pow(ki*kj,shape2/2-1);
    p3=exp(2*lgammafn(aa)-(2*lgammafn(shape1/2)+2*lgammafn(shape2/2)-aa*log1p(-rho2)));
    p2= appellF4(aa,aa,shape1/2,shape2/2,rho2*zi*zj,rho2*ki*kj);
  res=p1*p2*p3;
} else  {p1=R_pow(zi,shape1/2-1)*R_pow(ki,shape2/2-1)*exp(lgammafn(aa)-lgammafn(shape1/2)-lgammafn(shape2/2));
         p2=R_pow(zj,shape1/2-1)*R_pow(kj,shape2/2-1)*exp(lgammafn(aa)-lgammafn(shape1/2)-lgammafn(shape2/2));
         res=p1*p2;}
return(res/R_pow(dd,2));
}







double biv_Unif(double rho,double ui,double uj)
{
  double res,rho2;
if(fabs(rho)>1e-8){
  rho2=rho*rho;
  res= 4*(1-rho2)*appellF4(2,2,1,1,rho2*ui*uj,rho2*(1-ui)*(1-uj));}
 else  {res= 1.0;}
return(res);
}



double biv_Kumara(double rho,double zi,double zj,double ai,double aj,double shape1,double shape2,double min,double  max)
{
  double xx=0.0,yy=0.0,ki=0.0,kj=0.0,p1=0.0,p2=0.0,rho2=0.0,res=0.0;
    double dd=(max-min);
  zi=(zi-min)/dd;zj=(zj-min)/dd;

 ki=1-R_pow(zi,shape2); kj=1-R_pow(zj,shape2);
if(rho) {
  rho2=rho*rho;
  xx=rho2*R_pow(ki*kj,shape1);
  yy=rho2*(1-R_pow(ki,shape1))*(1-R_pow(kj,shape1));
  p1=R_pow(shape1*shape2,2)*R_pow(zi*zj,shape2-1)*R_pow(ki*kj,shape1-1)*R_pow(1-rho2,2);
  p2= appellF4(2,2,1,1,xx,yy);
  res=p1*p2;}
 else  {p1=shape1*shape2*R_pow(zi,shape2-1)* R_pow(ki,shape1-1);
         p2=shape1*shape2*R_pow(zj,shape2-1)* R_pow(kj,shape1-1);
         res=p1*p2;}
return(res/R_pow(dd,2));
}



double biv_Kumara2(double rho,double zi,double zj,double ai,double aj,double shape1,double shape2,double min,double  max)
{
  double xx=0.0,yy=0.0,ki=0.0,kj=0.0,p1=0.0,p2=0.0,rho2=0.0,res=0.0;
  double mi=1/(1+exp(-ai)), mj=1/(1+exp(-aj));
  double dd=(max-min), shapei=log(0.5)/log1p(-R_pow(mi,shape2)),shapej=log(0.5)/log1p(R_pow(mj,shape2));
  zi=(zi-min)/dd;zj=(zj-min)/dd;

 ki=1-R_pow(zi,shape2); kj=1-R_pow(zj,shape2);
if(rho) {
  rho2=rho*rho;
  xx=rho2*R_pow(ki,shapei)*R_pow(kj,shapej);
  yy=rho2*(1-R_pow(ki,shapei))*(1-R_pow(kj,shapej));
  p1=R_pow(sqrt(shapei*shapej)*shape2,2)*R_pow(zi*zj,shape2-1)*R_pow(ki,shapei-1)*R_pow(kj,shapej-1)*R_pow(1-rho2,2);
  p2= appellF4(2,2,1,1,xx,yy);
  res=p1*p2;}
else  {p1=shapei*shape2*R_pow(zi,shape2-1)*R_pow(ki,shapei-1);
       p2=shapej*shape2*R_pow(zj,shape2-1)*R_pow(kj,shapej-1);
         res=p1*p2;}
return(res/R_pow(dd,2));
}

/***********************************************************************************/
/************ functions for binomial or negative binomial  two piece *********/
/***********************************************************************************/

double pbnorm22(double lim1,double lim2,double corr)
{
    double  lowe[2]  = {0,0}, uppe[2] = {lim1,lim2}, corre[1] = {corr};
    double  value;
    int     infin[2] = {0,0};
    value            = F77_CALL(bvnmvn)(lowe,uppe,infin,corre); 
    return(value);
}

/**********************************************************************/
double pblogi22(double lim1,double lim2,double corr)
{

double value=0.0,sum1=0.0,sum2=0.0,term=0.0,term2=0.0,kk=0,bb=0.0;
double corr21=1-R_pow(corr,2);
int m=0,n=0,p1,p2,p3;
double e1=exp(lim1);double e2=exp(lim2);
while(n<=400){
bb=exp(2*n*log(corr)+n*(lim1+lim2)- 2*log(n+1));
sum1=0.0;m=0;
p1=1+n;p3=p1+1;
while(m<=200)
{ 
      p2=2+n+m;
      term=exp(2*m*log(corr)+log(hypergeo(p1,p2,p3,-e1))+
                           log(hypergeo(p1,p2,p3,-e2))-2*lbeta(p1,m+1));
     sum1=sum1+ term;
     if(term<1e-7) {break;}
     m=m+1;      

}
term2=bb*sum1; 
sum2=sum2+term2;
if(term2<1e-7) {break;}
n=n+1;
}
kk=exp(2*log(corr21)+lim1+lim2);
value= sum2*kk;
 return(value);
}

// cdf bivariate half-normal distribution
double pbhalf_gauss(double zi,double zj,double rho,double nugget)
{
  double dens = 0;
  dens = pbnorm22(zi,zj,(1-nugget)*rho) + pbnorm22(-zi,-zj,(1-nugget)*rho) -
              pbnorm22(-zi,zj,(1-nugget)*rho) - pbnorm22(zi,-zj,(1-nugget)*rho);
  return(dens); 
}
/****** cdf univariate two-piece gaussian distribution *****/
double pnorm_two_piece(double x, double eta)
{
    double cdf = 0;
    if (x <=  0) cdf = (1 + eta)*pnorm(x/(1 + eta),0,1,1,0);
    if (x > 0)   cdf = eta + (1 - eta)*pnorm(x/(1 - eta),0,1,1,0);
  return(cdf);
}
//***********************************************//
// cdf univariate half-gaussian distribution
double phalf_gauss (double z){
  double dens = 0;
  dens = 2 * pnorm(z,0,1,1,0) - 1;
  return(dens);
}

// cdf bivariate two_piece gaussian distribution
double pbnorm_two_piece(int *cormod, double h, double u, 
    double xi, double xj, double nugget, double var,double eta,double *par)
{
    double p11,g_eta=0.0,q_g_eta=0.0,dens=0.0,corr=0.0;
    double etamos,etamas;
    g_eta   = (1 - eta)/2;
    q_g_eta = qnorm(g_eta,0,1,1,0);
    corr    = CorFct(cormod,h,u,par,0,0);
    p11     = pbnorm22(q_g_eta,q_g_eta,(1-nugget)*corr);
    etamas = 1+eta;
    etamos = 1-eta;
    if((xi  <= 0) & (xj <= 0)) dens = (1 + pbhalf_gauss(-xi/etamas,-xj/etamas,corr,nugget) - phalf_gauss(-xi/etamas) - phalf_gauss(-xj/etamas)) * ( 1 + p11 - 2*pnorm(q_g_eta,0,1,1,0) );
    if((xi  > 0)  & (xj <= 0)) dens = (1 + pbhalf_gauss(0,-xj/etamas,corr,nugget) - phalf_gauss(-xj/etamas)) * ( 1 + p11 - 2*pnorm(q_g_eta,0,1,1,0) ) + 
      (phalf_gauss(xi/etamos) - pbhalf_gauss(xi/etamos,-xj/etamas,corr,nugget)) * (pnorm(q_g_eta,0,1,1,0) - p11);
    if((xi  <= 0) & (xj >  0)) dens =  (1 + pbhalf_gauss(-xi/etamas,0,corr,nugget) - phalf_gauss(-xi/etamas)) * ( 1 + p11 - 2*pnorm(q_g_eta,0,1,1,0) ) + 
      (phalf_gauss(xj/etamos) - pbhalf_gauss(-xi/etamas,xj/etamos,corr,nugget)) * (pnorm(q_g_eta,0,1,1,0) - p11);
    if((xi  > 0)  & (xj >  0)) dens =  ( 1 + p11 - 2*pnorm(q_g_eta,0,1,1,0) ) + (phalf_gauss(xi/etamos) - pbhalf_gauss(xi/etamos,0,corr,nugget))*(pnorm(q_g_eta,0,1,1,0)-p11) + 
     (phalf_gauss(xj/etamos) - pbhalf_gauss(0,xj/etamos,corr,nugget))* (pnorm(q_g_eta,0,1,1,0)-p11) + pbhalf_gauss(xi/etamos,xj/etamos,corr,nugget)* p11;
    return(dens);
}
/***********************************************************************************/
/***********************************************************************************/
/********** bivariate log-logistic **********/
double biv_LogLogistic(double corr,double zi,double zj,double mui, double muj, double shape)
{
    double c=gammafn(1+1/shape)*gammafn(1-1/shape);
    double A=0.0,res=0.0,B=0.0,C=0.0;
    double ci=exp(mui);double cj=exp(muj);
    double ki=R_pow(c*zi/ci,shape)+1;
    double kj=R_pow(c*zj/cj,shape)+1;
    double rho2=R_pow(corr,2);
    double kij=ki*kj;
   // if(corr)   {
        A=(R_pow(c*shape,2)/(ci*cj))*R_pow((c*c*zi*zj)/(ci*cj),shape-1)*R_pow(1-rho2,2);
        B=R_pow(kij,-2);
        C=appellF4(2,2,1,1,
        (rho2*R_pow(c*c*zi*zj,shape))/(R_pow(ci*cj,shape)*kij),
        rho2/(kij));
        res=A*B*C;
    //}else{    B=(c*shape/ci)*R_pow((c*zi/ci),shape-1)*R_pow(ki,-2);  C=(c*shape/cj)*R_pow((c*zj/cj),shape-1)*R_pow(kj,-2);
    //    res=B*C;}
    return(res);
}

/********** bivariate logistic **********/
double biv_Logistic(double corr,double zi,double zj,double mui, double muj, double sill)
{
    double a=0.0,A=0.0,res=0.0,B=0.0,C=0.0;
    double ci=mui;double cj=muj;
    double ki=exp((zi-ci)/sqrt(sill));
    double kj=exp((zj-cj)/sqrt(sill));
    double rho2=R_pow(corr,2);
    //if(corr)   {
        a=1-rho2;
        A=(ki*kj)/(R_pow(a,-2)*sill);
        B=R_pow((ki+1)*(kj+1),-2);
        C=appellF4(2,2,1,1,
        (rho2*ki*kj)/((ki+1)*(kj+1)),
        rho2/((ki+1)*(kj+1)));
        res=A*B*C;
    //} else{ B=ki*R_pow((ki+1),-2)/sqrt(sill);C=kj*R_pow((kj+1),-2)/sqrt(sill);res=B*C;}
    return(res);
}




//********** functions for bivariate tukey h *******************//

// compute lambert w function
double LambertW(double z) {
  const double eps = 4.0e-16;               // Precisione del calcolo
  const double em1 = 0.3678794411714423;    // e^(-1)
  double p, e, t, w;
  
  // Gestione dei casi speciali
  if (z < -em1 || isinf(z) || isnan(z)) {
    return R_NaN;  // Meglio restituire NAN piuttosto che 0.0 per input invalidi
  }
  
  if (z == 0.0) return 0.0;
  
  // Serie di approssimazione vicino a -e^(-1)
  if (z < -em1 + 1e-4) {
    double q = z + em1;
    double r = sqrt(q);
    double q2 = q * q;
    double q3 = q2 * q;
    return -1.0
           + 2.331643981597124203 * r
           - 1.812187885639363490 * q
           + 1.936631114492359755 * r * q
           - 2.353551201881614516 * q2
           + 3.066858901050631912 * r * q2
           - 4.175335600258177138 * q3
           + 5.858023729874774148 * r * q3
           - 8.401032217523977370 * q3 * q;
  }
  
  // Calcolo del valore iniziale di w per l'iterazione
  if (z < 1.0) {
    // Espansione in serie vicino a 0
    p = sqrt(2.0 * (M_E * z + 1.0)); // Uso di M_E da math.h
    w = -1.0 + p * (1.0 + p * (-1.0/3.0 + p * (1.0/72.0)));  // Migliorata precisione
  } else {
    // Approssimazione asintotica
    w = log(z);
    // Miglioramento per z > 3.0
    if (z > 3.0) {
      w -= log(w);
      // Ulteriore miglioramento per z molto grandi
      if (z > 30.0) {
        double l1 = log(w);
        w -= l1 * (1.0 - l1/(2.0 * w));
      }
    }
  }
  
  int i;
  for (i = 0; i < 6; i++) {  
    e = exp(w);
    t = w * e - z;
    p = w + 1.0;
    
    double f = e * p; 
    double f1 = (p + 1.0) * t / (2.0 * p);
    double d = t / (f - f1);
    
    w -= d;
    if (fabs(d) < eps * (1.0 + fabs(w))) {
      return w;
    }
  }
  return w;
}

// pdf bivariate gaussian distribution
double dbnorm(double x_i,double x_j,double mean_i,double mean_j,double sill,double corr)
{
    double  fraq = 1.0,dens = 0.0,aux1 = 1.0,z1 = 1.0,z2 = 1.0,z3 = 1.0,z = 1.0;
    fraq = 2*M_PI*sill*sqrt(1-corr*corr);
    z1   = (x_i - mean_i)*(x_i - mean_i);
    z2   = (x_j - mean_j)*(x_j - mean_j);
    z3   = 2*corr*(x_i - mean_i)*(x_j - mean_j);
    z    = (z1 + z2 - z3)/sill;
    aux1 = 2*(1-corr*corr); 
    dens = (1/fraq)*exp(-z/aux1);
    return(dens);
}



// pdf bivariate tukey h random field 
double biv_tukey_h(double corr,double data_i, double data_j, double mean_i, double mean_j, double tail, double sill)
{
  double dens = 0.0,x_i = 0.0,x_j = 0.0,est_mean_i = 0.0,est_mean_j = 0.0;
  double est_mean_ij = 1.0,extra = 1.0;

  est_mean_i = (data_i - mean_i)/sqrt(sill);
  est_mean_j = (data_j - mean_j)/sqrt(sill);
  
  x_i = inverse_lamb(est_mean_i,tail);
  x_j = inverse_lamb(est_mean_j,tail);

  est_mean_ij = 1/(est_mean_i*est_mean_j);
  extra       = 1/( (1 + LambertW(tail*est_mean_i*est_mean_i))*(1 + LambertW(tail*est_mean_j*est_mean_j)));
  dens = dbnorm(x_i,x_j,0,0,1,corr)*x_i*x_j * est_mean_ij * extra/sill;

  
  if((x_i==0.0)&&(x_j!=0.0))  dens = dbnorm(x_i,x_j,0,0,1,corr)*x_j/(est_mean_j*(1 + LambertW(tail*est_mean_j*est_mean_j)));
  if((x_j==0.0)&&(x_i!=0.0))  dens = dbnorm(x_i,x_j,0,0,1,corr)*x_i/(est_mean_i*(1 + LambertW(tail*est_mean_i*est_mean_i)));
  if((x_j==0.0)&&(x_i==0.0))  dbnorm(x_i,x_j,0,0,1,corr);
  return(dens);
}








/***** bivariate half tukey h ****/     
double biv_half_Tukeyh(double rho,double ti,double tj,double tail)
{
  double dens = 0.0;
  dens = biv_tukey_h(rho,ti,tj,0,0,tail,1) + biv_tukey_h(rho,-ti,-tj,0,0,tail,1) +
    biv_tukey_h(rho,-ti,tj,0,0,tail,1) +
    biv_tukey_h(rho,ti,-tj,0,0,tail,1);
  return(dens);
}
 

/***** bivariate two piece tukey h ****/ 
double biv_two_pieceTukeyh(double rho,double zi,double zj,double sill,double eta,double tail,
             double p11,double mui,double muj)
{
double res=0.0;  
double etamas=1+eta;
double etamos=1-eta;
double zistd=(zi-mui)/sqrt(sill);
double zjstd=(zj-muj)/sqrt(sill);
/*if(rho)   {*/
//if(zi>=mui&&zj>=muj)
    if(zistd>=0&&zjstd>=0)
{res=          (p11/R_pow(etamos,2))*biv_half_Tukeyh(rho,zistd/etamos,zjstd/etamos,tail);}
//if(zi>=mui&&zj<muj)
    if(zistd>=0&&zjstd<0)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*biv_half_Tukeyh(rho,zistd/etamos,-zjstd/etamas,tail);}
//if(zi<mui&&zj>=muj)
      if(zistd<0&&zjstd>=0)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*biv_half_Tukeyh(rho,-zistd/etamas,zjstd/etamos,tail);}
//if(zi<mui&&zj<muj)
    if(zistd<0&&zjstd<0)
{res=    ((p11+eta)/R_pow(etamas,2))*biv_half_Tukeyh(rho,-zistd/etamas,-zjstd/etamas,tail);}
/*}else{   if(zi>=mui)
         {res=dnorm(x_i,0,1,0)*x_i/((zistd/etamos)*(1+LambertW(tail*R_pow(zistd/etamos,2))));}
         if(zj<muj)
         {res=dnorm(x_j,0,1,0)*x_j/((zjstd/etamas)*(1+LambertW(tail*R_pow(zjstd/etamas,2))));}
      }*/
return(res/sill);
}

/*******************************************************************************/

double  binomialCoeff(int n, int k) 
{ 
    double res=0.0;
    res=lgammafn(n+1)-(lgammafn(k+1)+lgammafn(n-k+1));
    return(exp(res)); 
} 

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

double Prt(double corr, int r, int t, double mean_i, double mean_j) {
    const double rho2 = R_pow(corr, 2);
    const double one_minus_rho2 = 1.0 - rho2;
    const double log_rho2 = log(rho2);
    const double log_rho2_ratio = log_rho2 - log1p(-rho2);
    const double log_mean_i = log(mean_i);
    const double auxi = mean_i / one_minus_rho2;
    const double auxj = mean_j / one_minus_rho2;
    const double rho2auxi = rho2 * auxi;

    const int n = r - t;
    const int iter1 = 4000;
    const int iter2 = 3000;

    double sum = 0.0, sum1 = 0.0;
    double  prev_sum1 = 0.0;

    for (int m = 0; m <= iter1; m++) {
        int tm = t + m;
        double gamma_tm = lgammafn(tm);
        double gamma_m1 = lgammafn(m + 1);
        double gamma_t = lgammafn(t);

        // Prima somma (con k)
        double inner_sum = 0.0;
        double prev_inner_sum = 0.0;

        for (int k = 0; k <= iter2; k++) {
            int cc = t + m + n + k + 1;
            double gamma_cc = lgammafn(cc);
            double log_igam = log(igam(1 + k + tm, auxj));
            double hyper = hyperg(n, cc, rho2auxi);

            double log_q1 = R_finite(hyper) ? log(hyper) - gamma_cc
                                            : log(aprox_reg_1F1(n, cc, rho2auxi));

            double term = (k + m) * log_rho2_ratio
                          + gamma_tm
                          - gamma_m1 - gamma_t
                          + (cc - 1) * log_mean_i
                          + log_igam
                          + log_q1;

            double val = exp(term);
            if (!R_finite(val) || val < 1e-12) break;

            inner_sum += val;
            if (fabs(inner_sum - prev_inner_sum) < 1e-10) break;
            prev_inner_sum = inner_sum;
        }
        sum += inner_sum;
        // Seconda somma
        double q2_val;
        int cc2 = tm + n + 1;
        double hyper2 = hyperg(n + 1, cc2, rho2auxi);
        if (R_finite(hyper2))
            q2_val = exp(log(hyper2) - lgammafn(cc2));
        else
            q2_val = aprox_reg_1F1(n + 1, cc2, rho2auxi);
        double term1 = m * log_rho2_ratio
                       + gamma_tm
                       + (tm + n) * log_mean_i
                       - gamma_m1 - gamma_t
                       + log(q2_val)
                       + log(igam(tm, auxj));
        double val1 = exp(term1);
        if (!R_finite(val1) || val1 < 1e-12) break;

        sum1 += val1;
        if (fabs(sum1 - prev_sum1) < 1e-10) break;
        prev_sum1 = sum1;
    }

    double prt = exp(-auxi + log(sum1)) - exp(-auxi + log(sum));
    return prt;
}


//*****************************************************************************/
double Prr(double corr, int r, int t, double mean_i, double mean_j) {
    double rho2 = R_pow(corr, 2);    
    double prr = 0.0, term = 0.0, term1 = 0.0, term2 = 0.0, term3 = 0.0;
    double sum = 0.0, res00 = 0.0, res11 = 0.0, sum1 = 0.0, sum2 = 0.0, rho2k = 0.0;
    int k = 0, m = 0;

    // Precompute constants that do not change during loops
    double auxi = mean_i / (1 - rho2);
    double auxj = mean_j / (1 - rho2);
    double ei = exp(-mean_i);
    double ej = exp(-mean_j);
    int iter1 = 1500, iter2 = 1000;

    // Loop over k
    while (k < iter1) {
        rho2k = R_pow(rho2, k);  // Precompute powers of rho2
        // Inner loop over m
        m = 0;
        while (m < iter2) {
            term = (1 - rho2) * rho2k * R_pow(rho2, m) * 
                   exp(lgammafn(r + m) - lgammafn(r) - lgammafn(m + 1) + 
                        log(igam(r + k + m + 1, auxi)) + log(igam(r + k + m + 1, auxj)));

            // Early exit if term is small or non-finite
            if (fabs(term) < 1e-10 || !R_finite(term)) {
                break;
            }
            sum += term;
            m++;
        }
        // Recompute terms based on k
        double aa = lgammafn(k + 1) + lgammafn(r);
        double bb = lgammafn(r + k);
        double cc = igam(r + k, auxi);
        double dd = igam(r + k, auxj);
        term1 = rho2k * exp(bb + log(cc) + log(dd) - aa);
        term2 = ei * R_pow(1 / rho2, r) * exp(bb + log(igam(r + k, rho2 * auxi)) + log(dd) - aa);
        term3 = ej * R_pow(1 / rho2, r) * exp(bb + log(cc) + log(igam(r + k, rho2 * auxj)) - aa);
        // Early exit if any of the terms are non-finite
        if (!R_finite(term1) || !R_finite(term2) || !R_finite(term3)) {
            break;
        }
        sum1 += term1;
        sum2 += term2 + term3;
        // Check for convergence
        if (fabs(sum1 - res00) < 1e-10 && fabs(sum2 - res11) < 1e-10) {
            break;
        } else {
            res00 = sum1;
            res11 = sum2;
        }
        k++;
    }
    // Final computation for prr
    prr = pow1p(-rho2, r) * (-sum1 + sum2 + sum);
    return prr;
}



/***************************************************************/
double Pr0(double corr, int r, int t, double mean_i, double mean_j) {
    double rho2 = R_pow(corr, 2);
    double one_minus_rho2 = 1.0 - rho2;
    double log_rho2 = log(rho2);
    double log_one_minus_rho2 = log1p(-rho2);
    double log_rho_ratio = log_rho2 - log_one_minus_rho2;

    double auxi = mean_i / one_minus_rho2;
    double auxj = mean_j / one_minus_rho2;
    double log_mean_i = log(mean_i);
    
    double term = 0.0, q2 = 0.0, sum = 0.0, res0 = 0.0;
    double aux = 0.0, aux1 = 0.0;
    
    int n = r - t;
    int m = 0;
    const int iter = 4000;
    while (m <= iter) {
        aux = m * log_rho_ratio;
        aux1 = (m + n) * log_mean_i;
        double hyper = hyperg(n, m + n + 1, rho2 * auxi);
        if (!R_finite(hyper)) break;
        q2 = exp(log(hyper) - lgammafn(m + n + 1));
        double ig = igam(m + 1, auxj);
        if (!R_finite(q2) || !R_finite(ig)) break;
        term = exp(aux + aux1 + log(q2) + log(ig));
        if (term < 1e-12) break;
        sum += term;
        if (fabs(sum - res0) < 1e-10) break;
        res0 = sum;
        m++;
    }
    double part1 = exp(-mean_i + n * log_mean_i - lgammafn(n + 1));
    double part2 = exp(-auxi + log(sum));
    double pr0 = part1 - part2;

    return pr0;
}


/*******************************************************************************/
double P00(double corr, int r, int t, double mean_i, double mean_j) {
    const double rho2 = R_pow(corr, 2);
    const double one_minus_rho2 = 1.0 - rho2;
    const double log_rho2 = log(rho2);
    const double auxi = mean_i / one_minus_rho2;
    const double auxj = mean_j / one_minus_rho2;

    double sum = 0.0, prev_sum = 0.0;

    for (int k = 0; k < 4000; ++k) {
        int kp1 = k + 1;

        double igam_i = igam(kp1, auxi);
        double igam_j = igam(kp1, auxj);

        if (igam_i == 0 || igam_j == 0) break;

        double log_term = k * log_rho2 + log(igam_i) + log(igam_j);
        double term = exp(log_term);

        if (!R_finite(term) || term < 1e-12) break;

        sum += term;

        if (fabs(sum - prev_sum) < 1e-10) break;
        prev_sum = sum;
    }

    double p00 = -1 + exp(-mean_i) + exp(-mean_j) + one_minus_rho2 * sum;
    return p00;
}



void biv_pois_call(double *corr,int *r, int *t, double *mean_i, double *mean_j,double *res)
{
    *res = biv_Poisson(*corr,*r,*t,*mean_i,*mean_j);
}



double biv_Poisson(double corr,int r, int t, double mean_i, double mean_j)
{
double dens=0.0;
if(fabs(corr)>1e-6){
  if(r==t)
  {    if(r==0) dens=P00(corr,r,r,mean_i,mean_j);
       if(r>0)  dens=Prr(corr,r,r,mean_i,mean_j);
  }

else{
  if(r==0&&t>0) dens=Pr0(corr,t,r,mean_j,mean_i);
  if(r>0&&t==0) dens=Pr0(corr,r,t,mean_i,mean_j);

  if(r>0&&t>0)
   {  
   if(r>t) dens=Prt(corr,r,t,mean_i,mean_j);
   if(t>r) dens=Prt(corr,t,r,mean_j,mean_i);
   }
  }
}
else{
    //double lambda_i=exp(mean_i); double lambda_j=exp(mean_j);  
    double dens1= -mean_i +r*log(mean_i)-lgammafn(r+1);
    double dens2= -mean_j +t*log(mean_j)-lgammafn(t+1);
    dens=exp(dens1+dens2);
} 

return(dens);

}


double biv_PoissonZIP(double corr,int r, int t, double mean_i, double mean_j,double mup,double nugget1,double nugget2)
{
double dens=0.0,p,p00,p10,p01,p11;
p=pnorm(mup,0,1,1,0);
p00=pbnorm22(mup,mup,(1-nugget2)*corr);
p01=p-p00;
p10=p01;
p11=1-2*p+p00;

//Rprintf("%f %f \n",nugget1,nugget2);

if(r==0&&t==0)
     dens=p00  + p01*exp(-mean_i) + p10*exp(-mean_j)+p11*biv_Poisson((1-nugget1)*corr,0, 0, mean_i, mean_j);
if(r==0&&t>0)
      //dens=      p01*  exp(-mean_i+t*log(mean_i)-lgammafn(t+1))  + p11*biv_Poisson((1-nugget1)*corr,0, t, mean_i, mean_j);
      dens=      p01*  exp(-mean_j+t*log(mean_j)-lgammafn(t+1))  + p11*biv_Poisson((1-nugget1)*corr,0, t, mean_i, mean_j);
if(r>0&&t==0)
      //dens=      p10*  exp(-mean_j+r*log(mean_j)-lgammafn(r+1))  + p11*biv_Poisson((1-nugget1)*corr,r, 0, mean_i, mean_j);
      dens=      p10*  exp(-mean_i+r*log(mean_i)-lgammafn(r+1))  + p11*biv_Poisson((1-nugget1)*corr,r, 0, mean_i, mean_j);
if(r>0&&t>0)
      dens=      p11*biv_Poisson((1-nugget1)*corr,r, t, mean_i, mean_j);
return(dens);

}





double biv_Mis_PoissonZIP(double corr, double data_i, double data_j,
                          double mean_i, double mean_j, double mup,
                          double nugget1, double nugget2) {
    const int N = 2;
    double M[2][2];
    double dat[2];
    double dens = 0.0;

    double p = pnorm(mup, 0, 1, 1, 0);
    double one_minus_nug2_corr = (1.0 - nugget2) * corr;
    double one_minus_nug1_corr = (1.0 - nugget1) * corr;

    double p00 = pbnorm22(mup, mup, one_minus_nug2_corr);
    double p01 = p - p00;
    double p10 = p01;
    double p11 = 1.0 - 2.0 * p + p00;

    double corr1 = corr_pois(one_minus_nug1_corr, mean_i, mean_j);

    double sqrt_mean_i = sqrt(mean_i);
    double sqrt_mean_j = sqrt(mean_j);

    M[0][0] = mean_i;
    M[1][1] = mean_j;
    M[0][1] = M[1][0] = sqrt(mean_i * mean_j) * corr1;

    dat[0] = data_i - mean_i;
    dat[1] = data_j - mean_j;

    double pdf2 = dNnorm(N, (double **)M, dat); // Bivariate normal
    double pdf1i = dnorm(data_i, mean_i, sqrt_mean_i, 0); // Univariate marginals
    double pdf1j = dnorm(data_j, mean_j, sqrt_mean_j, 0);

    if (data_i > 0.0 && data_j > 0.0)
        dens = p11 * pdf2;
    else if (data_i > 0.0 && data_j == 0.0)
        dens = p11 * pdf2 + p10 * pdf1i;
    else if (data_i == 0.0 && data_j > 0.0)
        dens = p11 * pdf2 + p01 * pdf1j;
    else // data_i == 0.0 && data_j == 0.0
        dens = p11 * pdf2 + p01 * pdf1i + p10 * pdf1j + p00;

    return dens;
}


/*****/
double biv_binomnegZINB(int N,double corr,int r, int t, double mean_i, double mean_j,
    double nugget1,double nugget2,double mup)
{
double dens=0.0,ap,ap00,ap10,ap01,ap11;
ap=pnorm(mup,0,1,1,0);

ap00=pbnorm22(mup,mup,(1-nugget2)*corr);
ap01=ap-ap00;
ap10=ap01;
ap11=1-2*ap+ap00;

double p11=pbnorm22(mean_i,mean_j,(1-nugget1)*corr);
double p1=pnorm(mean_i,0,1,1,0);
double p2=pnorm(mean_j,0,1,1,0);


if(r==0&&t==0)
     dens=ap00  + ap01*R_pow(p1,N) + ap10*R_pow(p2,N)+ap11*biv_binomneg(N,0, 0 ,p1, p2, p11);
if(r==0&&t>0)
      dens=      ap01*  exp(lgammafn(N+t)-lgammafn(t+1)-lgammafn(N) +N*log(p2)+t*log1p(-p2) ) +
                               ap11*biv_binomneg(N,0, t, p1, p2, p11);
if(r>0&&t==0)
      dens=      ap10* exp(lgammafn(N+r)-lgammafn(r+1)-lgammafn(N) +N*log(p1)+r*log1p(-p1) ) +
                               ap11*biv_binomneg(N,r, 0, p1, p2, p11);
if(r>0&&t>0)
      dens=      ap11*biv_binomneg(N,r, t,p1, p2, p11);
return(dens);

}
/* binary misspecified negative binomial*/
double biv_binegbinary (int NN, int u, int v,double pu,double pv, double p11)
{
double res=0.0;double u_0,v_0,p00;

u_0=R_pow(pu,NN);
v_0=R_pow(pv,NN);
p00=R_pow(p11,NN);

if(u==0&&v==0)  res=p00;
if(u==0&&v>0)   res=u_0-p00;
if(v==0&&u>0)   res=v_0-p00;
if(u>0&&v>0)    res=1-(u_0-p00+v_0);
return(res);
}








/************* POIsson gammma ************/
double PG00(double corr, int r, int t, double mean_i, double mean_j, double a) {
    double rho2 = corr * corr;
    double beta_i = a / mean_i;
    double beta_j = a / mean_j;
    double beta_ij = beta_i * beta_j;

    double auxi = 1.0 / (1.0 + beta_i);
    double auxj = 1.0 / (1.0 + beta_j);
    double auxij = auxi * auxj;

    double log_beta_ij = log(beta_ij);
    double log_rho2 = log(rho2);
    double log_auxij = log(auxij);
    double log_p1mrho2 = log1p(-rho2);  // log(1 - rho2)

    double sum = 0.0, res0 = 0.0;
    const int iter1 = 600, iter2 = 600;

    for (int k = 0; k < iter1; ++k) {
        for (int l = 0; l < iter2; ++l) {
            double log_aa = (l + a - 1) * log_beta_ij + (k + l) * log_rho2 +
                            (k + l + a) * log_auxij + (a + 1) * log_p1mrho2;

            double bb = 2 * lgammafn(k + l + a + 1) - 2 * lgammafn(k + 2)
                        - lgammafn(l + 1) - lgammafn(a) - lgammafn(l + a);

            double hg1 = hypergeo(1, 1 - l - a, k + 2, -1.0 / beta_i);
            double hg2 = hypergeo(1, 1 - l - a, k + 2, -1.0 / beta_j);

            if (!R_finite(hg1) || !R_finite(hg2)) continue;

            double term = exp(log_aa + bb) * hg1 * hg2;

            if (!R_finite(term) || fabs(term) < 1e-30) break;

            sum += term;
        }

        if (fabs(sum - res0) < 1e-30) break;
        res0 = sum;
    }

    double p00 = -1.0 + R_pow(auxi * beta_i, a) + R_pow(auxj * beta_j, a) + sum;
    if (p00 < 1e-320) p00 = 1e-320;

    return p00;
}


double PGrr(double corr, int r, int t, double mean_i, double mean_j, double a) {
    const double rho2 = corr * corr;
    const double beta_i = a / mean_i;
    const double beta_j = a / mean_j;
    const double beta_ij = beta_i * beta_j;
    const double bmi = -1.0 / beta_i;
    const double bmj = -1.0 / beta_j;

    const double auxi = 1.0 / (1.0 + beta_i);
    const double auxj = 1.0 / (1.0 + beta_j);
    const double auxij = auxi * auxj;

    const double bri = 1.0 + beta_i - rho2;
    const double brj = 1.0 + beta_j - rho2;
    const double fji = brj * beta_i * auxij;
    const double fij = bri * beta_j * auxij;

    const double ra = r + a;
    const double log_rho2 = log(rho2);
    const double log_beta_ij = log(beta_ij);
    const double log_auxij = log(auxij);
    const double log1mrho2 = log1p(-rho2);
    const double r2br = -rho2 / bri;

    const int iter1 = 600, iter2 = 500, iter3 = 400;

    double sum = 0.0, sum1 = 0.0, sum2 = 0.0;
    double res00 = 0.0, res11 = 0.0;

    for (int k = 0; k < iter1; ++k) {
        int mm = r + k + 1;

        for (int l1 = 0; l1 < iter2; ++l1) {
            double ff = 1.0 - l1 - a;

            for (int l = 0; l < iter3; ++l) {
                int mml1 = mm + l + 1;
                int kll1 = k + l + l1;
                double log_aa = (l1 + a - 1) * log_beta_ij + kll1 * log_rho2 + (ra + kll1) * log_auxij + (ra + 1) * log1mrho2;
                double bb = lgammafn(r + l) + 2 * lgammafn(mm + l + l1 + a) - 2 * lgammafn(mml1) -
                            lgammafn(l + 1) - lgammafn(l1 + 1) - lgammafn(r) - lgammafn(a) - lgammafn(l1 + a);

                double hg1 = hypergeo(1, 1 - l1 - a, mml1, bmi);
                double hg2 = hypergeo(1, ff, mml1, bmj);
                if (!R_finite(hg1) || !R_finite(hg2)) break;

                double term = exp(log_aa + bb) * hg1 * hg2;
                if (fabs(term) < 1e-30 || !R_finite(term)) break;

                sum += term;
            }

            // Term1, Term2, Term3
            int kl1 = k + l1;
            double log_aa1 = (l1 + a) * log_beta_ij + kl1 * log_rho2 + (ra + kl1) * log_auxij + ra * log1mrho2;
            double bb1 = exp(lgammafn(r + k) + 2 * lgammafn(ra + kl1) - 2 * lgammafn(mm) -
                             lgammafn(k + 1) - lgammafn(l1 + 1) - lgammafn(r) - lgammafn(a) - lgammafn(l1 + a));

            double h1 = hypergeo(1, ff, mm, bmi);
            double h2 = hypergeo(1, ff, mm, bmj);
            double h3 = hypergeo(1, ff, mm, r2br);
            double h4 = hypergeo(1, ff, mm, -rho2 / brj);

            if (!R_finite(h1) || !R_finite(h2) || !R_finite(h3) || !R_finite(h4)) break;

            double aa1 = exp(log_aa1);
            double term1 = aa1 * h1 * h2 * bb1 / (beta_ij * auxij);
            double term2 = aa1 * h3 * h2 * bb1 / fij;
            double term3 = aa1 * h1 * h4 * bb1 / fji;

            if (fabs(term1) < 1e-30 || fabs(term2) < 1e-30 || fabs(term3) < 1e-30) break;

            sum1 += term1;
            sum2 += term2 + term3;
        }

        if (fabs(sum1 - res00) < 1e-30 && fabs(sum2 - res11) < 1e-30) break;

        res00 = sum1;
        res11 = sum2;
    }

    double prr = -sum1 + sum2 + sum;
    if (prr < 1e-320) prr = 1e-320;
    return prr;
}

double PGr0(double corr, int r, int t, double mean_i, double mean_j, double a) {
    const double rho2 = corr * corr;
    const double beta_i = a / mean_i;
    const double beta_j = a / mean_j;
    const double auxi = 1.0 / (1.0 + beta_i);
    const double auxj = 1.0 / (1.0 + beta_j);
    const double auxij = auxi * auxj;
    const double brho = beta_i - rho2;


    const int n = r - t;
    const double an = a + n;
    const double rb = -rho2 / (1.0 + brho);
    const double ib = -1.0 / beta_j;
    const double beta_i_pow_a = R_pow(beta_i, a);
    const double auxi_pow_na = R_pow(auxi, n + a);
    const double aux1 = beta_i_pow_a * auxi_pow_na *
                        exp(lgammafn(n + a) - lgammafn(n + 1) - lgammafn(a));

    double sum1 = 0.0;
    const int iter1 = 700;
    const int iter2 = 500;

    for (int l = 0; l < iter1; l++) {
        const int q = l + 1;
        const int nq = n + q;

        for (int l1 = 0; l1 < iter2; l1++) {
            const int ii = l + l1;
            const double l1a = l1 + a;
            const int cc = ii + (int)a;

            const double beta_i_l1a = R_pow(beta_i, l1a);
            const double beta_j_l1a_1 = R_pow(beta_j, l1a - 1.0);
            const double rho2_pow_ii = R_pow(rho2, ii);
            const double auxij_pow_cc = R_pow(auxij, cc);

            const double aa = beta_i_l1a * beta_j_l1a_1 * rho2_pow_ii *
                              pow1p(-rho2, an) * auxij_pow_cc * pow1p(brho, -n);

            const double bb = lgammafn(n + cc) + lgammafn(cc + 1) -
                              lgammafn(nq) - lgammafn(q + 1) -
                              lgammafn(l1 + 1) - lgammafn(a) - lgammafn(l1a);

            double hg1 = hypergeo(n, 1.0 - l1a, nq, rb);
            double hg2 = hypergeo(1.0, 1.0 - l1a, q + 1, ib);

            double term1 = aa * hg1 * hg2 * exp(bb);
            if (fabs(term1) < 1e-30 || !R_finite(term1))
                break;

            sum1 += term1;
        }
    }

    double pr0 = aux1 - sum1;
    if (pr0 < 1e-320) pr0 = 1e-320;
    return pr0;
}

double PGrt(double corr, int r, int t, double mean_i, double mean_j, double a) {
    double rho2 = R_pow(corr, 2);
    double beta_i = a / mean_i;
    double beta_j = a / mean_j;
    double brho = beta_i - rho2;

    double auxi = 1 / (beta_i + 1);
    double auxj = 1 / (beta_j + 1);
    double auxij = auxi * auxj;
    double rtrho = -rho2 / (1 + brho);
    double mbj = -1 / beta_j;

    double sum1 = 0.0, sum2 = 0.0;
    double res11 = 0.0;
    double prt = 0.0;

    int n = r - t;
    int s = n + 1;
    int tna = t + n + a;

    int iter1 = 700, iter2 = 500, iter3 = 400;

    for (int l = 0; l < iter1; l++) {
        int tl = t + l;
        int u = t + l;
        int u1 = u + 1;
        int us = u + s;
        int f = l + 1;

        for (int l1 = 0; l1 < iter2; l1++) {
            double jj = l1 + a;
            double jj1 = 1 - jj;
            int mm = u + jj;
            int ll1 = l + l1;
            int d = l1 + 1;

            for (int k = 0; k < iter3; k++) {
                int cc = u + k;
                int ii = cc + n;

                double aa = R_pow(beta_i, jj) * R_pow(beta_j, jj - 1) *
                            R_pow(rho2, k + ll1) * pow1p(-rho2, tna) *
                            R_pow(auxij, cc + jj) * pow1p(brho, -n);

                double bb = lgammafn(tl) + lgammafn(ii + jj) + lgammafn(cc + jj + 1) -
                            lgammafn(ii + 1) - lgammafn(cc + 2) - lgammafn(d) -
                            lgammafn(f) - lgammafn(t) - lgammafn(a) - lgammafn(jj);

                double term1 = aa * hypergeo(n, jj1, ii + 1, rtrho) *
                               hypergeo(1, jj1, cc + 2, mbj) * exp(bb);

                if (fabs(term1) < 1e-30 || !R_finite(term1)) break;
                sum1 += term1;
            }

            double aa1 = R_pow(beta_i, jj) * R_pow(beta_j, -jj1) *
                         R_pow(rho2, ll1) * pow1p(-rho2, tna) *
                         R_pow(auxij, mm - 1) * pow1p(brho, -s);

            double bb1 = lgammafn(tl) + lgammafn(n + mm) + lgammafn(mm) -
                         lgammafn(us) - lgammafn(u1) - lgammafn(d) -
                         lgammafn(f) - lgammafn(t) - lgammafn(a) - lgammafn(jj);

            double term2 = aa1 * hypergeo(s, jj1, us, rtrho) *
                           hypergeo(1, jj1, u1, mbj) * exp(bb1);

            if (fabs(term2) < 1e-30 || !R_finite(term2)) break;
            sum2 += term2;
        }

        if (fabs(sum2 - res11) < 1e-30) break;
        res11 = sum2;
    }

    prt = sum2 - sum1;
    if (prt < 1e-320) prt = 1e-320;
    return prt;
}


double biv_PoissonGamma(double corr,int r, int t, double mean_i, double mean_j, double a)
{
double dens=0.0;
if(fabs(corr)>1e-120){
  if(r==t)
  {    if(r==0) dens=PG00(corr,r,r,mean_i,mean_j,a);
       if(r>0)  dens=PGrr(corr,r,r,mean_i,mean_j,a);
  }

else{
  if(r==0&&t>0) dens=PGr0(corr,t,r,mean_j,mean_i,a);
  if(r>0&&t==0) dens=PGr0(corr,r,t,mean_i,mean_j,a);

  if(r>0&&t>0)
   {
   if(r>t) dens=PGrt(corr,r,t,mean_i,mean_j,a);
   if(t>r) dens=PGrt(corr,t,r,mean_j,mean_i,a);
   }
  }
}
else{
    double beta_i= a/mean_i;
    double beta_j= a/mean_j;
    double dens1= r*log(1/(1+beta_i))+a*log(beta_i/(1+beta_i))+lgammafn(a+r)-lgammafn(r+1)-lgammafn(a);
    double dens2= t*log(1/(1+beta_j))+a*log(beta_j/(1+beta_j))+lgammafn(a+t)-lgammafn(t+1)-lgammafn(a);
    dens=exp(dens1+dens2);
}
//Rprintf("%f \n",dens);
return(dens);

}



double biv_PoissonGammaZIP(double corr,int r, int t, double mean_i, double mean_j,double mup,double nugget1,double nugget2,double shape)
{
double dens=0.0,p,p00,p10,p01,p11;
p=pnorm(mup,0,1,1,0);
p00=pbnorm22(mup,mup,(1-nugget2)*corr);
p01=p-p00;
p10=p01;
p11=1-2*p+p00;


if(r==0&&t==0)
dens=p00  + p01* R_pow(shape/(mean_i+ shape ),shape) + p10*R_pow(shape/(mean_j+ shape ),shape)+p11*biv_PoissonGamma((1-nugget1)*corr,0, 0, mean_i, mean_j,shape);
if(r==0&&t>0)
dens=      p01*   dnbinom(0, shape,   mean_j /(mean_j+ shape ),0)  + p11*biv_PoissonGamma((1-nugget1)*corr,0, t, mean_i, mean_j,shape);
if(r>0&&t==0)
dens=      p10*   dnbinom(0, shape,   mean_i /(mean_i+ shape ),0)  + p11*biv_PoissonGamma((1-nugget1)*corr,r, 0, mean_i, mean_j,shape);
if(r>0&&t>0)
dens=      p11*biv_PoissonGamma((1-nugget1)*corr,r, t, mean_i, mean_j,shape);
return(dens);

}




/*####*/
double int_kuma(double x,double eta, double gam,double k,double m)
{
    double res=0.0;

    res=R_pow( 1-R_pow(x,1/eta),1/gam)*R_pow(x,k-m)*R_pow(1-x,m);
    return (res);
}

/*####*/
void integr_kuma(double *x, int n, void *ex){
    int i;double eta,gam,k,m;
    eta =    ((double*)ex)[0];  
    gam =    ((double*)ex)[1];  
    k=       ((double*)ex)[2];  
    m=       ((double*)ex)[3];  

    for (i=0;i<n;i++) {
        x[i]=int_kuma(x[i],eta,gam,k,m);}
    return;
}


double kumaintegral(double *param) {
    double ex[4], lower, upper, epsabs, epsrel, result, abserr, *work;
    int neval, ier, subdiv, lenw, last, *iwork;
    subdiv = 100;epsabs = R_pow(DBL_EPSILON, 0.25);epsrel = epsabs;    lenw = 4 * subdiv;           /* as instructed in WRE */
    iwork =   (int *) R_Calloc(subdiv, int);  /* idem */
    work = (double *) R_Calloc(lenw, double); /* idem */
    ex[0] = param[0]; //eta
    ex[1] = param[1]; //gam 
    ex[2] = param[2];  //k
    ex[3] = param[3];  //m
    lower=0;
    upper=1;
    Rdqags(integr_kuma, (void *) &ex,&lower, &upper, 
               &epsabs, &epsrel, &result,
               &abserr, &neval, &ier, &subdiv, &lenw, &last, iwork, work);
    R_Free(iwork);R_Free(work);
    return(result);
}







/*######*/
  double mean_kuma(double eta,double gam){
    double out=eta*beta(1+(1/gam),eta);
    return(out);
  }
/*******/
   double var_kuma(double eta, double gam){
    double mm=mean_kuma(eta,gam);
    double out=eta*beta(1+2*(1/gam),eta)-mm*mm;
    return(out);
  }

/******/
/******/
double corr_kuma(double rho,double eta,double gam){

  double corr=0.0, tol=1e-6; int iter=0;
  double rho2=R_pow(rho,2);
  int k=0;int m=0;

  if(fabs(rho)< tol) return(0.0);

  /*###############*/
if (eta==1.0&&gam==1.0){
    corr=((2*(rho2*(3*rho2-1)-R_pow(rho2-1,2)*log1p(-rho2)))/R_pow(rho2,2))-3;
    return(corr);
  }
else
{
  if (eta==1.0&&gam!=1.0){
    double res_K=0.0,res_M=0.0,sum_M=0.0, bb,aa;
    iter=10000;
    tol=1e-14;
    while(k<=iter){
      res_M=0;m=0;
      bb=     2*(log1p(-rho2) + k*log(rho));
      while(m<=k){

        //A=exp(lbeta(1+k-m,1+(1/gam)+m));
        aa=    -2*lbeta(k-m+1,m+1);
        sum_M= exp(aa + bb + 2*lbeta(1+k-m,1+(1/gam)+m));
        res_M=res_M+  sum_M;
        if ((sum_M<tol)|(sum_M>1e300)){break;}
       m=m+1;
      }
      res_K=res_K+res_M;
      if (res_M<tol){ break;}
      k=k+1;
    }
    double mm=mean_kuma(eta,gam);
    double vv=var_kuma(eta,gam);
    corr=(res_K-R_pow(mm,2))/vv;
  }
/******/
if (eta!=1.0&&gam==1.0){
   double res_K=0.0,res_M=0.0,sum_M=0.0,bb,aa,c1,c2,c3;
    iter=10000;
    tol=1e-14;
    while(k<=iter){
      res_M=0;m=0;
        bb= 2*(log1p(-rho2) + k*log(rho));
      while(m<=k){
        aa=-2*lbeta(k-m+1,m+1);
        c1=lgammafn(1+m)+lgammafn(1+k-m)-lgammafn(2+k);
        c2=lgammafn(1+m)+lgammafn(1+(1/eta)+k-m)-lgammafn(2+(1/eta)+k);
        c3=aa + bb;
        sum_M= exp(2*c1+c3)+exp(2*c2+c3)-2*exp(c1+c2+c3);
        res_M=res_M+  sum_M;
        if ((sum_M<tol)|(sum_M>1e300)){ break;}
        m=m+1;
      }
      res_K=res_K+res_M;
      if (res_M<tol){break;}
      k=k+1;
    }
    double mm=mean_kuma(eta,gam);
    double vv=var_kuma(eta,gam);
    corr=(res_K-R_pow(mm,2))/vv;
  }
  
 


  /**+*********************/
  if ((eta!=1.0)&&(gam!=1.0)){
   
       double *param;;
       param=(double *) R_Calloc(4,double);
       param[0]=eta;param[1]=gam;  //mu,alpha //beta
       
 
    double res_K=0.0,res_M=0.0,sum_M=0.0,aa,bb,p2,p1;
    k=0;res_K=0;iter=10000;tol=1e-14;
     while (k<=iter){
      res_M=0.0;m=0;
      bb= 2*(log1p(-rho2) + k*log(rho));
       param[2]=k;
      while(m<=k){
        param[3]=m;
        p1=kumaintegral(param);
        aa=-2*lbeta(k-m+1,m+1);
        p2=exp(aa+bb);
        sum_M=p2*p1*p1;
        res_M=res_M+  sum_M;
        if ((sum_M<tol) | (sum_M>1e300)){break;}
        m=m+1;
      }
      res_K=res_K+res_M;
      if (res_M<tol){ break;}
      k=k+1;
    }
    double mm=mean_kuma(eta,gam);
    double vv=var_kuma(eta,gam);
    corr=(res_K-R_pow(mm,2))/vv;
  
  }
    return(corr);
 }
 
}

  /*********************/
void corr_kuma_vec(double *rho,double *eta,double *gam,double *res, int *n)
{
int i=0;
for(i=0;i<=*n;i++)  res[i]=corr_kuma(rho[i],eta[0],gam[0]);
}

 /*********************/

void biv_unif_CopulaGauss_call(double *x,double *y,double *rho, double *res)
{
    *res = biv_unif_CopulaGauss(*x,*y,*rho);
}

// bivariatte density  gaussian copula
double biv_unif_CopulaGauss(double dat1, double dat2, double rho)
{
    if (fabs(dat1 - 1) < 0.0001) dat1 = 0.999;
    if (fabs(dat2 - 1) < 0.0001) dat2 = 0.998;

    // Calcola i quantili
    double a1 = qnorm(dat1, 0, 1, 1, 0);
    double a2 = qnorm(dat2, 0, 1, 1, 0);

    // Calcola le densità normali
    double g1 = dnorm(a1, 0, 1, 0);
    double g2 = dnorm(a2, 0, 1, 0);

    // Calcola e restituisci il risultato
    return biv_Norm(rho, a1, a2, 0, 0, 1, 1, 0) / (g1 * g2);
}


// bivariatte density  skewgaussian copula
double biv_unif_CopulaSkewGauss(double dat1, double dat2, double rho, double alpha)
{
    if (fabs(dat1 - 1) < 0.0001) dat1 = 0.999;
    if (fabs(dat2 - 1) < 0.0001) dat2 = 0.998;
    
    const double small = 1e-8;
    const double omega = sqrt(alpha * alpha + 1);
    
    // Calcola i quantili in una sola riga
    double a1 = qsn(dat1, omega, alpha, 0, small);
    double a2 = qsn(dat2, omega, alpha, 0, small);

    // Calcola le densità dei quantili
    double g1 = dsn(a1, omega, alpha, 0);
    double g2 = dsn(a2, omega, alpha, 0);

    // Calcola e restituisci il risultato
    return biv_skew(rho, a1, a2, 0, 0, 1, alpha, 0) / (g1 * g2);
}



 /*********************/
void biv_unif_CopulaClayton_call(double *x,double *y,double *rho, double *nu, double *res)
{
    *res = biv_unif_CopulaClayton(*x,*y,*rho,*nu);
}

// bivariate density clayton copula
double biv_unif_CopulaClayton(double dat1, double dat2, double rho, double nu)
{
    // Precalcola variabili
    double nu2 = nu / 2;
    double rho2 = rho * rho;
    double a = nu2 + 1;

    // Calcola i termini a1 e a2
    double a1 = R_pow(dat1, 1 / nu2);
    double a2 = R_pow(dat2, 1 / nu2);

    // Calcola e restituisce il risultato
    return a * log1p(-rho2) + log(appellF4(a, a, nu2, 1, rho2 * a1 * a2, rho2 * (1 - a1) * (1 - a2)));
}


double cdf_kuma(double y,double a, double b){
double res=1-R_pow(1-R_pow(y,a),b);
return(res);
}
double pdf_kuma(double y,double a, double b){
double res=(a*b)* R_pow(y,a-1)*R_pow(1-R_pow(y,a),b-1);
return(res);
}

/*********************************************************************/
double biv_cop(double rho,int type_cop,int cond,
             double z1,double z2,double mu1,double mu2,double *nuis,int model, int NN1,int NN2)
             {
double dens=0.0,rho1=0.0,nu=0.0;
double g1=0.0,g2=0.0,a1=0.0,a2=0.0,b1,b2=0.0;
switch(model) // Correlation functions are in alphabetical order
    {

/******* models on the real line ****************/
    case 1: // gaussian
      rho1=(1-nuis[0])*rho;
      b1=(z1-mu1)/sqrt(nuis[1]);b2=(z2-mu2)/sqrt(nuis[1]);
      a1=pnorm(b1,0,1,1,0);  a2=pnorm(b2,0,1,1,0); 
      g1=dnorm(b1,0,1,0)/sqrt(nuis[1]); //marginal 1
      g2=dnorm(b2,0,1,0)/sqrt(nuis[1]); //marginal 2
    break;
 case 34:  // Tukeyh
      rho1=(1-nuis[0])*rho;  
      b1=1;
      b2=1;
      a1=1;
      a2=1;
      g1=1;
      g2=1;
   break;
    case 25:  // Logistic
      rho1=(1-nuis[0])*rho;
      b1=1;b2=1;
      a1=plogis(z1, mu1, sqrt(nuis[1]),1,0);
      a2=plogis(z2, mu2, sqrt(nuis[1]),1,0);
      g1=dlogis(z1, mu1, sqrt(nuis[1]),0);
      g2=dlogis(z2, mu2, sqrt(nuis[1]),0);
   break;
   case 12: // t student
    //Rprintf("%f %f \n",1/nuis[0],nuis[1]);
      rho1=(1-nuis[1])*rho; 
      b1=(z1-mu1)/sqrt(nuis[2]);b2=(z2-mu2)/sqrt(nuis[2]);
      a1=pt(b1,1/nuis[0],1,0);  a2=pt(b2,1/nuis[0],1,0);
      g1=dt(b1,1/nuis[0],0)/sqrt(nuis[2]); //marginal 1
      g2=dt(b2,1/nuis[0],0)/sqrt(nuis[2]); //marginal 2
   break;
/******* models on the  positive real line ****************/
     case 24: // lognormal
      rho1=(1-nuis[0])*rho; 
      b1=mu1-nuis[1]/2; b2=mu2-nuis[1]/2;
      a1=plnorm(z1,b1, sqrt(nuis[1]),1,0);
      a2=plnorm(z2,b2, sqrt(nuis[1]),1,0);
      g1=dlnorm(z1,b1, sqrt(nuis[1]),0);
      g2=dlnorm(z2,b2, sqrt(nuis[1]),0);
  break;
     case 22: // loglogistic
      rho1=(1-nuis[0])*rho; 
      b1=1;
      b2=1;
      a1=1;
      a2=1;
      g1=1;
      g2=1;
  break;
   case 21: // gamma
      rho1=(1-nuis[0])*rho;
      b1=1;
      b2=1;
      a1=pgamma(z1,nuis[2]/2,R_pow(nuis[2]/(2*exp(mu1)),-1),1,0);  //revisar
      a2=pgamma(z2,nuis[2]/2,R_pow(nuis[2]/(2*exp(mu2)),-1),1,0);  //revisar
      g1=dgamma(z1,nuis[2]/2,R_pow(nuis[2]/(2*exp(mu1)),-1),0);  //revisar
      g2=dgamma(z2,nuis[2]/2,R_pow(nuis[2]/(2*exp(mu2)),-1),0);  //revisar
   break;
     case 26: // weibull
      rho1=(1-nuis[0])*rho;      
      b1=1;
      b2=1;
      //Rprintf("%f\n",nuis[2]);
      a1=pweibull(z1,nuis[2],exp(mu1)/(gammafn(1+1/nuis[2])),1,0);
      a2=pweibull(z2,nuis[2],exp(mu2)/(gammafn(1+1/nuis[2])),1,0);
      g1=dweibull(z1,nuis[2],exp(mu1)/(gammafn(1+1/nuis[2])),0);
      g2=dweibull(z2,nuis[2],exp(mu2)/(gammafn(1+1/nuis[2])),0);
   break;
/******* models on the a bounded support ****************/
   case 28: //Beta
      rho1=(1-nuis[0])*rho; 
      b1=(z1- nuis[4])/(nuis[5]-nuis[4]); b2=(z2- nuis[4])/(nuis[5]-nuis[4]);
      a1=pbeta(b1,nuis[2],nuis[3],0,0); a2=pbeta(b2,nuis[2],nuis[3],0,0);
      g1=dbeta(b1,nuis[2],nuis[3],0)/(nuis[5]-nuis[4]);//marginal 1
      g2=dbeta(b2,nuis[2],nuis[3],0)/(nuis[5]-nuis[4]);//marginal 2
   break;   
   case 50:  // Beta regression
    rho1=(1-nuis[0])*rho; 
      mu1=1/(1+exp(-mu1));
      mu2=1/(1+exp(-mu2));
    //  Rprintf("%f %f %f %f %f \n",mu1,rho,nuis[2],nuis[4],nuis[3]);
      b1=(z1- nuis[3])/(nuis[4]-nuis[3]); b2=(z2- nuis[3])/(nuis[4]-nuis[3]);
      a1=pbeta(b1, nuis[2]*mu1,(1-mu1)*nuis[2],0,0);
      a2=pbeta(b2, nuis[2]*mu2,(1-mu2)*nuis[2],0,0);
      g1=dbeta(b1, nuis[2]*mu1,(1-mu1)*nuis[2],0)/(nuis[4]-nuis[3]); //marginal 1
      g2=dbeta(b2, nuis[2]*mu2,(1-mu2)*nuis[2],0)/(nuis[4]-nuis[3]); //marginal 2
   break;
   case 33:  // kuma
    rho1=(1-nuis[0])*rho; 
      b1=(z1- nuis[4])/(nuis[5]-nuis[4]);
      b2=(z2- nuis[4])/(nuis[5]-nuis[4]);
      a1= cdf_kuma(b1,nuis[2],nuis[3]);
      a2= cdf_kuma(b2,nuis[2],nuis[3]);
      g1= pdf_kuma(b1,nuis[2],nuis[3])/(nuis[5]-nuis[4]);
      g2= pdf_kuma(b2,nuis[2],nuis[3])/(nuis[5]-nuis[4]);
   break;
   case 42:  // kuma regression
    rho1=(1-nuis[0])*rho; 
      mu1=1/(1+exp(-mu1));mu2=1/(1+exp(-mu2));
      b1=(z1- nuis[3])/(nuis[4]-nuis[3]);
      b2=(z2- nuis[3])/(nuis[4]-nuis[3]);
      mu1=log(1-R_pow(  mu1    ,nuis[2]))/log(0.5) ;
      mu2=log(1-R_pow(  mu2    ,nuis[2]))/log(0.5) ;
      a1= cdf_kuma(b1, nuis[2],1/mu1 );
      a2= cdf_kuma(b2, nuis[2],1/mu2 );
      g1= pdf_kuma(b1, nuis[2],1/mu1 )/(nuis[4]-nuis[3]);
      g2= pdf_kuma(b2, nuis[2],1/mu2 )/(nuis[4]-nuis[3]);
   break;
/************ discrete models *****************************/
case 30: // Poisson
      rho1=(1-nuis[0])*rho; 
      mu1=exp(mu1);mu2=exp(mu2);
      b2=dpois(z2,mu2,0);
      a1=qnorm(ppois(z1,  mu1,1,0),0,1,1,0);
      a2=qnorm(ppois(z1-1,mu1,1,0),0,1,1,0);
      g1=qnorm(ppois(z2,  mu2,1,0),0,1,1,0);      
      g2=qnorm(ppois(z2-1,mu2,1,0),0,1,1,0);  
      if(z1==0) {a2=-99;}
      if(z2==0) {g2=-99;}
   break;
      case 43: // Poisson inflated
      rho1=(1-nuis[0])*rho; 
      mu1=exp(mu1);mu2=exp(mu2);
      b2=dpois(z2,mu2,0);
      a1=qnorm(ppoisinflated(z1,  mu1,nuis[1]),0,1,1,0);
      a2=qnorm(ppoisinflated(z1-1,mu1,nuis[1]),0,1,1,0);
      g1=qnorm(ppoisinflated(z2,  mu2,nuis[1]),0,1,1,0);      
      g2=qnorm(ppoisinflated(z2-1,mu2,nuis[1]),0,1,1,0);  
      if(z1==0) {a2=-99;}
      if(z2==0) {g2=-99;}
   break;
     case 11: // binomial     /// 
      rho1=(1-nuis[0])*rho;
      mu1=pnorm(mu1,0,1,1,0);
      mu2=pnorm(mu2,0,1,1,0);
      b2=dbinom(z2,NN2,mu2,0);
      a1=qnorm(pbinom(z1,  NN1,mu1,1,0),0,1,1,0);
      a2=qnorm(pbinom(z1-1,NN1,mu1,1,0),0,1,1,0);
      g1=qnorm(pbinom(z2,  NN2,mu2,1,0),0,1,1,0);      
      g2=qnorm(pbinom(z2-1,NN2,mu2,1,0),0,1,1,0);
      if((z1<NN1+MAXERR)&&(z1>NN1-MAXERR)) {a1=99;}
      if(z1==0){a2=-99;}  
      if((z2<NN2+MAXERR)&&(z2>NN2-MAXERR)) {g1=99;}
      if(z2==0){g2=-99;}  
   break;
       case 16: // binomial neg     /// ok 
      rho1=(1-nuis[0])*rho;
      mu1=pnorm(mu1,0,1,1,0);
      mu2=pnorm(mu2,0,1,1,0);
      b2=dnbinom(z2,NN2,mu2,0);
      a1=qnorm(pnbinom(z1,  NN1,mu1,1,0),0,1,1,0);
      a2=qnorm(pnbinom(z1-1,NN1,mu1,1,0),0,1,1,0);
      g1=qnorm(pnbinom(z2,  NN2,mu2,1,0),0,1,1,0);      
      g2=qnorm(pnbinom(z2-1,NN2,mu2,1,0),0,1,1,0);
      //Rprintf("%f %f %d %d %f %f  %f\n",mu1,mu2,NN1,NN2,a1,a2,rho1);
      /*if(z1==0) {a2=-99;}
      if(z2==0) {g2=-99;}*/
        //if((z1<NN1+MAXERR)&&(z1>NN1-MAXERR)) {a1=99;}
      if(z1==0){a2=-99;}  
      //if((z2<NN2+MAXERR)&&(z2>NN2-MAXERR)) {g1=99;}
      if(z2==0){g2=-99;}  
   break;
    case 45: // binomial neg  inflated   /// ok 
      rho1=(1-nuis[0])*rho;
      mu1=pnorm(mu1,0,1,1,0);
      mu2=pnorm(mu2,0,1,1,0);
      b2=dnbinom(z2,NN2,mu2,0);
      a1=qnorm(pbneginflated(z1,  NN1,mu1,nuis[1]),0,1,1,0);
      a2=qnorm(pbneginflated(z1-1,NN1,mu1,nuis[1]),0,1,1,0);
      g1=qnorm(pbneginflated(z2,  NN2,mu2,nuis[1]),0,1,1,0);      
      g2=qnorm(pbneginflated(z2-1,NN2,mu2,nuis[1]),0,1,1,0);
      if(z1==0) {a2=-99;}
      if(z2==0) {g2=-99;}
   break;
       }
/*********************** end cases ***************/
/******************copula gaussiana*************************/
if(type_cop==1)  { 
   if(!(model==16||model==11||model==30))   //continous  models
     dens=log(biv_unif_CopulaGauss(a1,a2,rho1)) + log(g1) + log(g2);
    else                            // discrete                  
     dens= log( pbnorm22(a1,g1,rho1) -   pbnorm22(a2,g1,rho1)  - pbnorm22(a1,g2,rho1) + pbnorm22(a2,g2,rho1) );
}
/******************copula clayton and skew *************************/             
if(type_cop==2) 
{
      if(!(model==16||model==11||model==30)) { // continous  models
    if(model==50||model==42) nu=nuis[5];   // for beta2 regression
    if(model==21||model==22||model==26||model==12) nu=nuis[3];   // gammma weibull t  
    if(model==24||model==1||model==25||model==16||model==11||model==30)   nu=nuis[2]; // loggaussian gaussian logistic pois binom binomneg 
    dens= biv_unif_CopulaClayton (a1,a2,rho1,nu)+ log(g1)+log(g2);
               }
     else { Rf_error("not implemented."); } // discrete   models

}
/******************copula skew gaussian *************************/ 
if(type_cop==3) 
{
       if(!(model==16||model==11||model==30)) { // continous  models 
    if(model==50||model==42) nu=nuis[5];   // for beta2 regression
    if(model==21||model==22||model==26||model==12) nu=nuis[3];   // gammma weibull t  
    if(model==24||model==1||model==25||model==16||model==11||model==30)   nu=nuis[2];    // loggaussian gaussian logistic pois binom binomneg 
                    
    dens= log(biv_unif_CopulaSkewGauss(a1,a2,rho1,nu))+ log(g1)+log(g2);
     }
     else { Rf_error("not implemented."); } //   discrete   models
}

if(cond)  {
               if(!(model==16||model==11||model==30))     dens=dens-log(g2);
               else  dens=dens-log(b2);
          }
return(dens);
}

/***********************/
/***********************/
/******* some marginals (log)pdf  *****************/
/***********************/
/***********************/




static const double GL_x[20] = {
    -0.9931285991850949, -0.9639719272779138, -0.9122344282513259, -0.8391169718222188,
    -0.7463319064601508, -0.6360536807265150, -0.5108670019508271, -0.3737060887154195,
    -0.2277858511416451, -0.0765265211334973,  0.0765265211334973,  0.2277858511416451,
     0.3737060887154195,  0.5108670019508271,  0.6360536807265150,  0.7463319064601508,
     0.8391169718222188,  0.9122344282513259,  0.9639719272779138,  0.9931285991850949
};

static const double GL_w[20] = {
    0.0176140071391521, 0.0406014298003869, 0.0626720483341091, 0.0832767415767048,
    0.1019301198172404, 0.1181945319615184, 0.1316886384491766, 0.1420961093183820,
    0.1491729864726037, 0.1527533871307259, 0.1527533871307259, 0.1491729864726037,
    0.1420961093183820, 0.1316886384491766, 0.1181945319615184, 0.1019301198172404,
    0.0832767415767048, 0.0626720483341091, 0.0406014298003869, 0.0176140071391521
};


/*Gauss-Legendre semplice */
double owens_t(double h, double a) {
    if (a == 0.0) return 0.0;

    double sum = 0.0;
    double h_squared = h * h;
    double half_a = a / 2.0;

    for (int i = 0; i < 20; i++) {
        // Transform from [-1,1] to [0,a]
        double x = half_a * (GL_x[i] + 1.0);
        double denom = 1.0 + x * x;
        double fx = exp(-0.5 * h_squared * denom) / denom;
        sum += GL_w[i] * fx;
    }
    return (half_a * sum) / (2.0 * M_PI);
}




/*
double owens_t(double h, double a) {
    if (a == 0.0) return 0.0;
    int n = 900; 
    double sum = 0.0;
    double dx = a / n;
    double h_squared = h * h;
    for (int i = 0; i < n; i++) {
        double x1 = i * dx;
        double x2 = (i + 1) * dx;
        double fx1 = exp(-0.5 * h_squared * (1 + x1 * x1)) / (1 + x1 * x1);
        double fx2 = exp(-0.5 * h_squared * (1 + x2 * x2)) / (1 + x2 * x2);
        sum += (fx1 + fx2) * dx / 2.0;
    }
    return sum / (2 * M_PI);
}*/


double psn(double x,  double omega, double alpha, double tau) {
    double z = x / omega;
    return pnorm(z, 0.0, 1.0, 1, 0) - 2 * owens_t(z, alpha);
}
double dsn(double x, double omega, double alpha, double tau) {
    double z = x/ omega;
    double t = alpha * z + tau;
    return (2.0 / omega) * dnorm(z, 0.0, 1.0, 0) * pnorm(t, 0.0, 1.0, 1, 0);
}



double qsn(double p, double omega, double alpha, double tau, double tol) {
    if (omega <= 0) Rf_error("omega must be positive.");
    if (p <= 0.0) return -INFINITY;
    if (p >= 1.0) return INFINITY;

    int lower_tail = 1, log_p = 0;

    // Calcolo limiti iniziali basati sulla chi-quadro
    double max_q = sqrt(qchisq(p, 1, lower_tail, log_p)) + tau;
    double min_q = -sqrt(qchisq(1 - p, 1, lower_tail, log_p)) + tau;

    // Caso speciale tau == 0 e alpha infinito
    if (tau == 0 && isinf(alpha)) {
        return omega * (alpha > 0 ? max_q : min_q);
    }

    double abs_alpha = fabs(alpha);
    double p_adj = (alpha < 0) ? (1 - p) : p;

    // Inizializza estremi per la ricerca (bisezione + regula falsi alternata)
    double xa = qnorm(p_adj, 0.0, 1.0, 1, 0);
    double xb = sqrt(qchisq(p_adj, 1, 1, 0)) + fabs(tau);

    double fa = psn(xa, 1.0, abs_alpha, tau) - p_adj;
    double fb = psn(xb, 1.0, abs_alpha, tau) - p_adj;

    // Se i valori iniziali sono già abbastanza vicini, ritorna subito
    if (fabs(fa) < tol) {
        return omega * ((alpha < 0) ? -xa : xa);
    }
    if (fabs(fb) < tol) {
        return omega * ((alpha < 0) ? -xb : xb);
    }

    double xc = 0.0, fc = 0.0;
    int use_regula = 0;
    int iter = 0, max_iter = 100;

    while (iter < max_iter) {
        if (use_regula && (fb - fa) != 0.0) {
            // Regula falsi step
            xc = xb - fb * (xb - xa) / (fb - fa);
        } else {
            // Bisezione step
            xc = 0.5 * (xa + xb);
        }

        fc = psn(xc, 1.0, abs_alpha, tau) - p_adj;

        if (fabs(fc) < tol) {
            break;  // Convergenza raggiunta
        }

        // Aggiorna gli estremi dell'intervallo
        if (fc * fa < 0) {
            xb = xc;
            fb = fc;
        } else {
            xa = xc;
            fa = fc;
        }

        use_regula = !use_regula;
        iter++;
    }

    double x = (alpha < 0) ? -xc : xc;
    return omega * x;
}




double one_log_SkewGauss(double z, double m, double vari, double skew)
{
    const double LOG_2 = 0.6931471805599453; // log(2)
    double q = z - m;
    double skew2 = skew * skew;
    double omega = skew2 + vari;
    double sqrt_omega = sqrt(omega);
    double sqrt_var_omega = sqrt(vari * omega);

    double z1 = q / sqrt_omega;
    double z2 = skew * q / sqrt_var_omega;

    double res = LOG_2 - 0.5 * log(omega)
                 + dnorm(z1, 0.0, 1.0, 1)
                 + pnorm(z2, 0.0, 1.0, 1, 1);

    return res;
}




double one_log_tukeyh(double z, double m, double sill, double tail)
{
    const double sqrt_sill = sqrt(sill);
    const double q = (z - m) / sqrt_sill;
    const double q2 = q * q;
    const double x = inverse_lamb(q, tail);  // assume questa è indipendente da extra
    const double W = LambertW(tail * q2);
    const double log_extra = -log(1.0 + W);
    const double log_dnorm = -0.5 * M_LN_2PI - 0.5 * x * x;
    const double log_density = log_dnorm + log(fabs(x)) + log_extra - log(fabs(q)) - log(sqrt_sill);

    return log_density;
}


double one_log_tukeyhh(double z,double m, double sill, double h1,double h2)
{
  double  res=0.0;
 if(z>=m){
    res=one_log_tukeyh(z,m,sill,h2);
          }
  if(z<m){
    res=one_log_tukeyh(z,m,sill,h1);
         }
  return(res);
}
double one_log_T(double z,double m, double sill, double df)
{
  double  res=0.0;
  double q=(z-m)/sqrt(sill);
    res=lgammafn(0.5*(df+1))-(0.5*(df+1))*log1p(q*q/df)-log(sqrt(M_PI*df))-lgammafn(df/2)-0.5*log(sill);
  return(res);
}

double one_log_sas(double z,double m, double skew, double tail,  double vari)
{
  double  res=0.0,S,b;
  double x=(z-m)/(sqrt(vari));
    b=tail*asinh(x)-skew;
    S=sinh(b);
       res= tail*sqrt(1+S*S)*dnorm(S,0,1,0)/(sqrt(1+x*x)*sqrt(vari));
  return(log(res));
}

double one_log_beta(double z, double shape1,double shape2,double min,double  max)
{
  double  res=0.0;
  double q=(z-min)/(max-min);
  res=(shape1/2-1)*log(q)+(shape2/2-1)*log1p(-q)+lgammafn(0.5*(shape1+shape2))-lgammafn(shape1/2)-lgammafn(shape2/2)-log(max-min);
  return(res);
}

double one_log_kumma2(double z,double m, double shape1,double shape2,double min,double  max)
{
  double  res=0.0,k;
  double q=(z-min)/(max-min);k=1-R_pow(q,shape2);
  double m1=1/(1+exp(-m));
  double shapei=log(0.5)/log1p(-R_pow(m1,shape2));
  res=log(shapei)+log(shape2)+(shape2-1)*log(q)+(shapei-1)*log(k)-log(max-min);
  return(res);
}

double one_log_kumma(double z,double m, double shape1,double shape2,double min,double  max)
{
  double  res=0.0,k;
  double q=(z-min)/(max-min);k=1-R_pow(q,shape2);
  res=log(shape1)+log(shape2)+(shape2-1)*log(q)+(shape1-1)*log(k)-log(max-min);
  return(res);
}

double one_log_loggaussian(double z,double m, double sill)
{
  double  res=0.0;
  double q=z*exp(sill/2);
  res=-0.5*R_pow((log(q)-m),2)/sill-log(q)-log(sqrt(sill))-0.5*log(2*M_PI)+sill/2;
  return(res);
}
double one_log_weibull(double z,double m, double shape)
{
  double  res=0.0;
  double scale1=exp(m)/(gammafn(1+1/shape));
  res=log(shape) -log(scale1)  +  (shape-1)*(log(z)-log(scale1))   -    R_pow(z/scale1,shape);
  return(res);
}
double one_log_gamma(double z,double m, double shape)
{
  double  res=0.0;
  res=(shape/2)*log(shape/(2*exp(m)))+(shape/2-1)*log(z)-(shape/(2*exp(m)))*z-log(gammafn(shape/2));
  return(res);
}
/************************************************************************/
double one_log_two_pieceTukey(double z,double m, double sill,double tail, double eta)
{
  double  res=0.0;
  double y=(z-m)/sqrt(sill);
 if(y>=0)res=one_log_tukeyh(y/(1-eta),0,1,tail);       
 if(y<0) res=one_log_tukeyh(y/(1+eta),0,1,tail); 
  return(res-log(sqrt(sill)));
}
double one_log_two_pieceT(double z,double m, double sill, double df, double eta)
{
  double  res=0.0;
  double y=(z-m)/sqrt(sill);
 if(y>=0) res=one_log_T(y/(1-eta),0,1,df);       
 if(y<0)  res=one_log_T(y/(1+eta),0,1,df); 
  return(res-log(sqrt(sill)));
}
double one_log_two_pieceGauss(double z,double m, double sill, double eta)
{
  double  res=0.0;
  double y=(z-m)/sqrt(sill);
 if(y>=0) res=dnorm(y/(1-eta),0,1,1);       
 if(y<0)  res=dnorm(y/(1+eta),0,1,1); 
  return(res-log(sqrt(sill)));
}

double one_log_wrapped(double alfa,double u,double mi,double sill)
{
    double x,wrap_gauss=0.0; 
    x=u-2*atan(mi)-M_PI; 
    double k1=-alfa;

     while(k1<=alfa){
      wrap_gauss +=  dnorm((x+2*k1*M_PI),0,sqrt(sill),0); 
     
     k1 = k1+1;
      }
return(log(wrap_gauss));
}

/************************************************************************/
double one_log_gammagem(double z,double shape,double n)
{
  double  res=0.0;
  res=(shape/2)*log(n/2)+(shape/2-1)*log(z)-0.5*n*z-lgammafn(shape/2);
  return(res);
}

double one_log_bomidal(double z,double m, double sill,double nu,double delta, double eta)
{
  double  res=0.0;
  double q=(z-m)/sqrt(sill);
  double alpha=2*(delta+1)/nu;
  double nn=R_pow(2,1-alpha/2);
 if(z>=m){
    res=log(alpha)+(alpha-1)*log(q)-(alpha-1)*log1p(-eta)-log(2)+one_log_gammagem(R_pow(z/(1-eta),alpha),nu,nn)-0.5*log(sill);    
          }
  if(z<m){
    res=log(alpha)+(alpha-1)*log(-q)-(alpha-1)*log1p(eta)-log(2)+one_log_gammagem(R_pow(-z/(1+eta),alpha),nu,nn)-0.5*log(sill);
         }
  return(res);
}
double one_log_BinomnegZIP(int z,double n, double mu, double mup)
{
  double  res=0.0;
  double  p=pnorm(mup,0,1,1,0);
  double  pp=pnorm(mu,0,1,1,0);
 if(z==0){
    res=log(p+(1-p)*dnbinom(0,n,pp,0));    
          }
  if(z>0){
   res=log1p(-p)+ dnbinom(z,n,pp,1);
         }
  return(res);
}
double one_log_PoisZIP(int z,double lambda, double mup)
{
  double  res=0.0;
  double  p=pnorm(mup,0,1,1,0);
 if(z==0){
    res=log(p+(1-p)*dpois(0,lambda,0));    
          }
  if(z>0){
    res=log1p(-p)+dpois(z,lambda,1);   
         }
  return(res);
}

double one_log_PoisgammaZIP(int z,double lambda, double mup,double shape)
{
  double  res=0.0;
  double  p=pnorm(mup,0,1,1,0);
  double pp=lambda/(lambda+shape);
 if(z==0){
    res=log(p+(1-p)*dnbinom(0, shape, pp,0));
          }
  if(z>0){
    res=log1p(-p)+dnbinom(z, shape, pp,1);
         }
  return(res);
}

double one_log_dpoisgamma(int z,double lambda, double a)
{
  double  res;
  double beta= a/lambda;
   res= z*log(1/(1+beta))+a*log(beta/(1+beta))+lgammafn(a+z)-lgammafn(z+1)-lgammafn(a);
  return(res);
}


double one_log_negbinom_marg(int u,int N, double p)
{
  double  res;
    res=lgammafn(u+N)-(lgammafn(u+1)+lgammafn(N))+N*log(p)+u*log1p(-p);
  return(res);
}


double one_log_loglogistic(double z,double m, double shape)
{
  double  res;
  double c=gammafn(1+1/shape)*gammafn(1-1/shape);
  double k=R_pow(c*z/m,shape)+1;
  //res=log(c)+log(shape/m)+(shape-1)*log(c*z/m)-2*loh(k);
  res=log((c*shape/m)*R_pow((c*z/m),shape-1)*R_pow(k,-2));
  return(res);
}
double one_log_logistic(double z,double m, double sill)
{
  double  res;
  double k=exp((z-m)/sqrt(sill));
  res=log(k)-2*log(k+1)-0.5*log(sill);
  return(res);
}


