/*
 * Tricomi U for the real, positive-argument Kummer covariance in GeoModels.
 * Production entry point: kummer_ch, used by CorKummer.
 * Legacy U interfaces are preserved separately in the development archive.
 *
 * U(a,b,x) = int_0^inf exp(-x*t) t^(a-1) (1+t)^(b-a-1) dt / Gamma(a).
 * The covariance uses a=power2>0, b=1-smooth<1, x=lag^2/(2*scale^2).
 *
 * Evaluation paths:
 *   1. cancellation-checked small-x connection/logarithmic series;
 *   2. large-x expansion with first-omitted-term bound (a>0, a-b+1>0);
 *   3. saddle-scaled Gauss-Laguerre rules, cached by a (16/32/64 nodes);
 *   4. positive, log-scaled adaptive Gauss-Kronrod fallback.
 * No quantization/interpolation of parameters or distances. No clipping to [0,1].
 * No R_alloc/QUADPACK calls in the fallback. Fixed-size automatic workspaces.
 *
 * References: NIST DLMF 13.2.9, 13.2.11, 13.2.42, 13.4.4, 13.7.3.
 * All quadrature error estimates are numerical estimates, not certified bounds.
 * Compile without -ffast-math: finite checks and cancellation tests are essential.
 * Define KUMMER_STANDALONE for standalone libm tests (not for the R package).
 */
#ifdef KUMMER_STANDALONE
#include <math.h>
#define R_FINITE isfinite
#define R_NaN NAN
#define lgammafn lgamma
#else
#include "header.h"
#endif
#include <float.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif
#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#define KUMMER_TLS _Thread_local
#define KUMMER_HAVE_TLS 1
#elif defined(__GNUC__) || defined(__clang__)
#define KUMMER_TLS __thread
#define KUMMER_HAVE_TLS 1
#else
#define KUMMER_TLS
#define KUMMER_HAVE_TLS 0
#endif

typedef struct {
    double rtol;
    int max_series, max_asympt, max_quad;
} KummerConfig;
/* Same balanced settings used by kummer_ch before the cleanup. */
static const KummerConfig km_balanced_config = {1e-10, 256, 128, 384};
#define KM_QUAD_MAX 512
#define KM_PAR_CACHE 4

typedef struct { double x[64], lw[64]; int ready, n; } KummerRule;
typedef struct {
    int valid, series_ok, integer_m, norm_ready;
    double a,b,c,lga,logA,logB,signA,signB,coefferr;
    double int_logA,int_logB,int_h0,nu_norm,lognorm;
    KummerRule rules[3];
} KummerParams;
#if KUMMER_HAVE_TLS
static KUMMER_TLS KummerParams km_params[KM_PAR_CACHE];
static KUMMER_TLS unsigned km_par_next;
static KUMMER_TLS unsigned km_par_last;
#endif

/* A result contains log(U), never the unnormalized quadrature integral. */
typedef struct { double logu,relerr; int ok,method; } KummerResult;
static KummerResult km_bad(void) {
    KummerResult r={R_NaN,INFINITY,0,0}; return r;
}
static KummerResult km_result(double logu,double err,int method) {
    KummerResult r={logu,err, R_FINITE(logu),method}; return r;
}

static double km_sinpi(double x) { return sin(M_PI*remainder(x,2.0)); }
static double km_gsign(double x) { return x>0.0?1.0:copysign(1.0,km_sinpi(x)); }
static void km_add(double x,double *sum,double *correction) {
    double y=x-*correction, t=*sum+y;
    *correction=(t-*sum)-y; *sum=t;
}
static double km_softplus(double x) {
    return x>0.0?x+log1p(exp(-x)):log1p(exp(x));
}
#ifdef KUMMER_STANDALONE
/* Standalone-test adapter only. Production uses Rmath digamma(). */
static double km_psi(double x) {
    double s=0.0;
    while(x<16.0){s-=1.0/x;x+=1.0;}
    double v=1.0/(x*x);
    return s+log(x)-0.5/x-v*(1.0/12.0-v*(1.0/120.0-v*(1.0/252.0-
           v*(1.0/240.0-v*(1.0/132.0-v*(691.0/32760.0-v*(1.0/12.0)))))));
}
#else
#define km_psi digamma
#endif
static void km_init_params(KummerParams *p,double a,double b) {
    memset(p,0,sizeof(*p));
    p->valid=1; p->a=a; p->b=b; p->c=a+(1.0-b);
    p->lga=lgammafn(a);p->integer_m=-1;
    if(b<=1.0 && b==floor(b) && 1.0-b<=64.0) {
        int m=(int)(1.0-b);double hm=0.0;
        for(int j=1;j<=m;j++)hm+=1.0/j;
        p->integer_m=m;
        p->int_logA=m>0?lgammafn((double)m)-lgammafn(a+m):-INFINITY;
        p->int_logB=-lgammafn(m+1.0)-p->lga;
        p->int_h0=km_psi(a+m)+1.1544313298030657212130241801648-hm;
    }
    /* Near integer b the connection formula subtracts singular terms. */
    if (fabs(b-nearbyint(b))>1e-5 && p->c>0.0 && fabs(b)<1000.0) {
        p->logA=lgammafn(1.0-b)-lgammafn(p->c);
        p->logB=lgammafn(b-1.0)-p->lga;
        p->signA=km_gsign(1.0-b); p->signB=km_gsign(b-1.0);
        p->series_ok=R_FINITE(p->logA)&&R_FINITE(p->logB);
        p->coefferr=8.0*DBL_EPSILON*(1.0+fabs(p->logA)+fabs(p->logB)+2.0*fabs(p->lga));
    }
}
static KummerParams *km_get_params(double a,double b,KummerParams *local) {
#if KUMMER_HAVE_TLS
    KummerParams *p=&km_params[km_par_last];
    if(p->valid && p->a==a && p->b==b) return p;
    for(unsigned i=0;i<KM_PAR_CACHE;i++) {
        p=&km_params[i];
        if(p->valid && p->a==a && p->b==b) {km_par_last=i;return p;}
    }
    unsigned i=km_par_next; km_par_next=(i+1)%KM_PAR_CACHE;km_par_last=i;
    p=&km_params[i]; km_init_params(p,a,b); (void)local; return p;
#else
    km_init_params(local,a,b);return local;
#endif
}

/* DLMF 13.2.42, evaluated with scaled, compensated sums. */
static KummerResult km_series(const KummerParams *p,double x,const KummerConfig *cfg) {
    if(!p->series_ok || x>4.0 || (p->a+p->c)*x>30.0)return km_bad();
    double lb=p->logB+(1.0-p->b)*log(x), scale=fmax(p->logA,lb);
    if(!R_FINITE(scale))return km_bad();
    double t1=p->signA*exp(p->logA-scale),t2=p->signB*exp(lb-scale);
    double sum=0.0,comp=0.0,abssum=fabs(t1)+fabs(t2);
    km_add(t1,&sum,&comp);km_add(t2,&sum,&comp);
    for(int k=1;k<=cfg->max_series;k++) {
        double d1=p->b+k-1.0,d2=2.0-p->b+k-1.0;
        if(d1==0.0||d2==0.0)return km_bad();
        t1*=x*(p->a+k-1.0)/(k*d1);
        t2*=x*(p->c+k-1.0)/(k*d2);
        if(!R_FINITE(t1)||!R_FINITE(t2))return km_bad();
        km_add(t1,&sum,&comp);km_add(t2,&sum,&comp);
        abssum+=fabs(t1)+fabs(t2);
        if(p->b+k>0.0 && 2.0-p->b+k>0.0) {
            double q1=x/(k+1.0)*fmax(1.0,(p->a+k)/(p->b+k));
            double q2=x/(k+1.0)*fmax(1.0,(p->c+k)/(2.0-p->b+k));
            if(q1<0.8 && q2<0.8 && sum>0.0) {
                double tail=fabs(t1)*q1/(1.0-q1)+fabs(t2)*q2/(1.0-q2);
                double rounding=((16.0+8.0*k)*DBL_EPSILON+p->coefferr)*abssum;
                double err=(tail+rounding)/sum;
                if(err<=cfg->rtol*0.25)
                    return km_result(scale+log(sum),err,1);
                if(tail<rounding && rounding>cfg->rtol*0.25*sum)return km_bad();
            }
        }
    }
    return km_bad();
}

/* DLMF 13.2.9 and 13.2.11, specialized to b=1-m (m>=0).
 * A finite polynomial plus a logarithmic series. Digamma is cached once
 * per parameter pair; successive values use reciprocal recurrences. */
static KummerResult km_integer_series(const KummerParams *p,double x,const KummerConfig *cfg) {
    int m=p->integer_m;
    if(m<0 || x>4.0 || (p->a+m)*x>20.0)return km_bad();
    double lx=log(x),lb=p->int_logB+m*lx,scale=fmax(p->int_logA,lb);
    if(!R_FINITE(scale))return km_bad();
    double sum=0.0,comp=0.0,abssum=0.0;
    if(m>0) {
        double f=exp(p->int_logA-scale);km_add(f,&sum,&comp);abssum+=fabs(f);
        for(int j=1;j<m;j++) {
            f*=-x*(p->a+j-1.0)/(j*(m-j));
            km_add(f,&sum,&comp);abssum+=fabs(f);
        }
    }
    double t=exp(lb-scale)*((m&1)?1.0:-1.0),h=lx+p->int_h0,A=p->a+m;
    double ce=8.0*DBL_EPSILON*(1.0+fabs(p->lga)+fabs(lb)+fabs(p->int_h0));
    for(int k=0;k<=cfg->max_series;k++) {
        double v=t*h;
        if(!R_FINITE(v))return km_bad();
        km_add(v,&sum,&comp);abssum+=fabs(t)*(1.0+fabs(h));
        double rat=x/(k+1.0)*fmax(1.0,(A+k)/(m+1.0+k));
        if(rat<0.8 && sum>0.0) {
            double dh=1.0/(A+k)+1.0/(k+1.0)+1.0/(m+1.0+k);
            double tail=fabs(t)*rat*(fabs(h)/(1.0-rat)+dh/((1.0-rat)*(1.0-rat)));
            double rounding=((32.0+8.0*(k+m))*DBL_EPSILON+ce)*abssum;
            double err=(tail+rounding)/sum;
            if(err<=0.25*cfg->rtol)return km_result(scale+log(sum),err,3);
            if(tail<rounding && rounding>0.25*cfg->rtol*sum)return km_bad();
        }
        t*=x*(A+k)/((m+1.0+k)*(k+1.0));
        h+=1.0/(A+k)-1.0/(k+1.0)-1.0/(m+1.0+k);
    }
    return km_bad();
}

/* For a,c>0, Taylor's integral remainder in (1+t/x)^(-c)
 * bounds the scaled-U error by the first omitted 2F0 term. */
static KummerResult km_asympt(const KummerParams *p,double x,const KummerConfig *cfg) {
    double a=p->a,c=p->c;
    if(c<=0.0 || (a/x)*c>0.75)return km_bad();
    double sum=1.0,comp=0.0,term=1.0,abssum=1.0;
    for(int k=1;k<=cfg->max_asympt;k++) {
        double next=-term*((a+k-1.0)/x)*((c+k-1.0)/k);
        double rnd=(8.0+4.0*k)*DBL_EPSILON*abssum;
        if(sum>0.0 && fabs(next)+rnd<=0.25*cfg->rtol*sum)
            return km_result(-a*log(x)+log(sum),(fabs(next)+rnd)/sum,2);
        if(!R_FINITE(next)||fabs(next)>fabs(term))return km_bad();
        km_add(next,&sum,&comp);abssum+=fabs(next);term=next;
    }
    return km_bad();
}

/* Golub-Welsch by implicit QL on the Jacobi tridiagonal matrix. The
 * first eigenvector row gives NORMALIZED Gamma(a,1) quadrature weights.
 * Consequently Gamma(a) is never formed, even when a is large. */
static int km_laguerre_rule(KummerRule *r,int n,double a) {
    double d[64],e[64],z[64];
    for(int i=0;i<n;i++){d[i]=2.0*i+a;e[i]=sqrt((i+1.0)*(i+a));z[i]=0.0;}
    z[0]=1.0;e[n-1]=0.0;
    for(int l=0;l<n;l++) {
        int it=0;
        for(;;) {
            int m=l;
            while(m<n-1 && fabs(e[m])>DBL_EPSILON*(fabs(d[m])+fabs(d[m+1])))m++;
            if(m==l)break;
            if(++it>60)return 0;
            double p=d[l],g=(d[l+1]-p)/(2.0*e[l]),v=hypot(g,1.0);
            g=d[m]-p+e[l]/(g+copysign(v,g));
            double s=1.0,c=1.0; p=0.0;
            for(int i=m-1;i>=l;i--) {
                double f=s*e[i],bb=c*e[i];
                if(fabs(g)<=fabs(f)) {
                    c=g/f;v=hypot(c,1.0);e[i+1]=f*v;s=1.0/v;c*=s;
                } else {
                    s=f/g;v=hypot(s,1.0);e[i+1]=g*v;c=1.0/v;s*=c;
                }
                g=d[i+1]-p;v=(d[i]-g)*s+2.0*c*bb;
                p=s*v;d[i+1]=g+p;g=c*v-bb;
                f=z[i+1];z[i+1]=s*z[i]+c*f;z[i]=c*z[i]-s*f;
            }
            d[l]-=p;e[l]=g;e[m]=0.0;
        }
    }
    double wsum=0.0;
    for(int i=0;i<n;i++) {
        if(!(d[i]>0.0)||!R_FINITE(d[i]))return 0;
        wsum+=z[i]*z[i];
    }
    for(int i=0;i<n;i++){r->x[i]=d[i];r->lw[i]=z[i]==0.0?-INFINITY:2.0*log(fabs(z[i]))-log(wsum);}
    r->n=n;r->ready=1;return 1;
}

static double km_lag_eval(KummerRule *r,double a,double c,double x,double s) {
    double maxv=-INFINITY,sum=0.0,one_minus_xs=1.0-x*s;
    for(int i=0;i<r->n;i++) {
        double v=r->lw[i]+one_minus_xs*r->x[i]-c*log1p(s*r->x[i]);
        if(v>maxv){sum=sum*exp(maxv-v)+1.0;maxv=v;}
        else sum+=exp(v-maxv);
    }
    return a*log(s)+maxv+log(sum);
}
static double km_mode_t(double a,double b,double x) {
    double g=x+(1.0-b),disc=hypot(g,2.0*sqrt(a)*sqrt(x));
    if(g>=0.0)return 2.0*(a/disc)/(1.0+g/disc);
    return (0.5*disc-0.5*g)/x;
}
static KummerResult km_laguerre(KummerParams *p,double x,const KummerConfig *cfg) {
    if(p->c<=0.0 || p->a<0.01 || p->a>256.0)return km_bad();
    double t=km_mode_t(p->a,p->b,x),s=t/p->a;
    if(!(s>0.0)||!R_FINITE(s))return km_bad();
    double last=0.0;
    for(int k=0;k<3;k++) {
        KummerRule *r=&p->rules[k];
        if(!r->ready && !km_laguerre_rule(r,16<<k,p->a))return km_bad();
        double u=km_lag_eval(r,p->a,p->c,x,s);
        if(!R_FINITE(u))return km_bad();
        if(k>0) {
            double err=fabs(expm1(u-last));
            if(err<cfg->rtol*0.03125)
                return km_result(u,err*16.0+32.0*DBL_EPSILON,5);
        }
        last=u;
    }
    return km_bad();
}

/* Positive quadrature of the original integral on y=log(t), scaled at
 * its mode. This is also valid for integer/near-integer b: no subtraction
 * of two large hypergeometric-M terms, no inverse-Gamma underflow. */
typedef struct { double a,c,y0,q0,logq0,p0,soft0,h0; } KummerQuad;
static double km_logheight(double s,const KummerQuad *q) {
    double y=q->y0+s;
    double le=q->logq0+s;
    if(le>log(DBL_MAX))return -INFINITY;
    double xt=exp(le),dh;
    if(fabs(s)<0.5 && q->p0<1.0) {
        double es=expm1(s),v=log1p(q->p0*es);
        dh=q->a*s-q->q0*es-q->c*v;
    } else {
        dh=q->a*s-(xt-q->q0)-q->c*(km_softplus(y)-q->soft0);
    }
    return dh;
}
static double km_height(double s,const KummerQuad *q) {return exp(km_logheight(s,q));}
static double km_slope(double s,const KummerQuad *q) {
    double y=q->y0+s,xt=exp(q->logq0+s);
    double p=y>0?1.0/(1.0+exp(-y)):exp(y)/(1.0+exp(y));
    return q->a-xt-q->c*p;
}

/* 15-point Kronrod / embedded 7-point Gauss rule. */
static void km_gk15(double lo,double hi,const KummerQuad *q,double *val,double *err) {
    static const double x[8]={
        .991455371120812639206854697526329,.949107912342758524526189684047851,
        .864864423359769072789712788640926,.741531185599394439863864773280788,
        .586087235467691130294144838258730,.405845151377397166906606412076961,
        .207784955007898467600689403773245,0};
    static const double wk[8]={
        .022935322010529224963732008058970,.063092092629978553290700663189204,
        .104790010322250183839876322541518,.140653259715525918745189590510238,
        .169004726639267902826583426598550,.190350578064785409913256402421014,
        .204432940075298892414161999234649,.209482141084727828012999174891714};
    static const double wg[4]={.129484966168869693270611432679082,
        .279705391489276667901467771423780,.381830050505118944950369775488975,
        .417959183673469387755102040816327};
    double mid=0.5*(lo+hi),h=0.5*(hi-lo),fc=km_height(mid,q);
    double sk=wk[7]*fc,sg=wg[3]*fc, f1[7],f2[7];
    for(int i=0;i<7;i++) {
        f1[i]=km_height(mid-h*x[i],q);f2[i]=km_height(mid+h*x[i],q);
        double f=f1[i]+f2[i];sk+=wk[i]*f;if(i&1)sg+=wg[(i-1)/2]*f;
    }
    double avg=sk*0.5,asc=wk[7]*fabs(fc-avg);
    for(int i=0;i<7;i++)asc+=wk[i]*(fabs(f1[i]-avg)+fabs(f2[i]-avg));
    *val=sk*h;asc*=h;double e=fabs((sk-sg)*h);
    if(asc>0.0 && e>0.0)e=asc*fmin(1.0,pow(200.0*e/asc,1.5));
    *err=fmax(e,50.0*DBL_EPSILON*fabs(*val));
}
typedef struct {double lo,hi,v,e;} KummerInterval;
static KummerResult km_quad(const KummerParams *p,double x,const KummerConfig *cfg) {
    double t=km_mode_t(p->a,p->b,x);
    if(!(t>0.0)||!R_FINITE(t))return km_bad();
    KummerQuad q;
    q.a=p->a;q.c=p->c;q.y0=log(t);q.q0=x*t;q.logq0=log(x)+q.y0;q.p0=t/(1.0+t);q.soft0=log1p(t);
    q.h0=t>1.0?(p->b-1.0)*q.y0-p->c*log1p(1.0/t)-q.q0:
        p->a*q.y0-p->c*q.soft0-q.q0;
    double curvature=q.q0+q.c*q.p0*(1.0-q.p0);
    if(!(curvature>0.0)||!R_FINITE(curvature))return km_bad();
    double step=fmin(1.0,1.0/sqrt(curvature));
    if(!(step>0.0))return km_bad();
    double total=0.0,error=0.0,tailerr=0.0;
    KummerInterval iv[KM_QUAD_MAX];int n=0;
    /* Build geometrically sized panels from the mode in both directions,
     * so a narrow peak can never fall between all quadrature nodes. */
    for(int side=-1;side<=1;side+=2) {
        double prev=0.0,reach=step;int ended=0;
        for(int k=0;k<96;k++) {
            if(n>=cfg->max_quad)return km_bad();
            double lo=side<0?-reach:prev,hi=side<0?-prev:reach;
            iv[n].lo=lo;iv[n].hi=hi;
            km_gk15(lo,hi,&q,&iv[n].v,&iv[n].e);
            total+=iv[n].v;error+=iv[n].e;n++;
            double s=side*reach,h=km_logheight(s,&q),slope=km_slope(s,&q);
            /* Log-density is concave for c>=0. For c<0 a conservative
             * left-tail derivative lower bound is a-x*t; right-tail
             * negative derivative is bounded below by x*t-a-|c|. */
            double bound_slope;
            if(q.c>=0.0)bound_slope=side<0?slope:-slope;
            else {double xt=exp(q.logq0+s);bound_slope=side<0?q.a-xt:xt-q.a+q.c;}
            double tail=bound_slope>0.0?exp(h)/bound_slope:INFINITY;
            if(tail<=cfg->rtol*0.01*total){tailerr+=tail;ended=1;break;}
            prev=reach;reach*=2.0;
            if(!R_FINITE(reach))return km_bad();
        }
        if(!ended)return km_bad();
    }
    if(!(total>0.0)||!R_FINITE(total))return km_bad();
    while(error+tailerr>cfg->rtol*0.25*total) {
        if(n>=cfg->max_quad)return km_bad();
        int j=0;for(int i=1;i<n;i++)if(iv[i].e>iv[j].e)j=i;
        KummerInterval old=iv[j];double mid=(old.lo+old.hi)*0.5;
        if(mid==old.lo || mid==old.hi)return km_bad();
        iv[j].hi=mid;iv[n].lo=mid;iv[n].hi=old.hi;
        km_gk15(iv[j].lo,mid,&q,&iv[j].v,&iv[j].e);
        km_gk15(mid,iv[n].hi,&q,&iv[n].v,&iv[n].e);
        total+=iv[j].v+iv[n].v-old.v;error+=iv[j].e+iv[n].e-old.e;n++;
        if(!R_FINITE(total)||!R_FINITE(error))return km_bad();
    }
    return km_result(q.h0-p->lga+log(total),(error+tailerr)/total,4);
}

static KummerResult km_eval(KummerParams *p,double x,const KummerConfig *cfg) {
    KummerResult r;
    if(p->c==0.0)return km_result(-p->a*log(x),0.0,2);
    r=km_asympt(p,x,cfg);if(r.ok)return r;
    r=km_series(p,x,cfg);if(r.ok)return r;
    r=km_integer_series(p,x,cfg);if(r.ok)return r;
    r=km_laguerre(p,x,cfg);if(r.ok)return r;
    return km_quad(p,x,cfg);
}


/* Normalized CH covariance used by CorKummer. */
double kummer_ch(double alpha,double nu,double x) {
    if(!R_FINITE(alpha)||!R_FINITE(nu)||!R_FINITE(x)||alpha<=0.0||nu<=0.0||x<0.0)return R_NaN;
    if(x==0.0)return 1.0;
    KummerParams local,*p=km_get_params(alpha,1.0-nu,&local);
    if(!p->norm_ready || p->nu_norm!=nu) {
        p->lognorm=lgammafn(alpha+nu)-lgammafn(nu);
        p->nu_norm=nu;p->norm_ready=1;
    }
    KummerResult r=km_eval(p,x,&km_balanced_config);
    if(!r.ok)return R_NaN;
    double rho=exp(p->lognorm+r.logu);
    return R_FINITE(rho)&&rho>=0.0&&rho<=1.0+1e-9?rho:R_NaN;
}
