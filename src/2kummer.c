#include "header.h"

double hypU_wrap(double a, double b, double x) {
  double out;
  int md; /* method code --- not returned */
  int isfer = 0;

    F77_CALL(chgu)(&a, &b, &x, &out, &md, &isfer);
  if (out == 1e300) {
      out = -1;
  }
  if (isfer == 6) {
    out = -1;
  } else if (isfer != 0) {
    out = -1;
  }
  return out;
}

void hyperg_U_e_call( double *a,  double *b,  double *x, double *val)
{
    *val = hypU_wrap(*a, *b, *x);
}

double kummer(double a,double b, double c)
{
    double res=0.0;
    res=hypU_wrap(a,b,c);
    return (res);
}
