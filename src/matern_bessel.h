#ifndef GEOMODELS_MATERN_BESSEL_H
#define GEOMODELS_MATERN_BESSEL_H

/*
 * Thread-safe scaled modified Bessel K for the Matérn OpenMP kernel.
 *
 * work must have length at least 1 + floor(fabs(nu)).
 * On success, *value = exp(x) K_nu(x) and the function returns 1.
 * No R API, warning, allocation or global mutable state is used.
 */
int GeoModels_bessel_k_scaled(double x, double nu, double *work, double *value);

#endif
