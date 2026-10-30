#ifndef GEOMODELS_FIXED_MARGIN_ST_H
#define GEOMODELS_FIXED_MARGIN_ST_H
/* Models 85 and 89: normalized interaction, fixed marginal covariances. */
int GeoFixedMarginSTDimension(void);
int GeoFixedMarginSTValid(int code, const double *par, int spatial_dimension);
double GeoFixedMarginSTCorr(int code, double h, double u,
                            const double *par, int spatial_dimension);
#endif
