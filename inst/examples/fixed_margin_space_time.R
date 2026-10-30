## Examples for the replaced space-time families (85 and 89).
## Existing parameter names are unchanged; sep now denotes rho in [-1,1].
library(GeoModels)

param_gw <- list(mean=0, sill=1, nugget=0,
                 scale_s=1, scale_t=1,
                 power2_s=4, power2_t=4,
                 smooth_s=0, smooth_t=0, sep=0.7)
param_matern <- list(mean=0, sill=1, nugget=0,
                     scale_s=0.3, scale_t=0.5,
                     smooth_s=0.5, smooth_t=0.5, sep=0.7)

## No new arguments: standalone Euclidean lag-only curves use dimension 2.
## Coordinate-based calls infer dimension from their existing coordinates.
gw_curve <- GeoCorrFct(x=seq(0,1.4,length.out=101), t=c(0,.25,.5,.75),
                       corrmodel="Genwend_Genwend_nosep",param=param_gw)
mat_curve <- GeoCorrFct(x=seq(0,1.4,length.out=101), t=c(0,.25,.5,.75),
                        corrmodel="Matern_Matern_nosep",param=param_matern)

## Matrices, including the exact sparse support pattern of the GW family.
set.seed(42)
coords <- cbind(runif(20,0,2),runif(20,0,2))
coordt <- c(0,.2,.5,.9,1.4)
K_gw <- GeoCovmatrix(coordx=coords,coordt=coordt,
                     corrmodel="Genwend_Genwend_nosep",param=param_gw)
K_sparse <- GeoCovmatrix(coordx=coords,coordt=coordt,
                         corrmodel="Genwend_Genwend_nosep",param=param_gw,sparse=TRUE)
stopifnot(isTRUE(all.equal(as.matrix(K_sparse$covmatrix),as.matrix(K_gw$covmatrix),
                           tolerance=1e-8,check.attributes=FALSE)))

## Simulation and a minimal exact-ML example: estimate sep only first.
## The data layout is time rows x site columns.
y <- as.numeric(t(chol(as.matrix(K_gw$covmatrix))) %*%
                  rnorm(nrow(coords)*length(coordt)))
data <- matrix(y,nrow=length(coordt),byrow=TRUE)
fixed_gw <- param_gw;fixed_gw$sep <- NULL
fit <- GeoFit(data=data,coordx=coords,coordt=coordt,
              corrmodel="Genwend_Genwend_nosep",model="Gaussian",
              fixed=fixed_gw,start=list(sep=0),
              lower=list(sep=-1),upper=list(sep=1),
              likelihood="Full",type="Standard",optimizer="optimize",varest=FALSE)
print(fit)

## With all remaining parameters fixed, the plug-in predictor uses the same kernel.
pred <- GeoKrig(estobj=fit,loc=matrix(c(.6,.8),nrow=1),time=.65,mse=TRUE)
print(pred$pred)

## Matern: swap corrmodel and the corresponding fixed parameter list.
## In a first analysis keep smoothness (and GW power2_s/power2_t) fixed.
## GeoFit supports the existing joint optimizers; for GW the shape inequalities
## must hold at the initial point and all valid optimizer trial points.
