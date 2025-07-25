\name{GeoFit}
\alias{GeoFit}
\alias{print.GeoFit}
\encoding{UTF-8}
\title{Max-Likelihood-Based Fitting of Gaussian and non Gaussian random fields.}
\description{
  Maximum weighted composite-likelihood fitting for Gaussian  and some Non-Gaussian  
  univariate spatial, spatio-temporal and bivariate spatial random fieldss
  The function allows to fix any of the parameters and setting upper/lower bound in the optimization.
  Different methods of optimization can be used.

}
\usage{
GeoFit(data, coordx, coordy=NULL,coordz=NULL, coordt=NULL, coordx_dyn=NULL,copula=NULL,
      corrmodel=NULL,distance="Eucl",fixed=NULL,anisopars=NULL,
      est.aniso=c(FALSE,FALSE), grid=FALSE, likelihood='Marginal',
      lower=NULL,maxdist=Inf,neighb=NULL,
      maxtime=Inf, memdist=TRUE,method="cholesky", 
      model='Gaussian',n=1, onlyvar=FALSE ,
      optimizer='Nelder-Mead', 
      radius=1,  sensitivity=FALSE,sparse=FALSE, 
      start=NULL, taper=NULL, tapsep=NULL, 
      type='Pairwise', upper=NULL, varest=FALSE, 
      weighted=FALSE,X=NULL,nosym=FALSE,spobj=NULL,spdata=NULL)
}

\arguments{
  \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation)  or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid)
   or a
   (\eqn{t \times d}{t x d})-matrix (a single spatial-temporal realisation)   or an (\eqn{d \times d \times t \times n }{d x d x t})-array
   (a single spatial-temporal realisation on regular grid).
   For the description see the Section \bold{Details}.}
       \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix  or  (\eqn{d \times 3}{d x 3})-matrix 
     Coordinates on a sphere for a  fixed radius \code{radius} 
    are passed in lon/lat format expressed in decimal degrees.}
  \item{coordy}{A numeric vector giving 1-dimension of
    spatial coordinates; Optional argument, the default is \code{NULL}.}
      \item{coordz}{A numeric vector giving 1-dimension of
    spatial coordinates;  Optional argument, the default is \code{NULL}.}
  \item{coordt}{A numeric vector assigning 1-dimension of
    temporal coordinates.  Optional argument, the default is \code{NULL}
    then a spatial random fields is expected.}
  \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d x 2})-matrices
       containing dynamical (in time) spatial coordinates. Optional argument, the default is \code{NULL}}
  \item{copula}{String; the type of copula. It can be "Clayton" or "Gaussian"}
  \item{corrmodel}{String; the name of a correlation model, for the
    description see the Section \bold{Details}.}
  \item{distance}{String; the name of the spatial distance. The default is \code{Eucl}, the euclidean distance.
    See the Section \bold{Details}.}
  \item{fixed}{An optional named list giving the values of the parameters that
    will be considered as known values. The listed parameters for a
    given correlation function will be not estimated.}
  \item{anisopars}{ A list of two elements: "angle" and "ratio" i.e. the anisotropy angle and the anisotropy ratio, respectively.}
  \item{est.aniso}{ A bivariate logical vector  providing which  anisotropic parameters must be estimated.}

  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as spatial or spatial-temporal realisations on a set
     of non-equispaced spatial sites (irregular grid).}
  \item{likelihood}{String; the configuration of the composite
    likelihood. \code{Marginal} is the default, see the Section \bold{Details}.}
   \item{lower}{An optional named list giving the values  for the lower bound of the space parameter
    when the optimizer is  \code{L-BFGS-B} or \code{nlminb} or \code{bobyqa}  or \code{optimize}. The 
    names of the list must be the same of the names in the   \code{start} list.}
  \item{maxdist}{Numeric; an optional positive value indicating the maximum
    spatial distance considered in the composite 
    computation. See the Section \bold{Details} for more information.}
      \item{neighb}{Numeric; an optional positive integer indicating the 
    order of neighborhood  in the composite likelihood
    computation. See the Section \bold{Details} for more information.}
  \item{maxtime}{Numeric; an optional positive integer indicating the 
    order of temporal neighborhood  in the composite likelihood
    computation.}
        \item{memdist}{ Logical; if \code{TRUE} then  all  the distances  useful in the  composite likelihood 
      estimation are computed before  the optimization. \code{FALSE} is deprecated. }
   \item{method}{String; the type of matrix decomposition used in the simulation. Default is cholesky.
    The other possible choices is \code{svd}.}   
  \item{model}{String; the type of random fields and therefore the densities associated to the likelihood
    objects. \code{Gaussian} is the default, see the Section
    \bold{Details}.}
    \item{n}{Numeric; number of trials in a binomial random fields; number of successes in a negative binomial random fields }  
  \item{onlyvar}{Logical; if \code{TRUE}  (and varest is  \code{TRUE}) only the variance covariance matrix is computed
  without optimizing. \code{FALSE} is the default.}
  \item{optimizer}{String; the optimization algorithm
    (see \code{\link{optim}} for details). \code{Nelder-Mead} is the
    default. Other possible choices are \code{nlm}, \code{BFGS}, \code{SANN}, \code{L-BFGS-B} and \code{nlminb} and
    \code{bobyqa}.
    In these last three cases upper and lower bounds can be passed by the user. In the case
    of one-dimensional optimization, the function \code{optimize} is used.}   
   \item{radius}{Numeric; the radius of the sphere in the case of lon-lat coordinates.
   Default value is 1.}
   \item{sensitivity}{Logical; if \code{TRUE} then  the sensitivy  matrix is computed}
  \item{sparse}{Logical; if \code{TRUE} then  maximum likelihood is computed
  using sparse matrices algorithms (spam packake).It should be used
  with compactly supported covariance models.\code{FALSE} is the default.}
  \item{start}{An optional named list with the initial values of the
    parameters that are used by the numerical routines in maximization
    procedure. \code{NULL} is the default (see \bold{Details}).}
  \item{taper}{String; the name of the type of covariance matrix.
     It can be \code{Standard} (the default value) or
    \code{Tapering} for taperd covariance matrix.}
   \item{tapsep}{Numeric; an optional  value indicating the separabe parameter in the 
     space time adaptive taper (see \bold{Details}).}
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods (see \bold{Details}).}
  \item{upper}{An optional named list giving the values  for the upper bound 
      of the space parameter when the optimizer is or \code{L-BFGS-B} or  \code{bobyqa} or   \code{nlminb} or \code{optimize}. 
       The names of the list must be the same of the names in the   \code{start} list.}  
  \item{varest}{Logical; if \code{TRUE} the estimates' variances and
    standard errors are returned. For composite likelihood estimation it is deprecated.
    Use sensitivity TRUE and update the object using the function \code{GeoVarestbootstrap}
    \code{FALSE} is the default.}
  \item{weighted}{Logical; if \code{TRUE} the likelihood objects are
    weighted, see the Section \bold{Details}. If \code{FALSE} (the
    default) the composite likelihood is not weighted.}
   \item{X}{Numeric; Matrix of spatio(temporal)covariates in the linear mean specification.}
    \item{nosym}{Logical; if TRUE simmetric weights are not considered. This allows a faster but less efficient CL estimation. }
      \item{spobj}{An object of class sp or spacetime}
   \item{spdata}{Character:The name of data in the sp or spacetime object}
}

\details{
  \code{GeoFit} provides  weighted composite likelihood  estimation based on pairs and independence composite likelihood estimation
  for Gaussian and non Gaussian  random fields. Specifically, marginal and conditional  pairwise 
  likelihood are considered for each type of random field (Gaussian and not Gaussian).
  The optimization method is specified using \code{optimizer}. The default method is \code{Nelder-mead} and other available methods are  \code{nlm}, \code{BFGS}, \code{SANN}, \code{L-BFGS-B}, \code{bobyqa}
    and  \code{nlminb}. In the last three cases upper and lower bounds constraints in the optimization can be specified using \code{lower} and \code{upper} parameters.

  Depending on the dimension of \code{data} and on  the name of the correlation model, 
  the observations   are assumed as a realization of
  a spatial, spatio-temporal or bivariate random field.
  Specifically, with \code{data}, \code{coordx}, \code{coordy}, \code{coordt} parameters:
  \itemize{
    \item If \code{data} is a numeric \eqn{d}{d}-dimensional vector, \code{coordx} and \code{coordy} are two
    numeric \eqn{d}{d}-dimensional vectors (or \code{coordx} is (\eqn{d \times 2}{d x 2})-matrix and \code{coordy=NULL}), then the data are interpreted as a single spatial
    realisation observed on \eqn{d}{d} spatial sites;
    \item If \code{data} is a numeric (\eqn{t \times d}{t x d})-matrix, \code{coordx} and \code{coordy} are two
    numeric \eqn{d}{d}-dimensional vectors (or \code{coordx} is (\eqn{d \times 2}{d x 2})-matrix and \code{coordy=NULL}),
    \code{coordt} is a numeric \eqn{t}{t}-dimensional vector, then the data
    are interpreted as a single spatial-temporal realisation of a random fields observed on \code{d} spatial sites and for \code{t} times.
      \item If \code{data} is a numeric (\eqn{2 \times d}{2 x d})-matrix, \code{coordx} and \code{coordy} are two
    numeric \eqn{d}{d}-dimensional vectors (or \code{coordx} is (\eqn{d \times 2}{d x 2})-matrix and \code{coordy=NULL}), then the data
    are interpreted as a single spatial realisation of a bivariate random fields observed on \code{d} spatial sites. 
    \item If \code{data} is a list,  \code{coordxdyn}  is a list and   \code{coordt} is a numeric \eqn{t}{t}-dimensional vector, then the data
    are interpreted as a single spatial-temporal realisation of a random fields observed on dynamical spatial sites (different locations sites for each temporal instants) and for \code{t} times.
 }

Is is also possible to specify a matrix of covariates using \code{X}.
Specifically:
  \itemize{
    \item  In the spatial case  \code{X} must be a (\eqn{d \times k}{d x k}) covariates matrix associated to  \code{data}  a numeric \eqn{d}{d}-dimensional vector;
     \item  In the spatiotemporal case  \code{X} must be a (\eqn{N \times k}{N x k}) covariates matrix associated to  \code{data}  a numeric (\eqn{t \times d}{t x d})-matrix, where \eqn{N=t\times d};
     \item  In the spatiotemporal case  \code{X} must be a (\eqn{N \times k}{N x k}) covariates matrix associated to  \code{data}  a numeric (\eqn{t \times d}{t x d})-matrix, where \eqn{N=2\times d};
 }

  The \code{corrmodel} parameter allows to select a specific correlation
  function for the random fields. (See   \code{\link{GeoCovmatrix}} ).

  The \code{distance} parameter allows to consider differents kinds of spatial distances.
  The settings alternatives are:
  \enumerate{
    \item \code{Eucl}, the euclidean distance (default value);
    \item \code{Chor}, the chordal distance;
    \item \code{Geod}, the geodesic distance;
  }

  The \code{likelihood} parameter represents the composite-likelihood
  configurations. The settings alternatives are:
  \enumerate{
    \item \code{Conditional}, the composite-likelihood is formed by
    conditionals likelihoods;
    \item \code{Marginal}, the composite-likelihood is formed by
  marginals likelihoods (default value);
    \item \code{Full}, the composite-likelihood turns out to be the standard likelihood;
  }

  It must be coupled with the  \code{type} parameter that can be fixed to 
    \enumerate{
    \item \code{Pairwise}, the composite-likelihood is based on pairs;
    \item \code{Independence}, the composite-likelihood is based on indepedence;
    \item \code{Standard}, this is the option for  the standard likelihood;
  }

The possible combinations are:
 \enumerate{
    \item  \code{likelihood="Marginal"} and \code{type="Pairwise"} for  maximum marginal pairwise likelihood estimation (the default setting)
    \item  \code{likelihood="Conditional"} and \code{type="Pairwise"} for maximum  conditional pairwise likelihood estimation
    \item  \code{likelihood="Marginal"} and \code{type="Independence"} for  maximum independence composite likelihood estimation
   \item  \code{likelihood="Full"} and \code{type="Standard"} for  maximum stardard likelihood estimation 
}
The first three combinations can be used for any model. The standard likelihood can be used only for some specific model.

  The \code{model} parameter indicates the type of random field
  considered. The available options are:

  random fields with marginal symmetric distribution:
  \itemize{
    \item \code{Gaussian}, for a Gaussian random field.   
    \item \code{StudentT}, for a StudentT random field (see  Bevilacqua M., Caamaño C., Arellano Valle R.B., Morales-Oñate V., 2020).
    \item \code{Tukeyh}, for a Tukeyh random field.
    \item \code{Tukeyh2}, for a Tukeyhh random field. (see Caamaño et al., 2023)
    \item \code{Logistic}, for a Logistic random field.
 }   
   random fields with positive values and  right skewed marginal distribution:
  \itemize{
    \item \code{Gamma} for a Gamma  random fields (see  Bevilacqua M., Caamano C.,  Gaetan, 2020)
    \item \code{Weibull} for a Weibull  random fields (see  Bevilacqua M., Caamano C.,  Gaetan, 2020)
    \item \code{LogGaussian} for a  LogGaussian  random fields (see  Bevilacqua M., Caamano C.,  Gaetan, 2020)
    \item \code{LogLogistic} for a  LogLogistic random fields. 
}
random fields with with possibly asymmetric marginal distribution:
  \itemize{
    \item \code{SkewGaussian} for a skew Gaussian  random field (see   Alegrıa et al. (2017)).
    \item \code{SinhAsinh} for a Sinh-arcsinh  random field  (see Blasi et. al 2022).
    \item \code{TwopieceGaussian} for a Twopiece Gaussian random field   (see Bevilacqua et. al 2022).
    \item \code{TwopieceTukeyh} for a Twopiece Tukeyh random field   (see Bevilacqua et. al 2022).
 }
 random fields with for directional data
 \itemize{   
    \item \code{Wrapped} for a wrapped Gaussian  random field (see  Alegria A., Bevilacqua, M., Porcu, E. (2016))
    }
random fields with marginal counts data
 \itemize{  
     \item \code{Poisson} for a Poisson  random field (see Morales-Navarrete et. al 2021).
     \item \code{PoissonZIP} for a zero inflated Poisson  random field (see Morales-Navarrete et. al 2021).
     \item \code{Binomial} for a Binomial  random field.
     \item \code{BinomialNeg} for a negative Binomial random field.
     \item \code{BinomialNegZINB} for a zero inflated negative Binomial random field.
  }
random fields using Gaussian and Clayton copula (see Bevilacqua, Alvarado  and Caamaño, 2023) with the following marginal distribution:
 \itemize{   
    \item \code{Gaussian} for Gaussian random field.
    \item \code{Beta2}  for Beta   random field.
    }

  For  a given \code{model} the associated parameters are given by nuisance and correlation parameters. In order to obtain the nuisance parameters associated with a specific model use \code{\link{NuisParam}}. 
  In order to obtain the correlation parameters associated with a given correlation model use 
  \code{\link{CorrParam}}.

   All the nuisance and covariance parameters must be specified
   by the user using  the  \code{start}  and  the \code{fixed} parameter.
   Specifically:

   The option \code{start} sets the starting values  parameters in the optimization process for the parameters  to be estimated.
   The option \code{fixed} parameter allows to fix  some of the parameters.

   Regression parameters in the linear specification must be specified as \code{mean,mean1,..meank} (see \code{\link{NuisParam}}).
   In this case a matrix of covariates with suitable dimension must be specified using   \code{X}. 
   In the case of a single mean then \code{X} should not be specified and it is interpreted as a vector of ones.
   It is also possible to fix a vector of spatial or spatio-temporal means (estimated with other methods for instance).


  Two types  of binary weights can be used in  the weighted composite likelihood estimation based on pairs,
  one based on neighboords and one based on distances.
   
  In the first  case  binary weights are set to 1 or 0 depending if the pairs are neighboords of a certain order (1, 2, 3, ..) specified by the parameter (\code{neighb}). This weighting scheme is effecient for large-data sets since the computation of the 'useful' pairs  in based on  fast nearest neighbour 
  search (Caamaño et al., 2023).

  In the second case, binary weights are set to 1 or 0 depending if the pairs have distance lower than  (\code{maxdist}).
  This weighting scheme is less inefficient for large data. The same arguments of \code{neighb} applies  for \code{maxtime} that sets
   the  order (1, 2, 3, ..) of temporal neighboords in spatial-temporal field.

   The two approaces can be combined by considering  the pairs that have distances lower than  (\code{maxdist})
   in the set of pairs that are neighboords of a certain order.

  For estimation of the variance-covariance matrix of  the  weighted composite likelihood  estimates
   the option \code{sensitivity=TRUE} must be specified. Then the  \code{GeoFit}  object  must be updated   using the function \code{GeoVarestbootstrap}. This allows to estimate  the Godambe Information matrix  (see  Bevilacqua et. al. (2012) , Bevilacqua and Gaetan (2013)).
  Then standard error estimation, confidence intervals,  pvalues  and composite likelihood information critera can be obtained.
  
  
}

\value{
  Returns an object of class \code{GeoFit}.
  An object of class \code{GeoFit} is a list containing
  at most the following components:
  \item{bivariate}{Logical:\code{TRUE} if the Gaussian random fields is bivariate, otherwise  \code{FALSE};}
  \item{clic}{The composite information criterion, if the full
  likelihood is considered then it coincides with the Akaike information
  criterion;}
  \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates;}
  \item{coordx_dyn}{A  list of dynamical (in time) spatial coordinates;}
  \item{conf.int}{Confidence intervals for  standard maximum likelihood estimation;}
  \item{convergence}{A string that denotes if convergence is reached;}
  \item{copula}{The type of copula;}
  \item{corrmodel}{The correlation model;}
  \item{data}{The vector or matrix or array of data;}
  \item{distance}{The type of spatial distance;}
  \item{fixed}{A list of the  fixed parameters;}
  \item{iterations}{The number of iteration used by the numerical routine;}
  \item{likelihood}{The configuration of the composite likelihood;}
  \item{logCompLik}{The value of the log composite-likelihood at the maximum;}
  \item{maxdist}{The maximum spatial distance used in the weigthed composite likelihood. 
                       If no spatial distance is specified then it is NULL;}
  \item{maxtime}{The 
    order of temporal neighborhood  in the composite likelihood
    computation. }
  \item{message}{Extra message passed from the numerical routines;}
  \item{model}{The density associated to the likelihood objects;}
   \item{missp}{True if a misspecified  Gaussian model is ued in  the composite likelihhod;}
  \item{n}{The number of trials in a binominal random fields;the number of successes in a  negative Binomial random fieldss;} 
    \item{neighb}{The 
    order of spatial neighborhood  in the composite likelihood
    computation. }
  \item{ns}{The number of (different) location sites in the bivariate case;} 
   \item{nozero}{In the case of tapered likelihood  the percentage of non zero values in the covariance matrix.
  Otherwise is NULL.}
  \item{numcoord}{The number of spatial coordinates;}
  \item{numtime}{The number of the temporal realisations of the random fields;}
  \item{param}{A list of the   parameters' estimates;}
  \item{radius}{The radius of the sphere in the case of great  circle distance;}
  \item{stderr}{The vector of standard errors for  standard maximum likelihood estimation;}
  \item{sensmat}{The sensitivity matrix;}
  \item{varcov}{The matrix of the variance-covariance of the estimates;}
  \item{type}{The type of the likelihood objects.}
  \item{X}{The matrix of covariates;}
}


\references{


  \bold{General Composite-likelihood}:

  Varin, C., Reid, N. and Firth, D. (2011). An Overview of Composite
  Likelihood Methods. \emph{Statistica Sinica}, \bold{21}, 5--42.

  Varin, C. and Vidoni, P. (2005)
  A Note on Composite Likelihood Inference and Model Selection.
  \emph{Biometrika}, \bold{92}, 519--528.


  \bold{Non Gaussian random fields}: 

  Alegrıa A., Caro S., Bevilacqua M., Porcu E., Clarke J. (2017) \emph{Estimating covariance functions of multivariate skew-Gaussian random fields on the sphere}. Spatial Statistics \bold{22} 388-402

  Alegria A., Bevilacqua, M., Porcu, E. (2016) 
  Likelihood-based inference for multivariate space-time wrapped-Gaussian fields. 
    \emph{Journal of Statistical Computation and Simulation}. \bold{86(13)}, 2583--2597.

  Bevilacqua M., Caamano C.,  Gaetan C. (2020) On modeling positive continuous data with spatio-temporal dependence.  \emph{Environmetrics}   \bold{31(7)}

  Bevilacqua M., Caamaño C., Arellano Valle R.B., Morales-Oñate V.  (2020) Non-Gaussian Geostatistical Modeling using (skew) t Processes. \emph{Scandinavian Journal of Statistics}.

  Blasi F., Caamaño C., Bevilacqua M., Furrer R. (2022) A selective view of climatological data and likelihood estimation  \emph{Spatial Statistics}   10.1016/j.spasta.2022.100596

  Bevilacqua M., Caamaño C., Arellano-Valle R. B., Camilo Gomez C. (2022) A class of random fields with two-piece marginal distributions for modeling point-referenced data with spatial outliers. \emph{Test}  10.1007/s11749-021-00797-5
  
  Morales-Navarrete D., Bevilacqua M., Caamaño C., Castro  L.M. (2022) Modelling Point Referenced Spatial Count Data: A Poisson Process Approach   \emph{TJournal of the American Statistical Association} To appear

  Caamaño C., Bevilacqua M., López C., Morales-Oñate V. (2023) Nearest neighbours weighted composite likelihood based on pairs for (non-)Gaussian massive spatial data with an application to Tukey-hh random fields estimation  \emph{Computational Statistics and Data Analysis} To appear

  Bevilacqua M., Alvarado E., Caamaño C., (2023) A flexible Clayton-like spatial copula with application to bounded support data \emph{Journal of Multivariate Analysis} To appear 

  \bold{Weighted Composite-likelihood for (non-)Gaussian random fields}:

  Bevilacqua, M. Gaetan, C., Mateu, J. and Porcu, E. (2012)
  Estimating space and space-time covariance functions for large data
  sets: a weighted composite likelihood approach.
  \emph{Journal of the American Statistical Association, Theory &
    Methods}, \bold{107}, 268--280.

  Bevilacqua, M., Gaetan, C. (2015) Comparing composite likelihood methods based on pairs for spatial Gaussian random fields.  \emph{Statistics and Computing}, \bold{25(5)},  877-892.

  Caamaño C., Bevilacqua M., López C., Morales-Oñate V. (2023) Nearest neighbours weighted composite likelihood based on pairs for (non-)Gaussian massive spatial data with an application to Tukey-hh random fields estimation  \emph{Computational Statistics and Data Analysis} To appear


  \bold{Sub-sampling estimation}:

  Heagerty, P. J. and Lumley T. (2000) Window Subsampling of Estimating
  Functions with Application to Regression Models.
  \emph{Journal of the American Statistical Association, Theory &
    Methods}, \bold{95}, 197--211.
 
}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}
  
\examples{
library(GeoModels)

###############################################################
############ Examples of spatial Gaussian random fieldss ################
###############################################################


# Define the spatial-coordinates of the points:
set.seed(3)
N=300  # number of location sites
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)

# Define spatial matrix covariates and regression parameters
X=cbind(rep(1,N),runif(N))
mean <- 0.2
mean1 <- -0.5

# Set the covariance model's parameters:
corrmodel <- "Matern"
sill <- 1
nugget <- 0
scale <- 0.2/3
smooth=0.5


param<-list(mean=mean,mean1=mean1,sill=sill,nugget=nugget,scale=scale,smooth=smooth)

# Simulation of the spatial Gaussian random fields:
data <- GeoSim(coordx=coords,corrmodel=corrmodel, param=param,X=X)$data



################################################################
###
### Example 0. Maximum independence composite likelihood fitting of
### a Gaussian random fields (no dependence parameters)
### 
###############################################################
# setting starting parameters to be estimated
start<-list(mean=mean,mean1=mean1,sill=sill)

fit1 <- GeoFit(data=data,coordx=coords,likelihood="Marginal",
                    type="Independence", start=start,X=X)
print(fit1)


################################################################
###
### Example 1. Maximum conditional pairwise likelihood fitting of
### a Gaussian random fields using BFGS
### 
###############################################################
# setting fixed and starting parameters to be estimated
fixed<-list(nugget=nugget,smooth=smooth)
start<-list(mean=mean,mean1=mean1,scale=scale,sill=sill)

fit1 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, 
                    neighb=3,likelihood="Conditional",optimizer="BFGS",
                    type="Pairwise", start=start,fixed=fixed,X=X)
print(fit1)

################################################################
###
### Example 2. Standard Maximum likelihood fitting of
### a Gaussian random fields using nlminb
###
###############################################################
# Define the spatial-coordinates of the points:
set.seed(3)
N=250  # number of location sites
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)

param<-list(mean=mean,sill=sill,nugget=nugget,scale=scale,smooth=smooth)

data <- GeoSim(coordx=coords,corrmodel=corrmodel, param=param)$data

# setting fixed and parameters to be estimated
fixed<-list(nugget=nugget,smooth=smooth)
start<-list(mean=mean,scale=scale,sill=sill)

I=Inf
lower<-list(mean=-I,scale=0,sill=0)
upper<-list(mean=I,scale=I,sill=I)
fit2 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel,
                    optimizer="nlminb",upper=upper,lower=lower,
                    likelihood="Full",type="Standard", 
                    start=start,fixed=fixed)
print(fit2)


###############################################################
############ Examples of spatial non-Gaussian random fieldss #############
###############################################################


################################################################
###
### Example 3. Maximum pairwise likelihood fitting of a Weibull  random fields 
### with Generalized Wendland correlation with Nelder-Mead
### 
###############################################################
set.seed(524)
# Define the spatial-coordinates of the points:
N=300
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)
X=cbind(rep(1,N),runif(N))
mean=1; mean1=2 # regression parameters
nugget=0
shape=2
scale=0.2
smooth=0

model="Weibull"
corrmodel="GenWend"
param=list(mean=mean,mean1=mean1,scale=scale,
                     shape=shape,nugget=nugget,power2=4,smooth=smooth)
# Simulation of a  non stationary weibull random fields:
data <- GeoSim(coordx=coords, corrmodel=corrmodel,model=model,X=X,
           param=param)$data


fixed<-list(nugget=nugget,power2=4,smooth=smooth)
start<-list(mean=mean,mean1=mean1,scale=scale,shape=shape)

# Maximum independence likelihood:
fit <- GeoFit(data=data, coordx=coords, X=X, 
           likelihood="Marginal",type="Independence", corrmodel=corrmodel,
         ,model=model, start=start, fixed=fixed)
print(unlist(fit$param))

## estimating  dependence parameter fixing vector mean   parameter
Xb=as.numeric(X \%*\% c(mean,mean1))
fixed<-list(nugget=nugget,power2=4,smooth=smooth,mean=Xb)
start<-list(scale=scale,shape=shape)

# Maximum conditional composite-likelihood fitting of the random fields:
fit1 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
                    neighb=3,likelihood="Conditional",type="Pairwise",
                    optimizer="Nelder-Mead",
                    start=start,fixed=fixed)
print(unlist(fit1$param))



### joint estimation  of the dependence parameter and  mean   parameters
fixed<-list(nugget=nugget,power2=4,smooth=smooth)
start<-list(mean=mean,mean1=mean1,scale=scale,shape=shape)
fit2 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
                    neighb=3,likelihood="Conditional",type="Pairwise",X=X,
                    optimizer="Nelder-Mead",
                    start=start,fixed=fixed)
print(unlist(fit2$param))



################################################################
###
### Example 4. Maximum pairwise likelihood fitting of
### a Skew-Gaussian spatial  random fields with Wendland correlation
###
###############################################################
set.seed(261)
model="SkewGaussian"
# Define the spatial-coordinates of the points:
x <- runif(500, 0, 1)
y <- runif(500, 0, 1)
coords <- cbind(x,y)

corrmodel="Wend0"
mean=0;nugget=0
sill=1
skew=-4.5
power2=4
c_supp=0.2

# model parameters
param=list(power2=power2,skew=skew,
             mean=mean,sill=sill,scale=c_supp,nugget=nugget)
data <- GeoSim(coordx=coords, corrmodel=corrmodel,model=model, param=param)$data

plot(density(data))
fixed=list(power2=power2,nugget=nugget)
start=list(scale=c_supp,skew=skew,mean=mean,sill=sill)
lower=list(scale=0,skew=-I,mean=-I,sill=0)
upper=list(scale=I,skew=I,mean=I,sill=I)
# Maximum marginal pairwise likelihood:
fit1 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
                    neighb=3,likelihood="Marginal",type="Pairwise",
                    optimizer="bobyqa",lower=lower,upper=upper,
                    start=start,fixed=fixed)
print(unlist(fit1$param))


################################################################
###
### Example 5. Maximum pairwise likelihood fitting of 
### a Binomial random fields with exponential correlation 
###
###############################################################

set.seed(422)
N=250
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)
mean=0.1; mean1=0.8; mean2=-0.5 # regression parameters
X=cbind(rep(1,N),runif(N),runif(N)) # marix covariates
corrmodel <- "Wend0"
param=list(mean=mean,mean1=mean1,mean2=mean2,nugget=0,scale=0.2,power2=4)
# Simulation of the spatial Binomial-Gaussian random fields:
data <- GeoSim(coordx=coords, corrmodel=corrmodel, model="Binomial", n=10,X=X,
             param=param)$data


## estimating the marginal parameters using independence cl
fixed <- list(power2=4,scale=0.2,nugget=0)
start <- list(mean=mean,mean1=mean1,mean2=mean2)

# Maximum independence likelihood:
fit <- GeoFit(data=data, coordx=coords,n=10, X=X, 
           likelihood="Marginal",type="Independence",corrmodel=corrmodel,
         ,model="Binomial", start=start, fixed=fixed)
                  
print(fit)


## estimating  dependence parameter fixing vector mean   parameter
Xb=as.numeric(X \%*\% c(mean,mean1,mean2))
fixed <- list(nugget=0,power2=4,mean=Xb)
start <- list(scale=0.2)

# Maximum conditional pairwise likelihood:
fit1 <- GeoFit(data=data, coordx=coords, corrmodel=corrmodel,n=10, 
          likelihood="Conditional",type="Pairwise",  neighb=3
         ,model="Binomial", start=start, fixed=fixed)
                  
print(fit1)


## estimating jointly marginal   and dependnce parameters
fixed <- list(nugget=0,power2=4)
start <- list(mean=mean,mean1=mean1,mean2=mean2,scale=0.2)

# Maximum conditional pairwise likelihood:
fit2 <- GeoFit(data=data, coordx=coords, corrmodel=corrmodel,n=10, X=X, 
          likelihood="Conditional",type="Pairwise",  neighb=3
         ,model="Binomial", start=start, fixed=fixed)
                  
print(fit2)


###############################################################
######### Examples of Gaussian spatio-temporal random fieldss ###########
###############################################################
set.seed(52)
# Define the temporal sequence:
time <- seq(1, 9, 1)

# Define the spatial-coordinates of the points:
x <- runif(20, 0, 1)
y <- runif(20, 0, 1)
coords=cbind(x,y)

# Set the covariance model's parameters:
scale_s=0.2/3;scale_t=1
smooth_s=0.5;smooth_t=0.5
sill=1
nugget=0
mean=0

param<-list(mean=0,scale_s=scale_s,scale_t=scale_t,
 smooth_t=smooth_t, smooth_s=smooth_s ,sill=sill,nugget=nugget)

# Simulation of the spatial-temporal Gaussian random fields:
data <- GeoSim(coordx=coords,coordt=time,corrmodel="Matern_Matern",
              param=param)$data

################################################################
###
### Example 6. Maximum pairwise likelihood fitting of a
### space time Gaussian random fields with double-exponential correlation
###
###############################################################
# Fixed parameters
fixed<-list(nugget=nugget,smooth_s=smooth_s,smooth_t=smooth_t)
# Starting value for the estimated parameters
start<-list(mean=mean,scale_s=scale_s,scale_t=scale_t,sill=sill)

# Maximum composite-likelihood fitting of the random fields:
fit <- GeoFit(data=data,coordx=coords,coordt=time,
                    corrmodel="Matern_Matern",maxtime=1,neighb=3,
                    likelihood="Marginal",type="Pairwise",
                     start=start,fixed=fixed)
print(fit)



###############################################################
######### Examples of a bivariate Gaussian  random fields ###########
###############################################################

################################################################
### Example 7. Maximum pairwise  likelihood fitting of a
### bivariate Gaussian random fields with separable Bivariate  matern 
### (cross) correlation model 
###############################################################

# Define the spatial-coordinates of the points:
set.seed(89)
x <- runif(300, 0, 1)
y <- runif(300, 0, 1)
coords=cbind(x,y)
# parameters
param=list(mean_1=0,mean_2=0,scale=0.1,smooth=0.5,sill_1=1,sill_2=1,
           nugget_1=0,nugget_2=0,pcol=0.2)

# Simulation of a spatial bivariate Gaussian random fields:
data <- GeoSim(coordx=coords, corrmodel="Bi_Matern_sep", 
              param=param)$data

# selecting fixed and estimated parameters
fixed=list(mean_1=0,mean_2=0,nugget_1=0,nugget_2=0,smooth=0.5)
start=list(sill_1=var(data[1,]),sill_2=var(data[2,]),
           scale=0.1,pcol=cor(data[1,],data[2,]))


# Maximum marginal pairwise likelihood
fitcl<- GeoFit(data=data, coordx=coords, corrmodel="Bi_Matern_sep",
                     likelihood="Marginal",type="Pairwise",
                     optimizer="BFGS" , start=start,fixed=fixed,
                     neighb=c(4,4,4))
print(fitcl)

}
\keyword{Composite}
