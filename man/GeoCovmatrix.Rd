\name{GeoCovmatrix}
\alias{GeoCovmatrix}
\encoding{UTF-8}
\title{Spatial and Spatio-temporal Covariance Matrix of (non) Gaussian random fields}
\description{
  The function computes the  covariance matrix associated to 
  a spatial or spatio(-temporal) or a bivariate spatial Gaussian or non Gaussian randomm field
  with given  underlying covariance model and a  set of spatial location sites (and temporal instants).
}
\usage{
GeoCovmatrix(estobj=NULL,coordx, coordy=NULL,coordz=NULL, coordt=NULL,coordx_dyn=NULL, 
            corrmodel,distance="Eucl", grid=FALSE, maxdist=NULL, maxtime=NULL,
          model="Gaussian", n=1, param, anisopars=NULL, radius=1, sparse=FALSE,
          taper=NULL, tapsep=NULL, type="Standard",copula=NULL,X=NULL,spobj=NULL)
}

\arguments{
  \item{estobj}{An object of class Geofit that includes information about data, model and estimates.}
       \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix  or  (\eqn{d \times 3}{d x 3})-matrix 
     Coordinates on a sphere for a  fixed radius \code{radius} 
    are passed in lon/lat format expressed in decimal degrees.}
  \item{coordy}{A numeric vector giving 1-dimension of
    spatial coordinates; Optional argument, the default is \code{NULL}.}
      \item{coordz}{A numeric vector giving 1-dimension of
    spatial coordinates;  Optional argument, the default is \code{NULL}.}
  \item{coordt}{A numeric vector giving 1-dimension of
    temporal coordinates. At the moment implemented only for the
    Gaussian case. Optional argument, the default is \code{NULL}
    then a spatial random field is expected.}
   \item{coordx_dyn}{A list of \eqn{T} numeric (\eqn{d_t \times 2}{d x 2})-matrices
       containing dynamical (in time) coordinates. Optional argument, the default is \code{NULL}
    }
  \item{corrmodel}{String; the name of a correlation model, for the
    description see the Section \bold{Details}.}
  \item{distance}{String; the name of the spatial distance. The default
  is \code{Eucl}, the euclidean distance. See \code{\link{GeoFit}}.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as spatial or spatial-temporal realisations on a set
    of non-equispaced spatial sites (irregular grid). See
    \code{\link{GeoFit}}.}
   \item{maxdist}{Numeric; an optional positive value indicating the
  marginal spatial compact support in the case of tapered covariance matrix.
   See \code{\link{GeoFit}}.}
  \item{maxtime}{Numeric; an optional positive value indicating the
  marginal temporal compact support in the case of  spacetime tapered 
  covariance matrix.  See \code{\link{GeoFit}}.}
 \item{n}{Numeric; the number of trials in a binomial random fields.
         Default is \eqn{1}.} 
 \item{model}{String; the type of RF. See \code{\link{GeoFit}}.}
  \item{param}{A list of parameter values required  for the covariance
  model.}
      \item{anisopars}{ A list of two elements "angle" and "ratio" i.e. the anisotropy angle and the anisotropy ratio, respectively.}
  \item{radius}{Numeric; a value indicating  the radius of the sphere 
     when using covariance models valid using the great circle distance. 
     Default value is 1} 
  \item{sparse}{Logical; if \code{TRUE} the function return an object of class spam.
  This option  should be used when  a parametric compactly supporte covariance is 
  used. Default is FALSE. }
  \item{taper}{String; the name of the taper correlation
    function if type is  \code{Tapering}, see the Section \bold{Details}.}
   \item{tapsep}{Numeric; an optional  value indicating the separabe parameter in
    the space-time non separable taper or the colocated correlation parameter in a bivariate 
     spatial taper (see \bold{Details}).} 
   \item{type}{String; the type of covariance matrix
   \code{Standard} (the default) or \code{Tapering} for tapered covariance matrix. }
     \item{copula}{String; the type of copula. It can be "Clayton" or "Gaussian"}
   \item{X}{Numeric; Matrix of space-time covariates.}
         \item{spobj}{An object of class sp or spacetime}
   }


\value{
  Returns an object of class \code{GeoCovmatrix}.
  An object of class \code{GeoCovmatrix} is a list containing
  at most the following components:
  \item{bivariate}{Logical:\code{TRUE} if the Gaussian random field is bivariaete
       otherwise  \code{FALSE};}
  \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates;}
  \item{coordx_dyn}{A  list of \eqn{t}{t}  matrices of spatial coordinates;}
  \item{covmatrix}{The covariance matrix if \code{type} is\code{Standard}. An object of 
        class spam if \code{type} is \code{Tapering} or \code{Standard} and
        \code{sparse} is TRUE.}
  \item{corrmodel}{String: the correlation model;}
  \item{distance}{String: the type of spatial distance;}
  \item{grid}{Logical:\code{TRUE} if the spatial data are in a regular grid, otherwise  \code{FALSE};}
  \item{nozero}{In the case of tapered matrix  the percentage of non zero values in the 
        covariance matrix. Otherwise is NULL.}
  \item{maxdist}{Numeric: the marginal spatial compact support 
        if \code{type} is \code{Tapering};}
  \item{maxtime}{Numeric: the marginal temporal compact support 
        if \code{type} is \code{Tapering}; }
  \item{n}{The number of trial for Binomial RFs}       
  \item{namescorr}{String: The names of the correlation parameters; }
  \item{numcoord}{Numeric: the number of spatial coordinates;}
  \item{numtime}{Numeric: the number the temporal  coordinates;}
  \item{model}{The type of RF, see \code{\link{GeoFit}}.}
  \item{param}{Numeric: The covariance parameters;}
  \item{tapmod}{String: the taper model if \code{type} is \code{Tapering}. Otherwise is NULL.}
  \item{spacetime}{\code{TRUE} if spatio-temporal and \code{FALSE} if spatial covariance model;}
  \item{sparse}{Logical: is the  returned object of class spam? ;}

 
}

\details{
In the spatial case, the covariance matrix of the 
  random vector \deqn{[Z(s_1),\ldots,Z(s_n,)]^T}  
  with a specific spatial covariance model is computed. Here  \eqn{n} is the number of the spatial location sites.

 In the space-time case,  the covariance matrix of the 
  random vector \deqn{[Z(s_1,t_1),Z(s_2,t_1),\ldots,Z(s_n,t_1),\ldots,Z(s_n,t_m)]^T}
    with a specific  space time covariance model is computed.   Here \eqn{m}  is the number of temporal instants.

  In the bivariate case, the covariance matrix of the 
  random vector \deqn{[Z_1(s_1),Z_2(s_1),\ldots,Z_1(s_n),Z_2(s_n)]^T}
  with a specific spatial bivariate covariance model is computed.

  The location site \eqn{s_i} can be a point in the \eqn{d}-dimensional euclidean space with \eqn{d=2}
  or a point (given in lon/lat degree format) on a  sphere of arbitrary radius.


A list with all the implemented space and space-time and bivariate 
 correlation models is given below.
The argument \code{param} is a list including all the parameters of a given
correlation  model specified by  the argument \code{corrmodel}.
For each correlation model one can check the associated  parameters' names using \code{\link{CorrParam}}.
In what follows 
\eqn{\kappa>0}, \eqn{\beta>0}, \eqn{ \alpha, \alpha_s, \alpha_t \in (0,2] },  and \eqn{\gamma \in [0,1]}.
The associated parameters in the argument \code{param} are 
\code{smooth}, \code{power2}, \code{power}, \code{power_s}, \code{power_t} and \code{sep} respectively.
 Moreover let \eqn{1(A)=1}  when \eqn{A} is true and \eqn{0} otherwise.


\itemize{
  \item Spatial correlation models:
    \enumerate{

    \item \eqn{GenCauchy} (generalised \eqn{Cauchy} in Gneiting and Schlater (2004)) defined as:
    \deqn{R(h) = ( 1+h^{\alpha} )^{-\beta / \alpha}} 
      If \eqn{h} is the geodesic distance then \eqn{\alpha \in (0,1]}.

  \item \eqn{Matern} defined as:
    \deqn{R(h) = 2^{1-\kappa} \Gamma(\kappa)^{-1} h^\kappa K_\kappa(h)}{2^{1-\kappa} \Gamma(\kappa)^{-1} h^\kappa K_\kappa(h)}
     If \eqn{h} is the geodesic distance then \eqn{\kappa \in (0,0.5]}

   \item \eqn{Kummer} (Kummer hypergeometric in  Ma and Bhadra (2022)) defined as:
    \deqn{R(h) = \Gamma(\kappa+\alpha) U(\alpha,1-\kappa,0.5 h^2 ) / \Gamma(\kappa+\alpha)}
     \eqn{U(.,.,.)} is the Kummer hypergeometric function.
     If \eqn{h} is the geodesic distance then \eqn{\kappa \in (0,0.5]}

    \item \eqn{Kummer{\_}Matern} It is a rescaled version of the  \eqn{Kummer} model  
     that is \eqn{h} must be divided by \eqn{(2*(1+\alpha))^{0.5}}.  When  \eqn{\alpha} goes to infinity it is the Matern model. 

    \item \eqn{Wave} defined as:
    \deqn{R(h)=sin(h)/h} 
    This model is valid only for dimensions less than or equal to 3.

    \item \eqn{GenWend} (Generalized Wendland in  Bevilacqua et al.(2019)) defined as:
      \deqn{R(h) =  A (1-h^2)^{\beta+\kappa} F(\beta/2,(\beta+1)/2,2\beta+\kappa+1,1-h^2) 1(h \in [0,1]) }
      where   \eqn{\mu \ge 0.5(d+1)+\kappa}
       and \eqn{A=(\Gamma(\kappa)\Gamma(2\kappa+\beta+1))/(\Gamma(2\kappa)\Gamma(\beta+1-\kappa)2^{\beta+1})} and $F(.,.,.)$ 
      is the Gaussian hypergeometric function.
     The cases \eqn{\kappa=0,1,2} correspond to the \eqn{Wend0}, \eqn{Wend1} and \eqn{Wend2}  models respectively. 

    \item \eqn{GenWend_{\_}Matern} (Generalized Wendland Matern in  Bevilacqua et al. (2022)). It is defined as a rescaled version of the Generalized Wendland that is \eqn{h} must be divided by 
    \eqn{(\Gamma(\beta+2\kappa+1)/\Gamma(\beta))^{1/(1+2\kappa)}}. When  \eqn{\beta} goes to infinity it is the Matern model.

    \item \eqn{GenWend_{\_}Matern2} (Generalized Wendland Matern second parametrization. It is defined as a rescaled version of the Generalized Wendland that is \eqn{h} must be multiplied  by \eqn{\beta}  and the smoothness parameter is
    \eqn{\kappa-0.5}. When  \eqn{\beta} goes to infinity it is the Matern model.
    

        \item \eqn{Hypergeometric} (Hypegeometric model in Bevilacqua et al  (2025)). 


     \item \eqn{Hypergeometric_{\_}Matern} (Hypegeometric model first parametrization). 
     
     \item \eqn{Hypergeometric_{\_}Matern2} (Hypegeometric model second  parametrization). 

      \item \eqn{Multiquadric} defined as:
    \deqn{R(h) = (1-\alpha0.5)^{2\beta}/(1+(\alpha0.5)^{2}-\alpha cos(h))^{\beta}, \quad h \in [0,\pi]}
    This model is valid on the unit sphere and \eqn{h} is the geodesic distance.
     \item \eqn{Sinpower} defined as:
    \deqn{R(h) = 1-(sin(h/2))^{\alpha},\quad h \in [0,\pi]}
      This model is valid on the unit sphere and \eqn{h} is the geodesic distance.
  
    \item \eqn{F{\_}Sphere} (F family in   Alegria et al. (2021)) defined as:
    \deqn{R(h) =   K*F(1/{\alpha},1/{\alpha}+0.5,
                2/{\alpha}+0.5+{\kappa}),\quad h \in [0,\pi]}
      where      \eqn{K =(\Gamma(a)\Gamma(i))/\Gamma(i)\Gamma(o)) }.      
      This model is valid on the unit sphere and \eqn{h} is the geodesic distance.
      
  }

  \item Spatio-temporal correlation models. 

  \itemize{
    \item Non-separable models:
    \enumerate{
      \item \eqn{Gneiting} defined as:
    \deqn{R(h, u) = e^{ -h^{\alpha_s}/((1+u^{\alpha_t})^{0.5 \gamma \alpha_s })}/(1+u^{\alpha_t})}

  \item \eqn{Gneiting}\code{_}\eqn{GC} 
    \deqn{R(h, u) =  e^{  -u^{\alpha_t} /((1+h^{\alpha_s})^{0.5 \gamma \alpha_t}) }/( 1+h^{\alpha_s}) }
    where \eqn{h}  can be both the euclidean and the geodesic distance

    \item \eqn{Iacocesare} 
    \deqn{R(h, u) = (1+h^{\alpha_s}+u^\alpha_t)^{-\beta}}

       \item \eqn{Porcu} 
    \deqn{R(h, u) = (0.5 (1+h^{\alpha_s})^\gamma +0.5 (1+u^{\alpha_t})^\gamma)^{-\gamma^{-1}}}

       \item \eqn{Porcu1} 
    \deqn{R(h, u) =(e^{ -h^{\alpha_s} ( 1+u^{\alpha_t})^{0.5 \gamma \alpha_s}}) / ((1+u^{\alpha_t})^{1.5})}
    
     \item \eqn{Stein} 
  \deqn{R(h, u) = (h^{\psi(u)}K_{\psi(u)}(h))/(2^{\psi(u)}\Gamma(\psi(u)+1))}
  where \eqn{\psi(u)=\nu+u^{0.5\alpha_t} }




      \item \eqn{Gneiting}\code{_}\eqn{mat}\code{_}\eqn{S},  defined as:
         \deqn{R(h, u) =\phi(u)^{\tau_t}Mat(h\phi(u)^{-\beta},\nu_s)}
      where  \eqn{\phi(u)=(1+u^{0.5\alpha_t})},
       \eqn{\tau_t \ge  3.5+\nu_s}, \eqn{\beta \in [0,1]}
      
      \item \eqn{Gneiting}\code{_}\eqn{mat}\code{_}\eqn{T},  defined interchanging h with u 
      in \eqn{Gneiting}\code{_}\eqn{mat}\code{_}\eqn{S}


    \item \eqn{Gneiting}\code{_}\eqn{wen}\code{_}\eqn{S},  defined as:
      \deqn{R(h, u) =\phi(u)^{\tau_t}GenWend(h \phi(u)^{\beta},\nu_s,\mu_s)}
      where  \eqn{\phi(u)=(1+u^{0.5\alpha_t})},
      \eqn{\tau_t \geq 2.5+2\nu_s}, \eqn{\beta \in [0,1]}
      
      \item \eqn{Gneiting}\code{_}\eqn{wen}\code{_}\eqn{T},  defined interchanging h with u 
      in \eqn{Gneiting}\code{_}\eqn{wen}\code{_}\eqn{S}


        \item \eqn{Multiquadric}\code{_}\eqn{st} defined as:
  \deqn{R(h, u)= ((1-0.5\alpha_s)^2/(1+(0.5\alpha_s)^2-\alpha_s \psi(u) cos(h)))^{a_s} , \quad h \in [0,\pi]}
  where \eqn{\psi(u)=(1+(u/a_t)^{\alpha_t})^{-1}}. This model is valid on the unit
   sphere and \eqn{h} is the geodesic distance.


   \item \eqn{Sinpower}\code{_}\eqn{st} defined as:
  \deqn{R(h, u)=(e^{\alpha_s cos(h) \psi(u)/a_s} (1+\alpha_s cos(h) \psi(u) /a_s))/k}
   where \eqn{\psi(u)=(1+(u/a_t)^{\alpha_t})^{-1} }    and \eqn{k=(1+\alpha_s/a_s) exp(\alpha_s/a_s), \quad h \in [0,\pi]}
   This model is valid on the unit sphere and \eqn{h} is the geodesic distance.
  

      }  


  \item Separable models.\cr

      Space-time separable correlation models are easly obtained as the product of a spatial and a temporal correlation model, that is
 \deqn{R(h,u)=R(h) R(u)}
 Several combinations are possible:
    \enumerate{
    \item  \eqn{Exp}\code{_}\eqn{Exp}  defined as:
    \deqn{R(h, u) =Exp(h)Exp(u)}
    \item  \eqn{Matern}\code{_}\eqn{Matern} defined as:
      \deqn{R(h, u) =Matern(h;\kappa_s)Matern(u;\kappa_t)}
      \item  \eqn{GenWend}\code{_}\eqn{GenWend} defined as 
     \deqn{R(h, u) = GenWend(h;\kappa_s,\mu_s) GenWend(u;\kappa_t,\mu_t)}.
    \item  \eqn{Stable}\code{_}\eqn{Stable} defined as:
      \deqn{R(h, u) =Stable(h;\alpha_s)Stable(u;\alpha_t)}

  }

  Note that some models are nested. (The \eqn{Exp}\code{_}\eqn{Exp} with \eqn{Matern}\code{_}\eqn{Matern} for instance.)
  }
    \item Spatial  bivariate correlation models (see below):
  \enumerate{
    \item  \eqn{Bi}\code{_}\eqn{Matern} (Bivariate  full Matern model) 
  \item  \eqn{Bi}\code{_}\eqn{Matern}\code{_}\eqn{contr} (Bivariate  Matern model with contrainsts) 
   \item  \eqn{Bi}\code{_}\eqn{Matern}\code{_}\eqn{sep} (Bivariate  separable Matern model ) 
       \item  \eqn{Bi}\code{_}\eqn{LMC} (Bivariate  linear model of coregionalization)
    \item  \eqn{Bi}\code{_}\eqn{LMC}\code{_}\eqn{contr} (Bivariate  linear model of coregionalization with constraints ) 
     \item  \eqn{Bi}\code{_}\eqn{Wendx} (Bivariate  full Wendland model) 
       \item  \eqn{Bi}\code{_}\eqn{Wendx}\code{_}\eqn{contr} (Bivariate  Wendland model with contrainsts) 
          \item  \eqn{Bi}\code{_}\eqn{Wendx}\code{_}\eqn{sep} (Bivariate  separable Wendland model) 
          \item  \eqn{Bi}\code{_}\eqn{F{\_}Sphere} (Bivariate  full F_model model on the unit sphere) 
    }


  \item Spatial taper.\cr
  For spatial covariance tapering the taper functions are:
    \enumerate{

    \item \eqn{Bohman} defined as:
      \deqn{T(h)=(1-h)(sin(2\pi h)/(2 \pi h))+(1-cos(2\pi h))/(2\pi^{2}h) 1_{[0,1]}(h)}
    \item \eqn{Wendlandx, \quad x=0,1,2} defined as:
         \deqn{T(h)=Wendx(h;x+2),   x=0,1,2}.
      }

  \item Spatio-temporal tapers. \cr
     For spacetime covariance tapering the taper functions are:
    \enumerate{
     \item \eqn{Wendlandx}\code{_}\eqn{Wendlandy}  (Separable tapers) \eqn{x,y=0,1,2} defined as:
        \deqn{T(h,u)=Wendx(h;x+2) Wendy(h;y+2),    x,y=0,1,2.}

     \item \eqn{Wendlandx}\code{_}\eqn{time} (Non separable temporal taper)  \eqn{x=0,1,2} defined as:
         \eqn{Wenx}\code{_}\eqn{time}, \eqn{x=0,1,2}  assuming \eqn{\alpha_t=2}, \eqn{\mu_s=3.5+x}
         and \eqn{\gamma \in [0,1]} to be fixed using \code{tapsep}.

    \item \eqn{Wendlandx}\code{_}\eqn{space} (Non separable spatial taper)  \eqn{x=0,1,2} defined as:
         \eqn{Wenx}\code{_}\eqn{space}, \eqn{x=0,1,2}  assuming \eqn{\alpha_s=2}, \eqn{\mu_t=3.5+x}
         and \eqn{\gamma \in [0,1]} to be fixed using \code{tapsep}.
        }

      


  \item Spatial  bivariate taper (see below).
      \enumerate{
       \item \eqn{Bi}\code{_}\eqn{Wendlandx, \quad x=0,1,2}
    }

     
}

 \bold{Remarks}:\cr
 In what follows we assume \eqn{\sigma^2,\sigma_1^2,\sigma_2^2,\tau^2,\tau_1^2,\tau_2^2,
 a,a_s,a_t,a_{11},a_{22},a_{12},\kappa_{11},\kappa_{22},\kappa_{12},f_{11},f_{12},f_{21},f_{22}} positive.

 The associated names of the parameters in  \code{param} are 
\code{sill}, \code{sill_1},\code{sill_2},
\code{nugget}, \code{nugget_1},\code{nugget_2},
\code{scale},\code{scale_s},\code{scale_t}, \code{scale_1},\code{scale_2},\code{scale_12},
\code{smooth_1},\code{smooth_2},\code{smooth_12}, \code{a_1},\code{a_12},\code{a_21},\code{a_2}
 respectively.



   Let \eqn{R(h)} be a  spatial correlation model given in standard notation.
   Then the covariance model
   applied with arbitrary variance, nugget  and scale equals to \eqn{\sigma^2} if \eqn{h=0} and
  \deqn{C(h)=\sigma^2(1-\tau^2 ) R( h/a,,...), \quad     h > 0}
  with nugget parameter  \eqn{\tau^2} between 0 and 1.
  Similarly if \eqn{R(h,u)} is a spatio-temporal correlation model given in standard notation,
   then the covariance model is \eqn{\sigma^2} if \eqn{h=0} and \eqn{u=0} and 
  \deqn{C(h,u)=\sigma^2(1-\tau^2 )R(h/a_s ,u/a_t,...) \quad h>0, u>0 }

  Here \sQuote{...} stands for additional parameters.

  The bivariate models implemented are the following :
  \enumerate{
  \item \eqn{Bi}\code{_}\eqn{Matern}  defined as: 
  \deqn{C_{ij}(h)=\rho_{ij} (\sigma_i  \sigma_j+\tau_i^2 1(i=j,h=0))  Matern(h/a_{ij},\kappa_{ij}) \quad i,j=1,2.\quad h\ge 0}
 where \eqn{\rho=\rho_{12}=\rho_{21}} is the correlation colocated parameter and   \eqn{\rho_{ii}=1}.
 The model  \eqn{Bi}\code{_}\eqn{Matern}\code{_}\eqn{sep}  (separable matern) is a special case when
  \eqn{a=a_{11}=a_{12}=a_{22}} and \eqn{\kappa=\kappa_{11}=\kappa_{12}=\kappa_{22}}.
   The model  \eqn{Bi}\code{_}\eqn{Matern}\code{_}\eqn{contr}  (constrained matern)
   is a special case when \eqn{a_{12}=0.5  (a_{11} + a_{22})} and  \eqn{\kappa_{12}=0.5  (\kappa_{11} + \kappa_{22})}

 \item \eqn{Bi}\code{_}\eqn{GenWend}   defined as:
    \deqn{C_{ij}(h)=\rho_{ij} (\sigma_i  \sigma_j+\tau_i^2 1(i=j,h=0))   GenWend(h/a_{ij},\nu_{ij},\kappa_ij) \quad i,j=1,2.\quad h\ge 0}
     where \eqn{\rho=\rho_{12}=\rho_{21}} is the correlation colocated parameter and   \eqn{\rho_{ii}=1}.
   The model  \eqn{Bi}\code{_}\eqn{GenWend}\code{_}\eqn{sep}  (separable Genwendland) is a special case when
  \eqn{a=a_{11}=a_{12}=a_{22}} and \eqn{\mu=\mu_{11}=\mu_{12}=\mu_{22}}.
   The model  \eqn{Bi}\code{_}\eqn{GenWend}\code{_}\eqn{contr}  (constrained Genwendland)
   is a special case when \eqn{a_{12}=0.5  (a_{11} + a_{22})} and  \eqn{\mu_{12}=0.5  (\mu_{11} + \mu_{22})}

  \item  \eqn{Bi}\code{_}\eqn{LMC} defined as:
   \deqn{C_{ij}(h)=\sum_{k=1}^{2} (f_{ik}f_{jk}+\tau_i^2 1(i=j,h=0))   R(h/a_{k})       }
   where  \eqn{R(h)} is a correlation model.
    The model \eqn{Bi}\code{_}\eqn{LMC}\code{_}\eqn{contr}
    is a special case when  \eqn{f=f_{12}=f_{21}}.   
 Bivariate LMC models, in the current version of the package, is obtained with \eqn{R(h)}  equal to the exponential correlation  model.
 }

}
\references{

  Alegria, A.,Cuevas-Pacheco, F.,Diggle, P, Porcu E. (2021)
   The F-family of covariance functions: A Matérn analogue for modeling random fields on spheres.
   \emph{Spatial Statistics} 43, 100512

  Bevilacqua, M., Faouzi, T., Furrer, R., and Porcu, E. (2019). Estimation and prediction 
  using generalized Wendland functions under fixed domain asymptotics. 
  \emph{Annals of Statistics}, 47(2), 828--856.

  Bevilacqua, M., Caamano-Carrillo, C., and Porcu, E. (2022). Unifying compactly supported 
  and Matérn covariance functions in spatial statistics. 
  \emph{Journal of Multivariate Analysis}, 189, 104949.

  Daley J. D., Porcu E., Bevilacqua M. (2015) 
  Classes of compactly supported covariance functions for multivariate random fields. 
  \emph{Stochastic Environmental Research and Risk Assessment}. 29 (4), 1249--1263. 

  Emery, X. and Alegria, A. (2022). The gauss hypergeometric covariance kernel for modeling second-order stationary random fields in euclidean spaces: its compact support, properties and spectral representation. \emph{Stochastic Environmental Research and Risk Assessment}. 36 2819--2834.

  Gneiting, T. (2002). Nonseparable, stationary covariance functions for space-time data.
  \emph{Journal of the American Statistical Association}, 97, 590--600.

  Gneiting T, Kleiber W., Schlather M. 2010. Matern cross-covariance 
  functions for multivariate random fields.  
  \emph{Journal of the American Statistical Association}, 105, 1167--1177.

  Ma, P., Bhadra, A. (2022). 
  Beyond Matérn: on a class of interpretable confluent hypergeometric covariance functions.
  \emph{Journal of the American Statistical Association},1--14.

  Porcu, E.,Bevilacqua, M. and Genton M.  (2015) 
  Spatio-Temporal Covariance and Cross-Covariance Functions of the Great Circle Distance on a Sphere. 
   \emph{Journal of the American Statistical Association}. 
   DOI: 10.1080/01621459.2015.1072541


   Gneiting, T. and Schlater M. (2004) 
   Stochastic models that separate fractal dimension and the Hurst effect. 
    \emph{SSIAM Rev} 46, 269--282.

}

\seealso{
  \code{\link{GeoKrig}},
  \code{\link{GeoSim}},
   \code{\link{GeoFit}}}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{
library(GeoModels)
################################################################
###
### Example 1.Estimated  spatial covariance matrix associated to
### a WCL estimates using the Matern correlation model
###
###############################################################

set.seed(3)
N=300  # number of location sites
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)
# Set the covariance model's parameters:
corrmodel <- "Matern"
mean=0.5; 
sill <- 1
nugget <- 0
scale <- 0.2/3
smooth=0.5

param<-list(mean=mean,sill=sill,nugget=nugget,scale=scale,smooth=smooth)
data <- GeoSim(coordx=coords,corrmodel=corrmodel, param=param)$data
fixed<-list(nugget=nugget,smooth=smooth)
start<-list(mean=mean,scale=scale,sill=sill)

fit0 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, 
                    neighb=3,likelihood="Conditional",optimizer="BFGS",
                    type="Pairwise", start=start,fixed=fixed)
print(fit0)

#estimated covariance matrix using Geofit object
mm=GeoCovmatrix(fit0)$covmatrix
#estimated covariance matrix
mm1 = GeoCovmatrix(coordx=coords,corrmodel=corrmodel,
           param=c(fit0$param,fit0$fixed))$covmatrix
sum(mm-mm1)

################################################################
###
### Example 2. Spatial covariance matrix associated to
### the GeneralizedWendland-Matern correlation model
###
###############################################################

# Correlation Parameters for Gen Wendland model 
CorrParam("GenWend_Matern")
# Gen Wendland Parameters 
param=list(sill=1,scale=0.04,nugget=0,smooth=0,power2=1/1.5)

matrix2 = GeoCovmatrix(coordx=coords, corrmodel="GenWend_Matern", param=param,sparse=TRUE)

# Percentage of no zero values 
matrix2$nozero


################################################################
###
### Example 3. Spatial covariance matrix associated to
### the Kummer correlation model
###
###############################################################

# Correlation Parameters for kummer model 
CorrParam("Kummer")
param=list(sill=1,scale=0.2,nugget=0,smooth=0.5,power2=1)

matrix3 = GeoCovmatrix(coordx=coords, corrmodel="Kummer", param=param)

matrix3$covmatrix[1:4,1:4]


################################################################
###
### Example 4. Covariance matrix associated to
### the space-time double Matern correlation model
###
###############################################################

# Define the temporal-coordinates:
times = seq(1, 4, 1)

# Correlation Parameters for double Matern model
CorrParam("Matern_Matern")

# Define covariance parameters
param=list(scale_s=0.3,scale_t=0.5,sill=1,smooth_s=0.5,smooth_t=0.5)

# Simulation of a spatial Gaussian random field:
matrix4 = GeoCovmatrix(coordx=coords, coordt=times, corrmodel="Matern_Matern",
                     param=param)

dim(matrix4$covmatrix)

################################################################
###
### Example 5. Spatial Covariance matrix associated to
### a  skew gaussian RF with Matern correlation model
###
###############################################################

param=list(sill=1,scale=0.3/3,nugget=0,skew=4,smooth=0.5)
# Simulation of a spatial Gaussian random field:
matrix5 = GeoCovmatrix(coordx=coords,  corrmodel="Matern", param=param, 
                     model="SkewGaussian")

# covariance matrix
matrix5$covmatrix[1:4,1:4]

################################################################
###
### Example 6. Spatial Covariance matrix associated to
### a  Weibull RF with Genwend correlation model
###
###############################################################

param=list(scale=0.3,nugget=0,shape=4,mean=0,smooth=1,power2=5)
# Simulation of a spatial Gaussian random field:
matrix6 = GeoCovmatrix(coordx=coords,  corrmodel="GenWend", param=param, 
                     sparse=TRUE,model="Weibull")

# Percentage of no zero values 
matrix6$nozero

################################################################
###
### Example 7. Spatial Covariance matrix associated to
### a  binomial gaussian RF with Generalized Wendland correlation model
###
###############################################################

param=list(mean=0.2,scale=0.2,nugget=0,power2=4,smooth=0)
# Simulation of a spatial Gaussian random field:
matrix7 = GeoCovmatrix(coordx=coords,  corrmodel="GenWend", param=param,n=5, 
                     sparse=TRUE,
                     model="Binomial")

as.matrix(matrix7$covmatrix)[1:4,1:4]


################################################################
###
### Example 8.  Covariance matrix associated to
### a bivariate Matern exponential correlation model
###
###############################################################

set.seed(8)
# Define the spatial-coordinates of the points:
x = runif(4, 0, 1)
y = runif(4, 0, 1)
coords = cbind(x,y)

# Parameters 
param=list(mean_1=0,mean_2=0,sill_1=1,sill_2=2,
           scale_1=0.1,  scale_2=0.1,  scale_12=0.1,
           smooth_1=0.5, smooth_2=0.5, smooth_12=0.5,
            nugget_1=0,nugget_2=0,pcol=-0.25)

# Covariance matrix 
matrix8 = GeoCovmatrix(coordx=coords, corrmodel="Bi_matern", param=param)$covmatrix

matrix8

}
\keyword{Simulation}
