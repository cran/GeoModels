\name{plot.GeoVariogram}
\alias{plot.GeoVariogram}
\encoding{UTF-8}
\title{Plot empirical spatial, spatio-temporal and bivariate semivariograms}

\description{
 Plot an empirical spatial, spatio-temporal, or bivariate semivariogram
 returned by \code{\link{GeoVariogram}}.
}

\usage{
\method{plot}{GeoVariogram}(x, \dots)
}

\arguments{
 \item{x}{an object of the class \code{"GeoVariogram"} }
 \item{\dots}{other arguments to be passed to the function
 \code{\link{plot}} }
}

\details{
 Spatial objects are shown as a standard lag plot. Spatio-temporal objects
 show the empirical surface and its spatial and temporal margins. Bivariate
 objects show the two marginal semivariograms and the cross-semivariogram;
 missing empirical bins are handled as \code{NA}.
}

\value{
 Produces a plot.
 No values are returned.
}

\seealso{
 \code{\link{GeoVariogram}} for variogram computation and examples.
}

\keyword{Variogram}
