\name{plot.GeoVariogramDir}
\alias{plot.GeoVariogramDir}
\encoding{UTF-8}
\title{Plot empirical directional semi-variogram}

\description{
  Plots empirical directional semi-variograms for objects of class \code{"GeoVariogramDir"} as produced by \code{\link{GeoVariogramDir}}.
  All directions are displayed in a single plot, each with a different color and a legend indicating the direction (e.g., "0°", "45°", etc.).
}

\usage{
\method{plot}{GeoVariogramDir}(x, ..., main = "Directional Empirical Semivariograms",
  pch = 20, lwd = 1, col = 1:8, ylab = "Semivariogram", xlab = "Distance")
}

\arguments{
  \item{x}{An object of class \code{"GeoVariogramDir"} as produced by \code{\link{GeoVariogramDir}}.}
  \item{main}{A main title for the plot.}
  \item{pch}{Plotting character (point type) for the points (default: 20).}
  \item{lwd}{Line width for the lines connecting points (default: 1).}
  \item{col}{A vector of colors, one for each direction (default: 1:8).}
  \item{ylab}{Label for the y-axis (default: "Semivariogram").}
  \item{xlab}{Label for the x-axis (default: "Lag").}
  \item{\dots}{Additional graphical parameters passed to \code{\link{plot}}.}
}

\details{
  This function plots all empirical directional semi-variograms in a single graph, using different colors and a legend in the top left corner that indicates the direction (e.g., "0°", "45°", etc.).
  Each direction is represented by points connected by lines.
}

\value{
  Produces a plot.
  No values are returned.
}

\seealso{
  \code{\link{GeoVariogramDir}} for directional variogram computation and examples.
}



\keyword{Variogram}
