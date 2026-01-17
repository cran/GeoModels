\name{rainNLD}
\alias{rainNLD}
\docType{data}
\title{April Precipitation over the Netherlands with Distance-to-Coast Covariate}
\usage{data(rainNLD)}
\description{
A numeric matrix containing gridded April total precipitation (in mm) over the Netherlands,
together with projected coordinates (UTM, km) and a distance-to-coast covariate (km).
The precipitation field is derived from the WorldClim v2.1 monthly climatology for April
(i.e., long-term monthly averages rather than a specific year).
}
\format{
A numeric matrix with 71,401 rows and 4 columns:
\describe{
  \item{Column 1}{Easting coordinate (UTM zone 31N), in kilometers.}
  \item{Column 2}{Northing coordinate (UTM zone 31N), in kilometers.}
  \item{Column 3}{Distance to the coastline, in kilometers (computed from Natural Earth coastline geometry, transformed to UTM).}
  \item{Column 4}{April total precipitation, in millimeters.}
}
}
\details{
Coordinates were projected from longitude/latitude (WGS84) to UTM zone 31N and rescaled to kilometers.
The distance-to-coast covariate was computed as the Euclidean distance (in UTM) from each grid cell
to the nearest coastline.
}
\source{
Fick, S.E., Hijmans, R.J. (2017) WorldClim 2: new 1km spatial resolution climate surfaces for global land areas.
\emph{International Journal of Climatology}, \bold{37}, 4302--4315.
}
\keyword{datasets}
