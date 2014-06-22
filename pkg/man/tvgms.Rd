\name{tvgms}
\docType{data}
\encoding{utf8}
\alias{tvgms}
\title{Spatio-temporal variogram models for global daily air temperatures}
\description{
Variograms of residuals from multiple linear regression of mean, min., max. daily temperatures on  geometric temperature trend, MODIS LST, elevation, and topographic wetness index. The models is computed from GSOD and ECA&D and GHCN-Daily data. The obtained global models for mean, minimum, and maximum  temperature can be used to produce gridded images of daily temperatures at high spatial and temporal resolution.
}
\usage{data(tvgms)}
\format{
  A list of 8 space-time sum-metric models for daily air temperatures, units: space km, time days.
  \describe{
 \item{\code{tmeanGSODECAD}}{Variogram for residuals from multiple linear regression of mean daily temperature on  geometric temperature trend, MODIS LST, elevation, and topographic wetness index, data used: GSOD and ECA&D}
\item{\code{tmeanGSODECAD_noMODIS}}{Variogram for residuals from multiple linear regression of mean daily temperature on  geometric temperature trend,  elevation, and topographic wetness index, data used: GSOD and ECA&D}
\item{\code{tminGSODECAD}}{Variogram for residuals from multiple linear regression of min. daily temperature on  geometric temperature trend, MODIS LST, elevation, and topographic wetness index, data used: GSOD and ECA&D}
\item{\code{tminGHCND}}{Variogram for residuals from multiple linear regression of min. daily temperature on  geometric temperature trend, MODIS LST, elevation, and topographic wetness index, data used: GHCN-Daily}
\item{\code{tminGSODECAD_noMODIS}}{Variogram for residuals from multiple linear regression of min. daily temperature on  geometric temperature trend,  elevation, and topographic wetness index, data used: GSOD and ECA&D}
\item{\code{tmaxGSODECAD}}{Variogram for residuals from multiple linear regression of max. daily temperature on  geometric temperature trend, MODIS LST, elevation, and topographic wetness index, data used: GSOD and ECA&D}
\item{\code{tmaxGHCND}}{Variogram for residuals from multiple linear regression of max. daily temperature on  geometric temperature trend, MODIS LST, elevation, and topographic wetness index, data used: GHCN-Daily}
\item{\code{tmaxGSODECAD_noMODIS}}{Variogram for residuals from multiple linear regression of max. daily temperature on  geometric temperature trend,  elevation, and topographic wetness index, data used: GSOD and ECA&D}
   
  }
}
\author{Milan Kilibarda}
\references{
Kilibarda, M., T. Hengl, G. B. M. Heuvelink, B. Graeler, E. Pebesma, M. Percec Tadic, and B. Bajat (2014), Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution, J. Geophys. Res. Atmos., 119, 2294-2313, doi:10.1002/2013JD020803.
}
\examples{
data(tvgms)
tvgms[[1]] # model for mean daily temp.
}
\keyword{datasets}
