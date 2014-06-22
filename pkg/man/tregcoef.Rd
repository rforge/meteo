\name{tregcoef}
\docType{data}
\encoding{utf8}
\alias{tregcoef}
\title{Multiple linear regression coefficients for global daily air temperatures}
\description{
Multiple linear regression coefficients for mean, min., max. daily temperature on  geometric temperature trend, MODIS LST, elevation, and topographic wetness index. The models is computed from GSOD and ECA&D and GHCN-Daily data.
}
\usage{data(tregcoef)}
\format{
  A list of 8 multiple linear regression coefficients for daily air temperatures.
  \describe{
 \item{\code{tmeanGSODECAD}}{Multiple linear regression coefficients of mean daily temperature on  geometric temperature trend, MODIS LST, elevation, and topographic wetness index, data used: GSOD and ECA&D}
\item{\code{tmeanGSODECAD_noMODIS}}{Multiple linear regression coefficients of mean daily temperature on  geometric temperature trend,  elevation, and topographic wetness index, data used: GSOD and ECA&D}
\item{\code{tminGSODECAD}}{Multiple linear regression coefficients of min. daily temperature on  geometric temperature trend, MODIS LST, elevation, and topographic wetness index, data used: GSOD and ECA&D}
\item{\code{tminGHCND}}{Multiple linear regression coefficients of min. daily temperature on  geometric temperature trend, MODIS LST, elevation, and topographic wetness index, data used: GHCN-Daily}
\item{\code{tminGSODECAD_noMODIS}}{Multiple linear regression coefficients of min. daily temperature on  geometric temperature trend,  elevation, and topographic wetness index, data used: GSOD and ECA&D}
\item{\code{tmaxGSODECAD}}{Multiple linear regression coefficients of max. daily temperature on  geometric temperature trend, MODIS LST, elevation, and topographic wetness index, data used: GSOD and ECA&D}
\item{\code{tmaxGHCND}}{Multiple linear regression coefficients of max. daily temperature on  geometric temperature trend, MODIS LST, elevation, and topographic wetness index, data used: GHCN-Daily}
\item{\code{tmaxGSODECAD_noMODIS}}{Multiple linear regression coefficients of max. daily temperature on  geometric temperature trend,  elevation, and topographic wetness index, data used: GSOD and ECA&D}
   
  }
}
\author{Milan Kilibarda}
\references{
Kilibarda, M., T. Hengl, G. B. M. Heuvelink, B. Graeler, E. Pebesma, M. Percec Tadic, and B. Bajat (2014), Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution, J. Geophys. Res. Atmos., 119, 2294-2313, doi:10.1002/2013JD020803.
}
\examples{
data(tregcoef)
tregcoef[[1]] # model for mean daily temp.
}
\keyword{datasets}
