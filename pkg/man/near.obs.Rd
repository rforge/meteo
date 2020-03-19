\name{near.obs}
\alias{near.obs}
\title{Finds n nearest observations from given locations.}
\description{The function finds n nearest observations from given locations and creates an object of \link[base]{data.frame} class, where first n columns are Euclidean distances to n nearest locations and next n columns are observations at n nearest stations, and rows are given locations. Further more it can calculate averages in circles with different radiuses and can find nearest observation in quadrants. It is based on \link[nabor]{knn} function of package nabor.}
\usage{near.obs(locations, observations, zcol = 1, n.obs = 50,
  avg = FALSE,  direct = FALSE,  increment,  range, rm.dupl = TRUE
)}
\arguments{
  \item{locations}{ data.frame with longitude (or x) and latitude (or y) columns, or \link[sp]{SpatialPoints}, \link[sp]{SpatialPointsDataFrame} or \link[sp]{SpatialPixelsDataFrame} object.}
  \item{observations}{ data.frame with longitude (or x), latitude (or y) and observation columns, or \link[sp]{SpatialPoints} or \link[sp]{SpatialPointsDataFrame} object with observation column.}
  \item{zcol}{ numeric or character; give the observation column position or name in \code{observations}  object.}
  \item{n.obs}{ numeric; number of nearest observations to be found. Note that it cannot be larger than number of obsevrations.}
  \item{avg}{ boolean; Will averages in circles with different radiuses be calculated.}
  \item{direct}{ boolean; Will nearest observation in quadrants be calculated.}
  \item{increment}{ numeric; Increment of radiuses for calculation of averages in circles with different radiuses. Units depends on CRS of coordinates.}
  \item{range}{ numeric; Maximum radius for calculation of averages in circles with different radiuses. Units depends on CRS of coordinates.}
  \item{rm.dupl}{ boolean; Will duplicates, i.e. nearest observations where Euclidean distance is 0, be removed from the result.}
}
\value{ \link[base]{data.frame} object.  Rows are given locations. First n columns are Euclidean distances to n nearest locations. Next n columns are observations at n nearest stations. The following columns are averages in circles with different radiuses if avg is set to True. The following columns are nearest observation in quadrants if direct is set to True.
}
\author{ALeksandar Sekulic \email{asekulic@grf.bg.ac.rs}}
\note{The function can be used in any case if it is needed to find n nearest observations from given locations and distances to them.}
\seealso{
\code{\link[nabor]{knn}}
}
\examples{
# prepare data
# load observation - data.frame of mean temperatures
data(dtempc) 
str(dtempc)
data(stations)
#
str(stations)
lonmin=18 ;lonmax=22.5 ; latmin=40 ;latmax=46
library(sp)
library(spacetime)
serbia = point.in.polygon(stations$lon, stations$lat, c(lonmin,lonmax,lonmax,lonmin), 
                          c(latmin,latmin,latmax,latmax))
stations= stations[ serbia!=0, ] # stations in Serbia approx.
# create STFDF
temp <- meteo2STFDF(dtempc,st, crs= CRS('+proj=longlat +datum=WGS84'))
str(temp)
# get one day (SpatialPointsDataFrame)
temp <- temp[, 1]
# find 5 nearest observations and distances to them (remove duplicates)
nearest_obs <- near.obs(locations = temp,
                        observations = temp,
                        zcol = "tempc",
                        n.obs = 5,
                        rm.dupl = TRUE)
str(nearest_obs)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

