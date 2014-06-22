\name{rfillspgaps}
\alias{rfillspgaps}
\title{Close gaps of a grid or raster Layer data}
\description{The function close gaps of a grid or raster Layer data by using IDW}
\usage{ rfillspgaps(rasterLayer,maskPol=NULL,nmax =25, ...)
}
\arguments{
\item{rasterLayer}{raster Layer or SpatiaGrid(Pixels)DF containing NAs}
\item{maskPol}{SpatialPolygons or SpatialPolygonsDataFrame}
\item{nmax}{see \link[gstat]{idw} }
\item{...}{ arguments passed to \link[gstat]{idw} }
}
\value{ raster Layer or SpatiaGrid(Pixels)DF object with NA replaced using IDW}
\author{Milan Kilibarda \email{kili@grf.bg.ac.rs}}
\references{
Kilibarda, M., T. Hengl, G. B. M. Heuvelink, B. Graeler, E. Pebesma, M. Percec Tadic, and B. Bajat (2014), Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution, J. Geophys. Res. Atmos., 119, 2294-2313, doi:10.1002/2013JD020803.
}
\seealso{
\code{\link{rfilltimegaps}}
\code{\link{pred.strk}}
}
\examples{
   data(nlmodis20110712)
   data(NLpol)
   
   # spplot(nlmodis20110712, col.regions=bpy.colors())
   # fill spatial gaps
   r=rfillspgaps(nlmodis20110712,NLpol)
   # spplot(r, col.regions=bpy.colors())
  
}