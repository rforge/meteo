\name{rfilltimegaps}
\alias{rfilltimegaps}
\title{Disaggregation in the time dimension through the use of splines for each pixel}
\description{The function creates an object of \link[spacetime]{STFDF} class, spatio-temporal data with full space-time grid, from another \link[spacetime]{STFDF} and fills attribute data for missing dates using splines. }
\usage{ rfilltimegaps(stfdf,tunits="day",attrname=1, ...)
}
\arguments{
\item{stfdf}{\link[spacetime]{STFDF} object with time information of minimum length 2, and gap between}
\item{tunits}{increment of the sequence used to generete time infromation for temporal gap}
\item{attrname}{varible from \link[spacetime]{STFDF} to be splined}
\item{...}{ arguments passed to \link[stats]{spline} }
}
\value{ \link[spacetime]{STFDF} object}
\author{Milan Kilibarda \email{kili@grf.bg.ac.rs}}
\references{
Kilibarda, M., T. Hengl, G. B. M. Heuvelink, B. Graeler, E. Pebesma, M. Percec Tadic, and B. Bajat (2014), Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution, J. Geophys. Res. Atmos., 119, 2294-2313, doi:10.1002/2013JD020803.
}
\seealso{
\code{\link{rfillspgaps}}
\code{\link{pred.strk}}
}
\examples{
   data(nlmodis20110704)
   data(nlmodis20110712)
   data(NLpol)
   
   # fill spatial gaps
   r1=rfillspgaps(nlmodis20110704,NLpol)
   r2=rfillspgaps(nlmodis20110712,NLpol)
   
   r1 <- as(r1,"SpatialPixelsDataFrame")
   names(r1)='m1'
   r2 <- as(r2,"SpatialPixelsDataFrame")
   names(r2)='m2'

   r1@data <- cbind(r1@data, r2@data)
   
   df <- reshape(r1@data, varying=list(1:2), v.names="modis", direction="long", 
                 times=as.Date(c('2011-07-04','2011-07-12')), ids=1:dim(r1)[1])
   
   library(spacetime)
   stMODIS<- STFDF(as(r1, "SpatialPixels"), 
                  time= as.Date(c('2011-07-04','2011-07-12')), 
                  data.frame(modis=df[,'modis']))
   
   # stplot(stMODIS, col.regions=bpy.colors())
   stMODIS <- rfilltimegaps(stMODIS)
   # stplot(stMODIS, col.regions=bpy.colors())
}
