\name{dtemp}
\docType{data}
\encoding{latin1}
\alias{dtemp}
\title{Mean daily temperature in degrees Celsius for july 2011}
\description{Sample data set showing values of merged mean daily temperature measurements from the Global Surface Summary of day data (GSOD), for the month July 2011.}
\usage{data(dtemp)}
\format{
The \code{dtemp} contains the following columns:
  \describe{
  \item{\code{STNID}}{factor; station ID (World Meteorological Organization unique station number)}
  \item{\code{TimeSpan.begin}}{POSIXct; begin of the measurement period}
  \item{\code{TimeSpan.end}}{POSIXct; end of the measurement period}
  \item{\code{observedValue}}{numeric; mean daily temperature for the day in degrees Celsius to tenths}
  \item{\code{measurementError}}{numeric; measurement error estimated as 1 / number of measured temperatures per day}
}
}
\note{The data summaries provided here are based on data exchanged under the World Meteorological Organization (WMO) World Weather Watch Program. To prepare a point map, merge with the \code{\link{stations}} table containing stations' coordinates.}
\author{Milan Kilibarda and Tomislav Hengl}
\references{
\itemize{
\item Global Surface Summary of the day data (\url{ftp://ftp.ncdc.noaa.gov/pub/data/gsod/}) 
}
}
\examples{
## load data and convert to space-time class:
data(dtemp)
data(wmo)
}
\keyword{datasets}