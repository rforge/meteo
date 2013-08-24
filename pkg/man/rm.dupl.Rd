\name{rm.dupl}
\alias{rm.dupl}
\title{Find point pairs with equal spatial coordinates from \link[spacetime]{STFDF} object.
}
\description{ 
Find point pairs with equal spatial coordinates from \link[spacetime]{STFDF} object and remove stations with less observations.
}
\usage{rm.dupl(obj, zcol = 1, zero.tol = 0)}
\arguments{
  \item{obj}{ \link[spacetime]{STFDF} object}
  \item{zcol}{ variable column name, or column number, from obj@data }
  \item{zero.tol}{distance values less than or equal to this threshold value are considered as duplicates}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ \link[spacetime]{STFDF} object with removed duplicate stations. Stations with less observation is removed, if number of observation is the same for two stations the first is removed.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Milan Kilibarda \email{kili@grf.bg.ac.rs}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# load observation - data frame
data(dtempc) 
# load stations - data frame
data(stations)

str(dtempc)
str(stations)
# create STFDF
temp <- meteo2STFDF(dtempc,stations) # create STFDF from 2 data frames

nrow(temp@sp) # number of stations before removing dupl.

temp2 <-rm.dupl(temp)
nrow(temp2@sp) # number of stations after

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
