\name{gdistance-package}
\alias{gdistance-package}
\alias{gdistance}
\docType{package}
\title{
gdistance: geographic distance calculations
}
\description{
Calculate distances and routes on geographic grids.
}
\details{
\tabular{ll}{
Package: \tab gdistance\cr
Type: \tab Package\cr
Version: \tab 1.1-1\cr
Date: \tab 2011-01-04\cr
License: \tab GPL (>=3)\cr
LazyLoad: \tab yes\cr
}
Distances can be calculated following these steps. 

1. Read spatial grid data into the R environment, using the raster package. Function: \code{\link[raster:raster]{raster}}.

2. Construct an object of the class TransitionLayer or TransitionStack. Function: \code{\link[gdistance:transition]{transition}}.

3. Correct diagonal connections and projection distortions. Function: \code{\link[gdistance:geoCorrection]{geoCorrection}}.

4. Get coordinates for the starting and end points of routes. 

5. Calculate distances and routes. Functions: \code{\link[gdistance:accCost]{accCost}}, 
\code{\link[gdistance:costDistance]{costDistance}}, 
\code{\link[gdistance:commuteDistance]{commuteDistance}},
\code{\link[gdistance:rSPDistance]{rSPDistance}},
\code{\link[gdistance:shortestPath]{shortestPath}},
\code{\link[gdistance:passage]{passage}},
\code{\link[gdistance:pathInc]{pathInc}}.

}
\author{
Jacob van Etten \email{jacobvanetten@yahoo.com}
}
\references{
Cite this package as:
van Etten, J. 2017. R Package gdistance: Distances and routes on geographical grids. \emph{Journal of Statistical Software}, 76(13), 1-21. doi:10.18637/jss.v076.i13

For functions using random walks, cite:
van Etten, J., and Hijmans, R.J. 2010. A geospatial modelling approach integrating archaeobotany and genetics to trace the origin and dispersal of domesticated plants. \emph{PLoS ONE}, 5(8), p. e12060. doi:10.1371/journal.pone.0012060
}
\keyword{package}
\keyword{spatial}
\seealso{
\code{\link[raster:raster-package]{raster}}
}
