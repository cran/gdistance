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
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}
Distances can be calculated following these steps. 

1. Read spatial grid data into the R environment, using the raster package. Function: \code{\link[raster:raster]{raster}}.

2. Construct an object of the class TransitionLayer or TransitionStack. Function: \code{\link[gdistance:transition]{transition}}.

3. Correct diagonal connections and projection distortions. Function: \code{\link[gdistance:geoCorrection]{geoCorrection}}.

4. Get coordinates for the starting and end points of routes. 

5. Calculate distances and routes. Functions: \code{\link[gdistance:accCost]{accCost}}, 
\code{\link[gdistance:costDistance]{costDistance}}, 
\code{\link[gdistance:resistanceDistance]{resistanceDistance}},
\code{\link[gdistance:shortestPath]{shortestPath}},
\code{\link[gdistance:passage]{passage}},
\code{\link[gdistance:pathInc]{pathInc}}.

}
\author{
Jacob van Etten \email{jacobvanetten@yahoo.com}
}
\references{
  Ray, N. 2005. PATHMATRIX: a geographical information system tool to compute effective distances among samples. \emph{Molecular Ecology Notes} 5, 177 - 180.
  \url{http://cmpg.unibe.ch/software/pathmatrix/}

  McRae, B.H. 2006. Isolation by resistance. \emph{Evolution} 60(8), 1551 - 1561.
  
  McRae B.H., B.G. Dickson, and T. Keitt. 2008. Using circuit theory to model connectivity in ecology, evolution, and conservation. \emph{Ecology} 89:2712-2724.
  \url{http://www.circuitscape.org}
}
\keyword{package}
\keyword{spatial}
\seealso{
\code{\link[raster:raster-package]{raster}}
}
