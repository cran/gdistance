\name{shortestPath}
\alias{shortestPath}
\alias{shortestPath,TransitionLayer,Coords,Coords-method}

\title{ Shortest path }
\description{
  Calculates the shortest path from an origin to a goal. 
}
\usage{
shortestPath(transition, origin, goal, ...)
}
\arguments{
  \item{transition}{\code{TransitionLayer} object}
  \item{origin}{SpatialPoints, vector or matrix with coordinates}
  \item{goal}{SpatialPoints, vector or matrix with coordinates}
  \item{...}{Additional argument: output}
}
\details{
  As an additional argument output, the desired output object can be specified: either \dQuote{TransitionLayer} (default) or \dQuote{SpatialLines}.

  If there is more than one path either (1) transition values in the TransitionLayer get values of 1 / number of paths or (2) the SpatialLines object will contain more than one line.
}
\value{
  Transition object.
}

\author{ Jacob van Etten \email{jacobvanetten@yahoo.com}.}

\seealso{ \code{\link[gdistance]{costDistance}}, \code{\link[gdistance]{accCost}} }
\examples{
  #Create a new raster and set all its values to a random value.
  r <- raster(nrows=18, ncols=36) 
  r <- setValues(r,runif(ncell(r), 0.4, 0.6))

  #Create a Transition object from the raster
  tr <- transition(r,mean,8)
  trC <- geoCorrection(tr, type="c")

  #Create two sets of coordinates
  sP1 <- SpatialPoints(cbind(-105,55))
  sP2 <- SpatialPoints(cbind(105,-55))

  #Calculate the shortest path
  sPath1 <- shortestPath(tr, sP1, sP2)
  sPath2 <- shortestPath(tr, sP1, sP2, output="SpatialLines")
  plot(raster(sPath1))
  lines(sPath2)
}

\keyword{ spatial }
