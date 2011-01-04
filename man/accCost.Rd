\name{accCost}
\alias{accCost}
\alias{accCost,TransitionLayer,RasterLayer-method}
\alias{accCost,TransitionLayer,Coords-method}
\title{ Accumulated Cost Surface }
\description{
Calculates the accumulated cost surface from one or more origins. }
\usage{
accCost(transition, fromCoords)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{transition}{ object of class TransitionLayer }
  \item{fromCoords}{ origin point locations (SpatialPoints, matrix or numeric class)}
}
\details{
If more than one coordinate is supplied in fromCoords, the function calculates the minimum least-cost distance from any origin point.

The function uses Dijkstra's algorithm (as implemented in the igraph package).
}
\value{
RasterLayer
}
\references{
E.W. Dijkstra. 1959. A note on two problems in connexion with graphs. \emph{Numerische Mathematik} 1, 269 - 271.
}
\author{ Jacob van Etten \email{jacobvanetten@yahoo.com}}
\seealso{ \code{\link[igraph:shortest.paths]{shortest.paths}}, \code{\link[gdistance]{geoCorrection}}, \code{\link[gdistance]{costDistance}} }
\examples{
r <- raster(ncol=36,nrow=18)
r <- setValues(r,rep(1,times=ncell(r)))
tr <- transition(r, mean, directions=8)
tr <- geoCorrection(tr, type="c", multpl=FALSE)
cmd <- accCost(tr, c(0,0))
}
\keyword{ spatial }
