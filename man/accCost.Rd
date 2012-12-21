\name{accCost}
\alias{accCost}
\alias{accCost,TransitionLayer,RasterLayer-method}
\alias{accCost,TransitionLayer,Coords-method}
\title{ Accumulated Cost Surface }
\description{
Calculates the accumulated cost surface from one or more origins. }
\usage{
accCost(x, fromCoords)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ object of class TransitionLayer }
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
\seealso{ \code{\link[igraph]{shortest.paths}}, \code{\link[gdistance]{geoCorrection}}, \code{\link[gdistance]{costDistance}} }
\examples{
#example equivalent to that in the documentation on r.cost in GRASS
r <- raster(nrows=6, ncols=7, xmn=0, xmx=7, ymn=0, ymx=6, crs="+proj=utm +units=m")
 
r[] <- c(2, 2, 1, 1, 5, 5, 5,
         2, 2, 8, 8, 5, 2, 1,
         7, 1, 1, 8, 2, 2, 2,
         8, 7, 8, 8, 8, 8, 5,
         8, 8, 1, 1, 5, 3, 9,
         8, 1, 1, 2, 5, 3, 9)

T <- transition(r, function(x) 1/mean(x), 8) 
# 1/mean: reciprocal to get permeability
T <- geoCorrection(T)

c1 <- c(5.5,1.5) 
c2 <- c(1.5,5.5)

A <- accCost(T, c1)
plot(A)
text(A)

}
\keyword{ spatial }
