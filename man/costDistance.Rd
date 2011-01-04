\name{costDistance}
\alias{costDistance}
\alias{costDistance,TransitionLayer,Coords,missing-method}
\alias{costDistance,TransitionLayer,Coords,Coords-method}
\title{ Cost distance (least-cost distance) }
\description{
Calculate the least-cost distance between points. }
\usage{
costDistance(transition, fromCoords, toCoords)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{transition}{ object of class \code{TransitionLayer} }
  \item{fromCoords}{ first set of point locations (of \code{SpatialPoints}, matrix or numeric class)}  
  \item{toCoords}{ second, optional set of point locations (of \code{SpatialPoints}, matrix or numeric class)}
}
\details{
Cost units between cells are defined as the reciprocal of the values in the transition matrix.

The function uses Dijkstra's algorithm, as implemented in the igraph package.

A projection correction is needed for accuracy in the case of grid data for a longlat raster (see function \code{\link{geoCorrection}}). 
}
\value{
distance matrix (S3 class dist or matrix)
}
\references{ 
E.W. Dijkstra. 1959. A note on two problems in connexion with graphs. \emph{Numerische Mathematik} 1, 269-271.

}
\author{ Jacob van Etten \email{jacobvanetten@yahoo.com} }
\seealso{ \code{\link{accCost}}, \code{\link{geoCorrection}} }

\examples{
#create a new raster and set all its values to unity.
r <- raster(nrows=18, ncols=36) 
r <- setValues(r,runif(ncell(r),0,1))

#create a Transition object from the raster
tr <- transition(r,mean,8)

#asymmetric
ncf <- function(x) max(x) - x[1] + x[2]
tr2 <- transition(r,ncf,8, symm=FALSE)

#create two sets of coordinates
sP1 <- cbind(c(65,5,-65),c(55,35,-35))
sP2 <- cbind(c(50,15,-40),c(80,20,-5))

#from and to identical
costDistance(tr,sP1)
costDistance(tr2,sP1)

#from and to different
costDistance(tr,sP1,sP2)
costDistance(tr2,sP1,sP2)

}
\keyword{ spatial }

