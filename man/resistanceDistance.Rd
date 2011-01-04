\name{resistanceDistance}
\alias{ resistanceDistance }
\alias{resistanceDistance,TransitionLayer,Coords-method}

\title{ Resistance distance }
\description{
Calculates the resistance distance between points. }
\usage{
resistanceDistance (transition, coords)
}

\arguments{
  \item{transition}{ \code{TransitionLayer} object) }
  \item{coords}{ point locations coordinates (of SpatialPoints, matrix or numeric class)}  

}
\details{
  Resistance distance is equivalent to the expected random-walk commute time between nodes in a graph (Chandra et al. 1997).

  The function implements the algorithm given by Fouss et al. (2007).

  Before calculating resistance distances from a \code{TransitionLayer} object, see if you need to apply the function \code{\link{geoCorrection}}. 
}
\value{
distance matrix (S3 class dist or matrix)
}

\references{ 
Chandra, A.K., Raghavan, P., Ruzzo, W.L., Smolensy, R. & Tiwari, P. 1996. The electrical resistance of a graph captures its commute and cover times. \emph{Computational Complexity}, 6(4), 312-340.

Fouss, F., Pirotte, A., Renders, J.-M. & Saerens, M. 2007. Random-walk computation of similarities between nodes of a graph with application to collaborative recommendation. \emph{IEEE Transactions on Knowledge and Data Engineering}, 19(3), 355-369.

McRae, B.H. 2006. Isolation by resistance. \emph{Evolution} 60(8), 1551-1561. 
  
  \url{http://www.circuitscape.org}
}
\author{ Jacob van Etten \email{jacobvanetten@yahoo.com} }
\seealso{ \code{\link{geoCorrection}} }
\examples{
#Create a new raster and set all its values to unity.
raster <- raster(nrows=18, ncols=36) 
raster <- setValues(raster,rep(1,ncell(raster)))

#Create a Transition object from the raster
tr <- transition(raster,mean,4)

#Create two sets of coordinates
sP1 <- SpatialPoints(cbind(c(65,5,-65),c(55,35,-35)))
sP2 <- SpatialPoints(cbind(c(50,15,-40),c(80,20,-5)))

#Calculate the resistance distance between the points
resistanceDistance(tr,sP1)

}
\keyword{ spatial }