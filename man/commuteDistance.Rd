\name{commuteDistance}
\alias{ commuteDistance }
\alias{commuteDistance,TransitionLayer,Coords-method}

\title{ Commute-time distance }
\description{
Calculates the resistance distance between points. }
\usage{
commuteDistance (x, coords)
}

\arguments{
  \item{x}{ \code{TransitionLayer} object) }
  \item{coords}{ point locations coordinates (of SpatialPoints, matrix or numeric class)}

}
\details{
  This function calculates the expected random-walk commute time between nodes in a graph. It is defined as the effective distance (resistance distance) between the selected nodes multiplied by the volume of the graph, which is the sum of the conductance weights of all the edges in the graph (Chandra et al. 1997). The result represents the average number of steps that is needed to commute between the nodes during a random walk.

  The function implements the algorithm given by Fouss et al. (2007).

  Before calculating commute-time distances from a \code{TransitionLayer} object, see if you need to apply the function \code{\link{geoCorrection}}. 
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
r <- raster(nrows=18, ncols=36) 
r <- setValues(r,rep(1,ncell(raster)))

#Create a Transition object from the raster
tr <- transition(r,function(x) 1/mean(x),4)

#Create two sets of coordinates
sP1 <- SpatialPoints(cbind(c(65,5,-65),c(55,35,-35)))
sP2 <- SpatialPoints(cbind(c(50,15,-40),c(80,20,-5)))

#Calculate the resistance distance between the points
commuteDistance(tr,sP1)

}
\keyword{ spatial }