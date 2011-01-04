\name{pathInc}
\alias{pathInc}
\alias{pathInc,TransitionLayer,Coords,Coords,missing,character,missing,missing-method}
\alias{pathInc,TransitionLayer,Coords,Coords,missing,character,numeric,missing-method}
\alias{pathInc,TransitionLayer,Coords,Coords,missing,character,missing,TransitionLayer-method}
\alias{pathInc,TransitionLayer,Coords,Coords,missing,character,numeric,TransitionLayer-method}
\alias{pathInc,TransitionLayer,Coords,Coords,missing,character,missing,TransitionStack-method}
\alias{pathInc,TransitionLayer,Coords,Coords,missing,character,numeric,TransitionStack-method}

\title{Incidence of paths from a common origin: overlap and non-overlap}
\description{
Calculate the overlap and non-overlap of two paths departing from a common origin. Two algorithms to are available: random walk and randomised shortest paths.}
\usage{
pathInc(transition, origin, fromCoords, toCoords, type, theta, weight)
}
\arguments{
  \item{transition}{transition matrix (class Transition)}
  \item{origin}{coordinates of the origin (one point location, SpatialPoints, matrix or numeric class)}  
  \item{fromCoords}{coordinates of the destinations (SpatialPoints, matrix or numeric class)}
  \item{toCoords}{second set of coordinates of the destinations (only "missing" implemented)}
  \item{type}{a vector of one or two calculation types c("divergent","joint")}  
  \item{theta}{value > 0 and < 20 (randomised shortest paths) or missing (random walk)}
  \item{weight}{weight matrix -- Reciprocals of the non-zero values are used as weights. If missing, reciprocals of the transition matrix are used.}
}
\details{
This function calculates the overlap and non-overlap in the random walk trajectories or randomised shortest paths originating in a single point to sets of other locations.
First, probabilities of passage are calculated for each cell-to-cell connection. 
The shared path is calculated as:

Pjoint = Pi * Pj 

Divergence is calculated as: 

Pdivergent = max(0, max(Pi, Pj) * (1 - min(Pi, Pj)) - min(Pi, Pj))

Both are summed over all connections ij in the transition matrix.

}
\value{
dist object or list of dist objects
}
\references{
McRae B.H., B.G. Dickson, and T. Keitt. 2008. Using circuit theory to model connectivity in ecology, evolution, and conservation. \emph{Ecology} 89:2712-2724.
  
Saerens M., L. Yen, F. Fouss, and Y. Achbany. 2009. Randomized shortest-path problems: two related models. \emph{Neural Computation}, 21(8):2363-2404.

van Etten, J., and R.J. Hijmans. 2010. A geospatial modelling approach integrating archaeobotany and genetics to trace the origin and dispersal of domesticated plants. \emph{PLoS ONE} 5(8): e12060.
}
\author{ Jacob van Etten \email{jacobvanetten@yahoo.com}. Implementation of randomised shortest paths based on Matlab code by Marco Saerens. }
\examples{
#Create TransitionLayer
r <- raster(ncol=36,nrow=18)
r <- setValues(r,rep(1,times=ncell(r)))
tr <- transition(r,mean,directions=4)

#Two different types of correction are required
trR <- geoCorrection(tr, type="r", multpl=FALSE)
trC <- geoCorrection(tr, type="c", multpl=FALSE)

#Create TransitionStack
ts <- stack(trR, trR)

#Points for origin and coordinates between which to calculate path (non)overlaps
sP0 <- SpatialPoints(cbind(0,0))
sP1 <- SpatialPoints(cbind(c(65,5,-65),c(-55,35,-35)))

#Randomised shortest paths
#rescaling is needed: exp(-theta * trC) should give reasonable values
trC <- trC / median(transitionMatrix(trC)@x) #divide by median of the non-zero values
pathInc(trC, origin=sP0, fromCoords=sP1, type=c("divergent","joint"), theta=2)

#Random walk
pathInc(trR, origin=sP0, fromCoords=sP1, type=c("divergent","joint"))

#TransitionStack as weights
pathInc(trR, origin=sP0, fromCoords=sP1, type=c("divergent","joint"), weight=ts)

}
\keyword{ spatial }

