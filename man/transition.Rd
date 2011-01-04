\name{transition }
\alias{transition }
\alias{transition,RasterLayer-method}
\alias{transition,RasterBrick-method}
\title{ Create an object of the class Transition }
\description{
  Create a Transition object from a RasterLayer or RasterBrick object. Transition values are calculated with a user-defined function from the grid values.
}
\usage{
transition (object, transitionFunction, directions, ...)
}
\arguments{
  \item{object}{ \code{RasterLayer} or \code{RasterBrick} (raster package) }
  \item{transitionFunction}{ Function to calculate transition values from grid values}
  \item{directions}{ Directions in which cells are connected (4, 8, 16, or other), see \code{\link[raster]{adjacency}} }
  \item{...}{ Additional arguments, see below }
}
\details{
Users may use one of three methods to construct a Transition* object with this function.

1) \code{TransitionLayer} from \code{RasterLayer}

\code{transition(object, transisitonFunction, directions, symm)}

When a symmetric transition matrix is required, the user should supply a transitionFunction f that obeys f(i,j) = f(j,i) (a commutative function). 
The function \code{transition} does no commutativity check.
To obtain an asymmetric transition matrix, a non-commutative function should be supplied and an additional argument `symm' should be set to FALSE.

2) \code{TransitionLayer} from \code{RasterBrick}

\code{transition(object, transisitonFunction="mahal", directions)}

This method serves to summarize several layers of data in a single distance measure.
The distance between adjacent cells is the normalized reciprocal of the Mahalanobis distance (mean distance / (mean distance + distance ij).

3) \code{TransitionStack} from \code{RasterLayer}

\code{transition(object, transisitonFunction="barriers", directions, intervalBreaks)}

This method creates a \code{TransitionStack} with each layer containing a discrete boundary between areas in \code{object}.
Areas are defined by intervals in \code{object}.
The argument \code{intervalBreaks} is a vector of interval breaks corresponding to the values in \code{object}.

If between a pair of cells i and j, min(i,j) < break AND max(i,j) > break, then the value ij in the transition matrix becomes 1.
All other values in the transition matrix remain 0.
In contrast with the above methods, this method produces resistance matrices by default.
The package classInt offers several methods to define intervals.

}
\value{
Transition object
}

\author{ Jacob van Etten \email{jacobvanetten@yahoo.com}}
\examples{
  #Create a new raster and set all its values to unity.
  r <- raster(nrows=18, ncols=36) 
  r <- setValues(r, runif(ncell(r)))

  #Create a Transition object from the raster
  tr <- transition(r, transitionFunction=mean, directions=4)
  tr #the sparse matrix is of class dsCMatrix (symmetric)

  #Create an asymmetric transition matrix
  #first, an arbitrary non-commutative function 
  f <- function(x) max(x) - x[1] + x[2] 
  tr2 <- transition(r, f, 4, symm=FALSE)
  tr2 #the sparse matrix is of class dgCMatrix (=asymmetric)
}
\keyword{ spatial }

