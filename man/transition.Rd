\name{transition }
\alias{transition}
\alias{transition,RasterLayer-method}
\alias{transition,RasterBrick-method}
\title{ Create an object of the class Transition }
\description{
  Create a Transition object from a RasterLayer or RasterBrick object. Transition values are calculated with a user-defined function from the grid values.
}
\usage{
transition (x, transitionFunction, directions, ...)
}
\arguments{
  \item{x}{ \code{RasterLayer} or \code{RasterBrick} (raster package) }
  \item{transitionFunction}{ Function to calculate transition values from grid values}
  \item{directions}{ Directions in which cells are connected (4, 8, 16, or other), see \code{\link[raster:adjacent]{adjacent}} }
  \item{...}{ additional arguments, see methods 1 and 3 below }
}
\details{
Users may use one of three methods to construct a Transition* object with this function.

1) \code{TransitionLayer} from \code{RasterLayer}

\code{transition(x, transisitonFunction, directions, symm)}

When a symmetric transition matrix is required, the user should supply a transitionFunction f that obeys f(i,j) = f(j,i) (a commutative function). 
The function \code{transition} does no commutativity check.
To obtain an asymmetric transition matrix, a non-commutative function should be supplied and an additional argument `symm' should be set to FALSE.

2) \code{TransitionLayer} from \code{RasterBrick}

\code{transition(x, transitionFunction="mahal", directions)}

This method serves to summarize several layers of data in a single distance measure.
The distance between adjacent cells is the normalized reciprocal of the Mahalanobis distance (mean distance / (mean distance + distance ij).

3) \code{TransitionStack} from \code{RasterLayer}

In contrast with the above methods, this method produces resistance matrices by default.

a) Continuous variables - barriers

\code{transition(x, transitionFunction="barriers", directions, symm, intervalBreaks)}

This method creates a \code{TransitionStack} with each layer containing a discrete boundary between areas in \code{x}.
Areas are defined by intervals in \code{x}.
The argument \code{intervalBreaks} is a vector of interval breaks corresponding to the values in \code{x}.
If between a pair of cells i and j, min(i,j) < break AND max(i,j) > break, then the value ij in the transition matrix becomes 1.
All other values in the transition matrix remain 0.
The package classInt offers several methods to define intervals.
If symm is changed from the default (TRUE) to "up" or "down", it will give either only the upslope (symm="up") or downslope (symm="down") barriers.

b) Categorical variables - barriers 

\code{transition(x, transitionFunction="barriers", directions)}

In this case, areas are defined as categories in the input raster.
A raster with a categorical variable can be created with \code{asFactor()}.
The layers of the resulting TransitionStack contain all possible combinations of categories.
Which layer contains the combination of categories i and j out of n categories, can be determined with these formulae:

if \code{symm} is \code{TRUE}: layer(i,j) = n*(j-1) - j*(j-1)/2 + i-j.
if \code{symm} is \code{FALSE} and i>j: layer(i,j) = ((n*(j-1) - j*(j-1)/2 + i-j) * 2) - 1.
if \code{symm} is \code{FALSE} and i<j: layer(i,j) = (n*(j-1) - j*(j-1)/2 + i-j) * 2.

c) Categorical variables - areas

\code{transition(x, transitionFunction="areas", directions)}

Here, areas are also a categorical variable (see under 3b).
The layers in the resulting TransitionStack represent each one area.
Connections between two cells which are each inside the area are set to 1.
Connections between a cell inside and a cell outside the area are set to 0.5.
Connections between two cells outside the area are set to 0.


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
  
  #Barriers - interval breaks
  tr3 <- transition(r, "barriers", 8, intervalBreaks=c(0.25, 0.5, 0.75))
  nlayers(tr3)

  #Barriers - categories
  r <- round(r * 2)
  r <- asFactor(r)
  tr4 <- transition(r, "barriers", 8)
  nlayers(tr4)
  plot(raster(tr4[[2]]))
  
  #Areas
  r <- round(r * 2)
  r <- asFactor(r)
  tr5 <- transition(r, "areas", 8)
  nlayers(tr5)
  plot(raster(tr5[[2]]))
}
\keyword{ spatial }