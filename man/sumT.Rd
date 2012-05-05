\name{sumReciprocal}
\alias{sumReciprocal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Reciprocal of the sum of the reciprocals of conductance values in Transition* objects}
\description{
Reciprocal of the sum of the reciprocals of conductance Transition* objects }
\usage{
sumReciprocal(x1, x2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{ \code{TransitionLayer} object}
  \item{x2}{ \code{TransitionLayer} object}  
}
\details{
To calculate the total resistance of two resistors that are serially connected, we should add their resistance values. However, if we work with conductance values, we need to take the reciprocal of the summed reciprocals of the conductance values. This function does that when adding two TransitionLayers with conductance values (\code{matrixValues(tr) == "conductance"}). 

For a TransitionLayer with resistance values (\code{matrixValues(tr) == "resistance"}), the function will not take reciprocals for that object, but will still take a reciprocal for the final product (which will consequently have conductance values).
}
\value{
TransitionLayer object containing conductance values.
}
\author{ Jacob van Etten \email{jacobvanetten@yahoo.com} }
\examples{
  #Create a new raster and set all its values to unity.
  raster <- raster(nrows=18, ncols=36) 
  raster <- setValues(raster,rep(1,ncell(raster)))

  #Create TransitionLayer objects
  tr1 <- transition(raster,mean,4)
  tr2 <- tr1
  matrixValues(tr1)
  
  #Set one to resistance
  matrixValues(tr2) <- "resistance"
  
  #Sum the two objects
  sumReciprocal(tr1,tr2)
}
\keyword{ spatial }

