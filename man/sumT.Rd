\name{sumT}
\alias{sumT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summing transition objects}
\description{
Special methods to sum Transition objects }
\usage{
sumT(transition1, transition2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{transition1}{ \code{TransitionLayer} object}
  \item{transition2}{ \code{TransitionLayer} object}  
}
\details{
This function is meant to sum Transition objects taking the resistance values instead of the conductance values. In other words, the function calculates the reciprocal of the summed reciprocals. 

If \code{matrixValues(tr)} is equal to dQuote{conductance} the reciprocal of the values is taken (1 / conductance), if it is equal to dQuote{resistance}, the values are used. 

Zero values are ignored.

To treat a Transition object as containing resistance values instead of conductance values, use \code{matrixValues(tr)<-"resistance"}.
}
\value{
Transition object containing conductance values.
}
\note{
This function is useful to create barriers (with a certain resistance), with the rest having no resistance (value 0). This is difficult to realize with conductance values, as a matrix mostly filled with Inf values would be needed, which would be a dense matrix, occupying much memory.
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
  sumT(tr1,tr2)
}
\keyword{ spatial }

