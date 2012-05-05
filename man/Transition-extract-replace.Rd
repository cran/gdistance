\name{Transition* methods}
\docType{methods}
\alias{[,TransitionLayer,numeric,numeric,missing-method}
\alias{[,TransitionLayer,matrix,missing,missing-method}
\alias{[<-,TransitionLayer,matrix,missing,ANY-method}
\alias{[<-,TransitionLayer,numeric,numeric,ANY-method}
\alias{[[,TransitionStack,numeric,missing-method}
\alias{[[<-,TransitionStack,numeric,missing,TransitionData-method}

\title{ Extracting and replacing: class Transition}
\description{
Methods for functions \code{[} and \code{[<-} for object of the class TransitionLayer.
Methods for functions \code{[[} and \code{[[<-} for object of the class TransitionStack.

Also see \code{\link[gdistance]{adjacencyFromTransition}}.
}
\examples{
  #Create a new raster and set all its values to unity.
  r <- raster(nrows=18, ncols=36) 
  r <- setValues(r,rep(1,ncell(r)))

  #Create TransitionLayer objects
  tr1 <- transition(r,mean,4)
  tr2 <- tr1

  #Extracting and replacing
  tr1[cbind(1:9,1:9)] <- tr2[cbind(1:9,1:9)]
  tr1[1:9,1:9] <- tr2[1:9,1:9]
  tr1[1:5,1:5]
}
\keyword{methods}