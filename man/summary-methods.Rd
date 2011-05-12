\name{Summary-methods}

\docType{methods}

\alias{Summary-methods}
\alias{sum,TransitionStack-method}
\alias{mean,TransitionStack-method}

\title{ Summary methods }

\description{
 The following summary methods are available:
 
 \code{mean, Median, max, min, range, prod, sum, any, all}
}


\note{

 These methods compute a summary statistic based on cell values of layers in a TransitionStack. 
 The result of these methods is always a single TransitionLayer. 
}

\value{a TransitionLayer}

\author{Jacob van Etten}

\examples{
#Create a new raster and set all its values to unity.
raster <- raster(nrows=18, ncols=36) 
raster <- setValues(raster,rep(1,ncell(raster)))

#Create a Transition object from the raster
tr <- transition(raster,mean,4)

trS <- stack(tr, tr*2)

#Apply a Summary method
trSum <- sum(trS)

#plot(raster(trMean))
}

\keyword{methods}
\keyword{spatial}
