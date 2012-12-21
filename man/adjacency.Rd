\name{adjacencyFromTransition}

\alias{adjacencyFromTransition}

\title{Adjacent cells}

\description{
Identify pairs of cells that are adjacent.  
}

\usage{
adjacencyFromTransition(x) 
}

\arguments{
  \item{x}{TransitionLayer}
}

\details{
Extracts the indices of those cells that are connected (e.g. cells {i,j} that have a non-zero value in the transition matrix). 

Cell numbers are unique indices of cells in the original grid. 
By convention, cell numbers start with 1 in the upper-left corner of the grid and increase from left to right and from top to bottom.
}

\value{
A two column matrix with each row containing a pair of adjacent cells. 
}

\seealso{
\link[raster:adjacent]{adjacent}
}

\author{Jacob van Etten \email{jacobvanetten@yahoo.com}}


\examples{
r <- raster(nrows=6, ncols=7, xmn=0, xmx=7, ymn=0, ymx=6, crs="+proj=utm +units=m")
r[] <- runif(6*7)
T <- transition(r, function(x) 1/mean(x), 8) 
adjacencyFromTransition(T)
}

\keyword{spatial}
