\name{adjacencyFromTransition}

\alias{adjacencyFromTransition}

\title{Adjacent cells}

\description{
Identify pairs of cells that are adjacent.  
}

\usage{
adjacencyFromTransition(transition) 
}

\arguments{
  \item{transition}{TransitionLayer}
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
\link[raster:adjacency]{adjacency}
}

\author{Jacob van Etten \email{jacobvanetten@yahoo.com}}


\examples{

}

\keyword{spatial}
