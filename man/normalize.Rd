\name{normalize}
\alias{normalize}
\alias{normalize,TransitionLayer-method}
\title{ normalize }
\description{  Normalize the transition matrix.
}

\usage{
normalize(x, ...)
}
\arguments{
  \item{x}{object of class Transition*)}
  \item{...}{optional argument \code{method} (see below)}
 }
\details{
Normalization of the weighted adjacency matrix in the \code{Transition*} object.
Matrix values are divided by their respective row-sums, column-sums, or the product of the square-roots of both (symmetric normalization).
The default method is row-normalization. 
To use the other normalization methods, users can set the optional \code{method} argument to either "col" or "symm".
For random walk calculations a symmetric matrix is needed (\code{method = "symm"}).
}
\value{
  \code{TransitionLayer} object
}
\references{
von Luxburg, U. 2007. A tutorial on spectral clustering. \emph{Statistics and Computing} 17(4), 395-416.
\url{http://arxiv.org/PS_cache/arxiv/pdf/0711/0711.0189v1.pdf}

Chung, F. 1997. \emph{Spectral Graph Theory}. Conference Board of the Mathematical Sciences, Washington.
}
\author{ Jacob van Etten \email{jacobvanetten@yahoo.com} }
\examples{
r <- raster(ncol=36,nrow=18)
r <- setValues(r,rep(1,times=ncell(r)))
tr <- transition(r, mean, directions=8)

tr1 <- normalize(tr)
tr2 <- normalize(tr, method="symm")
}
\keyword{methods}
\keyword{spatial}