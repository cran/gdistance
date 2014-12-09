\name{Transition* classes}
\docType{class}
\alias{TransitionLayer-class}
\alias{TransitionStack-class}
\alias{TransitionData-class}
\alias{Transition-class}
\alias{coerce,TransitionLayer,sparseMatrix-method}
\alias{coerce,TransitionLayer,RasterLayer-method}
\alias{coerce,TransitionLayer,TransitionStack-method}
\alias{coerce,TransitionLayer,TransitionData-method}
\alias{coerce,RasterLayer,TransitionLayer-method}
\alias{coerce,TransitionData,sparseMatrix-method}
\alias{show,TransitionLayer-method}
\alias{show,TransitionStack-method}

\alias{==,TransitionLayer,TransitionLayer-method}
\alias{==,TransitionStack,TransitionStack-method}

\title{Transition* classes}
\description{
TransitionLayer and TransitionStack (or Transition*) are the core classes of the package gdistance. 
They are the main input into the functions to calculate distances and routes.

An object of the class TransitionLayer contains two main elements:
a. a transition matrix with transition values between connected cells in a raster - an object of class sparseMatrix (package Matrix); 
b. information on the extent, resolution and projection of the underlying raster - an object of class Raster (package raster).

All slots belong to these two elements from other package, except two additional slots:
1. slot transitionCells, which is only used internally in the package;
2. slot matrixValues indicates if the nonzero values of the transition matrix contains conductance or resistance values.

Class TransitionStack contains one or more transition matrices.

Class Transition is the union of TransitionLayer and TransitionStack.

}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Transition", nrows, ncols, xmin, xmax, ymin, ymax, projection)}.
}
\section{Slots}{
	\describe{
    \item{\code{transitionMatrix}:}{Object of class \code{"sparseMatrix"}}
    \item{\code{transitionCells}:}{Object of class \code{"integer"}}
    \item{\code{matrixValues}:}{Object of class \code{"character"}}
    \item{\code{ncols}:}{Object of class \code{"integer"}}
    \item{\code{nrows}:}{Object of class \code{"integer"}}
    \item{\code{crs}:}{Object of class \code{"CRS"} (sp package)}	
    \item{\code{extent}:}{Object of class \code{"Extent"}}
	\item{\code{layernames}:}{Object of class \code{"vector"}}
	
  }
}
\section{Extends}{
Class \code{"\linkS4class{Raster}"}.

}
%\section{Methods}{
%  \describe{
%    \item{Arith}{\code{signature(e1 = "Transition", e2 = "Transition")} }
%    \item{coerce}{\code{signature(from = "Transition", to = "sparseMatrix")}}
%    \item{initialize}{\code{signature(.Object = "Transition")}}
%    \item{Math}{\code{signature(x = "Transition")} }
%    \item{show}{\code{signature(object = "Transition")} }
%	 }
%}

\author{ Jacob van Etten <jacobvanetten@yahoo.com> }
\examples{
showClass("TransitionLayer")

tr <- new("TransitionLayer", nrows=as.integer(36), ncols=as.integer(18),
		extent=extent(c(xmin=-180,xmax=180, ymin=-90,ymax=90)),
		crs=CRS("+proj=longlat +datum=WGS84"))

tr <- new("TransitionLayer",nrows=as.integer(36),ncols=as.integer(18),
        extent=extent(c(xmin=-180, xmax=180, ymin=-90,ymax=90)),
		crs=CRS(""))
}
\keyword{classes}
