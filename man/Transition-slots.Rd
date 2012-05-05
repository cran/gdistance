\name{Transition slots}
\alias{transitionMatrix}
\alias{transitionMatrix<-}
\alias{matrixValues}
\alias{matrixValues<-}
\alias{transitionCells}
\alias{transitionData}
\alias{transitionMatrix,TransitionLayer,logical-method}
\alias{transitionMatrix,TransitionLayer,missing-method}
\alias{transitionMatrix<-,TransitionLayer,sparseMatrix-method}
\alias{transitionMatrix,TransitionData,missing-method}
\alias{transitionCells,TransitionLayer-method}
\alias{transitionCells,TransitionData-method}
\alias{matrixValues,TransitionLayer-method}
\alias{matrixValues<-,TransitionLayer,character-method}
\alias{matrixValues,TransitionStack-method}
\alias{matrixValues<-,TransitionStack,character-method}
\alias{transitionData,TransitionLayer-method}
\alias{transitionData,TransitionStack-method}



\title{Extract or change elements of Transition* objects}
\description{
These functions are to be used to access slots of Transition* objects.
}
\usage{
transitionMatrix(x, inflate)
transitionCells(x)
matrixValues(x)
}
\arguments{
  \item{x}{object of class Transition*}
  \item{inflate}{logical (default is TRUE)} 
}