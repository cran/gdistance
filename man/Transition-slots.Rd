\name{Transition slots}
\alias{transitionMatrix}
\alias{transitionMatrix<-}
\alias{matrixValues}
\alias{matrixValues<-}
\alias{transitionCells}
\alias{transitionData}
\alias{transitionMatrix,TransitionLayer-method}
\alias{transitionMatrix<-,TransitionLayer,sparseMatrix-method}
\alias{transitionCells,TransitionLayer-method}
\alias{matrixValues,TransitionLayer-method}
\alias{matrixValues<-,TransitionLayer,character-method}
\alias{matrixValues,TransitionStack-method}
\alias{matrixValues<-,TransitionStack,character-method}
\alias{transitionData,TransitionLayer-method}
\alias{transitionData,TransitionStack-method}
\alias{transitionMatrix,TransitionData-method}

\title{Extract or change elements of Transition* objects}
\description{
These functions are to be used to access slots of Transition* objects.
}
\usage{
transitionMatrix(transition)
transitionCells(transition)
matrixValues(transition)
}
\arguments{
  \item{transition}{object of class Transition*}
}