# Author: Jacob van Etten jacobvanetten@yahoo.com
# Date :  September 2011
# Version 1.2
# Licence GPL v3

adjacencyFromTransition <- function(transition)
{
	transitionMatr <- as(transition,"sparseMatrix")
	transition.dgT <- as(transitionMatr,"dgTMatrix")
	adjacency <- cbind(transition.dgT@i+1,transition.dgT@j+1)
	return(adjacency)
}

