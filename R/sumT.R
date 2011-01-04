sumT <- function(transition1, transition2)
{
	if(matrixValues(transition1) == "conductance") 
	{
		transition1@transitionMatrix@x <- 1 / transition1@transitionMatrix@x
		transition1@transitionMatrix@x[transition1@transitionMatrix@x == Inf] <- 0
	}
	if(matrixValues(transition2) == "conductance") 
	{
		transition2@transitionMatrix@x <- 1 / transition2@transitionMatrix@x
		transition2@transitionMatrix@x[transition2@transitionMatrix@x == Inf] <- 0
	}
	newTransition <- transition1 + transition2
	newTransition@transitionMatrix@x <- 1 / newTransition@transitionMatrix@x
	if(sum(newTransition@transitionMatrix@x == Inf) > 0) 
	{
		warning("Inf values introduced. Set to 0.")
		newTransition@transitionMatrix@x[newTransition@transitionMatrix@x == Inf] <- 0
	}
	matrixValues(newTransition) <- "conductance"
	return(newTransition)
}
