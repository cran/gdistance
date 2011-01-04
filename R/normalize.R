# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  March 2010
# Version beta
# Licence GPL v3

setGeneric("normalize", function(transition, ...) standardGeneric("normalize"))

setMethod("normalize", signature(transition = "TransitionLayer"), def = function(transition, method="row")
	{
		return(.normalize(transition, method))
	}
)

.normalize <- function(transition, method)
	{
		tr <- transitionMatrix(transition)
		if(!(method %in% c("row","col","symm"))){stop("invalid method argument")}
		if(method=="symm")
		{
			tr <- t(tr) 
			rs <- (rowSums(tr)^-.5)
			rs[rs == Inf] <- 0
			tr <- tr * rs
			tr <- t(tr)
			cs <- colSums(tr)^-.5
			cs[cs == Inf] <- 0
			tr <- tr * rs
		}
		if(method=="row")
		{
			rs <- 1 / rowSums(tr)
			rs[rs == Inf] <- 0
			tr <- tr * rs
		}
		if(method=="col")
		{
			rs <- 1 / colSums(tr)
			rs[rs == Inf] <- 0
			tr <- tr * rs
		}
		transitionMatrix(transition) <- tr
		return(transition)
	}