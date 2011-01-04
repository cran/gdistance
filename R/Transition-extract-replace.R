# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  January 2009
# Version 1.0
# Licence GPL v3

setGeneric("transitionMatrix", function(transition) standardGeneric("transitionMatrix"))

setMethod ("transitionMatrix", signature(transition = "TransitionLayer"),
	function(transition){
		transition@transitionMatrix
	}
)

setMethod ("transitionMatrix", signature(transition = "TransitionData"),
	function(transition){
		transition@transitionMatrix
	}
)

setGeneric("transitionMatrix<-", function(transition, value) standardGeneric("transitionMatrix<-"))

setReplaceMethod ("transitionMatrix", signature(transition = "TransitionLayer", value = "sparseMatrix"),
	function(transition, value){
		if(dim(value)[1] != dim(value)[2]){stop("sparse matrix has to be square")}
		transition@transitionMatrix <- value
		return(transition)
	}
)

setGeneric("transitionCells", function(transition = "TransitionLayer") standardGeneric("transitionCells"))

setMethod ("transitionCells", signature(transition = "TransitionLayer"),
	function(transition){
		if(transition@transitionCells[1] != 0){return(transition@transitionCells)}
		else{return(1:ncell(transition))}
	}
)

setGeneric("matrixValues", function(transition = "TransitionLayer") standardGeneric("matrixValues"))

setMethod ("matrixValues", signature(transition = "TransitionLayer"),
	function(transition){
		transition@matrixValues
	}
)

setGeneric("matrixValues<-", function(transition, value) standardGeneric("matrixValues<-"))

setReplaceMethod ("matrixValues", signature(transition = "TransitionLayer", value = "character"),
	function(transition, value){
		if (value == "resistance" | value == "conductance") 
		{
			transition@matrixValues <- value
			return(transition)
		}
		else {stop("matrixValues can only be set to resistance or conductance")}
	}
)

setMethod("[", signature(x = "TransitionLayer", i="numeric", j="numeric", drop="missing"), function(x,i,j)
	{
		tm <- transitionMatrix(x)
		tm <- tm[i,j]
		return(tm)
	}
)

setMethod("[", signature(x = "TransitionLayer", i="matrix", j="missing", drop="missing"), function(x,i)
	{
		tm <- transitionMatrix(x)
		tm <- tm[i]
		return(tm)
	}
)

setMethod("[<-", signature(x = "TransitionLayer", i="matrix", j="missing", value="ANY"),
		function(x, i, value){
			tm <- transitionMatrix(x)
			tm[i] <- value
			x@transitionMatrix <- tm
			return(x)
		}
)

setMethod("[<-", signature(x = "TransitionLayer", i="numeric", j="numeric", value="ANY"),
		function(x, i, j, value)
		{
			tm <- transitionMatrix(x)
			tm[i,j] <- value
			transitionMatrix(x) <- tm
			return(x)
		}
)

setGeneric("transitionMatrix<-", function(transition, value) standardGeneric("transitionMatrix<-"))

setReplaceMethod ("transitionMatrix", signature(transition = "TransitionLayer", value = "sparseMatrix"),
	function(transition, value){
		if(dim(value)[1] != dim(value)[2]){stop("sparse matrix has to be square")}
		transition@transitionMatrix <- value
		return(transition)
	}
)

setMethod('nlayers', signature(object='TransitionStack'), 
	function(object)
	{
		return(length(object@transition)) 
    }
)


setMethod("[[", signature(x = "TransitionStack", i="numeric", j="missing"), function(x,i)
	{
		if (!(all(i %in% 1:nlayers(x)))){stop("indices should correspond to layers")}
		else
		{
			if(length(i)==1)
			{
				result <- new("TransitionLayer", nrows=nrow(x),ncols = ncol(x),xmin = xmin(x),xmax = xmax(x),
				ymin = ymin(x), ymax = ymax(x), projection=projection(x))
				result@transitionMatrix <- x@transition[[i]]@transitionMatrix
				result@transitionCells <- x@transition[[i]]@transitionCells
				result@matrixValues <- x@transition[[i]]@matrixValues
			}			
			if(length(i)>1)
			{
				result <- x
				result@transition <- x@transition[i]
			}
		}
		return(result)
	}
)

setMethod("[[<-", signature(x = "TransitionStack", i="numeric", j="missing", value="TransitionData"), function(x,i, value)
	{
		x@transition[[i]] <- value
		return(x)
	}
)

setGeneric("transitionData", function(transition) standardGeneric("transitionData"))

setMethod ("transitionData", signature(transition = "TransitionLayer"),
	function(transition){
		as(transition, "TransitionData")
	}
)

setMethod ("transitionData", signature(transition = "TransitionStack"),
	function(transition){
		transition@transition
	}
)