# Author: Jacob van Etten
# Date: September 2010
# Version 1.0
# Licence GPL v3


if (!isGeneric("summary")) {
	setGeneric("summary", function(object, ...)
		standardGeneric("summary"))
}	

setMethod('summary', signature(object='TransitionLayer'), 
	function(object, ...) {
		stop("Not yet implemented")
	}
)

