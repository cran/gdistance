# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  January 2009
# Version 1.0
# Licence GPL v3

setGeneric("accCost", function(transition, fromCoords) standardGeneric("accCost"))

setMethod("accCost", signature(transition = "TransitionLayer", fromCoords = "Coords"), def = function(transition, fromCoords)
	{
		fromCoords <- .coordsToMatrix(fromCoords) 
		fromCells <- cellFromXY(transition, fromCoords)
		if(!all(!is.na(fromCells))){
			warning("some coordinates not found and omitted")
			fromCells <- fromCells[!is.na(fromCells)]
		}
			tr <- transitionMatrix(transition)
		tr <- rBind(tr,rep(0,nrow(tr)))
		tr <- cBind(tr,rep(0,nrow(tr)))
	
		startNode <- nrow(tr) #extra node to serve as origin
		adjP <- cbind(rep(startNode, times=length(fromCells)), fromCells)
	
		tr[adjP] <- Inf
	
		adjacencyGraph <- graph.adjacency(tr, mode="directed", weighted=TRUE)
		E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight		
	
		shortestPaths <- shortest.paths(adjacencyGraph, v=startNode-1)[-startNode]

		result <- as(transition, "RasterLayer")
		result <- setValues(result, shortestPaths)	
		return(result)
	}
)

