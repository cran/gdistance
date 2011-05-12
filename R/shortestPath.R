# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  January 2009, code added January 2010
# Version 1.0
# Licence GPL v3

#check if Transition and RasterLayers coincide, etc.

setGeneric("shortestPath", function(transition, origin, goal, ...) standardGeneric("shortestPath"))

setMethod("shortestPath", signature(transition = "TransitionLayer", origin = "Coords", goal = "Coords"), 
	def = function(transition, origin, goal, output="TransitionLayer")
	{
		origin <- .coordsToMatrix(origin)
		goal <- .coordsToMatrix(goal)
		return(.shortestPath(transition, origin, goal, output))		
	}
)

.shortestPath <- function(transition, origin, goal, output)
{
	originCells <- cellFromXY(transition, origin)
	goalCells <- cellFromXY(transition, goal)
	indexOrigin <- originCells - 1
	indexGoal <- goalCells - 1
	if(isSymmetric(transitionMatrix(transition))) {mode <- "undirected"} else {mode <- "directed"}
	adjacencyGraph <- graph.adjacency(transitionMatrix(transition), mode=mode, weighted=TRUE)
	E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight

	shortestPaths <- get.shortest.paths(adjacencyGraph, indexOrigin, indexGoal) #only first in indexOrigin is taken into account!
	
	if(output=="TransitionLayer")
	{
		
		result <- transition
		transitionMatrix(result) <- Matrix(0, ncol=ncell(transition), nrow=ncell(transition))			
		for(i in 1:length(shortestPaths))
		{
			sPVector <- (shortestPaths[[i]] + 1) #igraph 6.0 -> change
			adj <- cbind(sPVector[-(length(sPVector))], sPVector[-1])
			adj <- rbind(adj,cbind(adj[,2], adj[,1]))
			transitionMatrix(result)[adj] <- 1/length(shortestPaths)
		}
	}

	if(output=="TransitionStack")
	{
		result <- transition
		transitionMatrix(result) <- Matrix(0, ncol=ncell(transition), nrow=ncell(transition))			
		for(i in 1:length(shortestPaths))
		{
			resultNew <- result
			sPVector <- (shortestPaths[[i]] + 1) #igraph 6.0 -> change
			adj <- cbind(sPVector[-(length(sPVector))], sPVector[-1])
			adj <- rbind(adj,cbind(adj[,2], adj[,1]))
			transitionMatrix(resultNew)[adj] <- 1/length(shortestPaths)
			result <- stack(result, resultNew)
		}
		result <- result[[2:nlayers(result)]]	
	}

	if(output=="SpatialLines")
	{
		linesList <- vector(mode="list", length=length(shortestPaths))
				
		for(i in 1:length(shortestPaths))
		{
			sPVector <- (shortestPaths[[i]] + 1) #igraph 6.0 -> change
			coords <- xyFromCell(transition, sPVector)
			linesList[[i]] <- Line(coords)
		}
		
		LinesObject <- Lines(linesList, ID = as.character(1:length(shortestPaths)))
		result <- SpatialLines(list(LinesObject), proj4string = CRS(projection(transition)))
	}

	return(result)

}