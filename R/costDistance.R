# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  January 2009
# Version 1.0
# Licence GPL v3

#TODO check if coordinate systems are equal.
#TODO check if bounding box of coordinates falls inside bb of transition

setGeneric("costDistance", function(transition, fromCoords, toCoords) standardGeneric("costDistance"))

setMethod("costDistance", signature(transition = "TransitionLayer", fromCoords = "Coords", toCoords = "Coords"), def = function(transition, fromCoords, toCoords)
	{
		return(.cd(transition, fromCoords, toCoords))
	}
)

.cd <- function(transition, fromCoords, toCoords)
{
	fromCoords <- .coordsToMatrix(fromCoords)
	toCoords <- .coordsToMatrix(toCoords)

	fromCells <- cellFromXY(transition, fromCoords)
	toCells <- cellFromXY(transition, toCoords)
	
	if(!all(!is.na(fromCells))){
		warning("some coordinates not found and omitted")
		fromCells <- fromCells[!is.na(fromCells)]
	}
	
	if(!all(!is.na(toCells))){
		warning("some coordinates not found and omitted")
		toCells <- toCells[!is.na(toCells)]
	}
	
	costDist <- matrix(NA, nrow=length(fromCoords[,1]),ncol=length(toCoords[,1]))
	rownames(costDist) <- rownames(fromCoords)
	colnames(costDist) <- rownames(toCoords)
		
	if(isSymmetric(transitionMatrix(transition))) {m <- "undirected"} else {m <- "directed"}
	adjacencyGraph <- graph.adjacency(transitionMatrix(transition), mode=m, weighted=TRUE)
	
	E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight

	uniqueFromCells <- unique(fromCells)
	uniqueToCells <- unique(toCells)		
	shortestPaths <- matrix(nrow=length(uniqueFromCells),ncol=length(uniqueToCells))

	for (i in 1:length(uniqueFromCells))
	{
		shortestPaths[i,] <- shortest.paths(adjacencyGraph, uniqueFromCells[i]-1, mode="out")[,uniqueToCells]
	}

	index1 <- match(fromCells,uniqueFromCells)
	index2 <- match(toCells,uniqueToCells)
	costDist[] <- shortestPaths[index1,index2]
	return(costDist)
}

setMethod("costDistance", signature(transition = "TransitionLayer", fromCoords = "Coords", toCoords = "missing"), def = function(transition, fromCoords)
	{
		return(.cd2(transition, fromCoords))
	}
)

.cd2 <- function(transition, fromCoords)
{
		fromCoords <- .coordsToMatrix(fromCoords)
		fromCells <- cellFromXY(transition, fromCoords)

		if(!all(!is.na(fromCells))){
			warning("some coordinates not found and omitted")
			fromCells <- fromCells[!is.na(fromCells)]
		}
		
		costDist <- matrix(NA, nrow=length(fromCoords[,1]),ncol=length(fromCoords[,1]))
		rownames(costDist) <- rownames(fromCoords)
		colnames(costDist) <- rownames(fromCoords)
		
		if(isSymmetric(transitionMatrix(transition))) {m <- "undirected"} else {m <- "directed"}
		adjacencyGraph <- graph.adjacency(transitionMatrix(transition), mode=m, weighted=TRUE)
		
		E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight

		uniqueFromCells <- unique(fromCells)
		shortestPaths <- matrix(ncol=length(uniqueFromCells),nrow=length(uniqueFromCells))
		for (i in 1:length(uniqueFromCells))
		{
			shortestPaths[i,] <- shortest.paths(adjacencyGraph, uniqueFromCells[i]-1, mode="out")[,uniqueFromCells]
		}
		index <- match(fromCells,uniqueFromCells)
		costDist[] <- shortestPaths[index,index]
		if(m=="undirected") {costDist <- as.dist(costDist)}
		return(costDist)
}