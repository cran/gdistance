# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  January 2009
# Version beta
# Licence GPL v3

#TODO check if coordinate systems are equal.
#TODO check if bounding box of coordinates falls inside bb of transition
#TODO coordinates in same cell: distance = 0

setGeneric("resistanceDistance", function(transition, coords) standardGeneric("resistanceDistance"))

setMethod("resistanceDistance", signature(transition = "TransitionLayer", coords = "Coords"), def = function(transition, coords) 
	{
		return(.rD(transition, coords))
	}
)

.rD <- function(transition, coords){
		if(class(transitionMatrix(transition)) != "dsCMatrix"){stop("symmetric transition matrix required (dsCMatrix) in TransitionLayer object")}
		coords <- .coordsToMatrix(coords)

		rd <- matrix(NA,nrow=length(coords[,1]),ncol=length(coords[,1]))
		rownames(rd) <- rownames(coords)
		colnames(rd) <- rownames(coords)
		allFromCells <- cellFromXY(transition, coords)
		
		if(!all(!is.na(allFromCells))){
			warning("some coordinates not found and omitted")
			allFromCells <- allFromCells[!is.na(allFromCells)]
		}

		transition <- .transitionSolidify(transition)		
		fromCells <- allFromCells[allFromCells %in% transitionCells(transition)]
		if (length(fromCells) < length(allFromCells)) 
		{
			warning(length(fromCells)," out of ",length(allFromCells)," locations were found in the fully connected transition matrix. NAs introduced.")
		}
		else{}
		fromCells <- unique(allFromCells)
		Lr <- .Laplacian(transition)
		n <- max(Lr@Dim)
		Lr <- Lr[-n,-n]
		n <- max(Lr@Dim)
		C <- 1e-300 * (n + 1) #This should avoid too big floating points as "Voltage differences", but give a number that can still be divided by n+1
		Lplus <- matrix(ncol=length(fromCells),nrow=length(fromCells))
		index <- match(fromCells,transitionCells(transition))
		Lr <- Cholesky(Lr)
		for (i in 1:length(fromCells))
		{
			ei <- matrix((-C/(n+1)), ncol=1, nrow=n)
			ei[index[i],] <- C-(C/(n+1))
			xi <- solve(Lr,ei)
			#xi <- as.vector(xi)
			#Lplusallrows <- c(xi-sum(xi/(n+1)),(sum(xi)/(n+1))) This is not necessary and with big floating points it may break
			#Lplus[,i] <- Lplusallrows[index]
			Lplus[,i] <- c(as.vector(xi),0)[index]
		}
		Lplus <- Lplus / C
		rdSS <- (-2*Lplus + matrix(diag(Lplus),nrow=length(fromCells),ncol=length(fromCells)) 
			+ t(matrix(diag(Lplus),nrow=length(fromCells),ncol=length(fromCells)))) * sum(transitionMatrix(transition))
		index1 <- which(allFromCells %in% fromCells)
		index2 <- match(allFromCells[allFromCells %in% fromCells],fromCells)
		rd[index1,index1] <- rdSS[index2,index2]
		rd <- as.dist(rd)
		attr(rd, "method") <- "resistance"
		return(rd)
}