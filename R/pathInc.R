# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  January 2010
# Version 1.0
# Licence GPL v3

#transitionCells()

setGeneric("pathInc", function(transition, origin, fromCoords, toCoords, type, theta, weight) standardGeneric("pathInc"))

setMethod("pathInc", signature(transition = "TransitionLayer", origin = "Coords", fromCoords = "Coords", 
	toCoords = "missing", type="character", theta="missing", weight="missing"), 
	def = function(transition, origin, fromCoords, type)
	{
		origin <- .coordsToMatrix(origin)
		from <- .coordsToMatrix(fromCoords)
		preparedMatrix <- .preparationMatrix(transition, 0)
		preparedIndex <- .preparationIndex1(transition, origin, from, type)
		prepared <- c(preparedMatrix,preparedIndex,list(type=type))
		Intermediate <- .randomWalk(prepared)
		result <- .finishFlow(prepared, Intermediate)
		return(result)
	}
)

setMethod("pathInc", signature(transition = "TransitionLayer", origin = "Coords", fromCoords = "Coords", 
	toCoords = "missing", type="character", theta="numeric", weight="missing"), 
	def = function(transition, origin, fromCoords, type, theta)
	{
		if(theta < 0 | theta > 20 ) {stop("theta value out of range (between 0 and 20)")}
		origin <- .coordsToMatrix(origin)
		from <- .coordsToMatrix(fromCoords)
		preparedMatrix <- .preparationMatrix(transition, 0)
		preparedIndex <- .preparationIndex1(transition, origin, from, type)
		prepared <- c(preparedMatrix,preparedIndex,list(type=type))
		Intermediate <- .randomSP(prepared, theta)
		result <- .finishFlow(prepared, Intermediate)
		return(result)
	}
)

setMethod("pathInc", signature(transition = "TransitionLayer", origin = "Coords", fromCoords = "Coords", 
	toCoords = "missing", type="character", theta="missing", weight="TransitionLayer"), 
	def = function(transition, origin, fromCoords, type, weight)
	{
		origin <- .coordsToMatrix(origin)
		from <- .coordsToMatrix(fromCoords)
		preparedMatrix <- .preparationMatrix(transition, weight)
		preparedIndex <- .preparationIndex1(transition, origin, from, type)
		prepared <- c(preparedMatrix,preparedIndex,list(type=type))
		Intermediate <- .randomWalk(prepared)
		result <- .finishFlow(prepared, Intermediate)
		return(result)
	}
)

setMethod("pathInc", signature(transition = "TransitionLayer", origin = "Coords", fromCoords = "Coords", 
	toCoords = "missing", type="character", theta="numeric", weight="TransitionLayer"), 
	def = function(transition, origin, fromCoords, type, theta, weight)
	{
		if(theta < 0 | theta > 20 ) {stop("theta value out of range (between 0 and 20)")}
		origin <- .coordsToMatrix(origin)
		from <- .coordsToMatrix(fromCoords)
		preparedMatrix <- .preparationMatrix(transition, weight)
		preparedIndex <- .preparationIndex1(transition, origin, from, type)
		prepared <- c(preparedMatrix,preparedIndex,list(type=type))
		Intermediate <- .randomSP(prepared, theta)
		result <- .finishFlow(prepared, Intermediate)
		return(result)
	}
)

setMethod("pathInc", signature(transition = "TransitionLayer", origin = "Coords", fromCoords = "Coords", 
	toCoords = "missing", type="character", theta="missing", weight="TransitionStack"), 
	def = function(transition, origin, fromCoords, type, weight)
	{
		origin <- .coordsToMatrix(origin)
		from <- .coordsToMatrix(fromCoords)
		preparedMatrix <- .preparationMatrix(transition, weight)
		preparedIndex <- .preparationIndex1(transition, origin, from, type)
		prepared <- c(preparedMatrix,preparedIndex,list(type=type))
		Intermediate <- .randomWalk(prepared)
		result <- .finishFlowStack(prepared, Intermediate)
		return(result)
	}
)

setMethod("pathInc", signature(transition = "TransitionLayer", origin = "Coords", fromCoords = "Coords", 
	toCoords = "missing", type="character", theta="numeric", weight="TransitionStack"), 
	def = function(transition, origin, fromCoords, type, theta, weight)
	{
		if(theta < 0 | theta > 20 ) {stop("theta value out of range (between 0 and 20)")}
		origin <- .coordsToMatrix(origin)
		from <- .coordsToMatrix(fromCoords)
		preparedMatrix <- .preparationMatrix(transition, weight)
		preparedIndex <- .preparationIndex1(transition, origin, from, type)		
		prepared <- c(preparedMatrix,preparedIndex,list(type=type))
		Intermediate <- .randomSP(prepared, theta)
		result <- .finishFlowStack(prepared, Intermediate)
		return(result)
	}
)

.preparationIndex1 <- function(transition, origin, fromCoords, type)
{
		if(!all(type %in% c("divergent","joint"))) {stop("type can only have values \'joint\' and/or \'divergent\'")}

		originCell <- cellFromXY(transition, origin)

		if (!(originCell %in% transitionCells(transition))) {stop("the origin refers to a zero row/column in the transition matrix (unconnected)")} 

		allFromCells <- cellFromXY(transition, fromCoords)
		fromCells <- allFromCells[allFromCells %in% transitionCells(transition)]
		if (length(fromCells) < length(allFromCells)) 
		{
			warning(length(fromCells)," out of ",length(allFromCells[,1])," locations were found inside the transition matrix. NAs introduced.")
		}
		fromCells <- unique(fromCells)

		indexCoords <- match(fromCells,transitionCells(transition))
		indexOrigin <- match(originCell,transitionCells(transition))
		

		
		result <- list(transition=transition,
				fromCoords=fromCoords,
				allFromCells=allFromCells, 
				fromCells=fromCells,
				indexCoords=indexCoords, 
				indexOrigin=indexOrigin)
		return(result)
}
		
.preparationMatrix <- function(transition, weight)
{
		
		transition <- .transitionSolidify(transition)
		
		A <- as(transitionMatrix(transition),"lMatrix")
		A <- as(A,"dMatrix")
		AIndex <- as(A, "dgTMatrix")
		index1 <- cbind(transitionCells(transition)[as.integer(AIndex@i+1)],transitionCells(transition)[as.integer(AIndex@j+1)]) 
		index2 <- cbind(as.integer(AIndex@i+1),as.integer(AIndex@j+1))
		#if symmetric? index <- index[index[,1] < index[,2],]
		
		Size <- nrow(index1)
		
		if(class(weight) == "numeric")
		{
			R <- 1/transition[index2]
			R[R == Inf] <- 0
		}
		if(class(weight) == "TransitionLayer")
		{
			if(matrixValues(weight) == "conductance"){R <- 1/weight[index1]} else{R <- weight[index1]} 
			R[R == Inf] <- 0
		}
		if(class(weight) == "TransitionStack")
		{
			R <- matrix(nrow=nlayers(weight), ncol=length(index1[,1]))
			for(i in 1:nlayers(weight))
			{
				if(matrixValues(weight[[i]]) == "conductance"){R[i,] <- 1/weight[[i]][index1]} else{R[i,] <- weight[[i]][index1]}
			}
			R[R == Inf] <- 0
		}

		result <- list(transition=transition,
						index=index2,
						A=A,
						R=R,
						Size=Size)
		return(result)
}

.randomWalk <- function(prepared)
{
	transition <- prepared$transition
	indexCoords <- prepared$indexCoords
	indexOrigin <- prepared$indexOrigin
	fromCells <- prepared$fromCells
	index <- prepared$index
	Size <- prepared$Size
	A <- prepared$A

		
	L <- .Laplacian(transition)
	Lr <- L[-dim(L)[1],-dim(L)[1]]
	n <- max(Lr@Dim)
	Lr <- Cholesky(Lr)

	if(!(canProcessInMemory(transition, length(fromCells)*10))) #depending on memory availability, currents are calculated in a piecemeal fashion or all at once
	{
		filenm=rasterTmpFile()
		Flow <- raster(nrows=length(fromCells), ncols=Size)
		Flow <- writeStart(Flow, filenm, overwrite=TRUE)
		for(i in 1:(length(fromCells)))
		{
			matrixRow <- .currentM(L, Lr, A, n, indexOrigin, indexCoords[i], index)
			Flow <- writeValues(Flow, matrixRow, start=1)
		}
		Flow <- writeStop(Flow)
	}
	else
	{
		Flow <- matrix(nrow=Size,ncol=length(fromCells))
		for(i in 1:(length(fromCells)))
		{
			Flow[,i] <- .currentM(L, Lr, A, n, indexOrigin, indexCoords[i], index)
		}
	}
	return(Flow)
}


.randomSP <- function(prepared, theta)
{
	transition <- prepared$transition
	cj <- prepared$indexCoords
	ci <- prepared$indexOrigin
	index <- prepared$index
	Size <- prepared$Size
	fromCells <- prepared$fromCells
		
	tr <- transitionMatrix(transition)
	tc <- transitionCells(transition)
	
	trR <- tr
	trR@x <- 1 / tr@x

	nr <- dim(tr)[1] 
	Id <- Diagonal(nr) 
	rs <- rowSums(tr)
	rs[rs>0] <- 1/rs[rs>0]
	P <- tr * rs

	W <- trR
	W@x <- exp(-theta * trR@x) #zero values are not relevant because of next step exp(-theta * trR@x)
	W <- W * P 

	if(!(canProcessInMemory(transition, length(fromCells)*10))) 
	#this does not take into account the exact memory needed for matrix solving...
	{
		filenm=rasterTmpFile()
		Flow <- raster(nrows=length(cj), ncols=Size)
		Flow <- writeStart(Flow, filenm, overwrite=TRUE)
		for(i in 1:(length(cj)))
		{
			matrixRow <- transitionMatrix(.probPass(transition, Id, W, nr, ci, cj[i], tc, totalNet="net", output="TransitionLayer"))[index]
			Flow <- writeValues(Flow, matrixRow, 1)
		}
		Flow <- writeStop(Flow)
	}
	else
	{
		Flow <- matrix(nrow=Size,ncol=length(cj))
		for(i in 1:(length(cj)))
		{
			Flow[,i] <- transitionMatrix(.probPass(transition, Id, W, nr, ci, cj[i], tc, totalNet="net", output="TransitionLayer"))[index]
		}
	}
	return(Flow)
}	

.finishFlow <- function(prepared, Flow)
{
	fromCells <- prepared$fromCells
	allFromCells <- prepared$allFromCells
	fromCoords <- prepared$fromCoords
	type <- prepared$type
	
	Size <- prepared$Size
	R <- prepared$R
	if("divergent" %in% type) {divFlow <- matrix(ncol=length(fromCells),nrow=length(fromCells))}
	if("joint" %in% type) {jointFlow <- matrix(ncol=length(fromCells),nrow=length(fromCells))}	

	if(class(Flow) == "RasterLayer")
	{
		nr <- min(10,nrow(Flow))
		end <- ceiling(length(fromCells)/nr)

		
		nrows1 <- nr
		startrow1 <- 1
		dataRows1 <- getValues(Flow, row=startrow1, nrows=nrows1)
		dataRows1 <- matrix(dataRows1,nrow=ncol(Flow)) #rows are cell transitions, columns are locations

		for(j in 1:end)
		{
			if("divergent" %in% type) 
			{
				for(k in 1:nrows1)
				{
					index <- startrow1+k-1
					divFlow[startrow1:(startrow1+nrows1-1),index] <- colSums(matrix(pmax(pmax(dataRows1[,k],dataRows1) * 
					   (1-pmin(dataRows1[,k],dataRows1)) - pmin(dataRows1[,k],dataRows1), 0), nrow= Size) * R)
					#this fills the lower triangle only (plus some upper triangle blocks of size nr around the diagonal)
				}
			}
			if("joint" %in% type) 
			{
				for(l in 1:nrows1)
				{
					index <- startrow1+l-1
					jointFlow[startrow1:(startrow1+nrows1-1),index] <- colSums(matrix((dataRows1[,l] * dataRows1), nrow=Size) * R)

				}
			}
			if(j != end)
			{
				for(m in (j+1):end)
				{
					nrows2 <- min(nr, length(fromCells) - (m - 1) * nr)
					startrow2 <- (m-1)*nr+1
					dataRows2 <- getValues(Flow, row=startrow2, nrows=nrows2)
					dataRows2 <- matrix(dataRows2,nrow=ncol(Flow))
				
					if("divergent" %in% type) 
					{
						for(n in 1:nrows1)
						{
							index <- startrow1+n-1 
							divFlow[startrow2:(startrow2+nrows2-1),index] <- colSums(matrix(pmax(pmax(dataRows1[,n],dataRows2) * 
							(1-pmin(dataRows1[,n],dataRows2)) - pmin(dataRows1[,n],dataRows2), 0), nrow= Size) * R)
						}
					}
					if("joint" %in% type) 
					{
						for(o in 1:nrows1)
						{
							index <- startrow1+o-1
							jointFlow[startrow2:(startrow2+nrows2-1),index] <- colSums(matrix((dataRows1[,o] * dataRows2), nrow=Size) * R)
						}
					}
				}
				nrows1 <- min(nr, length(fromCells) - j * nr)
				startrow1 <- j*nr+1
				dataRows1 <- matrix(getValues(Flow, row=startrow1, nrows=nrows1),nrow=ncol(Flow))
			}
		}
	}
	if(class(Flow) == "matrix")
	{
		if("divergent" %in% type)
		{
			for(j in 1:(length(fromCells)))
			{
				divFlow[j,] <- colSums(matrix(pmax(pmax(Flow[,j],Flow) * (1-pmin(Flow[,j],Flow)) - pmin(Flow[,j],Flow), 0), nrow= Size) * R)
			}
		}
		if("joint" %in% type)
		{
			for(j in 1:(length(fromCells)))
			{
				jointFlow[j,] <- colSums((Flow[,j] * Flow) * R)
			}
		}
	}
	index1 <- which(allFromCells %in% fromCells)
	index2 <- match(allFromCells[allFromCells %in% fromCells], fromCells)
		
	if("divergent" %in% type)
	{
		divFl <- matrix(nrow=length(allFromCells),ncol=length(allFromCells))
		rownames(divFl) <- rownames(fromCoords)
		colnames(divFl) <- rownames(fromCoords)
		divFlow <- as.matrix(as.dist(divFlow, diag=TRUE))
		divFl[index1,index1] <- divFlow[index2,index2]
		divFl <- as.dist(divFl)
		attr(divFl, "method") <- "divergent path"
	}
	if("joint" %in% type)
	{
		jointFl <- matrix(nrow=length(allFromCells),ncol=length(allFromCells))
		rownames(jointFl) <- rownames(fromCoords)
		colnames(jointFl) <- rownames(fromCoords)
		jointFlow <- as.matrix(as.dist(jointFlow, diag=TRUE))
		jointFl[index1,index1] <- jointFlow[index2,index2]
		jointFl <- as.dist(jointFl)
		attr(jointFl, "method") <- "joint path"		
	}
	if(length(type) > 1) {return(list(divergent=divFl, joint=jointFl))}
	if(length(type) == 1 & type == "divergent") {return(divFl)}
	if(length(type) == 1 & type == "joint") {return(jointFl)}
} 

.finishFlowStack <- function(prepared, Flow) 
{
	fromCells <- prepared$fromCells
	allFromCells <- prepared$allFromCells
	fromCoords <- prepared$fromCoords
	type <- prepared$type
	
	Size <- prepared$Size
	R <- prepared$R
	if("divergent" %in% type) 
	{
		divFlow <- vector("list", nrow(R))
		divFlowi <- matrix(ncol=length(fromCells),nrow=length(fromCells))
		for(i in 1:nrow(R)){divFlow[[i]] <- divFlowi}
	}
	if("joint" %in% type) 
	{
		jointFlow <- vector("list", nrow(R))
		jointFlowi <- matrix(ncol=length(fromCells),nrow=length(fromCells))
		for(i in 1:nrow(R)){jointFlow[[i]] <- jointFlowi}
	}	
	if(class(Flow) == "RasterLayer")
	{
		nr <- min(10,nrow(Flow))
		end <- ceiling(length(fromCells)/nr)

		nrows1 <- nr
		startrow1 <- 1
		dataRows1 <- getValues(Flow, row=startrow1, nrows=nrows1)
		dataRows1 <- matrix(dataRows1,nrow=ncol(Flow)) #rows are cell transitions, columns are locations

		for(j in 1:end)
		{
			if("divergent" %in% type) 
			{
				for(k in 1:nrows1)
				{
					index <- startrow1+k-1
					for(i in 1:nrow(R))
					{
						divFlow[[i]][startrow1:(startrow1+nrows1-1),index] <- colSums(matrix(pmax(pmax(dataRows1[,k],dataRows1) * 
							(1-pmin(dataRows1[,k],dataRows1)) - pmin(dataRows1[,k],dataRows1), 0), nrow= Size) * R[i,])
					}
					#this fills the lower triangle only (plus some upper triangle blocks of size nr around the diagonal)
				}
			}
			if("joint" %in% type) 
			{
				for(l in 1:nrows1)
				{
					index <- startrow1+l-1
					for(i in 1:nrow(R))
					{
						jointFlow[[i]][startrow1:(startrow1+nrows1-1),index] <- colSums(matrix((dataRows1[,l] * dataRows1), nrow=Size) * R[i,])
					}
				}
			}
			if(j != end)
			{
				for(m in (j+1):end)
				{
					nrows2 <- min(nr, length(fromCells) - (m - 1) * nr)
					startrow2 <- (m-1)*nr+1
					dataRows2 <- getValues(Flow, row=startrow2, nrows=nrows2)
					dataRows2 <- matrix(dataRows2,nrow=ncol(Flow))
				
					if("divergent" %in% type) 
					{
						for(n in 1:nrows1)
						{
							index <- startrow1+n-1
							for(i in 1:nrow(R))
							{
								divFlow[[i]][startrow2:(startrow2+nrows2-1),index] <- colSums(matrix(pmax(pmax(dataRows1[,n],dataRows2) * 
									(1-pmin(dataRows1[,n],dataRows2)) - pmin(dataRows1[,n],dataRows2), 0), nrow= Size) * R[i,])
							}
						}
					}
					if("joint" %in% type) 
					{
						for(o in 1:nrows1)
						{
							index <- startrow1+o-1
							for(i in 1:nrow(R))
							{
								jointFlow[[i]][startrow2:(startrow2+nrows2-1),index] <- colSums(matrix((dataRows1[,o] * dataRows2), nrow=Size) * R[i,])
							}
						}
					}
				}
				nrows1 <- min(nr, length(fromCells) - j * nr)
				startrow1 <- j*nr+1
				dataRows1 <- matrix(getValues(Flow, row=startrow1, nrows=nrows1),nrow=ncol(Flow))
			}
		}
	}
	if(class(Flow) == "matrix")
	{
		if("divergent" %in% type)
		{
			for(j in 1:(length(fromCells)))
			{
				for(i in 1:nrow(R))
				{
					divFlow[[i]][j,] <- colSums(matrix(pmax(pmax(Flow[,j],Flow) * (1-pmin(Flow[,j],Flow)) - pmin(Flow[,j],Flow), 0), nrow= Size) * R[i,])
				}
			}
		}
		if("joint" %in% type)
		{
			for(j in 1:(length(fromCells)))
			{
				for(i in 1:nrow(R))
				{
					jointFlow[[i]][j,] <- colSums((Flow[,j] * Flow) * R[i,])
				}
			}
		}
	}
	index1 <- which(allFromCells %in% fromCells)
	index2 <- match(allFromCells[allFromCells %in% fromCells], fromCells)
	#convert to stack from here


	if("divergent" %in% type)
	{
		divFlnew <- matrix(nrow=length(allFromCells),ncol=length(allFromCells))
		rownames(divFlnew) <- rownames(fromCoords)
		colnames(divFlnew) <- rownames(fromCoords)
		
		for(i in 1:nrow(R))
		{
			divFl <- divFlnew
			divFlowi <- as.matrix(as.dist(divFlow[[i]], diag=TRUE))
			divFl[index1,index1] <- divFlowi[index2,index2]
			divFl <- as.dist(divFl)
			divFl <- as.vector(divFl)
			divFlow[[i]] <- divFl
		}
		divFlow <- as.data.frame(divFlow)
		colnames(divFlow) <- paste("div", 1:nrow(R), sep="")
	}
	if("joint" %in% type)
	{
		jointFlnew <- matrix(nrow=length(allFromCells),ncol=length(allFromCells))
		rownames(jointFlnew) <- rownames(fromCoords)
		colnames(jointFlnew) <- rownames(fromCoords)
		
		for(i in 1:nrow(R))
		{
			jointFl <- jointFlnew
			jointFlowi <- as.matrix(as.dist(jointFlow[[i]], diag=TRUE))
			jointFl[index1,index1] <- jointFlowi[index2,index2]
			jointFl <- as.dist(jointFl)
			jointFl <- as.vector(jointFl)
			jointFlow[[i]] <- jointFl
		}
		jointFlow <- as.data.frame(jointFlow)
		colnames(jointFlow) <- paste("joint", 1:nrow(R), sep="")
	}
	if(length(type) > 1) {return(cbind(divFlow, jointFlow))}
	if(length(type) == 1 & type == "divergent") {return(divFlow)}
	if(length(type) == 1 & type == "joint") {return(jointFlow)}
}