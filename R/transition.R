# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  January 2009
# Version 1.0
# Licence GPL v3

setGeneric("transition", function(object, transitionFunction, directions, ...) standardGeneric("transition"))

setMethod("transition", signature(object = "RasterLayer"), def = function(object, transitionFunction, directions, symm=TRUE, intervalBreaks=NULL)
		{
			if(class(transitionFunction)=="character") 
			{
				if(transitionFunction != "barriers" & transitionFunction != "areas")
				{
					stop("argument transitionFunction invalid")
				}
				if(transitionFunction=="barriers")
				{
					return(.barriers(object, directions, symm, intervalBreaks))
				}
				if(transitionFunction=="areas")
				{
					return(.areas(object, directions))
				}
			} else {
				return(.TfromR(object, transitionFunction, directions, symm))
			}
		}
)

.TfromR <- function(object, transitionFunction, directions, symm)
{
			transition <- new("TransitionLayer",nrows=nrow(object),
				ncols=ncol(object),xmin=xmin(object),xmax=xmax(object),ymin=ymin(object),ymax=ymax(object),
				projection=projection(object, asText=FALSE))
			transitionMatr <- transitionMatrix(transition)
			adj <- adjacency(object,which(!is.na(getValues(object))),which(!is.na(getValues(object))),directions=directions)
			if(symm){adj <- adj[adj[,1] < adj[,2],]}
			dataVals <- cbind(getValues(object)[adj[,1]],getValues(object)[adj[,2]])
			transition.values <- apply(dataVals,1,transitionFunction)
			if(!all(transition.values>=0)){warning("transition function gives negative values")}
			transitionMatr[adj] <- as.vector(transition.values)
			if(symm)
			{
				transitionMatr <- forceSymmetric(transitionMatr)
			}
			transitionMatrix(transition) <- transitionMatr
			matrixValues(transition) <- "conductance"
			return(transition)
}

.barriers <- function(x, directions, symm, intervalBreaks) {

	Xlayer <- new("TransitionLayer",nrows=nrow(x),
			ncols=ncol(x),xmin=xmin(x),xmax=xmax(x),ymin=ymin(x),ymax=ymax(x),
			projection=projection(x, asText=FALSE))
	matrixValues(Xlayer) <- "resistance"
	Xstack <- as(Xlayer, "TransitionStack") * 0
	#Xstack@transition <- vector(list,...)
	
	if(x@data@isfactor) {

		vals <- unlist(x@data@attributes[[1]])
		n <- length(vals)
		
		if(symm)
		{
			maxn <- (n^2 - n)/2
			for(i in 1:maxn)
			{
				j <- .matrIndex(i,n)
				XlayerNew <- Xlayer
				cells1 <- which(getValues(x) == vals[j[1]])
				cells2 <- which(getValues(x) == vals[j[2]])
				adj1 <- adjacency(x, cells1, cells2, directions)
				adj2 <- adjacency(x, cells2, cells1, directions)
				adj <- rbind(adj1,adj2)
				XlayerNew[adj] <- 1
				Xstack <- stack(Xstack, XlayerNew)
			}
		} else {
			maxn <- (n^2 - n)/2
			for(i in 1:maxn)
			{
				j <- .matrIndex(i,n)
				XlayerNew1 <- Xlayer
				XlayerNew2 <- Xlayer
				cells1 <- which(getValues(x) == vals[j[1]])
				cells2 <- which(getValues(x) == vals[j[2]])
				adj1 <- adjacency(x, cells1, cells2, directions)
				adj2 <- adjacency(x, cells2, cells1, directions)
				XlayerNew1[adj1] <- 1
				XlayerNew2[adj2] <- 1				
				Xstack <- stack(Xstack, XlayerNew1, XlayerNew2)
			}
		}
	
	} else {
	
		Xmin <- transition(x, min, directions)
		Xmax <- transition(x, max, directions)
		index1 <- adjacency(x, 1:ncell(x), 1:ncell(x), directions)
		XminVals <- Xmin[index1]
		XmaxVals <- Xmax[index1]

		if(symm == TRUE)
		{
			for(i in 1:length(intervalBreaks))
			{
				index2 <- index1[XminVals < intervalBreaks[i] & XmaxVals > intervalBreaks[i],]
				XlayerNew <- Xlayer
				XlayerNew[index2] <- 1
				Xstack <- stack(Xstack,XlayerNew)
			}
		}
		if(symm=="up" | symm=="down"){stop("not implemented yet")}
	}
	
	Xstack <- Xstack[[2:nlayers(Xstack)]]	
	return(Xstack)
}


.areas <- function(x, directions) {

	Xlayer <- new("TransitionLayer",nrows=nrow(x),
			ncols=ncol(x),xmin=xmin(x),xmax=xmax(x),ymin=ymin(x),ymax=ymax(x),
			projection=projection(x, asText=FALSE))
	matrixValues(Xlayer) <- "resistance"
	Xstack <- as(Xlayer, "TransitionStack") * 0
	#Xstack@transition <- vector(list,...)
	
	if(x@data@isfactor) {

		vals <- unlist(x@data@attributes[[1]])
		n <- length(vals)
		
			for(i in 1:n)
			{
				transitionFunction <- function(v) {return(sum(v == i) / 2)}
				XlayerNew <- .TfromR(x, transitionFunction, directions, symm=TRUE)
				Xstack <- stack(Xstack,XlayerNew)
			}
			
	} else {
		warning("not yet implemented for raster with non-factor variables. Contact author.")
	}
	Xstack <- Xstack[[2:nlayers(Xstack)]]	
	return(Xstack)
}



setMethod("transition", signature(object = "RasterBrick"), def = function(object, transitionFunction="mahal", directions)
		{
			if(transitionFunction != "mahal")
			{
				stop("only Mahalanobis distance method implemented for RasterBrick")
			}
			x <- cbind(1:ncell(object),getValues(object))
			x <- na.omit(x)
			dataCells <- x[,1]
			adj <- adjacency(object,dataCells,dataCells,directions=directions)
			x.minus.y <- x[match(adj[,1],x[,1]),-1]-x[match(adj[,2],x[,1]),-1]
			cov.inv <- solve(cov(x[,-1]))
			mahaldistance <- apply(x.minus.y,1,function(x){sqrt((x%*%cov.inv)%*%x)})
			mahaldistance <- mean(mahaldistance)/(mahaldistance+mean(mahaldistance))
			transitiondsC <- new("dsCMatrix", 
					p = as.integer(rep(0,ncell(object)+1)),
					Dim = as.integer(c(ncell(object),ncell(object))),
					Dimnames = list(as.character(1:ncell(object)),as.character(1:ncell(object)))
			)
			transitiondsC[adj] <- mahaldistance
			transition <- new("TransitionLayer",nrows=nrow(object),ncols=
			ncol(object),xmin=xmin(object),xmax=xmax(object),ymin=ymin(object),ymax=ymax(object),
			projection=projection(object, asText=FALSE), matrixValues="conductance")
			transitionMatrix(transition) <- transitiondsC
			return(transition)
		}
)