# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  January 2009
# Version 1.0
# Licence GPL v3

setGeneric("transition", function(object, transitionFunction, directions, ...) standardGeneric("transition"))

setMethod("transition", signature(object = "RasterLayer"), def = function(object, transitionFunction, directions, symm=TRUE, ...)
		{
			if(class(transitionFunction)=="character") 
			{
				if(transitionFunction=="barriers")
				{
					return(.barriers(object, directions, intervalBreaks))
				}
				else{return(.TfromR(object, transitionFunction, directions, symm))}
			}
			else{return(.TfromR(object, transitionFunction, directions, symm))}
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

.barriers <- function(x, directions, intervalBreaks) {
	maxVal <- max(getValues(x))
	minVal <- min(getValues(x))
	
	Xmin <- transition(x, min, 8)
	Xmax <- transition(x, max, 8)
	index1 <- adjacency(x, 1:ncell(x), 1:ncell(x), directions)
	XminVals <- Xmin[index1]
	XmaxVals <- Xmax[index1]
	Xstack <- as(Xmin, "TransitionStack") * 0
	Xlayer <- Xmin * 0
	matrixValues(Xlayer) <- "resistance"
	for(i in 1:length(intervalBreaks))
	{
		index2 <- index1[XminVals < intervalBreaks[i] & XmaxVals > intervalBreaks[i],]
		Xlayer[index2] <- 1
		Xstack <- stack(Xstack,Xlayer)
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