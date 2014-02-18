### R code from vignette source 'gdistance.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: figure0
###################################################
library("gdistance")
rex <- raster(matrix(1,4,4))

a <- rep(c(1.3333),times=5)
b <- c(-1.3333, -0.6666, 0, 0.6666, 1.3333)

x1 <- c(-a, b)
x2 <- c(a, b)
y1 <- c(b, -a)
y2 <- c(b, a)
x <- cbind(x1,x2)
y <- cbind(y1,y2)

par(mfrow=c(1,3), mar = c(5,5,4,2) + 0.1,
    oma = c(0,0,0,0) + 0.1)

x4 <- transition(rex, mean, 4)
g4 <- graph.adjacency(transitionMatrix(x4), mode="undirected")
gridLayout <- xyFromCell(x4, 1:ncell(x4))
plot(g4,layout=gridLayout, edge.color="black", vertex.color="black", vertex.label=NA, main="4 cells")
for(i in 1:dim(x)[1]){lines(x[i,],y[i,], col="lightgray")}
plot(g4, layout=gridLayout, add=TRUE, edge.color="black", vertex.color="black", vertex.label=NA)

x8 <- transition(rex, mean, 8)
g8 <- graph.adjacency(transitionMatrix(x8), mode="undirected")
plot(g8,layout=gridLayout, edge.color="black", vertex.color="black", vertex.label=NA, main="8 cells")
for(i in 1:dim(x)[1]){lines(x[i,],y[i,], col="lightgray")}
plot(g8, layout=gridLayout, add=TRUE, edge.color="black", vertex.color="black", vertex.label=NA)

x16 <- transition(rex, mean, 16)
g16 <- graph.adjacency(transitionMatrix(x16), mode="undirected")
plot(g16, layout=gridLayout, edge.color="black", vertex.color="black", vertex.label=NA, , main="16 cells")
for(i in 1:dim(x)[1]){lines(x[i,],y[i,], col="lightgray")}
plot(g16,layout=gridLayout, add=TRUE, edge.color="black", vertex.color="black", vertex.label=NA)



###################################################
### code chunk number 2: gdistanc00
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 3: gdistance-1
###################################################
r <- raster(ncol=3,nrow=3)
r[] <- 1:ncell(r)
r


###################################################
### code chunk number 4: figure1
###################################################
plot(r, main="r")
text(r)


###################################################
### code chunk number 5: gdistance-2
###################################################
plot(r, main="r")
text(r)


###################################################
### code chunk number 6: gdistance-3
###################################################
library("gdistance")
r[] <- 1
tr1 <- transition(r, transitionFunction=mean, directions=8)


###################################################
### code chunk number 7: gdistance-4
###################################################
tr1


###################################################
### code chunk number 8: gdistance-5
###################################################
r[] <- runif(9)
ncf <- function(x) max(x) - x[1] + x[2]
tr2 <- transition(r, ncf, 4, symm=FALSE)
tr2


###################################################
### code chunk number 9: gdistance-6
###################################################
tr3 <- tr1*tr2
tr3 <- tr1+tr2
tr3 <- tr1*3
tr3 <- sqrt(tr1)


###################################################
### code chunk number 10: gdistance-7
###################################################
tr3[cbind(1:9,1:9)] <- tr2[cbind(1:9,1:9)]
tr3[1:9,1:9] <- tr2[1:9,1:9]
tr3[1:5,1:5]


###################################################
### code chunk number 11: gdistance-8 (eval = FALSE)
###################################################
## image(transitionMatrix(tr1))


###################################################
### code chunk number 12: figure2
###################################################
print(image(transitionMatrix(tr1)))


###################################################
### code chunk number 13: figure3
###################################################
plot(raster(tr3), main="raster(tr3)")


###################################################
### code chunk number 14: gdistance-9
###################################################
tr1C <- geoCorrection(tr1, type="c", multpl=FALSE)
tr2C <- geoCorrection(tr2, type="c", multpl=FALSE)


###################################################
### code chunk number 15: gdistance-10
###################################################
r3 <- raster(ncol=18, nrow=9)
r3 <- setValues(r3, runif(18*9)+5)
tr3 <- transition(r3, mean, 4)
tr3C <- geoCorrection(tr3, type="c", multpl=FALSE, scl=TRUE)
tr3R <- geoCorrection(tr3, type="r", multpl=FALSE, scl=TRUE)


###################################################
### code chunk number 16: gdistance-11
###################################################
CorrMatrix <- geoCorrection(tr3, type="r", multpl=TRUE, scl=TRUE)
tr3R <- tr3 * CorrMatrix


###################################################
### code chunk number 17: gdistance-12
###################################################
sP <- cbind(c(-100, 100, -100), c(-50, 50, 50))


###################################################
### code chunk number 18: gdistance-13
###################################################
costDistance(tr3C, sP)
commuteDistance(tr3R, sP)
rSPDistance(tr3R, sP, sP, theta=1e-12, totalNet="total")


###################################################
### code chunk number 19: gdistance-14
###################################################
origin <- SpatialPoints(cbind(0, 0))
rSPraster <- passage(tr3C, origin, sP[3,], theta=3)


###################################################
### code chunk number 20: figure4
###################################################
plot(rSPraster, main="rSPraster")


###################################################
### code chunk number 21: gdistance-15
###################################################
r1 <- passage(tr3C, origin, sP[1,], theta=1)
r2 <- passage(tr3C, origin, sP[3,], theta=1)
rJoint <- min(r1, r2) #Figure 6
rDiv <- max(max(r1, r2) * (1 - min(r1, r2)) - min(r1, r2), 0) #Figure 7


###################################################
### code chunk number 22: figure5
###################################################
plot(rJoint, main="rJoint")


###################################################
### code chunk number 23: figure6
###################################################
plot(rDiv, main="rDiv")


###################################################
### code chunk number 24: gdistance-17
###################################################
pathInc(tr3C, origin, sP)


###################################################
### code chunk number 25: figure7
###################################################
plot(function(x)exp(-3.5 * abs(x + 0.05)), -1, 1, xlab="slope", ylab="speed (m/s)")
lines(cbind(c(0,0),c(0,3.5)), lty="longdash")


###################################################
### code chunk number 26: gdistance-18
###################################################
r <- raster(system.file("external/maungawhau.grd", 
  package="gdistance"))


###################################################
### code chunk number 27: gdistance-19
###################################################
heightDiff <- function(x){x[2] - x[1]}
hd <- transition(r,heightDiff,8,symm=FALSE)
slope <- geoCorrection(hd, scl=FALSE)


###################################################
### code chunk number 28: gdistance-20
###################################################
adj <- adjacent(r, cells=1:ncell(r), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- exp(-3.5 * abs(slope[adj] + 0.05))


###################################################
### code chunk number 29: gdistance-21
###################################################
x <- geoCorrection(speed, scl=FALSE) 


###################################################
### code chunk number 30: gdistance-22
###################################################
A <- c(2667670,6479000)
B <- c(2667800,6479400)
AtoB <- shortestPath(x, A, B, output="SpatialLines")
BtoA <- shortestPath(x, B, A, output="SpatialLines")


###################################################
### code chunk number 31: fig8plot
###################################################
plot(r, main="")
lines(AtoB, col="red", lwd=2)
lines(BtoA, col="blue")
text(A[1]-10,A[2]-10,"A")
text(B[1]+10,B[2]+10,"B")


###################################################
### code chunk number 32: fig8
###################################################
plot(r, main="")
lines(AtoB, col="red", lwd=2)
lines(BtoA, col="blue")
text(A[1]-10,A[2]-10,"A")
text(B[1]+10,B[2]+10,"B")


###################################################
### code chunk number 33: gdistance-24
###################################################
Europe <- raster(system.file("external/Europe.grd", 
  package="gdistance"))
Europe[is.na(Europe)] <- 0
data(genDist)
data(popCoord)
pC <- as.matrix(popCoord[c("x","y")])


###################################################
### code chunk number 34: gdistance.Rnw:724-726
###################################################
plot(Europe, main="")
text(pC[,1],pC[,2],unlist(popCoord["Population"]),cex=.7)


###################################################
### code chunk number 35: gdistance-25
###################################################
geoDist <- pointDistance(pC, longlat=TRUE)
geoDist <- as.dist(geoDist)
Europe <- aggregate(Europe,3)
tr <- transition(Europe, mean, directions=8)
trC <- geoCorrection(tr, "c", scl=TRUE)
trR <- geoCorrection(tr, "r", scl=TRUE)
cosDist <- costDistance(trC,pC)
resDist <- commuteDistance(trR, pC)
cor(genDist,geoDist)
cor(genDist,cosDist)
cor(genDist,resDist)


###################################################
### code chunk number 36: gdistance-26
###################################################
origin <- unlist(popCoord[22,c("x","y")])
pI <- pathInc(trC, origin=origin, from=pC, 
  functions=list(overlap))
cor(genDist,pI[[1]])


