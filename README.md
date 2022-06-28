# expdepth
R functions for computation of the bivariate expectile depth and regions and plot the BExPlot

Code created by I. Cascos and M. Ochoa. The theoretical description of the procedures is avilable at
https://e-archivo.uc3m.es/handle/10016/28434

Copy and paste the piece of code below, which uses the `ddalpha` package in order to plot some BExPlots, contours of bivariate expectile regions and a D-D chart for the expectile depth.

```{r}
source("https://raw.githubusercontent.com/icascos/expdepth/master/BExPlot.R")
source("https://raw.githubusercontent.com/icascos/expdepth/master/exactBExPlot.R")
source("https://raw.githubusercontent.com/icascos/expdepth/master/exactexp.R")
source("https://raw.githubusercontent.com/icascos/expdepth/master/expdepth.R")
source("https://raw.githubusercontent.com/icascos/expdepth/master/distexpdepth.R")
require(ddalpha)
data(hemophilia)
attach(hemophilia)
AHF.normal <-cbind(AHFactivity[gr=="normal"],AHFactivity.1[gr=="normal"])
AHF.carrier <-cbind(AHFactivity[gr=="carrier"],AHFactivity.1[gr=="carrier"])

# approximate BExPlots
par(mfrow=c(3,2))
BExPlot(AHF.normal)
BExPlot(AHF.carrier)

# exact BExPlot
exactBExPlot(AHF.normal)

# exact expectile regions
plot(AHF.carrier,pch=3)
points(x=mean(AHF.carrier[,1]),y=mean(AHF.carrier[,2]),pch=16)
lines(exactexp(AHF.carrier))
lines(exactexp(AHF.carrier,alpha=0.01))
lines(exactexp(AHF.carrier,alpha=0.05))
lines(exactexp(AHF.carrier,alpha=0.4))

# DD-plot for the expectile depth
X <- matrix(rnorm(100),ncol=2)
X <- rbind(X,matrix(rnorm(4,mean=5),ncol=2))
Y <- matrix(rnorm(100),ncol=2) 
Y <- rbind(Y,matrix(rnorm(4,mean=-5),ncol=2))
sim.data <-rbind(X,Y)
e1 <- vector(length=nrow(sim.data))
for(i in 1:length(e1)) e1[i]<-expdepth(x=sim.data[i,],data=X)
e2 <- vector(length=nrow(sim.data))
for(i in 1:length(e2)) e2[i]<-expdepth(x=sim.data[i,],data=Y)
plot(x=e1,y=e2,main="D-D plot, expectile depth",
     xlab="X",ylab="Y")
abline(0,1,lty=2)

# DD-plot for the distorted expectile depth (sigmoid distortion)
ed1 <- vector(length=nrow(sim.data))
for(i in 1:length(ed1)) ed1[i]<-distexpdepth(x=sim.data[i,],data=X,gtilde=sigmoid,sim=TRUE)
ed2 <- vector(length=nrow(sim.data))
for(i in 1:length(ed2)) ed2[i]<-distexpdepth(x=sim.data[i,],data=Y,gtilde=sigmoid,sim=TRUE)
plot(x=ed1,y=ed2,main="D-D plot, distorted expectile depth, sigmoid",
     xlab="X",ylab="Y")
abline(0,1,lty=2)
```
