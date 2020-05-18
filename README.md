# expdepth
R functions for computation of the bivariate expectile depth and regions and plot the BExPlot

Code created by I. Cascos and M. Ochoa. The theoretical description of the procedures is avilable at
https://e-archivo.uc3m.es/handle/10016/28434

Copy and paste the piece of code below, which uses the `ddalpha` package in order to plot some BExPlots, contours of bivariate expectile regions and a D-D chart for the expectile depth.

```{r}
source("https://raw.githubusercontent.com/icascos/expdepth/master/BExPlot.R")
source("https://raw.githubusercontent.com/icascos/expdepth/master/exactexp.R")
source("https://raw.githubusercontent.com/icascos/expdepth/master/expdepth.R")
require(ddalpha)
data(hemophilia)
attach(hemophilia)
AHF.normal <-cbind(AHFactivity[gr=="normal"],AHFactivity.1[gr=="normal"])

par(mfrow=c(2,2))
BExPlot(AHF.normal)
BExPlot(AHF.carrier)

AHF.carrier <-cbind(AHFactivity[gr=="carrier"],AHFactivity.1[gr=="carrier"])
plot(AHF.carrier,pch=3)
points(x=mean(AHF.carrier[,1]),y=mean(AHF.carrier[,2]),pch=16)
lines(exactexp(AHF.carrier))
lines(exactexp(AHF.carrier,alpha=0.01))
lines(exactexp(AHF.carrier,alpha=0.05))
lines(exactexp(AHF.carrier,alpha=0.4))

AHF <-rbind(AHF.normal,AHF.carrier)
e1 <- vector(length=nrow(AHF))
for(i in 1:nrow(AHF)) e1[i]<-expdepth(x=AHF[i,],data=AHF.normal)
e2 <- vector(length=nrow(AHF))
for(i in 1:nrow(AHF)) e2[i]<-expdepth(x=AHF[i,],data=AHF.carrier)
plot(x=e1,y=e2,main="D-D plot AHF activy and AHF antigen",
     xlab="Non-carrier depths",ylab="Carrier depths")
abline(0,1,lty=2)
```
