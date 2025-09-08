setwd("~/Dropbox/Home/Projects/Hemlock")
library( fields)
load("tornadoCountsSpace.rda")
source("possionSmoother.R")
source("logisticSmoother.R")

# synthetic test of the logistic smoother

n<- 1000
set.seed(222)
s<- sort(runif( n))
p<-  s*(1-s)*3.9

y0<- rbinom( n,1, p )
look0<- logisticSmoother(s, y0, df=6)
pHat<- exp(look0$fitted.values)/(1 + exp(look0$fitted.values))
plot( s, pHat, cex=.5)
lines(s,p, col="magenta" )
fit<- glm( y0~ s + I(s^2), family=binomial())

lines(s, fit$fitted.values, col="blue" )


look<- possionSmoother(loc, y, df=40)

# tornado data as 0,1 
yTest<- ifelse( y ==0, 0, 1)
look2<- logisticSmoother(loc, yTest, df=40)

set.panel( 1,2)
gHat1<- predictSurface( look, nx=80, ny=80)
gHat2<- predictSurface( look2, nx=80, ny=80)
imagePlot( gHat1$x, gHat1$y, exp(gHat1$z ),
           col=turbo(256))
imagePlot( gHat2$x, gHat2$y, exp(gHat2$z )/(1+ exp(gHat2$z)),
           col=turbo(256))
#points( loc[yTest==1,], col="white")


