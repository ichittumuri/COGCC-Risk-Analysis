

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/IC_Scripts/2024-2025/fields_example")
library( fields)
source("possionSmoother.R")
source("logisticSmoother.R")

# synthetic test of the logistic smoother

n<- 1000
set.seed(222)
s<- sort(runif( n))
p<-  s*(1-s)*3.9
s<- cbind( s)
y0<- rbinom( n,1, p )
look0<- logisticSmoother(s, y0, lambda= 1e-3)
pHat<- exp(look0$fitted.values)/(1 + exp(look0$fitted.values))
plot( s, p, cex=.5, type="l")
lines(s,pHat, col="magenta" )
fit<- glm( y0~ s + I(s^2), family=binomial())
points( s, y0, pch=16, cex=.5)
lines(s, fit$fitted.values, col="blue" )


#look<- possionSmoother(loc, y, df=40)


load("tornadoCountsSpace.rda")
ind<- y>0
bubblePlot(loc[ind,], y[ind])

# tornado data as 0,1 
yTest<- ifelse( y ==0, 0, 1)
ind<- yTest==1
plot(loc[ind,])

look2<- logisticSmoother(loc, yTest, lambda=1e-1)

gHat2<- predictSurface( look2, nx=80, ny=80)

imagePlot( gHat2$x, gHat2$y, exp(gHat2$z )/(1+ exp(gHat2$z)),
           col=turbo(256))
points( loc[yTest==1,], col="white", pch=".")
US(add=TRUE, col="magenta")


# grid search over lambda and saving approx profile loglikelihood

M<- 10
# this grid found by trail and error
lambdaGrid<- 10**seq( -2,2,length.out=M)
logLike<- rep( NA, M)
for( k in M:1){
  if( k>1){
    nuOld<- look2$fitted.values
  }
  else{
    nuOld<-NULL
  }
  # use previous fit as starting values. 
  look2<- logisticSmoother(loc, yTest, lambda=lambdaGrid[k],
                           nuOld=nuOld)
  logLike[k]<- look2$summary["lnProfileLike.FULL"]
  cat(lambdaGrid[k], logLike[k], fill=TRUE)
}

plot( log10(lambdaGrid),logLike )


# a bit of a tricky way to fill in the grid search
# to get a better idea of maximum
lGrid<- seq( -2,3,length.out=250)
look<- splint(log10(lambdaGrid),logLike, lGrid )
lines( lGrid, look, col="blue")
logLambdaHat<- lGrid[which.max(look)]
xline( logLambdaHat, col="blue", lty=2)

# spline fit with MLE for lambda 
MLEFit<- logisticSmoother(loc, yTest, lambda=10**logLambdaHat)
#not much features to this -- maybe not the best data set for illustration!

surface( MLEFit)
US(add=TRUE, col="magenta")
points( loc[yTest==1,], pch=16, cex=.2)
title("logit surface ")

gHat<- predictSurface( MLEFit)

imagePlot( gHat$x, gHat$y, exp( gHat$z)/ ( 1+exp( gHat$z) ))
US(add=TRUE, col="magenta")
points( loc[yTest==1,], pch=16, cex=.2)
title("logit surface ")



