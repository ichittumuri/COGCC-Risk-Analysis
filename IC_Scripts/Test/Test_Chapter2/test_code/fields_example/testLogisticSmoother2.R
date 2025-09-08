

setwd("~/Dropbox/Home/Desktop2/PhDProjects/FlowIs")
library( fields)

source("logisticSmoother.R")
source("logisticSmootherBayes.R")

# synthetic test of the logistic smoother

n<- 500
set.seed(226)
s<-  sort(runif( n))
s<- seq( 0,1,length.out=n)
p<-  s*(1-s)^3 +.001
p<- p / (max(p)*1.05)
s<- cbind( s)
y<- rbinom( n,1, p )
plot( s, y, cex=.5, type="p", pch=16)

look0<- logisticSmootherBayes(s, y, sigma=sqrt(1/4e-3), aRange=2)
look<- logisticSmoother(s, y,lambda=5e-4)
pHat0<- exp(look0$fitted.values)/(1 + exp(look0$fitted.values))
pHat<- exp(look$fitted.values)/(1 + exp(look$fitted.values))
lines(s,pHat0, col="magenta", lwd=3)
lines(s,pHat, col="grey", lwd=3, lty=2 ) 

lines( s, p, lty=2)
print(look0$logLike)
#lines(s,pHat, col="grey", lwd=3 )


set.seed(226)
n<- 50
s<-  sort(runif( n))
s<- seq( 0,1,length.out=n)
p<-  s*(1-s)^3 +.001
p<- p / (max(p)*1.05)
s<- cbind( s)
y<- rbinom( n,1, p )

L<- 25
sigmaGrid<- 10^seq( log10(.5), log10(15), length.out=L)
logLikeGrid<- rep( NA, L)
pOld<- mean(y)
nuStart<-  rep( log(pOld/(1-pOld)), length(y))
for( j in L:1){
  fit<- logisticSmootherBayes(s, y, 
                              nuStart=nuStart, 
                              sigma=sigmaGrid[j], 
                              aRange=2)
  logLikeGrid[j]<- fit$logLike
  nuStart<- fit$fitted.values
}

ind<- which.max(logLikeGrid)
print( sigmaGrid[ind])

fit<- logisticSmootherBayes(s, y, 
                            nuStart=nuStart, 
                            sigma=sigmaGrid[ind], 
                            aRange=2)
plot( s, y)
lines(s,fit$pHat, 
      col="grey", 
      lwd=3, lty=2 ) 







