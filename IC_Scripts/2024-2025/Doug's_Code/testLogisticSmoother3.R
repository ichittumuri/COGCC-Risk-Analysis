

setwd("~/Dropbox/Home/Desktop2/PhDProjects/FlowIs")
library( fields)

source("logisticSmoother.R")
source("logisticSmootherBayes.R")

set.seed(226)
n<- 50
s<-  sort(runif( n))
s<- seq( 0,1,length.out=n)
p<-  s*(1-s)^3 +.001
pTrue<- p / (max(p)*1.05)
s<- cbind( s)
y<- rbinom( n,1, pTrue )

L<- 25
sigmaGrid<- 10^seq( log10(.5), log10(15), length.out=L)
logLikeGrid<- rep( NA, L)
pOld<- mean(y)
nuStart<-  rep( log(pOld/(1-pOld)), length(y))
for( j in L:1){
  fit<- logisticSmootherBayes(s, y, 
                              nuStart=nuStart, 
                              sigma=sigmaGrid[j]) 
                              
  logLikeGrid[j]<- fit$logLike
  nuStart<- fit$fitted.values
}

ind<- which.max(logLikeGrid)
print( sigmaGrid[ind])

fit<- logisticSmootherBayes(s, y, 
                            nuStart=nuStart, 
                            sigma=sigmaGrid[ind])
fit1<- logisticSmootherBayes(s, y, 
                            nuStart=nuStart, 
                            sigma=40) 
                            
plot( s, y, pch=16, col="grey")
lines(s,fit$pHat, 
      lwd=3, lty=2, col="magenta" ) 
lines(s,fit1$pHat, 
      lwd=3, lty=2, col="magenta" ) 

lines(s,pTrue,col="orange3"  )







