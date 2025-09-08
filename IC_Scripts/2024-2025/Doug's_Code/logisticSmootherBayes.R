logisticSmootherBayes<- function(s, y, sigma, M=50,
                                 nuStart=NULL,
                                 aRange=NULL,
                                 nIter=50)
                                {
  # "nu" is in the logit/linear model space
  
  if( is.null(aRange)){
    aRange<- 1.5* (max(s)- min(s))
  }
  if( is.null(nuStart)){
  pOld<- mean(y)
  nuStart<- log( pOld/ (1-pOld))
  nuOld<- rep( nuStart, length(y))
  }
  else{
    nuOld<- nuStart
  }
  
  X<- cbind( 1, s)
  K<- sigma^2*stationary.cov(s,s, Covariance = Matern,
                             smoothness= 1.5, 
                            aRange= aRange)
  
  #K<- (sigma^2)*Tps.cov( s,s, cardinalX = c(.25, .75), m=2)
  colTable<- turbo(nIter)
  for( k in 1:nIter){
    pOld<- exp( nuOld)/ ( 1+ exp(nuOld))
    #lines(s, pOld, col= colTable[k])
    W<- c(pOld*(1-pOld))
    z<- nuOld  + (1/W)*( y- pOld)
    # in place of WLS -- a smoothing/curve fitting step 
    M<- diag(1/W) + K
    Minv<- solve(M)
    betaHat<- solve( t(X)%*%Minv%*%X, t(X)%*%Minv%*% z)
    #print( betaHat)
    predFixed<- X%*%betaHat
    tempObj<- K%*% ( Minv%*%(z-predFixed))
                    
    nuNew <- tempObj + predFixed
    
    testConv<-  mean( abs(nuNew- nuOld) )
    #cat( k, test, fill=TRUE)
    if( testConv< 1e-5){
      break
    }
    nuOld<- nuNew 
  }
  pHat<-exp(nuNew)/(1+ exp(nuNew))
  W<- c(pHat*(1-pHat))
  z<- nuNew  + (1/W)*( y- pHat)
  # in place of WLS -- a smoothing/curve fitting step 
  M<- diag(1/W) + K
  Minv<- solve( M)
  res<- z- X%*%betaHat
  SS<- sum(res*(Minv%*%res)) # (y - XbetaHat)^T Minv (y - XbetaHat)
  H<- chol(M)
  logDet<- 2*sum(log( diag(H))) # log |Minv|
  logLike<-  -.5*SS - .5*logDet # + on det because inverse 
  
  outObject<- list( 
      fitted.values=nuNew, 
      pHat= pHat,
      beta= betaHat, 
      gHat= tempObj, 
      logLike=logLike, 
      iter=k,
      aRange=aRange,
      sigma=sigma)

# compute log likelihood
# 

return(outObject)
}
