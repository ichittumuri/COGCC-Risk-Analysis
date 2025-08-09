logisticSmoother<- function(s, y, Z=NULL, lambda, nuOld=NULL){
  # assume y are 0,1 s
  if( is.null(nuOld)){
  pStart<- mean(y)
  
  pOld<- pStart
  nuStart<- log( pOld/ (1-pOld))  #modfiy so it doesn't reach infty 
  # nu is in the logit/linear model spacex
  nuOld<- rep( nuStart, length(y))
  }
  
  for( k in 1:20){
    pOld<- exp( nuOld)/ ( 1+ exp(nuOld)) 
    # Clamp pOld to [1e-6, 1-1e-6] so it never becomes exactly 0 or 1 (avoids zero weights)
    pOld <- pmin(pmax(pOld, 5e-4), 1 - 5e-4) #
    
    W<- c(pOld*(1-pOld))
    z<- nuOld  + (1/W)*( y- pOld)
    # in place of WLS -- a smoothing/curve fitting step 
    # note that smoothing found by default method ( CV)
    tempObj<- spatialProcess( s,z,
                              cov.function ="Tps.cov",
                   weights=W,
                   lambda=lambda,
                   Z=Z)
                    
    nuNew <- tempObj$fitted.values
    testConv<-  mean( abs(nuNew- nuOld) )
    #cat( k, testConv, fill=TRUE)
    if( testConv< 1e-5){
      break
    }
    nuOld<- nuNew
  }
  print( k )
  return(tempObj)
}