latkrig_logisticSmoother<- function(s, y, lambda, nuOld=NULL){
  # assume y are 0,1 s
  if( is.null(nuOld)){
    pStart<- mean(y)
    
    pOld<- pStart
    nuStart<- log( pOld/ (1-pOld))
    # nu is in the logit/linear model space
    nuOld<- rep( nuStart, length(y))
  }
  
  for( k in 1:20){
    cat("Starting iteration", k, "\n")
    eps  <- 1e-6
    pOld <- pmin(pmax(pOld, eps), 1 - eps) #clamp probabilities away from 0 and 1
    W    <- pOld * (1 - pOld)
    z<- nuOld  + (1/W)*( y- pOld)
    # in place of WLS -- a smoothing/curve fitting step 
    # note that smoothing found by default method ( CV)
    # fit thin‐plate‐spline via LatticeKrig
    tempObj <- LatticeKrig(
      x            = s,
      y            = z,
      cov.function = "Tps.cov",
      weights      = W,
      lambda       = lambda
    )
    
    # predict nuNew at the original sites
    nuNew <- predict(tempObj, xnew = s)
    
    testConv<-  mean( abs(nuNew- nuOld) )
    #cat( k, test, fill=TRUE)
    if( testConv< 1e-5){
      break
    }
    nuOld<- nuNew
  }
  print( k )
  return(tempObj)
}