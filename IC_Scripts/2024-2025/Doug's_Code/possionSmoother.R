possionSmoother<- function(s, y, df){
  nuOld<- rep( log( mean(y)), length(y))
  
  for( k in 1:10){
    muI<- exp( nuOld) 
    W<- c(muI)
    z<- nuOld  + (1/muI)*( y- muI)
    # in place of WLS -- a smoothing/curve fitting step 
    # note that smooting found by default method ( CV)
    tempObj<- Tps( s,z,
                   weights=W,
                   df=df, give.warnings=FALSE)
    nuNew <- tempObj$fitted.values
    test<-  mean( abs(nuNew- nuOld) )
    cat( k, test, fill=TRUE)
    if( test< 1e-5){
      break
    }
    nuOld<- nuNew
  }
  return(tempObj)
}