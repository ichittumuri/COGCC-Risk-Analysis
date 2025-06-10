makeKrigingWeights<- function(obj, sStar){
  n<- length(obj$y)
  IM<- diag( 1, n)
  
  nPred<- nrow( sStar)
  W<- matrix(NA,n,nPred)
  for( k in 1:n){
    W[k,]<- predict( obj,sStar, ynew=IM[,k])
  }
 # test.for.zero( t(W)%*%obj$y, obj$fitted.values[1:10])
  return(W)
}