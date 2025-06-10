logisticSmootherMonteCarlo<- function(s, y, sigma, M=50,
                                 nuStart=NULL,
                                 aRange=NULL,
                                 nIter=50,
                                 M =20)
                                {
  obj0<- logisticSmootherBayes(s, y, sigma=sigma, 
                               nuStart=nuStart, aRange=aRange)
  n<- length( y)
  X<- cbind( 1, s)
  #K<- sigma^2*stationary.cov(s,s, Covariance = Tps.cov,
  #                          lambda= 1/sigma^2)
  K<- sigma^2*stationary.cov(s,s, Covariance = Matern,
                             smoothness= 1.5, 
                             aRange= aRange)
  KC<- chol( K)
  logDetK<- 2*sum( log( diag(KC)))

  E<- matrix( rnorm(n*M ),n,M)
  gSim<- t( KC)%*%E
  # find log likelihoods with simulated g
  nuSim<- c(X%*%obj0$beta) + gSim
  tmp<- log( 1 + exp( nuSim))
  logL<- t(nuSim)%*%y - colSums( tmp) - .5* colSums( E^2) - logDetK
  
  
  
 

# compute log likelihood
# 

return(outObject)
}
