


# Fisher scoring
# starting value set  intercept to 
# p based on simple binomial and the rest to zero
logisticRegressionSimple<- function(X,y, 
                              baselineLogit=NULL,
                              betaStart=NULL, 
                              NIT=50, 
                              tol=1e-6,
                              verbose=FALSE){

m<- ncol(X)

if(is.null(betaStart)){
  pHat0<- mean(y)
  logit0<- exp(pHat0) / (1+exp(pHat0))
  betaOld<- c(logit0,rep(0,(m-1))) # this is only used for convergence test in first iteration
}
else{
  betaOld<- betaStart
}

if( verbose){
cat( "starting values", betaOld,fill=TRUE)
}

for( k in 1:NIT){
  
  # update pHat
  if( is.null(baselineLogit)){
    logit<- X%*%betaOld 
  }
  else{
    logit<- X%*%betaOld + baselineLogit
  }
  pHat<- exp( logit)/ ( 1+exp( logit) )
  
  W<- ( pHat*(1-pHat))
  # Fisher method of scoring -- variant of Newtons method
  #
  #S<- t(X)%*%c(y-pHat)
  #I<- t(X)%*% (diag(pHat*(1-pHat))%*%X)
  #betaNew<- betaOld + solve(I)%*%S
  # below is the identical computation but 
  # tricking the lm function to do it!
  # not obvious why this works
  z<- logit  + (1/W)*( y- pHat)
  fit<- lm( z~X-1, weights=W)
  betaNew<- fit$coefficients
  
  convTest<- mean( abs( betaNew - betaOld))
  if( verbose){
    cat( k, convTest,fill=TRUE)
  }
  if( convTest< tol){
    break
  }
  betaOld<- betaNew
}

# final fit 
lmfit<- lm( z~X-1, weights=W)
betaNew<- fit$coefficients
SE<- summary( fit)$coeffients[, 2]


return( list( beta = betaNew, 
      baseline = baselineLogit, 
      SE = SE,
      lmfit=lmfit,
       iter = k
        )
)
}