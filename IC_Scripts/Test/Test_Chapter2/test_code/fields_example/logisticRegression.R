


# Fisher scoring
# starting value set  intercept to 
# p based on simple binomial and the rest to zero
logisticRegression<- function(X,y, 
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
  S<- t(X)%*%c(y-pHat)
  
  # Fisher method of scoring -- variant of Newtons method
  #
  #I<- t(X)%*% (diag(pHat*(1-pHat))%*%X)
  #betaNew<- betaOld + solve(I)%*%S
  
  # more numericaly stable  method that avoid  computing X^T%*%W%*% X
  #
  qrX<- qr( X*c(sqrt(pHat*(1-pHat))))
  R<- qr.R(qrX)
  # Here I = t(R)%*%R  R is upper triangular
  # evaluate I^{-1}S by two triangular solves  
  #  t(R)S = v and then  R u = v. u is then I^{-1}S
  FisherStep<- forwardsolve(R,transpose=TRUE,
                            upper.tri=TRUE, x=S)
  FisherStep<- backsolve( R, FisherStep)
  betaNew<- betaOld + FisherStep                       
  
 
  convTest<- mean( abs( betaNew - betaOld))
  if( verbose){
    cat( k, convTest,fill=TRUE)
  }
  if( convTest< tol){
    break
  }
  betaOld<- betaNew
}


# invert  for final information matrix
qrX<- qr( X*c(sqrt(pHat*(1-pHat))))
R<- qr.R(qrX)
I<- t(R)%*% R
idMatrix<- diag( 1, ncol(R))
tmp<- forwardsolve(R,transpose=TRUE,
                          upper.tri=TRUE, x=idMatrix)
IInv<- backsolve( R, tmp)
SE<- sqrt(diag( IInv))

return( list( beta = betaNew, 
          baseline = baselineLogit, 
                SE = SE,
       iter = k,
       I= I,
       pHat= pHat
        )
)
}