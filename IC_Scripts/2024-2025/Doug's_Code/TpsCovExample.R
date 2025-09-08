library( fields)
epsilon <- .5

s0<- rbind( c(0,0),
            c(0,-epsilon),
            c(-epsilon,0),
            c(epsilon,0),
            c( 0,epsilon)
            )

contrast<- cbind(c( 1,-.25,-.25,-.25,-.25))
# important that contrast zeroes out a linear surface. 

cardinalX<- rbind(c(0,0),
                  c(1,0),
                  c(0,1)
)

M<- 100
delta<- seq( 0,3,length.out=M)
covTps<- rep( NA, M)

vars1<- rep( NA, M)
covTpsCenter<- rep( NA, M)

tmp<- Tps.cov(s0,s0, cardinalX= cardinalX)
vars0 <- c(t(contrast)%*%tmp%*%contrast)

for(k in 1:M ){
  s1<- s0 + delta[k]
   tmp<- Tps.cov(s0,s1, cardinalX= cardinalX)
   
   covTpsCenter[k]<- tmp[1,1] # raw covariance of two center points
   covTps[k] <- t(contrast)%*%tmp%*%contrast
   
   
   tmp<- Tps.cov(s1,s1, cardinalX= cardinalX)
   vars1[k] <- t(contrast)%*%tmp%*%contrast
}

plot( delta, covTps/sqrt( vars0*vars1), type="l",
      xlab="distance separation", ylab="correlation contrasts")
# variance is constant 
stats( vars1)


# here is the plot for just correlation among the center points
S0<- rbind( c(0,0))
S1<- matrix(s0, byrow=TRUE, nrow=M, ncol=2) +delta 
COV01<- Tps.cov(S0,S1, cardinalX= cardinalX)
VAR0<- c(Tps.cov(S0,S0, cardinalX= cardinalX))
VAR1<-  diag( Tps.cov(S1,S1, cardinalX= cardinalX))

# funny variation due to interaction of delta spacing with
# cardinal points.

plot( delta, VAR1)
title("weird - not stationary")

plot( delta, COV01/ sqrt( VAR1*VAR0), type="l")
title( "correlations")