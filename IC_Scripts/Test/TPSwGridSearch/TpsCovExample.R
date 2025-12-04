library( fields)

epsilon <- .5
s0<- rbind( c(0,0),
            c(0,-epsilon),
            c(-epsilon,0),
            c(epsilon,0),
            c( 0,epsilon)
            )
# Think of s0 as a center point with its four cardinal neighbors at distance ε.

contrast<- cbind(c( 1,-.25,-.25,-.25,-.25))
# important that contrast zeroes out a linear surface. 
#
# A 5-point **contrast** on the cross of locations in s0.
# Order in s0 is: center, south, west, east, north.
# Weights: +1 on the center, −1/4 on each of the 4 neighbors.
# This computes:   center − mean(neighbors)
#
# Why these weights?
# • They **sum to zero**:        1 + 4*(−1/4) = 0  → removes any constant level.
# • They give **zero first moments** in x and y:
#     sum(w_i * x_i) = 0,  sum(w_i * y_i) = 0     → removes any linear tilt (a + bx + cy).
#   So for any plane f(x,y)=a+bx+cy, this contrast evaluates to 0
#   (the center equals the average of its four symmetric neighbors).
#
# What does it measure then?
# • Only the local **curvature/bumps** (the part TPS actually models/penalizes).
# • It’s the 2-D analog of a **second difference** (discrete Laplacian stencil).
#
# Why is that important for TPS/IRF-2?
# • TPS is only defined up to a plane; raw variances depend on where you pin that plane.
# • Using a plane-killing contrast gives **stationary increments**:
#   the contrast’s variance is (effectively) constant across locations and
#   its correlation depends mainly on distance.
#
# Make a weight vector that does: center − average(neighbors).
# This cancels any flat/tilted plane so it only measures curvature.

cardinalX<- rbind(c(0,0),
                  c(1,0),
                  c(0,1)
)
# Think of a thin-plate spline (TPS) surface like a flexible sheet with little bumps on it. 
# You can lift the whole sheet up/down (change its level).
# You can tilt it like a tray (change its overall slope).
# Those moves don’t change the bumps. TPS only “cares” about the bumps.
# Because of that, the math can’t decide the sheet’s exact level/tilt unless we pin it down.
#
# This picks three pin locations. A plane (level + tilt) is determined by three points, so:
# 1st point fixes the height reference.
# 2nd point fixes tilt left–right.
# 3rd point fixes tilt up–down.
#
# Pick 3 anchor points to "pin" the plane (level & tilt),
# so TPS covariances are well-defined.

M<- 100
delta<- seq( 0,3,length.out=M)
#shift cross by delta, 100 different times of magnitude between 0-3 

covTps<- rep( NA, M)
vars1<- rep( NA, M)
covTpsCenter<- rep( NA, M)
# create empty vectors of M size

tmp<- Tps.cov(s0,s0, cardinalX= cardinalX)
# Get the covariance between every pair of the 5 points in s0 (TPS model)
vars0 <- c(t(contrast)%*%tmp%*%contrast)
# Var of a weighted sum = weights × covariance table × weights.

for(k in 1:M ){
# Try the k-th distance between our two plus-shaped stamps.
  
  s1<- s0 + delta[k]
  # Make a second copy of the plus-shape and slide it diagonally by delta[k].
  # (Every point moves by (delta[k], delta[k]).)
  
   tmp<- Tps.cov(s0,s1, cardinalX= cardinalX)
   # Ask the “blanket model” how the readings from HOME (s0) and MOVED (s1)
   # tend to wiggle together (a 5x5 cross-covariance matrix).
   
   covTpsCenter[k]<- tmp[1,1] # raw covariance of two center points
   # Pick out the togetherness of the two CENTER dots:
   # home center vs moved center.
   
   covTps[k] <- t(contrast)%*%tmp%*%contrast
   # Now measure togetherness of the “squish test” (center - average of arms)
   # at home vs the same squish test at the moved cross.
   
   tmp<- Tps.cov(s1,s1, cardinalX= cardinalX)
   # Ask how the MOVED cross’s own readings wiggle with themselves
   # (a 5x5 covariance/variance matrix for s1).
   
   vars1[k] <- t(contrast)%*%tmp%*%contrast
   # How wiggly the MOVED “squish test” is by itself — its variance.
}

plot( delta, # how far we slid the second cross
      covTps/sqrt( vars0*vars1), #covariance between the two contrasts 
      # / the variance of that contrast at the original and shifted crosses
      # sqrt turns covariance into a correlation
      type="l",
      xlab="distance separation", ylab="correlation contrasts")
# correlation of “curvature readings” vs how far apart the two crosses are.
# Interpretation: in the thin-plate-spline (TPS) model, 
# the local curvature at one place tells you less and less about the curvature at another place as separation grows; 
# after about that distance, they’re essentially unrelated.

stats( vars1) # variance is constant 


# here is the plot for just correlation among the center points
S0<- rbind( c(0,0))
# home center dot at spot (0,0
S1<- matrix(s0, byrow=TRUE, nrow=M, ncol=2) +delta 
# all the moved center dots
COV01<- Tps.cov(S0,S1, cardinalX= cardinalX)
# How much do the heights of the home center and each moved center wiggle together?
VAR0<- c(Tps.cov(S0,S0, cardinalX= cardinalX))
# How much does the home center wiggle by itself? 
VAR1<-  diag( Tps.cov(S1,S1, cardinalX= cardinalX))
# And how much does each moved center wiggle by itself?

# funny variation due to interaction of delta spacing with
# cardinal points.

plot( delta, VAR1)
# X (delta): how far the second dot was moved away.
# Y (VAR1): how much the height at that moved dot can wiggle all by itself (its variance).
title("weird - not stationary")

plot( delta, COV01/ sqrt( VAR1*VAR0), type="l")
# X (delta): same distance as above.
# Y: correlation between the home center height and the moved dot height
title( "correlations")

