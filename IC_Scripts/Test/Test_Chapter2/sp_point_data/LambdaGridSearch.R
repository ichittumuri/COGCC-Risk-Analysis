LambdaGridSearch <- function(coord_model,
                             risk_model,
                             Z = NULL,
                             M           = 10,
                             exp_range   = c(-2,  2),
                             fine_length = 250,
                             init_nuOld  = NULL,
                             warm_start  = TRUE) {

  exps       <- seq(exp_range[1], exp_range[2], length.out = M)
  lambdaGrid <- 10^exps
  logLike    <- numeric(M)
  look2      <- NULL
  
  for (k in M:1) {
    if (k < M && warm_start) {
      nuOld <- look2$fitted.values
    } else {
      nuOld <- init_nuOld
    }
    
    look2 <- logisticSmoother(
      coord_model,
      risk_model,
      lambda = lambdaGrid[k],
      nuOld  = nuOld,
      Z = Z
    )
    
    logLike[k] <- look2$summary["lnProfileLike.FULL"]
    cat("λ =", format(lambdaGrid[k], digits = 3),
        "→ logLik =", round(logLike[k], 3), "\n")
  }
  
  lGrid          <- seq(-2, 3, length.out = fine_length)
  profile_spline <- fields::splint(log10(lambdaGrid), logLike, lGrid)
  
  logLambdaHat <- lGrid[which.max(profile_spline)]
  bestLambda   <- 10^logLambdaHat
  cat("Estimated log10(λ) =", round(logLambdaHat, 3), "\n")
  cat("Ideal λ =", signif(bestLambda, 3), "\n")
  
  invisible(list(
    exps           = exps,
    lambdaGrid     = lambdaGrid,
    logLike        = logLike,
    lGrid          = lGrid,
    profile_spline = profile_spline,
    logLambdaHat   = logLambdaHat,
    bestLambda     = bestLambda
  ))
}

PlotLambdaProfile <- function(profile) {
  plot(
    log10(profile$lambdaGrid), profile$logLike, type = "b",
    xlab = "log10(λ)", ylab = "Profile log likelihood"
  )
  lines(profile$lGrid, profile$profile_spline, col = "blue")
  abline(v = profile$logLambdaHat, col = "blue", lty = 2)
}
