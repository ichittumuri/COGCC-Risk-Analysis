
# this is the code for the function:
#   # --- grid‐search over lambda on the TRAIN set to find MLE of λ ---
#   M <- 10
# lambdaGrid <- 10^seq(-2, 2, length.out = M)
# logLike    <- numeric(M)
# 
# for (k in M:1) {
#   # warm-start
#   if (k < M) {
#     nuOld <- look2$fitted.values
#   } else {
#     nuOld <- NULL
#   }
#   
#   # positional args: coord, response, then named lambda & nuOld
#   look2 <- logisticSmoother(
#     coord_model,
#     risk_model,
#     lambda = lambdaGrid[k],
#     nuOld  = nuOld
#   )
#   
#   logLike[k] <- look2$summary["lnProfileLike.FULL"]
#   cat("λ =", format(lambdaGrid[k], digits = 3),
#       "→ logLik =", round(logLike[k], 3), "\n")
# }
# 
# # plot coarse grid
# plot(
#   log10(lambdaGrid), logLike, type = "b",
#   xlab = "log10(λ)", ylab = "Profile log likelihood"
# )
# 
# # spline‐interpolate on a finer grid
# lGrid          <- seq(-2, 3, length.out = 250)
# profile_spline <- splint(log10(lambdaGrid), logLike, lGrid)
# lines(lGrid, profile_spline, col = "blue")
# 
# # find the maximizer
# logLambdaHat <- lGrid[which.max(profile_spline)]
# abline(v = logLambdaHat, col = "blue", lty = 2)
# cat("Estimated log10(λ) =", round(logLambdaHat, 3), "\n")
# 
# # re-fit at the best λ
# bestLambda <- 10^logLambdaHat
# cat("Ideal λ =", signif(bestLambda, 3), "\n")
# MLEFit     <- logisticSmoother(
#   coord_model,
#   risk_model,
#   lambda = bestLambda
# )