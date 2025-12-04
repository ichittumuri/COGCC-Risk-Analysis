# =============================================================================
# logisticSmoother: IRLS + spatial smoother + optional Z and weights
# =============================================================================
library(fields)

logisticSmoother <- function(
    s,
    y,
    lambda,
    Z       = NULL,
    betaHat = NULL,
    nuOld   = NULL,
    weights = NULL
) {
  # IRLS + spatial smoothing on coords s + optional covariates Z
  # 'weights' are observation weights v_i (e.g. distance-based)
  
  n <- length(y)
  
  if (is.null(weights)) {
    weights <- rep(1, n)
  }
  if (length(weights) != n) {
    stop("Length of 'weights' must match length(y).")
  }
  if (any(!is.finite(weights))) {
    stop("Non-finite values in 'weights'.")
  }
  
  # ---------------------------------------------------------------------------
  # Initial linear predictor
  # ---------------------------------------------------------------------------
  if (is.null(nuOld)) {
    if (!is.null(Z) && !is.null(betaHat)) {
      nuOld <- as.vector(Z %*% betaHat)
    } else {
      pStart  <- mean(y)
      pStart  <- min(max(pStart, 1e-6), 1 - 1e-6)  # avoid 0 or 1
      nuStart <- log(pStart / (1 - pStart))
      nuOld   <- rep(nuStart, n)
    }
  }
  
  eps_p <- 1e-6      # keep probabilities away from 0/1
  eps_w <- 1e-6      # avoid division by zero in W_core
  
  # ---------------------------------------------------------------------------
  # IRLS iterations
  # ---------------------------------------------------------------------------
  for (k in 1:20) {
    # Clip nuOld to avoid overflow in exp()
    nuOld <- pmin(pmax(nuOld, -20), 20)
    
    # Stable logistic
    pOld <- 1 / (1 + exp(-nuOld))
    pOld <- pmin(pmax(pOld, eps_p), 1 - eps_p)
    
    # IRLS weights
    W_core <- pOld * (1 - pOld)
    W_core <- pmax(W_core, eps_w)
    W      <- W_core * weights
    
    # Working response
    z <- nuOld + (y - pOld) / W_core
    
    if (!is.null(Z)) {
      z_spatial <- z - as.vector(Z %*% betaHat)
    } else {
      z_spatial <- z
    }
    
    # Spatial smoothing step
    tempObj <- spatialProcess(
      s,
      z_spatial,
      cov.function = "Tps.cov",
      weights      = W,
      lambda       = lambda
    )
    
    f_hat <- tempObj$fitted.values
    
    if (any(!is.finite(f_hat))) {
      stop("Non-finite fitted values from spatialProcess; check weights / lambda.")
    }
    
    # Updated linear predictor
    if (!is.null(Z)) {
      nuNew <- as.vector(Z %*% betaHat) + f_hat
    } else {
      nuNew <- f_hat
    }
    
    if (any(!is.finite(nuNew))) {
      stop("Non-finite values in nuNew (probability overflow/underflow).")
    }
    
    testConv <- mean(abs(nuNew - nuOld))
    if (!is.finite(testConv)) {
      stop("Non-finite testConv; IRLS divergence.")
    }
    if (testConv < 1e-5) break
    
    nuOld <- nuNew
  }
  
  cat("[IRLS Iterations]:", k, "\n")
  return(tempObj)
}
