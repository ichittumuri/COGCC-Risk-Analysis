# =============================================================================
# SigmaGridSearch_spatialZ: grid search over σ for distance weights
# =============================================================================

SigmaGridSearch_spatialZ <- function(
    sigma_grid,
    train_data,
    test_data,
    train_coord_model,
    test_coord_model,
    Z_train,
    Z_test,
    train_risk_model
) {
  # Helper: log-likelihood for Bernoulli
  val_loglik <- function(y, p) {
    eps <- 1e-15
    p <- pmin(pmax(p, eps), 1 - eps)
    sum(y * log(p) + (1 - y) * log(1 - p))
  }
  
  loglik_results <- numeric(length(sigma_grid))
  
  for (i in seq_along(sigma_grid)) {
    sigma_dist <- sigma_grid[i]
    cat("\n--- Sigma =", sigma_dist, "---\n")
    
    d_tr <- train_data$match_distance_m
    
    # Distance-based weights for training
    w_train <- ifelse(
      train_data$risk == 1 & !is.na(d_tr),
      exp(-d_tr / sigma_dist),
      1
    )
    
    # GLM on Z with these weights
    Z_df_train <- as.data.frame(Z_train)
    Z_df_train$risk <- train_data$risk
    
    glm_sel <- try(
      glm(
        risk ~ . - 1,
        data    = Z_df_train,
        family  = binomial(link = "logit"),
        weights = w_train
      ),
      silent = TRUE
    )
    
    if (inherits(glm_sel, "try-error")) {
      warning("GLM failed for sigma = ", sigma_dist, "; setting loglik = NA")
      loglik_results[i] <- NA
      next
    }
    
    betaHat <- coef(glm_sel)
    
    # Lambda grid search for this sigma
    profile <- try(
      LambdaGridSearch(
        coord_model = train_coord_model,
        risk_model  = train_risk_model,
        Z           = Z_train,
        betaHat     = betaHat,
        weights     = w_train,
        doPlot      = FALSE
      ),
      silent = TRUE
    )
    
    if (inherits(profile, "try-error")) {
      warning("LambdaGridSearch failed for sigma = ", sigma_dist, "; setting loglik = NA")
      loglik_results[i] <- NA
      next
    }
    
    bestLambda <- profile$bestLambda
    
    # Fit spatial smoother with the chosen lambda
    MLEFit <- try(
      logisticSmoother(
        s       = train_coord_model,
        y       = train_risk_model,
        lambda  = bestLambda,
        Z       = Z_train,
        betaHat = betaHat,
        weights = w_train
      ),
      silent = TRUE
    )
    
    if (inherits(MLEFit, "try-error")) {
      warning("logisticSmoother failed for sigma = ", sigma_dist, "; setting loglik = NA")
      loglik_results[i] <- NA
      next
    }
    
    # Predict spatial part on TEST coords
    f_hat_test   <- as.numeric(predict(MLEFit, test_coord_model))
    
    # GLM linear predictor at TEST
    eta_glm_test <- as.vector(Z_test %*% betaHat)
    
    # Combined linear predictor and probabilities
    nu_all <- eta_glm_test + f_hat_test
    p_all  <- plogis(nu_all)
    
    # Test log-likelihood
    loglik_results[i] <- val_loglik(test_data$risk, p_all)
    cat("  → Test log-likelihood =", round(loglik_results[i], 2), "\n")
  }
  
  results_df <- data.frame(
    sigma  = sigma_grid,
    loglik = loglik_results
  )
  
  # Handle case where some are NA
  if (all(is.na(loglik_results))) {
    warning("All sigmas failed; returning NA.")
    best_sigma <- NA
  } else {
    best_idx   <- which.max(loglik_results)
    best_sigma <- sigma_grid[best_idx]
    
    cat("\n====================================\n")
    cat("BEST σ =", best_sigma,
        "with log-likelihood =", round(loglik_results[best_idx], 2), "\n")
    cat("====================================\n")
  }
  
  # Optional quick plot
  if (!all(is.na(loglik_results))) {
    plot(
      results_df$sigma, results_df$loglik, type = "b",
      xlab = "sigma (distance scale, m)",
      ylab = "Test log-likelihood",
      main = "Sigma grid search (spatial + Z + distance weights)"
    )
  }
  
  invisible(list(
    results_df = results_df,
    best_sigma = best_sigma
  ))
}
