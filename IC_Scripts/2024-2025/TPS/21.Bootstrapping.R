# ============================================================
# 1) TRAIN fitted probabilities
#    used to generate synthetic bootstrap responses
# ============================================================

f_hat_train <- as.numeric(predict(MLEFit, train_coord_model))
eta_glm_train <- as.vector(Z_train %*% betaHat)
nu_train <- eta_glm_train + f_hat_train
pHat_train <- plogis(nu_train)

# ============================================================
# 2) TEST predictions from original fitted model
#    used to compute the original test log-likelihood
# ============================================================

f_hat_test <- as.numeric(predict(MLEFit, test_coord_model))
eta_glm_test <- as.vector(Z_test %*% betaHat)
nu_test <- eta_glm_test + f_hat_test
prob_test <- plogis(nu_test)

eps <- 1e-15
p_hat <- pmin(pmax(prob_test, eps), 1 - eps)
risk_model <- test_risk_model

loglik_i <- risk_model * log(p_hat) +
  (1 - risk_model) * log(1 - p_hat)

loglik_obs <- sum(loglik_i)

cat("Original test log-likelihood =", round(loglik_obs, 2), "\n")

# ============================================================
# 3) Enter bootstrap land
# ============================================================

nBoot <- 1000
loglikBoot <- rep(NA, nBoot)

set.seed(233)

for (k in 1:nBoot) {
  
  # synthetic TRAIN response
  YS <- rbinom(n = length(pHat_train), size = 1, prob = pHat_train)
  
  # refit weighted GLM
  Z_df_train_boot <- as.data.frame(Z_train)
  Z_df_train_boot$risk <- YS
  
  fitS <- glm(
    risk ~ . -1,
    data    = Z_df_train_boot,
    family  = binomial(link = "logit"),
    weights = w_train
  )
  
  betaTmp <- coef(fitS)
  
  # refit weighted spatial smoother
  MLEFitS <- logisticSmoother(
    s       = train_coord_model,
    y       = YS,
    lambda  = bestLambda,
    Z       = Z_train,
    betaHat = betaTmp,
    weights = w_train
  )
  
  # predict on TEST set
  f_hat_testS <- as.numeric(predict(MLEFitS, test_coord_model))
  eta_glm_testS <- as.vector(Z_test %*% betaTmp)
  nu_testS <- eta_glm_testS + f_hat_testS
  prob_testS <- plogis(nu_testS)
  
  # bootstrap test log-likelihood
  pS <- pmin(pmax(prob_testS, eps), 1 - eps)
  
  loglik_iS <- risk_model * log(pS) +
    (1 - risk_model) * log(1 - pS)
  
  loglikBoot[k] <- sum(loglik_iS)
  
  if (k %% 50 == 0) {
    cat("finished bootstrap", k, "of", nBoot, "\n")
  }
}

# ============================================================
# 4) Compare bootstrap distribution to original
# ============================================================

summary(loglikBoot)
quantile(loglikBoot, c(0.025, 0.975), na.rm = TRUE)

# --- Compute summaries ---
ci <- quantile(loglikBoot, c(0.025, 0.975), na.rm = TRUE)
boot_mean <- mean(loglikBoot, na.rm = TRUE)

cat("Original test log-likelihood =", round(loglik_obs, 2), "\n")
cat("Bootstrap mean log-likelihood =", round(boot_mean, 2), "\n")
cat("Bootstrap sd =", round(sd(loglikBoot, na.rm = TRUE), 2), "\n")
cat("95% CI =", round(ci, 2), "\n")

# --- Histogram ---
hist(loglikBoot,
     breaks = 30,
     xlim = quantile(loglikBoot, c(0.01, 0.99)),
     main = "Bootstrap test log-likelihood",
     xlab = "Test log-likelihood")

# --- Add vertical lines ---
abline(v = loglik_obs, col = "red", lwd = 2, lty = 2)   # observed
abline(v = boot_mean, col = "blue", lwd = 2, lty = 2)  # mean
abline(v = ci, col = "black", lwd = 2, lty = 2)          # 95% CI

# --- Legend ---
legend("topright",
       legend = c("Observed", "Bootstrap mean", "95% CI"),
       col = c("red", "blue", "black"),
       lwd = 2,
       lty = c(2, 2, 2))

