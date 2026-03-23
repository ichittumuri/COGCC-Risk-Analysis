# ============================================================
# STEP 1: Setup
# ============================================================

# spatial part at TRAIN locations
f_hat_train <- as.numeric(predict(MLEFit, train_coord_model))

# fixed-effects part at TRAIN locations
eta_glm_train <- as.vector(Z_train %*% betaHat)

# full linear predictor on train
nu_train <- eta_glm_train + f_hat_train

# fitted probabilities on train
pHat_train <- plogis(nu_train)

loglik_obs <- loglik_weighted_spatialZ

# make one synthetic dataset
YS <- rbinom(n = length(pHat_train), size = 1, prob = pHat_train)

# replace the training response with YS.
Z_df_train_boot <- as.data.frame(Z_train)
Z_df_train_boot$risk <- YS

# refit the weighted GLM:
fitS <- glm(
  risk ~ . -1,
  data    = Z_df_train_boot,
  family  = binomial(link = "logit"),
  weights = w_train
)

betaTmp <- coef(fitS)

# refit the spatial smoother:
MLEFitS <- logisticSmoother(
  s       = train_coord_model,
  y       = YS,
  lambda  = bestLambda,
  Z       = Z_train,
  betaHat = betaTmp,
  weights = w_train
)

f_hat_testS <- as.numeric(predict(MLEFitS, test_coord_model))
eta_glm_testS <- as.vector(Z_test %*% betaTmp)
nu_testS <- eta_glm_testS + f_hat_testS
prob_testS <- plogis(nu_testS)


eps <- 1e-15
pS <- pmin(pmax(prob_testS, eps), 1 - eps)

loglik_iS <- test_risk_model * log(pS) +
  (1 - test_risk_model) * log(1 - pS)

loglikBoot_k <- sum(loglik_iS)

# compare different loglikihoods
print(loglikBoot_k)
print(loglik_obs)

















