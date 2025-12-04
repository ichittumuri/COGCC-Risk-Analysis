# =============================================================================
# 1) Setup & Libraries
# =============================================================================
library(sf)
library(dplyr)
library(ggplot2)
library(fields)
library(caret)
library(readr)
library(tibble)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# =============================================================================
# 2) Sources & Data
# =============================================================================
source("OGlogisticSmoother.R")
source("OGLambdaGridSearch.R")

df_balanced <- read.csv("df_balanced.csv", stringsAsFactors = FALSE)

set.seed(42)

# caret prefers factor for stratification; convert back to 0/1 afterward
dat_fac <- df_balanced %>%
  mutate(risk = factor(risk, levels = c(0, 1)))

train_idx <- createDataPartition(dat_fac$risk, p = 0.70, list = FALSE)

train_data <- dat_fac[train_idx, ] %>%
  mutate(risk = as.integer(as.character(risk)))

test_data <- dat_fac[-train_idx, ] %>%
  mutate(risk = as.integer(as.character(risk)))

# *** NEW: distance column for training spills ***
train_data$dist_to_flowline <- as.numeric(train_data$match_distance_m)

cat("\nPrevalence check (should match overall ~5%/95%):\n")
cat("Overall p(Spill)= ", mean(as.integer(as.character(dat_fac$risk)) == 1), "\n")
cat("Train p(Spill)= ", mean(train_data$risk == 1), "\n")
cat("Test p(Spill)= ", mean(test_data$risk == 1), "\n")

# =============================================================================
# 3) Choose sigma: CV with simple 2-predictor GLM (same idea as before)
#    This is just a diagnostic step to get best_sigma + elbow_sigma.
# =============================================================================

# Helper to get log-likelihood
val_loglik <- function(y, p) {
  eps <- 1e-15
  p <- pmin(pmax(p, eps), 1 - eps)
  sum(y * log(p) + (1 - y) * log(1 - p))
}

# One-hot encode on TRAIN only, then keep fluidOther + elevation
dummies <- dummyVars(risk ~ ., data = train_data, fullRank = TRUE)
X_train <- predict(dummies, newdata = train_data) %>% as.data.frame()
y_train <- train_data$risk

Z_train <- X_train[, c("fluidOther", "elevation")]
train_sel <- cbind.data.frame(risk = y_train, Z_train)

set.seed(2025)
sigma_grid <- seq(5, 200, by = 5)
K <- 10
folds <- createFolds(train_sel$risk, k = K, list = TRUE)

cv_loglik <- numeric(length(sigma_grid))

for (s in seq_along(sigma_grid)) {
  sigma <- sigma_grid[s]
  fold_ll <- numeric(K)
  
  for (k in seq_len(K)) {
    val_idx <- folds[[k]]
    tr_idx  <- setdiff(seq_len(nrow(train_sel)), val_idx)
    
    tr <- train_sel[tr_idx, ]
    vl <- train_sel[val_idx, ]
    
    tr_dist <- train_data$dist_to_flowline[tr_idx]
    
    # exponential distance-based weights for spills
    w_tr <- ifelse(
      tr$risk == 1 & !is.na(tr_dist),
      exp(-tr_dist / sigma),
      1
    )
    
    fit <- glm(
      risk ~ fluidOther + elevation,
      data    = tr,
      family  = binomial(),
      weights = w_tr
    )
    
    p_hat <- predict(fit, newdata = vl, type = "response")
    fold_ll[k] <- val_loglik(vl$risk, p_hat)
  }
  
  cv_loglik[s] <- mean(fold_ll)
  cat("sigma =", sigma, "→ mean CV log-lik =", cv_loglik[s], "\n")
}

best_sigma <- sigma_grid[which.max(cv_loglik)]
cat("\nBest sigma (10-fold CV) =", best_sigma, "\n")

# Elbow: first sigma achieving 90% of improvement
ll_min <- cv_loglik[1]
ll_max <- max(cv_loglik)
total_gain <- ll_max - ll_min
f <- 0.9
threshold <- ll_min + f * total_gain

elbow_idx <- which(cv_loglik >= threshold)[1]
elbow_sigma <- sigma_grid[elbow_idx]
cat("Elbow sigma (90% gain) =", elbow_sigma, "\n")

# =============================================================================
# 4) Build weights for the 3 spatial models
# =============================================================================

dist_train <- train_data$dist_to_flowline

# Unweighted: all 1's
w_unweighted <- rep(1, nrow(train_data))

# Elbow sigma weights
w_elbow <- ifelse(
  train_data$risk == 1 & !is.na(dist_train),
  exp(-dist_train / elbow_sigma),
  1
)

# Best sigma weights
w_best <- ifelse(
  train_data$risk == 1 & !is.na(dist_train),
  exp(-dist_train / best_sigma),
  1
)

# Quick sanity check
summary(w_unweighted)
summary(w_elbow)
summary(w_best)

# =============================================================================
# 5) Smoother fit + lambda search for each weighting scheme
#    (Assumes LambdaGridSearch and logisticSmoother can accept weights=)
# =============================================================================

train_coord_model <- as.matrix(train_data[, c("lon","lat")])
train_risk_model  <- train_data$risk

# --- (a) Unweighted spatial model --------------------------------------------
profile_unw    <- LambdaGridSearch(train_coord_model, train_risk_model)  # no weights
bestLambda_unw <- profile_unw$bestLambda

fit_unw <- logisticSmoother(
  train_coord_model,
  train_risk_model,
  lambda  = bestLambda_unw
  # no weights argument
)

# --- (b) Elbow-weighted spatial model ----------------------------------------
profile_elbow <- LambdaGridSearch(
  train_coord_model,
  train_risk_model,
  weights = w_elbow              # *** NEW: pass weights ***
)
bestLambda_elbow <- profile_elbow$bestLambda

fit_elbow <- logisticSmoother(
  train_coord_model,
  train_risk_model,
  lambda  = bestLambda_elbow,
  weights = w_elbow              # *** NEW ***
)

# --- (c) Best-sigma-weighted spatial model -----------------------------------
profile_best <- LambdaGridSearch(
  train_coord_model,
  train_risk_model,
  weights = w_best               # *** NEW ***
)
bestLambda_best <- profile_best$bestLambda

fit_best <- logisticSmoother(
  train_coord_model,
  train_risk_model,
  lambda  = bestLambda_best,
  weights = w_best               # *** NEW ***
)

# =============================================================================
# 6) Predict on TEST for each model + log-likelihood
# =============================================================================

test_coord_model <- as.matrix(test_data[, c("lon","lat")])
test_risk_model  <- test_data$risk

evaluate_spatial <- function(fit, model_name, test_coords, test_data) {
  nu_hat <- predict(fit, test_coords)
  nu_hat <- as.numeric(nu_hat)
  p_hat  <- plogis(nu_hat)
  p_hat  <- pmin(pmax(p_hat, 1e-15), 1 - 1e-15)
  
  df <- test_data %>%
    mutate(predicted_prob = p_hat)
  
  y_true <- df$risk
  y_prob <- df$predicted_prob
  y_pred <- as.integer(y_prob >= 0.5)
  
  loglik_i <- y_true * log(y_prob) +
    (1 - y_true) * log(1 - y_prob)
  loglik   <- sum(loglik_i)
  rmse     <- sqrt(mean((y_prob - y_true)^2))
  
  tp <- sum(y_pred == 1 & y_true == 1)
  tn <- sum(y_pred == 0 & y_true == 0)
  fp <- sum(y_pred == 1 & y_true == 0)
  fn <- sum(y_pred == 0 & y_true == 1)
  
  accuracy  <- (tp + tn) / length(y_true)
  precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_)
  recall    <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
  f1        <- ifelse((precision + recall) > 0,
                      2 * precision * recall / (precision + recall),
                      NA_real_)
  sensitivity <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
  specificity <- ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_)
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  
  mcc <- ifelse(
    (tp + fp) > 0 & (tp + fn) > 0 & (tn + fp) > 0 & (tn + fn) > 0,
    (tp * tn - fp * fn) /
      sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)),
    NA_real_
  )
  
  metrics <- tibble(
    model             = model_name,
    set               = "TEST",
    n                 = length(y_true),
    log_likelihood    = loglik,
    rmse              = rmse,
    accuracy          = accuracy,
    precision         = precision,
    recall            = recall,
    f1                = f1,
    sensitivity       = sensitivity,
    specificity       = specificity,
    balanced_accuracy = balanced_accuracy,
    mcc               = mcc,
    threshold         = 0.5
  )
  
  list(pred_df = df, metrics = metrics)
}

res_unw   <- evaluate_spatial(fit_unw,   "Spatial_unweighted", test_coord_model, test_data)
res_elbow <- evaluate_spatial(fit_elbow, "Spatial_elbow",      test_coord_model, test_data)
res_best  <- evaluate_spatial(fit_best,  "Spatial_best",       test_coord_model, test_data)

all_metrics <- bind_rows(
  res_unw$metrics,
  res_elbow$metrics,
  res_best$metrics
)
print(all_metrics)

# =============================================================================
# 7) Save outputs (predictions + metrics)
# =============================================================================

write_csv(
  res_unw$pred_df %>%
    transmute(model = "Spatial_unweighted", lon, lat, risk, pred = predicted_prob),
  "predictions_spatial_unweighted.csv"
)
write_csv(
  res_elbow$pred_df %>%
    transmute(model = "Spatial_elbow", lon, lat, risk, pred = predicted_prob),
  "predictions_spatial_elbow.csv"
)
write_csv(
  res_best$pred_df %>%
    transmute(model = "Spatial_best", lon, lat, risk, pred = predicted_prob),
  "predictions_spatial_best.csv"
)

write_csv(all_metrics, "metrics_spatial_all.csv")

# =============================================================================
# 8) Plots – e.g., unweighted vs elbow vs best
# =============================================================================
# Example: predicted map for unweighted (you can clone this for elbow/best)

ggplot(res_unw$pred_df, aes(x = lon, y = lat)) +
  geom_point(aes(color = predicted_prob), size = 2) +
  scale_color_viridis_c(
    option = "turbo",
    name   = "Predicted Risk",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    oob    = scales::squish
  ) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "Spatial-only Predicted Spill Risk – Unweighted",
    x     = "Longitude",
    y     = "Latitude"
  )
