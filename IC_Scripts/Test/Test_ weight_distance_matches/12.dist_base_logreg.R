# =============================================================================
# 1) Setup & Libraries
# =============================================================================
library(dplyr)
library(ggplot2)
library(scales)
library(caret)
library(readr)
library(tibble)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# =============================================================================
# 2) Load Data & Define Formula
# =============================================================================
df_balanced <- read.csv("df_balanced.csv", stringsAsFactors = FALSE)

form <- risk ~ status + flowline_action + location_type + fluid + material +
  diameter_in + length_ft + max_operating_pressure + elevation + line_age_yr

# =============================================================================
# 3) Train-Test Split (70/30, Stratified by Risk)
# =============================================================================
set.seed(42)

dat_fac <- df_balanced %>%
  mutate(risk = factor(risk, levels = c(0, 1)))

train_idx <- createDataPartition(dat_fac$risk, p = 0.70, list = FALSE)

train_data <- dat_fac[train_idx, ] %>%
  mutate(risk = as.integer(as.character(risk)))

test_data <- dat_fac[-train_idx, ] %>%
  mutate(risk = as.integer(as.character(risk)))

# Distance column (only spills have non-NA)
train_data$dist_to_flowline <- as.numeric(train_data$match_distance_m)
test_data$dist_to_flowline <- as.numeric(test_data$match_distance_m)

# =============================================================================
# Prevalence check + sizes
# =============================================================================
cat("\nPrevalence check (should match overall ~5%/95%):\n")

summarize_split <- function(data, name) {
  n_total <- nrow(data)
  n_spill <- sum(data$risk == 1)
  n_nosp  <- sum(data$risk == 0)
  pct_spill <- mean(data$risk == 1)
  
  cat(sprintf("%s set: n = %d  |  spills = %d  |  no-spills = %d  |  p(Spill) = %.3f\n",
              name, n_total, n_spill, n_nosp, pct_spill))
}
summarize_split(dat_fac,   "Overall")
summarize_split(train_data, "Train")
summarize_split(test_data,  "Test")

# =============================================================================
# 4) One-hot encode categorical variables (TRAIN and TEST)
# =============================================================================
dummies <- dummyVars(risk ~ ., data = train_data, fullRank = TRUE)

X_train <- predict(dummies, newdata = train_data) %>% as.data.frame()
X_test  <- predict(dummies, newdata = test_data) %>% as.data.frame()

y_train <- train_data$risk
y_test  <- test_data$risk

# Keep only the two LASSO-selected predictors
Z_train <- X_train[, c("fluidOther", "elevation")]
Z_test  <- X_test[, c("fluidOther", "elevation")]

train_sel <- cbind.data.frame(risk = y_train, Z_train)
test_sel  <- cbind.data.frame(risk = y_test,  Z_test)

write.csv(train_sel, "train_sel.csv", row.names = FALSE)
write.csv(test_sel,  "test_sel.csv",  row.names = FALSE)

# =============================================================================
# 5) CV to choose best sigma and elbow sigma for distance weights
# =============================================================================
val_loglik <- function(y, p) {
  eps <- 1e-15
  p <- pmin(pmax(p, eps), 1 - eps)
  sum(y * log(p) + (1 - y) * log(1 - p))
}

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

# Best sigma (pure CV optimum)
best_sigma <- sigma_grid[which.max(cv_loglik)]
cat("\nBest sigma (10-fold CV) =", best_sigma, "\n")

# Elbow sigma: first sigma achieving 90% of total improvement
ll_min <- cv_loglik[1]
ll_max <- max(cv_loglik)
total_gain <- ll_max - ll_min
f <- 0.9
threshold <- ll_min + f * total_gain

elbow_idx <- which(cv_loglik >= threshold)[1]
elbow_sigma <- sigma_grid[elbow_idx]
cat("Elbow sigma (90% gain) =", elbow_sigma, "\n")

# =============================================================================
# 6) Fit three GLMs: unweighted, elbow-weighted, best-weighted
# =============================================================================
dist_train <- train_data$dist_to_flowline

# unweighted (baseline)
glm_unweighted <- glm(
  risk ~ fluidOther + elevation,
  data   = train_sel,
  family = binomial()
)

# elbow-weighted
w_elbow <- ifelse(
  train_sel$risk == 1 & !is.na(dist_train),
  exp(-dist_train / elbow_sigma),
  1
)
glm_elbow <- glm(
  risk ~ fluidOther + elevation,
  data    = train_sel,
  family  = binomial(),
  weights = w_elbow
)

# best-weighted (CV optimum)
w_best <- ifelse(
  train_sel$risk == 1 & !is.na(dist_train),
  exp(-dist_train / best_sigma),
  1
)
glm_best <- glm(
  risk ~ fluidOther + elevation,
  data    = train_sel,
  family  = binomial(),
  weights = w_best
)

cat("\nGLM summaries (train fit):\n")
print(summary(glm_unweighted))
print(summary(glm_elbow))
print(summary(glm_best))

# =============================================================================
# 7) Evaluation helper (applies to any fitted glm)
# =============================================================================
evaluate_model <- function(fit, model_name, test_sel, test_data) {
  # predictions
  p_hat <- predict(fit, newdata = test_sel, type = "response")
  eps  <- 1e-15
  p_hat <- pmin(pmax(p_hat, eps), 1 - eps)
  
  test_out <- test_sel %>%
    mutate(pred_prob = p_hat) %>%
    bind_cols(test_data %>% dplyr::select(lon, lat))
  
  y_true <- test_out$risk
  y_prob <- test_out$pred_prob
  y_pred <- as.integer(y_prob >= 0.5)
  
  loglik_i <- y_true * log(y_prob) + (1 - y_true) * log(1 - y_prob)
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
  
  list(test_out = test_out, metrics = metrics)
}

# =============================================================================
# 8) Evaluate all three models on TEST set
# =============================================================================
res_unw   <- evaluate_model(glm_unweighted, "GLM_unweighted", test_sel, test_data)
res_elbow <- evaluate_model(glm_elbow,      "GLM_elbow",      test_sel, test_data)
res_best  <- evaluate_model(glm_best,       "GLM_best",       test_sel, test_data)

# Combine metrics
all_metrics <- bind_rows(
  res_unw$metrics,
  res_elbow$metrics,
  res_best$metrics
)

print(all_metrics)
