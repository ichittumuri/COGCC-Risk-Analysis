# =============================================================================
# Exponential distance weighting diagnostic:
# Show that CV prefers very large sigma (i.e. almost no downweighting)
# =============================================================================

# 1) Setup & Libraries --------------------------------------------------------
library(dplyr)
library(caret)
library(ggplot2)
library(readr)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# 2) Load data ---------------------------------------------------------------
df_balanced <- read.csv("df_balanced.csv", stringsAsFactors = FALSE)

# (optional) formula kept for reference; we only use fluidOther + elevation
form <- risk ~ status + flowline_action + location_type + fluid + material +
  diameter_in + length_ft + max_operating_pressure + elevation + line_age_yr

# 3) Train–test split (same as main script) ----------------------------------
set.seed(42)

dat_fac <- df_balanced %>%
  mutate(risk = factor(risk, levels = c(0, 1)))

train_idx <- createDataPartition(dat_fac$risk, p = 0.70, list = FALSE)

train_data <- dat_fac[train_idx, ] %>%
  mutate(risk = as.integer(as.character(risk)))
test_data  <- dat_fac[-train_idx, ] %>%
  mutate(risk = as.integer(as.character(risk)))

# Distance column (only spills have non-NA)
train_data$dist_to_flowline <- as.numeric(train_data$match_distance_m)
test_data$dist_to_flowline <- as.numeric(test_data$match_distance_m)

# 4) One-hot encode on TRAIN only and keep fluidOther + elevation -----------
dummies <- dummyVars(risk ~ ., data = train_data, fullRank = TRUE)

X_train <- predict(dummies, newdata = train_data) %>% as.data.frame()
y_train <- train_data$risk

# choose the two predictors selected by LASSO
Z_train <- X_train[, c("fluidOther", "elevation")]

# bind into a clean training frame
train_sel <- cbind.data.frame(risk = y_train, Z_train)

# 5) Helper: validation log-likelihood ---------------------------------------
val_loglik <- function(y, p) {
  eps <- 1e-15
  p <- pmin(pmax(p, eps), 1 - eps)
  sum(y * log(p) + (1 - y) * log(1 - p))
}

# 6) K-fold CV over a big sigma grid (exponential weights) -------------------
set.seed(2025)

sigma_grid <- seq(5, 200, by = 5)   # "a ton of sigmas"; increase if you want
K <- 10                              # 5-fold CV (stratified on risk)

folds <- createFolds(train_sel$risk, k = K, list = TRUE)

cv_loglik <- numeric(length(sigma_grid))

for (s in seq_along(sigma_grid)) {
  sigma <- sigma_grid[s]
  fold_ll <- numeric(K)
  
  for (k in 1:K) {
    # split indices
    val_idx <- folds[[k]]
    tr_idx  <- setdiff(seq_len(nrow(train_sel)), val_idx)
    
    tr <- train_sel[tr_idx, ]
    vl <- train_sel[val_idx, ]
    
    # distances for these rows (from original train_data)
    tr_dist <- train_data$dist_to_flowline[tr_idx]
    
    # exponential distance-based weights for spills
    w_tr <- ifelse(
      tr$risk == 1 & !is.na(tr_dist),
      exp(-tr_dist / sigma),
      1
    )
    
    # fit weighted logistic regression
    fit <- glm(
      risk ~ fluidOther + elevation,
      data    = tr,
      family  = binomial(),
      weights = w_tr
    )
    
    # predict on validation fold
    p_hat <- predict(fit, newdata = vl, type = "response")
    
    # validation log-likelihood
    fold_ll[k] <- val_loglik(vl$risk, p_hat)
  }
  
  cv_loglik[s] <- mean(fold_ll)
  cat("sigma =", sigma, "→ mean CV log-lik =", cv_loglik[s], "\n")
}

# 7) Find best sigma and show that it's at the upper edge --------------------
best_sigma <- sigma_grid[which.max(cv_loglik)]
cat("\nBest sigma (5-fold CV, exponential weights) =", best_sigma, "\n")

# ----- simple elbow rule: 90% of max improvement -----
ll_min <- cv_loglik[1]
ll_max <- max(cv_loglik)
total_gain <- ll_max - ll_min

f <- 0.9  # fraction of total gain you want (try 0.9, 0.95, etc.)
threshold <- ll_min + f * total_gain

elbow_idx <- which(cv_loglik >= threshold)[1]
elbow_sigma <- sigma_grid[elbow_idx]
elbow_sigma

# 8) Plot: sigma vs mean CV log-likelihood -----------------------------------
cv_df <- data.frame(sigma = sigma_grid, cv_ll = cv_loglik)

ggplot(cv_df, aes(x = sigma, y = cv_ll)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = best_sigma, linetype = "dashed", color = "red") +
  geom_vline(xintercept = elbow_sigma, linetype = "dotted", color = "blue") +
  labs(
    title = "Exponential distance weighting:\nCV log-likelihood vs sigma",
    x = "Sigma (m)",
    y = "Mean CV log-likelihood (10-fold)"
  ) +
  theme_minimal()
