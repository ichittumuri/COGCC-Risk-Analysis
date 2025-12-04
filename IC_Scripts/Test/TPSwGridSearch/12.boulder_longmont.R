# New cell= command, option, I or alt, windows symbol, I 
# Run = windows symbol, enter 
# Run rscript = Command + Option + R

# =============================================================================
# 1) Setup & Libraries
# =============================================================================
library(sf)
library(dplyr)
library(ggplot2)
library(fields)
library(caret)

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

cat("\nPrevalence check (should match overall ~5%/95%):\n")
cat("Overall p(Spill)= ", mean(as.integer(as.character(dat_fac$risk)) == 1), "\n")
cat("Train p(Spill)= ", mean(train_data$risk == 1), "\n")
cat("Test p(Spill)= ", mean(test_data$risk == 1), "\n")

# =============================================================================
# 6) Smoother fit + lambda search
# =============================================================================
train_coord_model <- as.matrix(train_data[, c("lon","lat")])
train_risk_model  <- train_data$risk

profile    <- LambdaGridSearch(train_coord_model, train_risk_model)
bestLambda <- profile$bestLambda

train_MLEFit <- logisticSmoother(
  train_coord_model,
  train_risk_model,
  lambda = bestLambda
)

test_coord_model <- as.matrix(test_data[, c("lon","lat")])
test_risk_model  <- test_data$risk

# =============================================================================
# 7) Predict from fitted model
# =============================================================================
nu_all_mle   <- predict(train_MLEFit, test_coord_model)
nu_all_mle   <- as.numeric(nu_all_mle)               # drop 1-col matrix -> vector
prob_all_mle <- plogis(nu_all_mle)                   # same as exp(nu)/(1+exp(nu))
prob_all_mle <- as.numeric(prob_all_mle)             # ensure plain numeric

test_data <- test_data %>%
  mutate(predicted_prob = prob_all_mle)

# ---- Standardized outputs (Spatial-only) ------------------------------------
spatial_only_preds <- test_data %>%
  transmute(
    model = "SpatialOnly",
    idx   = row_number(),
    lon   = as.numeric(lon),
    lat   = as.numeric(lat),
    risk  = as.integer(risk),
    pred  = as.numeric(predicted_prob),
    lo    = NA_real_,
    hi    = NA_real_
  )

readr::write_csv(spatial_only_preds, "predictions_spatial_only.csv")

# =============================================================================
# 8) Total log-likelihood
# =============================================================================
nu_hat <- predict(train_MLEFit, test_coord_model)
p_hat  <- exp(nu_hat) / (1 + exp(nu_hat))  # same as plogis(nu_hat)
risk_model <- test_data$risk

loglik_i <- risk_model * log(p_hat) +
  (1 - risk_model) * log(1 - p_hat)

total_loglik <- sum(loglik_i)

cat("Total log-likelihood =", round(total_loglik, 2), "\n")

# =============================================================================
# 9) Observed Risk plot (categorical 0/1)
# =============================================================================
ggplot() +
  geom_point(
    data = test_data %>% filter(risk == 0),
    aes(x = lon, y = lat, color = factor(risk)),
    size = 1.7, alpha = 0.8
  ) +
  geom_point(
    data = test_data %>% filter(risk == 1),
    aes(x = lon, y = lat, color = factor(risk)),
    size = 1.7, alpha = 0.8
  ) +
  scale_color_manual(
    values = c("0" = "#42a5f5", "1" = "#e53935"),  # blue vs red
    name   = "Observed Risk",
    labels = c("0" = "No Spill", "1" = "Spill")
  ) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "Observed Spill Risk",
    x = "Longitude",
    y = "Latitude"
  )

# =============================================================================
# 10) Predicted Risk plot (Turbo continuous scale)
# =============================================================================
ggplot(test_data, aes(x = lon, y = lat)) +
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
    title = "Predicted Spill Risk – Without Z Matrix",
    x     = "Longitude",
    y     = "Latitude"
  )

# Metrics (same as GLM section)
y_true <- test_data$risk
y_prob <- test_data$predicted_prob
y_pred <- as.integer(y_prob >= 0.5)

loglik_i <- y_true * log(pmax(pmin(y_prob, 1 - 1e-15), 1e-15)) +
  (1 - y_true) * log(pmax(pmin(1 - y_prob, 1 - 1e-15), 1e-15))
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

# --- Add these: sensitivity, specificity, balanced accuracy, MCC -------------
sensitivity <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)  # recall for positives
specificity <- ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_)  # recall for negatives
balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
balanced_accuracy

mcc <- ifelse(
  (tp + fp) > 0 & (tp + fn) > 0 & (tn + fp) > 0 & (tn + fn) > 0,
  (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)),
  NA_real_
)

spatial_only_metrics <- tibble::tibble(
  model           = "SpatialOnly",
  set             = "TEST",
  n               = length(y_true),
  log_likelihood  = loglik,
  rmse            = rmse,
  accuracy        = accuracy,
  precision       = precision,
  recall          = recall,
  f1              = f1,
  sensitivity     = sensitivity,
  specificity     = specificity,
  balanced_accuracy = balanced_accuracy,
  mcc             = mcc,
  threshold       = 0.5
)
readr::write_csv(spatial_only_metrics, "metrics_spatial_only.csv")
