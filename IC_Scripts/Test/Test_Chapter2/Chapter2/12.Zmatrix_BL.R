# =============================================================================
# 1) Setup & Libraries
# =============================================================================
library(sf)
library(dplyr)
library(ggplot2)
library(caret)
library(fields)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# =============================================================================
# 2) Sources & Data
# =============================================================================
source("ZlogisticSmoother.R")       
source("ZLambdaGridSearch.R")

df_balanced <- read.csv("df_balanced.csv", stringsAsFactors = FALSE)

# =============================================================================
# 3) Stratified Train/Test Split
# =============================================================================
set.seed(42)
dat_fac <- df_balanced %>% mutate(risk = factor(risk, levels = c(0, 1)))

train_idx <- createDataPartition(dat_fac$risk, p = 0.70, list = FALSE)

train_data <- dat_fac[train_idx, ] %>%
  mutate(risk = as.integer(as.character(risk)))

test_data <- dat_fac[-train_idx, ] %>%
  mutate(risk = as.integer(as.character(risk)))

cat("\nPrevalence check:\n")
cat("Overall p(Spill) =", mean(as.integer(as.character(dat_fac$risk)) == 1), "\n")
cat("Train p(Spill)   =", mean(train_data$risk == 1), "\n")
cat("Test p(Spill)    =", mean(test_data$risk == 1), "\n")

# =============================================================================
# 4) Build X and Z matrices FOR EACH SPLIT
# =============================================================================
to.factor <- c("status","flowline_action","location_type","fluid","material")

## Training
train_data <- train_data %>% mutate(across(all_of(to.factor), as.factor))
X_train <- model.matrix(
  ~ status + flowline_action + location_type + fluid + material +
    diameter_in + length_ft + elevation + line_age_yr,
  data = train_data
)[, -1]

Z_train <- X_train[, c("fluidOther", "elevation")]
Z_df_train <- as.data.frame(Z_train)
Z_df_train$risk <- train_data$risk

train_risk_model   <- train_data$risk
train_coord_model  <- as.matrix(train_data[, c("lon","lat")])

glm_sel <- glm(risk ~ . -1, data = Z_df_train, family = binomial(link = "logit"))
betaHat <- coef(glm_sel)

## Testing
test_data <- test_data %>% mutate(across(all_of(to.factor), as.factor))
X_test <- model.matrix(
  ~ status + flowline_action + location_type + fluid + material +
    diameter_in + length_ft + elevation + line_age_yr,
  data = test_data
)[, -1]

Z_test <- X_test[, c("fluidOther", "elevation")]
Z_df_test <- as.data.frame(Z_test)
Z_df_test$risk <- test_data$risk  

test_risk_model   <- test_data$risk
test_coord_model  <- as.matrix(test_data[, c("lon","lat")])

# =============================================================================
# 7) Lambda search + final fit with Z and betaHat
# =============================================================================
profile <- LambdaGridSearch(
  train_coord_model,
  train_risk_model,
  Z       = Z_train,
  betaHat = betaHat
)

bestLambda <- profile$bestLambda

MLEFit <- logisticSmoother(
  train_coord_model,
  train_risk_model,
  lambda  = bestLambda,
  Z       = Z_train,
  betaHat = betaHat
)

# =============================================================================
# 7) Predict from fitted model (with Z)
# =============================================================================
nu_all_mle   <- predict(MLEFit, test_coord_model, Z = Z_test)
nu_all_mle   <- as.numeric(nu_all_mle)
prob_all_mle <- plogis(nu_all_mle)
prob_all_mle <- as.numeric(prob_all_mle)

test_data <- test_data %>%
  mutate(predicted_prob = prob_all_mle)

spatial_Z_preds <- test_data %>%
  transmute(
    model = "SpatialPlusZ",
    idx   = row_number(),
    lon   = as.numeric(lon),
    lat   = as.numeric(lat),
    risk  = as.integer(risk),
    pred  = as.numeric(predicted_prob),
    lo    = NA_real_,
    hi    = NA_real_
  )

readr::write_csv(spatial_Z_preds, "predictions_spatial_plus_z.csv")

# ============================================================================
# 8) Total log-likelihood
# =============================================================================
nu_hat <- predict(MLEFit, test_coord_model, Z = Z_test)
p_hat  <- exp(nu_hat) / (1 + exp(nu_hat))  # same as plogis(nu_hat)
risk_model <- test_data$risk

loglik_i <- risk_model * log(p_hat) +
  (1 - risk_model) * log(1 - p_hat)

total_loglik <- sum(loglik_i)

cat("Total log-likelihood =", round(total_loglik, 2), "\n")

# =============================================================================
# 9) Observed Risk plot
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
    values = c("0" = "skyblue", "1" = "firebrick"),
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
# 10) Predicted Risk plot
# =============================================================================
ggplot(test_data, aes(x = lon, y = lat)) +
  geom_point(aes(color = predicted_prob), size = 2) +
  scale_color_gradient(
    low   = "skyblue",
    high  = "firebrick",
    name  = "Predicted Risk",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    oob    = scales::squish
  ) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "Predicted Spill Risk – With Z Matrix",
    x     = "Longitude",
    y     = "Latitude"
  )


# Metrics (same recipe)
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

spatial_Z_metrics <- tibble::tibble(
  model         = "SpatialPlusZ",
  set           = "TEST",
  n             = length(y_true),
  log_likelihood = loglik,
  rmse          = rmse,
  accuracy      = accuracy,
  precision     = precision,
  recall        = recall,
  f1            = f1,
  threshold     = 0.5
)

readr::write_csv(spatial_Z_metrics, "metrics_spatial_plus_z.csv")