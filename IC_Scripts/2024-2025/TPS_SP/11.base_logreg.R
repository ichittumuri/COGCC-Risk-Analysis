# =============================================================================
# 1) Setup & Libraries
# =============================================================================
library(dplyr)
library(ggplot2)
library(scales)
library(caret)
library(readr)

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

# =============================================================================
# Prevalence check + sizes
# =============================================================================
cat("\nPrevalence check (should match overall ~5%/95%):\n")

# Helper function to summarize a dataset
summarize_split <- function(data, name) {
  n_total <- nrow(data)
  n_spill <- sum(data$risk == 1)
  n_nosp  <- sum(data$risk == 0)
  pct_spill <- mean(data$risk == 1)
  
  cat(sprintf("%s set: n = %d  |  spills = %d  |  no-spills = %d  |  p(Spill) = %.3f\n",
              name, n_total, n_spill, n_nosp, pct_spill))
}

# Print summaries
summarize_split(dat_fac,  "Overall")
summarize_split(train_data, "Train")
summarize_split(test_data,  "Test")

# =============================================================================
# 1) One-hot encode categorical variables
# =============================================================================
dummies <- dummyVars(risk ~ ., data = train_data, fullRank = TRUE)

X_train <- predict(dummies, newdata = train_data) %>% as.data.frame()
X_test  <- predict(dummies, newdata = test_data) %>% as.data.frame()

y_train <- train_data$risk
y_test  <- test_data$risk

Z_train <- X_train[, c("fluidOther", "elevation")]
Z_test  <- X_test[, c("fluidOther", "elevation")]

# Bind response for fitting
train_sel <- cbind.data.frame(risk = y_train, Z_train)
test_sel  <- cbind.data.frame(risk = y_test,  Z_test)

# Save the datasets as CSVs
write.csv(train_sel, "train_sel.csv", row.names = FALSE)
write.csv(test_sel,  "test_sel.csv",  row.names = FALSE)

# =============================================================================
# 4) Fit Baseline GLM (No LASSO)
# =============================================================================
glm_base <- glm(risk ~ fluidOther + elevation, data = train_sel, family = binomial())
cat("\nGLM summary (train fit):\n")
print(summary(glm_base))

# =============================================================================
# 5) Predictions on TEST Set + 95% CIs
# =============================================================================
nu_hat <-predict(glm_base, newdata = test_sel, type = "response")
p_hat  <- exp(nu_hat) / (1 + exp(nu_hat)) 
risk_model <- test_sel$risk

loglik_i <- risk_model * log(p_hat) +
  (1 - risk_model) * log(1 - p_hat)

total_loglik <- sum(loglik_i)

cat("Total log-likelihood =", round(total_loglik, 2), "\n")

# =============================================================================
# 5) Predictions on TEST Set (no CIs)
# =============================================================================
p_hat <- predict(glm_base, newdata = test_sel, type = "response")

# Numerical safety
eps  <- 1e-15
p_hat <- pmin(pmax(p_hat, eps), 1 - eps)

# Build a test_out for consistency
test_out <- test_sel %>%
  mutate(pred_prob = p_hat) %>%
  # attach lon/lat from the original test_data by row
  bind_cols(test_data %>% dplyr::select(lon, lat))

# =============================================================================
# 6) Evaluation Metrics (TEST)
# =============================================================================
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

conf_mat <- data.frame(
  Actual_0 = c(tn, fn),
  Actual_1 = c(fp, tp),
  row.names = c("Pred_0", "Pred_1")
)

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

metrics <- tibble::tibble(
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

# =============================================================================
# 7) Save Outputs (no CI columns)
# =============================================================================
readr::write_csv(
  test_out %>% transmute(lon, lat, risk, pred_prob),
  "test_predictions_glm70_30.csv"
)

readr::write_csv(metrics, "metrics_glm70_30.csv")

conf_mat_out <- tibble::tibble(
  Predicted = c("Pred_0", "Pred_1"),
  Actual_0  = c(tn, fn),
  Actual_1  = c(fp, tp)
)
readr::write_csv(conf_mat_out, "confusion_matrix_glm70_30.csv")

# Standardized outputs
glm_preds <- test_out %>%
  transmute(model = "GLM", idx = row_number(),
            lon, lat, risk = as.integer(risk),
            pred = pred_prob)
readr::write_csv(glm_preds, "predictions_glm.csv")

glm_metrics <- metrics %>%
  mutate(model = "GLM") %>%
  select(model, everything())
readr::write_csv(glm_metrics, "metrics_glm.csv")

# =============================================================================
# 8) Visualization (Turbo color scheme)
# =============================================================================
df_sel <- test_out %>%
  transmute(lon, lat, risk = as.integer(risk), predicted_prob = pred_prob)

# Observed Map
ggplot() +
  geom_point(data = df_sel %>% filter(risk == 0),
             aes(x = lon, y = lat, color = factor(risk)),
             size = 1.7, alpha = 0.8) +
  geom_point(data = df_sel %>% filter(risk == 1),
             aes(x = lon, y = lat, color = factor(risk)),
             size = 1.7, alpha = 0.8) +
  scale_color_manual(values = c("0" = "#42a5f5", "1" = "#e53935"),  # blue vs red
                     name = "Observed Risk",
                     labels = c("No Spill", "Spill")) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Observed Spill Risk (TEST set)",
       x = "Longitude", y = "Latitude")

# Predicted Map with Turbo scale
ggplot(df_sel %>% arrange(predicted_prob),
       aes(x = lon, y = lat)) +
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
  labs(title = "Predicted Spill Risk – GLM (TEST set)",
       x = "Longitude", y = "Latitude")
