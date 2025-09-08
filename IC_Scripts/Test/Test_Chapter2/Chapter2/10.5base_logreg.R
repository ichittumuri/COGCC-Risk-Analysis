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

cat("\nPrevalence check (should match overall ~5%/95%):\n")
cat("Overall p(Spill) =", mean(as.integer(as.character(dat_fac$risk)) == 1), "\n")
cat("Train p(Spill)   =", mean(train_data$risk == 1), "\n")
cat("Test p(Spill)    =", mean(test_data$risk == 1), "\n")

# =============================================================================
# 4) Fit Baseline GLM (No LASSO)
# =============================================================================
glm_base <- glm(form, data = train_data, family = binomial())
cat("\nGLM summary (train fit):\n")
print(summary(glm_base))

# =============================================================================
# 5) Predictions on TEST Set + 95% CIs
# =============================================================================
pred_link_test <- predict(glm_base, newdata = test_data, type = "link", se.fit = TRUE)
eta    <- pred_link_test$fit
se_eta <- pred_link_test$se.fit

p_hat <- plogis(eta)
p_lo  <- plogis(eta - 1.96 * se_eta)
p_hi  <- plogis(eta + 1.96 * se_eta)

# Numerical safety
eps <- 1e-15
p_hat <- pmin(pmax(p_hat, eps), 1 - eps)
p_lo  <- pmin(pmax(p_lo, eps), 1 - eps)
p_hi  <- pmin(pmax(p_hi, eps), 1 - eps)

test_out <- test_data %>%
  mutate(pred_prob = p_hat, ci_lo = p_lo, ci_hi = p_hi)

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

metrics <- tibble::tibble(
  set = "TEST",
  n = length(y_true),
  log_likelihood = loglik,
  rmse = rmse,
  accuracy = accuracy,
  precision = precision,
  recall = recall,
  f1 = f1,
  threshold = 0.5
)

cat("\nTEST metrics:\n")
print(metrics)
cat("\nConfusion matrix (rows=Pred, cols=Actual):\n")
print(conf_mat)

# =============================================================================
# 7) Save Outputs
# =============================================================================
# Predictions with CI
readr::write_csv(
  test_out %>% transmute(lon, lat, risk, pred_prob, ci_lo, ci_hi),
  "test_predictions_glm70_30.csv"
)

# Metrics
readr::write_csv(metrics, "metrics_glm70_30.csv")

# Confusion matrix
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
            pred = pred_prob, lo = ci_lo, hi = ci_hi)
readr::write_csv(glm_preds, "predictions_glm.csv")

glm_metrics <- metrics %>%
  mutate(model = "GLM") %>%
  select(model, everything())
readr::write_csv(glm_metrics, "metrics_glm.csv")

# =============================================================================
# 8) Visualization
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
  scale_color_manual(values = c("0" = "skyblue", "1" = "firebrick"),
                     name = "Observed Risk",
                     labels = c("No Spill", "Spill")) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Observed Spill Risk (TEST set)",
       x = "Longitude", y = "Latitude")

# Predicted Map
ggplot(df_sel %>% arrange(predicted_prob),  # Low → high so high overlays on top
       aes(x = lon, y = lat)) +
  geom_point(aes(color = predicted_prob), size = 2) +
  scale_color_gradient(low = "skyblue", high = "firebrick",
                       name = "Predicted Risk",
                       limits = c(0, 1),
                       breaks = seq(0, 1, 0.25),
                       oob = scales::squish) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Predicted Spill Risk – GLM (TEST set)",
       x = "Longitude", y = "Latitude")


