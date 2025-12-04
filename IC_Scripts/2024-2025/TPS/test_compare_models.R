# 16.compare_models.R
# =============================================================================
# Compare GLM, Spatial-only, Spatial+Z, and Spatial+Z+weights
# using log-likelihood, Brier score, ROC-AUC + plots
# =============================================================================

# 1) Setup --------------------------------------------------------------------
library(dplyr)
library(readr)
library(ggplot2)
library(pROC)
library(purrr)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

all_preds <- read_csv("all_predictions_combined.csv")

# 3) Helper functions ---------------------------------------------------------
loglik_bernoulli <- function(y, p) {
  eps <- 1e-15
  p <- pmin(pmax(p, eps), 1 - eps)
  sum(y * log(p) + (1 - y) * log(1 - p))
}

brier_score <- function(y, p) {
  mean((p - y)^2)
}

auc_roc <- function(y, p) {
  # y must be 0/1 numeric or factor; p = predicted prob
  r <- pROC::roc(response = y, predictor = p, quiet = TRUE)
  as.numeric(r$auc)
}

# 4) Compute metrics per model -----------------------------------------------
metrics_by_model <- all_preds %>%
  group_by(model) %>%
  summarise(
    n      = n(),
    loglik = loglik_bernoulli(risk, pred),
    brier  = brier_score(risk, pred),
    auc    = auc_roc(risk, pred),
    z_spieg = spiegelhalter_z(risk, pred),
    p_spieg = spiegelhalter_p(z_spieg),
    .groups = "drop"
  ) %>%
  arrange(desc(loglik))   # higher loglik (less negative) is better

print(metrics_by_model)

# 5) ROC curves for all models (ggplot) --------------------------------------
# Build ROC curve data for each model
roc_df_list <- all_preds %>%
  split(.$model) %>%
  imap(function(df, mname) {
    r <- pROC::roc(response = df$risk, predictor = df$pred, quiet = TRUE)
    tibble(
      model = mname,
      fpr   = 1 - r$specificities,
      tpr   = r$sensitivities
    )
  })

roc_df <- bind_rows(roc_df_list)

gg_roc <- ggplot(roc_df, aes(x = fpr, y = tpr, color = model)) +
  geom_line(linewidth = 1.1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "ROC Curves: GLM vs Spatial Models",
    x     = "False Positive Rate (1 - Specificity)",
    y     = "True Positive Rate (Sensitivity)",
    color = "Model"
  )

print(gg_roc)

# 6) Spatial maps of predicted risk ------------------------------------------
# You already have lon/lat in the prediction tables.

gg_spatial <- ggplot(all_preds, aes(x = lon, y = lat)) +
  geom_point(aes(color = pred), size = 2) +
  scale_color_viridis_c(
    option = "turbo",
    name   = "Predicted Risk",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    oob    = scales::squish,
    begin  = 0.12,
    end    = 0.9
  ) +
  coord_fixed() +
  facet_wrap(~ model, ncol = 2) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Predicted Spill Probability by Model (TEST set)",
    x     = "Longitude",
    y     = "Latitude"
  )

print(gg_spatial)

# 7) (Optional) Observed vs predicted probability calibration ----------------
# Quick calibration-style scatter or smooth for each model
gg_calib <-
  ggplot(all_preds, aes(x = pred, y = risk)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess", se = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = "longdash") +
  coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(~ model, ncol = 2) +
  theme_minimal() +
  labs(
    title = "Calibration Curves (Observed vs Predicted)",
    x = "Predicted spill probability",
    y = "Observed spill rate"
  )

print(gg_calib)

# 8) Save metrics table to CSV -----------------------------------------------
write_csv(metrics_by_model, "model_comparison_metrics.csv")

cat("\nSaved metrics to model_comparison_metrics.csv\n")
