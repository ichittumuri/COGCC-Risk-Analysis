# =============================================================================
# 16.plot_train_test_and_predictions.R
#   - Observed maps: Train + Test (2 panels)
#   - Predicted maps: 2x2 (GLM, SpatialOnly, SpatialPlusZ, SpatialPlusZ_wt)
#   Inputs:
#     - train_data_split.csv, test_data_split.csv
#     - all_predictions.csv
# =============================================================================
library(dplyr)
library(readr)
library(ggplot2)
library(scales)
library(tidyr)
library(patchwork)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")
theme_set(theme_grey())

# ---- Load data --------------------------------------------------------------
train_data <- read_csv("train_data_split.csv", show_col_types = FALSE)
test_data  <- read_csv("test_data_split.csv",  show_col_types = FALSE)
all_preds  <- read_csv("all_predictions.csv",  show_col_types = FALSE)

# =============================================================================
# 1) Observed Train/Test maps (2 panels, shared legend)
# =============================================================================
mk_sub <- function(df) {
  n <- nrow(df); s <- sum(df$risk == 1, na.rm = TRUE); p <- s / n
  sprintf("n=%d | spills=%d | p(Spill)=%.3f", n, s, p)
}

common_obs_map <- function(df, ttl) {
  ggplot(df %>% arrange(risk), aes(lon, lat)) +
    geom_point(aes(color = as.numeric(risk)), shape = 16, size = 2, alpha = 0.8) +
    scale_color_viridis_c(
      option = "turbo",
      name   = "Observed Outcome",
      limits = c(0, 1),
      breaks = c(0, 1),
      labels = c("No Spill", "Spill"),
      oob    = scales::squish,
      begin  = 0.12,
      end    = 0.9
    ) +
    coord_fixed() +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
    labs(
      title    = ttl,
      subtitle = mk_sub(df),
      x = "Longitude", y = "Latitude"
    )
}

p_train <- common_obs_map(train_data, "Observed — Train")
p_test  <- common_obs_map(test_data,  "Observed — Test")

p_obs_2panel <- (p_train | p_test) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

p_obs_2panel

# =============================================================================
# 2) Predicted maps (4 panels in 2x2, shared legend)
# =============================================================================
pred_long <- all_preds %>%
  pivot_longer(
    cols = c(pred_GLM, pred_SpatialOnly, pred_SpatialPlusZ, pred_SpatialPlusZ_wt),
    names_to  = "model",
    values_to = "pred"
  ) %>%
  mutate(
    model = factor(
      model,
      levels = c("pred_GLM", "pred_SpatialOnly", "pred_SpatialPlusZ", "pred_SpatialPlusZ_wt")
    )
  )

model_labels <- c(
  pred_GLM             = "Covariates Only",
  pred_SpatialOnly     = "Spatial Only",
  pred_SpatialPlusZ    = "Spatial + Covariates",
  pred_SpatialPlusZ_wt = "Spatial + Covariates (Weighted)"
)

p_pred_2x2 <- ggplot(pred_long, aes(lon, lat)) +
  geom_point(aes(color = pred), size = 2) +
  scale_color_viridis_c(
    option = "turbo",
    name   = "Probability",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    oob    = scales::squish,
    begin  = 0.12,
    end    = 0.9
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)) +
  labs(
    title = "Predicted Spill Probability",
    x = "Longitude", y = "Latitude"
  ) +
  facet_wrap(~ model, ncol = 2, labeller = as_labeller(model_labels))

p_pred_2x2
