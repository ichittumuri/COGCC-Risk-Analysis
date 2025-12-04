# =============================================================================
# Compare GLM, Spatial-only, Spatial+Z
# =============================================================================
library(dplyr)
library(readr)
library(ggplot2)
library(scales)
library(tidyr)
library(stringr)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# ---- 1) Load standardized predictions ---------------------------------------
pred_glm  <- read_csv("predictions_glm.csv", show_col_types = FALSE)
pred_so   <- read_csv("predictions_spatial_only.csv", show_col_types = FALSE)
pred_sz   <- read_csv("predictions_spatial_plus_z.csv", show_col_types = FALSE)

df_all <- bind_rows(pred_glm, pred_so, pred_sz) %>%
  mutate(
    model = factor(model, levels = c("GLM", "SpatialOnly", "SpatialPlusZ"))
  )

# Optional: save the combined dataset for downstream use
write_csv(df_all, "predictions_all_models.csv")

# ---- 2) Load metrics and combine --------------------------------------------
met_glm <- read_csv("metrics_glm.csv", show_col_types = FALSE)
met_so  <- read_csv("metrics_spatial_only.csv", show_col_types = FALSE)
met_sz  <- read_csv("metrics_spatial_plus_z.csv", show_col_types = FALSE)

metrics_all <- bind_rows(met_glm, met_so, met_sz) %>%
  mutate(
    model = factor(model, levels = c("GLM", "SpatialOnly", "SpatialPlusZ"))
  ) %>%
  arrange(model)

write_csv(metrics_all, "metrics_all_models.csv")
print(metrics_all)

# ---- 5) Metrics bar chart (optional) ----------------------------------------
metrics_table <- metrics_all %>%
  select(model, accuracy, precision, recall, f1) %>%
  arrange(model)

print(metrics_table)

# ---- Pull balanced accuracy from metrics_all ------------------------------
bal_acc <- metrics_all %>%
  select(model, balanced_accuracy) %>%
  mutate(
    balanced_accuracy = round(balanced_accuracy, 3)  # round for display
  )

# Build custom labels with balanced accuracy appended
# ---- Model labels without balanced accuracy ------------------------------
model_labels <- setNames(
  c("Covariates Only", "Spatial Only", "Spatial + Covariates"),
  levels(metrics_all$model)
)

# ---- Predicted risk maps --------------------------------------------------
df_map <- df_all %>%
  dplyr::group_by(model) %>%
  dplyr::arrange(pred, .by_group = TRUE) %>%
  dplyr::ungroup()

p_maps <- ggplot(df_map, aes(x = lon, y = lat)) +
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
  labs(
    title = "Predicted Spill Probability",
    x = "Longitude",
    y = "Latitude"
  ) +
  facet_wrap(~ model, ncol = 3, labeller = as_labeller(model_labels)) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)  # add box
  )

p_maps

# ---- Observed risk map with border -----------------------------------------
p_obs <- ggplot(df_all %>% arrange(risk), aes(x = lon, y = lat)) +
  geom_point(aes(color = as.numeric(risk)), size = 1.7, alpha = 0.8) +
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
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # add border
  ) +
  labs(
    # title = "Observed Spill Outcomes on Test Data",
    x = "Longitude",
    y = "Latitude"
  )

p_obs

# Round nicely
metrics_table <- metrics_all %>%
  select(model, accuracy, precision, recall, f1) %>%
  arrange(model) %>%
  mutate(
    across(c(accuracy, precision, recall, f1), ~ round(.x, 3))
  )

print(metrics_table)

library(dplyr)
library(ggplot2)
library(patchwork)

# 1) Reset theme
theme_set(theme_grey())

# 2) Helper: clean subtitles
mk_sub <- function(df) {
  n <- nrow(df); s <- sum(df$risk == 1, na.rm = TRUE); p <- s / n
  sprintf("n=%d | spills=%d | p(Spill)=%.3f", n, s, p)
}

sub_train <- mk_sub(train_data)
sub_test  <- mk_sub(test_data)

# 3) Common map with transparency, using same Turbo palette (no outlines)
common_map <- function(df, ttl, subttl) {
  ggplot(df %>% arrange(risk), aes(lon, lat)) +
    geom_point(
      aes(color = as.numeric(risk)), 
      shape = 16,   # solid circle, no outline
      size = 2,
      alpha = 0.8
    ) +
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
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    labs(
      title    = ttl,
      subtitle = subttl,
      x = "Longitude", y = "Latitude"
    )
}

p_train <- common_map(train_data, "Observed — Train", sub_train)
p_test  <- common_map(test_data,  "Observed — Test",  sub_test)

# 4) Combine
combined <- (p_train | p_test) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text()
  )

combined



# ---- Predicted maps (3 facets side by side) ----
ggsave(
  filename = "predicted_maps.png",
  plot = p_maps,
  width = 9, height = 3, dpi = 300
)

# ---- Observed map (single panel) ----
ggsave(
  filename = "observed_map.png",
  plot = p_obs,
  width = 4.5, height = 4, dpi = 300
)

# ---- Combined train/test observed (2 panels) ----
ggsave(
  filename = "observed_train_test.png",
  plot = combined,
  width = 8, height = 4, dpi = 300
)
