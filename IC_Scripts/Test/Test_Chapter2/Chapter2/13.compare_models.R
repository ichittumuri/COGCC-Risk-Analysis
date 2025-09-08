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

# ---- Predicted risk maps (GLM, Spatial-only, Spatial+Z) ---------------------
# Order within each model so high-risk points draw last (on top)
df_map <- df_all %>%
  group_by(model) %>%
  arrange(pred, .by_group = TRUE) %>%
  ungroup()

p_maps <- ggplot(df_map, aes(x = lon, y = lat)) +
  geom_point(aes(color = pred), size = 2) +
  scale_color_gradient(
    low    = "skyblue",
    high   = "firebrick",
    name   = "Predicted Risk",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    oob    = scales::squish
  ) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "Predicted Spill Risk by Model",
    x = "Longitude",
    y = "Latitude"
  ) +
  facet_wrap(~ model, ncol = 3)  # use ncol = 1 for a 3×1 layout

ggsave("plot_predicted_maps_by_model.png", p_maps, width = 12, height = 4.5, dpi = 300)
p_maps

