# compare_models_ci.R
# Reads model_out_glm_lasso.csv, model_out_spatial_only.csv, model_out_spatial_covs.csv
# and plots a single CI comparison for all three models.

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# -------------------------
# 1) Paths (edit if needed)
# -------------------------


setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

library(readr)
library(dplyr)
library(ggplot2)

# Read CSVs from your export step
glm_df <- read_csv("model_out_glm_lasso.csv", show_col_types = FALSE)
sp_df  <- read_csv("model_out_spatial_only.csv", show_col_types = FALSE)
spc_df <- read_csv("model_out_spatial_covs.csv", show_col_types = FALSE)

# Combine and order within each model
df_all <- bind_rows(glm_df, sp_df, spc_df) %>%
  group_by(model) %>%
  arrange(pred) %>%
  mutate(idx = row_number()) %>%
  ungroup()

# 3x1 facet plot, matching your style
ggplot(df_all, aes(x = idx, y = pred)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "red", alpha = 0.3) +
  geom_line(color = "black", size = 1) +
  facet_wrap(~model, ncol = 1, scales = "fixed") +
  labs(
    x = "Observation (sorted by predicted risk)",
    y = "Predicted Spill Risk",
    title = "Predicted Spill Risk ± 95% CI"
  ) +
  theme_minimal(base_size = 13)






# library(readr)
# library(dplyr)
# library(ggplot2)
# # Read all models
# glm_df <- read_csv("model_out_glm_lasso.csv", show_col_types = FALSE)
# sp_df  <- read_csv("model_out_spatial_only.csv", show_col_types = FALSE)
# spc_df <- read_csv("model_out_spatial_covs.csv", show_col_types = FALSE)
# 
# # Force consistent naming before plotting
# df_all <- df_all %>%
#   mutate(model = case_when(
#     grepl("glm", model, ignore.case = TRUE) ~ "GLM (LASSO)",
#     grepl("spatial\\s*only", model, ignore.case = TRUE) ~ "Spatial Only",
#     grepl("spatial.*cov", model, ignore.case = TRUE) ~ "Spatial + Covariates",
#     TRUE ~ model
#   ))
# 
# # Now define colors with exact matches
# model_colors <- c(
#   "GLM (LASSO)" = "#d73027",         # red
#   "Spatial Only" = "#1a9850",        # green
#   "Spatial + Covariates" = "#4575b4" # blue
# )
# 
# ggplot(df_all, aes(x = idx, y = pred, color = model, fill = model)) +
#   geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, colour = NA) +
#   geom_line(size = 1) +
#   facet_wrap(~model, ncol = 1, scales = "fixed") +
#   scale_color_manual(values = model_colors) +
#   scale_fill_manual(values = model_colors) +
#   labs(
#     x = "Observation (sorted by predicted risk)",
#     y = "Predicted Spill Risk",
#     title = "Predicted Spill Risk ± 95% CI"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     legend.position = "none",
#     panel.grid.minor = element_blank()
#   )
