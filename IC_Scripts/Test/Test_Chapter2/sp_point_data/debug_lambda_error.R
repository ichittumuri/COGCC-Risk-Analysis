#───────────────────────────────────────────────────────────────────────────────
# Minimal Example to Isolate Spatial Process Error at λ ≈ 4.64
#───────────────────────────────────────────────────────────────────────────────

library(readr)
library(dplyr)
library(fields)
library(MASS)


# model_df <- data.frame(
#   lon  = coord_model[, 1],
#   lat  = coord_model[, 2],
#   risk = risk_model
# )


# Set your working directory and source functions
setwd("~/Desktop/MINES/COGCC-Risk-Analysis/IC_Scripts/2024-2025/sp_point_data")
source("logisticSmoother.R")

# Load minimal model dataset
model_df <- read_csv("model_df.csv")
coord_model <- as.matrix(model_df[, c("lon", "lat")])
risk_model  <- model_df$risk

# Load full point data for Z matrix
usable_sf <- st_read("usable_point_data.geojson")

# Match the same spatial window as model_df
df <- usable_sf %>%
  st_drop_geometry() %>%
  select(-unique_id) %>%
  mutate(
    lon = as.numeric(sub("^\\((.*),.*$", "\\1", coords)),
    lat = as.numeric(sub("^.*,\\s*(.*)\\)$",  "\\1", coords))
  ) %>%
  filter(
    lon >= -105.85,
    lon <= -104.79 + 0.2,
    lat >=  39.42,
    lat <=  40.68
  )

# Ensure same ordering and row count (must match model_df exactly)
df_small <- df[1:nrow(model_df), ]  # truncate if needed

# Construct Z matrix
to.factor <- c("status", "flowline_action", "location_type", "fluid", "material")
df_mutated <- df_small %>%
  mutate(across(all_of(to.factor), as.factor)) %>%
  as.data.frame()

X <- model.matrix(~ status + flowline_action + location_type + fluid + material +
                    diameter_in + length_ft + line_age_yr, data = df_mutated)[, -1]

Z <- X[, c(
  "statusNew_Construction",
  "flowline_actionRealignment",
  "flowline_actionRegistration",
  "location_typeWell_Site",
  "fluidMultiphase",
  "fluidOther",
  "materialSteel",
  "diameter_in",
  "length_ft",
  "line_age_yr"
)]

# Use a linear model for starting values
glm_cov <- glm(model_df$risk ~ ., data = as.data.frame(cbind(Z, risk = model_df$risk)),
               family = binomial(link = "logit"))
init_nu <- predict(glm_cov, type = "link")

# Try calling logisticSmoother with the known failing λ
lambda <- 4.64159  # Approximate lambda that caused failure

fit_test <- logisticSmoother(
  s      = coord_model,
  y      = risk_model,
  Z      = Z,
  lambda = lambda,
  nuOld  = init_nu
)
