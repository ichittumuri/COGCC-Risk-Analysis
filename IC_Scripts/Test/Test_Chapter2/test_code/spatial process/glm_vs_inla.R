# 0. install.packages(c("sf","INLA","dplyr","leaflet"))
library(sf)
library(INLA)
library(dplyr)
library(leaflet)

# 1. Read & pare down
setwd("/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data")
df <- st_read("flowlines_pop_dems.geojson")
df_sf <- df %>%
  select(-Max_Elevation, -Min_Elevation, -operator_number,
         -flowline_id, -location_id, -root_cause)

# 2. Factors
df_sf <- df_sf %>%
  mutate(across(c(status, fluid, material, location_type), factor))

# 3. Centroids & drop NAs
cent         <- st_centroid(df_sf$geometry)
coords       <- st_coordinates(cent)
keep         <- complete.cases(coords)
df_model     <- df_sf[keep, ]
coords_model <- coords[keep, ]

# 3b. Transform to WGS84
df_model <- st_transform(df_model, crs = 4326)

# -------------------------------
# PART A: SPDE‐INLA
# -------------------------------

# 4. Build mesh
mesh <- inla.mesh.2d(
  loc      = coords_model,
  max.edge = c(0.1, 0.5),
  cutoff   = 0.05
)

# 5. Matern SPDE
spde <- inla.spde2.matern(mesh = mesh)

# 6. Projector
A <- inla.spde.make.A(mesh = mesh, loc = coords_model)

# 7. Stack
fixed_effects <- data.frame(
  intercept              = 1,
  status                 = df_model$status,
  fluid                  = df_model$fluid,
  material               = df_model$material,
  location_type          = df_model$location_type,
  diameter_in            = df_model$diameter_in,
  length_ft              = df_model$length_ft,
  max_operating_pressure = df_model$max_operating_pressure,
  line_age_yr            = df_model$line_age_yr,
  average_pop_density    = df_model$average_pop_density,
  Avg_Elevation          = df_model$Avg_Elevation
)

stk <- inla.stack(
  data    = list(y = df_model$risk),
  A       = list(1, A),
  effects = list(fixed_effects, spatial = 1:spde$n.spde),
  tag     = "est"
)

# 8. Formula with spatial random field
form_inla <- y ~ -1 +
  intercept + status + fluid + material + location_type +
  diameter_in + length_ft + max_operating_pressure +
  line_age_yr + average_pop_density + Avg_Elevation +
  f(spatial, model = spde)

# 9. Fit INLA
fit_inla <- inla(
  form_inla,
  family = "binomial",
  data   = inla.stack.data(stk),
  control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE),
  control.inla      = list(control.vb = list(enable = FALSE))
)

# 10. INLA predictions & class
idx <- inla.stack.index(stk, tag = "est")$data
df_model$inla_fitted   <- fit_inla$summary.fitted.values$mean[idx]
df_model$inla_class    <- as.integer(df_model$inla_fitted > 0.5)

# 11. INLA confusion & accuracy
orig_counts_inla <- table(Original=df_model$risk)
pred_counts_inla <- table(Predicted=df_model$inla_class)
conf_mat_inla    <- table(Observed=df_model$risk, Predicted=df_model$inla_class)

inla_accuracy <- sum(df_model$inla_class == df_model$risk) / nrow(df_model)

# -------------------------------
# PART B: STANDARD GLM
# -------------------------------

# 12. Fit standard logistic GLM (no spatial)
glm_fit <- glm(
  risk ~ status + fluid + material + location_type +
    diameter_in + length_ft + max_operating_pressure +
    line_age_yr + average_pop_density + Avg_Elevation,
  data   = df_model,
  family = binomial
)

# 13. GLM predictions & class
df_model$glm_fitted <- predict(glm_fit, type = "response")
df_model$glm_class  <- as.integer(df_model$glm_fitted > 0.5)

# 14. GLM confusion & accuracy
orig_counts_glm <- table(Original=df_model$risk)
pred_counts_glm <- table(Predicted=df_model$glm_class)
conf_mat_glm    <- table(Observed=df_model$risk, Predicted=df_model$glm_class)

glm_accuracy <- sum(df_model$glm_class == df_model$risk) / nrow(df_model)

# -------------------------------
# PART C: MODEL COMPARISON OUTPUT
# -------------------------------

cat("===== INLA SPDE MODEL =====\n")
cat("DIC:", fit_inla$dic$dic, "\n")
cat("WAIC:", fit_inla$waic$waic, "\n")
cat("Accuracy:", round(inla_accuracy, 3), "\n")
print(conf_mat_inla)

cat("\n===== STANDARD GLM =====\n")
cat("AIC:", AIC(glm_fit), "\n")
cat("Accuracy:", round(glm_accuracy, 3), "\n")
print(conf_mat_glm)

