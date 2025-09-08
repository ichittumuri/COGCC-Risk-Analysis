# 0. install.packages(c("sf","INLA","dplyr","leaflet"))
library(sf)
library(INLA)
library(dplyr)
library(leaflet)

# 1. Read in and pare down your data
setwd("/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data")
df <- st_read("flowlines_pop_dems.geojson")
df_sf <- df %>%
  select(-Max_Elevation, -Min_Elevation, -operator_number,
         -flowline_id, -location_id, -root_cause)

# 2. Convert categoricals to factors
df_sf <- df_sf %>%
  mutate(
    status        = factor(status),
    fluid         = factor(fluid),
    material      = factor(material),
    location_type = factor(location_type)
  )

# 3. Extract centroids and drop any NAs
cent         <- st_centroid(df_sf$geometry)
coords       <- st_coordinates(cent)
keep         <- !is.na(coords[,1]) & !is.na(coords[,2])
df_model     <- df_sf[keep, ]
coords_model <- coords[keep, ]

# 3b. Transform both sf objects to WGS84 for Leaflet
df_sf    <- st_transform(df_sf,    crs = 4326)
df_model <- st_transform(df_model, crs = 4326)

# 4. Build an SPDE mesh
mesh <- inla.mesh.2d(
  loc      = coords_model,
  max.edge = c(0.1, 0.5),
  cutoff   = 0.05
)

# 5. Define the Matern SPDE model
spde <- inla.spde2.matern(mesh = mesh)

# 6. Projector matrix
A <- inla.spde.make.A(mesh = mesh, loc = coords_model)

# 7. Stack fixed + spatial effects
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

# 8. Specify formula (no intercept in formula)
form <- y ~ -1 +
  intercept +
  status +
  fluid +
  material +
  location_type +
  diameter_in +
  length_ft +
  max_operating_pressure +
  line_age_yr +
  average_pop_density +
  Avg_Elevation +
  f(spatial, model = spde)

# 9. Fit the INLA model with VB‐correction disabled
fit <- inla(
  form,
  family = "binomial",
  data   = inla.stack.data(stk),
  control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE),
  control.inla      = list(control.vb = list(enable = FALSE))
)

# 10. Pull out the fitted probabilities
idx <- inla.stack.index(stk, tag = "est")$data
df_model$pred_risk <- fit$summary.fitted.values$mean[idx]

# 11. Post‐process predictions: apply 0.5 cutoff only
df_model <- df_model %>%
  mutate(
    pred_class = as.integer(pred_risk > 0.5)
  )

# 12. Quantify original vs. thresholded predictions
orig_counts <- table(Original = df_model$risk)
pred_counts <- table(Predicted = df_model$pred_class)
conf_mat    <- table(Observed = df_model$risk, Predicted = df_model$pred_class)

print(orig_counts)
print(pred_counts)
print(conf_mat)

# 13. Map the original 0/1 risk
leaflet(df_sf) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolylines(
    color   = ~ifelse(risk == 1, "red", "blue"),
    weight  = 2, opacity = 0.7
  ) %>%
  addLegend(
    "bottomright",
    colors = c("blue", "red"),
    labels = c("No Spill", "Spill"),
    title  = "Original Risk"
  )

# 14. Map the thresholded predictions
leaflet(df_model) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolylines(
    color   = ~ifelse(pred_class == 1, "red", "blue"),
    weight  = 2, opacity = 0.7
  ) %>%
  addLegend(
    "bottomright",
    colors = c("blue", "red"),
    labels = c("No Spill", "Spill"),
    title  = "Predicted Risk (0.5 cutoff)"
  )

# 15. Match vs. mismatch map
df_model <- df_model %>%
  mutate(match = ifelse(pred_class == risk, "Match", "Mismatch"))

leaflet(df_model) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolylines(
    color   = ~ifelse(match == "Match", "green", "orange"),
    weight  = 2, opacity = 0.7
  ) %>%
  addLegend(
    "bottomright",
    colors = c("green", "orange"),
    labels = c("Match", "Mismatch"),
    title  = "Model vs Observed"
  )

# 16. Uncertainty quantification (unchanged)
s <- fit$summary.fitted.values

df_model <- df_model %>%
  mutate(
    sd_risk  = s$sd[idx],
    lower95  = s$`0.025quant`[idx],
    upper95  = s$`0.975quant`[idx],
    ci_width = upper95 - lower95
  )

print(summary(df_model$sd_risk))
print(summary(df_model$ci_width))

print(table(
  Uncertainty = cut(
    df_model$sd_risk,
    breaks = c(0, 0.05, 0.1, 0.2, 1),
    labels = c("<0.05","0.05–0.1","0.1–0.2",">0.2"),
    include.lowest = TRUE
  )
))

# 17. Map the posterior SD
pal_sd <- colorNumeric("YlOrRd", domain = df_model$sd_risk)

leaflet(df_model) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolylines(
    color   = ~pal_sd(sd_risk),
    weight  = 2, opacity = 0.8
  ) %>%
  addLegend(
    "bottomright",
    pal    = pal_sd,
    values = ~sd_risk,
    title  = "Posterior SD\nof Predicted Risk"
  )
