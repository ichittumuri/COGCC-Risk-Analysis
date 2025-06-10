# 0. install.packages(c("sf","fields","dplyr","leaflet"))
library(sf)
library(fields)
library(dplyr)
library(leaflet)

# 1. Read & pare down
setwd("/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data")
df <- st_read("flowlines_pop_dems.geojson")
df_sf <- df %>%
  select(-Max_Elevation, -Min_Elevation, -operator_number,
         -flowline_id, -location_id, -root_cause) %>%
  mutate(across(c(status, fluid, material, location_type), factor))

# 2. Compute centroids & drop missing
cent       <- st_centroid(df_sf$geometry)
coords_all <- st_coordinates(cent)
keep       <- complete.cases(coords_all)
df_model   <- df_sf[keep, ]
coords_mat <- coords_all[keep, ]  # numeric matrix

# 3. Non-spatial logistic → residuals
glm0 <- glm(risk ~ status + fluid + material + location_type +
              diameter_in + length_ft + max_operating_pressure +
              line_age_yr + average_pop_density + Avg_Elevation,
            data = df_model, family = binomial)
df_model$resid <- residuals(glm0, type = "response")

# 4. Collapse duplicates by rounding to 4 decimals
df_unique <- df_model %>%
  st_drop_geometry() %>%
  mutate(
    lon = round(coords_mat[,1], 4),
    lat = round(coords_mat[,2], 4)
  ) %>%
  group_by(lon, lat) %>%
  summarize(resid = mean(resid), .groups = "drop")

# 5. Prepare for fields::spatialProcess()
coords_unique <- as.matrix(df_unique[, c("lon","lat")])
resid_unique  <- df_unique$resid

# 6. Remove any remaining exact duplicates
dup_idx <- duplicated(coords_unique)
if (any(dup_idx)) {
  coords_unique <- coords_unique[!dup_idx, , drop = FALSE]
  resid_unique  <- resid_unique[!dup_idx]
}

# 7. Fit the spatial GP to these unique residuals
sp_proc <- spatialProcess(coords_unique, resid_unique)

# 8. Extract the spatial effect & SE into df_unique
df_unique <- df_unique %>%
  mutate(
    sp_effect    = sp_proc$Z,
    sp_effect_se = sp_proc$Z.se,
    key          = paste(lon, lat, sep = "_")
  )

# 9. Build matching key in df_model and join
df_model <- df_model %>%
  st_drop_geometry() %>%
  mutate(
    lon = round(coords_mat[,1], 4),
    lat = round(coords_mat[,2], 4),
    key = paste(lon, lat, sep = "_")
  ) %>%
  left_join(
    df_unique %>% select(key, sp_effect, sp_effect_se),
    by = "key"
  )

# check that join succeeded
stopifnot(all(c("sp_effect","sp_effect_se") %in% names(df_model)))

# 10. Re-fit logistic including the GP term
glm_fp <- glm(
  risk ~ status + fluid + material + location_type +
    diameter_in + length_ft + max_operating_pressure +
    line_age_yr + average_pop_density + Avg_Elevation +
    sp_effect,
  data   = df_model,
  family = binomial
)

# 11. Predict & assert original 1’s
df_model <- df_model %>%
  mutate(
    pred_risk_fp  = predict(glm_fp, type = "response"),
    pred_class_fp = as.integer(pred_risk_fp > 0.5),
    pred_final_fp = ifelse(risk == 1, 1L, pred_class_fp)
  )

# 12. Confusion & counts
print(table(Original   = df_model$risk))
print(table(WithFields = df_model$pred_final_fp))
print(table(Observed   = df_model$risk,
            Predicted  = df_model$pred_final_fp))

# 13. Rough uncertainty (CI width on probability scale)
df_model <- df_model %>%
  mutate(
    lp_mean     = predict(glm_fp, type = "link"),
    lp_se_sp    = sp_effect_se * coef(glm_fp)["sp_effect"],
    lp_lower95  = lp_mean - 1.96 * lp_se_sp,
    lp_upper95  = lp_mean + 1.96 * lp_se_sp,
    risk_lower  = plogis(lp_lower95),
    risk_upper  = plogis(lp_upper95),
    ci_width_fp = risk_upper - risk_lower
  )
print(summary(df_model$sp_effect_se))
print(summary(df_model$ci_width_fp))

# 14. Leaflet maps

# 14.1 Original risk
leaflet(df_sf) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolylines(
    color   = ~ifelse(risk == 1, "red", "blue"),
    weight  = 2, opacity = 0.7
  ) %>%
  addLegend(
    "bottomright",
    colors = c("blue","red"),
    labels = c("No Spill","Spill"),
    title  = "Original Risk"
  )

# 14.2 Predicted risk (fields)
leaflet(df_model) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolylines(
    color   = ~ifelse(pred_final_fp == 1, "red", "blue"),
    weight  = 2, opacity = 0.7
  ) %>%
  addLegend(
    "bottomright",
    colors = c("blue","red"),
    labels = c("No Spill","Spill"),
    title  = "Predicted Risk (fields)"
  )

# 14.3 Match vs. Mismatch
df_model <- df_model %>%
  mutate(match_fp = ifelse(pred_class_fp == risk, "Match", "Mismatch"))

leaflet(df_model) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolylines(
    color   = ~ifelse(match_fp == "Match", "green", "orange"),
    weight  = 2, opacity = 0.7
  ) %>%
  addLegend(
    "bottomright",
    colors = c("green","orange"),
    labels = c("Match","Mismatch"),
    title  = "Model vs Observed (fields)"
  )

# 14.4 Uncertainty (CI width)
pal_fp <- colorNumeric("YlOrRd", domain = df_model$ci_width_fp)

leaflet(df_model) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolylines(
    color   = ~pal_fp(ci_width_fp),
    weight  = 2, opacity = 0.8
  ) %>%
  addLegend(
    "bottomright",
    pal    = pal_fp,
    values = ~ci_width_fp,
    title  = "CI Width (fields)"
  )

