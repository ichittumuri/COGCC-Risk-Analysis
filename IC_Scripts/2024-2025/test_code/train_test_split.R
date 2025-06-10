# 0. install.packages(c("sf","INLA","dplyr","leaflet"))
library(sf)
library(INLA)
library(dplyr)
library(leaflet)

set.seed(42)   # for reproducibility

# 1. Read & pare down
setwd("/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data")
df <- st_read("flowlines_pop_dems.geojson")
df_sf <- df %>%
  select(-Max_Elevation, -Min_Elevation, -operator_number,
         -flowline_id, -location_id, -root_cause)

# 2. Factors
df_sf <- df_sf %>%
  mutate(across(c(status, fluid, material, location_type), factor))

# 3. Centroids & drop NAs (for mesh construction only)
cent         <- st_centroid(df_sf$geometry)
coords       <- st_coordinates(cent)
keep         <- complete.cases(coords)
df_model     <- df_sf[keep, ]

# 4. Transform to WGS84
df_model <- st_transform(df_model, crs = 4326)

# 5. Stratified 80/20 train-test split
pos_idx <- which(df_model$risk == 1)
neg_idx <- which(df_model$risk == 0)

train_pos <- sample(pos_idx, size = floor(0.8 * length(pos_idx)))
train_neg <- sample(neg_idx, size = floor(0.8 * length(neg_idx)))
train_idx <- c(train_pos, train_neg)

df_train <- df_model[train_idx, ]
df_test  <- df_model[-train_idx, ]

# extract centroids for mesh only
coords_train <- st_coordinates(st_centroid(df_train$geometry))
coords_test  <- st_coordinates(st_centroid(df_test$geometry))

# ----------------------------------
# PART A: SPDE‐INLA on TRAIN / PREDICT on TEST
# ----------------------------------

# 6. Build mesh on TRAIN coords
mesh <- inla.mesh.2d(
  loc      = coords_train,
  max.edge = c(0.1, 0.5),
  cutoff   = 0.05
)
spde <- inla.spde2.matern(mesh = mesh)

# 7. Projector matrices
A_train <- inla.spde.make.A(mesh = mesh, loc = coords_train)
A_test  <- inla.spde.make.A(mesh = mesh, loc = coords_test)

# 8. Prepare stacks: TRAIN (y observed) + TEST (y = NA)
fixed_tr <- data.frame(
  intercept              = 1,
  status                 = df_train$status,
  fluid                  = df_train$fluid,
  material               = df_train$material,
  location_type          = df_train$location_type,
  diameter_in            = df_train$diameter_in,
  length_ft              = df_train$length_ft,
  max_operating_pressure = df_train$max_operating_pressure,
  line_age_yr            = df_train$line_age_yr,
  average_pop_density    = df_train$average_pop_density,
  Avg_Elevation          = df_train$Avg_Elevation
)
fixed_te <- fixed_tr[rep(1, nrow(df_test)), ]
rownames(fixed_te) <- NULL

stk_train <- inla.stack(
  data    = list(y = df_train$risk),
  A       = list(1, A_train),
  effects = list(fixed_tr, spatial = 1:spde$n.spde),
  tag     = "train"
)
stk_test <- inla.stack(
  data    = list(y = rep(NA, nrow(df_test))),
  A       = list(1, A_test),
  effects = list(fixed_te, spatial = 1:spde$n.spde),
  tag     = "test"
)
stk_all <- inla.stack(stk_train, stk_test)

# 9. Formula w/ spatial
form_inla <- y ~ -1 +
  intercept + status + fluid + material + location_type +
  diameter_in + length_ft + max_operating_pressure +
  line_age_yr + average_pop_density + Avg_Elevation +
  f(spatial, model = spde)

# 10. Fit INLA on TRAIN, compute predictions for TEST
fit_inla <- inla(
  form_inla,
  family = "binomial",
  data   = inla.stack.data(stk_all),
  control.predictor = list(A = inla.stack.A(stk_all), compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE),
  control.inla      = list(control.vb = list(enable = FALSE))
)

# extract indices
idx_train <- inla.stack.index(stk_all, tag="train")$data
idx_test  <- inla.stack.index(stk_all, tag="test")$data

df_train$inla_fitted <- fit_inla$summary.fitted.values$mean[idx_train]
df_train$inla_class  <- as.integer(df_train$inla_fitted > 0.5)

df_test$inla_pred  <- fit_inla$summary.fitted.values$mean[idx_test]
df_test$inla_class <- as.integer(df_test$inla_pred > 0.5)

# ----------------------------------
# PART B: STANDARD GLM on TRAIN / PREDICT on TEST
# ----------------------------------

glm_fit <- glm(
  risk ~ status + fluid + material + location_type +
    diameter_in + length_ft + max_operating_pressure +
    line_age_yr + average_pop_density + Avg_Elevation,
  data   = df_train,
  family = binomial
)

df_train$glm_fitted <- predict(glm_fit, type = "response")
df_train$glm_class  <- as.integer(df_train$glm_fitted > 0.5)

df_test$glm_pred  <- predict(glm_fit, newdata = df_test, type = "response")
df_test$glm_class <- as.integer(df_test$glm_pred > 0.5)

# ----------------------------------
# PART C: CLASS BALANCE & EVALUATION
# ----------------------------------

# 1. Class balance in train vs. test
train_counts <- table(df_train$risk)
test_counts  <- table(df_test$risk)

train_props  <- prop.table(train_counts)
test_props   <- prop.table(test_counts)

# 2. Balance in the test‐set predictions
inla_pred_counts <- table(df_test$inla_class)
glm_pred_counts  <- table(df_test$glm_class)

inla_pred_props  <- prop.table(inla_pred_counts)
glm_pred_props   <- prop.table(glm_pred_counts)

# 3. Print a nice summary
cat("=== TRUE RISK (train) ===\n")
print(train_counts); cat("→", round(train_props,3), "\n\n")

cat("=== TRUE RISK (test) ===\n")
print(test_counts);  cat("→", round(test_props,3),  "\n\n")

cat("=== INLA PREDICTIONS (test) ===\n")
print(inla_pred_counts); cat("→", round(inla_pred_props,3), "\n\n")

cat("===  GLM PREDICTIONS (test) ===\n")
print(glm_pred_counts);  cat("→", round(glm_pred_props,3),  "\n\n")

# 4. Confusion matrices & model criteria
conf_inla_test <- table(
  Observed  = df_test$risk,
  Predicted = df_test$inla_class
)
acc_inla_test <- mean(df_test$inla_class == df_test$risk)

conf_glm_test <- table(
  Observed  = df_test$risk,
  Predicted = df_test$glm_class
)
acc_glm_test <- mean(df_test$glm_class == df_test$risk)

cat("=== HOLD-OUT TEST RESULTS ===\n\n")
cat("INLA SPDE MODEL:\n")
cat("  DIC:", fit_inla$dic$dic, "\n")
cat("  WAIC:", fit_inla$waic$waic, "\n")
cat("  Test accuracy:", round(acc_inla_test,3), "\n")
print(conf_inla_test)

cat("\nSTANDARD GLM:\n")
cat("  AIC:", round(AIC(glm_fit),1), "\n")
cat("  Test accuracy:", round(acc_glm_test,3), "\n")
print(conf_glm_test)

# ----------------------------------
# PART D: INTERACTIVE MAP OF SPLIT & PREDICTIONS (using original LINESTRINGS)
# ----------------------------------

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  
  # Train (true risk)
  addPolylines(
    data    = df_train,
    color   = ~ifelse(risk       == 1, "red", "blue"),
    weight  = 2, opacity = 0.6,
    label   = ~paste0("Train – obs=", risk),
    group   = "Train: true risk"
  ) %>%
  
  # Test (true risk)
  addPolylines(
    data    = df_test,
    color   = ~ifelse(risk       == 1, "red", "blue"),
    weight  = 2, opacity = 0.6,
    label   = ~paste0("Test  – obs=", risk),
    group   = "Test: true risk"
  ) %>%
  
  # Test INLA predictions
  addPolylines(
    data    = df_test,
    color   = ~ifelse(inla_class == 1, "red", "blue"),
    weight  = 2, opacity = 0.6,
    label   = ~paste0("INLA  – pred=", inla_class),
    group   = "Test: INLA pred"
  ) %>%
  
  # Test GLM predictions
  addPolylines(
    data    = df_test,
    color   = ~ifelse(glm_class  == 1, "red", "blue"),
    weight  = 2, opacity = 0.6,
    label   = ~paste0("GLM   – pred=", glm_class),
    group   = "Test: GLM pred"
  ) %>%
  
  addLayersControl(
    overlayGroups = c(
      "Train: true risk",
      "Test: true risk",
      "Test: INLA pred",
      "Test: GLM pred"
    ),
    options = layersControlOptions(collapsed = FALSE)
  )

# # ----------------------------------
# # PART D: INTERACTIVE MAP OF SPLIT & PREDICTIONS (using centroids)
# # ----------------------------------
# 
# library(leaflet)
# library(dplyr)
# library(sf)
# 
# # 1. Turn each LINESTRING into its centroid POINT
# train_pts <- df_train %>% st_centroid()
# test_pts  <- df_test  %>% st_centroid()
# 
# # 2. Build the map with four overlay groups
# leaflet() %>%
#   addProviderTiles(providers$CartoDB.Positron) %>%
#   
#   # Train (true risk)
#   addCircleMarkers(
#     data        = train_pts,
#     color       = ~ifelse(risk == 1, "red", "blue"),
#     radius      = 5,
#     stroke      = FALSE,
#     fillOpacity = 0.8,
#     label       = ~paste0("Train • obs=", risk),
#     group       = "Train: true risk"
#   ) %>%
#   
#   # Test (true risk)
#   addCircleMarkers(
#     data        = test_pts,
#     color       = ~ifelse(risk == 1, "red", "blue"),
#     radius      = 5,
#     stroke      = FALSE,
#     fillOpacity = 0.8,
#     label       = ~paste0("Test • obs=", risk),
#     group       = "Test: true risk"
#   ) %>%
#   
#   # Test INLA predictions
#   addCircleMarkers(
#     data        = test_pts,
#     color       = ~ifelse(inla_class == 1, "red", "blue"),
#     radius      = 5,
#     stroke      = FALSE,
#     fillOpacity = 0.8,
#     label       = ~paste0("INLA • pred=", inla_class),
#     group       = "Test: INLA pred"
#   ) %>%
#   
#   # Test GLM predictions
#   addCircleMarkers(
#     data        = test_pts,
#     color       = ~ifelse(glm_class == 1, "red", "blue"),
#     radius      = 5,
#     stroke      = FALSE,
#     fillOpacity = 0.8,
#     label       = ~paste0("GLM  • pred=", glm_class),
#     group       = "Test: GLM pred"
#   ) %>%
#   
#   # Layer control
#   addLayersControl(
#     overlayGroups = c(
#       "Train: true risk",
#       "Test: true risk",
#       "Test: INLA pred",
#       "Test: GLM pred"
#     ),
#     options = layersControlOptions(collapsed = FALSE)
#   )
# 
