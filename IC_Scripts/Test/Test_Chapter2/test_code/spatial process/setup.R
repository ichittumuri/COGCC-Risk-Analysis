# ====================================================================
# R Script: COGCC Risk Analysis Pipeline (with Confusion Matrix)
# ====================================================================

# 1. Setup ----------------------------------------------------------------

# Set your working directory
setwd("/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data")

# Load required libraries
library(sf)        # spatial data I/O and transforms
library(dplyr)     # data manipulation
library(leaflet)   # interactive mapping
library(spdep)     # spatial weights & autocorrelation

# 2. Import & Initial Cleaning -------------------------------------------

# Read GeoJSON and force WGS84 to avoid datum mismatches
df <- st_read("flowlines_pop_dems.geojson") %>%
  st_transform(crs = 4326)

# Replace missing root_cause with "None"
df$root_cause[is.na(df$root_cause)] <- "None"

# Drop columns we won't use
df <- df %>%
  select(-Max_Elevation, -Min_Elevation,
         -operator_number, -flowline_id,
         -location_id, -root_cause)

# 3. Factor Encoding ------------------------------------------------------

df_encoded <- df %>%
  mutate(
    status        = factor(status),
    fluid         = factor(fluid),
    material      = factor(material),
    location_type = factor(location_type)
  )

# 4. Exploratory Mapping -------------------------------------------------

leaflet(df_encoded) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolylines(color = ~ifelse(risk == 1, "red", "blue"),
               weight = 2, opacity = 0.7) %>%
  addLegend("bottomright",
            colors = c("blue","red"),
            labels = c("No Spill","Spill"),
            title  = "Risk Level")

# 5. Preprocessing for Modeling ------------------------------------------

# Compute centroids and extract lon/lat
centroids    <- st_centroid(df_encoded$geometry)
coords       <- st_coordinates(centroids)
df_clean     <- df_encoded %>%
  mutate(lon = coords[,1],
         lat = coords[,2]) %>%
  filter(!is.na(lon) & !is.na(lat))

# Remove any exact duplicate points to avoid knearneigh errors
df_clean <- df_clean %>%
  distinct(lon, lat, .keep_all = TRUE)

# Final modeling dataframe
df_final <- df_clean %>%
  select(risk,
         status, fluid, material, location_type,
         diameter_in, length_ft, max_operating_pressure,
         line_age_yr, average_pop_density, Avg_Elevation,
         lon, lat)

# 6. Logistic Regression --------------------------------------------------

glm_model <- glm(risk ~ status + fluid + material + location_type +
                   diameter_in + max_operating_pressure + line_age_yr +
                   average_pop_density + Avg_Elevation,
                 data    = df_final,
                 family  = binomial,
                 control = glm.control(epsilon = 1e-8, maxit = 50))

# Review model output
summary(glm_model)

# 7. Confusion Matrix -----------------------------------------------------

# Predict probabilities and classes (threshold = 0.5)
df_final$pred_prob  <- predict(glm_model, type = "response")
df_final$pred_class <- ifelse(df_final$pred_prob > 0.5, 1, 0)

# Build and print confusion matrix
conf_mat <- table(Observed = df_final$risk, Predicted = df_final$pred_class)
print(conf_mat)

# 8. Spatial Autocorrelation of Residuals ---------------------------------

# 8a. Deviance residuals
df_final$residuals <- residuals(glm_model, type = "deviance")

# 8b. Jitter coordinates slightly to break any remaining exact ties
set.seed(42)
coords_mat <- cbind(df_final$lon, df_final$lat)
coords_jit <- apply(coords_mat, 2, function(x) jitter(x, factor = 1e-8))

# 8c. Build 5‐NN neighbor list on jittered coords
knn_nb <- knearneigh(coords_jit, k = 5)
nb     <- knn2nb(knn_nb, sym = TRUE)

# 8d. Convert to spatial weights, allowing isolates
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# 8e. Moran's I test
moran_test <- moran.test(df_final$residuals, lw, zero.policy = TRUE)
print(moran_test)

# 9. Feature Selection via Stepwise AIC/BIC -------------------------------

# Define full and null models
full_model <- glm(risk ~ status + fluid + material + location_type +
                    diameter_in + length_ft + max_operating_pressure +
                    line_age_yr + average_pop_density + Avg_Elevation,
                  data   = df_final,
                  family = binomial)

null_model <- glm(risk ~ 1, data = df_final, family = binomial)

# Run stepwise selection (BIC: k = log(n))
n_obs          <- nrow(df_final)
stepwise_model <- step(null_model,
                       scope     = list(lower = null_model, upper = full_model),
                       direction = "both",
                       k         = log(n_obs))

# Review selected model
summary(stepwise_model)


