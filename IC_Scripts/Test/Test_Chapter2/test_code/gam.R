# Load required packages
library(sf)
library(dplyr)
library(mgcv)

# Set working directory
setwd("/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data")

# Load the spatial dataset
flowlines <- st_read("final_dataset.geojson")

# Drop unnecessary metadata columns
drop.cols <- c("operator_number", "flowline_id", "location_id", 
               "construct_date", "spill_date", "root_cause")
flowlines <- flowlines %>% select(-any_of(drop.cols))

# Extract coordinates and drop geometry
flowlines <- flowlines %>%
  mutate(
    longitude = st_coordinates(.)[, 1],
    latitude  = st_coordinates(.)[, 2]
  ) %>% 
  st_drop_geometry()

# Convert character columns to factors
flowlines <- flowlines %>%
  mutate(across(
    c(status, flowline_action, location_type, fluid, material),
    as.factor
  )) %>%
  filter(complete.cases(.))  # Remove missing values

# Train/test split
set.seed(42)
train_idx <- sample(seq_len(nrow(flowlines)), 0.8 * nrow(flowlines))
train_data <- flowlines[train_idx, ]
test_data  <- flowlines[-train_idx, ]

# Fit spatial logistic regression model (GAM)
gam_model <- gam(risk ~ s(longitude, latitude), 
                 family = binomial(link = "logit"),
                 data = train_data)

# Predict logit values + standard errors on test data
pred <- predict(gam_model, newdata = test_data, type = "link", se.fit = TRUE)

# Convert logits to probabilities with 95% confidence intervals
test_data <- test_data %>%
  mutate(
    prob_risk  = plogis(pred$fit),
    prob_lower = plogis(pred$fit - 2 * pred$se.fit),
    prob_upper = plogis(pred$fit + 2 * pred$se.fit),
    ci_width   = prob_upper - prob_lower
  ) %>%
  mutate(across(c(prob_risk, prob_lower, prob_upper), ~pmin(pmax(., 0.001), 0.999)))

# Convert to sf and export to GeoJSON
test_sf <- st_as_sf(test_data, coords = c("longitude", "latitude"), crs = 4326)
st_write(test_sf, "predicted_risk_with_layers.geojson", delete_dsn = TRUE)

library(leaflet)
library(sf)

# Load GeoJSON with predictions
pred_sf <- st_read("predicted_risk_with_layers.geojson")

# Define color palettes for each layer
pal_mean <- colorNumeric("YlOrRd", pred_sf$prob_risk)
pal_upper <- colorNumeric("PuBu", pred_sf$prob_upper)
pal_ciwidth <- colorNumeric("Greys", pred_sf$ci_width)

# Create leaflet map with 3 overlay layers and a switcher
leaflet(pred_sf) %>%
  addProviderTiles(providers$CartoDB.Positron, group = "Base Map") %>%
  
  # Layer 1: Mean predicted risk
  addCircleMarkers(
    radius = 3,
    stroke = FALSE,
    fillOpacity = 0.8,
    color = ~pal_mean(prob_risk),
    group = "Predicted Risk",
    popup = ~paste0("Predicted: ", round(prob_risk, 3),
                    "<br>CI: [", round(prob_lower, 3), ", ", round(prob_upper, 3), "]")
  ) %>%
  
  # Layer 2: Upper 95% CI
  addCircleMarkers(
    radius = 3,
    stroke = FALSE,
    fillOpacity = 0.8,
    color = ~pal_upper(prob_upper),
    group = "Upper 95% CI",
    popup = ~paste0("Upper Bound: ", round(prob_upper, 3))
  ) %>%
  
  # Layer 3: CI Width (uncertainty)
  addCircleMarkers(
    radius = 3,
    stroke = FALSE,
    fillOpacity = 0.8,
    color = ~pal_ciwidth(ci_width),
    group = "CI Width (Uncertainty)",
    popup = ~paste0("CI Width: ", round(ci_width, 3))
  ) %>%
  
  # Add legends
  addLegend("bottomright", pal = pal_mean, values = ~prob_risk,
            title = "Predicted Risk", group = "Predicted Risk") %>%
  addLegend("bottomright", pal = pal_upper, values = ~prob_upper,
            title = "Upper 95% CI", group = "Upper 95% CI") %>%
  addLegend("bottomright", pal = pal_ciwidth, values = ~ci_width,
            title = "Uncertainty Range", group = "CI Width (Uncertainty)") %>%
  
  # Layer control panel
  addLayersControl(
    baseGroups = c("Base Map"),
    overlayGroups = c("Predicted Risk", "Upper 95% CI", "CI Width (Uncertainty)"),
    options = layersControlOptions(collapsed = FALSE)
  )