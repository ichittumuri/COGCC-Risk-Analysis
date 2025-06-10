# 1) Setup
setwd("/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data")
library(sf)
library(fields)
library(dplyr)
library(leaflet)

# 2) Source functions
source("./Doug_R_Functions/logisticSmoother.R")
source("./Doug_R_Functions/logisticRegressionSimple.R")

# 3) Load data and extract coords
flowlines <- st_read("final_dataset.geojson") %>% 
  st_transform(4326) %>% 
  mutate(X = st_coordinates(geometry)[,1], Y = st_coordinates(geometry)[,2]) %>% 
  st_drop_geometry()

# 4) Prepare modeling table
drop.cols <- c("unique_id", "operator_number", "flowline_id", "location_id",
               "construct_date", "spill_date", "root_cause")
id_coords_risk <- flowlines %>% select(unique_id, X, Y, risk)
model_df <- id_coords_risk %>% filter(!is.na(risk), !is.na(X), !is.na(Y))
coords_matrix <- as.matrix(model_df %>% select(X, Y))
y_binary <- model_df$risk

# 5) Run spatial IRLS smoother
lambda_chosen <- 1e-3
look <- logisticSmoother(s = coords_matrix, y = y_binary, lambda = lambda_chosen)
nu_hat <- look$fitted.values
p_hat <- exp(nu_hat) / (1 + exp(nu_hat))

# 6) Attach predictions and export
predictions_df <- data.frame(
  unique_id = model_df$unique_id,
  X = coords_matrix[,1],
  Y = coords_matrix[,2],
  p_hat = p_hat
)
predictions_sf <- st_as_sf(predictions_df, coords = c("X", "Y"), crs = 4326)
st_write(predictions_sf, "flowline_spatial_logistic_predictions.geojson", delete_dsn = TRUE)

# 7) Quick Leaflet preview
leaflet(predictions_sf) %>% 
  addProviderTiles(providers$CartoDB.Positron) %>% 
  addCircleMarkers(radius = 4, color = ~colorNumeric("RdYlBu", p_hat, reverse = TRUE)(p_hat),
                   stroke = FALSE, fillOpacity = 0.8,
                   popup = ~paste("ID:", unique_id, "<br> P(risk):", round(p_hat, 3))) %>% 
  addLegend("bottomright", pal = colorNumeric("RdYlBu", predictions_sf$p_hat, reverse = TRUE),
            values = predictions_sf$p_hat, title = "P(risk=1)")
