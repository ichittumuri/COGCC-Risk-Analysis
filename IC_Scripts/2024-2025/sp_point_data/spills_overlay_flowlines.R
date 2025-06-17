library(leaflet)
library(sf)

setwd("/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data")
 
# Matched (processed) files
flowlines_matched <- st_read("full_length_flowlines.geojson")
spills_matched <- st_read("updated_spills.geojson")

# Transform to WGS84
flowlines_matched <- st_transform(flowlines_matched, 4326)
spills_matched <- st_transform(spills_matched, 4326)

# Filter to only points
spills_matched <- spills_matched[st_geometry_type(spills_matched) %in% c("POINT", "MULTIPOINT"), ]

# Plot matched data
leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolylines(
    data = flowlines_matched,
    color = "blue",
    weight = 2,
    opacity = 0.6,
    group = "Flowlines"
  ) %>%
  addCircleMarkers(
    data = spills_matched,
    color = "red",
    radius = 2,
    stroke = FALSE,
    fillOpacity = 0.8,
    group = "Spills"
  ) %>%
  addLegend(
    position = "bottomright",
    colors = c("blue", "red"),
    labels = c("Flowlines", "Spills"),
    title = "Matched Data"
  )
