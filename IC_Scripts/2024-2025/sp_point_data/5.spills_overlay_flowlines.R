library(leaflet)
library(sf)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/IC_Scripts/2024-2025/sp_point_data")

line_data <- st_read("final_line_data.geojson")
spills <- st_read("updated_spills.geojson")
line_data <- st_transform(line_data, 4326)
spills <- st_transform(spills, 4326)
spills_matched <- spills[st_geometry_type(spills) %in% c("POINT", "MULTIPOINT"), ]

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolylines(
    data = line_data,
    color = "skyblue",
    weight = 2,
    opacity = 0.6,
    group = "Flowlines"
  ) %>%
  addCircleMarkers(
    data = spills,
    color = "firebrick",
    radius = 2,
    stroke = FALSE,
    fillOpacity = 0.8,
    group = "Spills"
  ) %>%
  leaflet::addLegend(
    position = "bottomright",
    colors   = c("skyblue", "firebrick"),
    labels   = c("Flowlines", "Spills"),
    title    = "Line Data"
  ) %>%
  addLayersControl(
    overlayGroups = c("Flowlines", "Spills"),
    options       = layersControlOptions(collapsed = FALSE)
  )