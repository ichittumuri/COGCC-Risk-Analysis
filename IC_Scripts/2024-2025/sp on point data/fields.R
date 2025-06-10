# 1) Prep -------------------------------------------------------------------

# install.packages(c("sf","fields","pROC"))  # if you haven’t already
setwd("~/Desktop/MINES/COGCC-Risk-Analysis/IC_Scripts/2024-2025/fields_example")
library(sf)              # for reading GeoJSON and handling geometry
library(fields)          # for spatialProcess()
library(pROC)            # to find the best ROC threshold
library(dplyr)

# 2) Source Doug’s functions :contentReference[oaicite:0]{index=0} ------------------
source("logisticSmoother.R")        # Doug’s spatial smoother :contentReference[oaicite:1]{index=1}:contentReference[oaicite:2]{index=2}
source("logisticRegressionSimple.R")# basic IWLS logistic (if you need it) :contentReference[oaicite:3]{index=3}:contentReference[oaicite:4]{index=4}

# 3) Load & clean your data -------------------------------------------------
flowlines <- st_read("final_dataset.geojson")

# drop the ID/date columns you don’t want
drop.cols <- c("unique_id","operator_number","flowline_id","location_id",
               "construct_date","spill_date","root_cause")
flowlines <- flowlines[, setdiff(names(flowlines), drop.cols)]

# factor your categoricals
to.factor <- c("status","flowline_action","location_type","fluid","material")
flowlines <- flowlines %>% mutate(across(all_of(to.factor), as.factor))

# 4) Extract planar coordinates ---------------------------------------------
#  – choose a projection (here Web Mercator) so distances are in meters
flow_proj <- st_transform(flowlines, 3857)
#  – take each line’s midpoint as its “location”
centroids <- st_centroid(flow_proj$geometry)
coords    <- st_coordinates(centroids)           # n×2 matrix of (X,Y)

# get your response vector
y_all <- flow_proj$risk

# remove any rows with missing coords or missing y
valid <- which(!is.na(coords[,1]) & !is.na(y_all))
s_all <- coords[valid, ]
y      <- y_all[valid]

# 5) Choose λ by profiling (optional) ---------------------------------------
lambda_grid <- 10^seq(-4, -1, length=10)
profile_ll  <- numeric(length(lambda_grid))

for(i in seq_along(lambda_grid)){
  fit_i <- logisticSmoother(s_all, y, lambda=lambda_grid[i])
  profile_ll[i] <- fit_i$summary["lnProfileLike.FULL"]
}
λ_best <- lambda_grid[which.max(profile_ll)]
cat("Best λ:", λ_best, "\n")

# 6) Bagging loop (1 000× 90/10 splits) -------------------------------------
set.seed(2025)
n       <- nrow(s_all)
pos_idx <- which(y==1)
neg_idx <- which(y==0)

# matrix to hold all prob predictions
prob_mat <- matrix(NA_real_, nrow=n, ncol=1000)

for(b in 1:1000){
  # stratified 90% train
  train_idx <- c(
    sample(pos_idx, length(pos_idx)*0.9),
    sample(neg_idx, length(neg_idx)*0.9)
  )
  # fit spatial smoother on train
  fit_b <- logisticSmoother(s_all[train_idx,], y[train_idx], lambda=λ_best)
  # predict log-odds at all points, then convert to probability
  nu_pred    <- predict(fit_b, s_all)
  prob_mat[,b] <- exp(nu_pred)/(1+exp(nu_pred))
}

# average across bags
prob_avg <- rowMeans(prob_mat)

# 7) Choose cutoff via ROC --------------------------------------------------
roc_obj    <- roc(y, prob_avg)
cutoff     <- coords(roc_obj, "best", ret="threshold")
cat("ROC-optimal threshold:", cutoff, "\n")

# 8) Write results back to GeoJSON ------------------------------------------
flow_proj$prob_risk <- NA_real_
flow_proj$prob_risk[valid] <- prob_avg

flow_proj$pred_class <- as.integer(flow_proj$prob_risk > cutoff)

st_write(flow_proj,
         "final_dataset_with_spatial_probs.geojson",
         delete_dsn = TRUE)