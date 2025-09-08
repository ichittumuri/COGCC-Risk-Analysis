# 0. install.packages(c("sf","spBayes","dplyr","leaflet"))
library(sf)
library(spBayes)
library(dplyr)
library(leaflet)

# 1. Read in & pare down data
setwd("/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data")
df <- st_read("flowlines_pop_dems.geojson")
df_sf <- df %>%
  select(-Max_Elevation, -Min_Elevation, -operator_number,
         -flowline_id, -location_id, -root_cause) %>%
  mutate(across(c(status, fluid, material, location_type), factor))

# 2. Extract centroids & drop any NAs
cent       <- st_centroid(df_sf$geometry)
coords      <- st_coordinates(cent)
keep        <- complete.cases(coords)
df_model    <- df_sf[keep, ]
coords_mat  <- coords[keep, ]

# 3. Scale coords to kilometers (avoids huge covariance entries)
coords_km   <- coords_mat / 1000

# 4. Add tiny jitter so no two rows are exactly identical
set.seed(42)
coords_km_jit <- coords_km + matrix(rnorm(nrow(coords_km)*2, sd = 1e-6), ncol = 2)

# 5. Build formula & design matrix
form_sb <- risk ~ status + fluid + material + location_type +
  diameter_in + length_ft + max_operating_pressure +
  line_age_yr + average_pop_density + Avg_Elevation
X_mat   <- model.matrix(form_sb, data = df_model)
p       <- ncol(X_mat)
n       <- nrow(df_model)

# 6. Fit Bayesian spatial GLM via MCMC
#    - We include 'w' in starting/tuning to satisfy spGLM API
#    - We bump up tau.sq (the nugget) so covariance matrix A + τ²I is PD
spmod <- spGLM(
  formula    = form_sb,
  data       = df_model,
  coords     = coords_km_jit,
  family     = "binomial",
  cov.model  = "exponential",
  n.samples  = 10000,
  
  # starting values
  starting   = list(
    beta     = rep(0,   p),      # flat starting for regression
    w        = rep(0,   n),      # spatial random effects
    phi      = 0.3,              # moderate range (in km⁻¹)
    sigma.sq = 1,                # process variance
    tau.sq   = 5                 # nugget variance ↑→ more diagonal dominance
  ),
  
  # tuning (Metropolis step sizes)
  tuning     = list(
    beta     = rep(0.01, p),
    w        = rep(0.01, n),
    phi      = 0.1,
    sigma.sq = 0.5,
    tau.sq   = 1
  ),
  
  # priors
  priors     = list(
    "beta.Flat",
    "phi.Unif"    = c(0.05, 2),  # allow φ between 0.05–2 km⁻¹
    "sigma.sq.IG" = c(2, 1),     # IG(2,1) mean=1
    "tau.sq.IG"   = c(2, 5)      # IG(2,5) mean=5 (stronger nugget)
  ),
  
  verbose = TRUE
)

# 7. Posterior predictive draws at original sites
sp_pred <- spPredict(
  spmod,
  pred.coords = coords_km_jit,   # use the same scaled+jittered coords
  pred.covars = X_mat,
  start       = 2000,            # burn-in
  thin        = 10
)

# 8. Summarize posterior predictive probabilities
ppd    <- sp_pred$p.y.predictive   # draws × locations
mean_p <- apply(ppd, 2, mean)
sd_p   <- apply(ppd, 2, sd)
lo95   <- apply(ppd, 2, quantile, 0.025)
hi95   <- apply(ppd, 2, quantile, 0.975)

df_model <- df_model %>%
  mutate(
    pred_risk_sb    = mean_p,
    pred_sd_sb      = sd_p,
    pred_lower95_sb = lo95,
    pred_upper95_sb = hi95,
    ci_width_sb     = pred_upper95_sb - pred_lower95_sb,
    pred_class_sb   = as.integer(pred_risk_sb > 0.5),
    pred_final_sb   = ifelse(risk == 1, 1L, pred_class_sb),
    match_sb        = ifelse(pred_class_sb == risk, "Match", "Mismatch")
  )

# 9. Print counts & confusion matrix
orig_counts    <- table(Original    = df_model$risk)
pred_counts_sb <- table(spBayesPred = df_model$pred_final_sb)
conf_mat_sb    <- table(Observed    = df_model$risk,
                        Predicted   = df_model$pred_final_sb)
print(orig_counts)
print(pred_counts_sb)
print(conf_mat_sb)

# 10. Summarize Bayesian uncertainty
print(summary(df_model$pred_sd_sb))
print(summary(df_model$ci_width_sb))

# 11. Leaflet maps

## 11.1 Original risk
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

## 11.2 Predicted risk (spBayes, asserted 1’s)
leaflet(df_model) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolylines(
    color   = ~ifelse(pred_final_sb == 1, "red", "blue"),
    weight  = 2, opacity = 0.7
  ) %>%
  addLegend(
    "bottomright",
    colors = c("blue","red"),
    labels = c("No Spill","Spill"),
    title  = "Predicted Risk (spBayes)"
  )

## 11.3 Match vs. Mismatch
leaflet(df_model) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolylines(
    color   = ~ifelse(match_sb == "Match", "green", "orange"),
    weight  = 2, opacity = 0.7
  ) %>%
  addLegend(
    "bottomright",
    colors = c("green","orange"),
    labels = c("Match","Mismatch"),
    title  = "Model vs Observed (spBayes)"
  )

## 11.4 Uncertainty (posterior SD)
pal_sb <- colorNumeric("YlOrRd", domain = df_model$pred_sd_sb)
leaflet(df_model) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolylines(
    color   = ~pal_sb(pred_sd_sb),
    weight  = 2, opacity = 0.8
  ) %>%
  addLegend(
    "bottomright",
    pal    = pal_sb,
    values = ~pred_sd_sb,
    title  = "Posterior SD (spBayes)"
  )

