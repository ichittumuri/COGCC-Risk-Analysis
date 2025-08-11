# New cell= command, option, I or alt, windows symbol, I 
# Run = windows symbol, enter 
# Run rscript = Command + Option + R

# =============================================================================
# 1) Setup & Libraries
# =============================================================================
# 1) Setup --------------------------------------------------------------------
library(sf)
library(dplyr)
library(ggplot2)
library(glmnet)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# =============================================================================
# 2) Sources & Data
# =============================================================================

usable_sf <- st_read("final_no_mkrig_dups.geojson")

# =============================================================================
# 3) Prep model data + quick class balance
# =============================================================================
coords_mat        <- st_coordinates(usable_sf)[, c("X","Y")]
usable_sf$lon     <- coords_mat[, 1]
usable_sf$lat     <- coords_mat[, 2]

df <- usable_sf |>
  st_drop_geometry() |>
  select(-unique_id) |>
  filter(
    lon >= -105.85,
    lon <= -104.79 + 0.2,
    lat >=  39.42,
    lat <=  40.68
  )

model_data <- df |>
  select(lon, lat, risk)

model_data %>%
  count(risk) %>%
  mutate(prop = n / sum(n)) %>%
  print()

# =============================================================================
# 4) Downsample 0-risk points to ~5% spills
# =============================================================================
set.seed(123) 

spills <- df %>% filter(risk == 1)
non_spills <- df %>% filter(risk == 0)

n_spills <- nrow(spills)
target_total <- ceiling(n_spills / 0.05)  # total rows for ~5% spills
n_non_spills_needed <- target_total - n_spills

non_spills_sample <- non_spills %>%
  sample_n(size = min(n_non_spills_needed, nrow(non_spills)))

df_balanced <- bind_rows(spills, non_spills_sample)

cat("\n--- NEW RATIO (~5% spills) ---\n")
df_balanced %>%
  count(risk) %>%
  mutate(prop = n / sum(n)) %>%
  print()

# =============================================================================
# 5) Observed Data Map
# =============================================================================
ggplot() +
  geom_point(
    data = df_balanced %>% filter(risk == 0),
    aes(x = lon, y = lat, color = factor(risk)),
    size = 1.7, alpha = 0.8
  ) +
  geom_point(
    data = df_balanced %>% filter(risk == 1),
    aes(x = lon, y = lat, color = factor(risk)),
    size = 1.7, alpha = 0.8
  ) +
  scale_color_manual(
    name = "Risk (0/1)",
    values = c("0" = "skyblue", "1" = "firebrick")
  ) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "Zoomed to Boulder–Greeley",
    x = "Longitude",
    y = "Latitude"
  )

# 5) Design matrix for LASSO ---------------------------------------------------
#    Use the SAME predictors you’ve been using
form <- risk ~ status + flowline_action + location_type + fluid + material +
  diameter_in + length_ft + max_operating_pressure + elevation + line_age_yr

y <- df_balanced$risk
X <- model.matrix(form, data = df_balanced)[, -1]  # drop intercept

# 6) CV-LASSO ------------------------------------------------------------------
set.seed(123)
cv_lasso <- cv.glmnet(
  x = X, y = y,
  family = "binomial",
  alpha  = 1,
  nfolds = 10
)

lambda_min <- cv_lasso$lambda.min
lambda_1se <- cv_lasso$lambda.1se
cat("\nSelected lambdas:\n",
    "  lambda.min =", signif(lambda_min, 3), "\n",
    "  lambda.1se =", signif(lambda_1se, 3), "\n")

# Extract nonzero at lambda.1se (more parsimonious)
lasso_coefs <- coef(cv_lasso, s = "lambda.1se")
nz_idx      <- which(as.numeric(lasso_coefs) != 0)
sel_vars    <- setdiff(rownames(lasso_coefs)[nz_idx], "(Intercept)")

cat("\nNonzero features at lambda.1se:\n")
print(sel_vars)

# Guard: if nothing selected, fall back to lambda.min
if (length(sel_vars) == 0) {
  message("No features at lambda.1se; using lambda.min instead.")
  lasso_coefs <- coef(cv_lasso, s = "lambda.min")
  nz_idx      <- which(as.numeric(lasso_coefs) != 0)
  sel_vars    <- setdiff(rownames(lasso_coefs)[nz_idx], "(Intercept)")
  cat("\nNonzero features at lambda.min:\n")
  print(sel_vars)
}

# 7) Refit baseline GLM on selected features ----------------------------------
X_sel          <- X[, sel_vars, drop = FALSE]
df_sel         <- as.data.frame(X_sel)
df_sel$risk    <- y
df_sel$lon     <- df_balanced$lon
df_sel$lat     <- df_balanced$lat

glm_lasso <- glm(risk ~ ., data = df_sel, family = binomial())
print(summary(glm_lasso))

# 8) Predictions + log-likelihood ---------------------------------------------
df_sel$predicted_prob <- predict(glm_lasso, type = "response")

p_hat <- df_sel$predicted_prob
yy    <- df_sel$risk
loglik_i     <- yy * log(p_hat) + (1 - yy) * log(1 - p_hat)
total_loglik <- sum(loglik_i)
cat("\nTotal log-likelihood (LASSO-selected GLM) =", round(total_loglik, 2), "\n")

# 9) Observed map --------------------------------------------------------------
ggplot() +
  geom_point(
    data = df_sel %>% filter(risk == 0),
    aes(x = lon, y = lat, color = factor(risk)),
    size = 1.7, alpha = 0.8
  ) +
  geom_point(
    data = df_sel %>% filter(risk == 1),
    aes(x = lon, y = lat, color = factor(risk)),
    size = 1.7, alpha = 0.8
  ) +
  scale_color_manual(
    values = c("0" = "skyblue", "1" = "firebrick"),
    name   = "Observed Risk",
    labels = c("0" = "No Spill", "1" = "Spill")
  ) +
  coord_fixed(xlim = xlims, ylim = ylims) +
  theme_minimal() +
  labs(
    title = "Observed Spill Risk – LASSO Subset",
    x = "Longitude", y = "Latitude"
  )

# 10) Predicted map (0–1 scale; high on top) ----------------------------------
ggplot(
  df_sel %>% arrange(predicted_prob),  # low → high so high overlays
  aes(x = lon, y = lat)
) +
  geom_point(aes(color = predicted_prob), size = 2) +
  scale_color_viridis_c(
    option = "D",
    name   = "Predicted Risk",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    oob    = scales::squish
  ) +
  coord_fixed(xlim = xlims, ylim = ylims) +
  theme_minimal() +
  labs(
    title = "Predicted Spill Risk – GLM with LASSO-selected Features",
    x = "Longitude", y = "Latitude"
  )

coords   <- as.matrix(df_balanced[, c("lon","lat")])
res_dev  <- residuals(glm_lasso, type = "deviance")
knn_nb   <- knearneigh(coords, k = 5)
nb       <- knn2nb(knn_nb, sym = TRUE)
lw       <- nb2listw(nb, style = "W", zero.policy = TRUE)
moran_out <- moran.test(res_dev, lw, zero.policy = TRUE)
print(moran_out)
