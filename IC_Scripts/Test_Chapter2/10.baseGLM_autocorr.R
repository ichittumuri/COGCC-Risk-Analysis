# New cell= command, option, I or alt, windows symbol, I 
# Run = windows symbol, enter 
# Run rscript = Command + Option + R

# =============================================================================
# 1) Setup & Libraries
# =============================================================================
library(sf)
library(dplyr)
library(ggplot2)
library(fields)

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

# =============================================================================
# 6) Fit baseline GLM (no spatial term)
# =============================================================================
glm_baseline <- glm(
  risk ~ status + flowline_action + location_type + fluid + material +
    diameter_in + length_ft + max_operating_pressure + elevation + line_age_yr,
  data   = df_balanced,
  family = binomial()
)
print(summary(glm_baseline))

# =============================================================================
# 7) Predict from fitted model
# =============================================================================
df_balanced$predicted_prob <- predict(glm_baseline, type = "response")

# =============================================================================
# 8) Total log-likelihood
# =============================================================================
p_hat <- df_balanced$predicted_prob
y     <- df_balanced$risk
loglik_i <- y * log(p_hat) + (1 - y) * log(1 - p_hat)
total_loglik <- sum(loglik_i)
cat("Total log-likelihood (baseline GLM) =", round(total_loglik, 2), "\n")

# =============================================================================
# 9) Observed Risk plot
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
    values = c("0" = "skyblue", "1" = "firebrick"),
    name   = "Observed Risk",
    labels = c("0" = "No Spill", "1" = "Spill")
  ) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "Observed Spill Risk",
    x = "Longitude",
    y = "Latitude"
  )

# =============================================================================
# 10) Predicted Risk plot
# =============================================================================
ggplot(
  df_balanced %>% arrange(predicted_prob),  # sort low → high
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
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "Predicted Spill Risk – Baseline GLM (No Spatial Smoother)",
    x = "Longitude", y = "Latitude"
  )

# =============================================================================
# 11) Predicted probabilities ± 95% CI  (BASELINE GLM ONLY)
# =============================================================================
# 1) predict on the LINK (logit) scale *from the glm*, with SEs
pred_link <- predict(glm_baseline, type = "link", se.fit = TRUE)

eta   <- pred_link$fit              # logit
se_eta<- pred_link$se.fit           # SE on logit
p_hat <- plogis(eta)                # prob

# 2) delta method: SE on probability scale
se_p  <- se_eta * p_hat * (1 - p_hat)

# 3) 95% CI on prob scale (logit CI, then transform is also fine)
p_lo <- pmax(0, pmin(1, plogis(eta - 1.96 * se_eta)))
p_hi <- pmax(0, pmin(1, plogis(eta + 1.96 * se_eta)))

df_ci <- df_balanced %>%
  mutate(predicted_prob = p_hat, p_lo = p_lo, p_hi = p_hi) %>%
  arrange(predicted_prob) %>%
  mutate(idx = row_number())

ggplot(df_ci, aes(x = idx, y = predicted_prob)) +
  geom_ribbon(aes(ymin = p_lo, ymax = p_hi), fill = "red", alpha = 0.3) +
  geom_line(color = "black", size = 1) +
  labs(
    x = "Observation (sorted by predicted risk)",
    y = "Predicted Spill Risk",
    title = "Predicted Spill Risk ± 95% CI (Baseline GLM)"
  ) +
  theme_minimal()

# =============================================================================
# 12) Pearson residuals + QQ plot
# =============================================================================
pearson_resid <- (risk_model - prob_all_mle) / sqrt(prob_all_mle * (1 - prob_all_mle))

qqnorm(pearson_resid,
       main = "QQ plot of Pearson Residuals",
       xlab = "Theoretical Quantiles",
       ylab = "Observed Residuals")
qqline(pearson_resid, lty = 2)

# =============================================================================
# 13) Dunn–Smyth residuals (randomized quantile residuals, for Bernoulli (n = 1))
# =============================================================================
set.seed(123)

eps <- 1e-12
p <- pmin(pmax(prob_all_mle, eps), 1 - eps)   # clamp away from 0/1
y <- as.integer(risk_model)

lower <- ifelse(y == 0, 0,       1 - p)  # F(y^-)
upper <- ifelse(y == 0, 1 - p,   1)      # F(y)
u     <- runif(length(y), lower, upper)
r_ds  <- qnorm(u)

cat("Dunn–Smyth residuals: mean =", round(mean(r_ds),4),
    "sd =", round(sd(r_ds),4), "\n")

qqnorm(r_ds,
       main = "QQ plot of Dunn–Smyth Residuals (Bernoulli)",
       xlab = "Theoretical Quantiles", ylab = "Observed Residuals")
qqline(r_ds, lty = 2)

df_balanced$resid_ds <- r_ds

plot(p, r_ds, pch = 19, cex = 0.6, # residual vs fitted
     xlab = "Fitted probability", ylab = "Dunn–Smyth residual")
abline(h = 0, lty = 2)

# =============================================================================
# 12) Moran's I of GLM residuals to confirm SAC
# =============================================================================
# Keep if you want to re-check spatial autocorrelation on the *balanced* set
# using the same 5-NN structure as before.
coords   <- as.matrix(df_balanced[, c("lon","lat")])
res_dev  <- residuals(glm_baseline, type = "deviance")
knn_nb   <- knearneigh(coords, k = 5)
nb       <- knn2nb(knn_nb, sym = TRUE)
lw       <- nb2listw(nb, style = "W", zero.policy = TRUE)
moran_out <- moran.test(res_dev, lw, zero.policy = TRUE)
print(moran_out)

