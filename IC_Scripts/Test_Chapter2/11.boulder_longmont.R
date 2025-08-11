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
source("OGlogisticSmoother.R")
source("OGLambdaGridSearch.R")

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
# 6) Smoother fit + lambda search
# =============================================================================
coord_model <- as.matrix(df_balanced[, c("lon","lat")])
risk_model  <- df_balanced$risk

profile    <- LambdaGridSearch(coord_model, risk_model)
bestLambda <- profile$bestLambda

MLEFit <- logisticSmoother(
  coord_model,
  risk_model,
  lambda = bestLambda
)

# =============================================================================
# 7) Predict from fitted model
# =============================================================================
nu_all_mle <- predict(MLEFit, coord_model)
prob_all_mle <- exp(nu_all_mle) / (1 + exp(nu_all_mle))

df_balanced <- df_balanced %>%
  mutate(predicted_prob = prob_all_mle)

# =============================================================================
# 8) Total log-likelihood
# =============================================================================
nu_hat <- predict(MLEFit, coord_model)
p_hat  <- exp(nu_hat) / (1 + exp(nu_hat))  # same as plogis(nu_hat)
risk_model <- df_balanced$risk

loglik_i <- risk_model * log(p_hat) +
  (1 - risk_model) * log(1 - p_hat)

total_loglik <- sum(loglik_i)

cat("Total log-likelihood =", round(total_loglik, 2), "\n")

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
ggplot(df_balanced, aes(x = lon, y = lat)) +
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
    title = "Predicted Spill Risk – Without Z Matrix",
    x     = "Longitude",
    y     = "Latitude"
  )

# =============================================================================
# 11) Predicted probabilities ± 95% CI
# =============================================================================
nu_hat <- predict(MLEFit, coord_model)  # logit scale
p_hat  <- plogis(nu_hat)                 # prob scale

nu_se <- predictSE(MLEFit, x = coord_model)

p_se <- nu_se * p_hat * (1 - p_hat)

nu_lo <- nu_hat - 1.96 * nu_se
nu_hi <- nu_hat + 1.96 * nu_se

p_lo <- plogis(nu_lo)
p_hi <- plogis(nu_hi)

df_ci <- df_balanced %>%
  mutate(
    predicted_prob = p_hat,
    p_lo = p_lo,
    p_hi = p_hi
  ) %>%
  arrange(predicted_prob) %>%
  mutate(idx = row_number())

ggplot(df_ci, aes(x = idx, y = predicted_prob)) +
  geom_ribbon(aes(ymin = p_lo, ymax = p_hi), fill = "red", alpha = 0.3) +
  geom_line(color = "black", size = 1) +
  labs(
    x = "Observation (sorted by predicted risk)",
    y = "Predicted Spill Risk",
    title = "Predicted Spill Risk ± 95% CI"
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
