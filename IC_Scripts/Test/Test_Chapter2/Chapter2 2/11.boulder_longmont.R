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

df_balanced <- read.csv("df_balanced.csv", stringsAsFactors = FALSE)

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
# 9) Observed Spill Outcomes (Risk 0 vs 1, consistent colors)
# =============================================================================
ggplot() +
  # First: plot the "no spill" points (risk = 0), ordered by predicted prob
  geom_point(
    data = df_balanced %>% arrange(predicted_prob) %>% filter(risk == 0),
    aes(x = lon, y = lat, color = factor(risk)),
    size = 1.7, alpha = 0.7
  ) +
  # Then: plot the "spill" points (risk = 1) so they always overlay
  geom_point(
    data = df_balanced %>% arrange(predicted_prob) %>% filter(risk == 1),
    aes(x = lon, y = lat, color = factor(risk)),
    size = 2.2, alpha = 0.8
  ) +
  scale_color_manual(
    values = c("0" = "#42a5f5", "1" = "#e53935"),  # same scheme as above
    name   = "Observed Risk",
    labels = c("0" = "No Spill", "1" = "Spill")
  ) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "Observed Spill Risk",
    x = "Longitude", y = "Latitude"
  )

# =============================================================================
# 10) Predicted Spill Probability (continuous turbo gradient)
# =============================================================================
ggplot(
  df_balanced %>% arrange(predicted_prob),  # low first so high overlays
  aes(x = lon, y = lat)
) +
  geom_point(aes(color = predicted_prob), size = 2) +
  scale_color_viridis_c(
    option = "turbo",
    begin  = 0.15,
    end    = 0.85,
    name   = "Predicted Probability",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = scales::percent_format(accuracy = 1),
    oob    = scales::squish
  ) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "Predicted Spill Probability – Without Z Matrix",
    x     = "Longitude", y = "Latitude"
  )
