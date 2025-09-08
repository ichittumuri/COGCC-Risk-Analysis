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
source("ZlogisticSmoother.R")       
source("ZLambdaGridSearch.R")

df_balanced <- read.csv("df_balanced.csv", stringsAsFactors = FALSE)
 
# =============================================================================
# 6) Fixed-effects, design Z matrix
# =============================================================================
to.factor <- c("status","flowline_action","location_type","fluid","material")
df_balanced <- df_balanced %>% mutate(across(all_of(to.factor), as.factor))

X <- model.matrix(
  ~ status + flowline_action + location_type + fluid + material +
    diameter_in + length_ft + elevation + line_age_yr,
  data = df_balanced
)[, -1]

# select Z columns (only those that exist)
Z <- X[, c(
  "statusNew Construction",
  "fluidOther",
  "elevation"
)]

risk_model   <- df_balanced$risk
coord_model  <- as.matrix(df_balanced[, c("lon","lat")])

Z_df <- as.data.frame(Z)
Z_df$risk <- df_balanced$risk           

glm_sel <- glm(risk ~ . -1, data = Z_df, family = binomial(link = "logit"))
betaHat <- coef(glm_sel)

# =============================================================================
# 7) Lambda search + final fit with Z and betaHat
# =============================================================================
profile <- LambdaGridSearch(
  coord_model,
  risk_model,
  Z       = Z,
  betaHat = betaHat
)

bestLambda <- profile$bestLambda

MLEFit <- logisticSmoother(
  coord_model,
  risk_model,
  lambda  = bestLambda,
  Z       = Z,
  betaHat = betaHat
)

# =============================================================================
# 7) Predict from fitted model (with Z)
# =============================================================================
nu_all_mle <- predict(MLEFit, coord_model, Z = Z)
prob_all_mle <- plogis(nu_all_mle)  # same as exp(nu)/(1+exp(nu))

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

# 9) Observed Spill Outcomes — spill points always on top
ggplot() +
  # First: plot the "no spill" points (risk = 0) ordered by predicted prob
  geom_point(
    data = df_balanced %>% arrange(predicted_prob) %>% filter(risk == 0),
    aes(x = lon, y = lat, color = factor(risk)),
    size = 1.7, alpha = 0.7
  ) +
  # Then: plot the "spill" points (risk = 1) so they always sit on top
  geom_point(
    data = df_balanced %>% arrange(predicted_prob) %>% filter(risk == 1),
    aes(x = lon, y = lat, color = factor(risk)),
    size = 2.2, alpha = 0.8
  ) +
  scale_color_manual(
    values = c("0" = "#42a5f5", "1" = "#e53935"),  # light blue vs strong red
    name   = "Observed Outcome",
    labels = c("0" = "No Spill", "1" = "Spill")
  ) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "Observed Spill Outcomes",
    x = "Longitude", y = "Latitude"
  )


# 10) Predicted Spill Probability — with covariates (low → high layering)
ggplot(
  df_balanced %>% arrange(predicted_prob),
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
    title = "Predicted Spill Probability (with covariates)",
    x = "Longitude", y = "Latitude"
  )
