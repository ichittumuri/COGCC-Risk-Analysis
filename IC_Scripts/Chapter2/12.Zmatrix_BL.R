# New cell= command, option, I or alt, windows symbol, I 
# Run = windows symbol, enter 
# Run rscript = Command + Option + R
# example https://simonbrewer.github.io/geog5160/GEOG_5160_6160_lab04.html

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
    title = "Predicted Spill Risk – With Z Matrix",
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


# =============================================================================
# 12) export CIs
# =============================================================================
library(dplyr)
library(readr)
# if p_lo/p_hi came from apply/simulations, drop any dims:
pred_num <- as.numeric(df_ci$predicted_prob)  # or df_balanced$predicted_prob, if that's yours
lo_num   <- as.numeric(p_lo)                  # e.g. p_lo[,1] also works
hi_num   <- as.numeric(p_hi)

# stable id
id <- with(df_balanced, paste0(round(lon, 6), "_", round(lat, 6)))

out <- tibble(
  id    = id,
  lon   = as.numeric(df_balanced$lon),
  lat   = as.numeric(df_balanced$lat),
  pred  = pred_num,
  lo    = lo_num,
  hi    = hi_num,
  model = "Spatial + Covariates (LASSO)"   # <- change per script
)

write_csv(out, "model_out_spatial_covs.csv")






