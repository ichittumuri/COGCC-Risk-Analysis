# =============================================================================
# 1) Setup & Libraries
# =============================================================================
library(dplyr)
library(caret)
library(readr)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# =============================================================================
# 2) Load Data & Define Formula
# =============================================================================
df_balanced <- read.csv("df_balanced.csv", stringsAsFactors = FALSE)

form <- risk ~ status + flowline_action + location_type + fluid + material +
  diameter_in + length_ft + max_operating_pressure + elevation + line_age_yr

# Make sure categoricals are factors
cat_cols <- c("status", "flowline_action", "location_type", "fluid", "material")
df_balanced[cat_cols] <- lapply(df_balanced[cat_cols], factor)

# Response as numeric 0/1
df_balanced$risk <- as.integer(df_balanced$risk)

# =============================================================================
# 3) Single Train-Test Split (90/10, Stratified by Risk)
# =============================================================================
set.seed(42)

dat_fac <- df_balanced %>%
  mutate(risk = factor(risk, levels = c(0, 1)))

idx_90 <- createDataPartition(dat_fac$risk, p = 0.90, list = FALSE)

train_data <- dat_fac[idx_90, ] %>%
  mutate(risk = as.integer(as.character(risk)))

test_data <- dat_fac[-idx_90, ] %>%
  mutate(risk = as.integer(as.character(risk)))

# quick prevalence check
summarize_split <- function(data, name) {
  n_total <- nrow(data)
  pct_spill <- mean(data$risk == 1)
  cat(sprintf("%s: n = %d | p(Spill) = %.3f\n", name, n_total, pct_spill))
}
summarize_split(train_data, "Train (90%)")
summarize_split(test_data,  "Test  (10%)")

# =============================================================================
# 4) Fit Base Logistic Model and Evaluate TEST Log-Likelihood
# =============================================================================
fit_glm <- glm(form, data = train_data, family = binomial)

# Predict probabilities on TEST set
p_hat <- predict(fit_glm, newdata = test_data, type = "response")
y     <- test_data$risk

loglik_i <- y * log(p_hat) + (1 - y) * log(1 - p_hat)
total_loglik <- sum(loglik_i)

cat("Test set log-likelihood =", round(total_loglik, 2), "\n")

summary(fit_glm)
confint(fit_glm)


