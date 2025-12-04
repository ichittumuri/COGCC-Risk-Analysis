# =============================================================================
# 1) Setup & Libraries
# =============================================================================
library(dplyr)
library(caret)
library(readr)
library(ggplot2)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# =============================================================================
# 2) Load Data & Define Formula
# =============================================================================
df_balanced <- read.csv("df_balanced.csv", stringsAsFactors = FALSE)

form <- risk ~ status + flowline_action + location_type + fluid + material +
  diameter_in + length_ft + max_operating_pressure + elevation + line_age_yr

# make sure cats are factors and risk is integer
cat_cols <- c("status","flowline_action","location_type","fluid","material")
df_balanced[cat_cols] <- lapply(df_balanced[cat_cols], factor)
df_balanced$risk <- as.integer(df_balanced$risk)

pred_names <- c("status","flowline_action","location_type","fluid","material",
                "diameter_in","length_ft","max_operating_pressure","elevation","line_age_yr")

# =============================================================================
# 3) One split -> return BOTH test log-likelihood and test predictions
# =============================================================================
cv_spilt <- function(df, p_train = 0.90, eps = 1e-15) {
  # stratified split
  df$risk_fac <- factor(df$risk, levels = c(0,1))
  idx <- createDataPartition(df$risk_fac, p = p_train, list = FALSE)
  train <- df[idx, ]
  test  <- df[-idx, ]
  
  # build dummies on TRAIN only
  x_train <- train[, pred_names]
  x_test  <- test[, pred_names]
  dv <- dummyVars(~ ., data = x_train, fullRank = TRUE)
  X_train <- as.data.frame(predict(dv, newdata = x_train))
  X_test  <- as.data.frame(predict(dv, newdata = x_test))
  
  # keep_vars <- c("fluidOther","elevation")
  # for (nm in setdiff(keep_vars, names(X_train))) X_train[[nm]] <- 0
  # for (nm in setdiff(keep_vars, names(X_test)))  X_test[[nm]]  <- 0
  # 
  # X_train <- X_train[, keep_vars, drop = FALSE]
  # X_test  <- X_test[,  keep_vars, drop = FALSE]
  
  # fit & predict on TEST
  fit  <- glm(risk ~ ., data = cbind(risk = train$risk, X_train), family = binomial())
  p_hat <- predict(fit, newdata = X_test, type = "response")
  
  # --- compute total log-likelihood on TEST set ---
  y <- test$risk
  loglik_i <- y * log(p_hat) + (1 - y) * log(1 - p_hat)
  total_loglik <- sum(loglik_i)
  
  # return BOTH: scalar loglik + row-level test predictions
  data.frame(
    row_id = test$row_id,      # make sure row_id exists (we'll add it below)
    p_hat  = p_hat,
    stringsAsFactors = FALSE
  ) -> test_preds
  
  return(list(loglik = total_loglik, test_preds = test_preds))
}

# add row_id once so we can aggregate later
df_balanced$row_id <- seq_len(nrow(df_balanced))

# =============================================================================
# 4) Run 200 splits: collect log-liks and test predictions
# =============================================================================
set.seed(2025)
res_list <- replicate(200, cv_spilt(df_balanced), simplify = FALSE)

# vector of total test log-likelihoods
loglik_vec <- sapply(res_list, `[[`, "loglik")

# bind all TEST predictions (only when a row was held out)
test_preds_all <- do.call(rbind, lapply(res_list, `[[`, "test_preds"))

# =============================================================================
# 5) Boxplot of log-likelihoods across 200 splits
# =============================================================================
loglik_df <- data.frame(LogLik = loglik_vec)
ggplot(loglik_df, aes(x = "", y = LogLik)) +
  geom_boxplot(alpha = 0.6, width = 0.3) +
  # geom_jitter(width = 0.1, alpha = 0.4) +
  labs(title = "Log-Likelihood across 200 CV splits",
       x = NULL, y = "Total Log-Likelihood (TEST set)") +
  theme_minimal(base_size = 14)

summary(loglik_df)

# =============================================================================
# 6) Aggregate TEST predictions to per-row MEAN (and SD, count)
# =============================================================================
p_mean <- tapply(test_preds_all$p_hat, test_preds_all$row_id, mean)
p_sd   <- tapply(test_preds_all$p_hat, test_preds_all$row_id, sd)
p_n    <- table(test_preds_all$row_id)

oof_summary <- data.frame(
  row_id = as.integer(names(p_mean)),
  p_mean = as.numeric(p_mean),
  p_sd   = as.numeric(p_sd),
  n_test = as.integer(p_n[match(names(p_mean), names(p_n))]),
  stringsAsFactors = FALSE
)

# join back lon/lat/risk for mapping
map_df <- merge(
  df_balanced[, c("row_id","lon","lat","risk")],
  oof_summary,
  by = "row_id",
  all.x = TRUE
)

# save
write.csv(map_df, "test_mean_predictions_200splits.csv", row.names = FALSE)

# how many times rows appeared as TEST
summary(map_df$n_test)

# =============================================================================
# 7) Maps: mean TEST predictions (and optional SD)
# =============================================================================

# 1) Observed Risk Map (true labels)
ggplot() +
  geom_point(data = subset(map_df, risk == 0),
             aes(x = lon, y = lat, color = factor(risk)),
             size = 1.7, alpha = 0.8) +
  geom_point(data = subset(map_df, risk == 1),
             aes(x = lon, y = lat, color = factor(risk)),
             size = 1.7, alpha = 0.8) +
  scale_color_manual(values = c("0" = "#42a5f5", "1" = "#e53935"),
                     name = "Observed Risk",
                     labels = c("No Spill", "Spill")) +
  coord_fixed() + theme_minimal() +
  labs(title = "Observed Spill Risk (all rows)",
       x = "Longitude", y = "Latitude")

# Mean predicted risk map (turbo)
ggplot(map_df[order(map_df$p_mean), ], aes(lon, lat)) +
  geom_point(aes(color = p_mean), size = 2) +
  scale_color_viridis_c(
    option = "turbo", limits = c(0,1),
    breaks = seq(0,1,0.25), oob = scales::squish,
    name = "Mean Test Risk"
  ) +
  coord_fixed() + theme_minimal() +
  labs(title = "Mean Predicted Risk (from TEST predictions only)",
       x = "Longitude", y = "Latitude")

# Uncertainty (SD across TEST predictions)
ggplot(map_df, aes(lon, lat)) +
  geom_point(aes(color = p_sd), size = 2) +
  scale_color_viridis_c(name = "SD of Test Preds") +
  coord_fixed() + theme_minimal() +
  labs(title = "Uncertainty: SD of Test Predictions (200 splits)",
       x = "Longitude", y = "Latitude")
