# New cell= command, option, I or alt, windows symbol, I 
# Run = windows symbol, enter 
# Run rscirpt =  Command + Option + R

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/IC_Scripts/2024-2025/sp_point_data")

library(sf)        
library(dplyr)     
library(leaflet)   
library(spdep)     
library(tidyr) 
library(readr)

df <- st_read("usable_point_data.geojson") %>%
  st_transform(crs = 4326)

df_clean <- df %>%
  filter(
    !is.na(status),
    !is.na(flowline_action),
    !is.na(location_type),
    !is.na(fluid),
    !is.na(material),
    !is.na(diameter_in),
    !is.na(length_ft),
    !is.na(max_operating_pressure),
    !is.na(avg_population),
    !is.na(avg_elevation),
    !is.na(line_age_yr)
  )

# 1a) Coerce your columns to factors
cats <- c("status", "flowline_action", "location_type", "fluid", "material")

df_clean_count <- df_clean %>%
  mutate(across(all_of(cats), as.factor))

# 1b) Now count levels
n_levels <- sapply(df_clean_count[cats], nlevels)
n_levels

# 1c) Number of dummy (one‐hot) columns = levels − 1
n_dummies <- n_levels - 1
n_dummies

# Print level names for each factor column
levels_list <- lapply(df_clean_count[cats], levels)

# Print each set of levels with the column name
for (col_name in names(levels_list)) {
  cat(paste0("\nLevels for ", col_name, ":\n"))
  print(levels_list[[col_name]])
}

# Logistic Regression --------------------------------------------------

glm_model <- glm(
  risk ~ status + flowline_action+ location_type + fluid + material +
    diameter_in + length_ft + max_operating_pressure + avg_elevation +
    line_age_yr,
  data   = df_clean,
  family = binomial,
  control = glm.control(epsilon = 1e-8, maxit = 50)
)

summary(glm_model)

# Confusion Matrix -----------------------------------------------------

df_clean$pred_prob  <- predict(glm_model, type = "response")
df_clean$pred_class <- ifelse(df_clean$pred_prob > 0.5, 1, 0)
conf_mat <- table(Observed = df_clean$risk, Predicted = df_clean$pred_class)
print(conf_mat)

# Spatial Autocorrelation of Residuals ---------------------------------
df_clean$residuals <- residuals(glm_model, type = "deviance")
df_clean <- df_clean %>% 
  st_drop_geometry()

coords_mat <- df_clean %>%
  mutate(coords = as.character(coords)) %>%
  separate(coords, into = c("X","Y"), sep = ",\\s*", remove = TRUE) %>%
  mutate(
    X = parse_number(X),
    Y = parse_number(Y)
  ) %>%
  select(X, Y) %>%
  as.matrix()

head(coords_mat)

set.seed(42)
coords_jit <- coords_mat
coords_jit[,1] <- jitter(coords_mat[,1], factor = 1e-8)
coords_jit[,2] <- jitter(coords_mat[,2], factor = 1e-8)

knn_nb     <- knearneigh(coords_jit, k = 5)
nb         <- knn2nb(knn_nb, sym = TRUE)
lw         <- nb2listw(nb, style = "W", zero.policy = TRUE)

moran_test <- moran.test(df_clean$residuals, lw, zero.policy = TRUE)
print(moran_test)

# Feature Selection via Stepwise AIC/BIC -------------------------------

full_model <- glm(risk ~ status + flowline_action+ location_type + fluid + material +
                    diameter_in + length_ft + max_operating_pressure + avg_elevation +
                    line_age_yr,
                  data   = df_clean,
                  family = binomial)

null_model <- glm(risk ~ 1, data = df_clean, family = binomial)

n_obs          <- nrow(df_clean)
stepwise_AIC <- step(null_model,
                       scope     = list(lower = null_model, upper = full_model),
                       direction = "both",
                       k=2) 
                      # 2=AIC, log(n_obs)= BIC 

summary(stepwise_AIC)

# Feature Selection via LASSO -------------------------------

library(glmnet)
y <- df_clean$risk
X <- model.matrix(
  risk ~ status + flowline_action+ location_type + fluid + material +
    diameter_in + length_ft + max_operating_pressure + avg_elevation +
    line_age_yr,
  data = df_clean
)[, -1]  

set.seed(123)
cv_lasso   <- cv.glmnet(x = X, y = y, family = "binomial", alpha = 1, nfolds = 10) #alpha 1=lasso, 0=Ridge, 0 < α < 1 = elastic net

lasso_coefs <- coef(cv_lasso, s = "lambda.1se") # largest lambda of mean CV within 1SE of the min
nz_idx <- which(lasso_coefs != 0)
sel_vars <- rownames(lasso_coefs)[nz_idx]
sel_vars <- setdiff(sel_vars, "(Intercept)")
sel_vars

X_sel  <- X[, sel_vars, drop = FALSE]
df_sel <- as.data.frame(X_sel)
df_sel$risk <- y

glm_lasso <- glm(risk ~ ., data = df_sel, family = binomial)
summary(glm_lasso)

