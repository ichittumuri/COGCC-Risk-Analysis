# New cell= command, option, I or alt, windows symbol, I 
# Run = windows symbol, enter 
# Run rscript = Command + Option + R

# =============================================================================
# 1) Setup & Libraries
# =============================================================================
library(sf)
library(dplyr)
library(spdep)
library(glmnet)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# =============================================================================
# 2) Load & Prepare Data
# =============================================================================
df <- st_read("final_no_mkrig_dups.geojson") |>
  st_transform(4326) |>
  mutate(across(c("status","flowline_action","location_type","fluid","material"), as.factor))

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
    !is.na(elevation),
    !is.na(line_age_yr)
  )

# =============================================================================
# 3) Print Categorical Levels
# =============================================================================
cats <- c("status", "flowline_action", "location_type", "fluid", "material")

n_levels  <- sapply(df_clean[cats], nlevels)
n_dummies <- n_levels - 1
print(n_levels)
print(n_dummies)

levels_list <- lapply(df_clean[cats], levels)
for (col_name in names(levels_list)) {
  cat(paste0("\nLevels for ", col_name, ":\n"))
  print(levels_list[[col_name]])
}

# =============================================================================
# 4) Logistic Regression
# =============================================================================
glm_model <- glm(
  risk ~ status + flowline_action + location_type + fluid + material +
    diameter_in + length_ft + max_operating_pressure + elevation + line_age_yr,
  data   = df_clean,
  family = binomial()
)
print(summary(glm_model))

# =============================================================================
# 5) Spatial Autocorrelation of Residuals (Moran's I)
# =============================================================================
df_clean$residuals <- residuals(glm_model, type = "deviance")

coords_mat <- st_coordinates(df_clean)[, c("X", "Y")]
df_clean   <- df_clean %>% st_drop_geometry()

knn_nb <- knearneigh(coords_mat, k = 5)
nb     <- knn2nb(knn_nb, sym = TRUE)
lw     <- nb2listw(nb, style = "W", zero.policy = TRUE)

moran_test <- moran.test(df_clean$residuals, lw, zero.policy = TRUE)
print(moran_test)

# =============================================================================
# 6) Feature Selection via LASSO (glmnet)
# =============================================================================
y <- df_clean$risk
X <- model.matrix(
  risk ~ status + flowline_action + location_type + fluid + material +
    diameter_in + length_ft + max_operating_pressure + elevation + line_age_yr,
  data = df_clean
)[, -1]

set.seed(123)
cv_lasso <- cv.glmnet(x = X, y = y, family = "binomial", alpha = 1, nfolds = 10)

lasso_coefs <- coef(cv_lasso, s = "lambda.1se")
nz_idx      <- which(lasso_coefs != 0)
sel_vars    <- setdiff(rownames(lasso_coefs)[nz_idx], "(Intercept)")
sel_vars

X_sel       <- X[, sel_vars, drop = FALSE]
df_sel      <- as.data.frame(X_sel)
df_sel$risk <- y

glm_lasso <- glm(risk ~ ., data = df_sel, family = binomial())
summary(glm_lasso)



