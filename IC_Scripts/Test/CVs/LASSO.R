# New cell= command, option, I or alt, windows symbol, I 
# Run = windows symbol, enter 
# Run rscript = Command + Option + R
# example https://simonbrewer.github.io/geog5160/GEOG_5160_6160_lab04.html

# =============================================================================
# 1) Setup & Libraries
# =============================================================================
# 1) Setup --------------------------------------------------------------------
library(sf)
library(dplyr)
library(ggplot2)
library(glmnet)
library(spdep)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# =============================================================================
# 2) Sources & Data
# =============================================================================

usable_sf <- st_read("final_dataset_subset.geojson")

usable_sf[] <- lapply(usable_sf, function(col) {
  if (is.factor(col) || is.character(col)) {
    col <- gsub(" ", "_", col)      # replace spaces in the values
    col <- factor(col)              # ensure factors stay factors
  }
  col
})

# =============================================================================
# 3) Prep model data + quick class balance
# =============================================================================
coords_mat        <- st_coordinates(usable_sf)[, c("X","Y")]
usable_sf$lon     <- coords_mat[, 1]
usable_sf$lat     <- coords_mat[, 2]

df <- usable_sf |>
  st_drop_geometry() |>
  # select(-unique_id) |>
  filter(
    lon >= -105.85,
    lon <= -104.79 + 0.2,
    lat >=  39.42,
    lat <=  40.68
  )

# =============================================================================
# 4) Downsample 0-risk points to ~5% spills
# =============================================================================
set.seed(42)

# Count positives/negatives
N1 <- sum(df$risk == 1, na.rm = TRUE)
N0 <- sum(df$risk == 0, na.rm = TRUE)

# Compute proportional keep-fraction for negatives
neg_frac <- min(1, (19 * N1) / N0)  # target ~5% positives overall

# Downsample NEGATIVES by the SAME fraction within each unique_id
neg_down <- df %>%
  filter(risk == 0) %>%
  group_by(unique_id) %>%
  sample_frac(neg_frac) %>%   # keeps the same proportion per line
  ungroup()

# Keep ALL positives
pos_keep <- df %>% filter(risk == 1)

# Combine back to point level
df_prop_down <- bind_rows(pos_keep, neg_down)

# Check achieved proportions
tab <- df_prop_down %>%
  count(risk) %>%
  mutate(prop = n / sum(n))

tab

write.csv(df_prop_down, "df_balanced.csv", row.names = FALSE)

# 5) Design matrix for LASSO ---------------------------------------------------
form <- risk ~ status + flowline_action + location_type + fluid + material +
  diameter_in + length_ft + max_operating_pressure + elevation + line_age_yr

df_balanced <- df_prop_down

y <- df_balanced$risk
X <- model.matrix(form, data = df_balanced)[, -1]  # drop intercept

# 6) CV-LASSO ------------------------------------------------------------------
set.seed(2025)
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
