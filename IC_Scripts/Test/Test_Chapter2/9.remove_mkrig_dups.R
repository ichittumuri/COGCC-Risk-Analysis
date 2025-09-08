library(sf)
library(fields)          

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# 1. Load data and extract coordinates
df <- st_read("final_dataset_subset.geojson")
coords_mat <- st_coordinates(df)[, c("X", "Y")]

# 2. Check duplicates (R default: rounded numeric comparison)
coord_rounded <- round(coords_mat, 6)
dupe_rows <- duplicated(as.data.frame(coord_rounded))
cat("R-level duplicates (rounded):", sum(dupe_rows), "\n")

# 3. Check duplicates (mKrig-style: convert to strings)
coord_strings <- cat.matrix(coords_mat)
dups_mkrig <- duplicated(coord_strings)
cat("mKrig-level duplicates:", sum(dups_mkrig), "\n") 

# 4. Inspect risk distribution for duplicates
risk_dups <- df$risk[dups_mkrig]
cat("Risk distribution for mKrig duplicates:\n")
print(table(risk_dups))

# 5. Drop mKrig-level duplicates from the sf
df_clean <- df[!dups_mkrig, ]

# 6. Final check
cat("Final number of points after removing mKrig-level duplicates:", nrow(df_clean), "\n")
cat("Risk distribution after cleaning:\n")
print(table(df_clean$risk))

# 7. Check for any NAs in the dataset
na_counts <- colSums(is.na(df))
cat("\n--- NA counts per column ---\n")
print(na_counts)

if (any(na_counts > 0)) {
  cat("\nColumns with NAs:\n")
  print(names(df)[na_counts > 0])
} else {
  cat("\nNo NAs found in the dataset.\n")
}

# 8. Drop NAs from dataset
df_no_na <- na.omit(df_clean)

cols_to_check <- c("flowline_action", "fluid", "material")
df_no_na <- df_clean[complete.cases(st_drop_geometry(df_clean)[, cols_to_check]), ]

cat("Rows before:", nrow(df_clean), "\n")
cat("Rows after removing NAs:", nrow(df_no_na), "\n")

# 9. Final check
cat("Final number of points after removing NAs:", nrow(df_no_na), "\n")
cat("Risk distribution after cleaning:\n")
print(table(df_no_na$risk))

# Optional: save cleaned dataset
st_write(df_no_na, "final_no_mkrig_dups.geojson", delete_dsn = TRUE)


