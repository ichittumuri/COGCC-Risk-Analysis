# ============================================
# Setup
# ============================================
library(sf)
library(fields)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

FLOWLINES_FP <- "flowline_points_50m_dedup.geojson"
SPILLS_FP    <- "spills_w_flowline_attributes.geojson"

FLOWLINES_CLEAN_FP <- "flowlines_dedup_nonNA.geojson"
SPILLS_CLEAN_FP    <- "spills_dedup_nonNA.geojson"
COMBINED_FP        <- "combined_flowlines_spills.geojson"

# columns to require non-NA (adjust as needed)
cols_to_check <- c("flowline_action","fluid","material")

# --- shim: ensure cat.matrix exists & returns a character vector key per row ---
if (!exists("cat.matrix")) {
  cat.matrix <- function(M) apply(M, 1, function(r) paste(r, collapse = " "))
}

# ============================================
# Step 2 — audit & clean (EXACT style) for FLOWLINES
# ============================================

cat("\n========== FLOWLINES ==========\n")
df <- flowlines

coords_mat <- st_coordinates(df)[, c("X","Y")]

# R-level dups (rounded)
coord_rounded <- round(coords_mat, 6)
dupe_rows <- duplicated(as.data.frame(coord_rounded))
cat("R-level duplicates (rounded):", sum(dupe_rows), "\n")

# mKrig-style dups (convert rows to strings)
coord_strings <- cat.matrix(coords_mat)
dups_mkrig <- duplicated(coord_strings)
cat("mKrig-level duplicates:", sum(dups_mkrig), "\n")

# risk distribution among mKrig dups (if present)
if ("risk" %in% names(df)) {
  risk_dups <- df$risk[dups_mkrig]
  cat("Risk distribution for mKrig duplicates:\n")
  print(table(risk_dups))
}

# drop mKrig-level dups
df_clean <- df[!dups_mkrig, ]
cat("Final number of points after removing mKrig-level duplicates:", nrow(df_clean), "\n")
if ("risk" %in% names(df_clean)) {
  cat("Risk distribution after cleaning:\n")
  print(table(df_clean$risk))
}

# NA audit on original df (to mirror your pattern)
na_counts <- colSums(is.na(st_drop_geometry(df)))
cat("\n--- NA counts per column (original df) ---\n")
print(na_counts)
if (any(na_counts > 0)) {
  cat("\nColumns with NAs:\n")
  print(names(df)[na_counts > 0])
} else {
  cat("\nNo NAs found in the dataset.\n")
}

# Drop NAs from cleaned df on specific columns
rows_before <- nrow(df_clean)
df_no_na <- df_clean[
  complete.cases(st_drop_geometry(df_clean)[, intersect(cols_to_check, names(df_clean)), drop = FALSE]),
]
cat("Rows before:", rows_before, "\n")
cat("Rows after removing NAs:", nrow(df_no_na), "\n")

if ("risk" %in% names(df_no_na)) {
  cat("Final risk distribution after NA removal:\n")
  print(table(df_no_na$risk))
}

flowlines_clean <- df_no_na
st_write(flowlines_clean, FLOWLINES_CLEAN_FP, delete_dsn = TRUE, quiet = TRUE)

# ============================================
# Step 3 — audit & clean (EXACT style) for SPILLS
# ============================================

cat("\n========== SPILLS ==========\n")
df <- spills

coords_mat <- st_coordinates(df)[, c("X","Y")]

# R-level dups (rounded)
coord_rounded <- round(coords_mat, 6)
dupe_rows <- duplicated(as.data.frame(coord_rounded))
cat("R-level duplicates (rounded):", sum(dupe_rows), "\n")

# mKrig-style dups (convert rows to strings)
coord_strings <- cat.matrix(coords_mat)
dups_mkrig <- duplicated(coord_strings)
cat("mKrig-level duplicates:", sum(dups_mkrig), "\n")

# risk distribution among mKrig dups (if present)
if ("risk" %in% names(df)) {
  risk_dups <- df$risk[dups_mkrig]
  cat("Risk distribution for mKrig duplicates:\n")
  print(table(risk_dups))
}

# drop mKrig-level dups
df_clean <- df[!dups_mkrig, ]
cat("Final number of points after removing mKrig-level duplicates:", nrow(df_clean), "\n")
if ("risk" %in% names(df_clean)) {
  cat("Risk distribution after cleaning:\n")
  print(table(df_clean$risk))
}

# NA audit on original df (to mirror your pattern)
na_counts <- colSums(is.na(st_drop_geometry(df)))
cat("\n--- NA counts per column (original df) ---\n")
print(na_counts)
if (any(na_counts > 0)) {
  cat("\nColumns with NAs:\n")
  print(names(df)[na_counts > 0])
} else {
  cat("\nNo NAs found in the dataset.\n")
}

# Drop NAs from cleaned df on specific columns
rows_before <- nrow(df)
df_no_na <- df_clean[
  complete.cases(st_drop_geometry(df_clean)[, intersect(cols_to_check, names(df_clean)), drop = FALSE]),
]
cat("Rows before:", rows_before, "\n")
cat("Rows after removing NAs:", nrow(df_no_na), "\n")

if ("risk" %in% names(df_no_na)) {
  cat("Final risk distribution after NA removal:\n")
  print(table(df_no_na$risk))
}

spills_clean <- df_no_na
st_write(spills_clean, SPILLS_CLEAN_FP, delete_dsn = TRUE, quiet = TRUE)


# ============================================
# Step 4 — remove overlaps (keep spill) + combine
# ============================================

# build rounded 6dp keys to find overlaps
fl_xy <- st_coordinates(flowlines_clean)[, c("X","Y")]
sp_xy <- st_coordinates(spills_clean)[,    c("X","Y")]

fl_key <- apply(round(fl_xy, 6), 1, paste, collapse = ",")
sp_key <- apply(round(sp_xy, 6), 1, paste, collapse = ",")

common_keys <- intersect(fl_key, sp_key)
cat("\n=== OVERLAP AUDIT (rounded 6 d.p.) ===\n")
cat("Flowlines BEFORE:", nrow(flowlines_clean), "\n")
cat("Spills:",           nrow(spills_clean),    "\n")
cat("Overlapping coord keys:", length(common_keys), "\n")

flowlines_no_overlap <- flowlines_clean[!(fl_key %in% common_keys), ]
cat("Flowlines AFTER removing overlaps:", nrow(flowlines_no_overlap), "\n")

# combine
if (st_crs(flowlines_no_overlap) != st_crs(spills_clean)) {
  spills_clean <- st_transform(spills_clean, st_crs(flowlines_no_overlap))
}

library(dplyr)

# 1) strip geom -> data.frames
fl_df  <- st_drop_geometry(flowlines_no_overlap)
sp_df  <- st_drop_geometry(spills_clean)
fl_g   <- st_geometry(flowlines_no_overlap)
sp_g   <- st_geometry(spills_clean)

# 2) union columns like pandas concat(ignore_index=True)
combined_df <- bind_rows(fl_df, sp_df)  # handles different columns

# 3) reattach geometry and CRS
combined <- st_sf(combined_df,
                  geometry = c(fl_g, sp_g),
                  crs      = st_crs(flowlines_no_overlap))

cat("Rows combined:", nrow(combined), "\n")



# ============================================
# Step 5 Checked combined mKrig dup counts
# ============================================

if (!exists("cat.matrix")) {
  cat.matrix <- function(M) apply(M, 1, function(r) paste(r, collapse = " "))
}

cat("\n========== COMBINED ==========\n")
df <- combined

coords_mat <- st_coordinates(df)[, c("X","Y")]

# R-level dups (rounded)
coord_rounded <- round(coords_mat, 6)
dupe_rows <- duplicated(as.data.frame(coord_rounded))
cat("R-level duplicates (rounded):", sum(dupe_rows), "\n")

# mKrig-style dups
coord_strings <- cat.matrix(coords_mat)
dups_mkrig <- duplicated(coord_strings)
cat("mKrig-level duplicates:", sum(dups_mkrig), "\n")

# risk distribution among mKrig dups (if present)
if ("risk" %in% names(df)) {
  risk_dups <- df$risk[dups_mkrig]
  cat("Risk distribution for mKrig duplicates:\n")
  print(table(risk_dups))
}

# drop mKrig-level dups
df_clean <- df[!dups_mkrig, ]
cat("Final number of points after removing mKrig-level duplicates:", nrow(df_clean), "\n")

if ("risk" %in% names(df_clean)) {
  cat("Risk distribution after cleaning:\n")
  print(table(df_clean$risk))
}

combined_clean <- df_clean
st_write(combined_clean, COMBINED_FP, delete_dsn = TRUE, quiet = TRUE)
cat("\nSaved:", COMBINED_FP, "\n")
