# Standalone end-to-end test harness for AthensDemo's glmnet_var_select block.
#
# Why this exists: the sandbox this was written in has no R interpreter and
# no network access to CRAN/conda, so the actual cv.glmnet run could not be
# executed there. This script lets you run the *real, unmodified*
# R/glmnet_var_select.R end-to-end on synthetic data, outside of spade,
# so you can confirm the PDF renders cleanly before July 29 without needing
# a live spade pipeline run first.
#
# Usage: run from the AthensDemo repo root (so "R/glmnet_var_select.R"
# resolves), with spadelib, arrow, dplyr, and glmnet already installed:
#
#   Rscript test_glmnet_var_select.R
#
# It writes glmnet_diagnostics.pdf to the current working directory --
# open it and check: coefficient table + recommendation (p1), intercept
# diagnostic (p2), predictor color key (p3), one CV/path page per alpha
# (p4-8 for the default 5-alpha grid), guidance page (last). 9 pages total
# for the default alpha grid of 5 values.

library(arrow)

setwd("/home/pradtke/source/AthensDemo")
block_file <- "R/glmnet_var_select.R"
if (!file.exists(block_file)) {
  stop("Run this from the AthensDemo repo root -- ", block_file, " not found.")
}

# Source everything in the block file except the trailing `run(handler)`
# call, so we get the real `handler` function (and its spade_types
# attribute) without needing a params.yaml / inputs/ spade harness.
lines <- readLines(block_file)
run_line <- grep("^run\\(handler\\)\\s*$", lines)
if (length(run_line) != 1) {
  stop("Expected exactly one top-level `run(handler)` line in ", block_file,
       " -- found ", length(run_line), ". Check the file hasn't changed shape.")
}
eval(parse(text = paste(lines[-run_line], collapse = "\n")), envir = .GlobalEnv)

# ---- Build synthetic data mimicking join_chm_stats.R's output schema ----
# 8 raw-km^2 canopy-height bins (HT0..HT35, per join_chm_stats.R's
# cut(breaks = c(-1,0,5,...,35)) binning), for a ~100-county run like the
# planned July 29 demo. Response is a zero-intercept linear combination of
# a few bins (mimicking Cao et al.'s zero-canopy/zero-volume physical
# relationship) plus noise, so the intercept diagnostic and predictor
# recommendation should both come back sensible.
set.seed(1)
n_counties <- 100
ht_cols <- paste0("HT", c(0, 5, 10, 15, 20, 25, 30, 35))

X <- sapply(ht_cols, function(nm) pmax(0, rgamma(n_counties, shape = 2, scale = 8)))
colnames(X) <- ht_cols

true_coef <- c(HT0 = 0, HT5 = 0, HT10 = 0.2, HT15 = 0.3,
                HT20 = 1.1, HT25 = 1.8, HT30 = 2.4, HT35 = 2.9)
y_true <- as.numeric(X %*% true_coef[ht_cols])
noise <- rnorm(n_counties, sd = 0.08 * sd(y_true))
ESTIMATE <- y_true + noise

df <- data.frame(
  STATECD = 51,
  COUNTYCD = seq_len(n_counties),
  CO_FIPS = 51000 + seq_len(n_counties),
  ESTIMATE = ESTIMATE,
  VARIANCE = (0.1 * abs(ESTIMATE))^2
)
df <- cbind(df, as.data.frame(X))

tmp_parquet <- tempfile(fileext = ".parquet")
arrow::write_parquet(df, tmp_parquet)
message("Synthetic joined_data written to ", tmp_parquet, " (", n_counties, " counties, ",
        length(ht_cols), " HT* bins)\n")

# ---- Call the handler exactly as spade would ----
joined_data <- File(path = tmp_parquet)
result <- handler(joined_data = joined_data, response = "ESTIMATE",
                   alphas = "0,0.25,0.5,0.75,1", seed = 42)

message("\nDone. Output file object: ", result@path, "\n")
message("Open it and confirm the page count/content matches the checklist above.\n")
