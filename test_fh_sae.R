# Standalone end-to-end test harness for AthensDemo's fh_sae block.
#
# Same rationale as test_glmnet_var_select.R: this was written in a
# sandbox with no R interpreter and no CRAN/conda access, so the real
# sae::mseFH() call could not be executed there. Run this yourself to
# confirm the block works before wiring it into a real spade run --
# this is especially important for fh_sae since its sae::mseFH() usage
# (field names on the returned fit object: est$eblup$eblup, est$fit$
# estcoef, est$fit$refvar, est$fit$goodness, and the top-level $mse
# vector) was written from documentation/memory, not tested against a
# live install of the `sae` package.
#
# Usage: run from the AthensDemo repo root (so "R/fh_sae.R" resolves),
# with spadelib, arrow, dplyr, and sae already installed:
#
#   Rscript test_fh_sae.R
#
# It writes fh_sae_diagnostics.pdf and fh_sae_results.csv to the current
# working directory. Check:
#   - No errors/warnings from the sae::mseFH() call or the coefficient/
#     goodness-of-fit extraction (page 1 should show real numbers, not
#     the "unavailable from this sae:: package version" fallback text --
#     if you see that fallback, the field names in fh_sae.R's tryCatch
#     blocks need updating for your installed sae:: version).
#   - fh_sae_results.csv has one row per county with FIA stats, FH_EBLUP,
#     FH_RMSE, FIA_SE_PCT_OF_EBLUP, FH_RMSE_PCT_OF_EBLUP, and SER columns.
#   - The console message reports most counties with SER < 1 (expected,
#     since this synthetic data has a real area-level random effect on
#     top of a real predictor signal -- see below).
#   - Page 3 (EBLUP vs. FIA direct) scatters roughly along the 1:1 line.
#   - Page 4 (RMSE% vs. SE%, both standardized by the EBLUP) shows most
#     points below the 1:1 line, consistent with most SER < 1.

library(arrow)

block_file <- "R/fh_sae.R"
if (!file.exists(block_file)) {
  stop("Run this from the AthensDemo repo root -- ", block_file, " not found.")
}

# Source everything in the block file except the trailing `run(handler)`
# call, so we get the real `handler` function without needing a full
# spade params.yaml / inputs/ harness.
lines <- readLines(block_file)
run_line <- grep("^run\\(handler\\)\\s*$", lines)
if (length(run_line) != 1) {
  stop("Expected exactly one top-level `run(handler)` line in ", block_file,
       " -- found ", length(run_line), ". Check the file hasn't changed shape.")
}
eval(parse(text = paste(lines[-run_line], collapse = "\n")), envir = .GlobalEnv)

# ---- Build synthetic data consistent with the Fay-Herriot generative
# model (unlike a plain regression test, this needs a real area-level
# random effect on top of sampling error, since that's what
# sae::mseFH() is actually modeling) ----
set.seed(2)
n_counties <- 90
ht_cols <- paste0("HT", c(0, 5, 10, 15, 20, 25, 30, 35))
predictor_cols <- c("HT20", "HT25", "HT30", "HT35")

X <- sapply(ht_cols, function(nm) pmax(0, rgamma(n_counties, shape = 2, scale = 8)))
colnames(X) <- ht_cols

true_coef <- c(HT20 = 1.0, HT25 = 1.6, HT30 = 2.2, HT35 = 2.8)
mu <- as.numeric(X[, predictor_cols] %*% true_coef[predictor_cols])

# Area-level (random-effect) variance: the true "SAE" signal, i.e. the
# part of each county's true value that predictors don't explain.
sigma2_u <- (0.06 * mean(mu))^2
u <- rnorm(n_counties, sd = sqrt(sigma2_u))

# Sampling variance shrinks with plot count, like a real direct estimator.
plot_count <- pmax(5, round(rnorm(n_counties, mean = 25, sd = 10)))
nonzero_plot_count <- pmax(1, round(plot_count * runif(n_counties, 0.5, 0.9)))
vardir <- (0.20 * mu / sqrt(plot_count / 25))^2
e <- rnorm(n_counties, sd = sqrt(vardir))

ESTIMATE <- mu + u + e
SE <- sqrt(vardir)

df <- data.frame(
  STATECD = 51,
  UNITCD = 1,
  UNITNM = "Unit 1",
  COUNTYCD = seq_len(n_counties),
  CO_FIPS = 51000 + seq_len(n_counties),
  COUNTY_NAME = paste0("County", seq_len(n_counties)),
  ESTIMATE = ESTIMATE,
  VARIANCE = vardir,
  SE = SE,
  SE_PERCENT = 100 * SE / ESTIMATE,
  PLOT_COUNT = plot_count,
  NON_ZERO_PLOTS = nonzero_plot_count
)
df <- cbind(df, as.data.frame(X))

tmp_parquet <- tempfile(fileext = ".parquet")
arrow::write_parquet(df, tmp_parquet)
message("Synthetic joined_data written to ", tmp_parquet, " (", n_counties,
        " counties, area-level random effect sigma2_u = ", round(sigma2_u, 4), ")\n")

# ---- Call the handler exactly as spade would ----
joined_data <- File(path = tmp_parquet)
result <- handler(
  joined_data = joined_data,
  response = "ESTIMATE",
  variance = "VARIANCE",
  predictors = paste(predictor_cols, collapse = ",")
)

message("\nDone. Output file objects:\n  PDF: ", result$result@path,
        "\n  CSV: ", result$table@path, "\n")
message("Open both and confirm against the checklist above.\n")
