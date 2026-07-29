# Block: fh_sae

message("Loading package libraries...\n")
library(spadelib)
library(arrow)
library(dplyr)
library(sae)

message("Finished loading package libraries...\n")

# Round numeric columns of a data.frame for compact printing in the PDF
# preview tables, without touching id/character columns.
round_df <- function(d, digits = 3) {
  num_cols <- vapply(d, is.numeric, logical(1))
  d[num_cols] <- lapply(d[num_cols], round, digits = digits)
  d
}

handler <- function(joined_data,
                    response             = "ESTIMATE",
                    variance             = "VARIANCE",
                    predictors           = "HT20,HT25,HT30,HT35",
                    id_cols              = "STATECD,UNITCD,UNITNM,COUNTYCD,CO_FIPS,COUNTY_NAME",
                    fia_se_col           = "SE",
                    fia_plot_count_col   = "PLOT_COUNT",
                    fia_nonzero_plot_col = "NON_ZERO_PLOTS",
                    min_plot_count       = 1,
                    intercept            = "TRUE",
                    method               = "REML"){

  SPADE = TRUE # set to FALSE to run outside of spade

  message("Reading joined FIA/GFCHM parquet file...\n")
  if(SPADE){
    df <- arrow::read_parquet(joined_data@path)
  }else{
    # for running outside of spade
    result_fn = "~/.spade/pipelines/PLACEHOLDER/outputs/result/result.parquet"
    df <- read_parquet(result_fn)
  }

  # sae::mseFH()/eblupFH() internally do `vardir <- data[, namevar]`
  # (single-bracket, no drop = FALSE). On a base data.frame that drops to
  # a plain atomic vector, which is what those functions expect -- but
  # arrow::read_parquet() returns a tibble, and tibbles' `[` never drops,
  # so `data[, namevar]` would instead return a one-column tibble and
  # break sae::'s internal median()/arithmetic on it ("need numeric
  # data"). Convert to a base data.frame up front so sae:: gets the plain
  # vector it's written to expect.
  df <- as.data.frame(df)

  id_cols_vec    <- trimws(strsplit(id_cols, ",")[[1]])
  predictor_cols <- trimws(strsplit(predictors, ",")[[1]])

  # Required for the model itself -- fail loudly if these are missing,
  # since there's no sensible fallback.
  required_cols <- c(response, variance, fia_se_col, predictor_cols)
  missing_required <- setdiff(required_cols, names(df))
  if (length(missing_required) > 0)
    stop("joined_data is missing required column(s): ",
         paste(missing_required, collapse = ", "))

  missing_id <- setdiff(id_cols_vec, names(df))
  if (length(missing_id) > 0) {
    message("Note: id_cols not found in joined_data (dropped): ",
            paste(missing_id, collapse = ", "), "\n")
    id_cols_vec <- setdiff(id_cols_vec, missing_id)
  }

  # Optional FIA sampling-effort columns. NON_ZERO_PLOTS in particular is
  # currently commented out of reformat_fia.R's dplyr::select(), so it
  # won't be present until that's uncommented upstream -- don't fail the
  # whole block over an optional descriptive column, just flag it.
  optional_cols <- c(fia_plot_count_col, fia_nonzero_plot_col)
  for (col in optional_cols) {
    if (!(col %in% names(df))) {
      message("Note: optional column '", col, "' not found in joined_data -- ",
              "filling with NA. If this is PLOT_COUNT or NON_ZERO_PLOTS, check ",
              "whether reformat_fia.R's column selection includes it.\n")
      df[[col]] <- NA_real_
    }
  }

  # Drop counties with unreliably few FIA field plots behind their direct
  # estimate (e.g. singleton-plot counties) before fitting the Fay-Herriot
  # model -- these estimates/variances are the least trustworthy inputs to
  # mseFH() and shouldn't be allowed to distort the fit. Uses the same
  # fia_plot_count_col/min_plot_count convention as glmnet_var_select.R;
  # point both blocks at the same threshold to fit/select on the same
  # county set. Rows where the column is genuinely absent (all-NA, i.e.
  # the "optional column not found" case just above) are left alone here --
  # there's nothing to filter on -- but a present column with per-row NA
  # plot counts still has those rows dropped, since an unknown plot count
  # isn't evidence the estimate is reliable.
  if (!all(is.na(df[[fia_plot_count_col]]))) {
    n_before <- nrow(df)
    low_plot <- is.na(df[[fia_plot_count_col]]) | df[[fia_plot_count_col]] <= min_plot_count
    if (any(low_plot)) {
      message("Dropping ", sum(low_plot), " of ", n_before, " county row(s) with ",
              fia_plot_count_col, " <= ", min_plot_count,
              " (or missing) -- unreliable direct estimates excluded from the FH fit.\n")
      df <- df[!low_plot, , drop = FALSE]
    }
  } else {
    message("Note: '", fia_plot_count_col, "' not available -- skipping ",
            "low-plot-count filter (min_plot_count = ", min_plot_count, " had no effect).\n")
  }

  # Fill in zeros for missing predictor values, consistent with
  # glmnet_var_select.R's convention for these HT* columns.
  df[predictor_cols] <- lapply(df[predictor_cols], function(x) {
    x[is.na(x)] <- 0
    x
  })

  keep <- stats::complete.cases(df[, c(response, variance, predictor_cols), drop = FALSE])
  if (any(!keep)) {
    message("Dropping ", sum(!keep), " county row(s) with missing response/variance/predictor values.\n")
    df <- df[keep, , drop = FALSE]
  }

  n_domains <- nrow(df)

  # Parse the intercept flag up front so a bad value (anything as.logical()
  # can't resolve to TRUE/FALSE) fails loudly here rather than silently
  # producing an NA-riddled formula downstream.
  intercept_flag <- as.logical(intercept)
  if (is.na(intercept_flag))
    stop("intercept must be \"TRUE\" or \"FALSE\" (got: '", intercept, "').")

  message("Fitting Fay-Herriot model for ", n_domains, " counties: ", response,
          " ~ ", paste(predictor_cols, collapse = " + "),
          if (!intercept_flag) " (no intercept)" else "", "\n")

  # intercept = FALSE fits the zero-intercept model (formula ~ predictors - 1),
  # matching glmnet_var_select.R's intercept = FALSE default/diagnostic --
  # review that block's diagnostic PDF (page 2, "Intercept Diagnostic") before
  # setting this. Its recommendation text there is the evidence to act on:
  # only switch to intercept = TRUE if it flags a meaningful CV-MSE gain
  # *and* a non-trivial fitted intercept; otherwise the zero-intercept model
  # (intercept = FALSE) is the one its own diagnostics support.
  fh_formula <- if (intercept_flag) {
    stats::as.formula(paste(response, "~", paste(predictor_cols, collapse = " + ")))
  } else {
    stats::as.formula(paste(response, "~", paste(predictor_cols, collapse = " + "), "- 1"))
  }

  # mseFH() fits the area-level Fay-Herriot model (via eblupFH() internally)
  # and returns both the EBLUP and its (analytic, Prasad-Rao/Yokum-type)
  # MSE in one call.
  #
  # sae::mseFH()/eblupFH() resolve `vardir` via `deparse(substitute(vardir))`
  # -- i.e. they capture the literal *text* of whatever expression is
  # passed for `vardir`, then re-look that text up as a column name in
  # `data`. That means passing an expression like `df[[variance]]` fails,
  # because the literal text "df[[variance]]" isn't a real column name --
  # it doesn't matter that the expression *evaluates* to the right vector,
  # since eblupFH() never uses that evaluated value when `data` is
  # supplied. `do.call()` sidesteps this: passing an actual symbol object
  # (`as.symbol(variance)`, e.g. the symbol `VARIANCE`) means the
  # constructed call's `vardir` argument *is* that symbol, so
  # `substitute(vardir)` returns it directly and deparses to the real
  # column name.
  fh_fit <- do.call(sae::mseFH, list(
    formula = fh_formula,
    vardir  = as.symbol(variance),
    method  = method,
    data    = df
  ))

  # mseFH() doesn't error on non-convergence -- it warns and returns
  # early with $eblup/$mse left as their NA initial values, which would
  # otherwise surface downstream as a confusing length-mismatch error
  # rather than a clear one.
  if (isFALSE(fh_fit$est$fit$convergence))
    stop("sae::mseFH did not converge (method = '", method, "'). Try a ",
         "different method, or review the predictor set for near-collinearity.")

  # eblupFH() returns $eblup as the point-estimate vector directly (in
  # practice a 1-column matrix, from X %*% beta arithmetic inside sae::),
  # not a data.frame with an "eblup" column -- as.numeric() flattens it
  # to a plain vector for arithmetic/plotting below.
  eblup   <- as.numeric(fh_fit$est$eblup)
  fh_mse  <- fh_fit$mse
  fh_rmse <- sqrt(fh_mse)

  # Standard-error-ratio (SER): RMSE of the FH EBLUP over the SE of the
  # FIA direct estimate, per county. SER < 1 means the model-based EBLUP
  # is more precise than the design-based direct estimate.
  se_direct <- df[[fia_se_col]]
  ser <- fh_rmse / se_direct

  # Percent-scale precision columns are both standardized by the FH EBLUP,
  # not by each series' own point estimate (i.e. not FIA SE% of ESTIMATE
  # and EBLUP RMSE% of EBLUP separately). Two reasons: FIA direct
  # estimates can be unstable in low-plot-count counties, which would
  # inflate the FIA side of the comparison for reasons unrelated to
  # actual precision; and using two different denominators makes the
  # pair harder to interpret side by side. Standardizing both by the more
  # stable EBLUP keeps them on one consistent, shared scale -- and, as a
  # side effect, makes the 1:1 line on that plot exactly the SER = 1
  # boundary, since the EBLUP denominator cancels in the ratio either way.
  fia_se_pct_of_eblup <- 100 * se_direct / eblup
  fh_rmse_pct_of_eblup <- 100 * fh_rmse / eblup

  results <- df[, id_cols_vec, drop = FALSE]
  results[[response]]             <- df[[response]]
  results[[fia_se_col]]           <- se_direct
  results[["SE_PERCENT"]]         <- if ("SE_PERCENT" %in% names(df)) {
    df[["SE_PERCENT"]]
  } else {
    100 * se_direct / df[[response]]
  }
  results[[fia_plot_count_col]]   <- df[[fia_plot_count_col]]
  results[[fia_nonzero_plot_col]] <- df[[fia_nonzero_plot_col]]
  results$FH_EBLUP                <- eblup
  results$FH_RMSE                 <- fh_rmse
  results$FIA_SE_PCT_OF_EBLUP     <- fia_se_pct_of_eblup
  results$FH_RMSE_PCT_OF_EBLUP    <- fh_rmse_pct_of_eblup
  results$SER                     <- ser

  n_improved <- sum(ser < 1, na.rm = TRUE)
  pct_improved <- 100 * n_improved / n_domains
  summary_msg <- sprintf(
    paste0("%d of %d counties (%.1f%%) have SER < 1 (FH EBLUP RMSE smaller ",
           "than the FIA direct SE) -- i.e. improved precision from small ",
           "area estimation."),
    n_improved, n_domains, pct_improved
  )
  message(summary_msg, "\n")

  csv_path <- "fh_sae_results.csv"
  write.csv(results, csv_path, row.names = FALSE)
  message("Results table written: ", csv_path, "\n")

  # ---- Model fit diagnostics: best-effort. sae::'s internal fit-object
  # field names have shifted across package versions historically, so
  # degrade gracefully to a placeholder message instead of hard-failing
  # the whole block over a diagnostics-page-only formatting issue.
  coef_text <- tryCatch({
    paste(utils::capture.output(print(fh_fit$est$fit$estcoef)), collapse = "\n")
  }, error = function(e) "Coefficient table unavailable from this sae:: package version.")

  refvar_text <- tryCatch({
    sprintf("Estimated area-level (random-effect) variance: %.4g", fh_fit$est$fit$refvar)
  }, error = function(e) "")

  goodness_text <- tryCatch({
    paste(utils::capture.output(print(fh_fit$est$fit$goodness)), collapse = "\n")
  }, error = function(e) "")

  message("Rendering diagnostic PDF...\n")
  pdf("fh_sae_diagnostics.pdf", width = 10, height = 8)

  # --- Page 1: model summary ---
  plot.new()
  title(paste0("Fay-Herriot EBLUP Model Summary - response: ", response))
  p1_text <- paste(c(
    paste0("Formula: ", deparse(fh_formula)),
    paste0("Method: ", method, "   Domains (counties): ", n_domains,
           "   Intercept: ", intercept_flag),
    paste0("min_plot_count filter: ", fia_plot_count_col, " > ", min_plot_count),
    "",
    "Coefficients:",
    coef_text,
    "",
    refvar_text,
    "",
    "Goodness of fit:",
    goodness_text
  ), collapse = "\n")
  text(0.02, 0.95, p1_text, adj = c(0, 1), family = "mono", cex = 0.7)

  # --- Page 2: SER summary + best/worst-county preview tables ---
  plot.new()
  title("Standard Error Ratio (SER) Summary")
  ord <- order(results$SER)
  preview_n <- min(15, nrow(results))
  best  <- results[ord[seq_len(preview_n)], ]
  worst <- results[rev(ord)[seq_len(preview_n)], ]
  preview_cols <- intersect(
    c(id_cols_vec, response, fia_se_col, "FH_EBLUP", "FH_RMSE", "SER"),
    names(results)
  )
  p2_text <- paste(c(
    summary_msg,
    "",
    sprintf("Median SER: %.3f   Min: %.3f   Max: %.3f",
            stats::median(results$SER, na.rm = TRUE),
            min(results$SER, na.rm = TRUE),
            max(results$SER, na.rm = TRUE)),
    "",
    "Most improved counties (lowest SER):",
    utils::capture.output(print(round_df(best[, preview_cols, drop = FALSE]), row.names = FALSE)),
    "",
    "Least improved counties (highest SER):",
    utils::capture.output(print(round_df(worst[, preview_cols, drop = FALSE]), row.names = FALSE)),
    "",
    "Full results for all counties are in fh_sae_results.csv."
  ), collapse = "\n")
  text(0.02, 0.95, p2_text, adj = c(0, 1), family = "mono", cex = 0.6)

  # --- Page 3: EBLUP vs. FIA direct estimate (bias check) ---
  plot(df[[response]], eblup,
       xlab = paste0("FIA direct estimate (", response, ")"),
       ylab = "Fay-Herriot EBLUP",
       main = "EBLUP vs. FIA Direct Estimate", pch = 19, col = "steelblue")
  abline(0, 1, lty = 2, col = "darkgrey")
  legend("topleft", legend = "1:1 line", lty = 2, col = "darkgrey", bty = "n")

  # --- Page 4: RMSE% vs. SE%, both standardized by the FH EBLUP ---
  plot(fia_se_pct_of_eblup, fh_rmse_pct_of_eblup,
       xlab = "FIA direct SE, % of FH EBLUP",
       ylab = "FH EBLUP RMSE, % of FH EBLUP",
       main = "Precision Comparison (both standardized by the FH EBLUP)",
       pch = 19, col = "darkorange")
  abline(0, 1, lty = 2, col = "darkgrey")
  legend("topleft",
         legend = c("1:1 line (SER = 1)", "points below: SER < 1 (FH more precise)"),
         lty = c(2, NA), col = c("darkgrey", NA), bty = "n")

  # --- Final page: guidance ---
  plot.new()
  title("How to Interpret This Diagnostic Output")
  guidance <- c(
    "1. Model summary (page 1): the fitted Fay-Herriot area-level model,",
    "   response ~ predictors, with vardir = the FIA direct estimate's",
    "   variance for each county. The coefficients and estimated",
    "   area-level (random-effect) variance describe how much variation",
    "   the predictors explain vs. what's left as area-level noise the",
    "   model borrows strength across counties to smooth over. The",
    "   Intercept line records whether this fit used a free intercept --",
    "   review glmnet_var_select's own intercept diagnostic (its PDF,",
    "   page 2) before choosing this rather than leaving it at the",
    "   default, and use the same choice consistently here.",
    "",
    "2. SER summary (page 2): SER = FH EBLUP RMSE / FIA direct SE, per",
    "   county. SER < 1 means the model-based EBLUP is more precise than",
    "   the design-based FIA direct estimate for that county -- the",
    "   headline result small area estimation is meant to deliver, and it",
    "   should hold for most Virginia counties, especially low-plot-count",
    "   ones where the direct estimate is least reliable. The preview",
    "   tables show the most- and least-improved counties; the full table",
    "   is in fh_sae_results.csv. Counties with fia_plot_count_col <=",
    "   min_plot_count were dropped before fitting (see page 1) since",
    "   their direct estimates/variances are the least trustworthy inputs",
    "   to the model.",
    "",
    "3. EBLUP vs. FIA direct estimate (page 3): checks for systematic",
    "   bias in the model-based estimates. Points should scatter evenly",
    "   around the 1:1 line; a consistent one-sided departure (e.g. EBLUPs",
    "   running high for low-estimate counties) would suggest the model",
    "   isn't capturing something the direct estimates are.",
    "",
    "4. Precision comparison (page 4): both series are standardized by",
    "   the FH EBLUP, not by each series' own point estimate -- so this",
    "   plot's 1:1 line is exactly the SER = 1 boundary, and points below",
    "   it are the counties with SER < 1. Standardizing both series by",
    "   the same, more stable denominator (rather than the FIA SE by the",
    "   FIA estimate and the EBLUP RMSE by the EBLUP separately) avoids",
    "   letting direct-estimate instability in small-sample counties",
    "   inflate one side of the comparison, and keeps both percentages on",
    "   a directly comparable footing.",
    "",
    "5. fh_sae_results.csv has the full per-county table -- FIA direct",
    "   stats (estimate, SE, SE%, plot count, nonzero plot count), the FH",
    "   EBLUP/RMSE, both EBLUP-standardized percent columns, and SER --",
    "   for sharing or use in downstream analysis."
  )
  text(0.02, 0.95, paste(guidance, collapse = "\n"),
       adj = c(0, 1), family = "mono", cex = 0.62)

  dev.off()

  message("Diagnostics written: fh_sae_diagnostics.pdf\n")

  list(
    result = File(path = "fh_sae_diagnostics.pdf"),
    table  = File(path = "fh_sae_results.csv")
  )
}

spade_types(handler) <- list(
  joined_data = "File",
  .return     = "File"
)

run(handler)
