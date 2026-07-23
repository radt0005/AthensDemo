# Block: glmnet_var_select

message("Loading package libraries...\n")
library(spade)
library(arrow)
library(dplyr)
library(glmnet)

message("Finished loading package libraries...\n")

handler <- function(joined_data,
                    response = "ESTIMATE",
                    alphas = "0,0.25,0.5,0.75,1",
                    seed = 42){

  SPADE = TRUE # set to FALSE to run outside of spade

  message("Reading joined FIA/GFCHM parquet file...\n")
  if(SPADE){
    df <- arrow::read_parquet(joined_data@path)
  }else{
    # for running outside of spade
    result_fn = "~/.spade/pipelines/019f7827-31fc-7bb3-be77-eb3daf03976d/019f636f-7e74-78af-8339-235ea5ec7393/outputs/result/result.parquet"
    df <- read_parquet(result_fn)
  }

  # fill in zeros for missing values in height bins
  height_cols <- grep("^HT", names(df))
  df[height_cols] <- lapply(df[height_cols], function(x) {
    x[is.na(x)] <- 0
    x
  })

  if (!(response %in% names(df)))
    stop(paste0("response column '", response, "' not found in joined_data"))

  ht_cols <- grep("^HT", names(df), value = TRUE)
  if (length(ht_cols) == 0)
    stop("No HT* (GFCHM height-bin) columns found in joined_data.")

  message("Predictors: ", paste(ht_cols, collapse = ", "), "\n")
  message("Response: ", response, "\n")

  X <- as.matrix(df[, ht_cols, drop = FALSE])
  y <- df[[response]]

  keep <- stats::complete.cases(X, y)
  if (any(!keep)) {
    message("Dropping ", sum(!keep), " row(s) with missing values.\n")
    X <- X[keep, , drop = FALSE]
    y <- y[keep]
  }

  alpha_vec <- as.numeric(trimws(strsplit(alphas, ",")[[1]]))
  set.seed(seed)

  message("Fitting cv.glmnet across alpha = ", paste(alpha_vec, collapse = ", "), "...\n")
  # Shared fold assignment across all alphas so CV-MSE comparisons are
  # apples-to-apples (otherwise each alpha gets its own random folds).
  n <- length(y)
  nfolds <- min(10, n)
  foldid <- sample(rep(seq_len(nfolds), length.out = n))

  cv_fits <- setNames(vector("list", length(alpha_vec)), paste0("alpha_", alpha_vec))
  for (i in seq_along(alpha_vec)) {
    cv_fits[[i]] <- cv.glmnet(X, y, alpha = alpha_vec[i], foldid = foldid, intercept = FALSE)
  }

  message("Fitting matching intercept = TRUE models (same folds) for intercept diagnostics...\n")
  # Same alpha grid, same foldid, only intercept differs, so any change in
  # CV-MSE below is attributable to the intercept choice and not to
  # different fold assignments or a different alpha search.
  cv_fits_int <- setNames(vector("list", length(alpha_vec)), paste0("alpha_", alpha_vec))
  for (i in seq_along(alpha_vec)) {
    cv_fits_int[[i]] <- cv.glmnet(X, y, alpha = alpha_vec[i], foldid = foldid, intercept = TRUE)
  }

  message("Building coefficient comparison table (lambda.min and lambda.1se per alpha)...\n")
  coef_cols <- list()
  for (i in seq_along(alpha_vec)) {
    a <- alpha_vec[i]
    fit <- cv_fits[[i]]
    coef_cols[[paste0("a", a, "_min")]] <- as.numeric(coef(fit, s = "lambda.min"))[-1]
    coef_cols[[paste0("a", a, "_1se")]] <- as.numeric(coef(fit, s = "lambda.1se"))[-1]
  }
  coef_table <- do.call(cbind, coef_cols)
  rownames(coef_table) <- ht_cols

  message("Ranking predictors by selection stability across the alpha grid...\n")
  # "Stable" = nonzero in every lambda.min/lambda.1se column across the
  # whole alpha grid (i.e. survives from ridge to lasso). Among stable
  # predictors, rank by mean standardized |coefficient| = |raw coef| *
  # sd(predictor) -- this puts bins with different natural variance on
  # a comparable footing for ranking purposes, without altering the fit
  # itself or the raw-km^2 coefficients shown in the table (glmnet's
  # own internal standardize = TRUE already fits fairly across scales;
  # this just makes the *reported ranking* reflect that same fairness).
  n_cols <- ncol(coef_table)
  sel_freq <- rowSums(coef_table != 0)
  predictor_sd <- apply(X, 2, sd)[rownames(coef_table)]
  std_coef_table <- coef_table * predictor_sd
  mean_abs_std_coef <- rowMeans(abs(std_coef_table))
  var_rank <- data.frame(
    predictor = rownames(coef_table),
    n_selected = sel_freq,
    n_total = n_cols,
    mean_abs_std_coef = mean_abs_std_coef
  )
  var_rank <- var_rank[order(-var_rank$n_selected, -var_rank$mean_abs_std_coef), ]

  stable_vars <- var_rank$predictor[var_rank$n_selected == n_cols]
  stable_vars_ranked <- stable_vars[order(-mean_abs_std_coef[stable_vars])]
  partial_vars <- var_rank$predictor[var_rank$n_selected > 0 & var_rank$n_selected < n_cols]

  var_recommendation <- if (length(stable_vars_ranked) > 0) {
    paste0(
      "Recommended predictors (nonzero across the full alpha grid, i.e. ",
      "ridge through lasso, at both lambda.min and lambda.1se), ranked by ",
      "mean standardized |coefficient| (coef x predictor SD): ",
      paste(stable_vars_ranked, collapse = ", "), ". ",
      if (length(partial_vars) > 0) {
        paste0("Weaker / selection-sensitive: ", paste(partial_vars, collapse = ", "), ".")
      } else {
        ""
      }
    )
  } else {
    "No predictor was retained across the full alpha grid -- selection is unstable; treat any single-alpha result with caution."
  }
  message(var_recommendation, "\n")

  message("Building intercept = TRUE vs FALSE diagnostic comparison...\n")
  # For each alpha, compare the CV-MSE achieved at that alpha's own
  # lambda.min under intercept = FALSE vs intercept = TRUE, and report
  # the fitted intercept term itself (as a raw value and as a % of the
  # mean response, so its practical size is easy to judge regardless of
  # the response's units).
  mean_y <- mean(y)
  int_diag_rows <- lapply(seq_along(alpha_vec), function(i) {
    a          <- alpha_vec[i]
    fit_false  <- cv_fits[[i]]
    fit_true   <- cv_fits_int[[i]]
    mse_false  <- fit_false$cvm[fit_false$lambda == fit_false$lambda.min]
    mse_true   <- fit_true$cvm[fit_true$lambda == fit_true$lambda.min]
    intercept_est <- as.numeric(coef(fit_true, s = "lambda.min"))[1]
    data.frame(
      alpha              = a,
      cvm_intercept_false = mse_false,
      cvm_intercept_true  = mse_true,
      pct_change_mse     = 100 * (mse_false - mse_true) / mse_false,
      intercept_est      = intercept_est,
      intercept_pct_mean_y = 100 * intercept_est / mean_y
    )
  })
  int_diag_table <- do.call(rbind, int_diag_rows)

  # Headline comparison: best achievable CV-MSE across the whole alpha
  # grid under each intercept regime (not just a single alpha), since
  # the best alpha can itself differ between the two regimes.
  best_false_idx <- which.min(sapply(cv_fits, function(f) min(f$cvm)))
  best_true_idx  <- which.min(sapply(cv_fits_int, function(f) min(f$cvm)))
  best_mse_false <- min(cv_fits[[best_false_idx]]$cvm)
  best_mse_true  <- min(cv_fits_int[[best_true_idx]]$cvm)
  overall_pct_change_mse <- 100 * (best_mse_false - best_mse_true) / best_mse_false
  best_true_intercept_est <- as.numeric(coef(cv_fits_int[[best_true_idx]], s = "lambda.min"))[1]
  best_true_intercept_pct_mean_y <- 100 * best_true_intercept_est / mean_y

  # Simple, transparent thresholds -- not a formal test, just a
  # heuristic flag to draw the eye toward what the numbers already show.
  # A meaningfully lower CV-MSE with intercept = TRUE *and* a
  # non-trivial fitted intercept both need to hold before we'd suggest
  # abandoning the zero-intercept (physically-motivated) model.
  mse_threshold  <- 2   # % improvement in CV-MSE
  int_threshold  <- 5   # % of mean response
  meaningful_mse_gain <- overall_pct_change_mse > mse_threshold
  meaningful_intercept <- abs(best_true_intercept_pct_mean_y) > int_threshold

  recommendation <- if (meaningful_mse_gain && meaningful_intercept) {
    paste0(
      "Recommendation: consider intercept = TRUE. Allowing an intercept ",
      "lowers the best achievable CV-MSE by ", sprintf("%.1f", overall_pct_change_mse),
      "%, and the fitted intercept (", sprintf("%.3g", best_true_intercept_est),
      ") is ", sprintf("%.1f", abs(best_true_intercept_pct_mean_y)),
      "% of the mean response -- large enough that the zero-intercept ",
      "assumption may not fit this dataset well."
    )
  } else if (meaningful_mse_gain && !meaningful_intercept) {
    paste0(
      "Recommendation: intercept = FALSE is still reasonable. CV-MSE improves ",
      "slightly with an intercept (", sprintf("%.1f", overall_pct_change_mse),
      "%), but the fitted intercept is small relative to the response (",
      sprintf("%.1f", abs(best_true_intercept_pct_mean_y)),
      "% of the mean) -- the gain is likely noise rather than a real ",
      "nonzero-intercept relationship."
    )
  } else {
    paste0(
      "Recommendation: intercept = FALSE (the zero-intercept model) is ",
      "supported by these results. CV-MSE does not improve meaningfully ",
      "when an intercept is added (", sprintf("%.1f", overall_pct_change_mse),
      "% change), consistent with the Cao et al. zero-canopy / zero-total ",
      "physical justification for this predictor set."
    )
  }
  message(recommendation, "\n")

  message("Rendering diagnostic PDF (selection table + CV curves + coefficient paths)...\n")
  pdf("glmnet_diagnostics.pdf", width = 10, height = 8)

  # --- Page 1: selected-coefficient table across alphas ---
  plot.new()
  title(paste0(
    "glmnet Variable Selection (lambda.min / lambda.1se) - response: ",
    response
  ))
  # title(paste0("glmnet Variable Selection (lambda.min / lambda.1se) \u2014 response: ", response))
  tbl_text <- capture.output(print(round(coef_table, 4)))
  full_p1_text <- paste(c(
    tbl_text,
    "",
    strwrap(var_recommendation, width = 100)
  ), collapse = "\n")
  text(0.02, 0.95, full_p1_text,
       adj = c(0, 1), family = "mono", cex = 0.7)

  # --- Page 2: intercept = TRUE vs FALSE diagnostic ---
  plot.new()
  title("Intercept Diagnostic: TRUE vs FALSE")
  int_tbl_text <- capture.output(print(round(int_diag_table, 4), row.names = FALSE))
  header_lines <- c(
    "Per-alpha cross-validated MSE (CV-MSE) at each fit's own lambda.min,",
    "intercept = FALSE vs TRUE (same folds both ways, so differences",
    "reflect the intercept choice only):",
    ""
  )
  footer_lines <- c(
    "",
    sprintf("Best achievable CV-MSE, intercept = FALSE: %.4g (alpha = %.2f)",
            best_mse_false, alpha_vec[best_false_idx]),
    sprintf("Best achievable CV-MSE, intercept = TRUE:  %.4g (alpha = %.2f, intercept = %.4g, %.1f%% of mean response)",
            best_mse_true, alpha_vec[best_true_idx], best_true_intercept_est,
            best_true_intercept_pct_mean_y),
    sprintf("Change in best CV-MSE from adding an intercept: %.1f%%", overall_pct_change_mse),
    "",
    strwrap(recommendation, width = 95)
  )
  full_text <- paste(c(header_lines, int_tbl_text, footer_lines), collapse = "\n")
  text(0.02, 0.95, full_text, adj = c(0, 1), family = "mono", cex = 0.65)

  # --- Page: predictor color key ---
  # Colors are assigned once and reused for every alpha's coefficient-path
  # plot below, so a single legend page covers all of them instead of a
  # cramped, easily-clipped legend on every panel.
  predictor_colors <- setNames(
    grDevices::palette.colors(n = length(ht_cols), palette = "Okabe-Ito"),
    ht_cols
  )
  plot.new()
  title("Predictor Color Key (used on all coefficient-path plots below)")
  legend("center", legend = ht_cols, col = predictor_colors, lty = 1, lwd = 3,
         bty = "n", cex = 1.1, ncol = 2)

  # --- Page(s): CV curve + coefficient path per alpha ---
  for (i in seq_along(alpha_vec)) {
    a <- alpha_vec[i]
    fit_cv <- cv_fits[[i]]
    fit_full <- glmnet(X, y, alpha = a, intercept = FALSE)

    op <- par(mfrow = c(1, 2))

    # Manual CV-MSE plot: same data as plot.cv.glmnet (error bars, red
    # points, lambda.min/lambda.1se dashed lines) but without the top
    # "# nonzero coefficients" axis, which plot.cv.glmnet always draws
    # and offers no argument to suppress.
    log_lambda <- log(fit_cv$lambda)
    plot(log_lambda, fit_cv$cvm, type = "n",
         xlab = "Log(Lambda)", ylab = "Mean-Squared Error",
         ylim = range(fit_cv$cvup, fit_cv$cvlo))
    segments(log_lambda, fit_cv$cvlo, log_lambda, fit_cv$cvup, col = "darkgrey")
    points(log_lambda, fit_cv$cvm, pch = 20, col = "red")
    abline(v = log(fit_cv$lambda.min), lty = 3)
    abline(v = log(fit_cv$lambda.1se), lty = 3)
    title(sprintf("CV-MSE, alpha = %.2f", a), line = 2.5)

    # Manual coefficient-path plot: same data as plot.glmnet(xvar =
    # "lambda") but colored/labeled by predictor name via the shared
    # legend page above, instead of glmnet's clipped numeric index labels.
    beta_mat <- as.matrix(fit_full$beta)          # rows = ht_cols, cols = lambda steps
    log_lambda_path <- log(fit_full$lambda)
    matplot(log_lambda_path, t(beta_mat), type = "l", lty = 1, lwd = 1.5,
            col = predictor_colors[rownames(beta_mat)],
            xlab = "Log(Lambda)", ylab = "Coefficients")
    abline(v = log(fit_cv$lambda.min), lty = 3)
    abline(v = log(fit_cv$lambda.1se), lty = 3)
    title(sprintf("Coefficient paths, alpha = %.2f", a), line = 2.5)

    par(op)
  }

  # --- Final page: guidance for interpreting this output ---
  plot.new()
  title("How to Interpret This Diagnostic Output")
  guidance <- c(
    "1. Coefficient table (page 1): compares which HT* bins are selected",
    "   (nonzero) at lambda.min (best CV-MSE) vs. lambda.1se (simplest model",
    "   within 1 SE of the best) across the alpha grid, from alpha = 0",
    "   (ridge, no selection) to alpha = 1 (lasso, sparse selection). A",
    "   predictor that stays nonzero across most alphas and both lambdas",
    "   is a stable, well-supported signal; one that only survives at",
    "   low alpha / lambda.min is a weaker, more selection-sensitive one.",
    "   The line below the table applies this rule automatically: it",
    "   first filters to predictors nonzero in every column (the",
    "   stability criterion), then ranks that filtered set by mean",
    "   standardized |coefficient| (coef x predictor SD) -- not raw",
    "   |coefficient| -- so a bin with more inherent variance across",
    "   counties isn't ranked higher purely because of its scale.",
    "",
    "2. Intercept diagnostic (page 2): this compares model fit with and",
    "   without a free intercept term, using identical folds so the only",
    "   thing that differs is the intercept choice.",
    "     - A small change in CV-MSE (roughly < 2%) means the intercept",
    "       is not doing meaningful work -- the zero-intercept model is",
    "       adequate and preferable for its physical interpretability",
    "       (zero canopy area implies zero volume).",
    "     - A meaningfully lower CV-MSE combined with a fitted intercept",
    "       that is large relative to the mean response (roughly > 5%)",
    "       suggests the zero-intercept assumption may not hold for this",
    "       dataset, and intercept = TRUE should be used instead.",
    "     - These thresholds are heuristic guides, not a formal",
    "       hypothesis test -- treat borderline cases as ambiguous and",
    "       prefer the physically-motivated model (intercept = FALSE)",
    "       unless the evidence for an intercept is clear-cut.",
    "",
    "3. CV-MSE curves and coefficient paths (following pages): the CV",
    "   curve's minimum marks lambda.min; the vertical dashed lines mark",
    "   lambda.min and lambda.1se. A CV curve that is flat over a wide",
    "   range of lambda suggests the model is insensitive to how hard it",
    "   is regularized -- a sign of either genuinely strong, redundant",
    "   predictors or limited sample size relative to predictor count.",
    "   The coefficient-path plot shows how each predictor's coefficient",
    "   shrinks toward zero as lambda increases; predictors that survive",
    "   to larger lambda (further right / higher on the shrinkage axis)",
    "   are the more dominant signals. Colors are consistent across every",
    "   alpha and match the predictor color key page just before these.",
    "",
    "4. Predictor scale: HT0...HT35 are raw per-county areas (km^2) in",
    "   each canopy-height class, not fractions of a county's total",
    "   forested area. This is what makes the zero-intercept model",
    "   physically sensible in the first place (zero canopy area implies",
    "   zero total volume) and matches the Cao et al. methodology --",
    "   the intercept diagnostic on page 2 is testing that assumption",
    "   directly against this predictor set, not against a compositional",
    "   (fraction-based) one."
  )
  text(0.02, 0.95, paste(guidance, collapse = "\n"),
       adj = c(0, 1), family = "mono", cex = 0.62)

  dev.off()

  message("Diagnostics written: glmnet_diagnostics.pdf\n")
  File(path = "glmnet_diagnostics.pdf")
}

spade_types(handler) <- list(
  joined_data = "File",
  .return     = "File"
)

run(handler)
