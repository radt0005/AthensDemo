# Block: fia_gfchm_diagnostics

message("Loading package libraries...\n")
library(spade)
library(arrow)
library(dplyr)

message("Finished loading package libraries...\n")

handler <- function(joined_data,
                     id_cols = "STATECD,COUNTYCD,CO_FIPS,COUNTYNM,STATEAB,STATENAME,UNITCD,UNITNM,EVAL_GRP,ATTRIBUTE_NBR",
                     fia_cols = ""){

  message("Reading joined FIA/GFCHM parquet file...\n")
  df <- arrow::read_parquet(joined_data@path)

  id_col_vec <- trimws(strsplit(id_cols, ",")[[1]])
  id_col_vec <- id_col_vec[id_col_vec %in% names(df)]

  ht_col_vec <- grep("^HT", names(df), value = TRUE)
  if (length(ht_col_vec) == 0)
    stop("No HT* (GFCHM height-bin) columns found in joined_data.")

  if (nzchar(fia_cols)) {
    fia_col_vec <- trimws(strsplit(fia_cols, ",")[[1]])
    missing <- setdiff(fia_col_vec, names(df))
    if (length(missing) > 0)
      stop(paste("fia_cols not found in joined_data:", paste(missing, collapse = ", ")))
  } else {
    is_numeric <- sapply(df, is.numeric)
    candidate_cols <- names(df)[is_numeric]
    fia_col_vec <- setdiff(candidate_cols, c(id_col_vec, ht_col_vec))
  }

  if (length(fia_col_vec) == 0)
    stop("No FIA estimate columns identified (after excluding id_cols and HT* columns). Pass fia_cols explicitly.")

  message("FIA estimate columns: ", paste(fia_col_vec, collapse = ", "), "\n")
  message("GFCHM height-bin columns: ", paste(ht_col_vec, collapse = ", "), "\n")

  fia_mat <- df[, fia_col_vec, drop = FALSE]
  ht_mat  <- df[, ht_col_vec, drop = FALSE]

  message("Computing FIA x GFCHM cross-correlation matrix...\n")
  cor_mat <- matrix(NA_real_, nrow = length(fia_col_vec), ncol = length(ht_col_vec),
                     dimnames = list(fia_col_vec, ht_col_vec))
  for (i in seq_along(fia_col_vec)) {
    for (j in seq_along(ht_col_vec)) {
      cor_mat[i, j] <- suppressWarnings(
        cor(fia_mat[[i]], ht_mat[[j]], use = "pairwise.complete.obs")
      )
    }
  }

  message("Rendering diagnostic PDF (table + heatmap + scatterplot grid)...\n")
  pdf("diagnostics.pdf", width = 10, height = 8)

  # --- Page 1: correlation matrix as a printed table ---
  plot.new()
  title("FIA x GFCHM Pearson Correlation Matrix")
  tbl_text <- capture.output(print(round(cor_mat, 3)))
  text(0.02, 0.95, paste(tbl_text, collapse = "\n"),
       adj = c(0, 1), family = "mono", cex = 0.8)

  # --- Page 2: correlation heatmap ---
  n_fia <- nrow(cor_mat); n_ht <- ncol(cor_mat)
  op <- par(mar = c(6, 8, 4, 2))
  image(1:n_ht, 1:n_fia, t(cor_mat[n_fia:1, , drop = FALSE]),
        col = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
        zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "",
        main = "Correlation: FIA estimates vs GFCHM height-bin area fractions")
  axis(1, at = 1:n_ht, labels = colnames(cor_mat), las = 2, cex.axis = 0.8)
  axis(2, at = 1:n_fia, labels = rev(rownames(cor_mat)), las = 2, cex.axis = 0.8)
  for (i in 1:n_fia) {
    for (j in 1:n_ht) {
      val <- cor_mat[n_fia - i + 1, j]
      text(j, i, sprintf("%.2f", val), cex = 0.7)
    }
  }
  par(op)

  # --- Page(s): scatterplot grid, FIA vars (rows) x HT bins (cols) ---
  max_panels_per_page <- 20
  ht_per_page <- max(1, floor(max_panels_per_page / n_fia))
  ht_chunks <- split(ht_col_vec, ceiling(seq_along(ht_col_vec) / ht_per_page))

  for (chunk in ht_chunks) {
    op <- par(mfrow = c(n_fia, length(chunk)), mar = c(4, 4, 2, 1))
    for (fcol in fia_col_vec) {
      for (hcol in chunk) {
        r <- cor_mat[fcol, hcol]
        plot(df[[hcol]], df[[fcol]],
             xlab = hcol, ylab = fcol,
             main = sprintf("r = %.2f", r),
             pch = 16, col = adjustcolor("steelblue", alpha.f = 0.6),
             cex.main = 0.9, cex.lab = 0.8)
        fit <- lm(df[[fcol]] ~ df[[hcol]])
        abline(fit, col = "firebrick", lwd = 1.5)
      }
    }
    par(op)
  }

  dev.off()

  message("Diagnostics written: diagnostics.pdf\n")
  File(path = "diagnostics.pdf")
}

spade_types(handler) <- list(
  joined_data = "File",
  .return     = "File"
)

run(handler)
