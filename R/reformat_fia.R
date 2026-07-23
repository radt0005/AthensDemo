# Block: reformat_fia
# Reformats GB_est output: selects columns, constructs CO_FIPS,
# joins UNITCD from reference CSV.

library(spade)
library(dplyr)

handler <- function(gb_est_result,
                    unit_ref) {

  # Read GB_est result
  result <- readRDS(gb_est_result@path)
  result <- result[, !duplicated(names(result)), drop = FALSE]
  # message("columns: ", paste(names(result), collapse = ", "))
  # message("duplicated: ", paste(names(result)[duplicated(names(result))], collapse = ", "))

  # Read survey unit reference
  unit_ref_df <- read.csv(unit_ref@path) %>% filter(UNITCD < 100) %>%
    dplyr::select(STATECD, COUNTYCD, UNITCD, UNITNM)

  # Construct CO_FIPS and select/reorder columns
  result <- result %>%
    dplyr::mutate(CO_FIPS = STATECD * 1000 + COUNTYCD) %>%
    dplyr::left_join(unit_ref_df, by = c("STATECD", "COUNTYCD")) %>%
    dplyr::select(
      # EVAL_GRP,
      STATECD,
      UNITCD,
      UNITNM,
      COUNTYCD,
      CO_FIPS,
      COUNTY_NAME,
      # ATTRIBUTE_NBR,
      ESTIMATE,
      VARIANCE,
      SE,
      SE_PERCENT,
      PLOT_COUNT
      # NON_ZERO_PLOTS,
      # TOT_POP_AC
    )

  saveRDS(result, "result.rds")
  File(path = "result.rds")
}

spade_types(handler) <- list(gb_est_result = "File", unit_ref = "File", .return = "File")

run(handler)
