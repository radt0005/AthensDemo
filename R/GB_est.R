# Block: GB_est

library(spade)
library(FIAapi)

handler <- function(EVAL_GRP = "512017",
                    ATTRIBUTE_NBR = 10,
                    GRP_BY_ATTRIB = "STATECD, COUNTYCD") {

  # Parse EVAL_GRP from comma-separated string to integer vector
  eval_grp_vec <- as.integer(trimws(strsplit(EVAL_GRP, ",")[[1]]))

  # Validate: no duplicate state codes (EVAL_GRP %/% 10000 = state FIPS)
  state_codes <- eval_grp_vec %/% 10000
  if (length(state_codes) != length(unique(state_codes))) {
    dup_states <- state_codes[duplicated(state_codes)]
    stop("Each state can appear at most once across EVAL_GRP values. ",
         "Duplicate state FIPS codes detected: ",
         paste(dup_states, collapse = ", "))
  }

  # Parse GRP_BY_ATTRIB
  grp_by_vec <- trimws(strsplit(GRP_BY_ATTRIB, ",")[[1]])

  # Call GB_api with vector of EVAL_GRP values
  result <- GB_api(
    EVAL_GRP      = eval_grp_vec,
    ATTRIBUTE_NBR = as.integer(ATTRIBUTE_NBR),
    GRP_BY_ATTRIB = grp_by_vec
  )

  saveRDS(result, "result.rds")
  File(path = "result.rds")
}

spade_types(handler) <- list(.return = "File")

run(handler)
