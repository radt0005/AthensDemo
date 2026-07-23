# Block: join_chm_stats

# test to see if the segfault goes away with this setting
Sys.setenv(ARROW_USER_SIMD_LEVEL = "AVX2")

message("Loading package libraries...\n")
library(spade)
library(terra)
library(sf)
library(exactextractr)
library(dplyr)
library(tidyr)
library(arrow)

arrow::set_cpu_count(1)
arrow::set_io_thread_count(2)

message("Finished loading package libraries...\n")

handler <- function(fia_data,
                    state_abbrev = "VA",
                    gfchm_file,
                    state_file,
                    county_file){

  message("Reading reformatted FIA estimates...\n")
  fia_df <- readRDS(fia_data@path)

  # API may return duplicate columns for STATECD...
  if(length(grep("STATECD",names(fia_df)))> 1){
    fia_df = fia_df[,-grep("STATECD",names(fia_df))[2]]
  }

  message("Reading county GFCHM file...\n")
  fn_gfchm <- gfchm_file@path
  message(paste("GFCHM raster filenames for",length(fn_gfchm),"counties in",state_abbrev))
  county_rast = rast(fn_gfchm)

  message("Reading state codes...\n")
  state_codes <- read.csv(state_file@path, stringsAsFactors = FALSE)
  state_codes <- state_codes %>%
    rename(STATECD = STATE,
           STATEAB = STUSAB,
           STATENAME = STATE_NAME) %>% arrange(STATECD) %>% distinct()
  state_fips_int <- state_codes$STATECD[which(state_codes$STATEAB == state_abbrev)]
  message("State codes read from file...\n")

  message("Reading county codes...\n")
  county_codes <- read.csv(county_file@path, stringsAsFactors = FALSE)
  county_codes <- county_codes %>% filter(STATECD == state_fips_int) %>%
    dplyr::select(STATECD,UNITCD,COUNTYCD,COUNTYNM) %>%
    arrange(COUNTYCD)
  message("County codes read from file...\n")

  message("State:", state_abbrev, " - Counties: ", nrow(county_codes), "\n")

  countynames = sort(county_codes$COUNTYNM)

  match_idx <- which(
    sapply(countynames, function(x)
      grepl(x, fn_gfchm, ignore.case = TRUE))
  )
  county_name <- countynames[match_idx]

  county_fips_int = county_codes$COUNTYCD[county_codes$COUNTYNM == county_name]

  county_name2 = strsplit(fn_gfchm,split='_')[[1]][2]
  message("State:", state_abbrev, " - County: ", county_name, "\n")
  message("Raster file for:", county_name2, " county read into memory.\n")

  if(!grepl(toupper(county_name),toupper(fn_gfchm)))
    stop("County name match error: join_chm_stats.R")
  county_rast$Layer_1[county_rast$Layer_1 > 100] <- NA

  bins <- (cut(values(county_rast$Layer_1),breaks=c(-1,0,5,10,15,20,25,30,35),labels=F)-1)*5
  county_bins = county_rast$Layer_1
  values(county_bins) <- bins

  chm_dist <- as.data.frame(table(bins))
  chm_dist$STATECD = state_fips_int
  chm_dist$COUNTYCD = county_fips_int
  chm_dist$COUNTYNM = county_name
  chm_dist$BIN_HT = as.numeric(as.character(chm_dist$bins))
  chm_dist$KM2 <- chm_dist$Freq * prod(res(county_rast))/1e6
  # chm_dist$AREA_FRAC = chm_dist$KM2 / sum(chm_dist$KM2)
  chm_dist$AREA = chm_dist$KM2

  chm_dist <- chm_dist %>% left_join(state_codes) %>%
    dplyr::select(-STATENS)

  message("Pivoting CHM height-bin distribution wide (one column per bin)...\n")
  chm_wide <- chm_dist %>%
    dplyr::mutate(CO_FIPS = STATECD * 1000 + COUNTYCD) %>%
    dplyr::select(STATECD, COUNTYCD, CO_FIPS, BIN_HT, AREA) %>%
    tidyr::pivot_wider(
      names_from = BIN_HT,
      values_from = AREA,
      names_prefix = "HT",
      values_fill = 0
    )

  message("Joining CHM height-bin columns onto FIA estimates by CO_FIPS...\n")
  result <- fia_df %>%
    dplyr::right_join(chm_wide, by = c("STATECD", "COUNTYCD", "CO_FIPS"))

  arrow::write_parquet(result, "result.parquet")
  File(path = "result.parquet")
}
spade_types(handler) <- list(
  fia_data    = "File",
  gfchm_file  = "File",
  state_file  = "File",
  county_file = "File",
  .return     = "File"
)

run(handler)
