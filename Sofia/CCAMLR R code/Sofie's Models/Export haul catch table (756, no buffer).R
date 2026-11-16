# -----------------------------------------------------------------
# Load CCAMLR catch/effort data from subareas 756 AND 759 (no spatial
# buffer applied -- every haul with data is kept) and export a
# haul-level lat/long + catch table to Excel.
# -----------------------------------------------------------------

library(dplyr)
library(purrr)
library(writexl)

# -----------------------------
# CONFIGURATION
# -----------------------------
base_dir <- "F:/CCAMLR Data"
subareas <- c("756", "759")

out_dir  <- "C:/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/Sofia/CCAMLR R code/Sofie's Models"
out_file <- file.path(out_dir, "CCAMLR_756_759_haul_catch_table.xlsx")

# -----------------------------
# LOAD + COMBINE ALL .Rds FILES FOR EACH SUBAREA
# -----------------------------
load_subarea <- function(subarea) {
  rds_dir <- file.path(base_dir, subarea, "R")
  rds_files <- list.files(rds_dir, pattern = "\\.Rds$", full.names = TRUE, ignore.case = TRUE)

  if (length(rds_files) == 0) {
    warning(paste("No .Rds files found in", rds_dir))
    return(list(C1 = tibble(), C1_CATCH = tibble()))
  }

  file_data <- map(rds_files, function(f) {
    dl <- readRDS(f)

    message("File: ", basename(f), " -- top-level list names: ",
            paste(names(dl), collapse = ", "))

    if (!all(c("C1", "C1_CATCH") %in% names(dl))) {
      warning("Skipping ", basename(f),
              " -- does not contain both 'C1' and 'C1_CATCH' ",
              "(found: ", paste(names(dl), collapse = ", "), ")")
      return(list(C1 = tibble(), C1_CATCH = tibble()))
    }

    list(
      C1       = dl[["C1"]]       %>% mutate(subarea = subarea, source_file = basename(f)),
      C1_CATCH = dl[["C1_CATCH"]] %>% mutate(subarea = subarea, source_file = basename(f))
    )
  })

  list(
    C1       = map_dfr(file_data, "C1"),
    C1_CATCH = map_dfr(file_data, "C1_CATCH")
  )
}

all_data <- map(subareas, load_subarea)

C1_all       <- map_dfr(all_data, "C1")
C1_CATCH_all <- map_dfr(all_data, "C1_CATCH")

# -----------------------------
# CHECK WHAT COLUMNS ARE ACTUALLY AVAILABLE
# -----------------------------
message("Columns available in C1: ",       paste(names(C1_all), collapse = ", "))
message("Columns available in C1_CATCH: ", paste(names(C1_CATCH_all), collapse = ", "))

# -----------------------------
# BUILD HAUL-LEVEL TABLE
# NOTE: no spatial buffer / site filtering -- every haul we have data
# for is retained, across both subareas.
# -----------------------------

# Haul-level fields we want, IF they exist in this pull of the data.
# (Matches the actual C1 columns confirmed via the message() output:
# 756_USA_2026-01-02.Rds. NOTE: 759_USA_2026-01-02.Rds turned out to be
# an EA8 (CEMP land-based predator/seabird monitoring) file, not a C1
# catch/effort file -- it has no haul/catch data and is skipped above.
# If you locate a real C1-format file for 759, this list will pick up
# its columns automatically as long as the names match.)
haul_cols_wanted <- c(
  "c1_id", "subarea", "season_ccamlr", "ship_code", "vessel_name",
  "vessel_nationality_code", "target_species", "asd_code",
  "fishing_purpose_code", "trawl_technique", "catchperiod_code",
  "date_catchperiod_start",
  "datetime_set_start", "datetime_set_end",
  "datetime_haul_start", "datetime_haul_end", "datetime_timezone",
  "latitude_set_end", "longitude_set_end",
  "latitude_haul_start", "longitude_haul_start",
  "depth_gear_set_end_m", "depth_gear_haul_start_m",
  "depth_bottom_set_end_m", "depth_bottom_haul_start_m",
  "gear_type_code", "gear_type", "mesh_size_mm", "trawl_speed_kn",
  "comments", "source_file"
)
haul_cols_present <- intersect(haul_cols_wanted, names(C1_all))

haul_info <- C1_all %>%
  select(all_of(haul_cols_present))

# Catch is kept per haul x taxon (all species/taxa retained, not just krill)
catch_info <- C1_CATCH_all %>%
  select(any_of(c("c1_id", "subarea", "taxon_code", "greenweight_caught_kg")))

haul_catch_table <- haul_info %>%
  left_join(catch_info, by = c("c1_id", "subarea"))

# -----------------------------
# EXPORT TO EXCEL
# -----------------------------
write_xlsx(haul_catch_table, out_file)

message("Saved haul catch table (", nrow(haul_catch_table), " rows) to: ", out_file)
