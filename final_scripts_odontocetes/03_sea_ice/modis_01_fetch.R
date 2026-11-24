# =============================================================================
# Script  : modis_01_fetch.R
# Section : 2.3 / Supplement - visual sea-ice imagery around the sites
# Purpose : Fetch MODIS true-colour images (Terra and Aqua corrected
#           reflectance, 250 m) from NASA GIBS (open WMS, no login) for days
#           with detections under substantial sea ice (box mean >= 30%), in a
#           120 x 120 km window centred on the site, so cloud-free dates can be
#           chosen for the figure.
# Inputs  : F_PRESENCE_DAILY, F_ENV_DAILY (box-mean ice)
#           https://gibs.earthdata.nasa.gov/wms/epsg3031/best/wms.cgi
# Outputs : DOWNLOAD_DIR/modis_truecolor/<site>_<yyyy-mm-dd>_<Terra|Aqua>.jpg
#           DOWNLOAD_DIR/modis_truecolor/candidate_dates.csv
#           outputs/03_sea_ice/modis_contact_sheet_<site>.png (for choosing dates)
# Author  : Natalie Posdaljian
# Date    : 2026-11-20
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Images are in Antarctic polar stereographic (EPSG:3031). We
#           acknowledge the use of imagery from NASA Worldview/GIBS, part of
#           NASA's Earth Science Data and Information System [CONFIRM wording].
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
suppressPackageStartupMessages({ library(tidyverse); library(terra) })

MIN_ICE <- 0.30
HALF_M <- 60000
PIXELS <- 600
GIBS <- "https://gibs.earthdata.nasa.gov/wms/epsg3031/best/wms.cgi"
LAYERS <- c(Terra = "MODIS_Terra_CorrectedReflectance_TrueColor",
            Aqua = "MODIS_Aqua_CorrectedReflectance_TrueColor")
out_dir <- file.path(DOWNLOAD_DIR, "modis_truecolor")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 120, HTTPUserAgent = "research-client")

presence <- read_csv(F_PRESENCE_DAILY, show_col_types = FALSE) |> mutate(date = as.Date(date))
ice <- read_csv(F_ENV_DAILY, show_col_types = FALSE) |> transmute(site, date = as.Date(date), ice_conc)
cand <- presence |>
  inner_join(ice, by = c("site", "date")) |>
  filter(ice_conc >= MIN_ICE) |>
  pivot_longer(any_of(c(MODELED_SPECIES, "BW58")), names_to = "species", values_to = "minutes") |>
  filter(minutes > 0) |>
  group_by(site, date, ice_pct = round(100 * ice_conc)) |>
  summarise(species = paste(species, collapse = ","), .groups = "drop")
write_csv(cand, file.path(out_dir, "candidate_dates.csv"))

centre <- function(st) {
  s <- SITES[SITES$site == st, ]
  crds(project(vect(data.frame(x = s$lon, y = s$lat), geom = c("x", "y"), crs = "EPSG:4326"),
               "EPSG:3031"))
}
for (i in seq_len(nrow(cand))) {
  xy <- centre(cand$site[i])
  bbox <- paste(xy[1] - HALF_M, xy[2] - HALF_M, xy[1] + HALF_M, xy[2] + HALF_M, sep = ",")
  for (sat in names(LAYERS)) {
    f <- file.path(out_dir, sprintf("%s_%s_%s.jpg", cand$site[i], cand$date[i], sat))
    if (file.exists(f)) next
    url <- sprintf(paste0("%s?SERVICE=WMS&REQUEST=GetMap&VERSION=1.3.0&LAYERS=%s&STYLES=",
                          "&CRS=EPSG:3031&BBOX=%s&WIDTH=%d&HEIGHT=%d&FORMAT=image/jpeg&TIME=%s"),
                   GIBS, LAYERS[[sat]], bbox, PIXELS, PIXELS, cand$date[i])
    tryCatch(download.file(url, f, mode = "wb", quiet = TRUE),
             error = function(e) message("failed: ", basename(f)))
  }
}

# Contact sheets (one per site) to choose cloud-free dates
for (st in unique(cand$site)) {
  files <- sort(list.files(out_dir, pattern = paste0("^", st, "_.*[.]jpg$"), full.names = TRUE))
  n <- length(files); nc <- 10; nr <- ceiling(n / nc)
  png(out_path("03_sea_ice", paste0("modis_contact_sheet_", st, ".png")),
      width = 220 * nc, height = 240 * nr)
  par(mfrow = c(nr, nc), mar = c(0, 0, 1.6, 0))
  for (f in files) {
    plot.new(); title(sub("[.]jpg$", "", basename(f)), cex.main = 1.1)
    rasterImage(jpeg::readJPEG(f), 0, 0, 1, 1)
    rect(0.5 - 20 / 120, 0.5 - 20 / 120, 0.5 + 20 / 120, 0.5 + 20 / 120, border = "red", lwd = 2)
  }
  dev.off()
}
message(nrow(cand), " candidate dates; images in ", out_dir)
