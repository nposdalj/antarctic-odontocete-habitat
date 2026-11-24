# =============================================================================
# Script  : bremen_02_compare_maps.R
# Section : Supplement - sea-ice maps around the sites (AMSR2 + MODIS)
# Purpose : Map AMSR2 ASI sea-ice concentration (University of Bremen,
#           6.25 km) around a site on selected detection days, next to the
#           MODIS true-colour image of the same day, with the +/-20 km box,
#           the recorder and the 15% ice edge, to show where the ice was
#           relative to the site (e.g., open water inside the box on days of
#           high box-mean concentration).
# Inputs  : DOWNLOAD_DIR/bremen_asi/asi_s6250_<year>_region.tif
#           DOWNLOAD_DIR/modis_truecolor/<site>_<date>_<sat>.jpg (modis_01_fetch.R)
#           outputs/03_sea_ice/amsr2_daily_site.csv (amsr2_site_series.R)
#           F_PRESENCE_DAILY
# Outputs : outputs/03_sea_ice/ice_maps_modis.png/.pdf
# Author  : Natalie Posdaljian
# Date    : 2026-11-20
# Project : Odontocete acoustic occurrence and habitat associations at three
#           HARP sites in the South Shetland Islands (EI, KGI, CI)
# Notes   : Dates were chosen from the MODIS contact sheets (cloud-free over
#           the box). KGI 2015-10-05 (killer whale under heavy ice) was fully
#           clouded in both MODIS images and is not mapped. The comparison
#           with the Copernicus reanalysis ice was removed (decision
#           N. Posdaljian 2026-11-20: Copernicus ice is not used).
# =============================================================================

BASE_DIR <- "."
source(file.path(BASE_DIR, "final_scripts_odontocetes", "config.R"))
suppressPackageStartupMessages({ library(tidyverse); library(terra); library(patchwork) })

SECTION <- "03_sea_ice"
EDGE_PCT <- 15
MAP_HALF_M <- 60000
MAP_DATES <- tribble(~site, ~date, ~sat,
                     "CI", "2016-05-11", "Terra",
                     "CI", "2016-08-17", "Terra",
                     "CI", "2016-08-21", "Aqua",
                     "CI", "2016-09-02", "Terra")

asi <- rast(sort(list.files(file.path(DOWNLOAD_DIR, "bremen_asi"), pattern = "_region[.]tif$",
                            full.names = TRUE)))
asi[asi > 100] <- NA
site_series <- read_csv(out_path(SECTION, "amsr2_daily_site.csv"), show_col_types = FALSE) |>
  mutate(date = as.Date(date))
presence <- read_csv(F_PRESENCE_DAILY, show_col_types = FALSE) |> mutate(date = as.Date(date))
boxes <- site_boxes(BOX_HALF_KM)

centre <- function(st) {
  s <- SITES[SITES$site == st, ]
  crds(project(vect(data.frame(x = s$lon, y = s$lat), geom = c("x", "y"), crs = "EPSG:4326"),
               "EPSG:3031"))
}
box_3031 <- function(st) {   # box outline, densified so its edges follow the projection
  b <- boxes[boxes$site == st, ]
  lon <- seq(b$lon_min, b$lon_max, length.out = 20); lat <- seq(b$lat_min, b$lat_max, length.out = 20)
  ring <- data.frame(x = c(lon, rep(b$lon_max, 20), rev(lon), rep(b$lon_min, 20)),
                     y = c(rep(b$lat_min, 20), lat, rep(b$lat_max, 20), rev(lat)))
  as.data.frame(crds(project(vect(ring, geom = c("x", "y"), crs = "EPSG:4326"), "EPSG:3031")))
}

panel_pair <- function(site, date, sat) {
  st <- site; d <- as.Date(date); xy <- centre(st)
  win <- ext(xy[1] - MAP_HALF_M, xy[1] + MAP_HALF_M, xy[2] - MAP_HALF_M, xy[2] + MAP_HALF_M)
  r <- crop(project(asi[[format(d, "%Y%m%d")]], "EPSG:3031", method = "near"), win)
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE) |> setNames(c("x", "y", "pct"))
  bx <- box_3031(st)
  day <- site_series[site_series$site == st & site_series$date == d, ]
  sp <- presence |> filter(site == st, date == d) |>
    pivot_longer(any_of(c(MODELED_SPECIES, "BW58"))) |> filter(value > 0) |> pull(name)
  edge_txt <- if (is.na(day$dist_ice_edge_km)) "" else if (day$dist_ice_edge_km >= 0)
    sprintf("Recorder %.0f km outside the ice edge", day$dist_ice_edge_km) else
    sprintf("Recorder %.0f km inside the ice edge", -day$dist_ice_edge_km)
  common <- list(geom_path(data = bx, aes(x, y), inherit.aes = FALSE, colour = "red", linewidth = 0.6),
                 annotate("point", x = xy[1], y = xy[2], shape = 21, fill = "yellow", size = 2),
                 coord_equal(xlim = c(win$xmin, win$xmax), ylim = c(win$ymin, win$ymax), expand = FALSE),
                 theme_void(base_size = 8),
                 theme(plot.title = element_text(size = 7.5, face = "bold"),
                       plot.subtitle = element_text(size = 6.8)))
  p_asi <- ggplot(df, aes(x, y, fill = pct)) + geom_raster() +
    geom_contour(aes(z = pct), breaks = EDGE_PCT, colour = "white", linewidth = 0.4, na.rm = TRUE) +
    scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 100),
                         na.value = "grey60", name = "AMSR2 sea-ice concentration (%)") +
    labs(title = sprintf("%s, %s (%s detected)", st, format(d, "%d %b %Y"),
                         paste(sp, collapse = ", ")),
         subtitle = sprintf("Box mean %.0f%%; open water %.0f%% of box\n%s",
                            100 * day$ice_conc, 100 * day$open_water_frac, edge_txt)) +
    common
  img <- jpeg::readJPEG(file.path(DOWNLOAD_DIR, "modis_truecolor",
                                  sprintf("%s_%s_%s.jpg", st, d, sat)))
  p_mod <- ggplot() +
    annotation_raster(img, win$xmin, win$xmax, win$ymin, win$ymax) +
    labs(title = " ", subtitle = " \n ") + common   # satellite (Terra/Aqua) given in the legend
  p_asi | p_mod
}
maps <- pmap(MAP_DATES, panel_pair)
fig <- wrap_plots(maps, ncol = 1) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
for (ext in c("png", "pdf")) ggsave(out_path(SECTION, paste0("ice_maps_modis.", ext)), fig,
                                    width = 6.5, height = 2.9 * nrow(MAP_DATES) + 0.6, dpi = 300)
