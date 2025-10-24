# ---- Packages ----
library(tidyverse)
library(scatterplot3d)
library(RColorBrewer)

# ---- Settings ----
species_vec <- c("BW29","BW37","Oo","Pm","Gm")
sites_keep  <- c("EI","KGI","CI")
out_dir_png <- "L:/Shared drives/Antarctic Marine Mammals/Data Visualizations/3D_T_S_0"

dir.create(out_dir_png, showWarnings = FALSE, recursive = TRUE)

# ---- Helper: readable names ----
name <- function(abbrev) {
  switch(abbrev,
         "CI"="Clarence Island",
         "KGI"="King George Island",
         "EI"="Elephant Island",
         "BW29"="Southern Bottlenose Whale",
         "BW37"="Gray's and Strap-toothed Whale BW37",
         "BW58"="Gray's and Strap-toothed Whale BW58",
         "Gm"="Long-finned Pilot Whale",
         "Oo"="Killer Whale",
         "Pm"="Sperm Whale",
         abbrev)
}

# ---- Load data ----
allData <- read.csv("/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/data/allData.csv") |>
  select(-X) |>
  mutate(date = as.Date(date))

# Columns for environmental variables (edit if yours differ)
col_O2  <- "o2_0"
col_SST <- "temperature_0"
col_SSS <- "salinity_0"

# Site colors
site_levels <- c("EI","KGI","CI")
site_cols   <- setNames(brewer.pal(3, "Set2"), site_levels)

# ---- Loop: make one PNG per species ----
for (sp in species_vec) {
  
  df <- allData |>
    filter(Site %in% sites_keep) |>
    mutate(
      Presence = ifelse(is.na(.data[[sp]]), NA_real_, as.numeric(.data[[sp]] > 0)),
      Site = factor(Site, levels = site_levels)
    ) |>
    filter(
      !is.na(.data[[col_O2]]),
      !is.na(.data[[col_SST]]),
      !is.na(.data[[col_SSS]]),
      !is.na(Presence),
      !is.na(ice_conc),
      ice_conc <= 15
    )
  
  if (nrow(df) == 0) next
  
  # Aesthetics
  pchs <- ifelse(df$Presence == 1, 16, 2)   # 16=filled circle, 2=open triangle
  cols <- site_cols[as.character(df$Site)]
  
  # Save PNG
  png(
    filename = file.path(out_dir_png, paste0("3D_O2_SST_SSS_v2", sp, ".png")),
    width = 2000, height = 1800, res = 240
  )
  par(mar = c(5.5, 5.5, 5.5, 3) + 0.1)
  
  s3d <- scatterplot3d(
    x = df[[col_SST]],     # SST (°C)
    y = df[[col_SSS]],     # SSS
    z = df[[col_O2]],      # Oxygen (mL/L)
    color = cols,
    pch   = pchs,
    angle = 55,            # azimuth for perspective
    grid  = TRUE,
    box   = TRUE,
    xlab  = "Sea Surface Temperature (°C)",
    ylab  = "Sea Surface Salinity",
    zlab  = "Oxygen (mL/L)",
    main  = paste0(name(sp), " — EI/KGI/CI  |  ● presence, △ absence  |  sea ice ≤ 15%")
  )
  
  # Legends
  legend("topright", inset = 0.02,
         title = "Site",
         legend = names(site_cols),
         col = site_cols, pch = 16, pt.cex = 1.2, bty = "n")
  legend("right", inset = 0.02,
         legend = c("Presence", "Absence"),
         pch = c(16, 2), pt.cex = 1.2, bty = "n")
  
  dev.off()
  message("Saved PNG for ", sp)
}
