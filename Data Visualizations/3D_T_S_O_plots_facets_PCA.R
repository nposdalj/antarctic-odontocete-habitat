# ---- Packages ----
library(tidyverse)
library(dplyr)
library(MASS)            # LDA
library(scatterplot3d)   # static 3D PNG
library(RColorBrewer)
library(viridis)         # color scale for O2 in T–S

# ---- Settings ----
species_vec <- c("BW29","BW37","Oo","Pm","Gm")
sites_keep  <- c("EI","KGI","CI")

out_dir <- "L:/Shared drives/Antarctic Marine Mammals/Data Visualizations/3D_T_S_0"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Pretty names ----
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
allData <- read.csv("/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/data/allData.csv")

# Safely drop column X if it exists (common when saving with row names)
if ("X" %in% colnames(allData)) {
  allData <- allData[, !(colnames(allData) == "X")]
}

# Convert date column to Date type
allData$date <- as.Date(allData$date)

# Columns (edit if yours differ)
col_O2  <- "o2_0"
col_SST <- "temperature_0"
col_SSS <- "salinity_0"

site_levels <- c("EI","KGI","CI")
site_cols   <- setNames(brewer.pal(3, "Set2"), site_levels)

# ---- Helpers ----
save_3d_png <- function(df_plot, main_title, axes_titles, file_png,
                        outline_by_site = TRUE, angle = 50) {
  
  pchs <- ifelse(df_plot$Presence == 1, 16, 2)   # ● presence, △ absence
  cols <- site_cols[as.character(df_plot$Site)]
  
  png(file_png, width = 2200, height = 2000, res = 260)
  par(mar = c(6,6,6,4) + 0.1)
  
  s3d <- scatterplot3d(
    x = df_plot$A1, y = df_plot$A2, z = df_plot$A3,
    color = cols, pch = pchs, angle = angle,
    grid = TRUE, box = TRUE,
    xlab = axes_titles[1], ylab = axes_titles[2], zlab = axes_titles[3],
    main = main_title
  )
  
  if (outline_by_site) {
    for (st in site_levels) {
      sub <- df_plot %>% filter(Site == st)
      if (nrow(sub) < 5) next
      proj <- s3d$xyz.convert(sub$A1, sub$A2, sub$A3)
      hh <- try(chull(proj$x, proj$y), silent = TRUE)
      if (!inherits(hh, "try-error") && length(hh) > 2) {
        polygon(proj$x[hh], proj$y[hh],
                border = site_cols[st],
                col = adjustcolor(site_cols[st], alpha.f = 0.12),
                lwd = 2)
      }
    }
  }
  
  legend("topright", inset = 0.02, title = "Site",
         legend = names(site_cols), col = site_cols,
         pch = 16, pt.cex = 1.2, bty = "n")
  legend("right", inset = 0.02,
         legend = c("Presence", "Absence"),
         pch = c(16, 2), pt.cex = 1.2, bty = "n")
  
  dev.off()
}

# ---- Loop per species ----
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
  
  # ---------- LDA view (max between-site separation) ----------
  # --- Build 3 plotting axes robustly (LDA + fallback to PCA) ---
  # Standardize env vars
  X <- scale(df[, c(col_O2, col_SST, col_SSS)])
  
  # If fewer than 2 unique sites, skip LDA entirely
  if (n_distinct(df$Site) < 2) next
  
  # LDA for site separation
  lda_fit <- MASS::lda(df$Site ~ X[,1] + X[,2] + X[,3])
  LD <- predict(lda_fit)$x
  LD <- as.data.frame(LD)           # n x k, where k <= min(3, #sites-1)
  
  # PCA for filler axes if needed
  pca <- prcomp(X, center = FALSE, scale. = FALSE)
  PC  <- as.data.frame(pca$x)
  
  # Assemble exactly three axes for plotting
  k <- ifelse(is.null(LD) || ncol(LD) == 0, 0, ncol(LD))
  if (k >= 3) {
    A <- LD[, 1:3, drop = FALSE]; axes_titles <- c("LD1","LD2","LD3")
  } else if (k == 2) {
    A <- cbind(LD[, 1:2, drop = FALSE], PC[, 1, drop = FALSE])
    axes_titles <- c("LD1","LD2","PC1")
  } else if (k == 1) {
    A <- cbind(LD[, 1, drop = FALSE], PC[, 1:2, drop = FALSE])
    axes_titles <- c("LD1","PC1","PC2")
  } else {  # k == 0 (edge case: LDA failed)
    A <- PC[, 1:3, drop = FALSE]
    axes_titles <- c("PC1","PC2","PC3")
  }
  
  names(A) <- c("A1","A2","A3")
  df_plot <- dplyr::bind_cols(df, A)
  
  
  # ---- Save the separation view (LDA + PCA fallback) ----
  sep_title <- sprintf(
    "%s — EI/KGI/CI | ● presence, △ absence | sea ice ≤ 15%%\nAxes: %s",
    name(sp), paste(axes_titles, collapse = ", ")
  )
  sep_file <- file.path(out_dir, paste0("SEPARATION3D_", sp, ".png"))
  save_3d_png(df_plot, sep_title, axes_titles, sep_file, outline_by_site = TRUE, angle = 50)
  
  # ---- Make a pure PCA view (PC1–PC3) and save ----
  pca_df <- dplyr::bind_cols(
    df,
    setNames(as.data.frame(PC[, 1:3, drop = FALSE]), c("A1","A2","A3"))
  )
  pca_title <- sprintf(
    "%s — EI/KGI/CI | ● presence, △ absence | sea ice ≤ 15%%\nPCA of O2/SST/SSS",
    name(sp)
  )
  pca_file <- file.path(out_dir, paste0("PCA3D_", sp, ".png"))
  save_3d_png(pca_df, pca_title, c("PC1","PC2","PC3"), pca_file, outline_by_site = TRUE, angle = 50)
  
  
  # ---------- T–S diagram with O2 as color (faceted by site) ----------
  ts_file <- file.path(out_dir, paste0("TS_O2_facets_", sp, ".png"))
  g <- df %>%
    ggplot(aes(x = .data[[col_SSS]], y = .data[[col_SST]])) +
    geom_point(aes(color = .data[[col_O2]],
                   shape = factor(Presence, levels = c(1,0), labels = c("Presence","Absence"))),
               alpha = 0.8, size = 1.5, stroke = 0.7) +
    scale_color_viridis(name = "Oxygen (mL/L)", option = "C", direction = -1) +
    scale_shape_manual(name = NULL, values = c(16, 2)) +
    facet_wrap(~ Site, nrow = 1) +
    labs(x = "Sea Surface Salinity", y = "Sea Surface Temperature (°C)",
         title = paste0(name(sp), " — T–S by site (O₂ as color)"),
         subtitle = "● presence, △ absence | sea ice ≤ 15%") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right", plot.title = element_text(face = "bold"))
  
  ggsave(ts_file, g, width = 11, height = 4.5, dpi = 300)
  
}
