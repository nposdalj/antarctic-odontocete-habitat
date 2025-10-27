# =======================================================
# Antarctic 3D & T–S Visualization Pipeline
# =======================================================
# Author: Natalie Posdaljian
# Description: Generates 3D and T–S plots across species and depths
# =======================================================

# ---- Packages ----
library(tidyverse)
library(plotly)
library(htmlwidgets)
library(scatterplot3d)
library(RColorBrewer)
library(MASS)            # LDA
library(viridis)         # color scale
library(oce)             # sigma-theta
library(reshape2)
library(tidyr)
library(stringr)
library(akima)      # for 2D interpolation (time x depth)
library(patchwork)  # to stack the 3 panels

# =======================================================
# ---- USER TOGGLES ----
# =======================================================
PLOT_INTERACTIVE_3D <- FALSE     # HTML interactive 3D plots
PLOT_STATIC_3D      <- FALSE     # Static 3D PNG (O2-SST-SSS)
PLOT_SEPARATION_3D  <- FALSE     # LDA/PCA 3D with convex hulls
PLOT_TS_BASIC       <- FALSE     # T–S plot (O2 as color)
PLOT_TS_SIGMA       <- FALSE     # T–S plot with σθ contours
PLOT_DEPTH_PROFILES_PER_SITE <- FALSE
PLOT_DEPTH_PROFILES_OVERLAY  <- FALSE   # overlay all sites per variable (optional)
MAKE_TIME_DEPTH_SECTIONS <- TRUE

# =======================================================
# ---- SETTINGS ----
# =======================================================
species_vec <- c("BW29","BW37","BW58","Oo","Pm","Gm")
sites_keep  <- c("EI","KGI","CI")
output_dir  <- "L:/Shared drives/Antarctic Marine Mammals/Data Visualizations/3D_T_S_0"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Site colors
site_levels <- c("EI","KGI","CI")
site_cols   <- setNames(brewer.pal(3, "Set2"), site_levels)

# =======================================================
# ---- HELPER FUNCTIONS ----
# =======================================================

# Site/species pretty names
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

depth_label <- function(dkey) paste0("z", dkey, " m")
depth_suffix <- function(dkey) paste0("_z", dkey, "m")

# =======================================================
# ---- LOAD DATA ----
# =======================================================
allData <- read.csv("/Users/nposd/Documents/GitHub/antarctic-odontocete-habitat/data/allData.csv")
if ("X" %in% colnames(allData)) {
  allData <- allData[, !(colnames(allData) == "X")]
}
if ("date" %in% names(allData)) {
  allData$date <- as.Date(allData$date)
}

# =======================================================
# ---- AUTO-DETECT DEPTHS ----
# =======================================================
temp_keys <- sub("^temperature_", "", grep("^temperature_", names(allData), value = TRUE))
salt_keys <- sub("^salinity_",    "", grep("^salinity_",    names(allData), value = TRUE))
o2_keys   <- sub("^o2_",          "", grep("^o2_",          names(allData), value = TRUE))

depth_keys <- Reduce(intersect, list(temp_keys, salt_keys, o2_keys))
depth_keys <- depth_keys[nchar(depth_keys) > 0]
if (length(depth_keys) == 0) stop("No matching temperature/salinity/o2 depth columns found.")

# =======================================================
# ---- PLOTTING FUNCTIONS ----
# =======================================================

# (A) Interactive 3D plot
save_plotly_3d <- function(df, sp, col_SST, col_SSS, col_O2, dkey) {
  df2 <- df %>% mutate(symbol = ifelse(Presence == 1, "circle", "triangle-up-open"))
  
  p <- plot_ly(
    data = df2,
    x = ~.data[[col_SST]],
    y = ~.data[[col_SSS]],
    z = ~.data[[col_O2]],
    color = ~Site,
    symbol = ~symbol,
    symbols = c("circle","triangle-up-open"),
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3, opacity = 0.7)
  ) |>
    layout(
      scene = list(
        xaxis = list(title = paste0("Temperature (°C) — ", depth_label(dkey))),
        yaxis = list(title = paste0("Salinity — ", depth_label(dkey))),
        zaxis = list(title = paste0("Oxygen (mL/L) — ", depth_label(dkey)))
      ),
      legend = list(title = list(text = "Site"))
    ) |>
    add_annotations(
      text = paste0("Species: ", name(sp),
                    " — Shapes: ● presence, △ absence — sea ice ≤ 15% — ", depth_label(dkey)),
      xref = "paper", yref = "paper", x = 0, y = 1.1, showarrow = FALSE
    )
  
  outfile <- file.path(output_dir, paste0("3D_O2_SST_SSS_", sp, depth_suffix(dkey), ".html"))
  saveWidget(p, outfile, selfcontained = TRUE)
  message("Saved HTML: ", outfile)
}

# (B) Static 3D PNG
save_s3d_png <- function(df, sp, col_SST, col_SSS, col_O2, dkey) {
  pchs <- ifelse(df$Presence == 1, 16, 2)
  cols <- site_cols[as.character(df$Site)]
  
  outfile <- file.path(output_dir, paste0("3D_O2_SST_SSS_", sp, depth_suffix(dkey), ".png"))
  png(outfile, width = 2000, height = 1800, res = 240)
  par(mar = c(5.5, 5.5, 5.5, 3) + 0.1)
  
  scatterplot3d(
    x = df[[col_SST]],
    y = df[[col_SSS]],
    z = df[[col_O2]],
    color = cols,
    pch   = pchs,
    angle = 55,
    grid  = TRUE,
    box   = TRUE,
    xlab  = paste0("Temperature (°C) — ", depth_label(dkey)),
    ylab  = paste0("Salinity — ", depth_label(dkey)),
    zlab  = paste0("Oxygen (mL/L) — ", depth_label(dkey)),
    main  = paste0(name(sp), " — EI/KGI/CI  |  ● presence, △ absence  |  sea ice ≤ 15%  |  ", depth_label(dkey))
  )
  
  legend("topright", inset = 0.02, title = "Site",
         legend = names(site_cols), col = site_cols, pch = 16, pt.cex = 1.2, bty = "n")
  legend("right", inset = 0.02,
         legend = c("Presence", "Absence"),
         pch = c(16, 2), pt.cex = 1.2, bty = "n")
  dev.off()
  message("Saved PNG (raw 3D): ", outfile)
}

# (C) Separation 3D (LDA/PCA fallback)
save_sep3d_png <- function(df, sp, col_SST, col_SSS, col_O2, dkey) {
  X <- scale(df[, c(col_O2, col_SST, col_SSS)])
  LD <- NULL
  if (dplyr::n_distinct(df$Site) >= 2) {
    lda_fit <- try(MASS::lda(df$Site ~ X[,1] + X[,2] + X[,3]), silent = TRUE)
    if (!inherits(lda_fit, "try-error")) {
      LD <- try(predict(lda_fit)$x, silent = TRUE)
      if (inherits(LD, "try-error")) LD <- NULL
    }
  }
  
  pca <- prcomp(X, center = FALSE, scale. = FALSE)
  PC  <- as.data.frame(pca$x)
  
  k <- ifelse(is.null(LD), 0, ncol(as.data.frame(LD)))
  if (k >= 3) {
    A <- as.data.frame(LD)[, 1:3, drop = FALSE]; axes_titles <- c("LD1","LD2","LD3")
  } else if (k == 2) {
    A <- cbind(as.data.frame(LD)[, 1:2, drop = FALSE], PC[,1, drop = FALSE])
    axes_titles <- c("LD1","LD2","PC1")
  } else if (k == 1) {
    A <- cbind(as.data.frame(LD)[, 1, drop = FALSE], PC[,1:2, drop = FALSE])
    axes_titles <- c("LD1","PC1","PC2")
  } else {
    A <- PC[, 1:3, drop = FALSE]; axes_titles <- c("PC1","PC2","PC3")
  }
  
  names(A) <- c("A1","A2","A3")
  df_plot <- dplyr::bind_cols(df, A)
  
  pchs <- ifelse(df_plot$Presence == 1, 16, 2)
  cols <- site_cols[as.character(df_plot$Site)]
  
  outfile <- file.path(output_dir, paste0("SEPARATION3D_", sp, depth_suffix(dkey), ".png"))
  png(outfile, width = 2200, height = 2000, res = 260)
  par(mar = c(6,6,6,4) + 0.1)
  
  s3d <- scatterplot3d(
    x = df_plot$A1, y = df_plot$A2, z = df_plot$A3,
    color = cols, pch = pchs, angle = 50,
    grid = TRUE, box = TRUE,
    xlab = axes_titles[1], ylab = axes_titles[2], zlab = axes_titles[3],
    main = sprintf("%s — EI/KGI/CI | ● presence, △ absence | sea ice ≤ 15%% | %s\nAxes: %s",
                   name(sp), depth_label(dkey), paste(axes_titles, collapse = ", "))
  )
  
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
  
  legend("topright", inset = 0.02, title = "Site",
         legend = names(site_cols), col = site_cols,
         pch = 16, pt.cex = 1.2, bty = "n")
  legend("right", inset = 0.02,
         legend = c("Presence", "Absence"),
         pch = c(16, 2), pt.cex = 1.2, bty = "n")
  dev.off()
  message("Saved PNG (separation 3D): ", outfile)
}

# (D) T–S plots (with optional sigma-theta contours)
save_TS_png <- function(df, sp, col_SST, col_SSS, col_O2, dkey, with_sigma = FALSE) {
  
  base <- ggplot(df, aes(x = .data[[col_SSS]], y = .data[[col_SST]])) +
    geom_point(aes(color = .data[[col_O2]],
                   shape = factor(Presence, levels = c(1,0),
                                  labels = c("Presence","Absence"))),
               alpha = 0.8, size = 1.5, stroke = 0.7) +
    scale_color_viridis(name = "Oxygen (mL/L)", option = "C", direction = -1) +
    scale_shape_manual(name = NULL, values = c(16, 2)) +
    facet_wrap(~ Site, nrow = 1) +
    labs(x = paste0("Salinity — ", depth_label(dkey)),
         y = paste0("Temperature (°C) — ", depth_label(dkey)),
         title = paste0(name(sp), " — T–S by site (O₂ as color) — ", depth_label(dkey)),
         subtitle = "● presence, △ absence | sea ice ≤ 15%") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right", plot.title = element_text(face = "bold"))
  
  if (!with_sigma) {
    fn <- file.path(output_dir, paste0("TS_O2_facets_", sp, depth_suffix(dkey), ".png"))
    ggsave(fn, base, width = 11, height = 4.5, dpi = 300)
    message("Saved T–S (O2 color): ", fn)
  } else {
    s_range <- range(df[[col_SSS]], na.rm = TRUE)
    t_range <- range(df[[col_SST]], na.rm = TRUE)
    s_seq <- seq(floor(s_range[1]*10)/10, ceiling(s_range[2]*10)/10, by = 0.05)
    t_seq <- seq(floor(t_range[1]),      ceiling(t_range[2]),      by = 0.25)
    grid <- expand.grid(S = s_seq, T = t_seq)
    grid$sigma_theta <- swSigmaTheta(salinity = grid$S, temperature = grid$T, pressure = 0)
    
    g <- base +
      geom_contour(data = grid,
                   aes(x = S, y = T, z = sigma_theta),
                   color = "grey35", size = 0.3, alpha = 0.85) +
      labs(subtitle = "● presence, △ absence | sea ice ≤ 15% | σθ contours (0 dbar)")
    
    fn <- file.path(output_dir, paste0("TS_O2_sigma_facets_", sp, depth_suffix(dkey), ".png"))
    ggsave(fn, g, width = 11.5, height = 4.8, dpi = 300)
    message("Saved T–S with σθ contours: ", fn)
  }
}

# =======================================================
# ---- MAIN LOOP (species × depth) ----
# =======================================================
for (sp in species_vec) {
  for (dkey in depth_keys) {
    
    col_O2  <- paste0("o2_",          dkey)
    col_SST <- paste0("temperature_", dkey)
    col_SSS <- paste0("salinity_",    dkey)
    
    if (!all(c(col_O2, col_SST, col_SSS) %in% names(allData))) next
    
    df <- allData %>%
      filter(Site %in% sites_keep) %>%
      mutate(
        Presence = ifelse(is.na(.data[[sp]]), NA_real_, as.numeric(.data[[sp]] > 0)),
        Site = factor(Site, levels = site_levels)
      ) %>%
      filter(
        !is.na(.data[[col_O2]]),
        !is.na(.data[[col_SST]]),
        !is.na(.data[[col_SSS]]),
        !is.na(Presence),
        !is.na(ice_conc),
        ice_conc <= 15
      )
    
    if (nrow(df) == 0) {
      message("No data after filtering for ", sp, " at ", depth_label(dkey), " — skipping.")
      next
    }
    
    if (PLOT_INTERACTIVE_3D) save_plotly_3d(df, sp, col_SST, col_SSS, col_O2, dkey)
    if (PLOT_STATIC_3D)      save_s3d_png(df, sp, col_SST, col_SSS, col_O2, dkey)
    if (PLOT_SEPARATION_3D)  save_sep3d_png(df, sp, col_SST, col_SSS, col_O2, dkey)
    if (PLOT_TS_BASIC)       save_TS_png(df, sp, col_SST, col_SSS, col_O2, dkey, with_sigma = FALSE)
    if (PLOT_TS_SIGMA)       save_TS_png(df, sp, col_SST, col_SSS, col_O2, dkey, with_sigma = TRUE)
  }
}

if (PLOT_DEPTH_PROFILES_PER_SITE || PLOT_DEPTH_PROFILES_OVERLAY) {
  profile_dir <- file.path(output_dir, "Depth_Profiles")
  dir.create(profile_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Helper: make long tidy table of temperature/salinity/o2 with numeric depth suffixes
  make_long_env <- function(df) {
    # find only columns like temperature_0, salinity_11, o2_635 (numeric suffix only)
    keep_cols <- grep("^(temperature|salinity|o2)_\\d+$", names(df), value = TRUE)
    if (length(keep_cols) == 0) return(data.frame(var = character(), depth = numeric(), value = numeric()))
    
    # subset safely using base R, then pivot_longer
    df_sub <- df[, keep_cols, drop = FALSE]
    long <- df_sub %>%
      pivot_longer(cols = everything(), names_to = "var_depth", values_to = "value") %>%
      separate(var_depth, into = c("var","depth"), sep = "_", remove = TRUE) %>%
      mutate(depth = suppressWarnings(as.numeric(depth))) %>%
      filter(!is.na(depth), var %in% c("temperature","salinity","o2"))
    return(long)
  }
  
  # colors and labels
  var_cols  <- c(temperature = "#E64B35", salinity = "#4DBBD5", o2 = "#00A087")
  var_labs  <- c(temperature = "Temperature (°C)", salinity = "Salinity", o2 = "Oxygen (mL/L)")
  
  # ---------- (A) Per-site profiles ----------
  if (PLOT_DEPTH_PROFILES_PER_SITE) {
    for (site in sites_keep) {
      df_site <- allData[allData$Site == site, ]
      long    <- make_long_env(df_site)
      if (nrow(long) == 0) {
        message("No depth columns for site ", site, " — skipping.")
        next
      }
      
      means <- aggregate(value ~ var + depth, data = long, FUN = mean, na.rm = TRUE)
      means <- means[order(means$var, means$depth), ]
      
      interp_list <- lapply(split(means, means$var), function(d) {
        d <- d[!is.na(d$value) & !is.na(d$depth), ]
        if (nrow(d) < 2) return(NULL)
        depth_grid <- seq(min(d$depth), max(d$depth), length.out = 300)
        data.frame(
          var   = unique(d$var),
          depth = depth_grid,
          value = approx(d$depth, d$value, xout = depth_grid, rule = 2)$y
        )
      })
      interp <- do.call(rbind, interp_list)
      if (is.null(interp) || nrow(interp) == 0) {
        message("Not enough depth points for site ", site, " — skipping.")
        next
      }
      
      p <- ggplot(interp, aes(x = value, y = depth, color = var)) +
        geom_path(size = 1.1) +
        scale_y_reverse(expand = c(0.02, 0)) +
        scale_color_manual(values = var_cols, labels = var_labs, name = NULL) +
        facet_wrap(~ var, ncol = 3, scales = "free_x",
                   labeller = as_labeller(var_labs)) +
        labs(
          title = paste("Depth Profiles —", name(site)),
          x = NULL, y = "Depth (m)"
        ) +
        theme_bw(base_size = 12) +
        theme(legend.position = "none",
              strip.text = element_text(face = "bold"),
              plot.title = element_text(face = "bold"))
      
      fpath <- file.path(profile_dir, paste0("DepthProfile_", site, ".png"))
      ggsave(fpath, p, width = 10, height = 4, dpi = 300)
      message("Saved depth profile: ", fpath)
    }
  }
  
  # ---------- (B) Overlay all sites per variable (optional) ----------
  if (PLOT_DEPTH_PROFILES_OVERLAY) {
    all_rows <- list()
    for (site in sites_keep) {
      df_site <- allData[allData$Site == site, ]
      long <- make_long_env(df_site)
      if (nrow(long) > 0) long$Site <- site
      all_rows[[site]] <- long
    }
    long_all <- do.call(rbind, all_rows)
    if (nrow(long_all) == 0) {
      message("No depth columns found for overlay — skipping.")
    } else {
      means_all <- aggregate(value ~ var + Site + depth, data = long_all, FUN = mean, na.rm = TRUE)
      means_all <- means_all[order(means_all$var, means_all$Site, means_all$depth), ]
      
      interp_all <- do.call(rbind, by(means_all, means_all[, c("var","Site")], function(d) {
        d <- d[!is.na(d$value) & !is.na(d$depth), ]
        if (nrow(d) < 2) return(NULL)
        depth_grid <- seq(min(d$depth), max(d$depth), length.out = 300)
        data.frame(
          var   = unique(d$var),
          Site  = unique(d$Site),
          depth = depth_grid,
          value = approx(d$depth, d$value, xout = depth_grid, rule = 2)$y
        )
      }))
      
      for (v in c("temperature","salinity","o2")) {
        sub <- interp_all[interp_all$var == v, ]
        if (nrow(sub) == 0) next
        
        p <- ggplot(sub, aes(x = value, y = depth, color = Site)) +
          geom_path(size = 1.1) +
          scale_y_reverse(expand = c(0.02, 0)) +
          scale_color_manual(values = site_cols, name = "Site") +
          labs(
            title = paste0("Depth Profile — ", var_labs[[v]]),
            x = NULL, y = "Depth (m)"
          ) +
          theme_bw(base_size = 12) +
          theme(legend.position = "right",
                plot.title = element_text(face = "bold"))
        
        fpath <- file.path(profile_dir, paste0("DepthProfile_ALLSITES_", v, ".png"))
        ggsave(fpath, p, width = 6.5, height = 5.5, dpi = 300)
        message("Saved overlay depth profile: ", fpath)
      }
    }
  }
}

if (MAKE_TIME_DEPTH_SECTIONS) {
  
  hov_dir <- file.path(output_dir, "TimeDepth_Sections")
  dir.create(hov_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Helper: build long tidy (date, Site, var, depth, value) from *_<depth> columns
  make_long_env_all <- function(df) {
    keep_cols <- grep("^(temperature|salinity|o2)_\\d+$", names(df), value = TRUE)
    if (length(keep_cols) == 0) {
      return(data.frame(date = as.Date(character()),
                        Site = character(), var = character(),
                        depth = numeric(), value = numeric()))
    }
    df_sub <- df[, c("date", "Site", keep_cols), drop = FALSE]
    long <- df_sub |>
      pivot_longer(cols = all_of(keep_cols), names_to = "var_depth", values_to = "value") |>
      separate(var_depth, into = c("var","depth"), sep = "_", remove = TRUE) |>
      mutate(depth = suppressWarnings(as.numeric(depth))) |>
      filter(!is.na(depth), var %in% c("temperature","salinity","o2"))
    long
  }
  
  var_labels <- c(temperature = "Temperature (°C)",
                  salinity    = "Salinity",
                  o2          = "Oxygen (mL/L)")
  
  # ---- Interpolate one variable at one site over a regular (time, depth) grid ----
  interp_time_depth <- function(df_var, extrapolate = FALSE) {
    df2 <- df_var[!is.na(df_var$value) & !is.na(df_var$depth) & !is.na(df_var$date), ]
    if (nrow(df2) < 5) return(NULL)
    
    tx <- as.numeric(df2$date)
    z  <- df2$value
    y  <- df2$depth
    
    xo <- seq(min(tx), max(tx), by = 1)                 # daily grid
    yo <- seq(min(y),  max(y),  length.out = 300)
    
    ii <- akima::interp(x = tx, y = y, z = z,
                        xo = xo, yo = yo,
                        duplicate = "mean",
                        linear = TRUE,
                        extrap = extrapolate)            # FALSE = no extrapolation
    
    out <- expand.grid(date_num = ii$x, depth = ii$y)
    out$value <- as.vector(ii$z)
    out$date  <- as.Date(out$date_num, origin = "1970-01-01")
    out
  }
  
  # ---- Pretty heatmap with explicit per-site limits ----
  hov_plot <- function(grid_df, title, fill_lab, lims) {
    ggplot(grid_df, aes(x = date, y = depth, fill = value)) +
      geom_raster(interpolate = TRUE, na.rm = TRUE) +
      scale_y_reverse(expand = c(0.001, 0)) +
      scale_fill_viridis(option = "C", name = fill_lab,
                         limits = lims, oob = scales::squish) +
      labs(x = NULL, y = "Depth (m)", title = title) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0),
        axis.title.x = element_blank(),
        panel.grid = element_blank()
      )
  }
  
  # ---- Build once: long tidy data frame ----
  long_all <- make_long_env_all(allData[allData$Site %in% sites_keep &
                                          !is.na(allData$ice_conc) &
                                          allData$ice_conc <= 15, ])
  
  # ---- One figure per site (stacked T/S/O) ----
  for (site in sites_keep) {
    site_long <- subset(long_all, Site == site)
    if (nrow(site_long) == 0) next
    
    plots <- list()
    for (v in c("temperature", "salinity", "o2")) {
      sub <- subset(site_long, var == v)
      if (nrow(sub) < 5) next
      
      # Per-site, per-variable color limits from observed data
      v_minmax <- range(sub$value, na.rm = TRUE)
      
      # Interpolate over time and depth (set extrapolate = FALSE to avoid artifacts)
      grid <- interp_time_depth(sub, extrapolate = FALSE)
      if (is.null(grid)) next
      
      plots[[v]] <- hov_plot(
        grid,
        paste0(name(site), " — ", var_labels[[v]]),
        var_labels[[v]],
        lims = v_minmax
      )
    }
    
    if (length(plots) == 0) {
      message("No data to plot for site ", site)
      next
    }
    
    # Stack vertically (shared x-axis visually)
    fig <- wrap_plots(plots, ncol = 1, guides = "collect") &
      theme(legend.position = "right")
    
    # Save as wide PNG
    fn <- file.path(hov_dir, paste0("TimeDepth_TSO_", site, ".png"))
    ggsave(fn, fig, width = 15, height = 5, dpi = 300)
    message("Saved time–depth section: ", fn)
  }
}


message("✅ All selected plots complete.")
