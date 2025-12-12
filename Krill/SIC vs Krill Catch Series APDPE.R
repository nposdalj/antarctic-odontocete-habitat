## ----------------------------------------------------------
## Setup
## ----------------------------------------------------------

# Load data
df <- read.csv(
  "L:/Shared drives/Antarctic Marine Mammals/Krill Data/Sofie's Analysis/APDPE_krill_env_monthly.csv",
  header = TRUE
)

df_year <- read.csv(
  "L:/Shared drives/Antarctic Marine Mammals/Krill Data/Sofie's Analysis/APDPE_krill_env_yearly.csv",
  header = TRUE
)

# ---- Output directory ----
out_dir <- "L:/Shared drives/Antarctic Marine Mammals/Krill Data/Sofie's Analysis/Figures"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ---- Krill catch units (EDIT IF NEEDED) ----
krill_units <- "t"   # e.g., "t", "kg", "t km^-2"


## ----------------------------------------------------------
## Generic plotting function (black points)
## ----------------------------------------------------------

plot_krill_lm_black <- function(df, xvar, yvar,
                                main, xlab, ylab,
                                pch = 19, cex = 1.5) {
  
  # Clean data
  df_clean <- df[!is.na(df[[xvar]]) & !is.na(df[[yvar]]), ]
  x <- df_clean[[xvar]]
  y <- df_clean[[yvar]]
  
  # Scatterplot
  plot(x, y,
       pch = pch,
       cex = cex,
       col = "black",
       xlab = xlab,
       ylab = ylab,
       main = main)
  
  # Linear model
  lm_model <- lm(y ~ x)
  abline(lm_model, col = "black", lwd = 3)
  
  # Stats
  summary_lm <- summary(lm_model)
  p_value <- coef(summary_lm)[2, 4]
  r2 <- summary_lm$r.squared
  
  legend("topright",
         bty = "n",
         legend = c(
           paste0("R² = ", round(r2, 3)),
           paste0("p = ", signif(p_value, 3))
         ))
}


###################
#### MONTHLY   ####
###################

# Sea ice % for monthly data
df$ice_pct <- df$sice_conc_mean_month_mean * 100

# Order by time
df <- df[order(df$Calendar_Year, df$Month), ]

# ---- Monthly year lags (12-month increments) ----
df$ice_lag0y <- df$ice_pct
df$ice_lag1y <- c(rep(NA, 12), head(df$ice_pct, -12))
df$ice_lag2y <- c(rep(NA, 24), head(df$ice_pct, -24))
df$ice_lag3y <- c(rep(NA, 36), head(df$ice_pct, -36))


## ---- MONTHLY: Same-year single plot ----
png(file.path(out_dir, "APDPE_monthly_same_year_SIC_vs_krill.png"),
    width = 2000, height = 1600, res = 300)

par(mfrow = c(1,1), mar = c(5,4,4,2) + 0.1)

plot_krill_lm_black(
  df,
  xvar = "ice_lag0y",
  yvar = "krill_weight",
  main = "APDPE Monthly: Same-Year SIC vs Krill Catch",
  xlab = "Sea Ice Concentration Same Year (%)",
  ylab = paste0("Krill Catch (", krill_units, ")")
)

dev.off()


## ---- MONTHLY: 2×2 lag panel ----
png(file.path(out_dir, "APDPE_monthly_SIC_vs_krill_lags_2x2.png"),
    width = 2400, height = 2400, res = 300)

par(mfrow = c(2,2), mar = c(5,4,4,2) + 0.1)

plot_krill_lm_black(df, "ice_lag0y", "krill_weight",
                    "Monthly: Same Year",
                    "SIC Same Year (%)",
                    paste0("Krill Catch (", krill_units, ")"))

plot_krill_lm_black(df, "ice_lag1y", "krill_weight",
                    "Monthly: 1-year Lag",
                    "SIC 1-year Lag (%)",
                    paste0("Krill Catch (", krill_units, ")"))

plot_krill_lm_black(df, "ice_lag2y", "krill_weight",
                    "Monthly: 2-year Lag",
                    "SIC 2-year Lag (%)",
                    paste0("Krill Catch (", krill_units, ")"))

plot_krill_lm_black(df, "ice_lag3y", "krill_weight",
                    "Monthly: 3-year Lag",
                    "SIC 3-year Lag (%)",
                    paste0("Krill Catch (", krill_units, ")"))

dev.off()


###################
#### YEARLY    ####
###################

# Sea ice % for yearly data
df_year$ice_pct <- df_year$mean_sice_conc * 100

# Order by year
df_year <- df_year[order(df_year$Calendar_Year), ]

# ---- Yearly lags ----
df_year$ice_lag1 <- c(NA, head(df_year$ice_pct, -1))
df_year$ice_lag2 <- c(NA, NA, head(df_year$ice_pct, -2))
df_year$ice_lag3 <- c(NA, NA, NA, head(df_year$ice_pct, -3))


## ---- YEARLY: Same-year single plot ----
png(file.path(out_dir, "APDPE_yearly_same_year_SIC_vs_krill.png"),
    width = 2000, height = 1600, res = 300)

par(mfrow = c(1,1), mar = c(5,4,4,2) + 0.1)

plot_krill_lm_black(
  df_year,
  xvar = "ice_pct",
  yvar = "total_krill",
  main = "APDPE Yearly: Same-Year SIC vs Krill Catch",
  xlab = "Sea Ice Concentration Same Year (%)",
  ylab = paste0("Krill Catch (", krill_units, ")")
)

dev.off()


## ---- YEARLY: 2×2 lag panel ----
png(file.path(out_dir, "APDPE_yearly_SIC_vs_krill_lags_2x2.png"),
    width = 2400, height = 2400, res = 300)

par(mfrow = c(2,2), mar = c(5,4,4,2) + 0.1)

plot_krill_lm_black(df_year, "ice_pct", "total_krill",
                    "Yearly: Same Year",
                    "SIC Same Year (%)",
                    paste0("Krill Catch (", krill_units, ")"))

plot_krill_lm_black(df_year, "ice_lag1", "total_krill",
                    "Yearly: 1-year Lag",
                    "SIC 1-year Lag (%)",
                    paste0("Krill Catch (", krill_units, ")"))

plot_krill_lm_black(df_year, "ice_lag2", "total_krill",
                    "Yearly: 2-year Lag",
                    "SIC 2-year Lag (%)",
                    paste0("Krill Catch (", krill_units, ")"))

plot_krill_lm_black(df_year, "ice_lag3", "total_krill",
                    "Yearly: 3-year Lag",
                    "SIC 3-year Lag (%)",
                    paste0("Krill Catch (", krill_units, ")"))

dev.off()
