# Load data
df <- read.csv("F://APDPE_krill_env_monthly.csv", header = TRUE)

# Convert sea ice concentration to percentage
df$ice_pct <- df$sice_conc_mean_month_mean * 100

# Keep only rows with krill data
df_clean <- subset(df, !is.na(krill_weight))

# Create a color gradient from blue to orange
n <- nrow(df_clean)
cols <- colorRampPalette(c("#40E0D0","#FF6F61"))(n)

# Sort colors according to sea ice concentration so gradient matches values
cols <- cols[rank(df_clean$ice_pct)]

# Scatterplot with larger points and gradient
plot(df_clean$ice_pct,
     df_clean$krill_weight,
     xlab = "Sea Ice Concentration (%)",
     ylab = "Krill Catch",
     pch = 19,
     col = cols,
     cex = 1.5,
     main = " APDPE annual SIC vs Krill Catch")

# Add smooth trend line (loess)
loess_model <- loess(krill_weight ~ ice_pct, data = df_clean)
ice_range <- seq(min(df_clean$ice_pct), max(df_clean$ice_pct), length.out = 200)
lines(ice_range, predict(loess_model, newdata = data.frame(ice_pct = ice_range)),
      col = "royalblue3", lwd = 3)


