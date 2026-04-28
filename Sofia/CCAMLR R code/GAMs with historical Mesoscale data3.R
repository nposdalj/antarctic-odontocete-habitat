#1 -- load environmental/krill data --

# -- visualization of data--
#   shape
#   normal distribution
#   fit according model

library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(readxl)
library(sf)
library(units)
library(car)
library(mgcv)
library(gratia)

# -----------------------------
# CONFIGURATIONS
# -----------------------------
env_file   <- file.path("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/Extended_PhysicsReanalysis_BiogeochemistryHindcast", "KGI_60km_mask_with_deseasoned.csv")
catch_file <- file.path("D:/CCAMLR Data/756/R", "756_USA_2026-01-02.Rds")


taxon_code <- "KRI"


start_date_env   <- as.Date("2011-01-05")
start_date_catch <- as.Date("2011-01-01")

buffer_km     <- 60
antarctic_crs <- 3031

sites <- data.frame(
  name = c("KGI"),
  lon  = c(-57.941917),
  lat  = c(-61.457817)
)

env_labels <- c(
  temperature_0  = "Temperature",
  salinity_0     = "Salinity",
  chla_0         = "Chlorophyll",
  productivity_0 = "Productivity",
  mlayer_0       = "Mixed Layer"
)

env_var <- "salinity_anom"
env_var <- "chla_anom"
env_var <- "productivity_anom"
env_var <- "mixed_layer_anom"
env_var <- "temp_anom"
env_var <- "o2_anom"

env_var <- "siconc_mean"
env_var <- "sithick_mean"
env_var <- "EKE_mean"
env_var <- "EKE_mad"
env_var <- "zos_mean"

env_label <- env_labels[[env_var]]

y_label <- env_labels[env_var]


data      <- read.csv(env_file)
data_list <- readRDS(catch_file)

C1       <- data_list[["C1"]]
C1_CATCH <- data_list[["C1_CATCH"]]


buffer_m <- as.numeric(set_units(buffer_km, "km") |> set_units("m"))

sites_sf <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = antarctic_crs)

sites_buffer <- st_buffer(sites_sf, dist = buffer_m)

C1_sf <- C1 %>%
  filter(!is.na(longitude_haul_start), !is.na(latitude_haul_start)) %>%
  st_as_sf(coords = c("longitude_haul_start", "latitude_haul_start"), crs = 4326) %>%
  st_transform(crs = antarctic_crs)

hauls_in_buffer <- st_join(C1_sf, sites_buffer, join = st_within)

haul_sites <- hauls_in_buffer %>%
  st_drop_geometry() %>%
  filter(!is.na(name)) %>%
  select(c1_id, Site = name)

# -----------------------------
# DAILY ENVIRONMENT DATA
# -----------------------------
daily_env_raw <- data %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= start_date_env) %>%
  group_by(date) %>%
  summarise(
    chla_anom  = mean(chla_anom,  na.rm = TRUE),
    o2_anom = mean(o2_anom,  na.rm = TRUE),
    productivity_anom = mean(productivity_anom, na.rm = TRUE),
    temp_anom = mean(temp_anom,  na.rm = TRUE),
    zos_mean = mean(zos_mean, na.rm = TRUE),
    salinity_anom = mean(salinity_anom, na.rm = TRUE),
    mixed_layer_anom = mean(mixed_layer_anom, na.rm = TRUE),
    siconc_mean = mean(siconc_mean, na.rm = TRUE),
    sithick_mean = mean(sithick_mean, na.rm = TRUE),
    EKE_mean = mean(EKE_mean, na.rm = TRUE),
    EKE_mad = mean(EKE_mad, na.rm = TRUE),
    .groups = "drop"
  )



date_range <- seq(min(daily_env_raw$date), max(daily_env_raw$date), by = "day")


daily_krill <- C1_CATCH %>%
  filter(taxon_code == taxon_code) %>%
  left_join(haul_sites, by = "c1_id") %>%
  filter(!is.na(Site)) %>%
  left_join(C1 %>% select(c1_id, date_catchperiod_start), by = "c1_id") %>%
  mutate(date = as.Date(date_catchperiod_start)) %>%
  filter(date >= start_date_catch) %>%
  group_by(date) %>%
  summarise(value = sum(greenweight_caught_kg, na.rm = TRUE), .groups = "drop")

daily_combined <- merge(daily_env_raw, daily_krill, by = "date", all.x = TRUE)


model_data = daily_combined

# how to reshape poison into normal dsitribution (log transform and a square root transform) 
# -------- model_data$value_sqrt - sqrt(model_data$value) / hist(model_data$value_sqrt) --------

# -----------------------------
# TRANSFORM KRILL DATA
# -----------------------------

# Square root transformation
model_data <- model_data %>%
  mutate(
    value_sqrt = sqrt(value),       # square root transform
    value_log  = log1p(value)       # log transform: log(1 + x) to handle zeros
  )

# Visualize the original vs transformed distributions
par(mfrow = c(1,3))  # 1 row, 3 columns
hist(model_data$value, main = "Original Krill Data", xlab = "kg", col = "lightblue", breaks = 30)
hist(model_data$value_sqrt, main = "Square Root Transform", xlab = "sqrt(kg)", col = "lightgreen", breaks = 30)
hist(model_data$value_log, main = "Log Transform", xlab = "log(kg+1)", col = "lightpink", breaks = 30)
par(mfrow = c(1,1))


#2 -- Check/adjust for autocorrelation --
#   Find how much autocorrelation -> adjust
#   How many days is data autocorrelated (find avg in autocorrelation bins)

gam_model <- gam(
  value_sqrt ~ s(EKE_mad, k = 4) + s(temp_anom, k = 4) + s(salinity_anom, k = 4),
  data = model_data,
  family = gaussian()
)

acf_res <- acf(residuals(gam_model), lag.max = 100, plot = TRUE)

CI = ggfortify:::confint.acf(acf_res)
ACFidx = which(acf_res[["acf"]] < CI, arr.ind = TRUE)
ACFval = ACFidx[1]

ACFval <- max(1, ACFval)  # prevent 0 or NA
my_unit <- paste(ACFval, "days")

print(paste("Estimated autocorrelation length (days):", ACFval))

model_data_binned_withNAs <- model_data %>%
  group_by(bin = floor_date(date, unit = my_unit)) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

model_data_binned <- model_data_binned_withNAs %>%
  filter(!is.na(value))

# -----------------------------
# CHECK AUTOCORRELATION
# -----------------------------
residuals_gam <- residuals(gam_model)

acf_res <- acf(residuals_gam, lag.max = 100, plot = TRUE)

# simple GAM model (need to know the shape of your data) with 'value' to the left of the formula and a few environmental variables to the right
# run function acf() on the residuals of your simple gam model 
# ACFval to average into bins


#3  -- Look for covariants (env_data)
#    Independent variables
#   Find similar variables -> remove similar

model_data_binned$siconc_mean[is.na(model_data_binned$siconc_mean)] <- 0





pred <- c("zos_mean", "mixed_layer_anom", "salinity_anom",  
          "temp_anom", 'chla_anom', 'EKE_mad', 'siconc_mean', 'o2_anom')
mod_formula <- paste('value_sqrt', "~", paste(pred, collapse = " + "))
vif_krill <- glm(as.formula(mod_formula),family=gaussian, data = model_data_binned)
vif(vif_krill)

#-----------------------------------------------------------------
#------ Removed zos_mean ----------------------------------------
#------ Removed sithick_mean and EKE_mean, unable to process -----
# ---------------------------------------------------------------


pred <- c("mixed_layer_anom", "salinity_anom",  
          "temp_anom", 'chla_anom', 'EKE_mad', 'siconc_mean', 'o2_anom')
mod_formula <- paste('value_sqrt', "~", paste(pred, collapse = " + "))
vif_krill <- glm(as.formula(mod_formula),family=gaussian, data = model_data_binned)
vif(vif_krill)

#4 -- Build GAM model ---
# side 1: dependent
#side 2: independent
# find p value: stepwise selection (take least significant and remove least significant variable until all variables are significant)

gam_model <- gam(
  value_sqrt ~ 
    s(mixed_layer_anom, k = 4) + 
    s(salinity_anom, k = 4) + 
    s(chla_anom, k = 4) + 
    s(EKE_mad, k = 4) + 
    s(siconc_mean, k = 4) + 
    s(o2_anom, k = 4),
  data = model_data_binned,
  family = gaussian()
)

summary(gam_model)

# -------- removed mixed layer --------

gam_model <- gam(
  value_sqrt ~ 
    s(salinity_anom, k = 4) + 
    s(chla_anom, k = 4) + 
    s(EKE_mad, k = 4) + 
    s(siconc_mean, k = 4) + 
    s(o2_anom, k = 4),
  data = model_data_binned,
  family = gaussian()
)

summary(gam_model)

# -------- removed salinity --------

gam_model <- gam(
  value_sqrt ~ 
    s(chla_anom, k = 4) + 
    s(EKE_mad, k = 4) + 
    s(siconc_mean, k = 4) + 
    s(o2_anom, k = 4),
  data = model_data_binned,
  family = gaussian()
)

summary(gam_model)

# --------- removed o2 --------

gam_model <- gam(
  value_sqrt ~ 
    s(chla_anom, k = 4) + 
    s(EKE_mad, k = 4) + 
    s(siconc_mean, k = 4),
  data = model_data_binned,
  family = gaussian()
)

summary(gam_model)

# ---------- removed chla

gam_model <- gam(
  value_sqrt ~  
    s(EKE_mad, k = 4) + 
    s(siconc_mean, k = 4),
  data = model_data_binned,
  family = gaussian()
)

summary(gam_model)
# ---------------------------------------
# ------------- FINAL GAM ---------------
# ---------------------------------------

model_data_binned$env_combined <- rowMeans(
  model_data_binned[, c("siconc_mean","EKE_mad")],
  na.rm = TRUE
)

draw(gam_model, scales = "fixed")  

summary(gam_model)
par(mfrow = c(2, 3))  
plot(gam_model, pages = 1, residuals = TRUE, pch = 1, shade = TRUE, all.terms = TRUE)

#5 -- Visualize model --

my_unit <- paste(ACFval, "days")

model_data_binned <- model_data %>%
  mutate(date = as.Date(date)) %>%
  group_by(bin = floor_date(date, unit = my_unit)) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

env_vars <- c(
  "siconc_mean",
  "EKE_mad"
)

model_data_binned <- model_data_binned %>%
  mutate(
    env_index = rowMeans(
      scale(across(all_of(env_vars))),
      na.rm = TRUE
    )
  )

gam_model <- gam(
  value_sqrt ~ s(env_index, k = 4),
  data = model_data_binned,
  family = gaussian(),
  method = "REML"
)

new_data <- data.frame(
  env_index = seq(
    min(model_data_binned$env_index, na.rm = TRUE),
    max(model_data_binned$env_index, na.rm = TRUE),
    length.out = 200
  )
)

pred <- predict(gam_model, newdata = new_data, se.fit = TRUE)

new_data <- new_data %>%
  mutate(
    fit = pred$fit,
    se = pred$se.fit,
    upper = fit + 1.96 * se,
    lower = fit - 1.96 * se
  )

gam_sum <- summary(gam_model)
dev_expl <- round(gam_sum$dev.expl * 100, 1)

p_label <- paste0(
  "p = ",
  formatC(gam_sum$s.table[1, "p-value"], format = "f", digits = 7)
)

ggplot(new_data, aes(x = env_index)) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "black", alpha = 0.25) +
  
  geom_line(aes(y = fit),
            color = "black", linewidth = 1.2) +
  
  geom_rug(data = model_data_binned,
           aes(x = env_index),
           inherit.aes = FALSE,
           sides = "b",
           alpha = 0.3) +
  
  annotate("label",
           x = Inf, y = Inf,
           label = p_label,
           hjust = 1.1, vjust = 1.2,
           size = 4,
           fill = "white") +
  
  labs(
    title = paste0(
      "Krill Catch at King George Island\n",
      "Deviance explained = ", dev_expl, "%"
    ),
    x = "Standarised Environmental Gradient",
    y = "Effect on Krill (sqrt)"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.3),
    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    
    plot.margin = margin(t = 10, r = 10, b = 30, l = 10)
  )





