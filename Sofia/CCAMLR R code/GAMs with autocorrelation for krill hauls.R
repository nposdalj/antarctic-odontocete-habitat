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

# -----------------------------
# CONFIGURATIONS
# -----------------------------
env_file   <- file.path("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/Deseasoning", "KGI_deseasoned.csv")
catch_file <- file.path("D:/CCAMLR Data/756/R", "756_USA_2026-01-02.Rds")

env_var    <- "chla_0"# temperature_0, salinity_0, chla_0, productivity_0, mlayer_0, SI_0
env_var2   <- "chla_0"
taxon_code <- "KRI"

start_date_env   <- as.Date("2011-01-05")
start_date_catch <- as.Date("2011-01-01")

buffer_km     <- 100
antarctic_crs <- 3031

sites <- data.frame(
  name = c("KGI"),
  lon  = c(-57.941917),
  lat  = c(-61.457817)
)


env_var <- "productivity_0"

env_labels <- c(
  temperature_0  = "Temperature",
  salinity_0     = "Salinity",
  chla_0         = "Chlorophyll",
  productivity_0 = "Productivity",
  mlayer_0       = "Mixed Layer"
)

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
    temperature_0  = mean(sst_none,  na.rm = TRUE),
    salinity_0     = mean(sss_none,  na.rm = TRUE),
    chla_0         = mean(chla_none, na.rm = TRUE),
    productivity_0 = mean(npp_none,  na.rm = TRUE),
    mlayer_0       = mean(mlayer_none, na.rm = TRUE),
    recording_effort = n(),
    .groups = "drop"
  )

date_range <- seq(min(daily_env_raw$date), max(daily_env_raw$date), by = "day")

daily_env_full <- expand.grid(date = date_range) %>%
  mutate(date = as.Date(date)) %>%
  left_join(daily_env_raw, by = "date") %>%
  mutate(recording_effort = ifelse(is.na(recording_effort), 0, recording_effort))


daily_krill <- C1_CATCH %>%
  filter(taxon_code == taxon_code) %>%
  left_join(haul_sites, by = "c1_id") %>%
  filter(!is.na(Site)) %>%
  left_join(C1 %>% select(c1_id, date_catchperiod_start), by = "c1_id") %>%
  mutate(date = as.Date(date_catchperiod_start)) %>%
  filter(date >= start_date_catch) %>%
  group_by(date) %>%
  summarise(value = sum(greenweight_caught_kg, na.rm = TRUE), .groups = "drop")

daily_combined <- merge(daily_env_full, daily_krill, by = "date", all.x = TRUE)

model_data <- daily_combined %>%
  filter(!is.na(value), !is.na(.data[[env_var]]))

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

library(mgcv)

# Build formula safely
gam_formula <- as.formula(
  paste0("value_sqrt ~ s(", env_var, ", k = 4) + s(", env_var2, ", k = 4)")
  )

# Fit GAM
gam_model <- gam(
  formula = gam_formula,
  data = model_data,
  family = gaussian()
)

# Check model ran
summary(gam_model)

# -----------------------------
# CHECK AUTOCORRELATION
# -----------------------------
residuals_gam <- residuals(gam_model)

acf_res <- acf(residuals_gam, lag.max = 100, plot = TRUE)

# -----------------------------
# ESTIMATE AUTOCORRELATION LENGTH
CI = ggfortify:::confint.acf(acf_res)
ACFidx = which(acf_res[["acf"]] < CI, arr.ind=TRUE)
ACFval = ACFidx[1]

print(paste("Estimated autocorrelation length (days):", ACFval))

# Assuming your dataframe is 'df' and date column is 'date_col'
my_unit = paste(ACFval, "days")
model_data_binned <- model_data %>%
  group_by(bin = floor_date(date, unit = my_unit)) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

# simple GAM model (need to know the shape of your data) with 'value' to the left of the formula and a few environmental variables to the right
# run function acf() on the residuals of your simple gam model 
# ACFval to average into bins


#3  -- Look for covariants (env_data)
  #    Independent variables
  #   Find similar variables -> remove similar
pred <- c("mlayer_0", "salinity_0", 
             "temperature_0", 'chla_0','productivity_0')
mod_formula <- paste('value_sqrt', "~", paste(pred, collapse = " + "))
vif_krill <- glm(as.formula(mod_formula),family=gaussian, data = model_data_binned)
vif(vif_krill)

pred <- c("mlayer_0", "salinity_0", 
          "temperature_0", 'chla_0')
mod_formula <- paste('value_sqrt', "~", paste(pred, collapse = " + "))
vif_krill <- glm(as.formula(mod_formula),family=gaussian, data = model_data_binned)
vif(vif_krill)

#4 -- Build GAM model ---
  # side 1: dependent
  #side 2: independent
  # find p value: stepwise selection (take least significant and remove least significant variable until all variables are significant)
gam_model <- gam(
  formula = value_sqrt ~ s(mlayer_0,k=4) + s(salinity_0,k = 4) + s(temperature_0, k = 4) + s(chla_0, k =4),
  data = model_data_binned,
  family = gaussian()
)
summary(gam_model)

gam_model <- gam(
  formula = value_sqrt ~ s(mlayer_0,k=4) + s(salinity_0,k = 4) + s(temperature_0, k = 4),
  data = model_data_binned,
  family = gaussian()
)
summary(gam_model)

gam_model <- gam(
  formula = value_sqrt ~ s(salinity_0,k = 4) + s(temperature_0, k = 4),
  data = model_data_binned,
  family = gaussian()
)
summary(gam_model)

gam_model <- gam(
  formula = value_sqrt ~ s(env_var, k = 4),
  data = model_data_binned,
  family = gaussian()
)
summary(gam_model)

gam_model <- gam(
  formula = value_sqrt ~ env_var,
  data = model_data_binned,
  family = gaussian()
)
summary(gam_model)
plot(gam_model,all.terms = TRUE)
plot(gam_model, residuals = TRUE, pch = 1, all.terms = TRUE)

#5 -- Visualize model -- 

x_var <- env_var
x_label <- env_labels[[x_var]]

pred_vars <- all.vars(delete.response(terms(gam_model)))

new_data <- model_data_binned %>%
  summarise(across(all_of(pred_vars), ~ mean(.x, na.rm = TRUE)))

new_data <- new_data[rep(1, 100), ]

new_data[[x_var]] <- seq(
  min(model_data_binned[[x_var]], na.rm = TRUE),
  max(model_data_binned[[x_var]], na.rm = TRUE),
  length.out = 100
)

pred <- predict(gam_model, newdata = new_data, se.fit = TRUE)

new_data$fit <- pred$fit
new_data$se  <- pred$se.fit

new_data$upper <- new_data$fit + 1.96 * new_data$se
new_data$lower <- new_data$fit - 1.96 * new_data$se

gam_sum <- summary(gam_model)

p_val <- if (!is.null(gam_sum$s.table)) {
  gam_sum$s.table[1, "p-value"]
} else {
  gam_sum$p.table[2, "Pr(>|t|)"]
}

p_label <- if (p_val < 0.001) {
  "p < 0.001"
} else {
  paste0("p = ", round(p_val, 3))
}

label_df <- data.frame(
  x = Inf,
  y = Inf,
  label = p_label
)

label_df$x <- label_df$x * 0.98
label_df$y <- label_df$y * 0.95


ggplot() +
  
  geom_ribbon(data = new_data,
              aes(x = .data[[x_var]], ymin = lower, ymax = upper),
              fill = "black", alpha = 0.3) +
  
  geom_line(data = new_data,
            aes(x = .data[[x_var]], y = fit),
            color = "black", size = 1.2) +
  
  geom_rug(data = model_data_binned,
           aes(x = .data[[x_var]]),
           sides = "b",
           color = "black",
           alpha = 0.8) +
  
  geom_label(
    data = label_df,
    aes(x = x, y = y, label = label),
    hjust = 1.1,
    vjust = 1.2,
    size = 5,
    fill = "white",
    color = "black",
    label.size = 0.5
  ) +
  coord_cartesian(clip = "off") +
  
  labs(
    title = paste("GAM: Krill (sqrt) vs", x_label),
    x = x_label,
    y = "Sqrt(Krill catch)"
  ) +
  
  theme_classic() +
  
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )
