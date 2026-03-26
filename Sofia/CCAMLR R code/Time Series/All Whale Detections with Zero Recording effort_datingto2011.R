library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(readxl)
library(sf)
library(units)

# -----------------------------
# CONFIGURATIONS
# -----------------------------
env_file   <- file.path("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Copernicus/Deseasoning", "KGI_deseasoned.csv")
catch_file <- file.path("data", "catch", "756_USA.Rds")

env_var    <- "mlayer_0" # temperature_0, salinity_0, chla_0, productivity_0, mlayer_0, SI_0
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

env_labels <- c(
  temperature_0  = "Temperature",
  salinity_0     = "Salinity",
  chla_0         = "Chlorophyll",
  productivity_0 = "Productivity",
  mlayer_0       = "Mixed Layer"
)

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

# ---------------------------------
# MODEL
# ---------------------------------
model <- lm(as.formula(paste(env_var, "~ value")), data = model_data)
summary_model <- summary(model)

r2   <- summary_model$r.squared
pval <- summary_model$coefficients[2,4]

model_line <- data.frame(
  value = seq(0, max(model_data$value, na.rm = TRUE), length.out = 500)
)
model_line[[env_var]] <- predict(model, newdata = model_line)

# -----------------------------
# PLOT
# -----------------------------
ggplot(model_data, aes(x = value, y = .data[[env_var]])) +
  geom_line(data = model_line, aes(x = value, y = .data[[env_var]]),
            color = "black", linewidth = 1.2) +
  geom_point(size = 2.2, color = "black", alpha = 1) +
  geom_smooth(method = "lm", color = "black", linewidth = 1.2,
              se = TRUE, fill = "grey", alpha = 0.2, fullrange = TRUE) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  labs(title = paste("Krill Catch versus", y_label, "at", sites$name),
       x = "Krill Catch (kg)", y = y_label) +
  annotate("label",
           x = Inf, y = Inf,
           label = paste0("R² = ", round(r2, 3), "\np = ", signif(pval, 3)),
           hjust = 1.1, vjust = 1.3,
           size = 5, label.size = 0.5, fill = "white", color = "black") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text  = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.3),
    axis.line  = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

