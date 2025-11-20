## ============================================================
## CCAMLR Statistical Bulletin – Krill, Fishery, Toothfish
## Subarea 48.1 / AP SSMUs
## Stacked bar plots
## ============================================================

library(tidyverse)
library(lubridate)

## ------------------ 1. Paths & file names -------------------

data_dir <- "L:/Shared drives/Antarctic Marine Mammals/Krill Data/CCAMLR Statistical Bulletin"

krill_file     <- file.path(data_dir, "AggregatedKrillCatch.csv")
fishery_file   <- file.path(data_dir, "AggregatedFisheryCatch.csv")
toothfish_file <- file.path(data_dir, "AggregatedToothfishLandings.csv")

## ------------------ 2. Load datasets ------------------------

krill_raw     <- readr::read_csv(krill_file, guess_max = 200000)
fishery_raw   <- readr::read_csv(fishery_file, guess_max = 200000)
toothfish_raw <- readr::read_csv(toothfish_file, guess_max = 200000)

## ------------------ 3. KRILL: headers & cleaning ------------

# Columns in AggregatedKrillCatch.csv:
# AKC_ID, Calendar_Year, Month, Group_SSMU_Code, SSMU_Code,
# Krill_Green_Weight, SSMU_Scale_Data_YN

krill <- krill_raw %>%
  rename(
    year         = Calendar_Year,
    month        = Month,
    group_ssmu   = Group_SSMU_Code,
    ssmu         = SSMU_Code,
    krill_weight = Krill_Green_Weight
  ) %>%
  mutate(
    year         = as.integer(year),
    month        = as.integer(month),
    group_ssmu   = as.character(group_ssmu),
    ssmu         = as.character(ssmu),
    krill_weight = as.numeric(krill_weight)
  )

# AP group SSMUs (from the file)
ap_ssmus <- c("APDPE", "APDPW", "APEI", "APPA", "APW", "APBSE", "APBSW", "APE")

krill_ap <- krill %>%
  filter(group_ssmu == "AP",
         ssmu %in% ap_ssmus,
         !is.na(krill_weight))

# Yearly totals by SSMU
krill_ap_year <- krill_ap %>%
  group_by(year, ssmu) %>%
  summarise(total_krill_weight = sum(krill_weight, na.rm = TRUE),
            .groups = "drop")

## ------------------ 3A. KRILL stacked bar (ALL YEARS) -------

p_krill_all_stack <- ggplot(krill_ap_year,
                            aes(x = year,
                                y = total_krill_weight,
                                fill = ssmu)) +
  geom_col(position = "stack") +
  labs(
    title    = "Krill catch (green weight) by AP SSMU",
    subtitle = "Aggregated krill catch – AP group SSMUs (all years)",
    x        = "Year",
    y        = "Krill catch (green weight, tonnes)",
    fill     = "SSMU"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

## ------------------ 3B. KRILL stacked bar (2014–2016) -------

krill_ap_year_1416 <- krill_ap_year %>%
  filter(year >= 2014, year <= 2016)

p_krill_1416_stack <- ggplot(krill_ap_year_1416,
                             aes(x = year,
                                 y = total_krill_weight,
                                 fill = ssmu)) +
  geom_col(position = "stack") +
  labs(
    title    = "Krill catch (green weight) by AP SSMU (2014–2016)",
    subtitle = "Aggregated krill catch – AP group SSMUs",
    x        = "Year",
    y        = "Krill catch (green weight, tonnes)",
    fill     = "SSMU"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

## ------------------ 3C. KRILL – APDPE & APEI only -----------

ap_pair <- c("APDPE", "APEI")

krill_ap_pair_year <- krill_ap_year %>%
  filter(ssmu %in% ap_pair)

# All years
p_krill_ap_pair_all <- ggplot(krill_ap_pair_year,
                              aes(x = year,
                                  y = total_krill_weight,
                                  fill = ssmu)) +
  geom_col(position = "stack") +
  labs(
    title    = "Krill catch – APDPE & APEI",
    subtitle = "Aggregated krill catch – AP group (all years)",
    x        = "Year",
    y        = "Krill catch (green weight, tonnes)",
    fill     = "SSMU"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

# 2014–2016
krill_ap_pair_year_1416 <- krill_ap_pair_year %>%
  filter(year >= 2014, year <= 2016)

p_krill_ap_pair_1416 <- ggplot(krill_ap_pair_year_1416,
                               aes(x = year,
                                   y = total_krill_weight,
                                   fill = ssmu)) +
  geom_col(position = "stack") +
  labs(
    title    = "Krill catch – APDPE & APEI (2014–2016)",
    subtitle = "Aggregated krill catch – AP group",
    x        = "Year",
    y        = "Krill catch (green weight, tonnes)",
    fill     = "SSMU"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

## ------------------ 4. FISHERY catch (Subarea 48.1) ---------
## IMPORTANT: this aggregated fishery file does NOT have SSMU
## codes, only GAR_Code (area). So we can filter to subarea 48.1
## (GAR_Code == '481'), but we cannot split by APDPE/APEI here.

# Columns in AggregatedFisheryCatch.csv:
# AFC_ID, AFE_ID, Flag_CTY_Code, Calendar_Year, Month,
# GAR_Code, Target_TXN_Code, GTY_Code, FAC_Code, VSZ_Code,
# Catch_TXN_Code, Green_Weight

fishery <- fishery_raw %>%
  rename(
    year        = Calendar_Year,
    month       = Month,
    area_code   = GAR_Code,
    catch_weight = Green_Weight
  ) %>%
  mutate(
    year         = as.integer(year),
    month        = as.integer(month),
    area_code    = as.character(area_code),
    catch_weight = as.numeric(catch_weight)
  )

# Subarea 48.1 corresponds to GAR_Code == "481"
fishery_481_year <- fishery %>%
  filter(area_code == "481") %>%
  group_by(year) %>%
  summarise(total_catch_weight = sum(catch_weight, na.rm = TRUE),
            .groups = "drop")

# Stacked bar with a single category (just overall catch)
p_fishery_481_all <- ggplot(fishery_481_year,
                            aes(x = year,
                                y = total_catch_weight)) +
  geom_col() +
  labs(
    title    = "Fishery catch – Subarea 48.1",
    subtitle = "Aggregated fishery catch (all years, GAR_Code = 481)",
    x        = "Year",
    y        = "Catch (green weight, tonnes)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )

fishery_481_year_1416 <- fishery_481_year %>%
  filter(year >= 2014, year <= 2016)

p_fishery_481_1416 <- ggplot(fishery_481_year_1416,
                             aes(x = year,
                                 y = total_catch_weight)) +
  geom_col() +
  labs(
    title    = "Fishery catch – Subarea 48.1 (2014–2016)",
    subtitle = "Aggregated fishery catch (GAR_Code = 481)",
    x        = "Year",
    y        = "Catch (green weight, tonnes)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )

## NOTE: There is no SSMU/APDPE/APEI info in this file, so we
## CANNOT make an APDPE/APEI-only fishery plot from this dataset.

## ------------------ 5. TOOTHFISH landings (Subarea 48.1) ----
## Same story as fishery: we have GAR_Code, but no SSMU codes.

# Columns in AggregatedToothfishLandings.csv:
# ATL_ID, Calendar_Year, GAR_Code, TXN_Code, Flag_CTY_Code,
# Landed_Product_Weight, Estimated_Green_Weight

#Only has one data point for 481 area

# toothfish <- toothfish_raw %>%
#   rename(
#     year                = Calendar_Year,
#     area_code           = GAR_Code,
#     toothfish_weight    = Estimated_Green_Weight  # use estimated green weight
#   ) %>%
#   mutate(
#     year             = as.integer(year),
#     area_code        = as.character(area_code),
#     toothfish_weight = as.numeric(toothfish_weight)
#   )
# 
# toothfish_481_year <- toothfish %>%
#   filter(area_code == "481") %>%
#   group_by(year) %>%
#   summarise(total_toothfish_weight = sum(toothfish_weight, na.rm = TRUE),
#             .groups = "drop")
# 
# p_toothfish_481_all <- ggplot(toothfish_481_year,
#                               aes(x = year,
#                                   y = total_toothfish_weight)) +
#   geom_col() +
#   labs(
#     title    = "Toothfish landings – Subarea 48.1",
#     subtitle = "Aggregated toothfish landings (all years, GAR_Code = 481)",
#     x        = "Year",
#     y        = "Estimated green weight (tonnes)"
#   ) +
#   theme_bw() +
#   theme(
#     plot.title = element_text(face = "bold")
#   )
# 
# toothfish_481_year_1416 <- toothfish_481_year %>%
#   filter(year >= 2014, year <= 2016)
# 
# p_toothfish_481_1416 <- ggplot(toothfish_481_year_1416,
#                                aes(x = year,
#                                    y = total_toothfish_weight)) +
#   geom_col() +
#   labs(
#     title    = "Toothfish landings – Subarea 48.1 (2014–2016)",
#     subtitle = "Aggregated toothfish landings (GAR_Code = 481)",
#     x        = "Year",
#     y        = "Estimated green weight (tonnes)"
#   ) +
#   theme_bw() +
#   theme(
#     plot.title = element_text(face = "bold")
#   )

## Again: no SSMU codes here, so APDPE & APEI cannot be isolated
## for the fishery/toothfish files you provided.

## ------------------ 6. Save all plots -----------------------

plot_dir <- data_dir  # same folder as data

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Krill – all SSMUs
ggsave(file.path(plot_dir, "Krill_AP_SSMU_all_years_stack.png"),
       p_krill_all_stack, width = 10, height = 6, dpi = 300)

ggsave(file.path(plot_dir, "Krill_AP_SSMU_2014_2016_stack.png"),
       p_krill_1416_stack, width = 10, height = 6, dpi = 300)

# Krill – APDPE & APEI
ggsave(file.path(plot_dir, "Krill_APDPE_APEI_all_years_stack.png"),
       p_krill_ap_pair_all, width = 10, height = 6, dpi = 300)

ggsave(file.path(plot_dir, "Krill_APDPE_APEI_2014_2016_stack.png"),
       p_krill_ap_pair_1416, width = 10, height = 6, dpi = 300)

# Fishery – Subarea 48.1
ggsave(file.path(plot_dir, "Fishery_48_1_all_years.png"),
       p_fishery_481_all, width = 10, height = 6, dpi = 300)

ggsave(file.path(plot_dir, "Fishery_48_1_2014_2016.png"),
       p_fishery_481_1416, width = 10, height = 6, dpi = 300)

# Toothfish – Subarea 48.1
# ggsave(file.path(plot_dir, "Toothfish_48_1_all_years.png"),
#        p_toothfish_481_all, width = 10, height = 6, dpi = 300)
# 
# ggsave(file.path(plot_dir, "Toothfish_48_1_2014_2016.png"),
#        p_toothfish_481_1416, width = 10, height = 6, dpi = 300)

## ============================================================
## End of script
## ============================================================
