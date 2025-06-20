#####################################################################################
# Create a dataset of envt predictors for each latlong from 1949-present
#
#####################################################################################

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(glue)
library(sf)
library(terra)

source("scripts/env_data_functions.R")

unique_data <- read_csv("data/unique_latlongs_time.csv") |>
  st_as_sf(coords = c("longitude", "latitude"),
           crs = 4326,
           remove = FALSE)

##
# Spring no3 & po4
##
cmems_rast <- rast("data/rasters/cmems_mod_glo_bgc_my_0.25deg_P1M-m_1750273483465.nc")
cmems_dates <- sort(unique(time(cmems_rast)))

# get values given NAs
cmems_extracted <- extract_with_nas(cmems_rast, unique_data)

cmems_extracted_long <- cmems_extracted |>
  pivot_longer(-c(cell,x,y, latitude, longitude), 
               names_to = c("variable", "depth", "date_idx"),
               names_sep = "_") |>
  select(-depth) |>
  pivot_wider(names_from = variable, values_from = value) |>
  mutate(date = cmems_dates[as.numeric(date_idx)],
         month = month(date),
         year = year(date),
         season = quarter(date, 
                          type = "year.quarter", 
                          fiscal_start = 12),
         season_name = case_when(
           grepl("\\.1", season) ~ "winter",
           grepl("\\.2", season) ~ "spring",
           grepl("\\.3", season) ~ "summer",
           grepl("\\.4", season) ~ "fall",
           .default = "other"
         )) |>
  as_tibble() 
  

 # check
 # ggplot(cmems_extracted_long,
 #        aes(y = si, x = year,
 #            color = latitude,
 #            group = paste(latitude, longitude))) +
 #   geom_line() +
 #   facet_wrap(vars(month)) +
 #   scale_color_viridis_c()

# get annual seasonal means and 
# then average seasonal means for the entire timeseries
cmems_seasonal <- cmems_extracted_long |>
  group_by(x, y, cell, longitude, latitude,
           season, season_name) |>
  summarize(across(no3:si, 
                   \(x) mean(x, na.rm = TRUE))) |>
  ungroup() |>
  mutate(year = gsub("\\.[1-4]", "", season)) |>
  select(-season) |>
  pivot_wider(names_from = season_name, 
              values_from = no3:si,
              names_sep = "_") |> 
  ungroup()

cmes_means <- cmems_seasonal |>
  group_by(x, y, cell, longitude, latitude) |>
  summarize(across(no3_winter:si_fall, 
                   \(x) mean(x, na.rm = TRUE),
                   .names = "{.col}_cellmean")) |> 
  ungroup()

## write them out
write_csv(cmems_seasonal, "data/env_data/cmems_seasonal.csv")
write_csv(cmes_means, "data/env_data/cmems_means.csv")

##
# lag summer mean sst
# summer mean sst
# spring mean sst
##

## with oisst
oisst_monthly_rast <- rast("data/rasters/oisst_monthly_mean_nma.nc")
oisst_times <- time(oisst_monthly_rast)

oisst_monthly <- extract_with_nas(oisst_monthly_rast,
                         unique_data) |>
  pivot_longer(-c(cell,x,y, latitude, longitude),
               values_to = "oisst") |>
  mutate(date_idx = gsub("sst_", "", name) |> as.numeric()) |>
  mutate(times = oisst_times[as.numeric(date_idx)],
         year = floor(times),
         month = round((times-year) * 12 + 1),
         #ok, now make the dates right
         date = ymd(paste(year, month, "01", sep = "-"))) |>
  make_date_to_year_mo_season() |> #see functions
  select(-name)

# get annual seasonal means and 
# then average seasonal means for the entire timeseries
oisst_seasonal <- oisst_monthly |>
  group_by(x, y, cell, longitude, latitude,
           season, season_name) |>
  summarize(across(oisst, mean)) |>
  ungroup() |>
  mutate(year = gsub("\\.[1-4]", "", season)) |>
  select(-season) |>
  pivot_wider(names_from = season_name, 
              values_from = oisst,
              names_prefix = "oisst_")|>
  group_by(x, y, cell, longitude, latitude) |>
  arrange(year) |>
  mutate(lag_oisst_summer = lag(oisst_summer)) |>
  ungroup()

oisst_means <- oisst_seasonal |>
  group_by(x, y, cell, longitude, latitude) |>
  summarize(across(oisst_winter:lag_oisst_summer, 
                   \(x) mean(x, na.rm = TRUE),
                   .names = "{.col}_cellmean"))

## write them out
write_csv(oisst_seasonal, "data/env_data/oisst_seasonal.csv")
write_csv(oisst_means, "data/env_data/oisst_means.csv")

## with hadsst 
had_rast <- rast("data/rasters/HadISST_sst_nwa.nc")
had_time <- time(had_rast)
names(had_rast) <- had_time


had_monthly <- extract_with_nas(had_rast,
                                  unique_data) |>
  pivot_longer(-c(cell,x,y, latitude, longitude),
               values_to = "hadsst",
               names_to = "date") |>
  mutate(date = ymd(date)) |>
  make_date_to_year_mo_season()  |>
  as_tibble()

# get annual seasonal means and 
# then average seasonal means for the entire timeseries
hadsst_seasonal <- 
  make_seasonal(had_monthly, "hadsst") |>
  group_by(x, y, cell, longitude, latitude) |>
  arrange(year) |>
  mutate(lag_hadsst_summer = lag(hadsst_summer)) |>
  ungroup()

hadsst_means <- hadsst_seasonal |>
  group_by(x, y, cell, longitude, latitude) |>
  summarize(across(hadsst_winter:lag_hadsst_summer, 
                   \(x) mean(x, na.rm = TRUE),
                   .names = "{.col}_cellmean")) |>
  ungroup()

## write them out
write_csv(hadsst_seasonal, "data/env_data/hadsst_seasonal.csv")
write_csv(hadsst_means, "data/env_data/hadsst_means.csv")

##
# spring turbidity "data/rasters/erdMH1kd490mday.nc"
# summer turbidity
# winter turbidity
##
turb_rast <- rast("data/rasters/erdMH1kd490mday.nc")

turb_monthly <- 
  make_monthly(turb_rast, 
               unique_data,
               "kd490")

turb_seasonal <- 
  make_seasonal(turb_monthly, "kd490") 

turb_means <- turb_seasonal |>
  group_by(x, y, cell, longitude, latitude) |>
  summarize(across(kd490_winter:kd490_fall, 
                   \(x) mean(x, na.rm = TRUE),
                   .names = "{.col}_cellmean")) |>
  ungroup()

## write out
write_csv(turb_seasonal, "data/env_data/turb_seasonal.csv")
write_csv(turb_means, "data/env_data/turb_means.csv")

##
# lag fall wave height
# winter wave height
##

wave_rast <- rast("data/rasters/reanalysis-era5-single-levels-monthly-means_swh.nc")

wave_monthly <- 
  make_monthly(wave_rast, 
               unique_data,
               "swh")

wave_seasonal <- 
  make_seasonal(wave_monthly, "swh") |>
  mutate(lag_swh_fall = lag(swh_fall)) |>
  ungroup()

wave_means <- wave_seasonal |>
  group_by(x, y, cell, longitude, latitude) |>
  summarize(across(swh_winter:swh_fall, 
                   \(x) mean(x, na.rm = TRUE),
                   .names = "{.col}_cellmean")) |>
  ungroup()

## write out
write_csv(wave_seasonal, "data/env_data/era5_waves_seasonal.csv")
write_csv(wave_means, "data/env_data/era5_waves_means.csv")
