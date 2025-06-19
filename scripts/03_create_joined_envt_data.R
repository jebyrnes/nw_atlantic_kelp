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
# what do we want to add?
# spring no3
# spring po4
# lag summer mean sst
# summer sst
# spring sst
# spring turbidity
# summer turbidity
# lag fall wave height
# winter wave height
# lag fall & winter wave height?
##

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
  summarize(across(no3:si, mean)) |>
  ungroup() |>
  mutate(year = gsub("\\.[1-4]", "", season)) |>
  select(-season) |>
  pivot_wider(names_from = season_name, 
              values_from = no3:si,
              names_sep = "_")

cmes_means <- cmems_seasonal |>
  group_by(x, y, cell, longitude, latitude) |>
  summarize(across(no3_winter:si_fall, 
                   mean,
                   .names = "{.col}_cellmean"))

## write them out
write_csv(cmems_seasonal, "data/env_data/cmems_seasonal.csv")
write_csv(cmes_means, "data/env_data/cmems_means.csv")

##
# lag summer mean sst
# summer mean sst
# spring mean sst
##
