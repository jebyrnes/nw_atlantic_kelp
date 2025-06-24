#####################################################################################
# Download or do some light reprocessing on data layers
#
#####################################################################################

library(readr)
library(dplyr)
library(ggplot2)
library(rerddap)
library(terra)
library(lubridate)
library(openeo)

library(ecmwfr)
library(ncdf4)
library(glue)
source("scripts/wf_keyset.R")

##
# unique latlongs
unique_data <- read_csv("data/unique_latlongs_time.csv")
##


##
# oisst
##
print("Getting OISST data from ERDDAP")
for(year in 1982:2024){
  print(glue("OISST in {year}"))
  
  sst <- griddap(
    datasetx = "ncdcOisst21Agg_LonPM180",
    fields = "sst",
    latitude = c(min(unique_data$latitude)-3, max(unique_data$latitude)+3),
    longitude = c(min(unique_data$longitude)-3, max(unique_data$longitude)+3),
    time = c(paste0(year, "-01-01"), paste0(year, "-12-31")), #get all sst values for month
    url = "https://coastwatch.pfeg.noaa.gov/erddap/"
    )

  file.copy(from = sst$summary$filename,
          to = glue("data/rasters/sst/ncdcOisst21Agg_LonPM180_nwa_{year}.nc"),
          overwrite = TRUE)
  
  cache_delete(sst)
  
  #check
 # a <- rast(glue("data/rasters/sst/ncdcOisst21Agg_LonPM180_nwa_{year}.nc"))
  #plot(a[[1]])
  
  # give errdap a break
  Sys.sleep(10)
  
}

# make OISST monthly mean and max
oisst_files <- list.files("data/rasters/sst/",
                          full.names = TRUE)
oisst_all <- rast(oisst_files)

oisst_monthly_mean <- tapp(oisst_all, "yearmonths", mean, na.rm = TRUE)
oisst_monthly_max <- tapp(oisst_all, "yearmonths", max, na.rm = TRUE)

writeCDF(oisst_monthly_mean, "data/rasters/oisst_monthly_mean_nma.nc", 
         unit = "C",
         varname = "sst",
         timename = "time",
         overwrite = TRUE)

writeCDF(oisst_monthly_max, "data/rasters/oisst_monthly_max_nma.nc", 
         unit = "C",
         varname = "sst",
         timename = "time",
         overwrite = TRUE)


#
# z <- as_tibble(oisst_monthly_mean) |> 
#   tidyr::pivot_longer(everything()) |>
#   mutate(date = gsub("ym_", "", name) |>
#            ym()) |>
#   group_by(date, year = year(date), month = month(date)) |>
#   summarize(mean_sst = mean(value),
#             min_sst = min(sst),
#             max_sst = max(sst)) |>
#   group_by(year) |>
#   summarize(max_sst = max(mean_sst))
# 
# plot(max_sst ~ year, data = z, type = "l")

##
# hadsst
##
had <- rast("data/rasters/HadISST_sst.nc")

had <- subset(had, time(had) > as.Date("1949-01-01"))
had <- subset(had, time(had) < as.Date("2023-01-01"))
NAflag(had) <- -1000.00000

nwa_had <- crop(had,
                ext(min(unique_data$longitude)-10, max(unique_data$longitude)+10,
                    min(unique_data$latitude)-10, max(unique_data$latitude)+10))

writeCDF(nwa_had, "data/rasters/HadISST_sst_nwa.nc", 
         unit = "C",
         varname = "sst",
         timename = "time",
         overwrite = TRUE)

##
# no3
##
# https://data.marine.copernicus.eu/product/GLOBAL_ANALYSISFORECAST_BGC_001_028/description
# https://data.marine.copernicus.eu/-/w8d93llpx6
# also has ch and po4 and ph and fe and si 
no3_rast <- rast("data/rasters/cmems_mod_glo_bgc_my_0.25deg_P1M-m_1750273483465.nc")

##
# turbidity
##

# erdMH1kd490mday 2003 - present monthly composite
# Dataset Title: 	Diffuse Attenuation K490, Aqua MODIS, NPP, L3SMI, Global, 4km, Science
# Quality, 2003-present (Monthly Composite)
#https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1kd490mday.graph

turb_dat <- griddap(
  datasetx = "erdMH1kd490mday",
  fields = "k490",
  latitude = c(min(unique_data$latitude)-3, max(unique_data$latitude)+3),
  longitude = c(min(unique_data$longitude)-3, max(unique_data$longitude)+3),
  time = c("2003-01-16T00:00:00Z", "2022-05-16T00:00:00Z"), 
  url = "https://coastwatch.pfeg.noaa.gov/erddap/"
)


file.copy(from = turb_dat$summary$filename,
          to = "data/rasters/erdMH1kd490mday.nc",
          overwrite = TRUE)

cache_delete(turb_dat)

##
# waves
# https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5
##
request <- list(
  "dataset_short_name" = "reanalysis-era5-single-levels-monthly-means",
  "product_type"   = "monthly_averaged_reanalysis",
  "variable"       = "significant_height_of_combined_wind_waves_and_swell",
  "year"           = 1949:2022,#substr(kelp_dat$sasdate, 1, 4),
  month            =1:12,
  "time"           = "00:00",
  "area" =    paste(min(unique_data$latitude)-3, min(unique_data$longitude)-3,
                    max(unique_data$latitude)-3, max(unique_data$longitude)+3, 
                    sep = "/"),

  "data_format"         = "netcdf",
  "target"         = "era5-wi.nc"
)

dat <- wf_request(
  request = request,
  user = "fdce2bf5-2fe6-4e87-bbd6-ec3f2f9074de"
)

file.copy(from = dat,
          to = "data/rasters/reanalysis-era5-single-levels-monthly-means_swh.nc",
          overwrite = TRUE)

# get time right
wi_rast <- rast("data/rasters/reanalysis-era5-single-levels-monthly-means_swh.nc")
#mos <- seq(ymd('1949-01-01'),ymd('2022-12-31'),by='months')
time(wi_rast) <- as.POSIXct(depth(wi_rast))
depth(wi_rast) <- NULL

# write with time info
writeCDF(wi_rast, 
         "data/rasters/reanalysis-era5-single-levels-monthly-means_swh.nc", 
         unit = "m",
         varname = "swh",
         timename = "time",
         overwrite = TRUE)
