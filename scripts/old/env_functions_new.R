library(ecmwfr)
library(terra)
library(sf)
library(lubridate)
library(dplyr)
library(tidyr)
source("scripts/wf_keyset.R")

get_WI_rast <- function(dat){
  
  # get year
  years <- seq(min(dat$year), max(dat$year), 1)
  aoi <- c(max(dat$latitude), min(dat$longitude), min(dat$latitude), max(dat$longitude))
  moi <- 1:12 #months of interest

  request <- list(
    "dataset_short_name" = "reanalysis-era5-single-levels-monthly-means",
    "product_type"   = "monthly_averaged_reanalysis",
    "variable"       = "significant_height_of_combined_wind_waves_and_swell",
    "year"           = years, 
    "month"          = moi,
    "time"           = "00:00",
    "area"           = aoi,# paste0(Data$latitude,"/",Data$longitude_new,"/",Data$latitude,"/",Data$longitude_new),
    "format"         = "netcdf",
    "target"         = "era5-wi.nc"
  )
  
  dl_dat <- wf_request(
    request = request,
    user = "fdce2bf5-2fe6-4e87-bbd6-ec3f2f9074de"
  )
  
  file.copy(dl_dat, "data/rasters/era5-wi.nc", overwrite = TRUE)
  # z <- rast("data/rasters/era5-wi.nc")
  # names(z)
  # 
  return(TRUE)
}

generate_cell_data <- function(dat, rfile = "data/rasters/era5-wi.nc"){
    # get unique coords then turn to an sf object
  dat_unique <- dat |>
    select(longitude, latitude) |>
    group_by(longitude, latitude) |>
    slice(1L) |>
    ungroup() |>
    st_as_sf(coords = c("longitude", "latitude"),
                  crs = 4326)
  
  r <- rast(rfile)
  
  # transform the points to the CRS of the raster
  dat_transformed <- st_transform(dat_unique, crs = crs(r)) |>
    st_buffer(dist = 55000) # 1 cell
  
#  plot(r[[1]])
#  plot(dat_transformed, add = TRUE)
  
  values <- extract(r, dat_transformed, 
                    xy = TRUE,
                    fun="mean", na.rm=TRUE#,
                  #  search_radius = 3 * 55000 #3 cells away, 55km per half degree
                    ) |>
    tidyr::pivot_longer(`swh_valid_time=-662688000`:`swh_valid_time=1669852800`)
  
  # turn names into dates and merge id with point 
  dat_unique_coords <- mutate(dat_unique, ID = 1:n()) |>
    mutate(longitude = st_coordinates(geometry)[,1],
           latitude = st_coordinates(geometry)[,2]) |>
    as_tibble() |> select(-geometry)

  
  ret <- left_join(dat_unique_coords, values) |>
    mutate(dt = gsub("swh_valid_time=", "", name) |> 
             as.numeric() |> as_datetime(),
           year = year(dt),
           month = month(dt)) |>
    select(-name, -dt)
  
  ret
    
}
