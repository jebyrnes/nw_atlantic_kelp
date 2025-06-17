#####################################################################################
# Functions to add average monthly sea surface temperature and wave intensity data to processed data sets
#
# Author: Amy Smith
# Script Last Updated: Jan 10, 2024
#
# Changelog
#####################################################################################

library(readr)
library(tidyverse)
library(rerddap)

library(ecmwfr)
library(ncdf4)


#function to add sst info to data frame row by row using year, month, lat, and lon
#change: skip rows that are already in final data set

get_SST_monthly <- function(kelp_dat){
  dat <- griddap(
    datasetx = "ncdcOisst21Agg_LonPM180",
    fields = "sst",
    latitude = c(kelp_dat$latitude, kelp_dat$latitude), 
    longitude = c(kelp_dat$longitude, kelp_dat$longitude),
    time = c(paste0(substr(kelp_dat$sasdate, 1, 7),"-01"), paste0(substr(kelp_dat$sasdate, 1, 7),"-31")), #get all sst values for month
    fmt = "csv"
  )
  return(dat)
}



SST_monthly <- function(kelp_dat){
  kelp_dat_processed <- as.data.frame(read_csv("data/env_data/sst"))
  kelp_dat_processed$sasdate <- as.Date(kelp_dat_processed$sasdate)
  
  kelp_dat <- kelp_dat |>
    filter(year(sasdate) > 1981) |>
    anti_join(kelp_dat_processed)
  
  Temp = data.frame(colnames = c("time", "zlev", "latitude", "longitude", "sst"))
  
  SST = data.frame(ncol = 6)
  SST <- kelp_dat_processed
  colnames(SST) = c("sst", "latitude", "longitude", "sasdate", "study", "site")
  
  Temp_kelp_dat = data.frame(colnames = colnames(kelp_dat))
  
  x=nrow(kelp_dat_processed)+1
  for (i in 1:nrow(kelp_dat))  { 
    print(x)
    Temp = get_SST_monthly(kelp_dat[i,]) #get sst for each row in kelp_dat and save as a row in temp
    mean_monthly_sst = mean(Temp$sst, na.rm=TRUE) #get mean monthly sst
    Temp_kelp_dat = kelp_dat[i,]
    n=0
    while(is.na(mean_monthly_sst) & n<5){ #if sst is NaN, move further offshore by 0.25 degrees longitude
      Temp_kelp_dat$longitude = Temp_kelp_dat$longitude + 0.25
      Temp = get_SST_monthly(Temp_kelp_dat)
      mean_monthly_sst = mean(Temp$sst, na.rm=TRUE)
      print(Temp_kelp_dat$longitude)
      n=n+1
    } 
    SST <- add_row(SST,
                   sst = mean_monthly_sst,
                   latitude = kelp_dat$latitude[i],
                   longitude = kelp_dat$longitude[i],
                   sasdate = kelp_dat$sasdate[i],
                   study = kelp_dat$study[i],
                   site = kelp_dat$site[i])
    write.csv(SST, "data/env_data/sst", row.names = FALSE)
    x=x+1
  }
  SST <- SST|>
    distinct()
  write.csv(SST, "data/env_data/sst", row.names = FALSE)
  return(SST)
}






get_SST_winter <- function(kelp_dat){
  dat <- griddap(
    datasetx = "ncdcOisst21Agg_LonPM180",
    fields = "sst",
    latitude = c(kelp_dat$latitude, kelp_dat$latitude), 
    longitude = c(kelp_dat$longitude, kelp_dat$longitude),
    time = c(paste0(kelp_dat$y, "-01-01"), paste0(kelp_dat$y, "-03-31")), #get all sst values for winter
    fmt = "csv"
  )
  return(dat)
}

SST_winter <- function(kelp_dat){
  kelp_dat_processed <- as.data.frame(read_csv("data/env_data/sst_winter"))
  colnames(kelp_dat_processed) <- c("sst", "latitude", "longitude", "y")
  
  kelp_dat <- kelp_dat |>
    mutate(y = year(sasdate)) |>
    distinct(latitude, longitude, y)
  
  kelp_dat <- kelp_dat |>
    filter(y > 1981) |>
    anti_join(kelp_dat_processed) 
  
  Temp = data.frame(colnames = c("time", "zlev", "latitude", "longitude", "sst"))
  
  SST = data.frame(sst = NA, latitude = NA, longitude = NA, y = NA)
  SST <- kelp_dat_processed
  
  Temp_kelp_dat = data.frame(colnames = colnames(kelp_dat))
  
  x=nrow(kelp_dat_processed)+1
  for (i in 1:nrow(kelp_dat))  { 
    print(x)
    Temp = get_SST_winter(kelp_dat[i,]) #get sst for each row in kelp_dat and save as a row in temp
    mean_winter_sst = mean(Temp$sst, na.rm=TRUE) #get mean seasonal sst
    Temp_kelp_dat = kelp_dat[i,]
    n=0
    while(is.na(mean_winter_sst) & n<5){ #if sst is NaN, move further offshore by 0.25 degrees longitude
      Temp_kelp_dat$longitude = Temp_kelp_dat$longitude + 0.25
      Temp = get_SST_winter(Temp_kelp_dat)
      mean_winter_sst = mean(Temp$sst, na.rm=TRUE)
      print(Temp_kelp_dat$longitude)
      n=n+1
    } 
    SST <- add_row(SST,
                   sst = mean_winter_sst,
                   latitude = kelp_dat$latitude[i],
                   longitude = kelp_dat$longitude[i],
                   y = kelp_dat$y[i])
    write.csv(SST, "data/env_data/sst_winter", row.names = FALSE)
    x=x+1
  }
  SST <- SST|>
    distinct()
  write.csv(SST, "data/env_data/sst_winter", row.names = FALSE)
  return(SST)
}








get_SST_spring <- function(kelp_dat){
  dat <- griddap(
    datasetx = "ncdcOisst21Agg_LonPM180",
    fields = "sst",
    latitude = c(kelp_dat$latitude, kelp_dat$latitude), 
    longitude = c(kelp_dat$longitude, kelp_dat$longitude),
    time = c(paste0(kelp_dat$y, "-04-01"), paste0(kelp_dat$y, "-06-31")), #get all sst values for spring
    fmt = "csv"
  )
  return(dat)
}

SST_spring <- function(kelp_dat){
  kelp_dat_processed <- as.data.frame(read_csv("data/env_data/sst_spring.csv", col_names=TRUE))
  colnames(kelp_dat_processed) <- c("sst", "latitude", "longitude", "y")
  
  kelp_dat <- kelp_dat |>
    mutate(y = year(sasdate)) |>
    distinct(latitude, longitude, y)
  
  kelp_dat <- kelp_dat |>
    filter(y > 1981) |>
    anti_join(kelp_dat_processed) 
  
  Temp = data.frame(colnames = c("time", "zlev", "latitude", "longitude", "sst"))
  
  SST = data.frame(sst = NA, latitude = NA, longitude = NA, y = NA)
  SST <- kelp_dat_processed
  
  Temp_kelp_dat = data.frame(colnames = colnames(kelp_dat))
  
  x=nrow(kelp_dat_processed)+1
  for (i in 1:nrow(kelp_dat))  { 
    print(x)
    #kelp_dat$y[i] <- ifelse(month(kelp_dat$sas_date[i])<04, year(kelp_dat$sasdate[i]), year(kelp_dat$sasdate[i])-1)
    Temp = get_SST_spring(kelp_dat[i,]) #get sst for each row in kelp_dat and save as a row in temp
    mean_spring_sst = mean(Temp$sst, na.rm=TRUE) #get mean seasonal sst
    Temp_kelp_dat = kelp_dat[i,]
    n=0
    while(is.na(mean_spring_sst) & n<5){ #if sst is NaN, move further offshore by 0.25 degrees longitude
      Temp_kelp_dat$longitude = Temp_kelp_dat$longitude + 0.25
      Temp = get_SST_spring(Temp_kelp_dat)
      mean_spring_sst = mean(Temp$sst, na.rm=TRUE)
      print(Temp_kelp_dat$longitude)
      n=n+1
    } 
    SST <- add_row(SST,
                   sst = mean_spring_sst,
                   latitude = kelp_dat$latitude[i],
                   longitude = kelp_dat$longitude[i],
                   y = kelp_dat$y[i])
    write.csv(SST, "data/env_data/sst_spring.csv", row.names = FALSE)
    x=x+1
  }
  SST <- SST|>
    distinct()
  write.csv(SST, "data/env_data/sst_spring.csv", row.names = FALSE)
  return(SST)
}


###### for grid calculations

get_SST_spring_grid <- function(kelp_dat){
  dat <- data.frame()
  dat <- griddap(
    datasetx = "ncdcOisst21Agg_LonPM180",
    fields = "sst",
    latitude = c(kelp_dat$latitude_min, kelp_dat$latitude_max), #### changed for grid cells
    longitude = c(kelp_dat$longitude_min, kelp_dat$longitude_max), #### changed for grid cells
    time = c("1981-10-01", "2024-02-01"), #get sst values within cell for entire data set time range
    fmt = "csv"
  )
  return(dat)
}

SST_spring_grid <- function(kelp_dat){
  kelp_dat_processed <- as.data.frame(read_csv("data/env_data/sst_spring_grid.csv", col_names=TRUE))
  colnames(kelp_dat_processed) <- c("sst", "cell_num")
  
  kelp_dat <- kelp_dat |>
    distinct(cell_num, .keep_all=T)
  
  print(kelp_dat)
  
  kelp_dat <- kelp_dat |>
    anti_join(kelp_dat_processed, by=join_by("cell_num" == "cell_num")) 
  
  Temp = data.frame(colnames = c("time", "zlev", "latitude", "longitude", "sst"))
  
  SST = data.frame(sst = NA, cell_num = NA)
  SST <- kelp_dat_processed
  
  Temp_kelp_dat = data.frame(colnames = colnames(kelp_dat))
  
  x=nrow(kelp_dat_processed)+1
  for (i in 1:nrow(kelp_dat))  { 
    print(x)
    Temp = get_SST_spring_grid(kelp_dat[i,]) 
    mean_spring_sst = mean(Temp$sst, na.rm=TRUE) 
    Temp_kelp_dat = kelp_dat[i,]
    
    SST <- add_row(SST,
                   sst = mean_spring_sst,
                   cell_num = kelp_dat$cell_num[i])
    write.csv(SST, "data/env_data/sst_spring_grid.csv", row.names = FALSE)
    x=x+1
  }
  SST <- SST|>
    distinct()
  write.csv(SST, "data/env_data/sst_spring_grid.csv", row.names = FALSE)
  return(SST)
}


########################################## Wave Intensity

get_WI_monthly <- function(kelp_dat){
  source("scripts/wf_keyset.R")
  request <- list(
    "dataset_short_name" = "reanalysis-era5-single-levels-monthly-means",
    "product_type"   = "monthly_averaged_reanalysis",
    "variable"       = "significant_height_of_combined_wind_waves_and_swell",
    "year"           = substr(kelp_dat$sasdate, 1, 4),
    "month"          = substr(kelp_dat$sasdate, 6, 7),
    "time"           = "00:00",
    "area"           = paste0(kelp_dat$latitude,"/",kelp_dat$longitude,"/",kelp_dat$latitude,"/",kelp_dat$longitude),
    "format"         = "netcdf",
    "target"         = "era5-wi.nc"
  )
  dat <- wf_request(
    request = request,
    user = "fdce2bf5-2fe6-4e87-bbd6-ec3f2f9074de"
  )
}



WI_monthly <- function(kelp_dat){
  kelp_dat_processed <- as.data.frame(read_csv("data/env_data/wi"))
  
  kelp_dat <- kelp_dat |>
    anti_join(kelp_dat_processed)
  
  Temp = data.frame(colnames = c("time", "depth", "latitude", "longitude", "wi"))
  
  WI = data.frame(ncol=6)
  WI <- kelp_dat_processed
  colnames(WI) = c("wi", "latitude", "longitude", "sasdate", "study", "site")
  
  Temp_kelp_dat = data.frame(colnames = colnames(kelp_dat))
  
  x=nrow(kelp_dat_processed)+1
  for (i in 1:nrow(kelp_dat))  { 
    print(x)
    dat = get_WI_monthly(kelp_dat[i,]) #get mean monthly wave intensity for each row in kelp_dat
    nc_data <- nc_open(dat)
    wi <- ncvar_get(nc_data, "swh")
    Temp_kelp_dat = kelp_dat[i,]
    n=0
    while(is.na(wi) & n<5){ #if wi is Na, move further offshore by 0.25 degrees longitude
      Temp_kelp_dat$longitude = Temp_kelp_dat$longitude + 0.25
      dat = get_WI_monthly(Temp_kelp_dat)
      nc_data <- nc_open(dat)
      wi <- ncvar_get(nc_data, "swh")
      print(Temp_kelp_dat$longitude)
      n=n+1
    } 
    WI <- add_row(WI,
                  wi = wi,
                  latitude = kelp_dat$latitude[i],
                  longitude = kelp_dat$longitude[i],
                  sasdate = kelp_dat$sasdate[i],
                  study = kelp_dat$study[i],
                  site = kelp_dat$site[i])
    write.csv(WI, "data/env_data/wi", row.names = FALSE)
    x=x+1
  }
  WI <- WI|>
    distinct()
  write.csv(WI, "data/env_data/wi", row.names = FALSE)
  return(WI)
}








get_WI_fall <- function(kelp_dat){
  source("scripts/wf_keyset.R")
  request <- list(
    "dataset_short_name" = "reanalysis-era5-single-levels-monthly-means",
    "product_type"   = "monthly_averaged_reanalysis",
    "variable"       = "significant_height_of_combined_wind_waves_and_swell",
    "year"           = as.numeric(kelp_dat$y) - 1,
    "month"          = kelp_dat$month,
    "time"           = "00:00",
    "area"           = paste0(kelp_dat$latitude,"/",kelp_dat$longitude_new,"/",kelp_dat$latitude,"/",kelp_dat$longitude_new),
    "format"         = "netcdf",
    "target"         = "era5-wi.nc"
  )
  dat <- wf_request(
    request = request,
    user = "fdce2bf5-2fe6-4e87-bbd6-ec3f2f9074de"
  )
}


get_WI_winter <- function(kelp_dat){
  source("scripts/wf_keyset.R")
  
  request <- list(
    "dataset_short_name" = "reanalysis-era5-single-levels-monthly-means",
    "product_type"   = "monthly_averaged_reanalysis",
    "variable"       = "significant_height_of_combined_wind_waves_and_swell",
    "year"           = kelp_dat$y,
    "month"          = kelp_dat$month,
    "time"           = "00:00",
    "area"           = paste0(kelp_dat$latitude,"/",kelp_dat$longitude_new,"/",kelp_dat$latitude,"/",kelp_dat$longitude_new),
    "format"         = "netcdf",
    "target"         = "era5-wi.nc"
  )
  dat <- wf_request(
    request = request,
    user = "fdce2bf5-2fe6-4e87-bbd6-ec3f2f9074de"
  )
}


WI_fall_winter <- function(kelp_dat){
  kelp_dat_processed <- as.data.frame(read_csv("data/env_data/wi_fall_winter.csv"))
  
  kelp_dat <- kelp_dat |>
    mutate(y = year(sasdate)) |>
    mutate(longitude_new = longitude) |>
    distinct(latitude, longitude, longitude_new, y)
  
  kelp_dat <- kelp_dat |>
    anti_join(kelp_dat_processed)
  print(nrow(kelp_dat))
  
  kelp_dat$month = NA
  
  WI = data.frame(ncol=6)
  WI <- kelp_dat_processed
  colnames(WI) = c("wi", "latitude", "longitude", "y")
  
  
  x=nrow(kelp_dat_processed)+1
  for (i in 1:nrow(kelp_dat))  { 
    print(x)
    wi = c()
    a = 1
    for (m in c(10, 11, 12)){ #get fall peak wi
      print(m)
      kelp_dat$month[i] = m
      dat = get_WI_fall(kelp_dat[i,]) #get mean monthly wave intensity for each row in kelp_dat
      nc_data <- nc_open(dat)
      wi[a] <- ncvar_get(nc_data, "swh")
      n=0
      while(is.na(wi[a]) & n<5){ #if wi is Na, move further offshore by 0.25 degrees longitude
        kelp_dat$longitude_new[i] = kelp_dat$longitude_new[i] + 0.25
        dat = get_WI_fall(kelp_dat[i,])
        nc_data <- nc_open(dat)
        wi[a] <- ncvar_get(nc_data, "swh")
        print(kelp_dat$longitude_new[i])
        n=n+1
      } 
      print(wi)
      a = a+1}
    for (m in c(01, 02, 03)){ #get winter peak wi
      print(m)
      kelp_dat$month[i] = m
      dat = get_WI_winter(kelp_dat[i,]) #get mean monthly wave intensity for each row in kelp_dat
      nc_data <- nc_open(dat)
      wi[a] <- ncvar_get(nc_data, "swh")
      n=0
      while(is.na(wi[a]) & n<5){ #if wi is Na, move further offshore by 0.25 degrees longitude
        kelp_dat$longitude_new[i] = as.numeric(kelp_dat$longitude_new[i]) + 0.25
        dat = get_WI_winter(kelp_dat[i,])
        nc_data <- nc_open(dat)
        wi[a] <- ncvar_get(nc_data, "swh")
        print(kelp_dat$longitude_new)
        n=n+1
      } 
      print(wi)
      a = a+1}
    mean_wi = mean(wi, na.rm = T)
    print(mean_wi)
    WI <- add_row(WI,
                  wi = mean_wi,
                  latitude = kelp_dat$latitude[i],
                  longitude = kelp_dat$longitude[i],
                  y = kelp_dat$y[i])
    write.csv(WI, "data/env_data/wi_fall_winter.csv", row.names = FALSE)
    x=x+1
  }
  WI <- WI|>
    distinct()
  write.csv(WI, "data/env_data/wi_fall_winter.csv", row.names = FALSE)
  return(WI)
}




WI_fall_winter_max <- function(kelp_dat){
  kelp_dat_processed <- as.data.frame(read_csv("data/env_data/wi_fall_winter_max.csv"))
  kelp_dat_processed$latitude <- as.numeric(kelp_dat_processed$latitude)
  kelp_dat_processed$longitude <- as.numeric(kelp_dat_processed$longitude)
  kelp_dat_processed$y <- as.numeric(kelp_dat_processed$y)
  
  kelp_dat <- kelp_dat |>
    mutate(y = year(sasdate)) |>
    mutate(longitude_new = longitude) |>
    distinct(latitude, longitude, longitude_new, y)
  
  kelp_dat <- kelp_dat |>
    anti_join(kelp_dat_processed)
  print(nrow(kelp_dat))
  
  kelp_dat$month = NA
  
  WI = data.frame(ncol=6)
  WI <- kelp_dat_processed
  colnames(WI) = c("wi", "latitude", "longitude", "y")
  
  
  x=nrow(kelp_dat_processed)+1
  for (i in 1:nrow(kelp_dat))  { 
    print(x)
    wi = c()
    a = 1
    for (m in c(10, 11, 12)){ #get fall peak wi
      print(m)
      kelp_dat$month[i] = m
      dat = get_WI_fall(kelp_dat[i,]) #get mean monthly wave intensity for each row in kelp_dat
      nc_data <- nc_open(dat)
      wi[a] <- ncvar_get(nc_data, "swh")
      n=0
      while(is.na(wi[a]) & n<5){ #if wi is Na, move further offshore by 0.25 degrees longitude
        kelp_dat$longitude_new[i] = kelp_dat$longitude_new[i] + 0.25
        dat = get_WI_fall(kelp_dat[i,])
        nc_data <- nc_open(dat)
        wi[a] <- ncvar_get(nc_data, "swh")
        print(kelp_dat$longitude_new[i])
        n=n+1
      } 
      print(wi)
      a = a+1}
    for (m in c(01, 02, 03)){ #get winter peak wi
      print(m)
      kelp_dat$month[i] = m
      dat = get_WI_winter(kelp_dat[i,]) #get mean monthly wave intensity for each row in kelp_dat
      nc_data <- nc_open(dat)
      wi[a] <- ncvar_get(nc_data, "swh")
      n=0
      while(is.na(wi[a]) & n<5){ #if wi is Na, move further offshore by 0.25 degrees longitude
        kelp_dat$longitude_new[i] = as.numeric(kelp_dat$longitude_new[i]) + 0.25
        dat = get_WI_winter(kelp_dat[i,])
        nc_data <- nc_open(dat)
        wi[a] <- ncvar_get(nc_data, "swh")
        print(kelp_dat$longitude_new)
        n=n+1
      } 
      print(wi)
      a = a+1}
    max_wi = max(wi)
    print(max_wi)
    WI <- add_row(WI,
                  wi = max_wi,
                  latitude = kelp_dat$latitude[i],
                  longitude = kelp_dat$longitude[i],
                  y = kelp_dat$y[i])
    write.csv(WI, "data/env_data/wi_fall_winter_max.csv", row.names = FALSE)
    x=x+1
  }
  WI <- WI|>
    distinct()
  write.csv(WI, "data/env_data/wi_fall_winter_max.csv", row.names = FALSE)
  return(WI)
}
####################### Turbidity

get_turb_monthly <- function(kelp_dat){
  dat <- griddap(
    datasetx = "erdMH1kd4901day",
    fields = "k490",
    latitude = c(kelp_dat$latitude-0.05, kelp_dat$latitude+0.05), 
    longitude = c(kelp_dat$longitude, kelp_dat$longitude+0.1),
    time = c(paste0(substr(kelp_dat$sasdate, 1, 7),"-01"), paste0(substr(kelp_dat$sasdate, 1, 7),"-31")), #get all sst values for month
    fmt = "csv"
  )
  return(dat)
}


turb_monthly <- function(kelp_dat){
  kelp_dat_processed <- as.data.frame(read_csv("data/env_data/turbidity.csv"))
  kelp_dat_processed$sasdate <- as.Date(kelp_dat_processed$sasdate)
  kelp_dat_processed$longitude <- as.numeric(kelp_dat_processed$longitude)
  kelp_dat_processed$latitude <- as.numeric(kelp_dat_processed$longitude)
  kelp_dat_processed$k490 <- as.numeric(kelp_dat_processed$k490)
  
  kelp_dat <- kelp_dat |>
    filter(year(sasdate) >= 2003) |>
    anti_join(kelp_dat_processed)
  
  Temp = data.frame(colnames = c("time", "latitude", "longitude", "k490"))
  
  Turbidity = data.frame(ncol = 6)
  Turbidity <- kelp_dat_processed
  colnames(Turbidity) = c("k490", "latitude", "longitude", "sasdate", "study", "site")
  
  Temp_kelp_dat = data.frame(colnames = colnames(kelp_dat))
  
  x=nrow(kelp_dat_processed)+1
  for (i in 1:nrow(kelp_dat))  { 
    print(x)
    Temp = get_turb_monthly(kelp_dat[i,]) #get sst for each row in kelp_dat and save as a row in temp
    mean_monthly_turb = mean(Temp$k490, na.rm=TRUE) #get mean monthly sst
    Temp_kelp_dat = kelp_dat[i,]
    n=0
    while(is.na(mean_monthly_turb) & n<5){ #if sst is NaN, move further offshore by 0.25 degrees longitude
      Temp_kelp_dat$longitude = Temp_kelp_dat$longitude + 0.1
      Temp = get_turb_monthly(Temp_kelp_dat)
      mean_monthly_turb = mean(Temp$k490, na.rm=TRUE)
      print(Temp_kelp_dat$longitude)
      n=n+1
    } 
    Turbidity <- add_row(Turbidity,
                         k490 = mean_monthly_turb,
                         latitude = kelp_dat$latitude[i],
                         longitude = kelp_dat$longitude[i],
                         sasdate = kelp_dat$sasdate[i],
                         study = kelp_dat$study[i],
                         site = kelp_dat$site[i])
    write.csv(Turbidity, "data/env_data/turbidity.csv", row.names = FALSE)
    x=x+1
  }
  Turbidity <- Turbidity|>
    distinct()
  write.csv(Turbidity, "data/env_data/turbidity.csv", row.names = FALSE)
  return(Turbidity)
}



get_turb_spring <- function(kelp_dat){
  dat <- griddap(
    datasetx = "erdMH1kd4901day", #4km resolution
    fields = "k490",
    latitude = c(kelp_dat$latitude-0.05, kelp_dat$latitude+0.05), 
    longitude = c(kelp_dat$longitude, kelp_dat$longitude+0.1),
    time = c(paste0(year(kelp_dat$sasdate),"-04-01"), paste0(year(kelp_dat$sasdate),"-06-30")), 
    fmt = "csv"
  )
  return(dat)
}


turb_spring <- function(kelp_dat){
  kelp_dat_processed <- as.data.frame(read_csv("data/env_data/turbidity_spring.csv"))|>
    mutate(y = as.numeric(y))
  
  kelp_dat <- kelp_dat |>
    mutate(y = year(sasdate))|>
    filter(y >= 2003) |>
    anti_join(kelp_dat_processed)
  
  Temp = data.frame(colnames = c("time", "latitude", "longitude", "k490"))
  
  Turbidity = data.frame(matrix(ncol = 4, nrow = 2))
  names(Turbidity) = c("k490", "latitude", "longitude", "y")
  #Turbidity[1,] <- NA
  Turbidity <- kelp_dat_processed
  
  Temp_kelp_dat <- data.frame(colnames = colnames(kelp_dat))
  
  x=nrow(kelp_dat_processed)+1
  for (i in 1:nrow(kelp_dat))  { 
    print(x)
    Temp <- get_turb_spring(kelp_dat[i,]) #get sst for each row in kelp_dat and save as a row in temp
    mean_spring_turb = mean(Temp$k490, na.rm=TRUE) #get mean monthly sst
    Temp_kelp_dat <- kelp_dat[i,]
    n=0
    while(is.na(mean_spring_turb) & n<5){ #if sst is NaN, move further offshore by 0.25 degrees longitude
      Temp_kelp_dat$longitude = Temp_kelp_dat$longitude + 0.1
      Temp = get_turb_monthly(Temp_kelp_dat)
      mean_monthly_turb = mean(Temp$k490, na.rm=TRUE)
      print(Temp_kelp_dat$longitude)
      n=n+1
    } 
    Turbidity <- add_row(Turbidity,
                         k490 = mean_spring_turb,
                         latitude = kelp_dat$latitude[i],
                         longitude = kelp_dat$longitude[i],
                         y = kelp_dat$y[i])
    write.csv(Turbidity, "data/env_data/turbidity_spring.csv", row.names = FALSE)
    x=x+1
  }
  Turbidity <- Turbidity|>
    distinct()
  write.csv(Turbidity, "data/env_data/turbidity_spring.csv", row.names = FALSE)
  return(Turbidity)
}



###### for grid calculations


get_turb_spring_grid <- function(kelp_dat){
  dat <- griddap(
    datasetx = "erdMH1kd4901day", #4km resolution
    fields = "k490",
    latitude = c(kelp_dat$latitude_min, kelp_dat$latitude_max), 
    longitude = c(kelp_dat$longitude_min, kelp_dat$longitude_max),
    time = c("2003-02-01", "2022-06-30"), 
    fmt = "csv"
  )
  return(dat)
}


turb_spring_grid <- function(kelp_dat){
  kelp_dat_processed <- as.data.frame(read_csv("data/env_data/turbidity_spring_grid.csv"))
  
  kelp_dat <- kelp_dat |>
    anti_join(kelp_dat_processed, by = join_by(cell_num == cell))
  
  Temp = data.frame(colnames = c("time", "latitude", "longitude",  "k490"))
  
  Turbidity = data.frame(matrix(ncol = 2, nrow = 2))
  names(Turbidity) = c("k490", "cell")
  Turbidity <- kelp_dat_processed
  
  x=nrow(kelp_dat_processed)+1
  for (i in 1:nrow(kelp_dat))  { 
    print(x)
    Temp <- get_turb_spring_grid(kelp_dat[i,]) #get sst for each row in kelp_dat and save as a row in temp
    mean_spring_turb = mean(Temp$k490, na.rm=TRUE) #get mean monthly sst
    
    Turbidity <- add_row(Turbidity,
                         k490 = mean_spring_turb,
                         cell = kelp_dat$cell[i])
    write.csv(Turbidity, "data/env_data/turbidity_spring_grid.csv", row.names = FALSE)
    x=x+1
  }
  Turbidity <- Turbidity|>
    distinct()
  write.csv(Turbidity, "data/env_data/turbidity_spring_grid.csv", row.names = FALSE)
  return(Turbidity)
}






####### Nitrates:

get_no3_spring <- function(kelp_dat){
  nitrate_data <- read_csv("data/env_data/GLODAPv2.2023_Atlantic_Ocean.csv")|>
    select(G2latitude, G2longitude, G2year, G2nitrate, G2month)
  nitrate_data$G2nitrate[nitrate_data["G2nitrate"] == -9999] <- NA ### replace -9999 with NA in NO3 dataset
  
  nitrate_data <- nitrate_data |>
    mutate(lat_between = between(kelp_dat$latitude, kelp_dat$latitude-0.025, kelp_dat$latitude+0.025),
           lon_between = between(kelp_dat$longitude, kelp_dat$longitude-0.025, kelp_dat$longitude+0.025))|>
    filter(lat_between == TRUE, 
           lon_between == TRUE, #### checks if there is an nitrate observation within 0.025 degrees of data point
           G2year == kelp_dat$y, 
           G2month %in% c(4,5,6))|>
    summarize(mean_nitrate = mean(G2nitrate, na.rm = TRUE))
  return(nitrate_data$mean_nitrate[1])
}

no3_spring <- function(kelp_dat){
  kelp_dat_processed <- as.data.frame(read_csv("data/env_data/no3_spring.csv"))
  kelp_dat_processed$longitude <- as.numeric(kelp_dat_processed$longitude)
  kelp_dat_processed$latitude <- as.numeric(kelp_dat_processed$latitude)
  kelp_dat_processed$y <- as.numeric(kelp_dat_processed$y)
  kelp_dat_processed$no3 <- as.numeric(kelp_dat_processed$no3)
  
  kelp_dat <- kelp_dat |>
    mutate(y = year(sasdate)) |>
    mutate(longitude_new = longitude) |>
    distinct(latitude, longitude, longitude_new, y)
  
  kelp_dat <- kelp_dat |>
    anti_join(kelp_dat_processed)
  print(nrow(kelp_dat))
  
  kelp_dat$month = NA
  
#  NO3 = data.frame(ncol=6)
  NO3 <- kelp_dat_processed
  colnames(NO3) = c("no3", "latitude", "longitude", "y")
  
  x=nrow(kelp_dat_processed)+1
  for (i in 1:nrow(kelp_dat))  { 
    print(x)
    dat = get_no3_spring(kelp_dat[i,]) #get mean monthly no3 for each row in kelp_dat
    n=0
    while(is.na(dat) & n<5){ #if dat is NA, move further offshore by 0.05 degrees longitude 
      kelp_dat$longitude_new[i] = kelp_dat$longitude_new[i] + 0.05
      dat = get_no3_spring(kelp_dat[i,])
      print(kelp_dat$longitude_new[i])
      n=n+1
    } 
    print(dat)
    NO3 <- add_row(NO3,
                   no3 = dat,
                   latitude = kelp_dat$latitude[i],
                   longitude = kelp_dat$longitude[i],
                   y = kelp_dat$y[i])
    write.csv(NO3, "data/env_data/no3_spring.csv", row.names = FALSE)
    x=x+1
  }
  NO3 <- NO3|>
    distinct()
  write.csv(NO3, "data/env_data/no3_spring.csv", row.names = FALSE)
  return(NO3)
}




get_no3_2010 <- function(kelp_dat){
  no3_data_2010 <- as.data.frame(rast("C:/Users/amyls/Downloads/no3_baseline_2000_2018_depthsurf_8486_b388_df7c_U1710889692967.nc"), xy=TRUE) |>
    rename(longitude = x, latitude = y)|>
    mutate(year = 2010) |>
    filter(longitude >= -75 & longitude <= -55 & latitude >= 38 & latitude <= 53)
  
  nitrate_data <- no3_data_2010 |>
    mutate(lat_between = between(kelp_dat$latitude, kelp_dat$latitude-0.025, kelp_dat$latitude+0.025),
           lon_between = between(kelp_dat$longitude, kelp_dat$longitude-0.025, kelp_dat$longitude+0.025))|>
    filter(lat_between == TRUE, 
           lon_between == TRUE)|> #### checks if there is an nitrate observation within 0.025 degrees of data point
    summarize(mean_nitrate = mean(no3_mean, na.rm = TRUE))
  return(nitrate_data$mean_nitrate[1])
}

get_no3_2020 <- function(kelp_dat){
  no3_data_2010 <- as.data.frame(rast("C:/Users/amyls/Downloads/no3_baseline_2000_2018_depthsurf_1e7f_655b_1964_U1710889680958.nc"), xy=TRUE) |>
    rename(longitude = x, latitude = y)|>
    mutate(year = 2010) |>
    filter(longitude >= -75 & longitude <= -55 & latitude >= 38 & latitude <= 53)
  
  nitrate_data <- no3_data_2010 |>
    mutate(lat_between = between(kelp_dat$latitude, kelp_dat$latitude-0.025, kelp_dat$latitude+0.025),
           lon_between = between(kelp_dat$longitude, kelp_dat$longitude-0.025, kelp_dat$longitude+0.025))|>
    filter(lat_between == TRUE, 
           lon_between == TRUE)|> #### checks if there is an nitrate observation within 0.025 degrees of data point
    summarize(mean_nitrate = mean(no3_mean, na.rm = TRUE))
  return(nitrate_data$mean_nitrate[1])
}

no3_decadal <- function(kelp_dat){
  kelp_dat_processed <- as.data.frame(read_csv("data/env_data/no3_decadal.csv"))
  kelp_dat_processed$longitude <- as.numeric(kelp_dat_processed$longitude)
  kelp_dat_processed$latitude <- as.numeric(kelp_dat_processed$latitude)
  kelp_dat_processed$y <- as.numeric(kelp_dat_processed$y)
  kelp_dat_processed$no3 <- as.numeric(kelp_dat_processed$no3)
  
  kelp_dat <- kelp_dat |>
    mutate(y = year(sasdate)) |>
    filter(y >= 2000) |>
    mutate(longitude_new = longitude) |>
    distinct(latitude, longitude, longitude_new, y)
  
  kelp_dat <- kelp_dat |>
    anti_join(kelp_dat_processed)
  print(nrow(kelp_dat))
  
  kelp_dat$month = NA
  
  NO3 = data.frame(ncol=6)
  NO3 <- kelp_dat_processed
  colnames(NO3) = c("no3", "latitude", "longitude", "y")
  
  x=nrow(kelp_dat_processed)+1
  for (i in 1:nrow(kelp_dat))  { 
    print(x)
    dat = ifelse(kelp_dat[i,]$y<= 2010, get_no3_2010(kelp_dat[i,]), get_no3_2020(kelp_dat[i,])) #get mean monthly wave intensity for each row in kelp_dat
    n=0
    while(is.na(dat) & n<5){ #if dat is NA, move further offshore by 0.05 degrees longitude 
      kelp_dat$longitude_new[i] = kelp_dat$longitude_new[i] + 0.05
      dat = ifelse(kelp_dat[i,]$y<= 2010, get_no3_2010(kelp_dat[i,]), get_no3_2020(kelp_dat[i,]))
      print(kelp_dat$longitude_new[i])
      n=n+1
    } 
    print(dat)
    NO3 <- add_row(NO3,
                   no3 = dat,
                   latitude = kelp_dat$latitude[i],
                   longitude = kelp_dat$longitude[i],
                   y = kelp_dat$y[i])
    write.csv(NO3, "data/env_data/no3_decadal.csv", row.names = FALSE)
    x=x+1
  }
  NO3 <- NO3|>
    distinct()
  write.csv(NO3, "data/env_data/no3_decadal.csv", row.names = FALSE)
  return(NO3)
}
######## For Grid Calculations

get_no3_spring_grid <- function(kelp_dat){
  nitrate_data <- read_csv("data/env_data/GLODAPv2.2023_Atlantic_Ocean.csv")|>
    select(G2latitude, G2longitude, G2nitrate)
  nitrate_data$G2nitrate[nitrate_data["G2nitrate"] == -9999] <- NA ### replace -9999 with NA in NO3 dataset
  
  nitrate_data <- nitrate_data |>
    filter(G2latitude >= kelp_dat[["latitude_min"]]
           & G2latitude <= kelp_dat[["latitude_max"]]
           & G2longitude >= kelp_dat[["longitude_min"]]
           & G2longitude <= kelp_dat[["longitude_max"]]) |>
    summarize(mean_nitrate = mean(G2nitrate, na.rm = TRUE))
  return(nitrate_data$mean_nitrate[1])
}


no3_spring_grid <- function(kelp_dat){
  kelp_dat_processed <- as.data.frame(read_csv("data/env_data/no3_spring_grid.csv"))
  kelp_dat_processed$cell <- as.numeric(kelp_dat_processed$cell)
  kelp_dat_processed$no3 <- as.numeric(kelp_dat_processed$no3)
  
  kelp_dat <- kelp_dat |>
    anti_join(kelp_dat_processed, by=join_by(cell_num == cell))
  print(nrow(kelp_dat))
  
  NO3 = data.frame()
  NO3 <- kelp_dat_processed
  #colnames(NO3) = c("no3", "cell")
  
  x=nrow(kelp_dat_processed)+1
  for (i in 1:nrow(kelp_dat))  { 
    print(x)
    dat = get_no3_spring_grid(kelp_dat[i,]) #get mean monthly nitrates for each row in kelp_dat
    print(dat)
    NO3 <- add_row(NO3,
                   no3 = dat,
                   cell = kelp_dat$cell[i])
    write.csv(NO3, "data/env_data/no3_spring_grid.csv", row.names = FALSE)
    x=x+1
  }
  NO3 <- NO3|>
    distinct()
  write.csv(NO3, "data/env_data/no3_spring_grid.csv", row.names = FALSE)
  return(NO3)
}

############# Salinity

get_salt_spring <- function(kelp_dat){
  dat <- griddap(
    datasetx = "NCOM_sfc8_agg",
    fields = "salinity",
    latitude = c(kelp_dat$latitude-0.05, kelp_dat$latitude+0.05), 
    longitude = c(kelp_dat$longitude, kelp_dat$longitude+0.1),
    time = c(paste0(year(kelp_dat$sasdate),"-04-01"), paste0(year(kelp_dat$sasdate),"-06-30")), 
    fmt = "csv"
  )
  return(dat)
}



salt_spring <- function(kelp_dat){
  kelp_dat_processed <- as.data.frame(read_csv("data/env_data/salinity_spring.csv"))|>
    mutate(y = as.numeric(y))
  
  kelp_dat <- kelp_dat |>
    mutate(y = year(sasdate))|>
    anti_join(kelp_dat_processed)
  
  Temp = data.frame(colnames = c("time", "latitude", "longitude", "salinity"))
  
  Salinity = data.frame(matrix(ncol = 4, nrow = 2))
  names(Salinity) = c("salinity", "latitude", "longitude", "y")
  Salinity <- kelp_dat_processed
  
  Temp_kelp_dat <- data.frame(colnames = colnames(kelp_dat))
  
  x=nrow(kelp_dat_processed)+1
  for (i in 1:nrow(kelp_dat))  { 
    print(x)
    Temp <- get_salt_spring(kelp_dat[i,]) #get sst for each row in kelp_dat and save as a row in temp
    mean_spring_salt = mean(Temp$salinity, na.rm=TRUE) #get mean monthly sst
    Temp_kelp_dat <- kelp_dat[i,]
    n=0
    while(is.na(mean_spring_salt) & n<5){ #if sst is NaN, move further offshore by 0.25 degrees longitude
      Temp_kelp_dat$longitude = Temp_kelp_dat$longitude + 0.1
      Temp = get_salt_monthly(Temp_kelp_dat)
      mean_monthly_salt = mean(Temp$salinity, na.rm=TRUE)
      print(Temp_kelp_dat$longitude)
      n=n+1
    } 
    Salinity <- add_row(Salinity,
                        salinity = mean_spring_salt,
                        latitude = kelp_dat$latitude[i],
                        longitude = kelp_dat$longitude[i],
                        y = kelp_dat$y[i])
    write.csv(Salinity, "data/env_data/salinity_spring.csv", row.names = FALSE)
    x=x+1
  }
  Salinity <- Salinity|>
    distinct()
  write.csv(Salinity, "data/env_data/salinity_spring.csv", row.names = FALSE)
  return(Salinity)
}
