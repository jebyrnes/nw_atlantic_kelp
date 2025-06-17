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

get_SST_monthly <- function(Data){
  dat <- griddap(
    datasetx = "ncdcOisst21Agg_LonPM180",
    fields = "sst",
    latitude = c(Data$latitude, Data$latitude), 
    longitude = c(Data$longitude, Data$longitude),
    time = c(paste0(substr(Data$sasdate, 1, 7),"-01"), paste0(substr(Data$sasdate, 1, 7),"-31")), #get all sst values for month
    fmt = "csv"
  )
  return(dat)
}



SST_monthly <- function(Data){
  Data_processed <- as.data.frame(read_csv("data/env_data/sst"))
  Data_processed$sasdate <- as.Date(Data_processed$sasdate)
  
  Data <- Data |>
    filter(year(sasdate) > 1981) |>
    anti_join(Data_processed)
  
  Temp = data.frame(colnames = c("time", "zlev", "latitude", "longitude", "sst"))
  
  SST = data.frame(ncol = 6)
  SST <- Data_processed
  colnames(SST) = c("sst", "latitude", "longitude", "sasdate", "study", "site")
  
  Temp_Data = data.frame(colnames = colnames(Data))
  
  x=nrow(Data_processed)+1
  for (i in 1:nrow(Data))  { 
    print(x)
    Temp = get_SST_monthly(Data[i,]) #get sst for each row in Data and save as a row in temp
    mean_monthly_sst = mean(Temp$sst, na.rm=TRUE) #get mean monthly sst
    Temp_Data = Data[i,]
    n=0
    while(is.na(mean_monthly_sst) & n<5){ #if sst is NaN, move further offshore by 0.25 degrees longitude
      Temp_Data$longitude = Temp_Data$longitude + 0.25
      Temp = get_SST_monthly(Temp_Data)
      mean_monthly_sst = mean(Temp$sst, na.rm=TRUE)
      print(Temp_Data$longitude)
      n=n+1
    } 
    SST <- add_row(SST,
                   sst = mean_monthly_sst,
                   latitude = Data$latitude[i],
                   longitude = Data$longitude[i],
                   sasdate = Data$sasdate[i],
                   study = Data$study[i],
                   site = Data$site[i])
    write.csv(SST, "data/env_data/sst", row.names = FALSE)
    x=x+1
  }
  SST <- SST|>
    distinct()
  write.csv(SST, "data/env_data/sst", row.names = FALSE)
  return(SST)
}






get_SST_winter <- function(Data){
  dat <- griddap(
    datasetx = "ncdcOisst21Agg_LonPM180",
    fields = "sst",
    latitude = c(Data$latitude, Data$latitude), 
    longitude = c(Data$longitude, Data$longitude),
    time = c(paste0(Data$y, "-01-01"), paste0(Data$y, "-03-31")), #get all sst values for winter
    fmt = "csv"
  )
  return(dat)
}

SST_winter <- function(Data){
  Data_processed <- as.data.frame(read_csv("data/env_data/sst_winter"))
  colnames(Data_processed) <- c("sst", "latitude", "longitude", "y")
  
  Data <- Data |>
    mutate(y = year(sasdate)) |>
    distinct(latitude, longitude, y)
  
  Data <- Data |>
    filter(y > 1981) |>
    anti_join(Data_processed) 
  
  Temp = data.frame(colnames = c("time", "zlev", "latitude", "longitude", "sst"))
  
  SST = data.frame(sst = NA, latitude = NA, longitude = NA, y = NA)
  SST <- Data_processed
  
  Temp_Data = data.frame(colnames = colnames(Data))
  
  x=nrow(Data_processed)+1
  for (i in 1:nrow(Data))  { 
    print(x)
    Temp = get_SST_winter(Data[i,]) #get sst for each row in Data and save as a row in temp
    mean_winter_sst = mean(Temp$sst, na.rm=TRUE) #get mean seasonal sst
    Temp_Data = Data[i,]
    n=0
    while(is.na(mean_winter_sst) & n<5){ #if sst is NaN, move further offshore by 0.25 degrees longitude
      Temp_Data$longitude = Temp_Data$longitude + 0.25
      Temp = get_SST_winter(Temp_Data)
      mean_winter_sst = mean(Temp$sst, na.rm=TRUE)
      print(Temp_Data$longitude)
      n=n+1
    } 
    SST <- add_row(SST,
                   sst = mean_winter_sst,
                   latitude = Data$latitude[i],
                   longitude = Data$longitude[i],
                   y = Data$y[i])
    write.csv(SST, "data/env_data/sst_winter", row.names = FALSE)
    x=x+1
  }
  SST <- SST|>
    distinct()
  write.csv(SST, "data/env_data/sst_winter", row.names = FALSE)
  return(SST)
}








get_SST_spring <- function(Data){
  dat <- griddap(
    datasetx = "ncdcOisst21Agg_LonPM180",
    fields = "sst",
    latitude = c(Data$latitude, Data$latitude), 
    longitude = c(Data$longitude, Data$longitude),
    time = c(paste0(Data$y, "-04-01"), paste0(Data$y, "-06-31")), #get all sst values for spring
    fmt = "csv"
  )
  return(dat)
}

SST_spring <- function(Data){
  Data_processed <- as.data.frame(read_csv("data/env_data/sst_spring.csv", col_names=TRUE))
  colnames(Data_processed) <- c("sst", "latitude", "longitude", "y")
  
  Data <- Data |>
    mutate(y = year(sasdate)) |>
    distinct(latitude, longitude, y)
  
  Data <- Data |>
    filter(y > 1981) |>
    anti_join(Data_processed) 
  
  Temp = data.frame(colnames = c("time", "zlev", "latitude", "longitude", "sst"))
  
  SST = data.frame(sst = NA, latitude = NA, longitude = NA, y = NA)
  SST <- Data_processed
  
  Temp_Data = data.frame(colnames = colnames(Data))
  
  x=nrow(Data_processed)+1
  for (i in 1:nrow(Data))  { 
    print(x)
    #Data$y[i] <- ifelse(month(Data$sas_date[i])<04, year(Data$sasdate[i]), year(Data$sasdate[i])-1)
    Temp = get_SST_spring(Data[i,]) #get sst for each row in Data and save as a row in temp
    mean_spring_sst = mean(Temp$sst, na.rm=TRUE) #get mean seasonal sst
    Temp_Data = Data[i,]
    n=0
    while(is.na(mean_spring_sst) & n<5){ #if sst is NaN, move further offshore by 0.25 degrees longitude
      Temp_Data$longitude = Temp_Data$longitude + 0.25
      Temp = get_SST_spring(Temp_Data)
      mean_spring_sst = mean(Temp$sst, na.rm=TRUE)
      print(Temp_Data$longitude)
      n=n+1
    } 
    SST <- add_row(SST,
                   sst = mean_spring_sst,
                   latitude = Data$latitude[i],
                   longitude = Data$longitude[i],
                   y = Data$y[i])
    write.csv(SST, "data/env_data/sst_spring.csv", row.names = FALSE)
    x=x+1
  }
  SST <- SST|>
    distinct()
  write.csv(SST, "data/env_data/sst_spring.csv", row.names = FALSE)
  return(SST)
}


###### for grid calculations

get_SST_spring_grid <- function(Data){
  dat <- data.frame()
  dat <- griddap(
    datasetx = "ncdcOisst21Agg_LonPM180",
    fields = "sst",
    latitude = c(Data$latitude_min, Data$latitude_max), #### changed for grid cells
    longitude = c(Data$longitude_min, Data$longitude_max), #### changed for grid cells
    time = c("1981-10-01", "2024-02-01"), #get sst values within cell for entire data set time range
    fmt = "csv"
  )
  return(dat)
}

SST_spring_grid <- function(Data){
  Data_processed <- as.data.frame(read_csv("data/env_data/sst_spring_grid.csv", col_names=TRUE))
  colnames(Data_processed) <- c("sst", "cell_num")
  
  Data <- Data |>
    distinct(cell_num, .keep_all=T)
  
  print(Data)
  
  Data <- Data |>
    anti_join(Data_processed, by=join_by("cell_num" == "cell_num")) 
  
  Temp = data.frame(colnames = c("time", "zlev", "latitude", "longitude", "sst"))
  
  SST = data.frame(sst = NA, cell_num = NA)
  SST <- Data_processed
  
  Temp_Data = data.frame(colnames = colnames(Data))
  
  x=nrow(Data_processed)+1
  for (i in 1:nrow(Data))  { 
    print(x)
    Temp = get_SST_spring_grid(Data[i,]) 
    mean_spring_sst = mean(Temp$sst, na.rm=TRUE) 
    Temp_Data = Data[i,]
    
    SST <- add_row(SST,
                   sst = mean_spring_sst,
                   cell_num = Data$cell_num[i])
    write.csv(SST, "data/env_data/sst_spring_grid.csv", row.names = FALSE)
    x=x+1
  }
  SST <- SST|>
    distinct()
  write.csv(SST, "data/env_data/sst_spring_grid.csv", row.names = FALSE)
  return(SST)
}


########################################## Wave Intensity

get_WI_monthly <- function(Data){
  source("scripts/wf_keyset.R")
  request <- list(
    "dataset_short_name" = "reanalysis-era5-single-levels-monthly-means",
    "product_type"   = "monthly_averaged_reanalysis",
    "variable"       = "significant_height_of_combined_wind_waves_and_swell",
    "year"           = substr(Data$sasdate, 1, 4),
    "month"          = substr(Data$sasdate, 6, 7),
    "time"           = "00:00",
    "area"           = paste0(Data$latitude,"/",Data$longitude,"/",Data$latitude,"/",Data$longitude),
    "format"         = "netcdf",
    "target"         = "era5-wi.nc"
  )
  dat <- wf_request(
    request = request,
    user = "fdce2bf5-2fe6-4e87-bbd6-ec3f2f9074de"
  )
}



WI_monthly <- function(Data){
  Data_processed <- as.data.frame(read_csv("data/env_data/wi"))
  
  Data <- Data |>
    anti_join(Data_processed)
  
  Temp = data.frame(colnames = c("time", "depth", "latitude", "longitude", "wi"))
  
  WI = data.frame(ncol=6)
  WI <- Data_processed
  colnames(WI) = c("wi", "latitude", "longitude", "sasdate", "study", "site")
  
  Temp_Data = data.frame(colnames = colnames(Data))
  
  x=nrow(Data_processed)+1
  for (i in 1:nrow(Data))  { 
    print(x)
    dat = get_WI_monthly(Data[i,]) #get mean monthly wave intensity for each row in Data
    nc_data <- nc_open(dat)
    wi <- ncvar_get(nc_data, "swh")
    Temp_Data = Data[i,]
    n=0
    while(is.na(wi) & n<5){ #if wi is Na, move further offshore by 0.25 degrees longitude
      Temp_Data$longitude = Temp_Data$longitude + 0.25
      dat = get_WI_monthly(Temp_Data)
      nc_data <- nc_open(dat)
      wi <- ncvar_get(nc_data, "swh")
      print(Temp_Data$longitude)
      n=n+1
    } 
    WI <- add_row(WI,
                  wi = wi,
                  latitude = Data$latitude[i],
                  longitude = Data$longitude[i],
                  sasdate = Data$sasdate[i],
                  study = Data$study[i],
                  site = Data$site[i])
    write.csv(WI, "data/env_data/wi", row.names = FALSE)
    x=x+1
  }
  WI <- WI|>
    distinct()
  write.csv(WI, "data/env_data/wi", row.names = FALSE)
  return(WI)
}








get_WI_fall <- function(Data){
  source("scripts/wf_keyset.R")
  request <- list(
    "dataset_short_name" = "reanalysis-era5-single-levels-monthly-means",
    "product_type"   = "monthly_averaged_reanalysis",
    "variable"       = "significant_height_of_combined_wind_waves_and_swell",
    "year"           = as.numeric(Data$y) - 1,
    "month"          = Data$month,
    "time"           = "00:00",
    "area"           = paste0(Data$latitude,"/",Data$longitude_new,"/",Data$latitude,"/",Data$longitude_new),
    "format"         = "netcdf",
    "target"         = "era5-wi.nc"
  )
  dat <- wf_request(
    request = request,
    user = "fdce2bf5-2fe6-4e87-bbd6-ec3f2f9074de"
  )
}


get_WI_winter <- function(Data){
  source("scripts/wf_keyset.R")
  
  request <- list(
    "dataset_short_name" = "reanalysis-era5-single-levels-monthly-means",
    "product_type"   = "monthly_averaged_reanalysis",
    "variable"       = "significant_height_of_combined_wind_waves_and_swell",
    "year"           = Data$y,
    "month"          = Data$month,
    "time"           = "00:00",
    "area"           = paste0(Data$latitude,"/",Data$longitude_new,"/",Data$latitude,"/",Data$longitude_new),
    "format"         = "netcdf",
    "target"         = "era5-wi.nc"
  )
  dat <- wf_request(
    request = request,
    user = "fdce2bf5-2fe6-4e87-bbd6-ec3f2f9074de"
  )
}


WI_fall_winter <- function(Data){
  Data_processed <- as.data.frame(read_csv("data/env_data/wi_fall_winter.csv"))
  
  Data <- Data |>
    mutate(y = year(sasdate)) |>
    mutate(longitude_new = longitude) |>
    distinct(latitude, longitude, longitude_new, y)
  
  Data <- Data |>
    anti_join(Data_processed)
  print(nrow(Data))
  
  Data$month = NA
  
  WI = data.frame(ncol=6)
  WI <- Data_processed
  colnames(WI) = c("wi", "latitude", "longitude", "y")
  
  
  x=nrow(Data_processed)+1
  for (i in 1:nrow(Data))  { 
    print(x)
    wi = c()
    a = 1
    for (m in c(10, 11, 12)){ #get fall peak wi
      print(m)
      Data$month[i] = m
      dat = get_WI_fall(Data[i,]) #get mean monthly wave intensity for each row in Data
      nc_data <- nc_open(dat)
      wi[a] <- ncvar_get(nc_data, "swh")
      n=0
      while(is.na(wi[a]) & n<5){ #if wi is Na, move further offshore by 0.25 degrees longitude
        Data$longitude_new[i] = Data$longitude_new[i] + 0.25
        dat = get_WI_fall(Data[i,])
        nc_data <- nc_open(dat)
        wi[a] <- ncvar_get(nc_data, "swh")
        print(Data$longitude_new[i])
        n=n+1
      } 
      print(wi)
      a = a+1}
    for (m in c(01, 02, 03)){ #get winter peak wi
      print(m)
      Data$month[i] = m
      dat = get_WI_winter(Data[i,]) #get mean monthly wave intensity for each row in Data
      nc_data <- nc_open(dat)
      wi[a] <- ncvar_get(nc_data, "swh")
      n=0
      while(is.na(wi[a]) & n<5){ #if wi is Na, move further offshore by 0.25 degrees longitude
        Data$longitude_new[i] = as.numeric(Data$longitude_new[i]) + 0.25
        dat = get_WI_winter(Data[i,])
        nc_data <- nc_open(dat)
        wi[a] <- ncvar_get(nc_data, "swh")
        print(Data$longitude_new)
        n=n+1
      } 
      print(wi)
      a = a+1}
    mean_wi = mean(wi, na.rm = T)
    print(mean_wi)
    WI <- add_row(WI,
                  wi = mean_wi,
                  latitude = Data$latitude[i],
                  longitude = Data$longitude[i],
                  y = Data$y[i])
    write.csv(WI, "data/env_data/wi_fall_winter.csv", row.names = FALSE)
    x=x+1
  }
  WI <- WI|>
    distinct()
  write.csv(WI, "data/env_data/wi_fall_winter.csv", row.names = FALSE)
  return(WI)
}




WI_fall_winter_max <- function(Data){
  Data_processed <- as.data.frame(read_csv("data/env_data/wi_fall_winter_max.csv"))
  Data_processed$latitude <- as.numeric(Data_processed$latitude)
  Data_processed$longitude <- as.numeric(Data_processed$longitude)
  Data_processed$y <- as.numeric(Data_processed$y)
  
  Data <- Data |>
    mutate(y = year(sasdate)) |>
    mutate(longitude_new = longitude) |>
    distinct(latitude, longitude, longitude_new, y)
  
  Data <- Data |>
    anti_join(Data_processed)
  print(nrow(Data))
  
  Data$month = NA
  
  WI = data.frame(ncol=6)
  WI <- Data_processed
  colnames(WI) = c("wi", "latitude", "longitude", "y")
  
  
  x=nrow(Data_processed)+1
  for (i in 1:nrow(Data))  { 
    print(x)
    wi = c()
    a = 1
    for (m in c(10, 11, 12)){ #get fall peak wi
      print(m)
      Data$month[i] = m
      dat = get_WI_fall(Data[i,]) #get mean monthly wave intensity for each row in Data
      nc_data <- nc_open(dat)
      wi[a] <- ncvar_get(nc_data, "swh")
      n=0
      while(is.na(wi[a]) & n<5){ #if wi is Na, move further offshore by 0.25 degrees longitude
        Data$longitude_new[i] = Data$longitude_new[i] + 0.25
        dat = get_WI_fall(Data[i,])
        nc_data <- nc_open(dat)
        wi[a] <- ncvar_get(nc_data, "swh")
        print(Data$longitude_new[i])
        n=n+1
      } 
      print(wi)
      a = a+1}
    for (m in c(01, 02, 03)){ #get winter peak wi
      print(m)
      Data$month[i] = m
      dat = get_WI_winter(Data[i,]) #get mean monthly wave intensity for each row in Data
      nc_data <- nc_open(dat)
      wi[a] <- ncvar_get(nc_data, "swh")
      n=0
      while(is.na(wi[a]) & n<5){ #if wi is Na, move further offshore by 0.25 degrees longitude
        Data$longitude_new[i] = as.numeric(Data$longitude_new[i]) + 0.25
        dat = get_WI_winter(Data[i,])
        nc_data <- nc_open(dat)
        wi[a] <- ncvar_get(nc_data, "swh")
        print(Data$longitude_new)
        n=n+1
      } 
      print(wi)
      a = a+1}
    max_wi = max(wi)
    print(max_wi)
    WI <- add_row(WI,
                  wi = max_wi,
                  latitude = Data$latitude[i],
                  longitude = Data$longitude[i],
                  y = Data$y[i])
    write.csv(WI, "data/env_data/wi_fall_winter_max.csv", row.names = FALSE)
    x=x+1
  }
  WI <- WI|>
    distinct()
  write.csv(WI, "data/env_data/wi_fall_winter_max.csv", row.names = FALSE)
  return(WI)
}
####################### Turbidity

get_turb_monthly <- function(Data){
  dat <- griddap(
    datasetx = "erdMH1kd4901day",
    fields = "k490",
    latitude = c(Data$latitude-0.05, Data$latitude+0.05), 
    longitude = c(Data$longitude, Data$longitude+0.1),
    time = c(paste0(substr(Data$sasdate, 1, 7),"-01"), paste0(substr(Data$sasdate, 1, 7),"-31")), #get all sst values for month
    fmt = "csv"
  )
  return(dat)
}


turb_monthly <- function(Data){
  Data_processed <- as.data.frame(read_csv("data/env_data/turbidity.csv"))
  Data_processed$sasdate <- as.Date(Data_processed$sasdate)
  Data_processed$longitude <- as.numeric(Data_processed$longitude)
  Data_processed$latitude <- as.numeric(Data_processed$longitude)
  Data_processed$k490 <- as.numeric(Data_processed$k490)
  
  Data <- Data |>
    filter(year(sasdate) >= 2003) |>
    anti_join(Data_processed)
  
  Temp = data.frame(colnames = c("time", "latitude", "longitude", "k490"))
  
  Turbidity = data.frame(ncol = 6)
  Turbidity <- Data_processed
  colnames(Turbidity) = c("k490", "latitude", "longitude", "sasdate", "study", "site")
  
  Temp_Data = data.frame(colnames = colnames(Data))
  
  x=nrow(Data_processed)+1
  for (i in 1:nrow(Data))  { 
    print(x)
    Temp = get_turb_monthly(Data[i,]) #get sst for each row in Data and save as a row in temp
    mean_monthly_turb = mean(Temp$k490, na.rm=TRUE) #get mean monthly sst
    Temp_Data = Data[i,]
    n=0
    while(is.na(mean_monthly_turb) & n<5){ #if sst is NaN, move further offshore by 0.25 degrees longitude
      Temp_Data$longitude = Temp_Data$longitude + 0.1
      Temp = get_turb_monthly(Temp_Data)
      mean_monthly_turb = mean(Temp$k490, na.rm=TRUE)
      print(Temp_Data$longitude)
      n=n+1
    } 
    Turbidity <- add_row(Turbidity,
                         k490 = mean_monthly_turb,
                         latitude = Data$latitude[i],
                         longitude = Data$longitude[i],
                         sasdate = Data$sasdate[i],
                         study = Data$study[i],
                         site = Data$site[i])
    write.csv(Turbidity, "data/env_data/turbidity.csv", row.names = FALSE)
    x=x+1
  }
  Turbidity <- Turbidity|>
    distinct()
  write.csv(Turbidity, "data/env_data/turbidity.csv", row.names = FALSE)
  return(Turbidity)
}



get_turb_spring <- function(Data){
  dat <- griddap(
    datasetx = "erdMH1kd4901day", #4km resolution
    fields = "k490",
    latitude = c(Data$latitude-0.05, Data$latitude+0.05), 
    longitude = c(Data$longitude, Data$longitude+0.1),
    time = c(paste0(year(Data$sasdate),"-04-01"), paste0(year(Data$sasdate),"-06-30")), 
    fmt = "csv"
  )
  return(dat)
}


turb_spring <- function(Data){
  Data_processed <- as.data.frame(read_csv("data/env_data/turbidity_spring.csv"))|>
    mutate(y = as.numeric(y))
  
  Data <- Data |>
    mutate(y = year(sasdate))|>
    filter(y >= 2003) |>
    anti_join(Data_processed)
  
  Temp = data.frame(colnames = c("time", "latitude", "longitude", "k490"))
  
  Turbidity = data.frame(matrix(ncol = 4, nrow = 2))
  names(Turbidity) = c("k490", "latitude", "longitude", "y")
  #Turbidity[1,] <- NA
  Turbidity <- Data_processed
  
  Temp_Data <- data.frame(colnames = colnames(Data))
  
  x=nrow(Data_processed)+1
  for (i in 1:nrow(Data))  { 
    print(x)
    Temp <- get_turb_spring(Data[i,]) #get sst for each row in Data and save as a row in temp
    mean_spring_turb = mean(Temp$k490, na.rm=TRUE) #get mean monthly sst
    Temp_Data <- Data[i,]
    n=0
    while(is.na(mean_spring_turb) & n<5){ #if sst is NaN, move further offshore by 0.25 degrees longitude
      Temp_Data$longitude = Temp_Data$longitude + 0.1
      Temp = get_turb_monthly(Temp_Data)
      mean_monthly_turb = mean(Temp$k490, na.rm=TRUE)
      print(Temp_Data$longitude)
      n=n+1
    } 
    Turbidity <- add_row(Turbidity,
                         k490 = mean_spring_turb,
                         latitude = Data$latitude[i],
                         longitude = Data$longitude[i],
                         y = Data$y[i])
    write.csv(Turbidity, "data/env_data/turbidity_spring.csv", row.names = FALSE)
    x=x+1
  }
  Turbidity <- Turbidity|>
    distinct()
  write.csv(Turbidity, "data/env_data/turbidity_spring.csv", row.names = FALSE)
  return(Turbidity)
}



###### for grid calculations


get_turb_spring_grid <- function(Data){
  dat <- griddap(
    datasetx = "erdMH1kd4901day", #4km resolution
    fields = "k490",
    latitude = c(Data$latitude_min, Data$latitude_max), 
    longitude = c(Data$longitude_min, Data$longitude_max),
    time = c("2003-02-01", "2022-06-30"), 
    fmt = "csv"
  )
  return(dat)
}


turb_spring_grid <- function(Data){
  Data_processed <- as.data.frame(read_csv("data/env_data/turbidity_spring_grid.csv"))
  
  Data <- Data |>
    anti_join(Data_processed, by = join_by(cell_num == cell))
  
  Temp = data.frame(colnames = c("time", "latitude", "longitude",  "k490"))
  
  Turbidity = data.frame(matrix(ncol = 2, nrow = 2))
  names(Turbidity) = c("k490", "cell")
  Turbidity <- Data_processed
  
  x=nrow(Data_processed)+1
  for (i in 1:nrow(Data))  { 
    print(x)
    Temp <- get_turb_spring_grid(Data[i,]) #get sst for each row in Data and save as a row in temp
    mean_spring_turb = mean(Temp$k490, na.rm=TRUE) #get mean monthly sst
    
    Turbidity <- add_row(Turbidity,
                         k490 = mean_spring_turb,
                         cell = Data$cell[i])
    write.csv(Turbidity, "data/env_data/turbidity_spring_grid.csv", row.names = FALSE)
    x=x+1
  }
  Turbidity <- Turbidity|>
    distinct()
  write.csv(Turbidity, "data/env_data/turbidity_spring_grid.csv", row.names = FALSE)
  return(Turbidity)
}






####### Nitrates:



get_no3_spring <- function(Data){
  nitrate_data <- read_csv("REU/GLODAPv2.2023_Atlantic_Ocean.csv")|>
    select(G2latitude, G2longitude, G2year, G2nitrate, G2month)
  nitrate_data$G2nitrate[nitrate_data["G2nitrate"] == -9999] <- NA ### replace -9999 with NA in NO3 dataset
  
  nitrate_data <- nitrate_data |>
    mutate(lat_between = between(Data$latitude, Data$latitude-0.025, Data$latitude+0.025),
           lon_between = between(Data$longitude, Data$longitude-0.025, Data$longitude+0.025))|>
    filter(lat_between == TRUE, 
           lon_between == TRUE, #### checks if there is an nitrate observation within 0.025 degrees of data point
           G2year == Data$y, 
           G2month %in% c(4,5,6))|>
    summarize(mean_nitrate = mean(G2nitrate, na.rm = TRUE))
  return(nitrate_data$mean_nitrate[1])
}

no3_spring <- function(Data){
  Data_processed <- as.data.frame(read_csv("data/env_data/no3_spring.csv"))
  Data_processed$longitude <- as.numeric(Data_processed$longitude)
  Data_processed$latitude <- as.numeric(Data_processed$latitude)
  Data_processed$y <- as.numeric(Data_processed$y)
  Data_processed$no3 <- as.numeric(Data_processed$no3)
  
  Data <- Data |>
    mutate(y = year(sasdate)) |>
    mutate(longitude_new = longitude) |>
    distinct(latitude, longitude, longitude_new, y)
  
  Data <- Data |>
    anti_join(Data_processed)
  print(nrow(Data))
  
  Data$month = NA
  
  NO3 = data.frame(ncol=6)
  NO3 <- Data_processed
  colnames(NO3) = c("no3", "latitude", "longitude", "y")
  
  x=nrow(Data_processed)+1
  for (i in 1:nrow(Data))  { 
    print(x)
    dat = get_no3_spring(Data[i,]) #get mean monthly wave intensity for each row in Data
    n=0
    while(is.na(dat) & n<5){ #if dat is NA, move further offshore by 0.05 degrees longitude 
      Data$longitude_new[i] = Data$longitude_new[i] + 0.05
      dat = get_no3_spring(Data[i,])
      print(Data$longitude_new[i])
      n=n+1
    } 
    print(dat)
    NO3 <- add_row(NO3,
                   no3 = dat,
                   latitude = Data$latitude[i],
                   longitude = Data$longitude[i],
                   y = Data$y[i])
    write.csv(NO3, "data/env_data/no3_spring.csv", row.names = FALSE)
    x=x+1
  }
  NO3 <- NO3|>
    distinct()
  write.csv(NO3, "data/env_data/no3_spring.csv", row.names = FALSE)
  return(NO3)
}




get_no3_2010 <- function(Data){
  no3_data_2010 <- as.data.frame(rast("C:/Users/amyls/Downloads/no3_baseline_2000_2018_depthsurf_8486_b388_df7c_U1710889692967.nc"), xy=TRUE) |>
    rename(longitude = x, latitude = y)|>
    mutate(year = 2010) |>
    filter(longitude >= -75 & longitude <= -55 & latitude >= 38 & latitude <= 53)
  
  nitrate_data <- no3_data_2010 |>
    mutate(lat_between = between(Data$latitude, Data$latitude-0.025, Data$latitude+0.025),
           lon_between = between(Data$longitude, Data$longitude-0.025, Data$longitude+0.025))|>
    filter(lat_between == TRUE, 
           lon_between == TRUE)|> #### checks if there is an nitrate observation within 0.025 degrees of data point
    summarize(mean_nitrate = mean(no3_mean, na.rm = TRUE))
  return(nitrate_data$mean_nitrate[1])
}

get_no3_2020 <- function(Data){
  no3_data_2010 <- as.data.frame(rast("C:/Users/amyls/Downloads/no3_baseline_2000_2018_depthsurf_1e7f_655b_1964_U1710889680958.nc"), xy=TRUE) |>
    rename(longitude = x, latitude = y)|>
    mutate(year = 2010) |>
    filter(longitude >= -75 & longitude <= -55 & latitude >= 38 & latitude <= 53)
  
  nitrate_data <- no3_data_2010 |>
    mutate(lat_between = between(Data$latitude, Data$latitude-0.025, Data$latitude+0.025),
           lon_between = between(Data$longitude, Data$longitude-0.025, Data$longitude+0.025))|>
    filter(lat_between == TRUE, 
           lon_between == TRUE)|> #### checks if there is an nitrate observation within 0.025 degrees of data point
    summarize(mean_nitrate = mean(no3_mean, na.rm = TRUE))
  return(nitrate_data$mean_nitrate[1])
}

no3_decadal <- function(Data){
  Data_processed <- as.data.frame(read_csv("data/env_data/no3_decadal.csv"))
  Data_processed$longitude <- as.numeric(Data_processed$longitude)
  Data_processed$latitude <- as.numeric(Data_processed$latitude)
  Data_processed$y <- as.numeric(Data_processed$y)
  Data_processed$no3 <- as.numeric(Data_processed$no3)
  
  Data <- Data |>
    mutate(y = year(sasdate)) |>
    filter(y >= 2000) |>
    mutate(longitude_new = longitude) |>
    distinct(latitude, longitude, longitude_new, y)
  
  Data <- Data |>
    anti_join(Data_processed)
  print(nrow(Data))
  
  Data$month = NA
  
  NO3 = data.frame(ncol=6)
  NO3 <- Data_processed
  colnames(NO3) = c("no3", "latitude", "longitude", "y")
  
  x=nrow(Data_processed)+1
  for (i in 1:nrow(Data))  { 
    print(x)
    dat = ifelse(Data[i,]$y<= 2010, get_no3_2010(Data[i,]), get_no3_2020(Data[i,])) #get mean monthly wave intensity for each row in Data
    n=0
    while(is.na(dat) & n<5){ #if dat is NA, move further offshore by 0.05 degrees longitude 
      Data$longitude_new[i] = Data$longitude_new[i] + 0.05
      dat = ifelse(Data[i,]$y<= 2010, get_no3_2010(Data[i,]), get_no3_2020(Data[i,]))
      print(Data$longitude_new[i])
      n=n+1
    } 
    print(dat)
    NO3 <- add_row(NO3,
                   no3 = dat,
                   latitude = Data$latitude[i],
                   longitude = Data$longitude[i],
                   y = Data$y[i])
    write.csv(NO3, "data/env_data/no3_decadal.csv", row.names = FALSE)
    x=x+1
  }
  NO3 <- NO3|>
    distinct()
  write.csv(NO3, "data/env_data/no3_decadal.csv", row.names = FALSE)
  return(NO3)
}
######## For Grid Calculations

get_no3_spring_grid <- function(Data){
  nitrate_data <- read_csv("data/env_data/GLODAPv2.2023_Atlantic_Ocean.csv")|>
    select(G2latitude, G2longitude, G2nitrate)
  nitrate_data$G2nitrate[nitrate_data["G2nitrate"] == -9999] <- NA ### replace -9999 with NA in NO3 dataset
  
  nitrate_data <- nitrate_data |>
    filter(G2latitude >= Data[["latitude_min"]]
           & G2latitude <= Data[["latitude_max"]]
           & G2longitude >= Data[["longitude_min"]]
           & G2longitude <= Data[["longitude_max"]]) |>
    summarize(mean_nitrate = mean(G2nitrate, na.rm = TRUE))
  return(nitrate_data$mean_nitrate[1])
}


no3_spring_grid <- function(Data){
  Data_processed <- as.data.frame(read_csv("data/env_data/no3_spring_grid.csv"))
  Data_processed$cell <- as.numeric(Data_processed$cell)
  Data_processed$no3 <- as.numeric(Data_processed$no3)
  
  Data <- Data |>
    anti_join(Data_processed, by=join_by(cell_num == cell))
  print(nrow(Data))
  
  NO3 = data.frame()
  NO3 <- Data_processed
  #colnames(NO3) = c("no3", "cell")
  
  x=nrow(Data_processed)+1
  for (i in 1:nrow(Data))  { 
    print(x)
    dat = get_no3_spring_grid(Data[i,]) #get mean monthly nitrates for each row in Data
    print(dat)
    NO3 <- add_row(NO3,
                   no3 = dat,
                   cell = Data$cell[i])
    write.csv(NO3, "data/env_data/no3_spring_grid.csv", row.names = FALSE)
    x=x+1
  }
  NO3 <- NO3|>
    distinct()
  write.csv(NO3, "data/env_data/no3_spring_grid.csv", row.names = FALSE)
  return(NO3)
}

############# Salinity

get_salt_spring <- function(Data){
  dat <- griddap(
    datasetx = "NCOM_sfc8_agg",
    fields = "salinity",
    latitude = c(Data$latitude-0.05, Data$latitude+0.05), 
    longitude = c(Data$longitude, Data$longitude+0.1),
    time = c(paste0(year(Data$sasdate),"-04-01"), paste0(year(Data$sasdate),"-06-30")), 
    fmt = "csv"
  )
  return(dat)
}



salt_spring <- function(Data){
  Data_processed <- as.data.frame(read_csv("data/env_data/salinity_spring.csv"))|>
    mutate(y = as.numeric(y))
  
  Data <- Data |>
    mutate(y = year(sasdate))|>
    anti_join(Data_processed)
  
  Temp = data.frame(colnames = c("time", "latitude", "longitude", "salinity"))
  
  Salinity = data.frame(matrix(ncol = 4, nrow = 2))
  names(Salinity) = c("salinity", "latitude", "longitude", "y")
  Salinity <- Data_processed
  
  Temp_Data <- data.frame(colnames = colnames(Data))
  
  x=nrow(Data_processed)+1
  for (i in 1:nrow(Data))  { 
    print(x)
    Temp <- get_salt_spring(Data[i,]) #get sst for each row in Data and save as a row in temp
    mean_spring_salt = mean(Temp$salinity, na.rm=TRUE) #get mean monthly sst
    Temp_Data <- Data[i,]
    n=0
    while(is.na(mean_spring_salt) & n<5){ #if sst is NaN, move further offshore by 0.25 degrees longitude
      Temp_Data$longitude = Temp_Data$longitude + 0.1
      Temp = get_salt_monthly(Temp_Data)
      mean_monthly_salt = mean(Temp$salinity, na.rm=TRUE)
      print(Temp_Data$longitude)
      n=n+1
    } 
    Salinity <- add_row(Salinity,
                        salinity = mean_spring_salt,
                        latitude = Data$latitude[i],
                        longitude = Data$longitude[i],
                        y = Data$y[i])
    write.csv(Salinity, "data/env_data/salinity_spring.csv", row.names = FALSE)
    x=x+1
  }
  Salinity <- Salinity|>
    distinct()
  write.csv(Salinity, "data/env_data/salinity_spring.csv", row.names = FALSE)
  return(Salinity)
}
