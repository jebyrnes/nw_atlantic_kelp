#' -----------------------------------------------
#' Analyze trends in drivers at sites
#'
#' -----------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(performance)
library(parameters)
library(lubridate)
library(modelbased)
library(glmmTMB)
library(broom.mixed)
library(car)
library(sf)
library(rnaturalearth)
library(terra)

theme_set(theme_bw(base_size = 12))
source("scripts/add_ecoregions_to_data.R")

## create coastline buffer

unique_latlong <- read_csv("data/unique_latlongs_time.csv") |>
  st_as_sf(crs = 4326,
           coords = c("longitude", "latitude"))

aoi <- unique_latlong |>
  st_bbox() + c(-3, -3, 3, 3) # some buffer

##
# Map of data
##

coastline <-  ne_states(country = c("United States of America", "Canada"), 
                        returnclass = "sf") |>
  st_crop(aoi)

coastline_buffer <- coastline |>
  summarize()  |>   
  st_buffer(dist = 5000)


coastline_buffer_2km <- coastline |>
  summarize()  |>   
  st_buffer(dist = 2000)

coastline_buffer_1km <- coastline |>
  summarize()  |>   
  st_buffer(dist = 1000)

ggplot(coastline_buffer) +
  geom_sf()

###
## function to extract
###
make_buffered_dat <- 
  function(a_rast,
           coastline_buffer,
           layer_name,
           return_buffered_rast = FALSE){
    
    # debug - oh, wait, useful!
    if(return_buffered_rast){
      return(a_rast |> mask(coastline_buffer))
    }
    
    # make a dataset
    names(a_rast) <- time(a_rast)
    
    dat_extracted <- extract(a_rast,
                             xy = TRUE,
                             coastline_buffer,
                             cells = TRUE) |>
      tidyr::pivot_longer(-c(ID, cell, x, y),
                          names_to = "date",
                          values_to = "values") |>
      filter(!is.na(values)) |>
      mutate(date = ymd(date),
             year = year(date),
             season = quarter(date, type = "year.quarter", fiscal_start = 12),
             season = case_when(
               grepl("\\.1", season) ~ "winter",
               grepl("\\.2", season) ~ "spring",
               grepl("\\.3", season) ~ "summer",
               grepl("\\.4", season) ~ "fall",
               .default = "other"
             )) |>
      group_by(ID, cell, x, y, year, season) |>
      summarize(mean_value = mean(values)) |>
      ungroup() |>
      tidyr::pivot_wider(names_from = season, 
                         values_from = mean_value,
                         names_glue = "layer_{season}")
    
    #return
    dat_extracted |>
      rename_with(~gsub("layer", layer_name, .x)) |>
      add_ecoregions_to_data()
  }

  
##
# hadsst
##

hadsst <- rast("data/rasters/HadISST_sst_nwa.nc") |>
  crop(aoi) |>
  mask(coastline_buffer)

ggplot() +
  geom_sf(data = coastline) +
  tidyterra::geom_spatraster(data = hadsst[[762]]) +
  scale_fill_distiller(palette = "RdBu", 
                       na.value = NA) +
  coord_sf(expand = FALSE, 
           xlim = c(aoi[1], aoi[3]),
           ylim = c(aoi[2], aoi[4])) +
  labs(fill = "ºC", subtitle = "June 2012 HADSST within 5km of coast")

ggsave("figures/hadsst_coastline.jpg", width = 4, height = 5)

# make a dataset
hadsst_extracted <- make_buffered_dat(hadsst,
                              coastline_buffer,
                              layer_name = "hadsst")

ggplot(hadsst_extracted, 
       aes(x = year, y = hadsst_summer, 
           color = y, group = cell)) +
  geom_line()


hadsst_summer_mod <- glmmTMB(
  hadsst_summer ~ 
    year  + (year|cell),
  data = hadsst_extracted
)

hadsst_spring_mod <- glmmTMB(
  hadsst_spring ~ 
    year + (1+year|cell),
  data = hadsst_extracted
)

##
# oisst
##
oisst <- rast("data/rasters/oisst_monthly_mean_nma.nc")
oisst_times <- time(oisst)
oisst_year <- floor(oisst_times)
oisst_month <- round((oisst_times-oisst_year) * 12 + 1)
#ok, now make the dates right
time(oisst) <- ymd(paste(oisst_year, oisst_month, "01", sep = "-"))

oisst_dat <- make_buffered_dat(oisst,
                               coastline_buffer,
                               layer_name = "oisst")


oisst_coast_rast <- make_buffered_dat(oisst,
                               coastline_buffer,
                               layer_name = "oisst",
                               return_buffered_rast=TRUE)

## Plots

ggplot() +
  geom_sf(data = coastline) +
  tidyterra::geom_spatraster(data = oisst_coast_rast[[366]]) +
  scale_fill_distiller(palette = "RdBu", 
                       na.value = NA) +
  coord_sf(expand = FALSE, 
           xlim = c(aoi[1], aoi[3]),
           ylim = c(aoi[2], aoi[4])) +
  labs(fill = "ºC", subtitle = "June 2012 OISST within 5km of coast")
ggsave("figures/oisst_coastline_mask.jpg", width = 5, height = 4)

ggplot(oisst_dat,
       aes(x = year, y = oisst_summer, 
           group = cell, color = y)) +
  geom_line(alpha = 0.4) +
 # geom_smooth(method = "lm", fill = NA, size = 0.4) +
  scale_color_viridis_b(option = "F", direction = -1) +
  guides(color = "none")+
  labs(x = "", y = "OISST ºC", subtitle = "Summer Average")

ggsave("figures/oisst_timeseries.jpg", width = 5, height = 4)

ggplot(oisst_dat,
       aes(x = year, y = oisst_spring, 
           group = cell, color = y)) +
  geom_line(alpha = 0.1) +
  geom_smooth(method = "lm", fill = NA, size = 0.4) +
  scale_color_viridis_b(option = "F", direction = -1) +
  guides(color = "none")+
  labs(x = "", y = "OISST ºC", subtitle = "Spring Average")


oisst_summer_mod <- glmmTMB(
  oisst_summer ~ 
    year + (1+year|cell),
  data = oisst_dat)

oisst_spring_mod <- glmmTMB(
  oisst_spring ~ 
    year + (1+year|cell),
  data = oisst_dat)

##
# waves
## 

wave_rast <- rast("data/rasters/reanalysis-era5-single-levels-monthly-means_swh.nc")



wave_dat <- make_buffered_dat(wave_rast,
                              coastline_buffer,
                              layer_name = "swh")


ggplot(wave_dat,
       aes(x = year, y = swh_fall, 
           group = cell, color = y)) +
  geom_line(alpha = 0.5) +
  #geom_smooth(method = "lm", fill = NA, size = 0.5) +
  scale_color_viridis_b(option = "H", direction = -1) +
  guides(color = "none")+
  labs(x = "", y = "Significant Wave Height (m)", subtitle = "Fall Average")
ggsave("figures/fall_wave_timeseries.jpg", width = 5, height = 4)


ggplot(wave_dat,
       aes(x = year, y = swh_winter, 
           group = cell, color = y)) +
  geom_line(alpha = 0.5) +
  #geom_smooth(method = "lm", fill = NA, size = 0.5) +
  scale_color_viridis_b(option = "H", direction = -1) +
  guides(color = "none")+
  labs(x = "", y = "Significant Wave Height (m)", subtitle = "Winter Average")


swh_fall_mod <- glmmTMB(
  swh_fall ~ 
    year + (year |cell),
  data = wave_dat
)

# swh in winter
swh_winter_mod <- glmmTMB(
  swh_winter ~ 
    year + (1 |cell),
  data = wave_dat)



##
# no3
##

cmems_rast <- rast("data/rasters/cmems_mod_glo_bgc_my_0.25deg_P1M-m_1750273483465.nc")

no3_rast <- cmems_rast["no3"]

no3_dat <- make_buffered_dat(no3_rast,
                             coastline_buffer,
                             layer_name = "no3")



ggplot() +
  geom_sf(data = coastline) +
  tidyterra::geom_spatraster(data = no3_rast |>
                               subset(time(no3_rast)==as.POSIXct("2012-04-01", tz = "UTC")) |>
                               mask(coastline_buffer)) +
  scale_fill_distiller(palette = "GnBu", 
                       na.value = NA,
                       transform = "log10") +
  geom_sf(data = coastline , fill = "lightgreen") +
  labs(subtitle = "April 2012 average CMEMS NO3",
       fill = expression(paste("mmol/", m^3)))+
  coord_sf(expand = FALSE, 
           xlim = c(aoi[1], aoi[3]),
           ylim = c(aoi[2], aoi[4])) 
ggsave("figures/no3_coastline_mask.jpg", width = 5, height = 4)


ggplot(no3_dat,
       aes(x = year, y = no3_spring, 
           group = cell, color = cell)) +
  geom_line(alpha = 0.5) +
  #geom_smooth(method = "lm", fill = NA, size = 0.5) +
  scale_color_viridis_b(option = "H", direction = -1) +
  guides(color = "none")+
  labs(x = "", y = "Nitrate", subtitle = "Spring Average")
ggsave("figures/spring_no3_timeseries.jpg", width = 5, height = 4)



no3_spring_mod <- glmmTMB(
  no3_spring ~ 
    year + (year| cell),
  data = no3_dat)

estimate_relation(no3_spring_mod, 
                  include_random = TRUE) |>
  plot(show_data = TRUE) + 
  guides(color = "none", fill = "none")

no3_winter_mod <- glmmTMB(
  no3_winter ~ 
    year + (year| cell),
  data = no3_dat)


##
# kd490
##

kd490_rast <- rast("data/rasters/erdMH1kd490mday.nc")

turb_dat <- make_buffered_dat(kd490_rast,
                              coastline_buffer,
                              layer_name = "kd490") 

ggplot(turb_dat,
       aes(x = year, y = kd490_spring, 
           group = cell, color = y)) +
  geom_line(alpha = 0.2) +
  #geom_smooth(method = "lm", fill = NA, size = 0.5) +
  scale_color_viridis_b(option = "H", direction = -1) +
  guides(color = "none")+
  labs(x = "", y = "Turbidity (kd490)", subtitle = "Spring Average")


kd490_spring_mod <- glmmTMB(
  kd490_spring ~ 
    year_c + (year |cell),
  REML = TRUE,
  data = turb_dat |> mutate(year_c = year-mean(year)))

###
# Put it all together
###
get_year_coef <- function(mod, 
                          coef = "year", 
                          name = NA,
                          add_ran_sd = TRUE){
  
  out <- tidy(mod) |>
    filter(grepl(coef, term)) |>
    filter(!grepl("cor_", term))
  
  if(add_ran_sd){
    # add them together
    out <- out |> mutate(
      std.error = ifelse(is.na(std.error), 
                         estimate, 
                         std.error),
      std.error = sum(std.error))
  }
  out |> filter(term == coef) |>
    select(estimate, std.error, p.value) |>
    mutate(model = name) |>
    relocate(model, estimate, std.error, p.value) |>
    rename(
      ` ` = model,
      coefficient = estimate,
      se = std.error,
      p = p.value
    )
}

envt_tab <- bind_rows(
  get_year_coef(hadsst_summer_mod, name = "HADSST summer"),
  get_year_coef(hadsst_spring_mod, name = "HADSST spring"),
  get_year_coef(oisst_summer_mod, name = "OISST summer"),
  get_year_coef(oisst_spring_mod, name = "OISST spring"),
  get_year_coef(swh_fall_mod, name = "SWH fall"),
  get_year_coef(swh_winter_mod, name = "SWH winter"),
  get_year_coef(no3_winter_mod, name = "NO3 winter"),
  get_year_coef(no3_spring_mod, name = "NO3 spring"),
  get_year_coef(kd490_spring_mod, name = "Turbidity spring", coef = "year_c")
)

write_csv(envt_tab, "tables/envt_coefs.csv")

gt::gt(envt_tab) |> 
  gt::fmt_number(decimals=3) |>
  gt::gtsave("tables/envt_coefs_gt.docx")

## plot

ggplot(envt_tab |>
         rename(Driver = ` `) |>
         mutate(Driver = forcats::fct_inorder(Driver)),
       aes(x = Driver, y = coefficient,
           ymin = coefficient - 2*se,
           ymax = coefficient + 2*se)) +
  geom_pointrange() +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "",
       y = "rate of change per year") +
  geom_hline(yintercept = 0,
             color = "red",
             lty = 2)

ggsave("figures/envt_tseries_change.jpg",
       width = 6, height = 5)
