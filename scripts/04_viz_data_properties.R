#' -------------------------------------------
#' Viz data properties of kelptime
#' and envt data
#' -------------------------------------------

library(readr)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(wesanderson)
library(sf)
library(terra)
library(tidyterra)

pal <- wes_palette("Zissou1", 4, type = "continuous")

theme_set(theme_bw(base_size = 14))

nwa_data <- read_csv("data/nwa_with_env.csv") |>
  mutate(ecoregion = factor(ecoregion, levels = 
                              c("Virginian",
                                "Gulf of Maine/Bay of Fundy",
                                "Scotian Shelf",
                                "Gulf of St. Lawrence - Eastern Scotian Shelf"
                              )))

unique_latlong <- read_csv("data/unique_latlongs_time.csv") |>
  st_as_sf(crs = 4326,
           coords = c("longitude", "latitude"))

aoi <- unique_latlong |>
  st_bbox() + c(-1, -1, 1, 1) # some buffer

##
# Map of data
##

coastline <- ne_states(country = c("United States of America", "Canada"), 
                       returnclass = "sf") |>
  st_crop(aoi)

basemap <- ggplot() +
  geom_sf(data = coastline , fill = "lightgreen") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "lightblue"))

basemap +
  geom_sf(data = unique_latlong, alpha = 0.8) +
  coord_sf(expand = FALSE)

ggsave("figures/sitemap.jpg", width = 6, height = 8)



##
# show data by ecoregion and measurement type
# number of trajectories
##

traj_dat <- nwa_data |>
  select(ecoregion,
         trajectory,
         has_bm, 
         has_sd,
         has_pc) |>
  pivot_longer(-c(ecoregion, trajectory)) |>
  filter(value>0) |>
  group_by(ecoregion, name) |>
  summarize(n_timeseries = n_distinct(trajectory)) |>
  mutate(name = case_when(
    name == "has_bm" ~ "Biomass",
    name == "has_sd" ~ "Stipe density",
    name == "has_pc" ~ "Percent cover"
  ))

ggplot(traj_dat,
       aes(y = ecoregion, x = name, 
           fill = n_timeseries,
           label = n_timeseries)) +
  geom_tile() +
  geom_label(fontface="bold") +
  labs(x = "", y = "") +
  scale_x_discrete(position = "top", 
                   expand = c(0,0))  +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_distiller(palette = "RdBu",
                       trans = "log10", 
                       guide = "none") +
  theme()

ggsave("figures/trajectories_by_type.jpg",
       width = 8, height = 4)

## table with # trajectories by ecoregion
nwa_data |>
  group_by(ecoregion) |>
  summarize(`Number of time series` = n_distinct(trajectory)) |>
  rename(Ecoregion = ecoregion) |>
  write_csv("tables/series_per_ecoregion.csv")

##
# show all timeseries collection times
##

nwa_data <- nwa_data |>
  mutate(trj = as.character(trajectory) |>
         forcats::fct_reorder(latitude, .desc = FALSE)) 

ggplot(nwa_data,
       aes(x = year, 
           y = trj |> as.numeric(),
           group = trj)) +
  geom_line(alpha = 0.7) +
  scale_y_continuous(guide = "none") +
  # scale_y_continuous(breaks = seq(0, 
  #                                 max(combined_clear$trj|>as.numeric()), 
  #                                 length.out=4),
  #                    labels = seq(min(combined_clear$latitude), 
  #                                 max(combined_clear$latitude), 
  #                                 length.out=4) |> round(2)) +
  labs(y = "\n\nN\n\n^\n|\nv\n\nS", x = "") +
  theme_bw() +
  facet_wrap(vars(focalUnit)) +
  theme(axis.title.y = element_text(angle = 0))


## 
# Show map of each envt predictor in June 2012
##
# oisst
oisst <- rast("data/rasters/oisst_monthly_mean_nma.nc")|>
  crop(aoi) 
basemap +
  geom_spatraster(data = oisst[[366]],
                  aes(fill = sst_366)) +
  scale_fill_distiller(palette = "RdBu", 
                       na.value = NA) +
  geom_sf(data = coastline , fill = "lightgreen") +
  coord_sf(expand = FALSE) +
  labs(title = "OISST average for June 2012",
       fill = "ºC")

ggsave("figures/oisst_june_2012.jpg",
        width = 6, height = 8)
# hadsst

hadsst <- rast("data/rasters/HadISST_sst_nwa.nc") |>
  crop(aoi)
basemap +
  geom_spatraster(data = hadsst[[762]],
                  aes(fill = sst_762)) +
  scale_fill_distiller(palette = "RdBu", 
                       na.value = NA) +
  geom_sf(data = coastline , fill = "lightgreen") +
  coord_sf(expand = FALSE) +
  labs(title = "HADSST average for June 2012",
       fill = "ºC")

ggsave("figures/hadsst_june_2012.jpg",
       width = 6, height = 8)

# no3
cmems <- rast("data/rasters/cmems_mod_glo_bgc_my_0.25deg_P1M-m_1750273483465.nc") |>
  crop(aoi) 
cmems <- cmems |>
  subset(time(cmems)==as.POSIXct("2012-06-01", tz = "UTC"))

basemap +
  geom_spatraster(data = cmems[[1]]) +
  scale_fill_distiller(palette = "GnBu", 
                       na.value = NA,
                       transform = "log10",
                       labels = comma) +
  geom_sf(data = coastline , fill = "lightgreen") +
  coord_sf(expand = FALSE) +
  labs(title = "CMEMS NO3 average for June 2012",
       fill = expression(paste("mmol/", m^3)))

ggsave("figures/no3_june_2012.jpg",
       width = 6, height = 8)

# kd490

turb_rast <- rast("data/rasters/erdMH1kd490mday.nc") |>
  crop(aoi) 
turb_rast <- turb_rast |>
  subset(time(turb_rast)==as.POSIXct("2012-06-16", tz = "UTC"))

basemap +
  geom_spatraster(data = turb_rast[[1]]) +
  scale_fill_distiller(palette = "YlOrBr", 
                       na.value = NA,
                       transform = "log10",
                       direction = 1) +
  geom_sf(data = coastline , fill = "lightgreen") +
  coord_sf(expand = FALSE) +
  labs(title = "June 2012 composite Kd490 average",
       fill = expression(m^-1))

ggsave("figures/turb_june_2012.jpg",
       width = 6, height = 8)

# swh
swh_rast <- rast("data/rasters/reanalysis-era5-single-levels-monthly-means_swh.nc") |>
  crop(aoi) 

swh_rast <- swh_rast |>
  subset(time(swh_rast)==as.POSIXct("2012-06-01", tz = "UTC"))

basemap +
  geom_spatraster(data = swh_rast[[1]]) +
  scale_fill_distiller(palette = "PuOr", 
                       na.value = NA,
                       direction = -1) +
  geom_sf(data = coastline , fill = "lightgreen") +
  coord_sf(expand = FALSE) +
  labs(title = "June 2012 Era5 SWH average",
       fill = "m")

ggsave("figures/swh_june_2012.jpg",
       width = 6, height = 8)


## 
# timeseries of envt parameters
##

## 
# show distance from sites to where data collected
# for different envt parameters
##