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


ggplot() +
  geom_sf(data = coastline , fill = "lightgreen") +
  theme_bw() +
  coord_sf(expand = FALSE) +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "lightblue")) +
  geom_sf(data = unique_latlong, alpha = 0.8)

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

## 
# timeseries of envt parameters
##

## 
# show distance from sites to where data collected
# for different envt parameters
##