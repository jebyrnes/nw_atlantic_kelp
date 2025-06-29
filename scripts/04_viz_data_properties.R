#' -------------------------------------------
#' Viz data properties of kelptime
#' and envt data
#' -------------------------------------------

library(readr)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(mregions)
library(wesanderson)
library(sf)
library(terra)
library(tidyterra)
library(patchwork)

pal <- wes_palette("Zissou1", 4, type = "continuous")

theme_set(theme_bw(base_size = 17))

#reads in nwa data
source("scripts/load_nwa_data.R")

ecoregions_shp <- mr_shp(key = "Ecoregions:ecoregions") |> 
  st_make_valid() |>
  dplyr::select(ecoregion, geometry) |>
  filter(ecoregion %in% levels(nwa_dat$ecoregion))


old_data <- read_csv("data/nwa_krumhansl_2016.csv")

old_unique_latlong <- old_data |>
  rename_all(tolower) |>
  group_by(latitude, longitude) |>
  slice(1L) |>
  ungroup()|>
  st_as_sf(crs = 4326,
           coords = c("longitude", "latitude"))

unique_latlong <- read_csv("data/unique_latlongs_time.csv") |>
  st_as_sf(crs = 4326,
           coords = c("longitude", "latitude"))

aoi <- unique_latlong |>
  st_bbox() + c(-3, -3, 3, 3) # some buffer

##
# Map of data
##

coastline <- ne_states(country = c("United States of America", "Canada"), 
                       returnclass = "sf") |>
  st_crop(aoi)

basemap <- ggplot() +
  geom_sf(data = coastline , fill = "seagreen3") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "lightblue")) +
  geom_sf(data = ecoregions_shp, fill = NA)

basemap +
  geom_sf(data = unique_latlong, alpha = 0.8, size = 2.5) +
  coord_sf(expand = FALSE, 
           xlim = c(aoi[1], aoi[3]),
           ylim = c(aoi[2], aoi[4])) 
  

ggsave("figures/sitemap.jpg", width = 6, height = 8)

basemap +
  geom_sf(data = old_unique_latlong, alpha = 0.8, size = 2.5) +
  coord_sf(expand = FALSE, 
           xlim = c(aoi[1], aoi[3]),
           ylim = c(aoi[2], aoi[4])) 

ggsave("figures/sitemap_krumhansl_2016.jpg", width = 6, height = 8)

##
# show data by ecoregion and measurement type
# number of trajectories
##

traj_dat <- nwa_dat |>
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

ggplot(traj_dat |> filter(!is.na(ecoregion)),
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
nwa_dat |>
  group_by(ecoregion) |>
  summarize(`Number of time series` = n_distinct(trajectory)) |>
  rename(Ecoregion = ecoregion) |>
  write_csv("tables/series_per_ecoregion.csv")

##
# show all timeseries collection times
##

nwa_dat <- nwa_dat |>
  mutate(trj = as.character(trajectory) |>
         forcats::fct_reorder(latitude, .desc = FALSE))

ggplot(nwa_dat |> filter(!is.na(ecoregion)),
       aes(x = year, 
           y = latitude,
           group = trj,
           color = ecoregion )) +
  geom_line(alpha = 0.7) +
  # scale_y_continuous(breaks = seq(0, 
  #                                 max(combined_clear$trj|>as.numeric()), 
  #                                 length.out=4),
  #                    labels = seq(min(combined_clear$latitude), 
  #                                 max(combined_clear$latitude), 
  #                                 length.out=4) |> round(2)) +
  labs(x = "", y = "", color = "") +
  scale_y_continuous(labels = ~ paste0(.x, "°")) +
  facet_wrap(vars(focalUnit))   +
  theme(legend.position = "bottom",
        legend.box="vertical")+ 
  theme_bw(base_size = 18) +
  theme()

ggsave("figures/length_by_focalUnit.jpg", width = 12, height = 4)


ggplot(nwa_dat,
       aes(x = year, 
           y = latitude,
           group = trj,
           color = ecoregion)) +
  geom_line(alpha = 0.7) +
  labs(x = "", y = "", color = "") +
  scale_y_continuous(labels = ~ paste0(.x, "°"))  +
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  theme_bw(base_size = 18)+
  theme(legend.position = "bottom",
        legend.box="vertical")


ggsave("figures/timeseries_with_length_latitude.jpg", width = 6, height = 4)

# whole shebang
ggplot(nwa_dat,
       aes(x = year, 
           y = trj |> as.numeric(),
           group = trj)) +
  geom_line(alpha = 0.7) +
  labs(x = "", y = "", fill = "") +
  theme_bw() +
  scale_y_continuous(guide = "none") +
  theme(legend.position = "bottom",
        legend.box="vertical") 

ggsave("figures/timeseries_with_length_trj.jpg", width = 8, height = 4)

##
# Histograms of properties
##

old_hist <- ggplot(data = old_data,
       aes(x = year)) + 
  geom_histogram(bins = 30) +
  xlim(c(1941, 2030)) +
  ylim(c(0,300)) +
  labs(x="", y = "# sites",
       subtitle = "Krumhansl et al. 2016")


new_hist <- ggplot(data = nwa_dat,
       aes(x = year)) + 
  geom_histogram(bins = 40)+
  xlim(c(1941, 2030)) +
  ylim(c(0,300)) +
  labs(x="", y = "# sites", 
       subtitle = "Current Dataset") 

old_hist+new_hist +
  plot_layout(axis_titles = "collect")

ggsave("figures/year_hist.jpg",
       width = 8, height = 4)

## Length
old_length <- old_data |>
  group_by(Site) |>
  summarize(study_length = max(year) - min(year)) |>
  ggplot(aes(x = study_length)) +
  geom_histogram(bins = 40) +
  xlim(c(-1, 80)) + ylim(c(-1, 55)) +
  labs(x = "time series duration (years)", y = "# sites",
       subtitle = "Krumhansl et al. 2016")

new_length <- nwa_dat |>
  group_by(site) |>
  summarize(study_length = max(year) - min(year)) |>
  ggplot(aes(x = study_length)) +
  geom_histogram(bins = 40) +
  xlim(c(-1, 80)) + ylim(c(-1, 55)) +
  labs(x = "time series duration (years)", y = "# sites",
       subtitle = "Current dataset")

old_length + new_length +
  plot_layout(axis_titles = "collect")


ggsave("figures/duration_hist.jpg",
       width = 8, height = 4)

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
  labs(title = "June 2012 average OISST",
       fill = "ºC")+
  coord_sf(expand = FALSE, 
           xlim = c(aoi[1], aoi[3]),
           ylim = c(aoi[2], aoi[4])) 

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
  labs(title = "June 2012 average HADSST",
       fill = "ºC")+
  coord_sf(expand = FALSE, 
           xlim = c(aoi[1], aoi[3]),
           ylim = c(aoi[2], aoi[4])) 

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
                       transform = "log10") +
  geom_sf(data = coastline , fill = "lightgreen") +
  labs(title = "June 2012 average CMEMS NO3",
       fill = expression(paste("mmol/", m^3)))+
  coord_sf(expand = FALSE, 
           xlim = c(aoi[1], aoi[3]),
           ylim = c(aoi[2]+2, aoi[4])) 

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
  labs(title = "June 2012 composite Kd490 average",
       fill = expression(m^-1))+
  coord_sf(expand = FALSE, 
           xlim = c(aoi[1], aoi[3]),
           ylim = c(aoi[2], aoi[4])) 

ggsave("figures/turb_june_2012.jpg",
       width = 6, height = 8)

# swh
swh_rast <- rast("data/rasters/reanalysis-era5-single-levels-monthly-means_swh.nc")|>
  crop(aoi) 


swh_rast <- swh_rast |>
  subset(time(swh_rast)==as.POSIXct("2012-06-01", tz = "UTC"))

basemap +
  geom_spatraster(data = swh_rast[[1]]) +
  scale_fill_distiller(palette = "PuOr", 
                       na.value = NA,
                       direction = -1) +
  geom_sf(data = coastline , fill = "lightgreen") +
  labs(title = "June 2012 Era5 SWH average",
       fill = "m")+
  coord_sf(expand = FALSE, 
           xlim = c(aoi[1], aoi[3]),
           ylim = c(aoi[2], aoi[4])) 

ggsave("figures/swh_june_2012.jpg",
       width = 6, height = 8)


## 
# timeseries of envt parameters
##

env_dat <- read_csv("data/merged_envt_data_all.csv") |>
  select(year, 
         oisst_summer, 
         hadsst_summer,
         no3_spring,
         swh_fall,
         kd490_summer
         ) |>
  pivot_longer(-year) |>
  filter(!is.na(value)) |>
  group_by(name) |>
  summarize(min_year = min(year),
            max_year = max(year)) |>
  mutate(proxy = 
           case_when(
             name == "hadsst_summer" ~"Temperature: HADSST",
             name == "oisst_summer" ~"Temperature: OISST",
             name == "kd490_summer" ~"Turbidity: kd490",
             name == "no3_spring" ~"Nitrate: CMEM Model",
             name == "swh_fall" ~"Waves: Era 5 Model"
             
           )) |>
  pivot_longer(c(min_year, max_year), 
               names_to = "type", 
               values_to = "year")

ggplot(env_dat,
       aes(x = year, y = proxy)) + 
  geom_line(size = 2) +
  labs(x = "", y = "",
       subtitle = "Temporal coverage of proxies")

ggsave("figures/proxy_timeline.jpg", width = 6, height = 4)

## 
# show distance from sites to where data collected
# for different envt parameters
##