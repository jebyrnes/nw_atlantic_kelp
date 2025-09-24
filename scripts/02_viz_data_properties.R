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
