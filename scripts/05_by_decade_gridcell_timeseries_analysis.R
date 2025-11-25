library(mregions)
library(sf)
library(glmmTMB)
library(ggplot2)
library(performance)
library(modelbased)
library(car)

source("scripts/load_nwa_data.R")
nwa_dat_sf <- nwa_dat |> 
  st_as_sf(crs = 4326, 
           coords = c("longitude", "latitude"),
           remove = FALSE)

# coastline 

ecoregions_shp <- mr_shp(key = "Ecoregions:ecoregions") |> 
  st_make_valid() |>
  dplyr::select(ecoregion, geometry) |>
  filter(ecoregion %in% levels(nwa_dat$ecoregion))

unique_latlong <- read_csv("data/unique_latlongs_time.csv") |>
  st_as_sf(crs = 4326,
           coords = c("longitude", "latitude"))

aoi <- unique_latlong |>
  st_bbox() + c(-3, -3, 3, 3) # some buffer

coastline <- rnaturalearth::ne_states(country = c("United States of America", "Canada"), 
                                      returnclass = "sf") |>
  st_crop(aoi) |>
  summarize()

# make a grid
ocean_grid <- st_make_grid(coastline, cellsize = 2) |>
  st_as_sf() |>
  mutate(gridcell = 1:n())

# crop it within 1km of the coastline
# create a coastline buffer within 1km of coast
coastline_buffer <- st_buffer(coastline, dist = 2000) |>
  st_difference(coastline)

plot(coastline)
plot(coastline_buffer, col = "red",  add = TRUE)

coastline_grid <- st_filter(ocean_grid, coastline_buffer)

st_intersection(coastline_grid, coastline_buffer) |> plot()

pred_grid <- st_coordinates(coastline_grid |> st_centroid()) |>
  as_tibble() |>
  rename(longitude = X, latitude = Y )

##
# determine what gridcell observations fall into
##
nwa_dat_gridded <- st_join(nwa_dat_sf, ocean_grid, join = st_within)
nwa_dat_gridded <- nwa_dat_gridded |>
  mutate(gridcell = as.factor(gridcell))


##
# see our gridcell sample sizes
##
nwa_grid_stats <- nwa_dat_gridded |>
  group_by(gridcell) |>
  summarize(n_sites = n_distinct(trajectory)) 

nwa_grid_stats |>
  pull(n_sites) |> hist(breaks = 100)

good_cells <- nwa_grid_stats |> 
  filter(n_sites>2) |>
  pull(gridcell)|>
  as.character() |> unique()

nwa_dat_gridded_filtered <- nwa_dat_gridded |>
  filter(as.character(gridcell) %in% 
           good_cells)

ggplot() +
  geom_sf(data = coastline) +
  geom_sf(data = nwa_dat_gridded_filtered) +
  geom_sf(data = coastline_grid, fill = NA)

## 
# fit a model where slope varies by gridcell
# and compare to one where slope varies by ecoregion
##


mod_gridcell <- glmmTMB(ln_focal_std_by_ecoregion ~ 
                          gridcell*year_c +
                           (1 + year_c |trajectory) + (1|study), 
                         dispformula =~focalUnit,
                         data = nwa_dat_gridded_filtered)


mod_ecoregion <- glmmTMB(ln_focal_std_by_ecoregion ~ 
                           eco_collapsed*year_c + 
                           (1 + year_c |trajectory) + (1|study), 
                         dispformula =~focalUnit,
                         data = nwa_dat_gridded)


mod_gridcell_decade <- glmmTMB(ln_focal_std_by_ecoregion ~ 
                          decade*gridcell*year_c +
                          (1 + year_c |trajectory) + (1|study), 
                        dispformula =~focalUnit,
                        data = nwa_dat_gridded_filtered)

estimate_slopes(mod_gridcell_decade,
                trend = "year_c",
                  by = c("decade", "gridcell"),
                backend = "emmeans",
                include_random = FALSE) |> 
  plot() +
  coord_flip()
