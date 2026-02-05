#' -----------------------------------------------
#' NWA Timeseries Modeling
#' looking at variation across space and spacetime 
#' with sdmTMB
#' -----------------------------------------------

library(sdmTMB)
library(sdmTMBextra)
library(sf)

####
## Load data and make into a spatially aware sf object
####
source("scripts/load_nwa_data.R")

# make sf, and also add UTM columns
# note will use UTM zone 19N; CRS = 32619.
# https://epsg.io/32619
nwa_dat_sf <- nwa_dat |> 
  st_as_sf(crs = 4326, 
           coords = c("longitude", "latitude"),
           remove = FALSE) |>
  add_utm_columns(units = "km")

# needed for some methods that don't like sf objects
nwa_dat <- nwa_dat|>
  add_utm_columns(units = "km")

###
## Create an AOI
## and get the coastline in that area
###

unique_latlong <- read_csv("data/unique_latlongs_time.csv") |>
  st_as_sf(crs = 4326,
           coords = c("longitude", "latitude"))

aoi <- unique_latlong |>
  st_bbox() + c(-3, -3, 3, 3) # some buffer

coastline <- rnaturalearth::ne_states(country = c("United States of America", "Canada"), 
                                      returnclass = "sf") |>
  st_crop(aoi) |>
  summarize()

# need it in km - using proj4 here
coastline_32619_km <- coastline |> 
  st_transform(crs = "+proj=utm +zone=19 +datum=WGS84 +units=km +no_defs")

ggplot() + 
  geom_sf(data= coastline|> st_transform(crs = 32619)) + 
  geom_sf(data = unique_latlong|> st_transform(crs = 32619))
###
## Make a mesh and 
## then make the coastline a barrier
###

# cutoff should be 1/3 of expected range over which the correlation
# drops to ~.13 in KM
basic_mesh <- make_mesh(nwa_dat_sf |> as_tibble(),
                        c("X", "Y"),
                        cutoff = 50)#, # minimum triangle edge length - distance between locations
#  max.edge = c(2, 80), # inner and outer max triangle lengths
# offset = c(1, 1)) # inner and outer border widths
plot(basic_mesh)



mesh_coastline <- sdmTMBextra::add_barrier_mesh(
  spde_obj = basic_mesh,
  barrier_sf = coastline |> st_transform(crs = 32619),
  plot = TRUE
)

plot(mesh_coastline)


####
## Fit a linear model
## where the coefficients and intercepts
## vary by space
####


mod_spatial <- sdmTMB(ln_focal_std_by_ecoregion ~ 
                        year_c + 
                        (1|study) +
                        (1+year_c | trajectory), 
                      spatial = "on",
                      spatial_varying = ~1 + year_c, 
                      mesh = mesh_coastline,
                      data = nwa_dat)

# predictions!
fitted_spatial <- predict(mod_spatial)

predict(mod_spatial) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  ggplot() +
  geom_sf(data = coastline) +
  geom_sf(aes(color = zeta_s_year_c + coef(mod_spatial)[2] )) +
  scale_color_viridis_c() +
  labs(color = "trend")


predict(mod_spatial) |>
  ggplot(aes(y = zeta_s_year_c + coef(mod_spatial)[2],
             x = latitude, color = ecoregion)) +
  geom_point()

predict(mod_spatial) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  ggplot() +
  geom_sf(data = coastline) +
  geom_sf(aes(color = (zeta_s_year_c + coef(mod_spatial)[2] )>0)) +
  labs(color = "direction")

