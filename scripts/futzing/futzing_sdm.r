library(sdmTMB)
library(sdmTMBextra)
library(mregions)
library(sf)

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

# make a 1x1 grid
ocean_grid <- st_make_grid(coastline, cellsize = 0.05) |>
  st_centroid() |>
  st_as_sf()

# crop it within 1km of the coastline
# create a coastline buffer within 1km of coast
coastline_buffer <- st_buffer(coastline, dist = 2000) |>
  st_difference(coastline)

plot(coastline)
plot(coastline_buffer, col = "red", lty = 0, add = TRUE)

coastline_grid <- st_filter(ocean_grid, coastline_buffer)

plot(coastline)
plot(coastline_buffer, col = "red", lty = 0, add = TRUE)
plot(coastline_grid, add = TRUE)

pred_grid <- st_coordinates(coastline_grid) |>
  as_tibble() |>
  rename(longitude = X, latitude = Y )

# make mesh
# first, create UTM coordinates in KM (remember for plotting!)
nwa_dat <- add_utm_columns(nwa_dat)

# cutoff should be 1/3 of expected range over which the correlation
# drops to ~.13 in KM
basic_mesh <- make_mesh(nwa_dat,
               c("X", "Y"),
#               fmesher_func = fmesher::fm_mesh_2d_inla,
               cutoff = 50)#, # minimum triangle edge length - distance between locations
             #  max.edge = c(2, 80), # inner and outer max triangle lengths
              # offset = c(1, 1)) # inner and outer border widths
plot(basic_mesh)

# coastline <- add_utm_columns(coastline)
# 
# mesh_coastline <- add_barrier_mesh(
#   spde_obj = basic_mesh,
#   barrier_sf = coastline
# )
# 
# plot(mesh_coastline)

# spatial variation in trend
# think about anisotropy
mod_spatial_0 <- sdmTMB(ln_focal_std_by_ecoregion ~ 
                          year_c + 
                       #   (1  + year_c|trajectory) + # do we need b/c space?
                        (1|study), 
                        spatial = "on",
                        #spatial_varying = ~1 + year_c,
                        mesh = basic_mesh,
                        data = nwa_dat)

mod_spatial <- sdmTMB(ln_focal_std_by_ecoregion ~ 
                          year_c + 
                          #   (1  + year_c|trajectory) + # do we need b/c space?
                          (1|study), 
                        spatial = "on",
                        spatial_varying = ~0 + year_c, # xy coord + study effect - don't need intercept?
                        mesh = basic_mesh,
                        data = nwa_dat)

predict(mod_spatial) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  ggplot() +
  geom_sf(data = coastline) +
  geom_sf(aes(color = zeta_s_year_c + coef(mod_spatial)[2] )) +
  scale_color_viridis_c() +
  labs(color = "trend")
#####

pred_dat <- replicate_df(nwa_dat, 
                         time_name = "year", 
                         time_values = 1:3)
# make a prediction grid that gets slopes 
# and then calculate the slope for each point
# OR calculate the estimate in two years and look at the difference




mod_spatiotemporal <- sdmTMB(ln_focal_std_by_ecoregion ~ 
                        year_c + 
                        #   (1  + year_c|trajectory) + # do we need b/c space?
                        (1|study), 
                      spatial = "on",
                      spatiotemporal = "rw",
                      time = "year",
                      spatial_varying = ~0 + year_c, # xy coord + study effect - don't need intercept?
                      time_varying = ~0 + year_c,
                      mesh = basic_mesh,
                      data = nwa_dat,
                     # silent = FALSE,
                      extra_time = c(1950, 1951, 1952, 1953, 1954, 1955,
                                     1956, 1957, 1958, 1959, 1960, 1961, 1962, 1963, 1964,
                                     1965, 1966, 1967, 1969, 1970, 1971, 1972, 1974, 1976))





mod_spatial2 <- sdmTMB(ln_focal_std_by_ecoregion ~ 
                        s(year_c) + 
                        (1  + year_c|trajectory) + 
                         (1|study), 
                      spatial = "on",
                      spatial_varying = ~1 + year_c,
                      mesh = mesh_coastline,
                      data = nwa_dat)

pred_grid2 <- tidyr::expand_grid(pred_grid,
                                 tibble(eco_collapsed = unique(nwa_dat$eco_collapsed)))
pred <- predict(mod_spatial,
                newdata = pred_grid2 |> mutate(year_c = 0),
                type = "response",
                re_form_iid = ~0)

pg <- pred_grid
pg$trend <- pred$zeta_s_year_c + coef(mod_spatial)[2]
ggplot(pg |> st_as_sf(coords = c("longitude", "latitude"),
                      crs = st_crs(coastline)),
       ) +
  geom_sf(data = coastline) +
    geom_sf(aes(fill = trend, color = trend)) +
  scale_fill_viridis_c() +
  scale_color_viridis_c()
#emmeans::emtrends(mod_spatial, ~eco_collapsed, "year_c") |>
#  plot()
# 
# mod_ecoregion <- sdmTMB(formula = ln_focal_std_by_ecoregion ~ 
#                            eco_collapsed*year_c + 
#                            (1  + year_c|trajectory) + (1|study), 
#                         # dispformula =~focalUnit,
#                         time_varying = ~ 0 + year_c,
#                         spatial_varying = ~1 + year_c,
#                          mesh = mesh_coastline,
#                          spatiotemporal = "ar1",
#                         spatial = "on",
#                         time = "year",
#                          data = nwa_dat)
# 

# mgcv
library(mgcv)

mod_spatiotemporal_gam <- gam(ln_focal_std_by_ecoregion ~ 
              s(longitude, latitude) + 
              s(year_c) + 
              ti(longitude, latitude,year_c, 
                 d=c(2,1)), 
           data=nwa_dat)

mod_spatial_gam <- gam(ln_focal_std_by_ecoregion ~ 
                 s(year_c, latitude, longitude) +
                   s(trajectory, bs = "re") +
                   s(study, bs = "re"),
                 data = nwa_dat)


mod_temporal_gam <- gam(ln_focal_std_by_ecoregion ~ 
                         s(year_c) +
                         s(latitude, longitude) +
                         (1  + year_c|trajectory) + 
                         (1|study), 
                       data = nwa_dat)


mod_spatiotemporal_gam <- gam(ln_focal_std_by_ecoregion ~ 
                         s(year_c) +
                         s(latitude, longitude) +
                         (1  + year_c|trajectory) + 
                         (1|study), 
                       data = nwa_dat)



# Fit GAM with spatiotemporal smooth
gam_fit <- gam(ln_focal_std_by_ecoregion ~ 
                 te(longitude, latitude, year, bs = c("tp", "tp", "tp")),
               data = nwa_dat,
               method = "REML",
               family = gaussian())  # Change to poisson() or nb() for counts

# This gives the local rate of change
derivs <- derivatives(gam_fit)

# Filter derivative to a specific year slice or time window
year_of_interest <- 2020
slice_deriv <- derivs %>%
  filter(abs(year - year_of_interest) < 1e-6)  # allow for floating point


# Create prediction grid
lon_seq <- seq(min(nwa_dat$longitude), max(nwa_dat$longitude), length.out = 100)
lat_seq <- seq(min(nwa_dat$latitude), max(nwa_dat$latitude), length.out = 100)
year_seq <- seq(1970, 2020, by = 2)

grid <- expand.grid(longitude = lon_seq, latitude = lat_seq, year = year_seq)
grid$pred_abund <- predict(gam_fit, newdata = grid, type = "response")



grid_change <- grid |>
  group_by(latitude, longitude) |>
  arrange(year) |>
  mutate(pred_change = pred_abund - lag(pred_abund),
         direction = sign(pred_change))

grid_change <- grid_change |>
  st_as_sf(coords = c("longitude", "latitude"),
           remove = FALSE,
           crs = 4326) |>
  st_intersection(coastline_buffer)

ggplot(grid_change |> filter(year %in% c(1975, 1985, 1995, 2005, 2015, 2020)),
       aes(x = longitude, y = latitude, 
           fill = pred_change, color = pred_change)) +
  geom_sf() +
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  facet_wrap(vars(year)) +
  theme_minimal()
  

ggplot(grid_change,
       aes(x = longitude, y = latitude, 
            color = as.factor(direction))) +
  geom_sf(size = 0.1) +
  facet_wrap(vars(year)) +
  theme_minimal()

