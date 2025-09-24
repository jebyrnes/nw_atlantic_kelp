library(sdmTMB)
library(sdmTMBextra)
library(mregions)
library(sf)

source("scripts/load_nwa_data.R")

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
basic_mesh <- make_mesh(nwa_dat,
               c("longitude", "latitude"),
               n_knots = 50)

mesh_coastline <- add_barrier_mesh(
  basic_mesh,
  coastline
)

plot(mesh_coastline)


mod_spatial <- sdmTMB(ln_focal_std_by_ecoregion ~ 
                          year_c + 
                          (1  + year_c|trajectory) + (1|study), 
                        spatial = "on",
                        spatial_varying = ~1 + year_c,
                        mesh = mesh_coastline,
                        data = nwa_dat)


mod_spatial <- sdmTMB(ln_focal_std_by_ecoregion ~ 
                        year_c + 
                        (1  + year_c|trajectory) + (1|study), 
                      spatiotemporal = "on",
                      spatial_varying = ~1 + year_c,
                      mesh = mesh_coastline,
                      data = nwa_dat)

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