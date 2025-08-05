library(sf)
library(mregions)

ecoregions_shp <- mr_shp(key = "Ecoregions:ecoregions") |> 
  st_make_valid() |>
  filter(lat > 30, long < -30, long > -80) |>
  dplyr::select(ecoregion, geometry) |>
  st_buffer(4*5500)

add_ecoregions_to_data <- function(dat, 
                                   latitude = "y",
                                   longitude = "x",
                                   eco = ecoregions_shp){
  dat |>
    st_as_sf(crs = 4326, 
             coords = c(longitude, latitude),
             remove = FALSE) |>
    st_join(eco) |>
    as_tibble() |>
    select(-geometry)
    
}
