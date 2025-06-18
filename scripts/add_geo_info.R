library(mregions) #devtools::install_github("ropensci/mregions")
library(sf)
library(dplyr)

ecoregions_shp <- mr_shp(key = "Ecoregions:ecoregions") |> 
  st_make_valid() |>
  dplyr::select(ecoregion, geometry)

provinces_shp <- mr_shp(key = "Ecoregions:provinces") |> 
  st_make_valid() |>
  dplyr::select(province, geometry)

realm_shp <- mr_shp(key = "Ecoregions:realm") |> 
  st_make_valid() |>
  dplyr::select(realm, geometry)

add_geo_info <- function(dat){
  dat |>
    st_as_sf(coords = c("longitude", "latitude"),
             crs = st_crs(ecoregions_shp),
             remove = FALSE)|>
    st_join(ecoregions_shp, 
            join = st_nearest_feature)|>
    st_join(provinces_shp, 
            join = st_nearest_feature)|>
    st_join(realm_shp, 
            join = st_nearest_feature)
}
