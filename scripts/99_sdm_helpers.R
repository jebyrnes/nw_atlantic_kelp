#' --------------------------------
#' Working with sdmTMB requires some
#' extra things that are shared beween several scripts
#' so here we set them up ONCE so that multiple scripts
#' can use meshes, coastlines, maps, and relevant
#' helper functions
#' ----------------------------------
#
library(sdmTMB)
library(sdmTMBextra)
library(sf)
source("scripts/constants_and_helpers.R")
source("scripts/load_nwa_data.R")


###
# Useful spatial objects
###


# make sf, and also add UTM columns
# note will use UTM zone 19N; CRS = 32619.
# https://epsg.io/32619
unique_latlong <- nwa_dat |> 
  st_as_sf(crs = 4326, 
           coords = c("longitude", "latitude"),
           remove = FALSE) |>
  add_utm_columns(units = "km", 
                  ll_names = c("longitude", "latitude"),
                  ll_crs = 4326)

# needed for some methods that don't like sf objects
nwa_dat <- nwa_dat|>
  add_utm_columns(units = "km")

###
## Create an AOI for analysis
## and get the coastline in that area
###

aoi <- unique_latlong |>
  st_bbox() + c(-1, -1, 1, 1) # some buffer

coastline <- rnaturalearth::ne_states(country = c("United States of America", "Canada"), 
                                      returnclass = "sf") |>
  st_crop(aoi) |>
  summarize()

# need it in km - using proj4 here
coastline_32619_km <- coastline |> 
  st_transform(crs = "+proj=utm +zone=19 +datum=WGS84 +units=km +no_defs")

# Check 
# ggplot() + 
#   geom_sf(data= coastline|> st_transform(crs = 32619)) + 
#   geom_sf(data = unique_latlong|> st_transform(crs = 32619))



###
## Function to make a mesh with a  
## coastline barrier
###

make_kelp_mesh <- function(cutoff = 30, 
                           buffer = -2,
                           debug_plot_mesh = FALSE,
                           debug_show_land = FALSE,
                           ...) {

    basic_mesh <- make_mesh(nwa_dat,
                            xy_cols = c("X", "Y"),
                            cutoff = cutoff, ...) # minimum triangle edge length - distance between locations
  
    
    mesh_coastline <- sdmTMBextra::add_barrier_mesh(
      spde_obj = basic_mesh,
      barrier_sf = coastline_32619_km |> st_buffer(buffer)
    )
      
    if(debug_plot_mesh){
      plot(mesh_coastline)
    }      
    
    if(debug_show_land){
        
        # show mesh using code from workshop
        # https://github.com/sdmTMB/sdmTMB-teaching/blob/main/dfo-tesa-2023/11-barrier-models.Rmd
        mesh_df_water <- mesh_coastline$mesh_sf[mesh_coastline$normal_triangles, ]
        mesh_df_land <- mesh_coastline$mesh_sf[mesh_coastline$barrier_triangles, ]
        ggplot(coastline_32619_km) +
          geom_sf() +
          geom_sf(data = mesh_df_water, size = 2, colour = "blue") +
          geom_sf(data = mesh_df_land, size = 2, colour = "darkgreen")
    }
    
    return(mesh_coastline)
    
}


####
# function to create valid areas
# for an interpolation plot
# so far away interpolations are not 
# included
####

make_valid_areas <- function(centers = 60, buffer = 40, debug_plot = FALSE){
  
  km <- kmeans(st_coordinates(unique_latlong), centers = centers)

  valid_areas <- unique_latlong |>
    mutate(cluster = km$cluster) |>
    group_by(cluster) |>
    summarize() |>
    st_concave_hull(ratio = 1, allow_holes = TRUE)  |>
    st_transform(st_crs(coastline_32619_km)) |>
    st_buffer(buffer) |> #km from sample areas
    summarize() |>
    st_difference(coastline_32619_km)
  
  if(debug_plot){
    ggplot() +
      geom_sf(data = coastline_32619_km) +
      geom_sf(data = valid_areas , fill = "blue") 
  }

  return(valid_areas)

}

valid_areas <- make_valid_areas()


####
## Function to make a grid for predictions
## Can take a while, so make one, save it out, and load
####

make_prediction_grid <- function(
    coastline_buffer = 10,
    grid_size = 200,
    debug_plot_coastline_buffer = FALSE,
    debug_plot_grid = FALSE,
    debug_plot_overlap_grid = FALSE,
    debug_plot_final_grid = FALSE,
    ...
    ){

  # crop it within 1km of the coastline
  # create a coastline buffer within 10km of coast
  coastline_buffer <- st_buffer(coastline_32619_km, dist = coastline_buffer) |>
    st_difference(coastline_32619_km)
  
  if (debug_plot_coastline_buffer) {
    ggplot() +
      geom_sf(data = coastline_32619_km, size = 1.5) +
      geom_sf(data = coastline_buffer,
              color = "red",
              fill = NA)
  }
  
  # create grid over the bounding box of the polygon
  full_grid <- st_make_grid(coastline_buffer, n = grid_size) |>
    st_sf()
  
  if (debug_plot_grid) {
    ggplot(full_grid) + geom_sf() +
      geom_sf(data = coastline_buffer)
  }

  # get the overlap
  overlap <- st_intersection(full_grid, coastline_buffer)

  if(debug_plot_overlap_grid) ggplot(overlap) + geom_sf()
  
  # create new points at the centroid of the grid cells
  prediction_points <- st_centroid(overlap)
  names(prediction_points) <- "geometry"
  st_geometry(prediction_points) <- "geometry"

  prediction_points <- prediction_points |>
    bind_cols(st_coordinates(x = prediction_points)) |>
    mutate(cell_id = 1:n())  
  
  
  if (debug_plot_final_grid) {
    ggplot() +
      geom_sf(data = coastline_32619_km) +
      geom_sf(data = valid_areas , fill = "lightblue") +
      geom_sf(data = prediction_points, size = 0.01)
  }
  
  ## add ecoregions and eco_collapsed
  ecor <- mregions::mr_shp(key = "Ecoregions:ecoregions") |>
    st_make_valid() |>
    dplyr::select(ecoregion, geometry)
  
  prediction_points <- prediction_points |>
    st_join(ecor |> st_transform(crs = st_crs(prediction_points)), join = st_nearest_feature) |>
    mutate(
      ecoregion = ifelse(ecoregion == "Sunda Shelf/Java Sea", "Scotian Shelf", ecoregion),
      eco_collapsed =
        forcats::fct_collapse(
          ecoregion,
          "Gulf of St. Lawrence - Newfoundland" =
            c(
              "Gulf of St. Lawrence - Eastern Scotian Shelf",
              "Southern Grand Banks - South Newfoundland",
              "Northern Grand Banks - Southern Labrador"
            )
        )
    )
  
  return(prediction_points)

}

# generate and save a grid so we don't have the startup
# time when we load this script
# full_prediction_points <- make_prediction_grid()
# 
# ggplot() +
#   geom_sf(data = coastline_32619_km) +
#   geom_sf(data = valid_areas , fill = "lightblue") +
#   geom_sf(data = full_prediction_points|> st_intersection(valid_areas),
#           size = 0.01, aes(color = eco_collapsed)) +
#   scale_color_brewer(palette = "Dark2")
# 
# saveRDS(full_prediction_points, "data/full_coastline_prediction_grid.rds")
# saveRDS(full_prediction_points |> st_intersection(valid_areas), "data/valid_coastline_prediction_grid.rds")

full_prediction_points <- readRDS("data/full_coastline_prediction_grid.rds")
prediction_points <- readRDS("data/valid_coastline_prediction_grid.rds")


####
# Useful helper functions
####

prediction_check_density <- function(mod, seed = 31415, 
                                     trans = identity, 
                                     return_data = FALSE, ...){
  m <- simulate(mod, 
                nsim = 100,
                type = "mle-mvn",
                seed = seed,
                mle_mvn_samples = "multiple") |>
    as_tibble() |>
    bind_cols(mod$data) |>
    tidyr::pivot_longer(matches("^V\\d"))
  
  if(return_data) return(m)

  
  ggplot() +
    geom_density(data = m ,  
                 mapping = aes(x = trans(value), group = name),
                 alpha = 0.1,
                 fill = NA, color = "lightgrey")+
    geom_density(data =bind_cols(tibble(x = mod$tmb_data$y_i),
                                 mod$data),
                 mapping = aes(x = trans(x)),
                 alpha = 0.1,
                 fill = NA, color = "darkblue")
}

dharma_plot <- function(mod, nsim = 300, 
                        type = "mle-mvn", 
                        seed = 2718282){
  set.seed(seed)
  simulate(mod, nsim = nsim, type = type) |>
    dharma_residuals(mod)
}

# Get predicted values from a model including spatial
# predictions and error and a yi variable for whatever the
# model predicted if this is on the model data

get_predicted_sdm_data <- function(mod, dat = NULL, nsim = 1e3,
                                   slopevar = "year_c",
                                   lwr = 0.025, upr = 0.975, ...){
  if(is.null(dat)){
    dat <- mod$data
    dat$yi <- mod$tmb_data$y_i
  }
  
  # for combined variable    
  zeta_s <- predict(mod, newdata = dat,
                    nsim = nsim, 
                    sims_var = "zeta_s", ...)
  
  sims <- spread_sims(mod, nsim = nsim)
  
  combined <- sims[[slopevar]] + t(zeta_s)
  dat$combined_slope_median <- apply(combined, 2, median)
  dat$combined_slope_mean <- apply(combined, 2, mean)
  dat$combined_slope_sd <- apply(combined, 2, sd)
  dat$combined_slope_lwr <- apply(combined, 2, quantile, probs = lwr)
  dat$combined_slope_upr <- apply(combined, 2, quantile, probs = upr)
  
  predict(mod, se=TRUE, newdata = dat, ...) |>
    mutate(slope_fe = coef(mod)[2],
           slope_fe_se = tidy(mod)$std.error[2])
}

plot_fit_sdm_model <- function(mod, dat = NULL){
  if(is.null(dat)) dat <- get_predicted_sdm_data(mod)
  
  ggplot(dat,
         aes(x = year,
             y = exp(est)-0.01,
             color = eco_collapsed,
             group = trajectory,
         )) +
    geom_line(linewidth = 1.1, alpha = 0.9) +
    labs(y = "Standardized Kelp Abundance", x = "",
         color = "") +
    guides(color=guide_legend(nrow=2,
                              byrow=TRUE))+
    theme_bw(base_size = 14)+
    theme(legend.position = "bottom",
          legend.box="vertical") +
    scale_color_brewer(palette = "Dark2") +
    geom_point(aes(y = yi), alpha = 0.1) +
    facet_wrap(vars(eco_collapsed)) +
    guides(color = "none")
}

###
# for the plots
###

slope_triptych <- function(dat,
                           color_lab = "Porportion Change\nPer Year",
                           title = "",
                           mean_title = "Mean Slope",
                           lwr_title = "Lower CI",
                           upr_title = "Upper CI",
                           size = 0.8,
                           limits = c(-0.07, 0.27),
                           n.breaks = 5,
                           ...){
  
  
  # show study slopes with FE + spatial variation
  mean_map <- ggplot(dat) +
    geom_sf(data = coastline) +
    geom_sf(aes(color = combined_slope_mean), size = size, ...) +
    scale_color_slope_b(limits = limits, n.breaks = n.breaks) +
    labs(color = color_lab, 
         title = title,
         subtitle = mean_title)
  
  
  lwr_map <- ggplot(dat) +
    geom_sf(data = coastline) +
    geom_sf(aes(color = combined_slope_lwr), size = size, ...) +
    scale_color_slope_b(limits = limits, n.breaks = n.breaks) +
    labs(color = color_lab, 
         subtitle = lwr_title)
  
  upr_map <- ggplot(dat) +
    geom_sf(data = coastline) +
    geom_sf(aes(color = combined_slope_upr), size = size, ...) +
    scale_color_slope_b(limits = limits, n.breaks = n.breaks) +
    labs(color = color_lab, 
         subtitle = upr_title) 
  
  
  mean_map /( lwr_map + upr_map)+
    plot_layout(guides = 'collect',
                width = c(1,2))&theme_light(base_size = 8)
  
}
