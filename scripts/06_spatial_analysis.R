#' -----------------------------------------------
#' NWA Timeseries Modeling
#' looking at variation across space 
#' with sdmTMB
#' -----------------------------------------------

library(sdmTMB)
library(sdmTMBextra)
library(sf)
library(ggplot2)
library(patchwork)
suppressWarnings({
source("scripts/99_sdm_helpers.R") #also loads data and constants/helpers
})

###
## Make objects we will use
##
####
mesh_coastline <- make_kelp_mesh(cutoff = 10)
valid_areas <- make_valid_areas()
# prediction_points <- full_prediction_points |> 
#   st_intersection(valid_areas)

ggplot() +
  inlabru::gg(mesh_coastline$mesh,
              edge.color = "black",
              ext.color = "darkgrey"
  ) +
  geom_point(data = nwa_dat, aes(x = X, y = Y)) +
  theme_void()

ggsave("figures/mesh_10.jpg", width = 7, height = 5)
####
## Fit a linear model
## where the coefficients and intercepts
## vary by space
####

mod_spatial <- sdmTMB(focal_std_by_all ~ 
                        year_c  +
                        (1 | study) +
                        (1  + year_c | trajectory), 
                      dispformula = ~focalUnit,
                      family = tweedie(link = "log"),
                      spatial = "on",
                      spatial_varying = ~0 + year_c, 
                      mesh = mesh_coastline,
                      data = nwa_dat)
##
# model checks
##
sanity(mod_spatial)

# check predictions
prediction_check_density(mod_spatial)
prediction_check_density(mod_spatial, trans = \(x) log(x+0.01))

# randomized quantile residuals
dharma_plot(mod_spatial)


# coefs
tidy(mod_spatial,  conf.int = TRUE)
tidy(mod_spatial, effects = "ran_pars", conf.int = TRUE)
spatial_mod_re <- tidy(mod_spatial, effects = "ran_vals")

perc_change_per_year <- ((exp(coef(mod_spatial)[2]) -1 )*100) |> round(2)
perc_change_per_year

##
# Visualize model outputs
##

# get some fitted data and some derived quantities
fitted_spatial <- get_predicted_sdm_data(mod_spatial,
                                         lwr = 0.05, upr = 0.95) |>
  st_as_sf(coords = c("X", "Y"), crs = st_crs(coastline_32619_km))

# Show the curves resulting from this model
plot_fit_sdm_model(mod_spatial, fitted_spatial)


# show study slopes with FE + spatial variation
mean_map <- ggplot(fitted_spatial) +
  geom_sf(data = coastline) +
  geom_sf(aes(color = combined_slope_mean), size = 0.8) +
  scale_color_slope_b(limits = c(-0.07, 0.027)) +
  labs(color = "Rate of Change\nper year", 
       title = "For Included Studies",
       subtitle = "Mean")


lwr_map <- ggplot(fitted_spatial) +
  geom_sf(data = coastline) +
  geom_sf(aes(color = combined_slope_lwr), size = 0.8) +
  scale_color_slope_b(limits = c(-0.07, 0.027)) +
  labs(color = "Rate of Change\nper year", 
       subtitle = "Lower 95% CI")

upr_map <- ggplot(fitted_spatial) +
  geom_sf(data = coastline) +
  geom_sf(aes(color = combined_slope_upr), size = 0.8) +
  scale_color_slope_b(limits = c(-0.07, 0.027)) +
  labs(color = "Rate of Change\nper year", 
       subtitle = "Upper 95% CI") 


mean_map /( lwr_map + upr_map)+
  plot_layout(guides = 'collect',
              width = c(1,2))&theme_light(base_size = 8)

# range -0.07 to 0.03 limits = c(-0.07, 0.03)

ggsave("figures/spatial_kelp_change_fitted.jpg", width = 10, height = 7)

##
# Interpolated Slopes
##

pred_data <- prediction_points |>  
  as_tibble() |> 
  mutate(year_c = 16,
         trajectory = nwa_dat$trajectory[1],
         study = nwa_dat$study[1],
         focalUnit = nwa_dat$focalUnit[1]
  )

spatial_interpolate <- get_predicted_sdm_data(mod_spatial,
                                              dat = pred_data,
                                         lwr = 0.05, upr = 0.95) |>
  st_as_sf(crs = st_crs(coastline_32619_km)) |>
  mutate(in_ci = ifelse(sign(combined_slope_mean) == sign(combined_slope_lwr)&
                          sign(combined_slope_mean) == sign(combined_slope_upr),
                        TRUE, FALSE))

# map of means
ggplot(spatial_interpolate) + 
  geom_sf(data = coastline_32619_km, color = NA, fill = "darkgrey") +
  geom_sf(aes(color = combined_slope_mean, alpha = in_ci), 
             size = 2) + 
  scale_color_distiller(palette = "RdYlBu", direction =2,
                        limits = c(-0.07, 0.07)) +
  labs(color = "Rate of Change\nper Year") +
  scale_alpha_manual(values = c(0.2, 1)) +
  guides(alpha = "none") +
  theme_light(base_size = 16)



ggsave("figures/spatial_kelp_change_interpolated.jpg", 
       width = 8, height = 6)

# all three - use function from sdmPrep
slope_triptych(spatial_interpolate, 
               limits = c(-0.07, 0.07),
               n.breaks = 7)

ggsave("figures/spatial_kelp_change_interpolated_triptych_std.jpg", 
       width = 10, height = 8)


# curves
spatial_slope <- predict(mod_spatial, 
                         type = "response",
                         newdata = prediction_points |>  
                           as_tibble() |> 
                           tidyr::crossing(
                             tibble(year_c = seq(min(nwa_dat$year_c), max(nwa_dat$year_c), length.out =100),
                                    year = seq(min(nwa_dat$year), max(nwa_dat$year), length.out =100),
                                    trajectory = rep(nwa_dat$trajectory[1],100),
                                    study = rep(nwa_dat$study[1],100),
                                    focalUnit = rep(nwa_dat$focalUnit[1],100)
                             )))

ggplot(nwa_dat, aes(x = year, color = Y,
                    y = focal_std_by_all)) +
  geom_point(alpha = 0.5) +
  geom_line(data = spatial_slope,
            aes(y = est, 
                color = Y,
                group = Y),
            alpha = 0.2)  +
  scale_color_viridis_c() +
  labs(y="", x = "Standardized Kelp Abundance", 
       color = "UTM Northing")

ggsave("figures/spatial_kelp_change_interpolated_curves.jpg", 
       width = 10, height = 8)
