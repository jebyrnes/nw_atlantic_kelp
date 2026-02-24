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
source("scripts/sdmPrep.R") #also loads data and constants/helpers


####
## Fit a linear model
## where the coefficients and intercepts
## vary by space
####

mod_spatial <- sdmTMB(ln_focal_std_by_realm ~ 
                        year_c  +
                        (1|study) +
                        (1+year_c | trajectory), 
                      dispformula = ~focalUnit,
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

# randomized quantile residuals
dharma_plot(mod_spatial)


# coefs
tidy(mod_spatial,  conf.int = TRUE)
tidy(mod_spatial, effects = "ran_pars", conf.int = TRUE)
spatial_mod_re <- tidy(mod_spatial, effects = "ran_vals")

##
# Visualize model outputs
##

# get some fitted data and some derived quantities
fitted_spatial <- get_predicted_sdm_data(mod_spatial,
                                         lwr = 0.05, upr = 0.95) 

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
  mutate(year_c = 0,
         trajectory = nwa_dat$trajectory[1],
         study = nwa_dat$study[1],
         focalUnit = nwa_dat$focalUnit[1]
  )

spatial_interpolate <- get_predicted_sdm_data(mod_spatial,
                                              dat = pred_data,
                                         lwr = 0.05, upr = 0.95) |>
  st_as_sf(crs = st_crs(coastline_32619_km))

# map of means
ggplot(spatial_interpolate) + 
  geom_sf(data = coastline_32619_km, color = NA) +
  geom_sf(aes(color = combined_slope_mean), 
             size = 1) + 
  scale_color_distiller(palette = "RdYlBu", direction =2) +
  labs(color = "Rate of Change\nper Year")


ggsave("figures/loss_by_space_interpolated_linear_realm_std.jpg", 
       width = 10, height = 8)

# all three - use function from sdmPrep
slope_triptych(spatial_interpolate, limits = c(-0.07, 0.026))

ggsave("figures/loss_by_space_interpolated_triptych_realm_std.jpg", 
       width = 10, height = 8)
