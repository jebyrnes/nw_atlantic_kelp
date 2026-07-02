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
mesh_nocoast <- make_mesh(nwa_dat,
                          xy_cols = c("X", "Y"),
                          cutoff = 10)
valid_areas <- make_valid_areas()
# prediction_points <- full_prediction_points |> 
#   st_intersection(valid_areas)

####
## Load the model
####
mod_spatial <- readRDS("models/mod_spatial.rds")


##
# Visualize model outputs
##

# get some fitted data and some derived quantities
mod_coefs <- tidy(mod_spatial)

fitted_spatial <- get_predicted_sdm_data(mod_spatial,
                                         lwr = 0.05, upr = 0.95) |>
  st_as_sf(coords = c("X", "Y"), crs = st_crs(coastline_32619_km))

slopes_only <- fitted_spatial |>
  group_by(trajectory) |>
  slice(1L) |>
  ungroup() |>
  select(trajectory, latitude, eco_collapsed, combined_slope_median:combined_slope_upr) |>
  arrange(combined_slope_mean) |>
  mutate(samp = 1:n())

# plot the slopes
ggplot(slopes_only,
       aes(x = latitude, y = combined_slope_mean,
           ymin = combined_slope_lwr, ymax = combined_slope_upr)) +
  geom_point() +
  geom_linerange(alpha = 0.2) +
  geom_hline(yintercept = 0, lty = 2, color = "red") +
  coord_flip() +
  facet_wrap(~eco_collapsed)

# Show the curves resulting from this model
plot_fit_sdm_model(mod_spatial, fitted_spatial)

ggsave("figures/spatial_kelp_change_fitted_curves.jpg", 
       width = 8, height = 5)

# show study slopes with FE + spatial variation
col_lims <- rep(max(abs(c(fitted_spatial$combined_slope_lwr, fitted_spatial$combined_slope_upr))),
                2) * c(-1,1)

slope_map <- ggplot(fitted_spatial) +
  geom_sf(data = coastline) +
  geom_sf(aes(color = combined_slope_mean), size = 2) +
  labs(color = "Rate of Change\nper year", 
       title = "For Included Studies",
       subtitle = "Mean")

slope_map +
  scale_color_slope_b(limits =c(-0.075, 0.075), n.breaks = 7) 
ggsave("figures/spatial_kelp_change_fitted_mean.jpg", width = 7, height = 7)


mean_map <- slope_map + scale_color_slope_b(limits =col_lims, n.breaks = 7) +


lwr_map <- ggplot(fitted_spatial) +
  geom_sf(data = coastline) +
  geom_sf(aes(color = combined_slope_lwr), size = 2) +
  scale_color_slope_b(limits = col_lims, n.breaks = 7) +
  labs(color = "Rate of Change\nper year", 
       subtitle = "Lower 90% CI")

upr_map <- ggplot(fitted_spatial) +
  geom_sf(data = coastline) +
  geom_sf(aes(color = combined_slope_upr), size = 2) +
  scale_color_slope_b(limits = col_lims, n.breaks = 7) +
  labs(color = "Rate of Change\nper year", 
       subtitle = "Upper 90% CI") 


mean_map /( lwr_map + upr_map)+
  plot_layout(guides = 'collect',
              width = c(1,2))&theme_light(base_size = 8)


ggsave("figures/spatial_kelp_change_fitted_triptych.jpg", width = 8, height = 5)

ggplot(fitted_spatial |>
         group_by(site) |>
         filter(sample_year == max(sample_year)),
       aes(x = latitude, y = combined_slope_mean,
           ymin = combined_slope_lwr, ymax = combined_slope_upr)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, color = "red", lty = 2) +
  facet_wrap(vars(eco_collapsed)) +
  coord_flip() 

##
# Interpolated Slopes
##

pred_data <- prediction_points |>  
  as_tibble() |> 
  mutate(year_c = 16,
         focalUnit = nwa_dat$focalUnit[1]
  )

spatial_interpolate <- get_predicted_sdm_data(mod_spatial,
                                              dat = pred_data,
                                         lwr = 0.05, upr = 0.95,
                                         re_form_iid = NA) |>
  st_as_sf(crs = st_crs(coastline_32619_km)) |>
  mutate(in_ci = ifelse(sign(combined_slope_mean) == sign(combined_slope_lwr)&
                          sign(combined_slope_mean) == sign(combined_slope_upr),
                        TRUE, FALSE))

# map of means
ggplot(spatial_interpolate) + 
  geom_sf(data = coastline_32619_km, color = NA, fill = "darkgrey") +
  geom_sf(aes(color = combined_slope_mean), 
             size = 1, shape = 19) + 
  scale_color_distiller(palette = "RdYlBu", direction =2,
                        limits = c(-0.075, 0.075),
                        breaks = seq(-0.06, 0.06, by = 0.02)) +
  #scale_color_slope_b(limits = c(-0.075, 0.075), n.breaks=7) +
  labs(color = "Porportion Change\nper Year") +
  scale_alpha_manual(values = c(0.2, 1)) +
  guides(alpha = "none") +
  theme_light(base_size = 16)


ggsave("figures/spatial_kelp_change_interpolated.jpg", 
       width = 8, height = 6)

# all three - use function from sdmPrep
slope_triptych(spatial_interpolate, 
               limits = c(-0.12, 0.12),
               n.breaks = 7)

ggsave("figures/spatial_kelp_change_interpolated_triptych_std.jpg", 
       width = 7, height = 5)


# curves
spatial_slope <- predict(mod_spatial, 
                         type = "response",
                         re_form_iid = NA,
                         newdata = prediction_points |>  
                           as_tibble() |> 
                           tidyr::crossing(
                             tibble(year_c = seq(min(nwa_dat$year_c), max(nwa_dat$year_c), length.out =100),
                                    year = seq(min(nwa_dat$year), max(nwa_dat$year), length.out =100),
                                    focalUnit = rep(nwa_dat$focalUnit[3],100)
                             )))

ggplot(nwa_dat, aes(x = year, color = Y,
                    y = focal_std_by_all)) +
  #geom_point(alpha = 0.5) +
  geom_line(data = spatial_slope,
            aes(y = est, 
                color = Y,
                group = Y),
            alpha = 0.2)  +
  scale_color_viridis_c() +
  labs(y="", x = "Standardized Kelp Abundance", 
       color = "UTM Northing")

ggsave("figures/spatial_kelp_change_interpolated_curves.jpg")



##
# emtrends(mod_spatial_eco, ~eco_collapsed, "year_c") |> plot()
# emtrends(mod_spatial_eco, ~eco_collapsed, "year_c") |> 
#   contrast(method = "pairwise", adjust= "none") |> plot()
