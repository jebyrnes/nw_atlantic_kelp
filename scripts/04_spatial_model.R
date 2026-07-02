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

ggplot() +
  geom_sf(data = coastline_32619_km) +
  inlabru::gg(mesh_coastline$mesh,
              edge.color = "black",
              ext.color = "darkgrey"
  ) +
  #geom_point(data = nwa_dat, aes(x = X, y = Y)) +
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

saveRDS(mod_spatial, "models/mod_spatial.rds")
# 
# mod_spatial_eco_add <- sdmTMB(focal_std_by_all ~ 
#                         year_c  + eco_collapsed+
#                         (1 | study) +
#                         (1  + year_c | trajectory), 
#                       dispformula = ~focalUnit,
#                       family = tweedie(link = "log"),
#                       spatial = "on",
#                       spatial_varying = ~0 + year_c, 
#                       mesh = mesh_coastline,
#                       data = nwa_dat)
# 
# saveRDS(mod_spatial_eco_add, "models/mod_spatial_eco_add.rds")
# 
# 
# mod_spatial_eco <- sdmTMB(focal_std_by_all ~
#                         year_c*eco_collapsed  +
#                         (1 | study) +
#                         (1  + year_c | trajectory),
#                       dispformula = ~focalUnit,
#                       family = tweedie(link = "log"),
#                       spatial = "on",
#                       spatial_varying = ~0 + year_c,
#                       mesh = mesh_coastline,
#                       data = nwa_dat)
# 
# saveRDS(mod_spatial_eco, "models/mod_spatial_eco.rds")

##
# model checks
##
sanity(mod_spatial)

# check predictions
prediction_check_density(mod_spatial)
prediction_check_density(mod_spatial, trans = \(x) log(x+0.01)) +
  labs(x = "std. kelp abundance", subtitle = "grey = simulated, blue = observed")
ggsave("figures/spatial_fit.jpg")


# where are we misaligned?
prediction_check_density(mod_spatial, trans = \(x) log(x+0.01)) +
  facet_wrap(vars(eco_collapsed)) +
  labs(x = "std. kelp abundance", subtitle = "grey = simulated, blue = observed")
ggsave("figures/spatial_fit_eco.jpg")

prediction_check_density(mod_spatial, trans = \(x) log(x+0.01)) +
  facet_wrap(vars(decade)) #good

prediction_check_density(mod_spatial, trans = \(x) log(x+0.01)) +
  facet_wrap(vars(focalUnit)) #good


prediction_check_density(mod_spatial, trans = \(x) log(x+0.01)) +
  facet_grid(vars(decade), vars(eco_collapsed), scale = "free_y") #huh

# randomized quantile residuals
dharma_plot(mod_spatial)


# coefs
tidy(mod_spatial,  conf.int = TRUE)
tidy(mod_spatial, effects = "ran_pars", conf.int = TRUE)
spatial_mod_re <- tidy(mod_spatial, effects = "ran_vals")

perc_change_per_year <- ((exp(coef(mod_spatial)[2]) -1 )*100) |> round(2)
perc_change_per_year
