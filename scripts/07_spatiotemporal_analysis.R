#' -----------------------------------------------
#' NWA Timeseries Modeling
#' looking at variation across spacetime 
#' with sdmTMB
#' -----------------------------------------------

library(sdmTMB)
library(sdmTMBextra)
library(sf)
library(ggplot2)
library(patchwork)
source("scripts/99_sdm_helpers.R") #also loads data and constants/helpers

###
## Make objects we will use
##
####
mesh_coastline <- make_kelp_mesh()
valid_areas <- make_valid_areas(buffer = 20)
prediction_points <- full_prediction_points |> 
  st_intersection(valid_areas)

####
## Create basis set and add to data
## See extensive conversation at https://github.com/sdmTMB/sdmTMB/issues/509
####

# Build smoother basis and split into:
# Xs: unpenalized component
# Zs: penalized component
smooth_formula <- focal_std_by_all ~ s(year,  bs = "tp", k = 10)

make_basis <- function(smooth_formula, dat, basis_prev = NULL){
  sm <- sdmTMB:::parse_smoothers(smooth_formula, data = dat, basis_prev = basis_prev)
  
  sx <- as.data.frame(sm$Xs)
  colnames(sx) <- paste0("SX", seq_len(ncol(sx)))
  sz <- as.data.frame(do.call(cbind, sm$Zs))
  colnames(sz) <- paste0("SZ", seq_len(ncol(sz)))
  dat <- cbind(dat, sx, sz)
  
  return(list(sx = sx, sz = sz, dat = dat, sm = sm))
}

nwa_basis <- make_basis(smooth_formula, nwa_dat)
nwa_dat_with_basis <- nwa_basis$dat

# for prediction
prediction_points_yr <- tidyr::expand_grid(prediction_points, 
                                           year = seq(min(nwa_dat_with_basis$year), 
                                                      max(nwa_dat_with_basis$year), 
                                                      by = 0.1))

pred_basis <- make_basis(smooth_formula, prediction_points_yr, basis_prev = nwa_basis$sm)

pred_data <- pred_basis$dat |>  
  as_tibble() |> 
  mutate(trajectory = nwa_dat_with_basis$trajectory[1],
         study = nwa_dat_with_basis$study[1],
         focalUnit = nwa_dat_with_basis$focalUnit[1]
  )


####
## Fit a GAM model
## where the coefficients and intercepts
## vary by space and time
####

##
fe_formula_penalized <- as.formula(paste("focal_std_by_all ~ 1 + ", 
                                         paste(colnames(nwa_basis$sx), collapse = " + "),
                                         "+ (1|study) + (1 +",
                                         paste(colnames(nwa_basis$sx), collapse = " + "),
                                         "| trajectory)"))

svc_formula_penalized <- as.formula(paste("~ 0 +", paste(colnames(nwa_basis$sz), collapse = " + ")))

mod_spatiotemporal<- sdmTMB(fe_formula_penalized, 
                      dispformula = ~focalUnit,
                      family = tweedie(link = "log"),
                      spatial = "on",
                      spatial_varying = svc_formula_penalized, 
                      mesh = mesh_coastline,
                      data = nwa_dat_with_basis,
                      control = 
                        sdmTMBcontrol(map = list(ln_tau_Z = factor(rep(1L, ncol(nwa_basis$sz)))))
)


##
# model checks
##
sanity(mod_spatiotemporal)

# check predictions
prediction_check_density(mod_spatiotemporal)
prediction_check_density(mod_spatiotemporal, trans = \(x) exp(x)+0.01)
#prediction_check_density(mod_spatiotemporal, trans = \(x) log(x+0.01))

# randomized quantile residuals
dharma_plot(mod_spatiotemporal)

##
# pretty plots
##

fitted <- predict(mod_spatiotemporal,
                  type = "response")

ggplot(fitted,
       aes(x = year, y = est, group = trajectory, color = eco_collapsed)) +
  geom_line() +
  scale_color_brewer(palette = "Dark2")+
  facet_wrap(vars(eco_collapsed), scale = "free_y") +
  guides(color="none") +
  geom_point(data = mod_spatiotemporal$data,
             aes(y = focal_std_by_all),
             alpha = 0.1) +
  labs(y = "Log Standardized Kelp Abundance", x = "")

## 
# Show a plot of predictions to 
# emphasize the GAM-nature of the data
##

predicted <- predict(mod_spatiotemporal,
                     type = "response",
                     newdata = pred_data)


ggplot(predicted|> 
         filter_time_eco_collapsed() ,
       aes(x = year, y = est, group = cell_id, color = Y)) +
  geom_line(alpha = 0.1)  +
  scale_color_viridis_c() +
  facet_wrap(vars(eco_collapsed))+
  geom_point(data = mod_spatiotemporal$data,
             aes(y = focal_std_by_all),
             alpha = 0.1, group = 1) +
  labs(y = "Standardized Kelp Abundance", x = "") +
  xlim(c(1975, 2025))

##
# Map predictions over time
##
yrs <- c(1975, 1985, 1995, 2005, 2015, 2020) + 5

predict_decades <-
  predicted |>
  filter(year %in% yrs) |>
  mutate(ind = match(year, yrs)) |>
  group_by(year) |> 
 # mutate(X = (X) + (ind-1)*2e3) |>
  ungroup() |>
  st_as_sf(coords = c("X", "Y"), crs = st_crs(coastline_32619_km))

ggplot() +
  geom_sf(data = coastline_32619_km) +
  geom_sf(data = predict_decades, aes(color = log(est)), size = 0.4) +
  scale_color_viridis_c() +
  facet_wrap(vars(year))

## 
# Get the first derivatives (x-lag(x))
# to show change in rate of change through time
##
