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
mesh_coastline <- make_kelp_mesh(cutoff = 10)
valid_areas <- make_valid_areas(buffer = 10)
prediction_points <- full_prediction_points |> 
  st_intersection(valid_areas)

####
## Create basis set and add to data
## See extensive conversation at https://github.com/sdmTMB/sdmTMB/issues/509
####

# Build smoother basis and split into:
# Xs: unpenalized component
# Zs: penalized component

make_basis <- function(smooth_formula, dat, basis_prev = NULL){
  sm <- sdmTMB:::parse_smoothers(smooth_formula, data = dat, basis_prev = basis_prev)
  
  sx <- as.data.frame(sm$Xs)
  colnames(sx) <- paste0("SX", seq_len(ncol(sx)))
  sz <- as.data.frame(do.call(cbind, sm$Zs))
  colnames(sz) <- paste0("SZ", seq_len(ncol(sz)))
  v <- as_tibble(sm$basis_out[[1]][[1]]$X)
  dat <- cbind(dat, sx, sz, v)
  
  return(list(sx = sx, sz = sz, dat = dat, sm = sm, v = v))
}

# the smoother
smooth_formula <- focal_std_by_all ~ s(year,  bs = "tp", k = 5)


nwa_basis <- make_basis(smooth_formula, nwa_dat)
nwa_dat_with_basis <- nwa_basis$dat

# for prediction
prediction_points_yr <-
  tidyr::expand_grid(prediction_points, year = seq(
    min(nwa_dat_with_basis$year),
    max(nwa_dat_with_basis$year),
    by = 0.1
  ))

pred_basis <- make_basis(smooth_formula, prediction_points_yr, basis_prev = nwa_basis$sm)

pred_data <- pred_basis$dat |>
  as_tibble() |>
  mutate(
    trajectory = nwa_dat_with_basis$trajectory[1],
    study = nwa_dat_with_basis$study[1],
    focalUnit = nwa_dat_with_basis$focalUnit[1]
  )


####
## Fit a GAM model
## where the coefficients and intercepts
## vary by space and time
####

##
fe_formula_unpenalized <- as.formula(paste("focal_std_by_all ~ 1 + ", 
                                           paste(colnames(nwa_basis$v), collapse = " + "),
                                           "+ (1|study) + (1  |trajectory)"))

svc_formula_unpenalized <- as.formula(paste("~ 0 +", paste(colnames(nwa_basis$v), collapse = " + ")))

fe_formula_penalized <- as.formula(paste("focal_std_by_all ~ 1 + ", 
                                         paste(colnames(nwa_basis$sx), collapse = " + "),
                                         "+ (1|study) + (1 +",
                                         paste(colnames(nwa_basis$sx), collapse = " + "),
                                         "| trajectory)"))

svc_formula_penalized <- as.formula(paste("~ 0 +", paste(colnames(nwa_basis$sz), collapse = " + ")))

mod_spatiotemporal<- sdmTMB(fe_formula_unpenalized, 
                      dispformula = ~focalUnit,
                      family = tweedie(link = "log"),
                      spatial = "on",
                      spatial_varying = svc_formula_unpenalized, 
                      mesh = mesh_coastline,
                      data = nwa_dat_with_basis,
                      control = 
                        sdmTMBcontrol(map = list(ln_tau_Z = factor(rep(1L, ncol(nwa_basis$v)))))
)


##
# model checks
##
sanity(mod_spatiotemporal)

# check predictions
prediction_check_density(mod_spatiotemporal)
prediction_check_density(mod_spatiotemporal, trans = \(x) log(x+0.01))

# randomized quantile residuals
dharma_plot(mod_spatiotemporal)

##
# Coefficients
##
tidy(mod_spatiotemporal)
tidy(mod_spatiotemporal, "ran_pars")
tidy(mod_spatiotemporal, "dispersion")

##
# pretty plots
##

fitted <- predict(mod_spatiotemporal,
                  type = "response")

ggplot(fitted,
       aes(x = year, y = est, group = trajectory, color = eco_collapsed)) +
  geom_line() +
  scale_color_brewer(palette = "Dark2")+
  facet_wrap(vars(eco_collapsed)) +
  guides(color="none") +
  geom_point(data = mod_spatiotemporal$data,
             aes(y = focal_std_by_all),
             alpha = 0.1) +
  labs(y = "Standardized Kelp Abundance", x = "")

ggsave("figures/spatiotemporal_fitted.jpg", 
       width = 7, height = 6)

## 
# Show a plot of predictions to 
# emphasize the GAM-nature of the data
##

predicted <- predict(mod_spatiotemporal,
                     type = "response",
                     newdata = pred_data)

predicted <- predicted |>
  mutate(eco_collapsed = 
           forcats::fct_relevel(eco_collapsed,
                                c("Virginian" ,
                                  "Gulf of Maine/Bay of Fundy" ,
                                  "Scotian Shelf" ,
                                  "Gulf of St. Lawrence - Newfoundland")),
         ord = as.numeric(eco_collapsed)*1e5 + Y,
         id = rank(ord)
  )

ggplot(predicted|> 
         filter_time_eco_collapsed()  ,
       aes(x = year, y = est, group = cell_id, color = Y)) +
  geom_line(alpha = 0.1)  +
  scale_color_viridis_c() +
  # geom_point(data = mod_spatiotemporal$data,
  #            aes(y = focal_std_by_all),
  #            alpha = 0.1, group = 1) +
   labs(y = "Standardized Kelp Abundance", x = "",
        color = "Northing") +
  theme_light(base_size = 18)

ggsave("figures/spatiotemporal_curves.jpg", 
       width = 8, height = 5)

# hovmoller
ggplot(predicted,
       aes(x = year, y = id, fill = est, color = est)) +
  geom_tile() +
  #scale_fill_viridis_c(option = "D") +
  colorspace::scale_fill_continuous_divergingx(palette = "BrBG", mid = 0.2, rev = FALSE) +
  colorspace::scale_color_continuous_divergingx(palette = "BrBG", mid = 0.2, rev = FALSE) +
  facet_wrap(vars(forcats::fct_rev(eco_collapsed)), 
              scale = "free_y", ncol = 1)+
  guides(y = "none") +
  labs(y = "", x = "", fill = "Standardized\nKelp Abundance") +
  theme_classic() 

##
# Map predictions over time
##
yrs <- c(1975, 1985, 1995, 2005, 2015) + 5

predict_decades <-
  predicted |>
  filter(year %in% yrs) |>
  mutate(ind = match(year, yrs)) |>
  group_by(year) |> 
 # mutate(X = (X) + (ind-1)*2e3) |>
  ungroup() |>
  st_as_sf(coords = c("X", "Y"), crs = st_crs(coastline_32619_km))

yrs_plot <- list(c(1980, 1990), c(1990, 2000),
            c(2000, 2010), c(2010, 2020))
ggplot() +
  geom_sf(data = coastline_32619_km) +
  geom_sf(data = predict_decades , 
          aes(color = est), size = 1) +
  scale_color_viridis_c(option = "A", transform = "logit") +
  facet_wrap(vars(year)) +
  labs(color = "Standardized\nKelp Abundance")


ggsave("figures/spatiotemporal_maps.jpg", 
       width = 8, height = 5)
## 
# Get the first derivatives (x-lag(x))
# to show change in rate of change through time
##

get_derivs <- function(dat){
  dat |>
    group_by(X, Y, ecoregion, trajectory, study, geometry) |>
    arrange(year) |>
    mutate(est_change = lead(est) - lag(est),
           log_est_change = log(lead(est)+1e-7) - log(lag(est)+1e-7),
           est_accel = lead(est_change) - lag(est_change) ) |>
    ungroup()
}



ggplot(predicted|> get_derivs(),
       aes(x = year, y = est_change, group = cell_id, color = Y)) +
  geom_line(alpha = 0.1)  +
  scale_color_viridis_c() +
 # facet_wrap(vars(eco_collapsed)) +
  labs(y = "Rate of Change", x = "") +
  geom_hline(yintercept = 0, lty = 2, lwd = 0.3) +
  xlim(c(1970, 2025)) +
  scale_y_continuous(transform = "pseudo_log")



ggplot(predicted|> get_derivs(),
       aes(x = year, fill = est_change, y = id)) +
  geom_tile()  +
  colorspace::scale_fill_continuous_divergingx(palette = 'RdYlBu', mid = 0) +
  facet_wrap(vars(forcats::fct_rev(eco_collapsed)), 
             scale = "free_y", ncol = 1)+
  labs(y = "", x = "", fill = "Rate of Change") +
  #xlim(c(1970, 2020)) +
  guides(y = "none")

ggsave("figures/spatiotemporal_rate_of_change_hovmoller.jpg", 
       width = 5, height = 7)
