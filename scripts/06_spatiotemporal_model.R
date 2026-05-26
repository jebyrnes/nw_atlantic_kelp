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
  v <- as.data.frame(sm$basis_out[[1]][[1]]$X)
  colnames(v) <- paste0("SV", seq_len(ncol(v)))
  dat <- cbind(dat, sx, sz, v)
  
  return(list(sx = sx, sz = sz, dat = dat, sm = sm, v = v))
}

# the smoother
smooth_formula <- focal_std_by_all ~ s(year,  bs = "tp", k = 5)


nwa_basis <- make_basis(smooth_formula, nwa_dat)
nwa_dat_with_basis <- nwa_basis$dat


####
## Fit a GAM model
## where the coefficients and intercepts
## vary by space and time
####

## Option 1 from https://github.com/sdmTMB/sdmTMB/issues/509
## Penalized smoother, then deviations from the smoother
fe_formula <- as.formula(paste("focal_std_by_all ~ 1 + ", 
                                         paste(colnames(nwa_basis$v), collapse = " + "),
                                         "+ (1|study) + (1 | trajectory)"))

svc_formula <- as.formula(paste("~ 0 +", 
                                paste(colnames(nwa_basis$sx), collapse = " + "), "+",
                                paste(colnames(nwa_basis$sz), collapse = " + ")))

mod_spatiotemporal<- sdmTMB(fe_formula,  
                      dispformula = ~focalUnit,
                      family = tweedie(link = "log"),
                      spatial = "on",
                      spatial_varying = svc_formula, 
                      mesh = mesh_coastline,
                      data = nwa_dat_with_basis,
                      control = 
                        sdmTMBcontrol(map = list(ln_tau_Z = factor(c(1L, rep(2L, ncol(nwa_basis$sz))))))
)

saveRDS(mod_spatiotemporal, "models/mod_spatiotemporal.rds")

##
# model checks
##
sanity(mod_spatiotemporal)

# check predictions
prediction_check_density(mod_spatiotemporal)

prediction_check_density(mod_spatiotemporal, trans = \(x) log(x+0.01)) +
  labs(x = "std. kelp abundance", subtitle = "grey = simulated, blue = observed")
ggsave("figures/spatiotemporal_fit_diagnostic.jpg")


prediction_check_density(mod_spatiotemporal, trans = \(x) log(x+0.01)) +
  facet_wrap(vars(eco_collapsed))+
  labs(x = "std. kelp abundance", subtitle = "grey = simulated, blue = observed")
ggsave("figures/spatiotemporal_fit_diagnostic_eco.jpg")

# randomized quantile residuals
dharma_plot(mod_spatiotemporal)

##
# Coefficients
##
tidy(mod_spatiotemporal)
tidy(mod_spatiotemporal, "ran_pars")
tidy(mod_spatiotemporal, "dispersion")

