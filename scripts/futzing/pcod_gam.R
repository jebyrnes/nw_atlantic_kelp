library(sdmTMB)
library(mgcv) #needed for parse_smoothers() to work
library(ggplot2)

## Make our mesh
mesh <- make_mesh(pcod, xy_cols = c("X", "Y"), cutoff = 10)

## fit the basic model to see what it looks like
fit <- sdmTMB(
  density ~ s(year, bs = "tp", k = 5),
  data = pcod,
  mesh = mesh,
  family = tweedie(link = "log"),
  spatial = "on"
)

ggeffects::ggpredict(fit, "year [2003:2017, by=0.1]") |> plot()

##
# Explore a spatiotemporal GAM
##

## add a basis set
source("https://raw.githubusercontent.com/sdmTMB/sdmTMB/106d8129168df17b6d94961882c6200281f49db2/R/smoothers.R")

sm <- sdmTMB:::parse_smoothers(density ~ s(year, bs = "tp", k = 5), data = pcod)

s1 <- bind_cols(sm$Xs, sm$Zs)
s2 <- as_tibble(sm$basis_out[[1]][[1]]$X)

pcod <- cbind(pcod, sm$basis_out[[1]][[1]]$X |> as.data.frame())

## Fit the S-T model
fit_gam <- sdmTMB(
  density ~ 1,
  spatial_varying = ~ 0 + V1 + V2 + V3 + V4,
  data = pcod,
  mesh = mesh,
  family = tweedie(link = "log"),
  spatial = "on"
)

##
# visualize
##
# create a year-basis set dataset for prediction
year_df <- gratia::basis(s(year, bs = "tp", k = 5), 
                  data = pcod, 
                  at = data.frame(year = seq(2003, 2017, by = 0.1)),
                  constraint = TRUE) |>
  tidyr::pivot_wider(names_from = ".bf", 
                     values_from = ".value", 
                     names_prefix = "V")

# combine the years and basis set above with all grid values
grid_yrs <- tidyr::expand_grid(qcs_grid,
                               year_df)

# predictions
pred <- predict(fit_gam, newdata = grid_yrs)

# show variation in trajectories while exploring a N-S gradient
ggplot(pred,
       aes(x=year, y = est, 
           group = paste(X, Y))) +
  geom_line(alpha = 0.01) +
  scale_color_viridis_c() #+
  #facet_wrap(vars(cut_interval(Y, 9)))

# nah, not much here

# show the actual years on a map
ggplot(pred |> filter(year == round(year)),
       aes(x=X, y = Y, fill = exp(est))) +
  geom_raster() +
  facet_wrap(vars(year)) +
  scale_fill_viridis_c(
    trans = "sqrt",
    # trim extreme high values to make spatial variation more visible:
    na.value = "yellow", limits = c(0, quantile(exp(pred$est), 0.995))
  )



######
# pcod bayes
#######

library(sdmTMB)
pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 5)
plot(pcod_spde)
d <- pcod
d$scaled_year <- (pcod$year - mean(pcod$year)) / 10

fit <- sdmTMB(
  density ~ scaled_year, 
  data = d,
  mesh = pcod_spde, 
  family = tweedie(link = "log"),
  spatial = "on",
  spatial_varying = ~ 0 + scaled_year
)

# grab the internal parameter list at estimated values:
pars <- sdmTMB::get_pars(fit)
# create a 'map' vector for TMB
# factor NA values cause TMB to fix or map the parameter at the starting value:
kappa_map <- factor(rep(NA, length(pars$ln_kappa)))
fit_mle <- update(
  fit,
  control = sdmTMBcontrol(
    start = list(
      ln_kappa = pars$ln_kappa #<
    ),
    map = list(
      ln_kappa = kappa_map #<
    )
  ),
  do_fit = FALSE #<
)

fitl_bayes <- tmbstan::tmbstan(fit$tmb_obj,
                               iter = 1e3,
                               chains = 2
)

pars_plot <- c("b_j[1]", "b_j[2]", "ln_tau_Z")
bayesplot::mcmc_trace(fitl_bayes, pars = pars_plot)
bayesplot::mcmc_pairs(fitl_bayes, pars = pars_plot)

a <- as.data.frame(fitl_bayes)

cor(a$`b_j[2]`, a$`zeta_s[87]`)


