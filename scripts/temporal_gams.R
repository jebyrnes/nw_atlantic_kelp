#' -----------------------------------------------
#' NWA Timeseries Analysis and Viz
#' to look at change by decade or other temporal variation
#' -----------------------------------------------

library(dplyr)
library(ggplot2)
library(lubridate)
library(glmmTMB)
library(mgcv)
library(gratia)
library(broom.mixed)
library(emmeans)
library(forcats)
library(modelbased)
library(car)

theme_set(theme_bw(base_size = 12))

#reads in nwa data
source("scripts/load_nwa_data.R")
(nwa_dat$year - nwa_dat$year_c)[1] -> mean_year

# where do we have what timeframe?
nwa_dat |>
  group_by(decade, eco_collapsed) |>
  summarize(sites = n_distinct(trajectory)) |>
  ggplot(aes(x = decade, y = eco_collapsed,
             fill = sites)) +
  geom_tile() +
  scale_fill_viridis_c(transform = "log10")

#####
## Model 1. Global trend, no regional variation, trajectories have own wiggliness
#####
## linear
mod_decadal_global <- glmmTMB(ln_focal_std_by_ecoregion ~ 
                         eco_collapsed +
                         year_c*decade + 
                         (1 + year_c |trajectory) + (1|study), 
                       dispformula =~focalUnit,
                       data = nwa_dat )

estimate_slopes(mod_decadal_global, by = "decade", backend = "emmeans") |> plot()

## GAM
mod_time_global <- bam(ln_focal_std_by_ecoregion ~
                         eco_collapsed +
                         s(year_c, k = 10, m = 2, bs = "tp") +
                         s(trajectory, bs = "re") +
                         s(trajectory, year_c, bs = "fs", m = 1, k = 3),
                       method = "REML",
                       family = "gaussian",
                       data = nwa_dat
                         )

estimate_relation(mod_time_global, 
                              by = c("year_c", "eco_collapsed"),
                              length = 300,
                  preserve_range = TRUE) |> 
  plot(show_data = TRUE,
       ribbon = "none")

# show RE
estimate_relation(mod_time_global,
                  data = nwa_dat,
                    # by = c("year_c", "eco_collapsed", "trajectory"),
                     preserve_range = TRUE) |> 
  ggplot(aes(x = year_c, y = Predicted, 
             color = eco_collapsed, group = trajectory)) +
  geom_line() +
  facet_wrap(vars(eco_collapsed))

#####
## Model 2. Global Trend, Regional variation with same wigliness
#####
## linear
mod_decadal_regional <- glmmTMB(ln_focal_std_by_ecoregion ~ 
                                eco_collapsed *
                                year_c*decade + 
                                (1 + year_c |trajectory) + (1|study), 
                              dispformula =~focalUnit,
                              data = nwa_dat )
estimate_slopes(mod_decadal_regional, by = c("decade", "eco_collapsed"), backend = "emmeans") |> 
  plot() +
  facet_wrap(vars(eco_collapsed), scale = "free_y") +
  guides(color = "none")

## GAM
mod_time_regional_global <- bam(ln_focal_std_by_ecoregion ~
                         eco_collapsed +
                           s(year_c, k = 10, m = 2, bs = "tp") +
                           s(year_c, eco_collapsed, k = 10, m = 2, bs = "fs") +
                           s(trajectory, bs = "re") +
                           s(trajectory, year_c, bs = "fs", m = 1, k = 3),
                         method = "REML",
                       family = "gaussian",
                       data = nwa_dat)

estimate_relation(mod_time_regional_global, 
                  by = c("year_c", "eco_collapsed"),
                  length = 300) |> 
  plot(ribbon = "none") 


#####
## Model 3. Regional trends only but with own wigliness, trajectories have own wiggliness
#####
## linear
## GAM
mod_time_regional_global_wiggly <- bam(ln_focal_std_by_ecoregion ~
                                  eco_collapsed +
                                  s(year_c, k = 10, m = 2, bs = "tp") +
                                  s(year_c, by = eco_collapsed, k = 10, m = 1, bs = "tp") +
                                  s(trajectory, bs = "re") +
                                  s(study, bs = "re") +
                                  s(trajectory, year_c, bs = "fs", m = 1, k = 3),
                                method = "REML",
                                family = "gaussian",
                                data = nwa_dat)

estimate_relation(mod_time_regional_global_wiggly, 
                  by = c("year_c", "eco_collapsed"),
                  length = 300) |> 
  plot(ribbon = "none",
       show_data = TRUE) +
  facet_wrap(vars(eco_collapsed)) +
  guides(color = "none") +
  labs(y = "Log Standardized Kelp Abundance") 

#####
## Model 4. Regional trends only with same wigliness, 
#####
## linear
## GAM
mod_time_regional <- bam(ln_focal_std_by_ecoregion ~
                                  eco_collapsed +
                                  s(year_c, eco_collapsed, 
                                    k = 10, m = 2, bs = "fs") +
                                  s(trajectory, bs = "re") +
                           s(trajectory, year_c, bs = "fs", m = 1, k = 3),
                         method = "REML",
                                family = "gaussian",
                                data = nwa_dat )

estimate_relation(mod_time_regional, 
                  by = c("year_c", "eco_collapsed"),
                  length = 300) |> 
  plot(ribbon = "none") 

fitted_values(mod_time_regional) |>
  ggplot(aes(x = year_c, y = .fitted, color = eco_collapsed, group = trajectory)) +
  geom_line() +
  facet_wrap(vars(eco_collapsed))

#####
## Model 5. Regional trends only with unique wigliness, 
#####
## linear
## GAM
mod_time_regional_wiggly <- bam(ln_focal_std_by_ecoregion ~
                                  eco_collapsed +
                                  s(year_c, by = eco_collapsed, 
                                    k = 10, m = 1, bs = "tp") +
                                  s(trajectory, bs = "re") +
                                  s(study, bs = "re") +
                                  s(trajectory, year_c, bs = "fs", m = 1, k = 3),
                                method = "REML",
                                family = "gaussian",
                                data = nwa_dat)

estimate_relation(mod_time_regional_wiggly, 
                  by = c("year_c", "eco_collapsed"),
                  length = 300) |> 
  plot(ribbon = "none") 

#####
## Model Comparisons
#####
aic_table <- tibble(
  modnames = c("One Global Trend",
              # "Global Trend, Regional Variation",
               "Global Trend, Wiggly Regional Variation",
            #   "Regional Trends Only",
               "Regional Trends Varying in Wiggliness"),
  
  mods = list(mod_time_global,
      # mod_time_regional_global,
       mod_time_regional_global_wiggly,
      # mod_time_regional,
       mod_time_regional_wiggly),
  
  df = purrr::map_dbl(mods, ~.x$df.null - .x$df.residual),
  
  AIC = purrr::map_dbl(mods, AIC),
  
  delta_AIC = AIC - min(AIC)
 
) |>
  arrange(delta_AIC) |>
  select(-mods)

aic_table

#####
## Plotting Comparisons
#####

estimate_slopes(mod_decadal_regional, by = c("decade", "eco_collapsed"), backend = "emmeans") |> 
  plot() +
  facet_wrap(vars(eco_collapsed), scale = "free_y") +
  guides(color = "none")

estimate_slopes(mod_decadal_global, by = "decade", backend = "emmeans") |> plot()
