#' -----------------------------------------------
#' NWA Timeseries Modeling
#' to look at change by decade or other temporal variation
#' -----------------------------------------------

library(dplyr)
library(ggplot2)
library(glmmTMB)
library(mgcv)
library(emmeans)
library(modelbased)
library(performance)

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
## Model where we split things up by decade
## But continue to assume linearity
## Note - will have some issues with the specification due to
## not every region being in every decade
#####
mod_decadal_regional <- glmmTMB(ln_focal_std_by_ecoregion ~ 
                                  (eco_collapsed:year_c)*decade + 
                                  (1 + year_c |trajectory) + 
                                  (1|study), 
                                dispformula =~focalUnit,
                                data = nwa_dat  )

check_model(mod_decadal_regional)
estimate_slopes(mod_decadal_regional, trend = "year_c", 
                by = c("decade", "eco_collapsed"), backend = "emmeans") |> plot() +
  facet_wrap(vars(eco_collapsed), scale = "free_y") +
  scale_color_brewer(palette = "Dark2")

# estimate_relation(mod_decadal_regional,
#                   backend = "emmeans",
#                   by = c("year_c", "decade", "eco_collapsed")) |>
#   plot() +
#   scale_color_brewer(palette = "Dark2")

#####
## GAM to see continuous time
#####
mod_time_regional <- bam(ln_focal_std_by_ecoregion ~
                                  eco_collapsed +
                           s(year_c, k = 10, bs = "tp") +
                                  s(year_c, by = eco_collapsed, 
                                    k = 10, bs = "sz") +
                                  s(trajectory, bs = "re") +
                                  s(study, bs = "re") +
                                  s(trajectory, year_c, bs = "fs", k = 3, m = 1),
                                method = "REML",
                                family = "gaussian",
                                data = nwa_dat)


mod_time_region_province <- bam(ln_focal_std_by_province ~
                           eco_collapsed +
                           s(year_c, k = 10, bs = "tp") +
                           s(year_c, by = eco_collapsed, 
                             k = 10, bs = "sz") +
                           s(trajectory, bs = "re") +
                           s(study, bs = "re") +
                           s(trajectory, year_c, bs = "fs", k = 3, m = 1),
                         method = "REML",
                         family = "gaussian",
                         data = nwa_dat)

check_model(mod_time_regional)

estimate_slopes(mod_time_regional, 
                trend = "year_c", 
                by = c("year_c", "eco_collapsed"), backend = "emmeans") |> plot() +
  facet_wrap(vars(eco_collapsed), scale = "free_y") +
  scale_color_brewer(palette = "Dark2")


estimate_relation(mod_time_regional, 
                  by = c("year_c", "eco_collapsed"),
                  length = 300) |> 
  plot(ribbon = "none",
       show_data = TRUE) +
  facet_wrap(vars(eco_collapsed)) +
  guides(color = "none") +
  labs(y = "Log Standardized Kelp Abundance") +
  scale_color_brewer(palette = "Dark2")


##
# Save the models
##

saveRDS(mod_decadal_regional, "models/mod_decadal_regional.rds")
saveRDS(mod_time_regional, "models/mod_time_regional.rds")
