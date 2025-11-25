#' -----------------------------------------------
#' NWA Timeseries Analysis and Viz
#' to look at change by decade or other temporal variation
#' -----------------------------------------------

library(dplyr)
library(ggplot2)
library(lubridate)
library(glmmTMB)
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

# plot some smooths
ggplot(nwa_dat, 
       aes(x = year, y = ln_focal_std_by_ecoregion,
           color = eco_collapsed)) +
  geom_point(alpha = 0.1) +
  stat_smooth(formula = y ~ s(x, bs = "cs", k = 5)) +
   theme(legend.position = "bottom")

##
#  glmmTMB fits
##

mod_decadal <- glmmTMB(ln_focal_std_by_ecoregion ~ 
                             year_c*decade + 
                             (1 + year_c |trajectory) + (1|study), 
                           dispformula =~focalUnit,
                           data = nwa_dat )

mod_decadal_eco <- glmmTMB(ln_focal_std_by_ecoregion ~ 
                         eco_collapsed*year_c*decade + 
                         (1 + year_c |trajectory) + (1|study), 
                       dispformula =~focalUnit,
                       data = nwa_dat |> 
                         filter(eco_collapsed != "Virginian") |>
                         filter(decade != "pre-1980s"))

Anova(mod_decadal_eco)


estimate_slopes(mod_decadal,
                trend = "year_c",
                by = c("decade"),
                backend = "emmeans",
                include_random = FALSE,
                ci = 0.66) |> 
  plot() +
  coord_flip()


estimate_slopes(mod_decadal_eco,
                trend = "year_c",
                by = c("decade", "eco_collapsed"),
                backend = "emmeans",
                include_random = FALSE) |> 
  plot() +
  coord_flip() +
  labs(x = "", y = "slope")


modelbased::estimate_relation(mod_decadal,
                              by = c("year_c", "decade"),
                              length = 200) |>
  plot(show_data = TRUE)

emtrends(mod_decadal_eco, ~ decade | eco_collapsed, "year_c")

modelbased::estimate_relation(mod_decadal_eco,
                              by = c("year_c", "decade", "eco_collapsed"
                                     ),
                              range = "grid",
                              preserve_range = FALSE) |>
  plot(show_data = TRUE,
       ribbon = list(fill = NA)) +
  ylim(c(-6,5))


##
# Let's try a gam....
##

mod_gam_decadal <- mgcv::bam(ln_focal_std_by_ecoregion ~ 
                               eco_collapsed + 
                                   s(year_c, k = 15, m = 2, bs = 'tp') +
                                   s(trajectory, bs = "re") +
                                   s(trajectory, year_c, bs = "re") + #random slope
                                   s(study, bs = "re"),
                                 data = nwa_dat)

mod_gam_decadal_eco <- mgcv::bam(ln_focal_std_by_ecoregion ~ 
                     s(year_c, by = eco_collapsed, k = 15, m = 2, bs = 'tp') +
                       s(trajectory, bs = "re") +
                       s(trajectory, year_c, bs = "re") + #random slope
                       s(study, bs = "re"),
                       data = nwa_dat)





gratia::fitted_values(mod_gam_decadal, 
                      scale = "response",
                      terms = c("(Intercept)", "eco_collapsed", "s(year_c)"))|>
  mutate(year = year_c + mean_year) |>
  ggplot(aes(x = year, y = .fitted, group = eco_collapsed)) +
  geom_line() +
  geom_ribbon(aes(group = eco_collapsed, 
                  ymin = .lower_ci,
                  ymax = .upper_ci),
              alpha = 0.4) +
  geom_point(data = nwa_dat, aes(y = ln_focal_std_by_ecoregion),
             alpha = 0.2)



gratia::fitted_values(mod_gam_decadal_eco, 
                      scale = "response",
                      terms = c("(Intercept)", "eco_collapsed", "s(year_c)"))|>
  mutate(year = year_c + mean_year) |>
  ggplot(aes(x = year, y = .fitted, group = eco_collapsed)) +
  geom_line() +
  geom_ribbon(aes(group = eco_collapsed, 
                  ymin = .lower_ci,
                  ymax = .upper_ci),
              alpha = 0.4) +
  geom_point(data = nwa_dat, aes(y = ln_focal_std_by_ecoregion),
             alpha = 0.2) +
  facet_wrap(vars(eco_collapsed))

change_rates <- gratia::derivatives(mod_gam_decadal) |>
  mutate(.smooth = gsub("s\\(year_c\\)", "", .smooth))

change_rates_eco <- gratia::derivatives(mod_gam_decadal_eco) |>
  mutate(.smooth = gsub("s\\(year_c\\):eco_collapsed", "", .smooth))


ggplot(change_rates,
       aes(x = year_c+mean_year, fill = .smooth, y = .derivative,
           ymin = .lower_ci, ymax = .upper_ci)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, color = "red", lty = 2) +
  theme(legend.position = "bottom") +
  xlim(c(1980, 2020)) +
  labs(color = "", fill = "")

ggplot(change_rates_eco,
       aes(x = year_c+mean_year, color = .smooth,
           fill = .smooth, y = .derivative,
           ymin = .lower_ci, ymax = .upper_ci)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, color = "red", lty = 2) +
  theme(legend.position = "bottom") +
  xlim(c(1980, 2020)) +
  labs(color = "", fill = "") +
  facet_wrap(vars(.smooth), scales = "free_y") 


####
decadal_slopes <- estimate_slopes(mod_decadal,
                                  trend = "year",
                                  by = "decade",
                                  backend = "emmeans",
                                  ci = 0.9) 

plot(decadal_slopes) +
  coord_flip() +
  scale_y_continuous(transform = coef_to_prec) +
  labs(x="", y = "% change/yr in standardized kelp")

ggsave("figures/decadal_slopes.jpg", width = 5, height = 5)


## Bayesian
library(brms)
mod_decadal_brm <- brm(brmsformula(ln_focal_std_by_ecoregion ~ 
                                   eco_collapsed*year_c*decade + 
                                   (1 + year_c |trajectory) + 
                                   (1|study),
                                 sigma ~ focalUnit),
                     data = nwa_dat)
