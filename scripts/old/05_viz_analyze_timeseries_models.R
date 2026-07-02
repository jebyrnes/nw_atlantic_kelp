#' -----------------------------------------------
#' NWA Timeseries Analysis and Viz
#' to look at change by decade or other temporal variation
#' -----------------------------------------------

library(dplyr)
library(ggplot2)
library(glmmTMB)
library(mgcv)
library(gratia)
library(broom.mixed)
library(emmeans)
library(modelbased)
library(car)

theme_set(theme_bw(base_size = 12))

#reads in nwa data
source("scripts/load_nwa_data.R")
(nwa_dat$year - nwa_dat$year_c)[1] -> mean_year

# Read in fit models
mod_decadal_global <- readRDS("models/mod_decadal_global.rds")
mod_time_global <- readRDS("models/mod_time_global.rds")
mod_decadal_regional_global <- readRDS("models/mod_decadal_regional_global.rds")
mod_time_regional_global_wiggly <- readRDS("models/mod_time_regional_global_wiggly.rds")
mod_decadal_regional <- readRDS("models/mod_decadal_regional.rds")
mod_time_regional_wiggly<- readRDS("models/mod_time_regional_wiggly.rds")

#####
## Model Comparisons
#####

## Decadal
aic_table_decadal <- tibble(
  modnames = c("One Global Trend",
               "Global Trend, Regional Variation",
               "Regional Trends Only"),
  
  mods = list(mod_decadal_global,
              mod_decadal_regional_global,
              mod_decadal_regional),
  
  #df = purrr::map_dbl(mods, ~.x$df.null - .x$df.residual),
  
  AIC = purrr::map_dbl(mods, AIC),
  
  delta_AIC = AIC - min(AIC)
  
) |>
  arrange(delta_AIC) |>
  select(-mods)

aic_table_decadal

## GAMS
aic_table_gams <- tibble(
  modnames = c("One Global Trend",
               "Global Trend, Wiggly Regional Variation",
               "Regional Trends Varying in Wiggliness"),
  
  mods = list(mod_time_global,
       mod_time_regional_global_wiggly,
       mod_time_regional_wiggly),
  
  df = purrr::map_dbl(mods, ~.x$df.null - .x$df.residual),
  
  AIC = purrr::map_dbl(mods, AIC),
  
  delta_AIC = AIC - min(AIC)
 
) |>
  arrange(delta_AIC) |>
  select(-mods)

aic_table_gams

#####
## Plotting Comparisons
#####

## How different are the two GAM predictions
modelbased::estimate_relation(mod_time_regional_global_wiggly,
                              by = c("year_c", "eco_collapsed"),
                              length = 200) |>
  plot(show_data = TRUE,
       point = list(alpha = 0.1),
       line = list(size = 1.5)) +
  facet_wrap(vars(eco_collapsed), scale = "free_y") +
  guides(color = "none", fill = "none") +
  labs(y = "Log Standardized Kelp Abundance", x = "") +
  scale_x_continuous(labels = \(x) round(x+mean_year))

modelbased::estimate_relation(mod_time_regional_wiggly,
                              by = c("year_c", "eco_collapsed"),
                              length = 200) |>
  plot(show_data = TRUE,
       point = list(alpha = 0.1),
       line = list(size = 1.5)) +
  facet_wrap(vars(eco_collapsed), scale = "free_y")+
  guides(color = "none", fill = "none") +
  labs(y = "Log Standardized Kelp Abundance", x = "") +
  scale_x_continuous(labels = \(x) round(x+mean_year))



# plot slopes
#term_names(mod_time_regional_global_wiggly)
nd <- nwa_dat |>
  group_by(eco_collapsed) |>
  reframe(year_c = 
            seq(from = min(year_c), 
                to = max(year_c), 
                length.out = 1000),
          study = study[1],
          trajectory = trajectory[1])

kelp_change <- response_derivatives(mod_time_regional_global_wiggly,
                                    data = nd,
                                   type = "forward",
                                   focal = "year_c",
                                   eps = 1e-03,
                                   level = 0.9,
                                   n_sim = 100,
                                   exclude = c("trajectory", "study"))

ggplot(kelp_change, 
       aes(x = year_c + mean_year, y = .derivative,
           ymin = .lower_ci, ymax = .upper_ci,
           color = eco_collapsed, fill = eco_collapsed)) +
  facet_wrap(vars(eco_collapsed), scale = "free_y") +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  geom_hline(yintercept = 0, lty = 2) +
  guides(color = "none", fill = "none") +
  labs(color = "", fill = "", y = "Change in Ln Std. Kelp", x = "") 

## Modelbased

kelp_change <- estimate_slopes(mod_time_regional_global_wiggly,
                               trend = "year_c",
                               by = c("year_c", "eco_collapsed"),
                               estimate = "average",
                               backend = "emmeans",
                               length = 1e2,
                               nuissance = c("trajectory", "study"))

plot(kelp_change) +
  facet_wrap(vars(eco_collapsed), scale = "free_y") +
  guides(color = "none", fill = "none") +
  geom_hline(yintercept = 0,  lty = 2) +
  labs(color = "", fill = "", y = "Change in Ln Std. Kelp", x = "") +
  scale_x_continuous(labels = \(x) round(x+mean_year),
                     limits = c(-40, 22))

kelp_change <- modelbased::estimate_slopes(mod_time_regional_global_wiggly, 
                            trend = "year_c",
                            by =  "eco_collapsed") 
kelp_change |> plot()

rate_of_change <- gratia::derivatives(mod_time_regional_global_wiggly, 
                                      type = "central") |>
  filter(!grepl("trajectory", .smooth))

ggplot(rate_of_change,
       aes(x = year_c+mean_year,  y = .derivative,
           ymin = .lower_ci, ymax = .upper_ci)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, color = "red", lty = 2) +
  theme(legend.position = "bottom") +
  xlim(c(1975, 2020)) +
  labs(color = "", fill = "", y = "Change in Ln Std. Kelp", x = "") +
  facet_wrap(vars(.smooth), scale = "free_y")


a <- estimate_relation(mod_time_regional_global_wiggly, data = nwa_dat)
ggplot(a, aes(x = year_c + mean_year, y = Predicted,
              group = trajectory, color = eco_collapsed)) +
  geom_line(alpha = 0.5) +
  facet_wrap(vars(eco_collapsed)) +
  guides(color = "none")

## Global

modelbased::estimate_relation(mod_decadal_global,
                              by = c("year_c", "decade"),
                              length = 200) |>
  plot(show_data = TRUE)

estimate_slopes(mod_decadal_global,
                trend = "year_c",
                by = c("decade"),
                backend = "emmeans",
                include_random = FALSE) |> 
  plot() +
  coord_flip()



# GAM
estimate_relation(mod_time_global, 
                  by = c("year_c", "eco_collapsed"),
                  length = 300) |> 
  plot(#ribbon = "none",
       show_data = TRUE,
       point = list(alpha = 0.1),
       line = list(size = 1.5)) +
  facet_wrap(vars(eco_collapsed)) +
  guides(color = "none", fill = "none") +
  labs(y = "Log Standardized Kelp Abundance") 


# plot slopes
global_change <- gratia::derivatives(mod_time_global, n = 200) |>
  filter(.smooth == "s(year_c)")


ggplot(global_change,
       aes(x = year_c+mean_year,  y = .derivative,
           ymin = .lower_ci, ymax = .upper_ci)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, color = "red", lty = 2) +
  theme(legend.position = "bottom") +
  xlim(c(1975, 2020)) +
  labs(color = "", fill = "", y = "Change in Ln Std. Kelp", x = "")



## Regional
estimate_slopes(mod_decadal_regional, trend = "year_c", 
                by = c("decade", "eco_collapsed"), 
                backend = "emmeans") |> 
  plot() +
  facet_wrap(vars(eco_collapsed), scale = "free_y") +
  guides(color = "none") +
  labs(y = "Slope of Trend") +
  coord_flip()

estimate_slopes(mod_decadal_global, by = "decade", backend = "emmeans") |> plot()
