#' -----------------------------------------------
#' NWA Timeseries Analysis and Viz
#'
#' -----------------------------------------------

library(dplyr)
library(ggplot2)
library(lubridate)
library(glmmTMB)
library(broom.mixed)
library(emmeans)
library(forcats)

theme_set(theme_bw(base_size = 12))

nwa <- read.csv("data/nwa_with_env.csv")

# initial plot
ggplot(nwa,
       aes(x = year, y = ln_focal_std_by_ecoregion,
           group = trajectory, color = study)) +
  geom_point(alpha = 0.1) +
  facet_wrap(vars(ecoregion)) +
  scale_color_discrete(guide = "none") +
  stat_smooth(method = "lm", fill = NA) +
  ylim(c(0,1.1))


ggplot(nwa,
       aes(x = year, y = ln_focal_std_by_ecoregion)) +
  geom_point(alpha = 0.1, aes(color = study)) +
  scale_color_discrete(guide = "none") +
  stat_smooth(method = "glm", fill = NA, 
              method.args = list(family = gaussian(link = "identity")),
              aes(group = ecoregion),
              color = "black") +
  ylim(c(0,1.1))



# Model
mod_ecoregion <- glmmTMB(lfse ~ 
                 ecoregion*year + 
                 (1 + year |trajectory) + (1|study), 
               dispformula =~focalUnit,
               family = ordbeta,
               data = nwa |>
                 mutate(lfse = scales::rescale(focal_std_by_ecoregion, c(0,1))))

mod_common_slope <- glmmTMB(lfse ~ 
                              ecoregion + year + 
                              (1 + year |trajectory) + (1|study), 
                            dispformula =~focalUnit,
                            family = ordbeta,
                            data = nwa |>
                              mutate(lfse = scales::rescale(focal_std_by_ecoregion, c(0,1))))

# check model
performance::check_convergence(mod_ecoregion)
performance::check_convergence(mod_common_slope)

performance::check_predictions(mod_ecoregion)
performance::check_predictions(mod_common_slope)

performance::check_residuals(mod_ecoregion) |> plot()
performance::check_residuals(mod_common_slope)


# What's the model tell us?
AIC(mod_no_eco)
AIC(mod) # it's this one

modelbased::estimate_grouplevel(mod_no_eco,
                                type = "total",
                                group = "trajectory:year") |>
  filter(!grepl("Intercept", Parameter))|> 
  mutate(Level = as.character(Coefficient)) |>
  plot() +
  guides(y = "none") +
  geom_hline(yintercept = 0, lty = 2, color = "red")

car::Anova(mod)
car::Anova(mod_no_eco)
tidy(mod, effects = "fixed")

# Plot slopes
slopes <- emtrends(mod, specs = ~ecoregion, var = "year")

plot(slopes)

# Show the model results
modelbased::estimate_relation(mod,
                              by = c("year", "ecoregion")) |>
  plot(show_data = TRUE)

##
# Percent Cover Only with Ordbetareg
##
perc_mod <- glmmTMB(p_cover ~ 
                      ecoregion*year + 
                      (1  |trajectory) + (1|study), 
                    family = ordbeta,
                    data = nwa |>
                      mutate(p_cover = percent_cover/100))

modelbased::estimate_expectation(perc_mod,
                              by = c("year"),
                              showdata = TRUE) |>
  plot(show_data = TRUE)

##
# brms
##
library(brms)
library(ordbetareg)

mod_ecoregion <- ordbetareg(
  bf(fse ~ ecoregion*year_c + 
       (1 + ecoregion*year_c |trajectory) + (1|study)),
  data = nwa |>
    mutate(lfse = scales::rescale(focal_std_by_ecoregion, c(0,1))),
                     cores = 4)


mod_ecoregion <- brm(
  bf(ln_focal_std_by_ecoregion ~ ecoregion*year_c + 
       (1 + ecoregion*year_c |trajectory) + (1|study),
     sigma ~ focalUnit),
  data = nwa,
  cores = 4,
  file = "models/mod_ecoregion.rds")

                                                       