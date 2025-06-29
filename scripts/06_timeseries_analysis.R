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

#reads in nwa data
source("scripts/load_nwa_data.R")

# initial plot
ggplot(nwa_dat,
       aes(x = year, y = ln_focal_std_by_ecoregion,
           group = trajectory, color = study)) +
  geom_point(alpha = 0.1) +
  facet_wrap(vars(ecoregion)) +
  scale_color_discrete(guide = "none") +
  stat_smooth(method = "lm", fill = NA, size = 0.5) +
  ylim(c(0,1.1))


ggplot(nwa_dat,
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
               REML = TRUE,
               data = nwa_dat |>
                 mutate(lfse = scales::rescale(focal_std_by_ecoregion, c(0,1))))

mod_common_slope <- glmmTMB(lfse ~ 
                              ecoregion + year + 
                              (1 + year |trajectory) + (1|study), 
                            dispformula =~focalUnit,
                            family = ordbeta,
                            REML = TRUE,
                            data = nwa_dat |>
                              mutate(lfse = scales::rescale(focal_std_by_ecoregion, c(0,1))))

# check model
performance::check_convergence(mod_ecoregion)
performance::check_convergence(mod_common_slope)

performance::check_predictions(mod_ecoregion)
performance::check_predictions(mod_common_slope)

performance::check_residuals(mod_ecoregion) |> plot()
performance::check_residuals(mod_common_slope) |> plot()


# What's the model tell us?
AIC(mod_ecoregion)
AIC(mod_common_slope) 

car::Anova(mod_ecoregion)

modelbased::estimate_grouplevel(mod_common_slope,
                                type = "total",
                                group = "trajectory:year") |>
  filter(!grepl("Intercept", Parameter))|> 
  mutate(Level = forcats::fct_reorder(Coefficient |> as.factor(), 
                                      Coefficient)) |>
  plot() +
  guides(y = "none") +
  geom_hline(yintercept = 0, lty = 2, color = "red")

car::Anova(mod_ecoregion)
car::Anova(mod_common_slope)
tidy(mod_ecoregion, effects = "fixed")

# Plot slopes
slopes <- emtrends(mod_ecoregion, specs = ~ecoregion, var = "year")

plot(slopes) +
  labs(y = "", x = "Logit Change")

# Show the model results
modelbased::estimate_relation(mod_ecoregion,
                              length = 100,
                              by = c("year", "ecoregion")) |>
  plot(show_data = TRUE) +
  labs(y = "Porpotion of Max Kelp", x = "",
       color = "", fill = "") +
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  theme_bw(base_size = 14)+
  theme(legend.position = "bottom",
        legend.box="vertical")

ggsave("figures/timeseries_with_curves.jpg",
       width = 7, height = 6)

# show with RE

modelbased::estimate_relation(mod_common_slope,
                              include_random = TRUE,
                              length = 200,
                              by = c("year", "trajectory")) |>
  plot(show_data = TRUE, 
       color = "grey",
       ribbon = list(alpha = 0)) +
  labs(y = "Porpotion of Max Kelp", x = "",
       color = "", fill = "") +
  guides(color = "none", fill = "none") +
  # guides(color=guide_legend(nrow=2,byrow=TRUE),
  #        fill=guide_legend(nrow=2,byrow=TRUE))+
  theme_bw(base_size = 14)+
  theme(legend.position = "bottom",
        legend.box="vertical")

##
# Percent Cover Only with Ordbetareg
##
perc_mod <- glmmTMB(p_cover ~ 
                      ecoregion*year_c + 
                      (1  + year_c|trajectory) + (1|study), 
                    family = ordbeta,
                    REML = FALSE,
                    data = nwa_dat)

check_convergence(perc_mod)
car::Anova(perc_mod)


# Show the model results
modelbased::estimate_relation(perc_mod,
                              by = c("year_c", "ecoregion")) |>
  plot(show_data = TRUE) +
  labs(y = "Porpotion Cover", x = "",
       color = "", fill = "") +
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  theme_bw(base_size = 14)+
  theme(legend.position = "bottom",
        legend.box="vertical") +
  scale_x_continuous(breaks = seq(-31, 11, 10),
                     labels = seq(-31, 11, 10)+
                       mean(nwa$year) |> round())

##
# brms
##
library(brms)
library(ordbetareg)

mod_ecoregion <- ordbetareg(
  bf(fse ~ ecoregion*year_c + 
       (1 + ecoregion*year_c |trajectory) + (1|study)),
  data = nwa_dat |>
    mutate(lfse = scales::rescale(focal_std_by_ecoregion, c(0,1))),
                     cores = 4)


mod_ecoregion <- brm(
  bf(ln_focal_std_by_ecoregion ~ ecoregion*year_c + 
       (1 + ecoregion*year_c |trajectory) + (1|study),
     sigma ~ focalUnit),
  data = nwa_dat,
  cores = 4,
  file = "models/mod_ecoregion.rds")

##
# Map viz?
##

