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
library(modelbased)
library(car)

theme_set(theme_bw(base_size = 12))


# for better coefficient understanding
coef_to_prec<-scales::new_transform(
  "coef_to_perc",
  transform = \(.x) exp(.x)-1,
  inverse = \(.y) log(.y+1) 
)
#reads in nwa data
source("scripts/load_nwa_data.R")

# initial plot
ggplot(nwa_dat,
       aes(x = year, y = focal_std_by_ecoregion,
           group = trajectory, color = study)) +
  geom_point(alpha = 0.1) +
  facet_wrap(vars(eco_collapsed)) +
  scale_color_discrete(guide = "none") +
  stat_smooth(method = "lm", fill = NA, size = 0.5) +
  ylim(c(0,1.1))


ggplot(nwa_dat,
       aes(x = year, y = focal_std_by_ecoregion)) +
  geom_point(alpha = 0.1, aes(color = study)) +
  scale_color_discrete(guide = "none") +
  stat_smooth(method = "glm", fill = NA, 
              method.args = list(family = gaussian(link = "identity")),
              aes(group = eco_collapsed),
              color = "black") +
  ylim(c(0,1.1))



# Model
mod_ecoregion <- glmmTMB(ln_focal_std_by_ecoregion ~ 
                  eco_collapsed*year + 
                 (1 + year |trajectory) + (1|study), 
               dispformula =~focalUnit,
               data = nwa_dat)

mod_common_slope <- glmmTMB(ln_focal_std_by_ecoregion ~ 
                              eco_collapsed+year + 
                              (1 + year |trajectory) + (1|study), 
                            dispformula =~focalUnit,
                            data = nwa_dat )


mod_ecoregion_ord <- glmmTMB(rescaled_std_by_ecoregion ~ 
                           eco_collapsed*year + 
                           (1 + year |trajectory) + (1|study), 
                         dispformula =~focalUnit,
                         family = ordbeta,
                         data = nwa_dat)

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

change_est <- tidy(mod_common_slope) |> filter(term == "year") |> pull(estimate)
perc_change_per_year <- ((exp(change_est) -1 )*100) |> round(2)

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
slopes <- emtrends(mod_ecoregion, specs = ~eco_collapsed, var = "year")

plot(slopes) +
  labs(y = "", x = "standardized change")

# Show the model results
modelbased::estimate_relation(mod_ecoregion,
                              length = 100,
                              by = c("year", "eco_collapsed")) |>
  as_tibble() |>
  mutate(Predicted = exp(Predicted)-0.01,
         CI_low = exp(CI_low)-0.01,
         CI_high = exp(CI_high)-0.01,
         CI_high = ifelse(CI_high>1.1, 1.1, CI_high)
  ) |>
  ggplot(aes(x = year)) +
  geom_line(aes(y = Predicted, 
                color = eco_collapsed),
            size = 1.5) +
  geom_ribbon(aes(y = Predicted, 
                  ymin = CI_low,
                  ymax = CI_high,
                  fill = eco_collapsed),
              alpha = 0.1) +
  geom_point(data = nwa_dat, 
             aes(y = focal_std_by_ecoregion, 
                 color = eco_collapsed),
             alpha = 0.2) +
  ylim(c(0,1.1)) +
  labs(y = "Standardized Kelp Abundance", x = "",
       color = "", fill = "") +
  guides(color=guide_legend(nrow=2,
                            byrow=TRUE))+
  theme_bw(base_size = 14)+
  theme(legend.position = "bottom",
        legend.box="vertical") +
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2") +
  annotate(x = 1960, y = .78, 
           geom = "label",
           label.size = 0,
           size = 4.5,
           label = paste0(perc_change_per_year,"% change per year"))

ggsave("figures/timeseries_with_curves.jpg",
       width = 7, height = 6)

# ordbeta plot?
modelbased::estimate_relation(mod_ecoregion_ord,
                                length = 100,
                                by = c("year", "eco_collapsed")) |>
  plot(show_data = TRUE)+
  ylim(c(0,1.1)) +
  labs(y = "Standardized Kelp Abundance", x = "",
       color = "", fill = "") +
  guides(color=guide_legend(nrow=2,
                            byrow=TRUE))+
  theme_bw(base_size = 14)+
  theme(legend.position = "bottom",
        legend.box="vertical") +
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")

ggsave("figures/timeseries_with_curves_ordbeta.jpg",
       width = 7, height = 6)
# show with RE

modelbased::estimate_relation(mod_common_slope,
                              include_random = TRUE,
                              length = 200,
                              by = c("year", "trajectory")) |>
  plot(show_data = TRUE, 
       color = "grey",
       ribbon = list(alpha = 0)) +
  labs(y = "Standardized Kelp Abundance", x = "",
       color = "", fill = "") +
  guides(color = "none", fill = "none") +
  theme_bw(base_size = 14)+
  theme(legend.position = "bottom",
        legend.box="vertical")

ggsave("figures/timeseries_all.jpg",
       width = 7, height = 6)

##
# Decadal change?
##
nwa_dat <- nwa_dat |>
  mutate(decade = (floor(year/10)*10) |> 
           paste0("s") |>
           as.factor(),
         decade = fct_collapse(decade,
                               "pre-1980s" = c(
                                 "1970s",
                                 "1960s",
                                 "1940s"
                               )),
         decade = fct_collapse(decade,
                               "post-2010" = c(
                                 "2010s",
                                 "2020s"
                               )))

mod_decadal <- glmmTMB(ln_focal_std_by_ecoregion ~ 
                              eco_collapsed+year*decade + 
                              (1 + year |trajectory) + (1|study), 
                            dispformula =~focalUnit,
                            data = nwa_dat )

Anova(mod_decadal)

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

# 
# parameters::model_parameters(mod_decadal,
#                              effects = "fixed",
#                              component = "conditional",
#                              include_sigma = FALSE,
#                              drop = "(eco_collapsed)",
#                              ci = 0.8) |> plot() +
#   scale_color_manual(values = c("black", "black")) +
#   labs(x = "% change per 1 unit change in x")
##
# Percent Cover Only with Ordbetareg
##
perc_mod <- glmmTMB(rescaled_std_by_ecoregion ~ 
                      ecoregion*year_c + 
                      (1  + year_c|trajectory) + (1|study), 
                    family = ordbeta,
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
                       mean(nwa_dat$year) |> round())

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

