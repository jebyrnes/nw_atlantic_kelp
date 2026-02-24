#' -----------------------------------------------
#' NWA Timeseries Analysis and Viz
#' in the style of the previous Krumhansl et al 2016
#' analysis where trend can vary by ecoregion
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
           group = trajectory, color = eco_collapsed)) +
  geom_point(alpha = 0.3) +
 # facet_wrap(vars(eco_collapsed)) +
  scale_color_discrete(guide = "none") +
  stat_smooth(method = "lm", fill = NA, 
              size = 0.5, alpha = 0.5) +
  ylim(c(0,1.1)) +
  ylim(c(0,1.1)) +
  labs(y = "Kelp Abundance Standardized\nBy Ecoregion", x = "",
       color = "") +
  guides(color=guide_legend(nrow=2,
                            byrow=TRUE))+
  theme_bw(base_size = 14)+
  theme(legend.position = "bottom",
        legend.box="vertical") +
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")


ggsave("figures/timeseries_raw_linear.jpg",
       width = 7, height = 6)

##
# Model by ecoregion
##
mod_ecoregion <- glmmTMB(ln_focal_std_by_ecoregion ~ 
                  eco_collapsed*year + 
                 (1 + year |trajectory) + (1|study), 
               dispformula =~focalUnit,
               data = nwa_dat)

saveRDS(mod_ecoregion, "models/mod_ecoregion.rds")

# check model
performance::check_convergence(mod_ecoregion)

performance::check_predictions(mod_ecoregion)

performance::check_residuals(mod_ecoregion) |> plot()


# What's the model tell us?

car::Anova(mod_ecoregion)

change_est <- 
  emtrends(mod_ecoregion, ~1, var = "year") |> tidy()|> pull(year.trend)

perc_change_per_year <- ((exp(change_est) -1 )*100) |> round(2)

modelbased::estimate_grouplevel(mod_ecoregion,
                                type = "total",
                                group = "trajectory:year") |>
  filter(!grepl("Intercept", Parameter))|> 
  mutate(Level = forcats::fct_reorder(Coefficient |> as.factor(), 
                                      Coefficient)) |>
  plot() +
  guides(y = "none") +
  geom_hline(yintercept = 0, lty = 2, color = "red")

tidy(mod_ecoregion, effects = "fixed")

# Plot slopes
slopes <- emtrends(mod_ecoregion, specs = ~eco_collapsed, var = "year")

plot(slopes) +
  labs(y = "", x = "standardized change")

estimate_slopes(mod_ecoregion, trend = "year", 
                by = "eco_collapsed",
                backend = "emmeans") |>
  plot()

ggsave("figures/ecoregion_slopes.jpg", 
       width = 7, height = 6)

contrast(slopes,
         method = "pairwise") 

# Show the model results

plot_mod <- function(a_mod){
  estimate_relation(a_mod,
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
    scale_fill_brewer(palette = "Dark2")
}

plot_mod(mod_ecoregion)

ggsave("figures/timeseries_ecoint_with_curves.jpg",
       width = 7, height = 6)


estimate_relation(mod_ecoregion,
                  include_random = TRUE,
                  length = 200,
                  data = NULL)|>
  as_tibble() |>
  mutate(Predicted = exp(Predicted)-0.01,
         CI_low = exp(CI_low)-0.01,
         CI_high = exp(CI_high)-0.01
  ) |>
  ggplot(aes(x = year)) +
  geom_line(aes(y = Predicted, 
                color = eco_collapsed,
                group = trajectory),
            size = 1.1,
            alpha = 0.5) +
  labs(y = "Standardized Kelp Abundance", x = "",
       color = "", fill = "") +
  guides(color=guide_legend(nrow=2,
                            byrow=TRUE))+
  theme_bw(base_size = 14)+
  theme(legend.position = "bottom",
        legend.box="vertical") +
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")

ggsave("figures/timeseries_all_modeled.jpg",
       width = 7, height = 6)
