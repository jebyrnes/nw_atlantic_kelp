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

theme_set(theme_bw(base_size = 12))

nwa <- read.csv("data/kelptime_nwa_data.csv") |> as_tibble() |>
  mutate(date = ymd(sasdate),
         year = year(date),
         year_f = numFactor(year),
         year_c = year - mean(year),
         trajectory = paste(study, site)) |>
  filter(!is.na(ln_focal_std_by_ecoregion)) |>
  group_by(trajectory) |>
  mutate(n_per_trajectory = n()) |>
  ungroup() |>
  filter(n_per_trajectory >= 2) # at least 3 data points per trajectory

# how many? And how have we improved?
nwa |>
  pull(trajectory) |>
  unique() |>
  length()

nwa_old <- read.csv("https://github.com/kelpecosystems/global_kelp_time_series/raw/refs/heads/master/05_HLM_analysis_code/formatted_data_3years.csv")

nwa_old |>
  filter(PROVINCE == ("Cold Temperate Northwest Atlantic")) |>
  filter(!is.na(ln_stdByECOREGION)) |>
  pull(trajectory_ID) |>
  unique() |>
  length()

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
mod <- glmmTMB(ln_focal_std_by_ecoregion ~ 0 + ecoregion*year + 
                 (1  |trajectory) ,#+ 
                # ou(0 + year_f | trajectory), 
               dispformula =~focalUnit,
               #using numfactor(year) for OU as a car structure
               data = nwa)

mod_no_eco <- glmmTMB(ln_focal_std_by_ecoregion ~ 0 + year + 
                 (1  |trajectory),# + 
              #   ou(0 + year_f | trajectory), 
               dispformula =~focalUnit,
               #using numfactor(year) for OU as a car structure
               data = nwa)

# check model
performance::check_model(mod)

# What's the model tell us?
AIC(mod_no_eco)
AIC(mod) # it's this one

car::Anova(mod)
tidy(mod, effects = "fixed")

# Plot slopes
slopes <- emtrends(mod, specs = ~ecoregion, var = "year")

plot(slopes)
