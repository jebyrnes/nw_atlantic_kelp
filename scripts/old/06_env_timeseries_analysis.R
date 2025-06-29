#' -----------------------------------------------
#' Analyze trends in drivers at sites
#'
#' -----------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(performance)
library(parameters)
library(modelbased)
library(glmmTMB)
library(broom.mixed)
library(car)
library(mregions)
library(sf)

theme_set(theme_bw(base_size = 12))

env_dat <- read_csv("data/merged_envt_data_all.csv") |>
  mutate(id = paste(latitude, longitude) |> as.factor() |> as.numeric() |> as.character(),
         year_c = year - mean(year),
         year_z = year_c/sd(year),
         year_f = as.factor(year))


ecoregions_shp <- mr_shp(key = "Ecoregions:ecoregions") |> 
  st_make_valid() |>
  dplyr::select(ecoregion, geometry)


add_geo_info <- function(dat){
  dat |>
    st_as_sf(coords = c("longitude", "latitude"),
             crs = st_crs(ecoregions_shp),
             remove = FALSE)|>
    st_join(ecoregions_shp, 
            join = st_nearest_feature)
}

env_dat <- env_dat |>
  add_geo_info()

get_unique_cells <- function(dat, cell_var){
  dat |>
    group_by(year, {{cell_var}}) |>
    slice(1L) |>
    ungroup() |>
    mutate(id = as.factor({{cell_var}}))
}

##
# hadsst
##

## Plots
ggplot(env_dat |> get_unique_cells(hadsst_seasonal_cell_number),
       aes(x = year, y = hadsst_summer, 
           group = id, color = id)) +
  geom_line() +
  scale_color_viridis_d(option = "F", direction = -1) +
  guides(color = "none") +
  labs(x = "", y = "HADSST ºC", subtitle = "Summer Average")

ggplot(env_dat|> get_unique_cells(hadsst_seasonal_cell_number),
       aes(x = year, y = hadsst_spring, 
           group = id, color = id)) +
  geom_line(alpha = 0.01) +
  geom_smooth(method = "lm", fill = NA) +
  scale_color_viridis_d(option = "F", direction = -1) +
  guides(color = "none")+
  labs(x = "", y = "HADSST ºC", subtitle = "Spring Average")


hadsst_summer_lm <- glmmTMB(
  hadsst_summer ~ 
    year  + (year|id),
  data = env_dat |> 
    get_unique_cells(hadsst_seasonal_cell_number) # subset so not huge
)

#had_slopes <- estimate_slopes(hadsst_lm, trend = "year")
emmeans::emtrends(hadsst_summer_lm,  var = "year") |> 
  as_tibble() |>
ggplot(
       aes(x = year.trend)) +
  geom_histogram(fill = "pink") +
  geom_vline(xintercept = 0, lty = 2, color = "red")

# same as initial ggplot
estimate_relation(hadsst_summer_lm, ci = NA) |>
  plot(show_data = TRUE,
       line = list(size = 1)) +
  guides(color = "none", fill = "none") +
  scale_color_viridis_d(option = "F", direction = -1) 
  
hadsst_spring_lm <- lm(
  hadsst_spring ~ 
    year * id,
  data = env_dat |>
    get_unique_cells(hadsst_seasonal_cell_number) # subset so not huge
)

# dist of slopes
emmeans::emtrends(hadsst_spring_lm, ~id, var = "year") |> 
  as_tibble() |>
ggplot(
       aes(x = year.trend)) +
  geom_histogram(fill = "pink") +
  geom_vline(xintercept = 0, lty = 2, color = "red")


# Smooth to deal with autocorrelation
  
hadsst_gam_mod <- mgcv::gam(
  hadsst_summer ~ 
    year * id +
    s(year, by = id), 
  data = env_dat |>
    get_unique_cells(hadsst_seasonal_cell_number)
)

estimate_relation(hadsst_gam_mod, ci = NA) |> 
  plot(show_data = TRUE) +
  guides(color = "none", fill = "none") +
  scale_color_viridis_d(option = "F", direction = -1) 

##
# oisst
##

## Plots
ggplot(env_dat|> get_unique_cells(oisst_seasonal_cell_number) |>
         filter(!is.na(oisst_summer)),
       aes(x = year, y = oisst_summer, 
           group = id, color = id)) +
  geom_line(alpha = 0.1) +
  geom_smooth(method = "lm", fill = NA, size = 0.4) +
  scale_color_viridis_d(option = "F", direction = -1) +
  guides(color = "none")+
  labs(x = "", y = "OISST ºC", subtitle = "Summer Average")

ggplot(env_dat|> get_unique_cells(oisst_seasonal_cell_number) |>
         filter(!is.na(oisst_spring)),
       aes(x = year, y = oisst_spring, 
           group = id, color = id)) +
  geom_line(alpha = 0.1) +
  geom_smooth(method = "lm", fill = NA, size = 0.4) +
  scale_color_viridis_d(option = "F", direction = -1) +
  guides(color = "none")+
  labs(x = "", y = "OISST ºC", subtitle = "Spring Average")

##
# swh
##

ggplot(env_dat|> get_unique_cells(era5_waves_seasonal_cell_number) |>
         filter(!is.na(swh_fall)),
       aes(x = year, y = swh_fall, 
           group = id, color = id)) +
  geom_line(alpha = 0.1) +
  geom_smooth(method = "lm", fill = NA, size = 0.5) +
  scale_color_viridis_d(option = "H", direction = -1) +
  guides(color = "none")+
  labs(x = "", y = "Significant Wave Height (m)", subtitle = "Fall Average")

ggplot(env_dat|> get_unique_cells(era5_waves_seasonal_cell_number) |>
         filter(!is.na(swh_winter)),
       aes(x = year, y = swh_winter, 
           group = id, color = id)) +
  geom_line(alpha = 0.1) +
  geom_smooth(method = "lm", fill = NA, size = 0.5) +
  scale_color_viridis_d(option = "H", direction = -1) +
  guides(color = "none")+
  labs(x = "", y = "Significant Wave Height (m)", subtitle = "Winter Average")



swh_fall_lm <- lm(
  swh_fall ~ 
    year * id,
  data = env_dat |>
    get_unique_cells(era5_waves_seasonal_cell_number) # subset so not huge
)

# Analysis
Anova(swh_fall_lm)
tidy(swh_fall_lm)[1:5,]
emmeans::emtrends(swh_fall_lm, ~1, "year")

# dist of slopes
emmeans::emtrends(swh_fall_lm, ~id, var = "year") |> 
  as_tibble() |>
  ggplot(
    aes(x = year.trend)) +
  geom_histogram(fill = "pink") +
  geom_vline(xintercept = 0, lty = 2, color = "red")

# swh in winter
swh_winter_lm <- lm(
  swh_winter ~ 
    year * id,
  data = env_dat |>
    get_unique_cells(era5_waves_seasonal_cell_number) # subset so not huge
)

# Analysis
Anova(swh_winter_lm)
tidy(swh_winter_lm)[1:5,]
emmeans::emtrends(swh_winter_lm, ~1, "year")

# dist of slopes
emmeans::emtrends(swh_winter_lm, ~id, var = "year") |> 
  as_tibble() |>
  ggplot(
    aes(x = year.trend)) +
  geom_histogram(fill = "pink") +
  geom_vline(xintercept = 0, lty = 2, color = "red")

##
# no3
##

ggplot(env_dat|> get_unique_cells(cmems_seasonal_cell_number) |>
         filter(!is.na(no3_spring)),
       aes(x = year, y = no3_spring, 
           group = id, color = id)) +
  geom_line(alpha = 0.5) +
  #geom_smooth(method = "lm", fill = NA, size = 0.5) +
  scale_color_viridis_d(option = "H", direction = -1) +
  guides(color = "none")+
  labs(x = "", y = "Nitrate", subtitle = "Spring Average")



no3_spring_lm <- lm(
  no3_spring ~ 
    year * id,
  data = env_dat |>
    get_unique_cells(cmems_seasonal_cell_number) # subset so not huge
)

# Analysis
Anova(no3_spring_lm)
emmeans::emtrends(no3_spring_lm, ~1, "year")

# dist of slopes
emmeans::emtrends(no3_spring_lm, ~id, var = "year") |> 
  as_tibble() |>
  ggplot(
    aes(x = year.trend)) +
  geom_histogram(fill = "pink", bins = 100) +
  geom_vline(xintercept = 0, lty = 2, color = "red",
             size = 0.3) 

## winter

ggplot(env_dat|> get_unique_cells(cmems_seasonal_cell_number) |>
         filter(!is.na(no3_winter)),
       aes(x = year, y = no3_winter, 
           group = id, color = id)) +
  geom_line(alpha = 0.5) +
  #geom_smooth(method = "lm", fill = NA, size = 0.5) +
  scale_color_viridis_d(option = "H", direction = -1) +
  guides(color = "none")+
  labs(x = "", y = "Nitrate", subtitle = "Winter Average")


no3_winter_lm <- lm(
  no3_winter ~ 
    year * id,
  data = env_dat |>
    get_unique_cells(cmems_seasonal_cell_number) # subset so not huge
)

# Analysis
Anova(no3_winter_lm)
emmeans::emtrends(no3_winter_lm, ~1, "year")

# dist of slopes
emmeans::emtrends(no3_winter_lm, ~id, var = "year") |> 
  as_tibble() |>
  ggplot(
    aes(x = year.trend)) +
  geom_histogram(fill = "pink", bins = 100) +
  geom_vline(xintercept = 0, lty = 2, color = "red",
             size = 0.3) 

##
# kd490
##

ggplot(env_dat|> get_unique_cells(turb_seasonal_cell_number) |>
         filter(!is.na(kd490_spring)),
       aes(x = year, y = kd490_spring, 
           group = id, color = id)) +
  geom_line(alpha = 0.5) +
  #geom_smooth(method = "lm", fill = NA, size = 0.5) +
  scale_color_viridis_d(option = "H", direction = -1) +
  guides(color = "none")+
  labs(x = "", y = "Turbidity (kd490)", subtitle = "Spring Average")


kd490_spring_lm <- lm(
  kd490_spring ~ 
    year * id,
  data = env_dat |>
    get_unique_cells(turb_seasonal_cell_number) # subset so not huge
)

# Analysis
Anova(kd490_spring_lm)
emmeans::emtrends(kd490_spring_lm, ~1, "year")

# dist of slopes
emmeans::emtrends(kd490_spring_lm, ~id, var = "year") |> 
  as_tibble() |>
  ggplot(
    aes(x = year.trend)) +
  geom_histogram(fill = "pink") +
  geom_vline(xintercept = 0, lty = 2, color = "red")
