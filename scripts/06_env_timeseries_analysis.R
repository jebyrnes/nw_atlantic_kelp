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
library(sf)
library(blme)

theme_set(theme_bw(base_size = 12))

env_dat <- read_csv("data/merged_envt_data_all.csv") |>
  mutate(id = paste(latitude, longitude) |> as.factor() |> as.numeric() |> as.character(),
         year_c = year - mean(year),
         year_z = year_c/sd(year),
         year_f = as.factor(year))

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

hadsst_summer_mod <- glmmTMB(
  hadsst_summer ~ 
    year  + (year|id),
  data = env_dat |> 
    get_unique_cells(hadsst_seasonal_cell_number) # subset so not huge
)

# same as initial ggplot
estimate_relation(hadsst_summer_mod, 
                  include_random = TRUE,
                  ci = NA) |>
  plot(show_data = TRUE,
       line = list(size = 1)) +
  guides(color = "none", fill = "none") +
  scale_color_viridis_d(option = "F", direction = -1) 

hadsst_spring_mod <- glmmTMB(
  hadsst_spring ~ 
    year + (1+year|id),
  data = env_dat |>
    get_unique_cells(hadsst_seasonal_cell_number) # subset so not huge
)

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


oisst_summer_mod <- glmmTMB(
  oisst_summer ~ 
    year + (1+year|id),
  data = env_dat |>
    get_unique_cells(oisst_seasonal_cell_number) # subset so not huge
)

oisst_spring_mod <- glmmTMB(
  oisst_spring ~ 
    year + (1+year|id),
  data = env_dat |>
    get_unique_cells(oisst_seasonal_cell_number) # subset so not huge
)

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


swh_fall_mod <- glmmTMB(
  swh_fall ~ 
    year + (1 |id),
  
  data = env_dat |>
    get_unique_cells(era5_waves_seasonal_cell_number) # subset so not huge
)

# swh in winter
swh_winter_mod <- glmmTMB(
  swh_winter ~ 
    year + (1 |id),
  data = env_dat |>
    get_unique_cells(era5_waves_seasonal_cell_number) # subset so not huge
)

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



no3_spring_mod <- glmmTMB(
  no3_spring ~ 
    year_c + (1 + year_c| id),
  data = env_dat |>
    get_unique_cells(cmems_seasonal_cell_number) # subset so not huge
)

estimate_relation(no3_spring_mod, include_random = TRUE) |>
  plot(show_data = TRUE) + 
  guides(color = "none", fill = "none")

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


no3_winter_mod <- glmmTMB(
  no3_winter ~ 
    year_c + (1+year_c|id),
  data = env_dat |>
    get_unique_cells(cmems_seasonal_cell_number) # subset so not huge
)


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


kd490_spring_mod <- glmmTMB(
  kd490_spring ~ 
    year + (1 |id),
  REML = TRUE,
  data = env_dat |>
    get_unique_cells(turb_seasonal_cell_number) # subset so not huge
)

###
# Put it all together
###
get_year_coef <- function(mod, coef = "year", name = NA){
  tidy(mod) |>
    filter(term == coef) |>
    select(estimate, std.error, p.value) |>
    mutate(model = name) |>
    relocate(model, estimate, std.error, p.value) |>
    rename(
      ` ` = model,
      coefficient = estimate,
      se = std.error,
      p = p.value
    )
}

envt_tab <- bind_rows(
  get_year_coef(hadsst_summer_mod, name = "HADSST summer"),
  get_year_coef(hadsst_spring_mod, name = "HADSST spring"),
  get_year_coef(oisst_summer_mod, name = "OISST summer"),
  get_year_coef(oisst_spring_mod, name = "OISST spring"),
  get_year_coef(swh_fall_mod, name = "SWH fall"),
  get_year_coef(swh_winter_mod, name = "SWH winter"),
  get_year_coef(no3_winter_mod, name = "NO3 winter", coef = "year_c"),
  get_year_coef(no3_spring_mod, name = "NO3 spring", coef = "year_c"),
  get_year_coef(kd490_spring_mod, name = "Turbidity spring")
)

write_csv(envt_tab, "tables/envt_coefs.csv")

gt::gt(envt_tab) |> 
  gt::fmt_number(decimals=3) |>
  gt::gtsave("tables/envt_coefs_gt.docx")
