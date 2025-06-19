#' -------------------------------------------
#' Merge env data sets and then merge
#' with the kelptime nwa data 
#' -------------------------------------------

library(dplyr)
library(purrr)
library(tidyr)
library(lubridate)

###
# envt merge
###

## Load all envt data sets
env_dat <- list.files("data/env_data",
                      full.names = TRUE) |>
  rev() |>
  ## strip off cell, x, and y
  map(~read_csv(.x) |> select(-c(x, y, cell))) |>
  ## merge them together by accumulating left joins
  reduce(left_join)

# checks
#sapply(env_dat, class)
# GGally::ggpairs(env_dat |> 
#                   select(contains("cellmean")) |>
#                   select(contains("summer")))



###
# merge envt and kelptime
###

## load kelptime
nwa_kelp <- read_csv("data/kelptime_nwa_data.csv")

## one big left join!
nwa_with_env <- left_join(nwa_kelp, env_dat)

## check and be done
visdat::vis_dat(nwa_with_env)

ggplot(nwa_with_env,
       aes(y = biomass_kg_wet_per_sq_m,
           x = oisst_summer,
           color = as.character(trajectory))) +
  geom_point() +
  scale_color_discrete(guide = "none")

write_csv(nwa_with_env, "data/nwa_with_env.csv")
