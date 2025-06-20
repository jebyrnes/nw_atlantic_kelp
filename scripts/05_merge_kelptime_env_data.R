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
## remember to start with oldest first to get all years
## which would be hadsst
get_had_first <- function(file_vec){
  idx <- which(grepl("hadsst_seasonal", file_vec))
  file_vec[c(idx, c(1:length(file_vec))[-idx])]
}

env_dat <- list.files("data/env_data",
                      full.names = TRUE) |>
  get_had_first() |>
  ## strip off cell, x, and y
  map(~read_csv(.x) |> select(-c(x, y, cell))) |>
  ## merge them together by accumulating left joins
  reduce(left_join)

# checks
visdat::vis_dat(env_dat)

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
           x = hadsst_summer,
           color = as.character(trajectory))) +
  geom_point() +
  scale_color_discrete(guide = "none")

write_csv(nwa_with_env, "data/nwa_with_env.csv")
