#' -------------------------------------------
#' Merge env data sets and then merge
#' with the kelptime nwa data 
#' -------------------------------------------

library(dplyr)
library(tidyr)
library(lubridate)
env_dat <- read_csv("data/merged_envt_data_all.csv")

###
# merge envt and kelptime
###

## load kelptime
nwa_kelp <- read_csv("data/kelptime_nwa_data.csv")

## one big left join!
nwa_with_env <- left_join(nwa_kelp, env_dat)

write_csv(nwa_with_env, "data/nwa_with_env.csv")


## check and be done
visdat::vis_dat(nwa_with_env)

ggplot(nwa_with_env,
       aes(y = percent_cover,
           x = swh_winter,
           color = as.character(trajectory))) +
  geom_point() +
  scale_color_discrete(guide = "none")

