#' -----------------------------------
#' download environmental predictors for data
#' and merge
#' -----------------------------------
source("scripts/add_env_driver_functions.R")
source("scripts/wf_keyset.R")
has_env_data <- TRUE #set to false if first time running or for new data

nwa_data <- read.csv("data/kelptime_nwa_data.csv") |> as_tibble() |>
  mutate(date = ymd(sasdate),
         year = year(date),
         year_c = year - mean(year),
         month = month(date),
         trajectory = paste(study, site))

# download the environmental data - only run if needed
if(!has_env_data){
  no3_spring(nwa_data)
  turb_spring(nwa_data)
  SST_spring(nwa_data)
  SST_winter(nwa_data)
  no3_decadal(nwa_data)

  WI_fall_winter(nwa_data)
  WI_fall_winter_max(nwa_data)
}
 
## Add to the dataset
nitrates <- read_csv("data/env_data/no3_spring.csv")|>
  distinct(y, latitude, longitude, .keep_all=T)
turbidity <- read_csv("data/env_data/turbidity_spring.csv")|>
  distinct(y, latitude, longitude, .keep_all=T)|>
  mutate(year = as.numeric(y))
SST <- read_csv("data/env_data/sst_spring.csv")|>
  distinct(y, latitude, longitude, .keep_all=T)
WI <- read_csv("data/env_data/wi_fall_winter.csv")|>
  distinct(y, latitude, longitude, .keep_all=T)


nwa_data_env <- nwa_data |>
  left_join(nitrates, by=join_by(year==y, latitude==latitude, longitude==longitude))|> #still some missing data
  left_join(WI, by=join_by(year==y, latitude==latitude, longitude==longitude))|>
  left_join(turbidity, by=join_by(year==year, latitude==latitude, longitude==longitude))|>
  left_join(SST, by=join_by(year==y, latitude==latitude, longitude==longitude))


write_csv(nwa_data_env, "data/nwa_with_env.csv")
