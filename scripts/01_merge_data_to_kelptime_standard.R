#' -----------------------------------------------
#' Merge datasets into Kelptime format
#' Author: Jarrett Byrnes jarrett.byrnes@umb.edu
#' Dataset Last Updated:  
#' Script Last Updated: June 8, 2023
#'
#' -----------------------------------------------

# load helpers
library(purrr)
source("scripts/data_checking_functions.R")
source("scripts/combine_data_functions.R")
source("scripts/add_geo_info.R")

# run the data check scripts
raw_data_files <- list.files("data/clean_data/timeseries/",
                             pattern = "csv")

map(raw_data_files, read_kelp_data) |>
  walk(check_problems)


# combine in old data as well
#moved GOM files to new clean data

# combine the data files - new and old
new_data <- 
  map_df(raw_data_files, process_one_study_to_combine)

# create combined data with geospatial info
# then get focal kelp standardized by each region
log_add_one <- function(x) log(x+1)
log_add_01 <- function(x) log(x+.01)

combined_data <- new_data |>
  filter(!is.na(study)) |>
  filter(!is.na(latitude)) |>
  add_geo_info() |>
  
  group_by(ecoregion, focalUnit) |>
  mutate(focal_std_by_ecoregion = standardize_by_max(focalKelp))|>
  ungroup() |>
  
  group_by(province, focalUnit) |>
  mutate(focal_std_by_province = standardize_by_max(focalKelp))|>
  ungroup() |>
  
  group_by(realm, focalUnit) |>
  mutate(focal_std_by_realm = standardize_by_max(focalKelp)) |>
  ungroup() |>
  
  group_by(focalUnit) |>
  mutate(focal_std_by_all = standardize_by_max(focalKelp)) |>
  ungroup() |>
  
  mutate(across(.cols = focal_std_by_ecoregion:focal_std_by_all,
                .fns = log_add_01,
                .names = "ln_{.col}")
  )

# filter out some duplicate studies 
# some of which we now have raw data for
combined_data <- combined_data |>
  filter(!(study %in% c("Attridge_et_al._2022", 
                        "Feehan_et_al._2019",
                        "Filbee-Dexteretal_et_al._2016",
                        "Egan_&_Yarish__1990"
                        )))

#write out data
write_csv(combined_data, "data/kelptime_nwa_all_data.csv")


#write out data that meets time requirements for analysis
combined_clear <- combined_data |>
  as_tibble() |>
  mutate(date = ymd(sasdate),
         month = month(date),
         year = year(date),
         year_f = as.character(year),
         year_c = year - mean(year),
         trajectory = paste(study, site)) |>
  filter(!is.na(ln_focal_std_by_ecoregion)) |>
  group_by(trajectory) |>
  mutate(n_per_trajectory = n()) |>
  ungroup() |>
  filter(n_per_trajectory >= 3) 

write_csv(combined_clear, 
           "data/kelptime_nwa_data.csv")

# unique lat/longs and date ranges
combined_clear |>
  group_by(latitude, longitude) |>
  summarize(min_year = min(year),
            max_year = max(year),
            .groups = "drop") |>
  write_csv("data/unique_latlongs_time.csv")

