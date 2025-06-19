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

combined_data <- new_data |>
  filter(!is.na(study)) |>
  add_geo_info() |>
  
  group_by(ecoregion, focalUnit) |>
  mutate(focal_std_by_ecoregion = standardize_by_max(focalKelp))|>
  ungroup() |>
  
  group_by(province) |>
  mutate(focal_std_by_province = standardize_by_max(focalKelp))|>
  ungroup() |>
  
  group_by(realm) |>
  mutate(focal_std_by_realm = standardize_by_max(focalKelp)) |>
  ungroup() |>
  
  mutate(across(.cols = focal_std_by_ecoregion:focal_std_by_realm,
                .fns = log_add_one,
                .names = "ln_{.col}")
  )


#write out data
write_csv(combined_data, "data/kelptime_nwa_all_data.csv")


#write out data that meets time requirements for analysis
combined_clear <- combined_data |>
  as_tibble() |>
  mutate(date = ymd(sasdate),
         year = year(date),
         year_f = as.character(year),
         year_c = year - mean(year),
         trajectory = paste(study, site)) |>
  filter(!is.na(ln_focal_std_by_ecoregion)) |>
  group_by(trajectory) |>
  mutate(n_per_trajectory = n()) |>
  ungroup() |>
  filter(n_per_trajectory >= 2)# at least 3 data points per trajectory

write_csv(combined_clear, 
           "data/kelptime_nwa_data.csv")

# unique lat/longs and date ranges
combined_clear |>
  group_by(latitude, longitude) |>
  summarize(min_year = min(year),
            max_year = max(year)) |>
  ungroup() |>
  write_csv("data/unique_latlongs_time.csv")

# plot timeseries
library(ggplot2)
combined_clear <- combined_clear |>
  mutate(trj = as.character(trajectory) |>
           forcats::fct_reorder(latitude, .desc = FALSE))

ggplot(combined_clear,
       aes(x = year, 
           y = latitude,#trj |> as.numeric(),
           group = trj,
           color = ecoregion)) +
  geom_line(alpha = 1) +
  scale_color_discrete(guide = "none") +
  # scale_y_continuous(breaks = seq(0, 
  #                                 max(combined_clear$trj|>as.numeric()), 
  #                                 length.out=4),
  #                    labels = seq(min(combined_clear$latitude), 
  #                                 max(combined_clear$latitude), 
  #                                 length.out=4) |> round(2)) +
   labs(y = "Latitude", x = "") +
  theme_bw()


