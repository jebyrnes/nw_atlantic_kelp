library(terra)

extract_with_nas <- function(a_rast_stack, 
                             lat_long_data,
                             search_radius = 3*55000){
  
  # na strategy from https://github.com/rspatial/terra/issues/1733
  # create e to deal with pixels too close to land
  # gets nearest non-NA cell
  e <- extract(a_rast_stack[[1]], 
               lat_long_data,
               xy=TRUE, 
               search_radius = search_radius,
               cells = TRUE, raw = TRUE) |> as_tibble()
  
  dat_extracted <- terra::extract(a_rast_stack, 
                                    e$cell,
                                    xy = TRUE) |> 
    as_tibble() |>
    mutate(cell = e$cell,
           longitude = lat_long_data$longitude,
           latitude = lat_long_data$latitude,
    )
  
  return(dat_extracted)
}

make_monthly <- function(dat_rast, 
                         lat_long_dat,
                         varname,
                         search_radius = 3*55000){
  names(dat_rast) <- time(dat_rast)
  
  extract_with_nas(dat_rast,
                   lat_long_dat,
                   search_radius = search_radius) |>
    pivot_longer(-c(cell,x,y, latitude, longitude),
                 values_to = varname,
                 names_to = "date") |>
    mutate(date = ymd(date)) |>
    make_date_to_year_mo_season()  |>
    as_tibble()
  
}

make_date_to_year_mo_season <- function(dat){
  dat |>
    mutate(
      month = month(date),
      year = year(date),
      season = quarter(date, type = "year.quarter", fiscal_start = 12),
      season_name = case_when(
        grepl("\\.1", season) ~ "winter",
        grepl("\\.2", season) ~ "spring",
        grepl("\\.3", season) ~ "summer",
        grepl("\\.4", season) ~ "fall",
        .default = "other"
      )
    ) 
}


make_seasonal <- function(monthly_dat, varname){
  monthly_dat |>
    group_by(x, y, cell, longitude, latitude,
             season, season_name) |>
    summarize(across({{varname}}, mean)) |>
    ungroup() |>
    mutate(year = gsub("\\.[1-4]", "", season)) |>
    select(-season) |>
    pivot_wider(names_from = season_name, 
                values_from = {{varname}},
                names_prefix = glue("{varname}_"))
}
