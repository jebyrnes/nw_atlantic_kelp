source("scripts/data_checking_functions.R")
library(lubridate)

mean_na <- function(x) mean(x, na.rm = TRUE)
max_na <- function(x) max(x, na.rm = TRUE)

process_one_study_to_combine <- function(datasetname, path = "data/clean_data/timeseries/"){
  print(glue("Working on {datasetname}"))
  
  # read in data
  kelp_data <- read_kelp_data(datasetname, path = path) 
  
  # reprocess
  ret <- kelp_data |>
    
    # for NWA, if no stipes or individuals, use the other
    # as stipes = individuals
    mutate(
      Stipe_Density_num_per_sq_m = ifelse(is.na(Stipe_Density_num_per_sq_m), 
                                          Individual_Density_num_per_sq_m,
                                          Stipe_Density_num_per_sq_m),
      Individual_Density_num_per_sq_m = ifelse(is.na(Individual_Density_num_per_sq_m),
                                               Stipe_Density_num_per_sq_m,
                                               Individual_Density_num_per_sq_m)
      
    ) |>
    
    mutate(trajectory_ID = paste(Study, Site) |> as.factor() |> as.numeric(),
           kelpPresent = ifelse(
             Percent_Cover > 0 | 
               Individual_Density_num_per_sq_m > 0 | 
               Biomass_kg_wet_per_sq_m > 0 | 
               Stipe_Density_num_per_sq_m > 0,
             1,
             NA_real_
           ),
           sasDate = paste(Sample_Year, Sample_Month, Sample_Day, sep = "-") |> ymd(),
           rawDate = interval( ymd("19600101"), ymd(sasDate)) %/% days(1)) |> #days since jan 1 1960
    arrange(Study, Site, sasDate) |>
    group_by(Study, Site) |>
    
    # add sample period
    mutate(
      period = sasDate |> as.factor() |> as.numeric()
    ) |>
    
    # Average each taxon within a sample_ID, Month, and Year (used to be period)
    # This is to ensure we don't get duplicate measurements, or if the same
    # permanent plot was sampled multiple times in the same month which happens
    # in some surveys
    group_by(Study, Site, Sample_Year, Sample_Month, Sample_ID, Taxon) |>
    summarize(
      Latitude = mean_na(Latitude),
      Longitude = mean_na(Longitude),
      Sample_Unit_Size_sq_m = mean_na(Sample_Unit_Size_sq_m),
      kelpPresent = max(kelpPresent),
      max_Depth = max_na(Depth_m),
      mean_Depth = mean_na(Depth_m),
      min_Depth = min(Depth_m, na.rm = TRUE),
      sasDate = mean_na(sasDate),
      Biomass_kg_wet_per_sq_m = mean_na(Biomass_kg_wet_per_sq_m),
      Stipe_Density_num_per_sq_m = mean_na(Stipe_Density_num_per_sq_m),
      Individual_Density_num_per_sq_m = mean_na(Individual_Density_num_per_sq_m),
      Percent_Cover = mean_na(Percent_Cover),
      .groups = "drop"
    ) |>
    
  # Sum across taxa for each sample and Sample.Year, Sample.Month - used to use period
  # in individual sample units
  group_by(Study, Site, Sample_Year, Sample_Month, Sample_ID) |>
    summarize(
      Latitude = mean_na(Latitude),
      Longitude = mean_na(Longitude),
      Sample_Unit_Size_sq_m = mean_na(Sample_Unit_Size_sq_m),
      kelpPresent = max(kelpPresent),
      max_Depth = max(max_Depth),
      mean_Depth = mean_na(mean_Depth),
      min_Depth = min(min_Depth),
      sasDate = mean_na(sasDate),
      Biomass_kg_wet_per_sq_m = sum(Biomass_kg_wet_per_sq_m),
      Stipe_Density_num_per_sq_m = sum(Stipe_Density_num_per_sq_m),
      Individual_Density_num_per_sq_m = sum(Individual_Density_num_per_sq_m),
      Percent_Cover = sum(Percent_Cover),
      .groups = "drop"
    ) |>
    
    # aggregate measurements by period - so getting average kelp cover
    # per site per period - 1 measurement - use sample year and month as period
   group_by(Study, Site, Sample_Year, Sample_Month) |>
    summarize(
      Latitude = mean_na(Latitude),
      Longitude = mean_na(Longitude),
      sample_unit_size_sq_m = mean_na(Sample_Unit_Size_sq_m),
      sasDate = mean_na(sasDate),
      max_Depth = max(max_Depth),
      mean_Depth = mean_na(mean_Depth),
      min_Depth = min(min_Depth),
      biomass_kg_wet_per_sq_m = mean_na(Biomass_kg_wet_per_sq_m),
      stipe_density_num_per_sq_m = mean_na(Stipe_Density_num_per_sq_m),
      individual_density_num_per_sq_m = mean_na(Individual_Density_num_per_sq_m),
      percent_cover = mean_na(Percent_Cover),
      # not sure why this was added - remove
    # biomass_KG_wet_per_sq_m_STD = sd(Biomass_kg_wet_per_sq_m),
    # Stipe_Density_num_per_sq_m_STD = sd(Stipe_Density_num_per_sq_m),
    # individual_per_sq_m_STD = sd(Individual_Density_num_per_sq_m),
    # percent_cover_STD = sd(Percent_Cover),
      .groups = "drop"
    ) |>
    rename_with(tolower) |>
    group_by(study) |>
    mutate(  has_bm = sum(!is.na(biomass_kg_wet_per_sq_m)),
             has_sd = sum(!is.na(stipe_density_num_per_sq_m)),
             has_id = sum(!is.na(individual_density_num_per_sq_m)),
             has_pc = sum(!is.na(percent_cover)),
             use_which = which.max(c(has_bm[1], has_sd[1], has_id[1], has_pc[1]))
          ) |>
    ungroup() |>
    
    
    # add focal kelp with preference order biomass, stipes, individuals, percent
    mutate(
      focalKelp = case_when(
        use_which==1  ~ biomass_kg_wet_per_sq_m,
        use_which==2  ~ stipe_density_num_per_sq_m,
        use_which==3  ~ individual_density_num_per_sq_m,
        use_which==4  ~ percent_cover,
        TRUE ~ NA_real_
      ),
      # focalKelp_STD = case_when(
      #   use_which==1  ~ biomass_kg_wet_per_sq_m_std,
      #   use_which==2  ~ stipe_density_num_per_sq_m_std,
      #   use_which==3  ~ individual_density_num_per_sq_m_std,
      #   use_which==4  ~ percent_cover_std,
      #   TRUE ~ NA_real_
      # ),
      focalUnit = case_when(
        use_which==1  ~ "biomass",
        use_which==2  ~ "stipe_density",
        use_which==3  ~ "individual_density",
        use_which==4  ~ "percent_cover",
        TRUE ~ NA_character_
      )
      )
  
  ret

}


# get the top 3 values, average them, and standardize by them
standardize_by_max <- function(x, samp=3) {
  m <- mean(sort(x, decreasing=T, na.last=T)[1:samp])
  x/m
}

ln_stdMax <- function(x, samp=3) {
  m <- mean(sort(x, decreasing=T, na.last=T)[1:samp])
  log(x/m+1/m)
}
