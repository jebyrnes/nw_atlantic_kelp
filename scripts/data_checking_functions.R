#####################################################################################
#' Functions to check processed datasets
#'
#' Author: Jarrett Byrnes jarrett.byrnes@umb.edu fixed up from Chat-GPT 3.5 translation of old SAS code
#' Dataset Last Updated:  
#' Script Last Updated: June 27, 2023
#'
#' Changelog
#'  20230727 - added in reading data frames and checking for duplicates and missing info
#####################################################################################

pacman::p_load(dplyr, tidyr, readr, glue, janitor, sqldf)

# debug/testing
#filename <- "temporal_data_hubline.csv"

# read in kelp data with correct column types
read_kelp_data <- function(filename, path = "data/clean_data/timeseries/"){
  read_csv(glue("{path}{filename}"), 
                        col_types = cols(
                          Entry = col_double(),
                          Study = col_character(),
                          Region = col_character(),
                          Site = col_character(),
                          Depth.m = col_double(),
                          Sample.ID = col_character(),
                          Sample.Year = col_double(),
                          Sample.Month = col_double(),
                          Sample.Day = col_double(),
                          Latitude = col_double(),
                          Longitude = col_double(),
                          Taxon = col_character(),
                          Fixed.or.Random = col_character(),
                          Sample.Unit.Size.sq.m = col_double(),
                          Percent.Cover = col_double(),
                          Stipe.Density.num.per.sq.m = col_double(),
                          Individual.Density.num.per.sq.m = col_double(),
                          Biomass.kg.wet.per.sq.m = col_double()
                        )
  ) |> clean_names(case = "none") |>
    #cleaning
    mutate(Sample_ID = ifelse(is.na(Sample_ID), Site, Sample_ID),
           Sample_Month = ifelse(is.na(Sample_Month), 6, Sample_Month),
           Sample_Day = ifelse(is.na(Sample_Day), 15, Sample_Day),
           Depth_m = ifelse(is.na(Depth_m), 9999, Depth_m),
           Sample_Unit_Size_sq_m = ifelse(is.na(Sample_Unit_Size_sq_m), 9999, Sample_Unit_Size_sq_m)
    )
}




# function to check for problems with dups, etc
check_problems <- function(kelp_data) {
  outName <- kelp_data$Study[1]
  
  # is there any missing info about the sample?
  missing_info <- kelp_data %>%
    filter(
      is.na(Site) |
        is.na(Depth_m) |
        is.na(Sample_Year) |
        is.na(Sample_Month) |
        is.na(Latitude) |
        is.na(Longitude) |
        is.na(Taxon) |
        is.na(Sample_ID) |
        is.na(Fixed_or_Random) |
        is.na(Study)
    )
  
  # are any of the data rows NA?
  missing_data <- kelp_data %>%
    filter(
      is.na(Percent_Cover) &
        is.na(Biomass_kg_wet_per_sq_m) &
        is.na(Individual_Density_num_per_sq_m) &
        is.na(Stipe_Density_num_per_sq_m)
    )
  
  
  
  # re-arrange to something standard
  kelp_data <- kelp_data %>%
    arrange(
      Site,
      Depth_m,
      Sample_Year,
      Sample_Month,
      Sample_Day,
      Latitude,
      Longitude,
      Taxon,
      Sample_ID
    )
  
  # check to see if there are duplicate entries
  pre_duplicates <- kelp_data %>%
    group_by(
      Site,
      Depth_m,
      Sample_Year,
      Sample_Month,
      Sample_Day,
      Latitude,
      Longitude,
      Taxon,
      Sample_ID
    ) %>%
    summarise(smallestEntry = min(Entry),
              largestEntry = max(Entry),
              .groups = "keep")
  
  duplicates <- pre_duplicates %>%
    filter(n() > 1)
  
  checkDouble <- kelp_data %>%
    group_by(
      Site,
      Depth_m,
      Sample_Year,
      Sample_Month,
      Sample_Day,
      Latitude,
      Longitude,
      Taxon,
      Sample_ID
    ) %>%
    filter(n() > 1)
  
  # function for writing out errors
  
  
  # function to write out problems
  exportIfRows <- function(dataName, out, errorName, isFatal) {
    query <- paste0("SELECT COUNT(*) AS nobs FROM ", dataName)
    mvar <- sqldf(query)
    if (mvar$nobs > 0) {
      cat(paste(out, errorName, "\n", sep = " "))
      
      NewDirectory <- paste("clean_data/dataQA/", out, sep = "")
      dir.create(NewDirectory, showWarnings = FALSE, recursive = TRUE)
      
      outfile <- paste("clean_data/dataQA/", out, "/", dataName, ".csv", sep = "")
      df <- eval(parse(text = dataName))
      write.csv(df, file = outfile, row.names = FALSE)
    }
    
    query2 <- paste0("SELECT COUNT(*) AS nobs FROM ", dataName)
    nobs <- sqldf(query2)
    
    if (nobs$nobs > 0) {
      Dataset_ErrorsT <- data.frame(
        study = out
      ) %>%
        mutate(
          #    study_ID = setStudyID,
          Error_type = errorName,
          Fatal = isFatal
        )
      
      Dataset_errors <- Dataset_ErrorsT
    }
  }
  
  
  # write out any errors
  exportIfRows("missing_info", outName, "missing_info", 1)
  exportIfRows("missing_data", outName, "No_kelp_data", 1)
  exportIfRows("checkDouble", outName, "Duplicate rows", 1)
  
  
}