library(quarto)
library(readr)
library(dplyr)
library(purrr)
setwd(here::here())

## load in kelptime
source("scripts/load_nwa_data.R")
studies <- unique(nwa_dat$study) |> sort()

## what are the studies in clean_data?
clean_dat <- list.files("data/clean_data/timeseries",
                        full.names = TRUE)


file_to_study <- read_csv(clean_dat, id = "path") |>
  rename_all(tolower) |>
  select(path, study) |>
  group_by(path, study) |>
  slice(1L) |>
  ungroup() |>
  filter(study %in% studies)

#####
# make the quarto outputs
#####

# reset wd for report generation
setwd("markdown/dataset_reports/")

for (i in 1:nrow(file_to_study)) {
  
  this_study <- file_to_study$study[i]
  this_path <- file_to_study$path[i]
  
  print(this_study)
  
  
  quarto_render(
    input = "dataset_report_template.qmd",
    output_file = paste0(this_study, "_timeseries_report.html"),
    execute_params = list(study = this_study, path = this_path),
    execute_dir = here::here()
  )
  
}


