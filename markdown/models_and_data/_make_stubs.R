library(here)
library(glue)
library(purrr)
library(janitor)
setwd(here::here())

## load in kelptime
source("scripts/load_nwa_data.R")
studies <- unique(nwa_dat$study) |> sort()

## what are the studies in clean_data?
clean_dat <- list.files("data/clean_data/timeseries",
                        full.names = TRUE)

analyzed_studies <- read_csv(clean_dat, id = "path") |>
  rename_all(tolower) |>
  select(path, study) |>
  group_by(path, study) |>
  slice(1L) |>
  ungroup() |>
  filter(study %in% studies) |>
  mutate(file = make_clean_names(study))

# 
# sites <- tibble::tribble(
#   ~file,     ~name,     ~dataset,
#   "site-a",  "Site A",  "site-a.csv",
#   "site-b",  "Site B",  "site-b.csv",
#   "site-c",  "Site C",  "site-c.csv"
#   # site-d.csv removed -> its stub will now get cleaned up automatically
# )

chapter_dir <- here("markdown", "models_and_data", "chapters")

# --- 1. write/update a stub for every current dataset ---
pwalk(analyzed_studies[1:3], function(path, study, file) {
  glue(
    "---\n",
    "title: \"{study}\"\n",
    "params:\n",
    "  path: \"{path}\"\n",
    "  dataset_name: \"{study}\"\n",
    "---\n\n",
    "{{{{< include ../_analysis-template.qmd >}}}}\n"
  ) |> writeLines(file.path(chapter_dir, glue("{file}.qmd")))
})

# --- 2. remove any stub that no longer has a matching dataset entry ---
expected_files <- glue("{analyzed_studies$file}.qmd")
existing_files <- list.files(chapter_dir, pattern = "\\.qmd$")

stale_files <- setdiff(existing_files, expected_files)

if (length(stale_files) > 0) {
  message("Removing stale chapter stubs: ", paste(stale_files, collapse = ", "))
  file.remove(file.path(chapter_dir, stale_files))
}