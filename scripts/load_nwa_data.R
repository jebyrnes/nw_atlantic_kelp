library(readr)
library(dplyr)
library(forcats)

nwa_dat <- read_csv("data/kelptime_nwa_data.csv", show_col_types = FALSE) |>
  # Rescaling
  mutate(rescaled_std_by_ecoregion = scales::rescale(focal_std_by_ecoregion, c(0,1)),
         year_z = (year - mean(year))/sd(year)) |>
  # for plotting, get N-S
  mutate(ecoregion = factor(ecoregion, levels = 
                              c("Virginian",
                                "Gulf of Maine/Bay of Fundy",
                                "Scotian Shelf",
                                "Gulf of St. Lawrence - Eastern Scotian Shelf",
                                "Southern Grand Banks - South Newfoundland",  
                                "Northern Grand Banks - Southern Labrador"
                              ))) |>
  # Some nicer names
  mutate(focalUnit = case_when(
    focalUnit == "biomass" ~ "Biomass",
    focalUnit == "percent_cover" ~ "Percent cover",
    focalUnit == "stipe_density" ~ "Stipe density"
  )) |>
  # for ordbeta of only percent cover
  mutate(p_cover = percent_cover/100,
         p_cover = p_cover>1, 1, p_cover)|>
  # to combine regions without much data
  mutate(eco_collapsed = 
           fct_collapse(ecoregion,
                        "Gulf of St. Lawrence - Newfoundland" = 
                          c("Gulf of St. Lawrence - Eastern Scotian Shelf",
                            "Southern Grand Banks - South Newfoundland",
                            "Northern Grand Banks - Southern Labrador"))) |>
  # french characters giving trouble
  mutate(site = ifelse(site == "I_le aux\nGoe\x89lands", "I_le aux Goelands", site),
         trajectory = ifelse(trajectory == "Gagnon_et_al._2005 I_le aux\nGoe\x89lands", "Gagnon_et_al._2005 I_le aux Goelands", trajectory)) |>
  # for decadal analysis
  mutate(decade = (floor(year/10)*10) |> 
           paste0("s") |>
           as.factor(),
         decade = fct_collapse(decade,
                               "pre-1980s" = c(
                                 "1970s",
                                 "1960s",
                                 "1940s"
                               )),
         decade = fct_collapse(decade,
                               "post-2010" = c(
                                 "2010s",
                                 "2020s"
                               )),
         # for GAMs
         eco_collapsed = as.factor(eco_collapsed),
         trajectory = as.factor(trajectory),
         study = as.factor(study))

maine_scotia <- nwa_dat |>
  filter(eco_collapsed %in% c("Gulf of Maine/Bay of Fundy",
                              "Scotian Shelf"))
