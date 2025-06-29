library(readr)
library(dplyr)

nwa_dat <- read_csv("data/nwa_with_env.csv") |>
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
         p_cover = p_cover>1, 1, p_cover)
