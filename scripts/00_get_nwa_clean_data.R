#' -----------------------------------------------
#' Move files from global analysis to WNA project
#' so that we don't have to drag and drop
#' -----------------------------------------------
pacman::p_load(dplyr, readr)

# get ye olde data
nwa_old <- read.csv("https://github.com/kelpecosystems/global_kelp_time_series/raw/refs/heads/master/05_HLM_analysis_code/formatted_data_3years.csv")

nwa_old |>
  filter(PROVINCE == ("Cold Temperate Northwest Atlantic")) |>
  filter(!is.na(ln_stdByECOREGION)) |>
  write_csv("data/nwa_krumhansl_2016.csv")

setwd(here::here())
# move files that work
f <- c("dmr_rasher_steneck.csv", "filbee_dexter_arctic_kelp.csv", 
       "Godbout_StLawrence.csv", "gulf_of_maine_ios_southwestern_gom.csv", 
       "hubline.csv", "keen_one.csv", "Mined_Data.csv", "nova_scotia_kelp.csv", 
       "Rhode_Island_Feehan_Grace_Narvaez.csv", "seabrook.csv",
       "Old_Mined_Data.csv")


file.copy(from = paste0("../clean_data/timeseries/", f),
          to = paste0("data/clean_data/timeseries/", f),
          overwrite = TRUE)

fold <- c("temporal_data_steneck_gom.csv", "temporal_data_witman_gom.csv")
  
file.copy(from = paste0("../old_clean_data/", fold),
          to = paste0("data/clean_data/timeseries/", fold),
          overwrite = TRUE)

# filter mined data and write
mined <- read_csv("data/clean_data/timeseries/Mined_Data.csv", na = c("", "NA")) |>
  bind_rows(read.csv("data/clean_data/timeseries/Old_Mined_Data.csv")) |>
  filter(Study %in% c("Attridge_et_al._2022", "Witman_et_al._2018",
                      "O'Brien_et_al._2018", "Feehan_et_al._2019", 
                      "Filbee-Dexteretal_et_al._2016", "Gagnon_et_al._2005",
                      "Gagnon_et_al_2004", "Steneck_et_al._2013" , "Steneck_1997",
                      "Witman_1987", "Siddon_&_Witman_2003", "Egan_&_Yarish__1990"))

write_csv(mined, "data/clean_data/timeseries/GOM_mined.csv")

# remove files no longer needed
file.remove("data/clean_data/timeseries/Old_Mined_Data.csv")
file.remove("data/clean_data/timeseries/Mined_Data.csv")
