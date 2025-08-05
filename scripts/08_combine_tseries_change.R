#' -----------------------------------------
#' Combine Trends and drivers
#' -----------------------------------------


library(dplyr)
library(readr)
library(ggplot2)
library(performance)
library(parameters)
library(modelbased)
library(glmmTMB)
library(broom.mixed)

had_wave_mod <- readRDS("models/hadsst_wave_mod.rds")
had_wave_coefs <- had_wave_mod |> tidy()
all_envt_mod <- readRDS("models/all_envt_mod.rds")
all_envt_coefs <- tidy(all_envt_mod)
envt_coefs <- read_csv("tables/envt_coefs.csv") |>
  rename(term = `...1`)

se_prod <- function(mx, my, sx, sy){
  sqrt((sx^2 + mx^2)*(sy^2+my^2) - mx^2*my^2)
}

###########
## Make a table of change based on had_wave_mode
###########

had_wave_coefs <- had_wave_coefs |>
  select(term, estimate, std.error) |>
  filter(term %in% 
           c("lag_hadsst_summer",
             "hadsst_spring",
             "swh_winter",
             "lag_swh_fall"))
  
had_change <- envt_coefs |>
  filter(term %in% 
           c("HADSST summer",
             "HADSST spring",
             "SWH winter",
             "SWH fall")) 

had_change <- had_change[c(1,2,4,3),]

had_change$`% change per year` <- 
  (had_wave_coefs[,2] * had_change[,2])[,1]

had_change$`se % change per year` <- se_prod(
  had_wave_coefs[,2],
  had_change[,2],
  had_wave_coefs[,3],
  had_change[,3])[,1]

had_est_change <- had_change |>
  select(-c(coefficient, se, p)) |>
  rename(`driver` = term) |>
  mutate(`% change per year` =
           (exp(`% change per year`)-1)*100,
         `se % change per year` =
           (exp(`se % change per year`)-1)*100,
         driver = 
           factor(driver,
                  levels = rev(c(
                    "HADSST summer",
                    "HADSST spring",
                    "SWH fall",
                    "SWH winter"
                  )))
           
         )
   
  
write_csv(had_est_change, "tables/had_change_est.csv")


###########
## Make a table of change based on all_envt_mod
###########

all_envt <- all_envt_coefs|>
  select(term, estimate, std.error) |>
  filter(term %in% 
           c("lag_oisst_summer",
             "oisst_spring",
             "swh_winter",
             "lag_swh_fall",
             "no3_spring",
             "no3_winter",
             "kd490_spring"))

all_envt_change <- envt_coefs |>
  filter(!(grepl("HADSST", term)))

all_envt_change <- all_envt_change[c(1,2,4,3,6,5,7),]

all_est_change <- all_envt_change |>
  select(-c(coefficient, se, p)) |>
  rename(`driver` = term) |>
  mutate(`% change per year` = 
           (all_envt[,2] * all_envt_change[,2])[,1],
         `se % change per year` =se_prod(
           all_envt[,2],
           all_envt_change[,2],
           all_envt[,3],
           all_envt_change[,3])[,1])|>
  mutate(`% change per year` =
           (exp(`% change per year`)-1)*100,
         `se % change per year` =
           (exp(`se % change per year`)-1)*100,
         driver = 
           factor(driver,
                  levels = rev(c(c("OISST summer", "OISST spring", "SWH winter", "SWH fall", "NO3 spring", 
                                   "NO3 winter", "Turbidity spring"))))
  ) 

write_csv(all_est_change, "tables/all_change_est.csv")

##
# Plots
##

ggplot(data= had_est_change,
       aes(y = driver,
           x = `% change per year`,
           xmin = `% change per year` - 1.96*`se % change per year`,
           xmax = `% change per year` + 1.96*`se % change per year`
       )) +
  geom_linerange() +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, 
             color = "red",
             lty = 2) +
  labs(y = "",
       x = "% change in kelp abundance per year due to driver")


ggsave("figures/had_wave_change.jpg", width = 5.4, height = 3)


# all

ggplot(data= all_est_change,
       aes(y = driver,
           x = `% change per year`,
           xmin = `% change per year` - 1.96*`se % change per year`,
           xmax = `% change per year` + 1.96*`se % change per year`
       )) +
  geom_point(size = 2) +
  geom_linerange() +
  geom_vline(xintercept = 0, 
             color = "red",
             lty = 2) +
  labs(y = "",
       x = "% change in kelp abundance per year due to driver")

ggsave("figures/all_change_est.jpg", width = 5.2, height = 3)


100*exp(70*log(1+.01*sum(all_est_change$`% change per year`)))
