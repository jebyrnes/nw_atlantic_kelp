#' -----------------------------------------------
#' NWA Timeseries Analysis and Viz
#'
#' -----------------------------------------------

library(dplyr)
library(ggplot2)
library(performance)
library(modelbased)
library(glmmTMB)
library(broom.mixed)
library(car)

theme_set(theme_bw(base_size = 14))

# for better coefficient understanding
coef_to_prec<-scales::new_transform(
  "coef_to_perc",
  transform = \(.x) exp(.x)-1,
  inverse = \(.y) log(.y+1) 
)

param_plot_ln <- function(mod){
  parameters::model_parameters(mod,
                               effects = "fixed",
                               component = "conditional",
                               include_sigma = FALSE,
                               drop = "(cellmean|ecoregion)",
                               ci = 0.8) |> plot() +
    scale_color_manual(values = c("black", "black")) +
    scale_x_continuous(transform = coef_to_prec) +
    labs(x = "% change per 1 unit change in x")
}

#reads in nwa data
source("scripts/load_nwa_data.R")

oisst_only_mod <- glmmTMB(ln_focal_std_by_ecoregion ~ 
                     ecoregion + year + 
                     lag_oisst_summer + 
                     oisst_spring +
                     lag_oisst_summer_cellmean +
                     oisst_spring_cellmean +
                     (1  +year|trajectory) + (1|study), 
                   dispformula =~focalUnit,
                   data = nwa)

oisst_only_mod_ordbera <- glmmTMB(rescaled_std_by_ecoregion ~ 
                            ecoregion + year + 
                            lag_oisst_summer + 
                            oisst_spring +
                            lag_oisst_summer_cellmean +
                            oisst_spring_cellmean +
                            (1 + year|trajectory) + (1|study), 
                          dispformula =~focalUnit,
                          family = ordbeta,
                          data = nwa)

hadsst_only_mod <- glmmTMB(ln_focal_std_by_ecoregion ~ 
                            ecoregion + year + 
                            lag_hadsst_summer + 
                            hadsst_spring +
                            lag_hadsst_summer_cellmean +
                            hadsst_spring_cellmean +
                            (1  |trajectory) + (1|study), 
                          dispformula =~focalUnit,
                          data = nwa)
car::Anova(oisst_only_mod)
car::Anova(hadsst_only_mod)

##
# Era 5 and Hadsst
##

hadsst_wave_mod <- glmmTMB(ln_focal_std_by_ecoregion ~ 
                             ecoregion + year_c +
                             #temp
                             lag_hadsst_summer + 
                             lag_hadsst_summer_cellmean +
                             
                             hadsst_spring +
                             hadsst_spring_cellmean +
                             
                             # waves
                             swh_winter + 
                             swh_winter_cellmean +
                             
                             lag_swh_fall +
                             lag_swh_fall_cellmean +
                             
                             (1+year  |trajectory) + (1|study), 
                           dispformula =~focalUnit,
                           data = nwa)

param_plot_ln(hadsst_wave_mod)


##
# Era 5 and Hadsst
##

oisst_wave_mod <- glmmTMB(ln_focal_std_by_ecoregion ~ 
                             0 + ecoregion + 
                             #temp
                             lag_oisst_summer + 
                             lag_oisst_summer_cellmean +
                             
                             oisst_spring +
                             oisst_spring_cellmean +
                             
                             # waves
                             swh_winter + 
                             swh_winter_cellmean +
                             
                             lag_swh_fall +
                             lag_swh_fall_cellmean +
                             
                             (1+year  |trajectory) + (1|study), 
                           dispformula =~0 + focalUnit,
                           data = nwa)



oisst_wave_mod_ordbeta <- glmmTMB(rescaled_std_by_ecoregion ~ 
                            0 + ecoregion + 
                            #temp
                            lag_oisst_summer + 
                            lag_oisst_summer_cellmean +
                            
                            oisst_spring +
                            oisst_spring_cellmean +
                            
                            # waves
                            swh_winter + 
                            swh_winter_cellmean +
                            
                            lag_swh_fall +
                            lag_swh_fall_cellmean +
                            
                            (1+year  |trajectory) + (1|study), 
                          dispformula =~0 + focalUnit,
                          family = ordbeta,
                          data = nwa)

parameters::model_parameters(oisst_wave_mod,
                             effects = "fixed",
                             component = "conditional",
                             include_sigma = FALSE,
                             drop = c("cellmean", "ecoregion"),
                             ci = 0.8) |> plot()

##
# Era 5,  OISST, CMEM for no3
##

temp_wave_nut_mod <- glmmTMB(
  rescaled_std_by_ecoregion ~
    0 + ecoregion + year_c +
    #temp
    lag_oisst_summer +
    lag_oisst_summer_cellmean +
    
    oisst_spring +
    oisst_spring_cellmean +
    
    # waves
    swh_winter +
    swh_winter_cellmean +
    
    lag_swh_fall +
    lag_swh_fall_cellmean +
    
    # nutrients
    no3_spring +
    no3_spring_cellmean +
    
    no3_winter +
    no3_winter_cellmean + 
    
    (1 + year  | trajectory) + (1 | study),
  dispformula =  ~ 0 + focalUnit,
  data = nwa,
  REML = FALSE
)



parameters::model_parameters(temp_wave_nut_mod,
                             effects = "fixed",
                             component = "conditional",
                             include_sigma = FALSE,
                             drop = c("cellmean", "ecoregion"),
                             ci = 0.9) |> plot() +
  scale_color_manual(values = c("black", "black"))

##
# Era 5,  OISST, CMEM, kd490
##
##



all_mod <- glmmTMB(
  rescaled_std_by_ecoregion ~
    0 + ecoregion + year_c +
    #temp
    lag_oisst_summer +
    lag_oisst_summer_cellmean +
    
    oisst_spring +
    oisst_spring_cellmean +
    
    # waves
    swh_winter +
    swh_winter_cellmean +
    
    lag_swh_fall +
    lag_swh_fall_cellmean +
    
    # nutrients
    no3_spring +
    no3_spring_cellmean +
    
    no3_winter +
    no3_winter_cellmean + 
    
    # turbidity
    kd490_spring +
    kd490_spring_cellmean +
    
    (1 + year  | trajectory) + (1 | study),
  dispformula =  ~ 0 + focalUnit,
  data = nwa,
  REML = FALSE
)

Anova(all_mod)

param_plot_ln(all_mod)

ggsave("figures/all_envt_params.jpg", width = 4.5, height = 4)
