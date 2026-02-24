library(ggplot2)
theme_set(theme_light(base_size = 12))

# for better coefficient understanding
coef_to_prec<-scales::new_transform(
  "coef_to_perc",
  transform = \(.x) exp(.x)-1,
  inverse = \(.y) log(.y+1) 
)

###
# we often need to filter predictions to years
# we actually have in each ecoregion
###
filter_time_eco_collapsed <- function(dat, reference = nwa_dat){
  reference_timing <- reference |>
    group_by(eco_collapsed) |>
    summarize(min_year_reference = min(year),
              max_year_reference = max(year))
  
  left_join(dat, reference_timing) |>
    filter(year > min_year_reference & year < max_year_reference)
  
}

##
# some ggplot2 scales
##
scale_color_slope_c <- function(...){
  scale_color_distiller(palette = "RdYlBu", direction =2, ...) 
}

scale_color_slope_b <- function(n.breaks = 5, ...){
  scale_color_fermenter(palette = "RdYlBu", direction =2, 
                        n.breaks=n.breaks,
                        ...) 
}

scale_fill_slope_c <- function(...){
  scale_fill_distiller(palette = "RdYlBu", direction =2, ...) 
}

scale_fill_slope_b <- function(n.breaks = 5, ...){
  scale_fill_fermenter(palette = "RdYlBu", direction =2, 
                        n.breaks=n.breaks,
                        ...) 
}
