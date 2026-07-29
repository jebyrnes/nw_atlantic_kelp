#' -----------------------------------------------
#' NWA Timeseries Modeling
#' looking at variation across spacetime 
#' with sdmTMB
#' -----------------------------------------------

library(sdmTMB)
library(sdmTMBextra)
library(sf)
library(ggplot2)
library(glue)
library(purrr)
source("scripts/99_sdm_helpers.R") #also loads data and constants/helpers
source("scripts/99_gam_sdmtmb_helpers.R")

###
## Make objects we will use
##
####
valid_areas <- make_valid_areas(buffer = 10)
prediction_points <- full_prediction_points |> 
  st_intersection(valid_areas)

cutoffs <- c(5, 10, 15, 20, 30, 50, 80, 100)
k_vals <- 3:12

outlist <- vector("list", 
                  length(cutoffs)*length(k_vals))
out_idx <- 0

## for each mesh cutoff - 
## minimum triangle edge length - distance between locations
for(a_cutoff in cutoffs){
  # make the mesh
  mesh_coastline <- make_kelp_mesh(cutoff = a_cutoff)
  
 ## for each number of splines
  for(a_k in k_vals){
    
    print(glue("k = {a_k}, cutoff = {a_cutoff}"))
    out_idx <- out_idx+1
    # fill in the outlist
    outlist[[out_idx]]$cutoff <- a_cutoff
    outlist[[out_idx]]$k <- a_k
    
    

    # make the basis set
    
    # the smoother
    smooth_formula <- focal_std_by_all ~ s(year,  bs = "tp", k = a_k)
    #smooth_formula <- focal_std_by_all ~ s(year,  bs = "gp", k = a_k, m=2) # evaluate if GP is a good idea
    
    nwa_basis <- make_basis(smooth_formula, nwa_dat)
    
    fe_formula <- as.formula(paste("focal_std_by_all ~ 1 + ", 
                                   paste(colnames(nwa_basis$v), collapse = " + "),
                                   "+ (1|study) + (1 | trajectory)"))
    
    svc_formula <- as.formula(paste("~ 0 +", 
                                    paste(colnames(nwa_basis$sx), collapse = " + "), "+",
                                    paste(colnames(nwa_basis$sz), collapse = " + ")))
    
    ## fit the model
    mod_spatiotemporal<- sdmTMB(fe_formula,  
                                dispformula = ~focalUnit,
                                family = tweedie(link = "log"),
                                spatial = "on",
                                spatial_varying = svc_formula, 
                                mesh = mesh_coastline,
                                data = nwa_basis$dat,
                                control = 
                                  sdmTMBcontrol(map = list(ln_tau_Z = factor(c(1L, rep(2L, ncol(nwa_basis$sz))))))
    )
    
    # get the sanity check
    outlist[[out_idx]]$sanity <- sanity(mod_spatiotemporal, silent = TRUE)
    
    # get the model coefs
    outlist[[out_idx]]$coefs <- tidy(mod_spatiotemporal,  "ran_pars")
    
    # get the AIC
    outlist[[out_idx]]$aic <- AIC(mod_spatiotemporal)
    
    # get the predicted values for curves with CIs for extrapolation
    # for prediction
    prediction_points_yr <-
      tidyr::expand_grid(prediction_points, year = seq(
        min(nwa_dat_with_basis$year),
        max(nwa_dat_with_basis$year),
        by = 0.1
      ))
    
    pred_basis <- make_basis(smooth_formula, prediction_points_yr, basis_prev = nwa_basis$sm)
    
    pred_data <- pred_basis$dat |>
      as_tibble() |>
      mutate(
        trajectory = nwa_dat_with_basis$trajectory[1],
        study = nwa_dat_with_basis$study[1],
        focalUnit = nwa_dat_with_basis$focalUnit[1]
      )
    
    predicted <- predict(mod_spatiotemporal,
                         type = "response",
                         newdata = pred_data)
    
    predicted <- predicted |>
      mutate(eco_collapsed = 
               forcats::fct_relevel(eco_collapsed,
                                    c("Virginian" ,
                                      "Gulf of Maine/Bay of Fundy" ,
                                      "Scotian Shelf" ,
                                      "Gulf of St. Lawrence - Newfoundland")),
             ord = as.numeric(eco_collapsed)*1e5 + Y,
             id = rank(ord)
      )
    
    outlist[[out_idx]]$dat <- 
      predicted |> select(year, est, cell_id, Y, eco_collapsed) |>
      mutate(k = a_k, cutoff = a_cutoff)

    # output a hovmoller plot
    ggplot(predicted,
           aes(x = year, y = id, fill = est)) +
      geom_tile() +
      #scale_fill_viridis_c(option = "D") +
      colorspace::scale_fill_continuous_divergingx(palette = "BrBG", mid = 0.2, rev = FALSE) +
      facet_wrap(vars(forcats::fct_rev(eco_collapsed)), 
                 scale = "free_y", ncol = 1)+
      guides(y = "none") +
      labs(y = "", x = "", fill = "Standardized\nKelp Abundance") +
      theme_classic() 
    
    ggsave(glue("figures/robustness_check/spatiotemporal_hovmoller_cutoff={a_cutoff}_k={a_k}.jpg"), 
           width = 8, height = 5)
    
    # output a curve plot
    ggplot(predicted|> 
             filter_time_eco_collapsed()  ,
           aes(x = year, y = est, group = cell_id, color = Y)) +
      geom_line(alpha = 0.1)  +
      scale_color_viridis_c() +
      labs(y = "Standardized Kelp Abundance", x = "",
           color = "Northing") +
      theme_light(base_size = 18)
    
    ggsave(glue("figures/robustness_check/spatiotemporal_curves_cutoff={a_cutoff}_k={a_k}.jpg"), 
                      width = 8, height = 5)
  }
}

saveRDS(outlist, "figures/robustness_check/outlist_robust.rds")

# compare sanity
sanity_check <- map(outlist,
    ~tibble(k = .x$k, cutoff = .x$cutoff, 
            sigmas_ok = .x$sanity$sigmas_ok,
            all_ok = .x$sanity$all_ok)
    ) |> list_rbind()

readr::write_csv(sanity_check, "figures/robustness_check/sanity.csv")

# compare AICs
aictab <- map(outlist,
                    ~tibble(k = .x$k, cutoff = .x$cutoff, 
                            aic = .x$aic)
) |> list_rbind() |>
  mutate(
    delta_aic = aic - min(aic, na.rm = TRUE)
  ) |>
  arrange(aic)

readr::write_csv(aictab, "figures/robustness_check/aictab.csv")


# make unified data
outdat <- map(outlist,
              ~.x$dat) |> list_rbind()


## For each mesh cutoff, make a panel plot by number of knots showing
## curves
for(a_cutoff in cutoffs){
  outdat |>
    filter(cutoff == a_cutoff) |>
    filter_time_eco_collapsed() |>
    ggplot(aes(x = year, y = est, group = cell_id, color = Y)) +
    geom_line(alpha = 0.1)  +
    scale_color_viridis_c() +
    labs(y = "Standardized Kelp Abundance", x = "",
         color = "Northing") +
    theme_light(base_size = 18) +
    facet_wrap(~k, scale = "free_y")
  
  ggsave(glue("figures/robustness_check/effect_of_k_cutoff={a_cutoff}.jpg"),
         width = 12,
         height = 8)
}
  

## For each number of knots, make a panel plot by mesh size showing
## curves
for(a_k in k_vals){
  outdat |>
    filter(k == a_k) |>
    filter_time_eco_collapsed() |>
    ggplot(aes(x = year, y = est, group = cell_id, color = Y)) +
    geom_line(alpha = 0.1)  +
    scale_color_viridis_c() +
    labs(y = "Standardized Kelp Abundance", x = "",
         color = "Northing") +
    theme_light(base_size = 18) +
    facet_wrap(~cutoff, scale = "free_y")
  
  ggsave(glue("figures/robustness_check/effect_of_cutoff_k={a_k}.jpg"),
         width = 12,
         height = 8)
}
