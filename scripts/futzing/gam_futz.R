library(mgcv)
library(modelbased)
library(patchwork)

# https://stats.stackexchange.com/questions/637423/conceptual-interpretation-of-bs-fs-and-by-term-in-gam/637540#637540
# https://stats.stackexchange.com/questions/657140/smooth-factor-interactions-in-generalized-additive-models/657164#657164

p_by <- gam(bill_len ~ 
              species +
              s(flipper_len, k = 5) +
              s(flipper_len, by = species, k = 5, m = 1),
            data = penguins,
            method = "REML") 

p_re <- gam(bill_len ~ 
           s(flipper_len, k = 5) +
           s(species, bs = "re")+ 
         s(species, flipper_len, bs = "re", k = 5), 
         data = penguins,
         method = "REML") 


p_fs <- gam(bill_len ~ 
              s(flipper_len, k = 5) +
              s(species, bs = "re") + 
              s(species, flipper_len, m = 1, bs = "fs", k = 5), 
            data = penguins,
            method = "REML")

p_sz <- gam(bill_len ~ 
              s(flipper_len, k = 5) +
              s(flipper_len, species, bs = "sz", k = 5, id = 1), 
            data = penguins,
            method = "REML")

g_by <- estimate_relation(p_by, by = c("flipper_len", "species"), length = 200) |> 
  plot()+ labs(title = "BY")


g_re <- estimate_relation(p_re, by = c("flipper_len", "species"), length = 200) |> 
  plot()+ labs(title = "RE")


g_fs <- estimate_relation(p_fs, by = c("flipper_len", "species"), length = 200) |> 
  plot()+ labs(title = "FS")


g_sz <- estimate_relation(p_sz, by = c("flipper_len", "species"), length = 200) |> 
  plot()+ labs(title = "SZ")


 g_re + g_fs + 
   g_by + g_sz + 
  plot_layout(guides = "collect",
              axis_titles = "collect",
              ncol = 2) &
  theme(legend.position='bottom')
 
 )
###
library(mgcv)
library(gratia)
library(ggplot2)

p_by <- gam(bill_len ~ 
              species +
              s(flipper_len, k = 5) +
              s(flipper_len, by = species, k = 5, m = 1),
            data = penguins) 

modelbased::estimate_relation(p_by,
                              by = c("flipper_len", "species"),
                              length = 400) |>
  plot(show_data = TRUE)

change <- response_derivatives(
  p_by,
  1focal = "flipper_len"
)

change_by_species <- response_derivatives(
  p_by,
  focal = "flipper_len",
  terms = c("flipper_len", "species")
)

modelbased::estimate_slopes(p_by, 
                            trend = "flipper_len",
                            by = c("flipper_len", "species"),
                            length = 400) |> plot()

head(change_by_species)


###

insight::get_datagrid(nwa_dat |> 
                        select(year_c, eco_collapsed),
                      by = c("year_c", "eco_collapsed"),
                      length = 1e4,
                      preserve_range = TRUE)
nd <- nwa_dat |>
  group_by(eco_collapsed) |>
  reframe(year_c = 
              seq(from = min(year_c), 
                  to = max(year_c), 
                  length.out = 1e4),
          study = study[1],
          trajectory = trajectory[1])
