source(here::here("markdown", "models_and_data", "_make_stubs.R"))
source(here::here("markdown", "models_and_data", "_make_quarto_yaml.R"))
quarto::quarto_render(here::here("markdown", "models_and_data"), as_job = FALSE)
