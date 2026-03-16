library(quarto)
countries <- c("Afghanistan", "Belgium", "India", "United Kingdom")
setwd(paste0(here::here(), "/markdown/test/"))

for (country in countries) {
  quarto_render(
    input = "life_expectancy_report.qmd",
    output_file = paste0("life_expectancy_", country, ".html"),
    execute_params = list(country = country),
    execute_dir = here::here()
  )
}
