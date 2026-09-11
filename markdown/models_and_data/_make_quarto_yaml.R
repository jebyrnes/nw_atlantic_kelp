library(here)
library(glue)

chapter_dir <- here("markdown", "models_and_data", "chapters")

# discover dataset chapters automatically, sorted for stable ordering
dynamic_chapters <- list.files(chapter_dir, pattern = "\\.qmd$") |>
  sort() |>
  file.path("chapters", ... = _)  # -> "chapters/site-a.qmd", etc.

# static chapters, in the order you want them to appear
static_front <- c("index.qmd")#, "methods.qmd")
#static_back  <- c("discussion.qmd", "references.qmd")

all_chapters <- c(static_front, dynamic_chapters)#, static_back)

template <- readLines(here("markdown", "models_and_data", "_quarto-template.yml"))

chapter_yaml <- c("  chapters:", glue("    - {all_chapters}"))

# splice the chapters block onto the end of the book: section
final_yml <- c(template, chapter_yaml)

writeLines(final_yml, here("markdown", "models_and_data", "_quarto.yml"))