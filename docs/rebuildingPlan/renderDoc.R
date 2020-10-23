# Render Rebuilding Plan document for HG Herring
# See more info on csasdown at:
# https://github.com/pbs-assess/csasdown

library(bookdown)

# setwd(here::here('docs/rebuildingPlan'))

# csasdown::draft("sr") # creates Science Response Template in working dir
bookdown::render_book(input = "index.Rmd", clean = TRUE, output_format = "csasdown::sr_pdf", config_file = "_bookdown.yml")
bookdown::render_book(input = "index.Rmd", clean = TRUE, output_format = "csasdown::sr_word", config_file = "_bookdown_docx.yml")

