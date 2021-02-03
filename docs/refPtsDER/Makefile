export PATH := /Library/TeX/texbin:$(PATH)    # add LaTeX path
export PATH := /usr/local/bin:~/.local/bin:$(PATH) # add pandoc-citeproc-preamble

# Cluster targets
pdf: index.pdf tandf.pdf appendices.pdf supplemental.pdf
docx: index.docx tandf.pdf appendices.pdf

# coauthor renders
coauth: intro.pdf intro.docx tandf.pdf

all: pdf docx

## pdf
index.pdf: index.Rmd Makefile
	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, output_format = "bookdown::pdf_document2", config_file = "_bookdown.yml")'

# docx
index.docx: index.Rmd Makefile
	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, output_format = "bookdown::word_document2", config_file = "_bookdownDocx.yml")'


index.html: index.Rmd Makefile
	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, output_format = "bookdown::html_document2", config_file = "_bookdown.yml")'	

intro.pdf: index.Rmd Makefile
	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, output_format = "bookdown::pdf_document2", config_file = "_intro.yml")'

intro.docx: index.Rmd Makefile
	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, output_format = "bookdown::word_document2", config_file = "_introDocx.yml")'

tandf.pdf: index.Rmd Makefile
	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, output_format = "bookdown::pdf_document2", config_file = "_tandf.yml")'


appendices.pdf: index.Rmd Makefile
	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, output_format = "bookdown::pdf_document2", config_file = "_appendices.yml")'


# docx
appendices.docx: index.Rmd Makefile
	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, output_format = "bookdown::word_document2", config_file = "_appendices.yml")'


supplemental.pdf: index.Rmd Makefile
	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, output_format = "bookdown::pdf_document2", config_file = "_supplemental.yml")'