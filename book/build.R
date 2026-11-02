# Builds the guide into book/_book. Run from the book folder:
#
#   Rscript build.R
#
# Requirements: R with bookdown, pkgload and ggplot2; Python 3.9 or later with
# NumPy (interpreter in QNCA_PYTHON, default "python", arguments after "|");
# Julia 1.6 or later (command in QNCA_JULIA, default "julia"); an internet
# connection for the datasets.
bookdown::render_book("index.Rmd", quiet = TRUE)
