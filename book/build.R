# Builds the guide into book/_book. Run from the book folder:
#
#   Rscript build.R
#
# Requirements: R with bookdown, pkgload and ggplot2; Python 3.9 or later with
# NumPy (interpreter in QNCA_PYTHON, default "python", arguments after "|");
# Julia 1.6 or later (command in QNCA_JULIA, default "julia"); an internet
# connection for the datasets.
bookdown::render_book("index.Rmd", quiet = TRUE)

# The code links of qnca functions point to rdrr.io, which lists packages on
# CRAN only, so the function names stay plain text until qnca is on CRAN.
for (f in list.files("_book", pattern = "\\.html$", full.names = TRUE)) {
  h <- readLines(f, encoding = "UTF-8", warn = FALSE)
  h <- gsub('<a href="https://rdrr\\.io/pkg/qnca/man/[^"]*">([^<]*)</a>', "\\1", h)
  writeLines(h, f, useBytes = TRUE)
}
