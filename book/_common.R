# Shared setup of every chapter: chunk options, the qnca R package of the
# repository, and the Python and Julia engines that run their chunks with the
# packages of the repository.

knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  fig.align = "center",
  out.width = "85%",
  fig.width = 6.5,
  fig.asp = 0.62,
  dpi = 110
)

repo_root <- normalizePath("..", winslash = "/")
suppressPackageStartupMessages(
  pkgload::load_all(file.path(repo_root, "r"), export_all = FALSE, quiet = TRUE)
)
library(ggplot2)
theme_set(theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey20")))

run_external <- function(cmd, args, code, ext) {
  file <- tempfile(fileext = ext)
  writeLines(code, file, useBytes = TRUE)
  if (.Platform$OS.type == "windows" && file.exists(cmd)) cmd <- utils::shortPathName(cmd)
  out <- suppressWarnings(system2(cmd, c(args, shQuote(file)), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop("The chunk failed:\n", paste(out, collapse = "\n"), call. = FALSE)
  }
  out
}

python_cmd <- strsplit(Sys.getenv("QNCA_PYTHON", "python"), "|", fixed = TRUE)[[1]]
Sys.setenv(PYTHONPATH = file.path(repo_root, "python", "src"))
julia_cmd <- Sys.getenv("QNCA_JULIA", "julia")

knitr::knit_engines$set(
  python = function(options) {
    out <- if (isTRUE(options$eval)) {
      run_external(python_cmd[1], python_cmd[-1], options$code, ".py")
    } else character()
    knitr::engine_output(options, options$code, out)
  },
  julia = function(options) {
    out <- if (isTRUE(options$eval)) {
      run_external(julia_cmd, c("--startup-file=no",
                                paste0("--project=", file.path(repo_root, "julia", "QuantileNecessity.jl"))),
                   options$code, ".jl")
    } else character()
    knitr::engine_output(options, options$code, out)
  }
)
