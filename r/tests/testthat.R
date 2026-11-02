library(testthat)
library(qnca)

args <- commandArgs(FALSE)
script_arg <- args[grepl("^--file=", args)]
if (length(script_arg) > 0L) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1L])))
} else {
  script_dir <- getwd()
}

test_dir(file.path(script_dir, "testthat"))
