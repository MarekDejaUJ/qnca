args <- commandArgs(FALSE)
script_arg <- args[grepl("^--file=", args)]
if (length(script_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1L]))
  root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
} else {
  root <- normalizePath(getwd())
}
setwd(root)

scripts <- c(
  "frontier_geometry.R",
  "pi_reduction.R",
  "outlier_fragility.R",
  "power_curve.R"
)

for (script in scripts) {
  source(file.path(root, "replication", "R", "figures", script), chdir = FALSE)
}

message("Figures written to ", file.path(root, "replication", "figures"))
