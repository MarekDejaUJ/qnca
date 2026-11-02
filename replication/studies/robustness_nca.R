# Robustness evaluation, NCA ceilings and outlier screen (design in README.md).
#
# Rebuilds every sample of studies S1 and S2 from the base data written by
# robustness.jl and computes the NCA CE-FDH, CR-FDH and QR ceilings with
# NCA 5.0.2. On the first 20 replicates of S1 (n <= 300) it also records the
# first case listed by nca_outliers for comparison with the native screen.
#
#   Rscript replication/studies/robustness_nca.R [cores]

suppressPackageStartupMessages({
  library(NCA)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 1L) as.integer(args[[1L]]) else 8L
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
root <- normalizePath(file.path(script_dir, ".."))
out_dir <- file.path(root, "results")
data_dir <- file.path(out_dir, "data")
if (packageVersion("NCA") != "5.0.2") stop("NCA 5.0.2 is required")

SCOPE <- c(0, 100, 0, 100)
M_LEVELS <- c(1, 2, 3, 5)
SIGMAS <- c(0, 2, 5, 10)
CHECK_REPS <- 20L

line_bottleneck <- function(slope, intercept, y = 50) {
  if (!is.finite(slope) || !is.finite(intercept)) return(NA_real_)
  if (slope <= 0) return(if (intercept >= y) SCOPE[1] else NA_real_)
  min(SCOPE[2], max(SCOPE[1], (y - intercept) / slope))
}

nca_rows <- function(X, Y) {
  d <- data.frame(X = X, Y = Y)
  m <- nca_analysis(d, "X", "Y", ceilings = c("ce_fdh", "cr_fdh", "qr"), scope = SCOPE)
  p <- m$summaries$X$params
  eff <- function(ceil) as.numeric(p["Effect size", ceil])
  hit <- Y >= 50
  data.frame(
    method = c("NCA CE-FDH", "CR-FDH", "QR"),
    d = c(eff("ce_fdh"), eff("cr_fdh"), eff("qr")),
    b50 = c(if (any(hit)) min(X[hit]) else SCOPE[2],
            line_bottleneck(as.numeric(p["Slope", "cr_fdh"]), as.numeric(p["Intercept", "cr_fdh"])),
            line_bottleneck(as.numeric(p["Slope", "qr"]), as.numeric(p["Intercept", "qr"]))),
    stringsAsFactors = FALSE
  )
}

first_outlier <- function(X, Y) {
  o <- nca_outliers(data.frame(X = X, Y = Y), "X", "Y", ceiling = "ce_fdh", scope = SCOPE)
  if (is.null(o) || nrow(o) == 0L) return(c(NA_integer_, NA_real_))
  c(as.integer(o$outliers[1]), as.numeric(o$dif.abs[1]))
}

run_s1 <- function(gen, n) {
  clean <- read.csv(file.path(data_dir, sprintf("s1_clean_%s_n%d.csv", gen, n)))
  planted <- read.csv(file.path(data_dir, sprintf("s1_planted_%s_n%d.csv", gen, n)))
  reps <- sort(unique(clean$rep))
  res <- mclapply(reps, function(rep) {
    base <- clean[clean$rep == rep, ]
    rows <- list(); checks <- list()
    r <- nca_rows(base$X, base$Y)
    rows[[length(rows) + 1L]] <- cbind(band = "clean", level = 0, r)
    for (band in c("mid", "top")) {
      pts <- planted[planted$rep == rep & planted$band == band, ]
      pts <- pts[order(pts$j), ]
      for (m in M_LEVELS) {
        X <- c(base$X, pts$x[seq_len(m)]); Y <- c(base$Y, pts$y[seq_len(m)])
        r <- nca_rows(X, Y)
        rows[[length(rows) + 1L]] <- cbind(band = band, level = m, r)
        if (n <= 300 && rep <= CHECK_REPS) {
          fo <- first_outlier(X, Y)
          checks[[length(checks) + 1L]] <- data.frame(rep = rep, band = band, m = m,
                                                      nca_case = fo[1], nca_dif_abs = fo[2])
        }
      }
    }
    list(rows = cbind(rep = rep, do.call(rbind, rows)),
         checks = if (length(checks)) do.call(rbind, checks) else NULL)
  }, mc.cores = cores)
  errs <- vapply(res, inherits, logical(1), "try-error")
  if (any(errs)) stop(sprintf("S1 %s n=%d: %d replicate(s) failed", gen, n, sum(errs)))
  rows <- do.call(rbind, lapply(res, `[[`, "rows"))
  checks <- do.call(rbind, lapply(res, `[[`, "checks"))
  if (!is.null(checks)) {
    write.csv(checks, file.path(data_dir, sprintf("s1_screencheck_nca_%s_n%d.csv", gen, n)),
              row.names = FALSE)
  }
  cbind(study = "S1", gen = gen, n = n, rows)
}

run_s2 <- function(gen, n) {
  base_all <- read.csv(file.path(data_dir, sprintf("s2_base_%s_n%d.csv", gen, n)))
  reps <- sort(unique(base_all$rep))
  res <- mclapply(reps, function(rep) {
    b <- base_all[base_all$rep == rep, ]
    rows <- lapply(SIGMAS, function(s) {
      X <- pmin(pmax(b$X + s * b$zx, 0), 100)
      Y <- pmin(pmax(b$Y + s * b$zy, 0), 100)
      cbind(band = "noise", level = s, nca_rows(X, Y))
    })
    cbind(rep = rep, do.call(rbind, rows))
  }, mc.cores = cores)
  errs <- vapply(res, inherits, logical(1), "try-error")
  if (any(errs)) stop(sprintf("S2 %s n=%d: %d replicate(s) failed", gen, n, sum(errs)))
  cbind(study = "S2", gen = gen, n = n, do.call(rbind, res))
}

t0 <- Sys.time()
out <- list()
for (gen in c("G1", "G2")) {
  for (n in c(100, 300, 1000)) {
    out[[length(out) + 1L]] <- run_s1(gen, n)
    cat(sprintf("S1 %s n=%d done (%.0f s)\n", gen, n, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  for (n in c(100, 300)) {
    out[[length(out) + 1L]] <- run_s2(gen, n)
    cat(sprintf("S2 %s n=%d done (%.0f s)\n", gen, n, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}
all_rows <- do.call(rbind, out)
all_rows <- all_rows[, c("study", "gen", "n", "rep", "band", "level", "method", "d", "b50")]
write.csv(all_rows, file.path(out_dir, "robustness_nca.csv"), row.names = FALSE)
writeLines(c(
  sprintf("NCA = %s", as.character(packageVersion("NCA"))),
  sprintf("R = %s", R.version.string),
  sprintf("cores = %d", cores),
  sprintf("rows = %d", nrow(all_rows)),
  sprintf("seconds = %.0f", as.numeric(difftime(Sys.time(), t0, units = "secs")))
), file.path(out_dir, "robustness_nca_provenance.txt"))
cat(sprintf("wrote %d rows\n", nrow(all_rows)))
