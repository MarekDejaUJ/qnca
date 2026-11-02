# Technology-acceptance application (settings in README.md).
#
# Rebuilds the rescaled latent scores with the published cIPMA workflow and, for
# the five direct antecedents of technology use, computes the NCA 5.0.2
# ceilings and outlier screen and the QNCA effects, permutation p-values,
# reference band, resolution table and deletion screen.
#
#   Rscript replication/studies/application.R

suppressPackageStartupMessages({
  library(seminr)
  library(cIPMA)
  library(NCA)
})

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
root <- normalizePath(file.path(script_dir, ".."))
out_dir <- file.path(root, "results")
for (f in list.files(file.path(root, "..", "r", "R"), full.names = TRUE)) source(f)
if (packageVersion("NCA") != "5.0.2") stop("NCA 5.0.2 is required")

PI_GRID <- c(1, 0.95, 0.90)
N_GRID <- 400L
B <- 1999L
M <- 999L
PAIRS <- c("AD", "EMV", "EOU", "PU", "CO")
PUBLISHED_D <- c(AD = 0.29, EMV = 0.33, EOU = 0.24, PU = 0.24, CO = 0.21)

# ---------------------------------------------------------------- scores

# Richter, Hauff, Kolev and Schubring (2023), Data in Brief 48, 109190;
# file "Extended TAM.csv" of Mendeley Data, doi:10.17632/pd5dp3phx2.4, CC BY 4.0.
tam_file <- file.path(root, "data", "tam.csv")
if (!file.exists(tam_file)) {
  utils::download.file(
    "https://data.mendeley.com/public-files/datasets/pd5dp3phx2/files/47b2cc21-ccde-48a1-a2af-9bcdad4e0e2f/file_downloaded",
    tam_file, mode = "wb", quiet = TRUE)
}
if (unname(tools::sha256sum(tam_file)) != "48930f478d415c2724135493c77b101a4aa1954795d8dcc035ec835bd600c7de") {
  stop("data/tam.csv differs from version 4 of the Mendeley Data file")
}
data <- read.csv(tam_file, fileEncoding = "UTF-8-BOM", check.names = FALSE)
mm <- constructs(
  reflective("PU", multi_items("PU_0", 1:3)),
  reflective("CO", multi_items("CO_0", 1:3)),
  reflective("EOU", multi_items("EOU_0", 1:3)),
  reflective("EMV", multi_items("EMV_0", 1:3)),
  reflective("AD", multi_items("AD_0", 1:3)),
  reflective("USE", single_item("USE_01"))
)
sm <- relationships(
  paths(from = "EOU", to = c("PU", "AD")),
  paths(from = c("PU", "EMV", "CO"), to = "AD"),
  paths(from = "AD", to = "USE")
)
model <- estimate_pls(data = data, measurement_model = mm, structural_model = sm)

unscaled_scores <- function(model, data) {
  w <- as.data.frame(as.table(model$outer_weights))
  colnames(w) <- c("Indicator", "Latent", "Std")
  w <- w[w$Std != 0, , drop = FALSE]
  w$Unstd <- w$Std / model$sdData[as.character(w$Indicator)]
  w$Norm <- w$Unstd / tapply(w$Unstd, w$Latent, sum)[as.character(w$Latent)]
  out <- data.frame(row.names = seq_len(nrow(data)))
  for (lv in unique(as.character(w$Latent))) {
    items <- model$mmMatrix[model$mmMatrix[, "construct"] == lv, "measurement"]
    out[[lv]] <- if (length(items) == 1L) data[[items]] else
      as.vector(as.matrix(data[, as.character(w$Indicator[w$Latent == lv]), drop = FALSE]) %*%
                  w$Norm[w$Latent == lv])
  }
  out
}
model$measurement_model <- model$mmMatrix
model$construct_scores <- as.matrix(unscaled_scores(model, data))
scales <- setNames(lapply(names(data), function(x) c(1, 5)), names(data))
scales[["USE_01"]] <- c(1, 7)
scores <- as.data.frame(cIPMA::get_rescaled_scores(model, scales))
for (nm in names(scores)) scores[[nm]] <- pmax(0, pmin(100, scores[[nm]]))
write.csv(scores, file.path(root, "data", "tam_scores_cipma.csv"), row.names = FALSE)

# ---------------------------------------------------------------- pairs

perm_matrix <- function(n, seed) {
  set.seed(seed)
  t(replicate(B, sample.int(n, n, replace = FALSE)))
}

rows <- list(); nca_out <- list(); qnca_out <- list(); band_out <- list(); res_out <- list()
for (idx in seq_along(PAIRS)) {
  xname <- PAIRS[idx]
  X <- scores[[xname]]; Y <- scores[["USE"]]
  scope <- c(min(X), max(X), min(Y), max(Y))
  d_frame <- data.frame(X = X, Y = Y)

  nm <- nca_analysis(d_frame, "X", "Y", ceilings = c("ce_fdh", "cr_fdh", "qr"))
  p <- nm$summaries$X$params
  eff <- function(ceil) as.numeric(p["Effect size", ceil])

  perms <- perm_matrix(length(Y), 20260930L + idx)
  q_env <- lapply(PI_GRID, function(pi) qnca(X, Y, pi = pi, n_grid = N_GRID, B = NA,
                                               scope = scope, permutations = perms))
  q_ls <- lapply(PI_GRID, function(pi) qnca(X, Y, pi = pi, n_grid = N_GRID, B = NA,
                                              scope = scope, monotone = "isotonic"))
  band <- spuriousness_band(X, Y, pi_grid = PI_GRID, n_grid = N_GRID, M = M,
                            scope = scope, seed = 20261030L + idx)
  b <- band$band
  b95 <- b[b$pi == 0.95, ]

  rows[[idx]] <- data.frame(
    pair = paste0(xname, " -> USE"),
    published_d = PUBLISHED_D[[xname]],
    rank_corr = band$rank_correlation,
    nca_ce_fdh = eff("ce_fdh"), nca_cr_fdh = eff("cr_fdh"), nca_qr = eff("qr"),
    d_100 = q_env[[1]]$d_pi, d_095 = q_env[[2]]$d_pi, d_090 = q_env[[3]]$d_pi,
    p_100 = q_env[[1]]$p_pi, p_095 = q_env[[2]]$p_pi, p_090 = q_env[[3]]$p_pi,
    d_ls_095 = q_ls[[2]]$d_pi, d_ls_090 = q_ls[[3]]$d_pi,
    U_095 = b95$upper, delta_095 = b95$excess, p_null_095 = b95$p_null,
    x_min = scope[1], x_max = scope[2], y_min = scope[3], y_max = scope[4]
  )
  band_out[[idx]] <- cbind(pair = xname, b)

  res <- qnca_resolution(X, Y, pi = 0.95, n_grid = N_GRID, scope = scope)
  res_out[[idx]] <- cbind(pair = xname, res)

  for (k in 1:2) {
    o <- nca_outliers(d_frame, "X", "Y", ceiling = "ce_fdh", k = k)
    if (!is.null(o) && nrow(o) > 0L) {
      nca_out[[length(nca_out) + 1L]] <- data.frame(
        pair = xname, k = k, rank = seq_len(nrow(o)), cases = as.character(o$outliers),
        eff_or = o$eff.or, eff_nw = o$eff.nw, dif_abs = o$dif.abs, dif_rel = o$dif.rel)
    }
    qo <- qnca_outliers(X, Y, pi_grid = PI_GRID, k = k, n_grid = N_GRID, scope = scope,
                        max_candidates = 12L)
    top <- qo$table[seq_len(min(10L, nrow(qo$table))), ]
    qnca_out[[length(qnca_out) + 1L]] <- cbind(pair = xname, k = k, rank = seq_len(nrow(top)), top,
                                               max_abs_100 = max(abs(qo$dif_abs[, 1])),
                                               max_abs_095 = max(abs(qo$dif_abs[, 2])),
                                               max_abs_090 = max(abs(qo$dif_abs[, 3])),
                                               d_full_100 = qo$d[1], d_full_095 = qo$d[2],
                                               d_full_090 = qo$d[3])
  }
  cat(sprintf("%s done\n", xname))
}

write.csv(do.call(rbind, rows), file.path(out_dir, "application_table.csv"), row.names = FALSE)
write.csv(do.call(rbind, band_out), file.path(out_dir, "application_band.csv"), row.names = FALSE)
write.csv(do.call(rbind, res_out), file.path(out_dir, "application_resolution.csv"), row.names = FALSE)
write.csv(do.call(rbind, nca_out), file.path(out_dir, "application_outliers_nca.csv"), row.names = FALSE)
write.csv(do.call(rbind, qnca_out), file.path(out_dir, "application_outliers_qnca.csv"), row.names = FALSE)
writeLines(c(
  sprintf("NCA = %s", as.character(packageVersion("NCA"))),
  sprintf("seminr = %s", as.character(packageVersion("seminr"))),
  sprintf("cIPMA = %s", as.character(packageVersion("cIPMA"))),
  sprintf("R = %s", R.version.string),
  sprintf("tam.csv sha256 = %s", "48930f478d415c2724135493c77b101a4aa1954795d8dcc035ec835bd600c7de")
), file.path(out_dir, "application_provenance.txt"))
