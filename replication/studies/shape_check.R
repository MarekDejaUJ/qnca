# Required condition at three outcome targets on the clean samples of S1,
# first 100 replicates per cell: NCA QR and CE-FDH ceilings, QNCA envelope,
# and the noise-free boundary. Output: results/shape_check.csv
suppressPackageStartupMessages(library(NCA))
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
root <- normalizePath(file.path(script_dir, ".."))
for (f in list.files(file.path(root, "..", "r", "R"), full.names = TRUE)) source(f)
S <- c(0, 100, 0, 100); YS <- c(20, 50, 80)
truth <- list(G1 = function(y) pmax(0, (y - 15) / 0.85), G2 = function(y) y^2 / 100)
rows <- list()
for (gen in c("G1", "G2")) for (n in c(300, 1000)) {
  cl <- read.csv(file.path(root, "results", "data", sprintf("s1_clean_%s_n%d.csv", gen, n)))
  for (rep in 1:100) {
    d <- cl[cl$rep == rep, c("X", "Y")]
    p <- nca_analysis(d, "X", "Y", ceilings = "qr", scope = S)$summaries$X$params
    a <- as.numeric(p["Intercept", "qr"]); b <- as.numeric(p["Slope", "qr"])
    qr_req <- if (is.finite(b) && b > 0) pmin(100, pmax(0, (YS - a) / b)) else rep(NA, 3)
    ce_req <- vapply(YS, function(y) min(d$X[d$Y >= y]), numeric(1))
    r95 <- qnca(d$X, d$Y, pi = 0.95, n_grid = 201L, B = NA, scope = S)$bottleneck
    q_req <- r95$X_required[match(YS, r95$outcome_level)]
    rows[[length(rows) + 1L]] <- data.frame(gen = gen, n = n, rep = rep, y = YS,
      truth = truth[[gen]](YS), ce_fdh = ce_req, qr = qr_req, qnca095 = q_req)
  }
}
out <- do.call(rbind, rows)
write.csv(out, file.path(root, "results", "shape_check.csv"), row.names = FALSE)
print(aggregate(cbind(truth, ce_fdh, qr, qnca095) ~ gen + n + y, data = out, FUN = mean), digits = 3)
