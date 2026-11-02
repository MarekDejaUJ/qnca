# Figures for the planted-case example and the robustness comparison.
#
#   Rscript replication/studies/figures.R <output-directory>

suppressPackageStartupMessages(library(ggplot2))

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
root <- normalizePath(file.path(script_dir, ".."))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else file.path(root, "results", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (f in list.files(file.path(root, "..", "r", "R"), full.names = TRUE)) source(f)

base_theme <- function() {
  theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          panel.spacing = grid::unit(1.4, "lines"),
          strip.text = element_text(face = "bold"))
}

save_pdf <- function(plot, path, width, height) {
  if (file.exists(path)) unlink(path)
  ok <- tryCatch({
    suppressWarnings(ggsave(filename = path, plot = plot, device = grDevices::cairo_pdf,
                            width = width, height = height, units = "in"))
    file.exists(path) && file.info(path)$size > 0
  }, error = function(e) FALSE)
  if (!ok) ggsave(filename = path, plot = plot, device = "pdf",
                  width = width, height = height, units = "in")
}

# ---------------------------------------------------------------- planted case

ex <- read.csv(file.path(root, "data", "example1.csv"))
SCOPE <- c(0, 100, 0, 100)
panels <- list("Clean sample" = ex[ex$planted == 0, ], "One efficient case added" = ex)
frontier_rows <- list(); point_rows <- list()
for (nm in names(panels)) {
  d <- panels[[nm]]
  spec <- list(
    list(label = "CE-FDH", pi = 1, monotone = "envelope"),
    list(label = "QNCA 0.95, envelope", pi = 0.95, monotone = "envelope"),
    list(label = "QNCA 0.95, least squares", pi = 0.95, monotone = "isotonic")
  )
  for (s in spec) {
    r <- qnca(d$X, d$Y, pi = s$pi, n_grid = 200L, B = NA, scope = SCOPE, monotone = s$monotone)
    frontier_rows[[length(frontier_rows) + 1L]] <- data.frame(
      panel = nm, frontier = s$label,
      outcome = r$bottleneck$outcome_level, condition = r$bottleneck$X_required)
  }
  point_rows[[nm]] <- cbind(panel = nm, d)
}
fr <- do.call(rbind, frontier_rows)
pts <- do.call(rbind, point_rows)
fr$panel <- factor(fr$panel, levels = names(panels))
pts$panel <- factor(pts$panel, levels = names(panels))
fr$frontier <- factor(fr$frontier, levels = c("CE-FDH", "QNCA 0.95, envelope", "QNCA 0.95, least squares"))
write.csv(fr, file.path(out_dir, "planted_case_frontiers.csv"), row.names = FALSE)

p1 <- ggplot() +
  geom_point(data = pts[pts$planted == 0, ], aes(X, Y), size = 0.6, alpha = 0.35, colour = "#4D4D4D") +
  geom_point(data = pts[pts$planted == 1, ], aes(X, Y), shape = 17, size = 2.2, colour = "#111111") +
  geom_path(data = fr, aes(condition, outcome, colour = frontier, linetype = frontier), linewidth = 0.7) +
  scale_colour_manual(name = NULL, values = c("CE-FDH" = "#111111",
                                              "QNCA 0.95, envelope" = "#1B9E77",
                                              "QNCA 0.95, least squares" = "#D95F02")) +
  scale_linetype_manual(name = NULL, values = c("CE-FDH" = "dashed",
                                                "QNCA 0.95, envelope" = "solid",
                                                "QNCA 0.95, least squares" = "dotted")) +
  facet_wrap(~ panel, nrow = 1) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) +
  labs(x = "Condition X", y = "Outcome Y") +
  base_theme()
save_pdf(p1, file.path(out_dir, "qnca_planted_case.pdf"), width = 6.2, height = 3.3)

# ---------------------------------------------------------------- retention

ret_file <- file.path(root, "results", "summary_retention.csv")
if (!file.exists(ret_file)) {
  cat("summary_retention.csv not found; retention figure skipped\n")
  quit(save = "no")
}
ret <- read.csv(ret_file)
keep <- c("CE-FDH" = "CE-FDH", "CR-FDH" = "CR-FDH", "QR" = "QR ceiling",
          "Screen 1% rule" = "Screen and delete",
          "Order-m 19" = "Order-m, m = 19",
          "QNCA least squares 0.95" = "QNCA 0.95, least squares",
          "QNCA envelope 0.95" = "QNCA 0.95, envelope",
          "QNCA envelope 0.9" = "QNCA 0.90, envelope")
ret <- ret[ret$n == 300 & ret$label %in% names(keep), ]
clean_rows <- ret[ret$band == "clean", ]
ret <- ret[ret$band != "clean", ]
zero <- unique(ret[, c("gen", "band", "label")]); zero$level <- 0; zero$mean_retention <- 1
for (i in seq_len(nrow(zero))) {
  hit <- clean_rows$gen == zero$gen[i] & clean_rows$label == zero$label[i]
  if (any(hit)) zero$mean_retention[i] <- clean_rows$mean_retention[hit][1]
}
ret <- rbind(ret[, c("gen", "band", "label", "level", "mean_retention")],
             zero[, c("gen", "band", "label", "level", "mean_retention")])
ret$method <- factor(keep[ret$label], levels = unname(keep))
ret$panel <- factor(paste(ifelse(ret$gen == "G1", "Linear ceiling", "Concave ceiling"),
                          ifelse(ret$band == "mid", "middle band", "top band"), sep = ", "),
                    levels = c("Linear ceiling, middle band", "Linear ceiling, top band",
                               "Concave ceiling, middle band", "Concave ceiling, top band"))
cols <- c("CE-FDH" = "#111111", "CR-FDH" = "#7F7F7F", "QR ceiling" = "#7570B3",
          "Screen and delete" = "#E7298A", "Order-m, m = 19" = "#A6761D",
          "QNCA 0.95, least squares" = "#D95F02", "QNCA 0.95, envelope" = "#1B9E77",
          "QNCA 0.90, envelope" = "#66A61E")
ltys <- c("CE-FDH" = "dashed", "CR-FDH" = "dotdash", "QR ceiling" = "longdash",
          "Screen and delete" = "twodash", "Order-m, m = 19" = "dotted",
          "QNCA 0.95, least squares" = "dotted", "QNCA 0.95, envelope" = "solid",
          "QNCA 0.90, envelope" = "solid")
p2 <- ggplot(ret, aes(level, mean_retention, colour = method, linetype = method)) +
  geom_hline(yintercept = 1, colour = "#BBBBBB", linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.1) +
  scale_colour_manual(name = NULL, values = cols) +
  scale_linetype_manual(name = NULL, values = ltys) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 5)) +
  facet_wrap(~ panel, nrow = 2) +
  coord_cartesian(ylim = c(0, 1.1)) +
  labs(x = "Added efficient cases", y = "Mean share of the clean effect size") +
  base_theme() +
  guides(colour = guide_legend(nrow = 3), linetype = guide_legend(nrow = 3))
save_pdf(p2, file.path(out_dir, "qnca_robustness.pdf"), width = 6.2, height = 5.2)
cat("figures written to", out_dir, "\n")
