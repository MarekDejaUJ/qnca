root <- normalizePath(getwd())
source(file.path(root, "replication", "R", "figures", "helpers.R"))
ensure_figure_dirs(root)
load_qnca(root)

library(ggplot2)

d <- read_dataset(root, 1000L)
frontier <- frontier_frame(d$X, d$Y, c(1.00), n_grid = 101L)
ce <- ce_fdh_ceiling_frame(d$X, d$Y)
poly <- empty_zone_polygon(frontier, pi_value = "1.00", x_min = 0)

p <- ggplot(d, aes(x = X, y = Y)) +
  geom_polygon(
    data = poly,
    aes(x = condition, y = outcome),
    inherit.aes = FALSE,
    fill = "#B2DF8A",
    alpha = 0.28
  ) +
  geom_point(size = 0.7, alpha = 0.28, color = "#4D4D4D") +
  geom_step(
    data = ce,
    aes(x = condition, y = outcome_ceiling, linetype = "Dul CE-FDH"),
    inherit.aes = FALSE,
    color = "#111111",
    linewidth = 0.55
  ) +
  geom_step(
    data = frontier,
    aes(x = X_required, y = outcome_level, linetype = "QNCA pi = 1"),
    color = "#1F78B4",
    linewidth = 0.65
  ) +
  scale_linetype_manual(name = NULL, values = c("Dul CE-FDH" = "dashed", "QNCA pi = 1" = "solid")) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) +
  labs(
    x = "Condition X",
    y = "Outcome Y"
  ) +
  base_theme()

save_pdf(p, file.path(root, "replication", "figures", "qnca_pi_reduction.pdf"), width = 3.4, height = 3.1)
