root <- normalizePath(getwd())
source(file.path(root, "replication", "R", "figures", "helpers.R"))
ensure_figure_dirs(root)
load_qnca(root)

library(ggplot2)

d <- read_dataset(root, 1000L)
pis <- c(1.00, 0.95, 0.90)
frontiers <- frontier_frame(d$X, d$Y, pis, n_grid = 101L)
write.csv(frontiers, file.path(root, "data", "figure_frontiers_n1000.csv"), row.names = FALSE)

poly <- empty_zone_polygon(frontiers, pi_value = "0.95", x_min = 0)
ce <- ce_fdh_ceiling_frame(d$X, d$Y)

p <- ggplot(d, aes(x = X, y = Y)) +
  geom_polygon(
    data = poly,
    aes(x = condition, y = outcome),
    inherit.aes = FALSE,
    fill = "#FDBF6F",
    alpha = 0.35
  ) +
  geom_point(size = 0.8, alpha = 0.35, color = "#4D4D4D") +
  geom_step(
    data = ce,
    aes(x = condition, y = outcome_ceiling),
    inherit.aes = FALSE,
    color = "#111111",
    linewidth = 0.45,
    linetype = "dashed"
  ) +
  geom_path(
    data = frontiers,
    aes(x = X_required, y = outcome_level, color = pi),
    linewidth = 0.85
  ) +
  scale_color_brewer(palette = "Dark2", name = "Tolerance pi") +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) +
  labs(
    x = "Condition X",
    y = "Outcome Y"
  ) +
  base_theme()

save_pdf(p, file.path(root, "replication", "figures", "qnca_frontier_geometry.pdf"), width = 3.4, height = 3.1)
