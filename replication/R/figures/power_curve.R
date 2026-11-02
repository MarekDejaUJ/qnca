root <- normalizePath(getwd())
source(file.path(root, "replication", "R", "figures", "helpers.R"))
ensure_figure_dirs(root)

library(ggplot2)

d <- read.csv(file.path(root, "data", "power_curve.csv"))
long <- rbind(
  data.frame(ceiling_strength = d$ceiling_strength, measure = "Permutation power", value = d$power),
  data.frame(ceiling_strength = d$ceiling_strength, measure = "Mean d_pi(0.95)", value = d$mean_d_pi_0.95)
)

p <- ggplot(long, aes(x = ceiling_strength, y = value, color = measure)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_wrap(~ measure, scales = "free_y", ncol = 1) +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  labs(
    x = "Ceiling strength",
    y = NULL
  ) +
  base_theme()

save_pdf(p, file.path(root, "replication", "figures", "qnca_power.pdf"), width = 3.4, height = 3.4)
