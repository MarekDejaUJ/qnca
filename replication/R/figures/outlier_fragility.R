root <- normalizePath(getwd())
source(file.path(root, "replication", "R", "figures", "helpers.R"))
ensure_figure_dirs(root)

library(ggplot2)

d <- read.csv(file.path(root, "data", "outlier_diagnostic.csv"))
long <- rbind(
  data.frame(case = d$case, pi = "1.00", d_pi = d$d_pi_1.00),
  data.frame(case = d$case, pi = "0.95", d_pi = d$d_pi_0.95)
)
long$case <- factor(long$case, levels = c("clean", "with_efficient_outlier"),
                    labels = c("Clean sample", "Efficient case added"))

p <- ggplot(long, aes(x = case, y = d_pi, fill = pi)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, color = "white", linewidth = 0.2) +
  scale_fill_brewer(palette = "Set2", name = "Tolerance pi") +
  labs(
    x = NULL,
    y = "Effect size d_pi"
  ) +
  ylim(0, 0.6) +
  base_theme()

save_pdf(p, file.path(root, "replication", "figures", "qnca_outlier_fragility.pdf"), width = 3.4, height = 2.8)
