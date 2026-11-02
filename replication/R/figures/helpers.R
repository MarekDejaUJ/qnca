ensure_figure_dirs <- function(root) {
  dir.create(file.path(root, "replication", "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "data"), recursive = TRUE, showWarnings = FALSE)
}

load_qnca <- function(root) {
  source(file.path(root, "r", "R", "qnca.R"))
}

save_pdf <- function(plot, path, width = 3.4, height = 2.8) {
  try(
    ggplot2::ggsave(
      filename = path,
      plot = plot,
      device = grDevices::cairo_pdf,
      width = width,
      height = height,
      units = "in"
    ),
    silent = TRUE
  )
  if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")
  if (!file.exists(path) || file.info(path)$size == 0L) {
    ggplot2::ggsave(
      filename = path,
      plot = plot,
      device = "pdf",
      width = width,
      height = height,
      units = "in"
    )
  }
}

frontier_frame <- function(X, Y, pis, n_grid = 101L, scope = c(0, 100, 0, 100)) {
  rows <- list()
  idx <- 1L
  for (pi_value in pis) {
    res <- qnca(X, Y, pi = pi_value, n_grid = n_grid, B = NA, scope = scope)
    rows[[idx]] <- data.frame(
      pi = sprintf("%.2f", pi_value),
      outcome_level = res$bottleneck$outcome_level,
      X_required = res$bottleneck$X_required
    )
    idx <- idx + 1L
  }
  do.call(rbind, rows)
}

ce_fdh_ceiling_frame <- function(X, Y, n_grid = 400L, scope = c(0, 100, 0, 100)) {
  x_grid <- seq(scope[1], scope[2], length.out = n_grid)
  y_ceiling <- vapply(x_grid, function(x) {
    hit <- X <= x
    if (!any(hit)) return(scope[3])
    max(Y[hit])
  }, numeric(1))
  data.frame(condition = x_grid, outcome_ceiling = y_ceiling)
}

empty_zone_polygon <- function(frontier, pi_value = "0.95", x_min = 0) {
  f <- frontier[frontier$pi == pi_value, ]
  rbind(
    data.frame(condition = x_min, outcome = f$outcome_level),
    data.frame(condition = rev(f$X_required), outcome = rev(f$outcome_level))
  )
}

read_dataset <- function(root, n = 1000L) {
  read.csv(file.path(root, "data", paste0("qnca_dataset_n", n, ".csv")))
}

base_theme <- function() {
  ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}
