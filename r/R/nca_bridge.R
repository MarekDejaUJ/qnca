#' Compare QNCA(pi = 1) against Dul's NCA package on one data frame.
#'
#' Requires the \pkg{NCA} package. Runs \code{NCA::nca_analysis} with the CE-FDH
#' and CR-FDH ceilings on a fixed scope, extracts the effect size and p-value,
#' and returns them alongside QNCA's d_pi at pi = 1 and the native CE-FDH value.
#' The three CE-FDH numbers should agree up to grid discretization, which is the
#' validation invariant reported in the paper.
#'
#' @param d data frame with columns X and Y.
#' @param scope length-4 c(x_min, x_max, y_min, y_max).
#' @param test.rep permutation repetitions for NCA's own test.
#' @param n_grid QNCA frontier grid resolution.
#' @return one-row data frame of comparison statistics.
#' @export
compare_qnca_nca <- function(d, scope = c(0, 100, 0, 100), test.rep = 1999L,
                             n_grid = 200L) {
  if (!requireNamespace("NCA", quietly = TRUE)) {
    stop("R package 'NCA' is not installed. Run install.packages('NCA').")
  }
  model <- NCA::nca_analysis(
    data = d, x = "X", y = "Y",
    ceilings = c("ce_fdh", "cr_fdh"),
    scope = scope,
    bottleneck.x = "actual", bottleneck.y = "actual",
    steps = seq(scope[1], scope[2], length.out = 11L),
    test.rep = test.rep
  )
  nca_d <- tryCatch(
    NCA::nca_extract(model, x = "X", ceiling = "ce_fdh", param = "Effect size"),
    error = function(e) NA_real_)
  nca_p <- tryCatch(
    NCA::nca_extract(model, x = "X", ceiling = "ce_fdh", param = "p-value"),
    error = function(e) NA_real_)

  q1 <- qnca(d$X, d$Y, pi = 1, n_grid = n_grid, B = NA, scope = scope)
  native <- nca_ce_fdh_d(d$X, d$Y, scope = scope)

  data.frame(
    n = nrow(d),
    qnca_d_pi1 = q1$d_pi,
    native_ce_fdh_d = native,
    nca_ce_fdh_d = as.numeric(nca_d),
    nca_ce_fdh_p = as.numeric(nca_p),
    abs_diff_qnca_vs_nca = abs(q1$d_pi - as.numeric(nca_d)),
    abs_diff_native_vs_nca = abs(native - as.numeric(nca_d)),
    test.rep = test.rep,
    qnca_n_grid = n_grid,
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    nca_version = as.character(utils::packageVersion("NCA"))
  )
}
