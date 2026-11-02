#' qnca: Quantile Necessary Condition Analysis
#'
#' Quantile Necessary Condition Analysis (QNCA) reports, for each outcome
#' target, the required level of a condition: a low type-1 quantile of the
#' condition among the cases that attain the target. The tolerance \code{pi} in
#' (0, 1] sets the share \code{1 - pi} of attaining cases allowed below the
#' requirement. At \code{pi = 1} the frontier is the CE-FDH ceiling of
#' Necessary Condition Analysis (Dul 2016).
#'
#' \code{\link{qnca}} gives the effect size, the permutation p-value and the
#' bottleneck table at one tolerance. \code{\link{qnca_resolution}} gives the
#' target-by-target account of a frontier. \code{\link{qnca_outliers}} and
#' \code{\link{qnca_sensitivity}} give the deletion and addition screens.
#'
#' @references Dul, J. (2016). Necessary Condition Analysis (NCA): Logic and
#'   methodology of "necessary but not sufficient" causality. Organizational
#'   Research Methods, 19(1), 10-52. \doi{10.1177/1094428115584005}
#' @keywords internal
#' @importFrom stats runif rnorm
#' @importFrom utils head tail combn
"_PACKAGE"
