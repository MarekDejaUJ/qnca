#' Quantile Necessary Condition Analysis (QNCA)
#'
#' QNCA generalises Dul's Necessary Condition Analysis. The deterministic NCA
#' ceiling is the \code{pi = 1} limit of the QNCA frontier family; lowering
#' \code{pi} discards the most efficient cases and raises the necessity floor,
#' trading a controlled amount of crispness for robustness to outliers.
#'
#' This file is intentionally dependency-free (base R + \pkg{stats} only) so it
#' agrees, to floating point, with the Python and Julia siblings. Tolerances are
#' read as exact decimals before the type-1 rank is formed.

#' Exact decimal parts of a tolerance.
#'
#' Returns the tail probability \code{1 - pi} as a reduced fraction
#' \code{num / den} of whole numbers, read from the shortest decimal that
#' round-trips to \code{pi}. Tolerances may carry at most nine decimal places.
qnca_tail_fraction <- function(pi) {
  if (length(pi) != 1L || !is.numeric(pi) || !is.finite(pi)) {
    stop("pi must be a single finite number")
  }
  if (!(pi > 0 && pi <= 1)) stop("pi must be in (0, 1]")
  txt <- NA_character_
  for (digits in 1:17) {
    cand <- sprintf("%.*g", digits, pi)
    if (as.numeric(cand) == pi) { txt <- cand; break }
  }
  parts <- regmatches(txt, regexec("^([0-9]+)(\\.([0-9]*))?([eE]([+-]?[0-9]+))?$", txt))[[1]]
  if (length(parts) == 0L) stop("pi must have a decimal representation")
  int_part <- parts[2]; frac_part <- parts[4]
  expo <- if (nzchar(parts[6])) as.integer(parts[6]) else 0L
  scale <- nchar(frac_part) - expo
  digits_txt <- sub("^0+(?=.)", "", paste0(int_part, frac_part), perl = TRUE)
  if (scale <= 0L) {
    # A whole-number tolerance inside (0, 1] can only be 1.
    return(c(num = 0, den = 1))
  }
  if (scale > 9L) stop("pi may carry at most nine decimal places")
  den <- 10^scale
  num <- den - as.numeric(digits_txt)
  if (num == 0) return(c(num = 0, den = 1))
  a <- num; b <- den
  while (b > 0) { r <- a %% b; a <- b; b <- r }
  c(num = num / a, den = den / a)
}

#' Selected type-1 order-statistic rank \code{max(1, ceiling(k * (1 - pi)))}.
#'
#' The ceiling is evaluated in exact integer arithmetic from the decimal value
#' of \code{pi}, so \code{pi = 0.95} and \code{k = 100} select rank 5.
#' @param k conditioning-set size (positive whole number).
#' @param pi tolerance in (0, 1].
#' @export
qnca_rank <- function(k, pi) {
  if (length(k) != 1L || !is.finite(k) || k < 1 || k != floor(k)) {
    stop("conditioning-set size must be a positive whole number")
  }
  if (k >= 2^31) stop("conditioning-set size is too large")
  qnca_rank_from_fraction(k, qnca_tail_fraction(pi))
}

qnca_rank_from_fraction <- function(k, fr) {
  a <- fr[["num"]]; D <- fr[["den"]]
  if (a == 0) return(1L)
  # Exact floor(k * a / D) with k < 2^31 and a < D <= 1e9, split so that every
  # intermediate product stays below 2^53.
  k_hi <- k %/% 65536; k_lo <- k %% 65536
  t1 <- k_hi * a
  q1 <- t1 %/% D; r1 <- t1 %% D
  t2 <- r1 * 65536 + k_lo * a
  q2 <- t2 %/% D; r2 <- t2 %% D
  rank <- q1 * 65536 + q2 + (r2 > 0)
  as.integer(min(k, max(1, rank)))
}

#' Type-1 quantile parameterised by the QNCA tolerance \code{pi}.
#' @export
qnca_quantile_type1_pi <- function(values, pi) {
  if (length(values) == 0L) return(NA_real_)
  x <- sort(as.numeric(values))
  x[qnca_rank(length(x), pi)]
}

qnca_validate_xy <- function(X, Y) {
  if (length(X) != length(Y)) stop("X and Y must have equal lengths")
  if (length(X) < 2L) stop("X and Y must contain at least two observations")
  if (!all(is.finite(X))) stop("X must contain only finite values")
  if (!all(is.finite(Y))) stop("Y must contain only finite values")
  if (!(min(X) < max(X))) stop("X must be non-degenerate")
  if (!(min(Y) < max(Y))) stop("Y must be non-degenerate")
  invisible(NULL)
}

qnca_validated_scope <- function(X, Y, scope = NULL) {
  bounds <- if (is.null(scope)) c(min(X), max(X), min(Y), max(Y)) else as.numeric(scope)
  if (length(bounds) != 4L) stop("scope must contain x_min, x_max, y_min, y_max")
  if (!all(is.finite(bounds))) stop("scope bounds must be finite")
  if (!(bounds[1] < bounds[2])) stop("scope requires x_min < x_max")
  if (!(bounds[3] < bounds[4])) stop("scope requires y_min < y_max")
  if (min(X) < bounds[1] || max(X) > bounds[2]) stop("X observations must lie inside scope")
  if (min(Y) < bounds[3] || max(Y) > bounds[4]) stop("Y observations must lie inside scope")
  bounds
}

qnca_validate_grid <- function(y_grid) {
  if (length(y_grid) == 0L) stop("y_grid must not be empty")
  if (!all(is.finite(y_grid))) stop("y_grid must contain only finite values")
  if (is.unsorted(y_grid)) stop("y_grid must be non-decreasing")
  invisible(NULL)
}

#' Quantile necessity frontier phi_pi(y), isotonic-fit and carried forward to the
#' scope ceiling so a fixed scope reproduces NCA's effect size.
#'
#' @param X,Y numeric vectors (condition, outcome).
#' @param pi tolerance in (0, 1].
#' @param y_grid numeric grid of outcome levels.
#' @param x_max optional fixed-scope condition ceiling.
#' @return numeric vector of required-X values, one per grid level.
#' @export
qnca_frontier <- function(X, Y, pi, y_grid, x_max = NULL) {
  X <- as.numeric(X); Y <- as.numeric(Y); y_grid <- as.numeric(y_grid)
  qnca_validate_xy(X, Y)
  qnca_validate_grid(y_grid)
  qnca_tail_fraction(pi)
  qnca_frontier_validated(X, Y, pi, y_grid, x_max = x_max)
}

qnca_frontier_validated <- function(X, Y, pi, y_grid, x_max = NULL) {
  fr <- qnca_tail_fraction(pi)
  phi <- vapply(y_grid, function(y) {
    hit <- Y >= y
    if (!any(hit)) return(NA_real_)
    x <- sort(X[hit])
    x[qnca_rank_from_fraction(length(x), fr)]
  }, numeric(1))
  fin <- which(is.finite(phi))
  if (length(fin) >= 2L) {
    phi[fin] <- stats::isoreg(y_grid[fin], phi[fin])$yf
  }
  # Under a fixed scope Dul's CE-FDH counts the whole high-outcome band as empty.
  # Frontier-only calls without x_max keep the last finite value.
  if (length(fin) >= 1L) {
    last <- fin[length(fin)]
    if (last < length(phi)) {
      phi[(last + 1L):length(phi)] <- if (is.null(x_max)) phi[last] else x_max
    }
  }
  phi
}

#' Effect size d_pi: trapezoidal empty-zone area left of the frontier / scope.
#' @export
qnca_d <- function(phi, y_grid, x_min, scope) {
  w <- pmax(0, phi - x_min)
  ok <- which(is.finite(w))
  if (length(ok) < 2L || scope <= 0) return(0)
  yg <- y_grid[ok]; wv <- w[ok]
  area <- sum(diff(yg) * (utils::head(wv, -1) + utils::tail(wv, -1)) / 2)
  max(0, min(1, area / scope))
}

qnca_validate_permutations <- function(permutations, n) {
  if (is.vector(permutations) && !is.matrix(permutations)) {
    perm <- matrix(permutations, nrow = 1L)
  } else {
    perm <- as.matrix(permutations)
  }
  storage.mode(perm) <- "integer"
  if (ncol(perm) != n) stop("permutations must be a B x n matrix")
  expected <- seq_len(n)
  for (i in seq_len(nrow(perm))) {
    if (!all(sort(as.integer(perm[i, ])) == expected)) {
      stop("each permutation row must contain 1:n exactly once")
    }
  }
  perm
}

#' Fit QNCA at a single tolerance.
#'
#' @param X,Y numeric vectors.
#' @param pi tolerance in (0, 1]; pi = 1 reproduces Dul's CE-FDH ceiling.
#' @param n_grid number of outcome grid levels.
#' @param B permutation repetitions; set 0 (or NA) to skip the p-value.
#' @param scope optional length-4 vector c(x_min, x_max, y_min, y_max).
#' @param permutations optional one-based B x n permutation matrix.
#' @return list with pi, d_pi, p_pi, scope and the bottleneck data frame.
#' @export
qnca <- function(X, Y, pi = 1, n_grid = 50L, B = 1999L, scope = NULL,
                 permutations = NULL) {
  X <- as.numeric(X); Y <- as.numeric(Y)
  qnca_validate_xy(X, Y)
  if (length(n_grid) != 1L || !is.finite(n_grid) || n_grid < 2) stop("n_grid must be at least 2")
  if (length(B) != 1L || (!is.na(B) && B < 0)) stop("B must be non-negative or NA")
  qnca_tail_fraction(pi)
  bounds <- qnca_validated_scope(X, Y, scope)
  x_min <- bounds[1]; x_max <- bounds[2]; y_min <- bounds[3]; y_max <- bounds[4]
  scope_area <- (x_max - x_min) * (y_max - y_min)
  y_grid <- seq(y_min, y_max, length.out = n_grid)
  phi_obs <- qnca_frontier_validated(X, Y, pi, y_grid, x_max = x_max)
  d_obs <- qnca_d(phi_obs, y_grid, x_min, scope_area)

  p_pi <- NA_real_
  if (!is.null(permutations)) {
    perm <- qnca_validate_permutations(permutations, length(Y))
    count <- 0L
    for (b in seq_len(nrow(perm))) {
      phi_p <- qnca_frontier_validated(X, Y[perm[b, ]], pi, y_grid, x_max = x_max)
      d_p <- qnca_d(phi_p, y_grid, x_min, scope_area)
      if (is.finite(d_p) && d_p >= d_obs) count <- count + 1L
    }
    p_pi <- (1 + count) / (nrow(perm) + 1)
  } else if (!is.na(B) && B > 0L) {
    count <- 0L
    for (b in seq_len(B)) {
      phi_p <- qnca_frontier_validated(X, sample(Y), pi, y_grid, x_max = x_max)
      d_p <- qnca_d(phi_p, y_grid, x_min, scope_area)
      if (is.finite(d_p) && d_p >= d_obs) count <- count + 1L
    }
    p_pi <- (1 + count) / (B + 1)
  }

  list(pi = pi, d_pi = d_obs, p_pi = p_pi, scope = scope_area,
       bottleneck = data.frame(outcome_level = y_grid, X_required = phi_obs))
}

#' Native CE-FDH effect size, computed the way Dul's NCA defines it.
#'
#' Builds the upper-left free-disposal hull directly from the data and integrates
#' the empty zone exactly. Lets the validation invariant d_pi(pi = 1) == this be
#' checked without the \pkg{NCA} package. Equivalent to \code{qnca}'s d at pi = 1
#' up to grid discretization.
#' @export
nca_ce_fdh_d <- function(X, Y, scope = NULL) {
  X <- as.numeric(X); Y <- as.numeric(Y)
  qnca_validate_xy(X, Y)
  bounds <- qnca_validated_scope(X, Y, scope)
  x_min <- bounds[1]; x_max <- bounds[2]; y_min <- bounds[3]; y_max <- bounds[4]
  scope_area <- (x_max - x_min) * (y_max - y_min)
  if (scope_area <= 0) return(0)

  ord <- order(-Y)
  Ys <- Y[ord]; Xs <- X[ord]
  run_min <- cummin(Xs)               # run_min[j] = min(Xs[1..j])
  clamp_y <- function(v) min(max(v, y_min), y_max)

  area <- 0
  top_lo <- clamp_y(Ys[1])
  if (y_max > top_lo) area <- area + (y_max - top_lo) * max(0, x_max - x_min)
  if (length(Ys) >= 2L) {
    for (k in 2:length(Ys)) {
      hi <- clamp_y(Ys[k - 1]); lo <- clamp_y(Ys[k])
      if (hi > lo) area <- area + (hi - lo) * max(0, run_min[k - 1] - x_min)
    }
  }
  bot_hi <- clamp_y(Ys[length(Ys)])
  if (bot_hi > y_min) area <- area + (bot_hi - y_min) * max(0, run_min[length(run_min)] - x_min)

  max(0, min(1, area / scope_area))
}

qnca_scope_tuple <- function(X, Y, scope = NULL) {
  qnca_validated_scope(X, Y, scope)
}

qnca_d_for_pi <- function(X, Y, pi, scope, n_grid) {
  qnca(X, Y, pi = pi, n_grid = n_grid, B = NA, scope = scope)$d_pi
}

qnca_empirical_quantile_values <- function(values, probs) {
  x <- sort(as.numeric(values))
  p <- pmin(1, pmax(0, as.numeric(probs)))
  idx <- ceiling(length(x) * p)
  idx <- pmin(length(x), pmax(1L, idx))
  x[idx]
}

qnca_quantile_type7 <- function(values, prob) {
  x <- sort(as.numeric(values))
  if (length(x) == 0L) return(NA_real_)
  if (length(x) == 1L || prob <= 0) return(x[1])
  if (prob >= 1) return(x[length(x)])
  h <- 1 + (length(x) - 1) * prob
  j <- floor(h)
  gamma <- h - j
  if (j >= length(x)) return(x[length(x)])
  (1 - gamma) * x[j] + gamma * x[j + 1L]
}

qnca_rank_normal_correlation <- function(X, Y) {
  n <- length(X)
  ax <- stats::qnorm(rank(X, ties.method = "average") / (n + 1))
  ay <- stats::qnorm(rank(Y, ties.method = "average") / (n + 1))
  if (stats::sd(ax) == 0 || stats::sd(ay) == 0) return(0)
  r <- suppressWarnings(stats::cor(ax, ay))
  if (!is.finite(r)) return(0)
  max(-0.999999, min(0.999999, r))
}

#' Proposed pi-resolved spuriousness band.
#'
#' Compares the observed d_pi curve with a Gaussian-copula null that preserves
#' the empirical marginals and observed normal-score rank correlation. This is a
#' diagnostic with unestablished operating characteristics; it does not modify
#' the QNCA estimator.
#'
#' @param X,Y numeric vectors.
#' @param pi_grid tolerance values.
#' @param n_grid number of outcome grid levels.
#' @param M number of null datasets.
#' @param scope optional length-4 vector c(x_min, x_max, y_min, y_max).
#' @param seed optional random seed.
#' @param normal_draws optional normal array with dimensions M x n x 2.
#' @param uniform_draws optional uniform array with dimensions M x n x 2.
#' @return list with a band data frame, null effect matrix and rank correlation.
#' @export
spuriousness_band <- function(X, Y, pi_grid = c(1, 0.95, 0.90), n_grid = 50L,
                              M = 999L, scope = NULL, seed = NULL,
                              normal_draws = NULL, uniform_draws = NULL) {
  X <- as.numeric(X); Y <- as.numeric(Y); pi_grid <- as.numeric(pi_grid)
  qnca_validate_xy(X, Y)
  if (length(n_grid) != 1L || !is.finite(n_grid) || n_grid < 2) stop("n_grid must be at least 2")
  if (length(pi_grid) == 0L) stop("pi_grid must not be empty")
  for (pi in pi_grid) qnca_tail_fraction(pi)
  scope_tuple <- qnca_scope_tuple(X, Y, scope)
  n <- length(X)
  if (!is.null(normal_draws) && !is.null(uniform_draws)) {
    stop("supply at most one of normal_draws and uniform_draws")
  }
  if (!is.null(uniform_draws)) {
    uniforms <- uniform_draws
    if (length(dim(uniforms)) != 3L || dim(uniforms)[2L] != n || dim(uniforms)[3L] != 2L) {
      stop("uniform_draws must have dimensions M x n x 2")
    }
    M <- dim(uniforms)[1L]
    draws <- NULL
  } else if (is.null(normal_draws)) {
    if (!is.null(seed)) set.seed(seed)
    draws <- array(stats::rnorm(M * n * 2L), dim = c(M, n, 2L))
    uniforms <- NULL
  } else {
    draws <- normal_draws
    if (length(dim(draws)) != 3L || dim(draws)[2L] != n || dim(draws)[3L] != 2L) {
      stop("normal_draws must have dimensions M x n x 2")
    }
    M <- dim(draws)[1L]
    uniforms <- NULL
  }
  if (M <= 0L) stop("M must be positive")

  d_obs <- vapply(pi_grid, function(pi) {
    qnca_d_for_pi(X, Y, pi, scope_tuple, n_grid)
  }, numeric(1))
  d_null <- matrix(NA_real_, nrow = M, ncol = length(pi_grid))
  r <- qnca_rank_normal_correlation(X, Y)
  scale <- sqrt(max(0, 1 - r * r))

  for (j in seq_len(M)) {
    if (!is.null(uniforms)) {
      Xj <- qnca_empirical_quantile_values(X, uniforms[j, , 1L])
      Yj <- qnca_empirical_quantile_values(Y, uniforms[j, , 2L])
    } else {
      zx <- draws[j, , 1L]
      zy <- r * zx + scale * draws[j, , 2L]
      Xj <- qnca_empirical_quantile_values(X, stats::pnorm(zx))
      Yj <- qnca_empirical_quantile_values(Y, stats::pnorm(zy))
    }
    for (k in seq_along(pi_grid)) {
      d_null[j, k] <- qnca_d_for_pi(Xj, Yj, pi_grid[k], scope_tuple, n_grid)
    }
  }

  lower <- vapply(seq_along(pi_grid), function(k) qnca_quantile_type7(d_null[, k], 0.025), numeric(1))
  median <- vapply(seq_along(pi_grid), function(k) qnca_quantile_type7(d_null[, k], 0.500), numeric(1))
  upper <- vapply(seq_along(pi_grid), function(k) qnca_quantile_type7(d_null[, k], 0.975), numeric(1))
  excess <- d_obs - upper
  p_null <- (1 + colSums(sweep(d_null, 2L, d_obs, FUN = ">="))) / (M + 1)

  list(
    band = data.frame(pi = pi_grid, d_obs = d_obs, lower = lower,
                      median = median, upper = upper, excess = excess,
                      p_null = p_null),
    d_null = d_null,
    rank_correlation = r
  )
}

#' Proposed extreme-vs-consistent consistency probe.
#'
#' Tracks d_pi across increasing subsample sizes for pi = 1 and a lower
#' tolerance. This is a descriptive diagnostic and can be confounded by hard
#' support bounds or marginal skewness.
#'
#' @param X,Y numeric vectors.
#' @param pi_pair two tolerance values, usually c(1, 0.90).
#' @param sizes integer subsample sizes.
#' @param reps number of subsamples per size.
#' @param n_grid number of outcome grid levels.
#' @param scope optional length-4 vector c(x_min, x_max, y_min, y_max).
#' @param seed optional random seed.
#' @param subsamples optional list of one-based reps x size matrices.
#' @return list with per-size summaries, effects, drift slopes and divergence.
#' @export
consistency_probe <- function(X, Y, pi_pair = c(1, 0.90), sizes = NULL,
                              reps = 100L, n_grid = 50L, scope = NULL,
                              seed = NULL, subsamples = NULL) {
  X <- as.numeric(X); Y <- as.numeric(Y); pi_pair <- as.numeric(pi_pair)
  if (length(pi_pair) != 2L) stop("pi_pair must contain exactly two tolerances")
  n <- length(X)
  if (is.null(sizes)) {
    sizes <- unique(pmin(n, pmax(2L, c(n %/% 4L, n %/% 2L, 3L * n %/% 4L, n))))
  }
  sizes <- as.integer(sizes)
  if (length(sizes) < 2L) stop("at least two subsample sizes are required")
  if (any(sizes < 2L) || any(sizes > n)) stop("subsample sizes must be between 2 and n")
  if (reps <= 0L) stop("reps must be positive")
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(subsamples) && length(subsamples) != length(sizes)) {
    stop("subsamples must contain one matrix per subsample size")
  }

  scope_tuple <- qnca_scope_tuple(X, Y, scope)
  effects <- array(NA_real_, dim = c(length(sizes), reps, length(pi_pair)))

  for (k in seq_along(sizes)) {
    m <- sizes[k]
    if (is.null(subsamples)) {
      rows <- t(replicate(reps, sample.int(n, m, replace = FALSE)))
    } else {
      rows <- as.matrix(subsamples[[k]])
      storage.mode(rows) <- "integer"
      if (!all(dim(rows) == c(reps, m))) stop("each subsample matrix must be reps x size")
      if (any(rows < 1L) || any(rows > n)) {
        stop("subsample indices must be one-based and within 1:n")
      }
    }
    for (rr in seq_len(reps)) {
      idx <- rows[rr, ]
      for (p in seq_along(pi_pair)) {
        effects[k, rr, p] <- qnca_d_for_pi(X[idx], Y[idx], pi_pair[p], scope_tuple, n_grid)
      }
    }
  }

  mean_d <- apply(effects, c(1, 3), mean)
  sd_d <- if (reps > 1L) apply(effects, c(1, 3), stats::sd) else matrix(0, nrow = length(sizes), ncol = length(pi_pair))
  g <- 1 / as.numeric(sizes)
  denom <- sum((g - mean(g))^2)
  beta <- rep(NA_real_, length(pi_pair))
  if (denom > 0) {
    for (p in seq_along(pi_pair)) {
      ybar <- mean_d[, p]
      beta[p] <- sum((g - mean(g)) * (ybar - mean(ybar))) / denom
    }
  }
  divergence <- if (all(is.finite(beta))) abs(beta[1]) - abs(beta[2]) else NA_real_
  summary <- data.frame(
    size = rep(sizes, times = length(pi_pair)),
    pi = rep(pi_pair, each = length(sizes)),
    mean_d = as.vector(mean_d),
    sd_d = as.vector(sd_d)
  )
  list(summary = summary, effects = effects, beta_drift = beta,
       divergence = divergence, sizes = sizes, pi_values = pi_pair)
}
