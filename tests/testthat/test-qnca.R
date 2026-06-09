test_that("validation invariant: d_pi(pi = 1) equals native CE-FDH d", {
  set.seed(42)
  for (n in c(100, 200, 500, 1000)) {
    d <- generate_reverse_L(n)
    for (scope in list(c(0, 100, 0, 100), NULL)) {
      q1 <- qnca(d$X, d$Y, pi = 1, n_grid = 400L, B = NA, scope = scope)
      ref <- nca_ce_fdh_d(d$X, d$Y, scope = scope)
      expect_lt(abs(q1$d_pi - ref), 2e-3)
    }
  }
})

test_that("frontier is non-decreasing", {
  set.seed(1)
  d <- generate_reverse_L(300)
  y_grid <- seq(0, 100, length.out = 50)
  phi <- qnca_frontier(d$X, d$Y, pi = 0.95, y_grid = y_grid)
  fin <- phi[is.finite(phi)]
  expect_true(all(diff(fin) >= -1e-9))
})

test_that("effect size lies in [0, 1] and grows as pi falls", {
  set.seed(2)
  d <- generate_reverse_L(500)
  d1  <- qnca(d$X, d$Y, pi = 1.00, n_grid = 200L, B = NA, scope = c(0,100,0,100))$d_pi
  d95 <- qnca(d$X, d$Y, pi = 0.95, n_grid = 200L, B = NA, scope = c(0,100,0,100))$d_pi
  d90 <- qnca(d$X, d$Y, pi = 0.90, n_grid = 200L, B = NA, scope = c(0,100,0,100))$d_pi
  expect_true(d1 >= 0 && d90 <= 1)
  expect_lte(d1, d95 + 1e-9)
  expect_lte(d95, d90 + 1e-9)
})

test_that("permutation p-value respects the 1/(B+1) floor", {
  set.seed(7)
  d <- generate_reverse_L(200)
  B <- 199L
  res <- qnca(d$X, d$Y, pi = 0.95, B = B, scope = c(0,100,0,100))
  expect_gte(res$p_pi, 1 / (B + 1) - 1e-12)
  expect_lte(res$p_pi, 1)
})

test_that("supplied permutations control the p-value", {
  set.seed(11)
  d <- generate_reverse_L(50)
  permutations <- matrix(rep(seq_len(nrow(d)), 5L), nrow = 5L, byrow = TRUE)
  res <- qnca(d$X, d$Y, pi = 0.95, B = NA, scope = c(0,100,0,100),
              permutations = permutations)
  expect_equal(res$p_pi, 1)
})

test_that("no necessity signal gives a small effect size", {
  set.seed(0)
  X <- runif(500, 0, 100); Y <- runif(500, 0, 100)
  res <- qnca(X, Y, pi = 1, n_grid = 200L, B = NA, scope = c(0,100,0,100))
  expect_lt(res$d_pi, 0.10)
})

test_that("spuriousness band has expected shape and bounds", {
  set.seed(21)
  d <- generate_reverse_L(80)
  set.seed(3)
  draws <- array(rnorm(12L * nrow(d) * 2L), dim = c(12L, nrow(d), 2L))
  res <- spuriousness_band(d$X, d$Y, pi_grid = c(1, 0.95), n_grid = 40L,
                           M = 12L, scope = c(0,100,0,100),
                           normal_draws = draws)
  expect_equal(dim(res$band), c(2L, 7L))
  expect_equal(dim(res$d_null), c(12L, 2L))
  expect_true(all(is.finite(as.matrix(res$band))))
  expect_true(all(res$band$lower <= res$band$median + 1e-12))
  expect_true(all(res$band$median <= res$band$upper + 1e-12))
  expect_true(all(res$band$p_null >= 1 / 13 - 1e-12))
  expect_true(all(res$band$p_null <= 1))
})

test_that("spuriousness band accepts uniform draws", {
  set.seed(31)
  d <- generate_reverse_L(60)
  set.seed(5)
  uniforms <- array(runif(9L * nrow(d) * 2L), dim = c(9L, nrow(d), 2L))
  res <- spuriousness_band(d$X, d$Y, pi_grid = c(1, 0.95), n_grid = 35L,
                           M = 9L, scope = c(0,100,0,100),
                           uniform_draws = uniforms)
  expect_equal(dim(res$d_null), c(9L, 2L))
  expect_true(all(is.finite(as.matrix(res$band))))
})

test_that("spuriousness band rejects two supplied streams", {
  set.seed(32)
  d <- generate_reverse_L(30)
  uniforms <- array(runif(4L * nrow(d) * 2L), dim = c(4L, nrow(d), 2L))
  normals <- array(rnorm(4L * nrow(d) * 2L), dim = c(4L, nrow(d), 2L))
  expect_error(
    spuriousness_band(d$X, d$Y, M = 4L, scope = c(0,100,0,100),
                      uniform_draws = uniforms, normal_draws = normals),
    "at most one"
  )
})

test_that("consistency probe returns finite summaries", {
  set.seed(22)
  d <- generate_reverse_L(90)
  res <- consistency_probe(d$X, d$Y, pi_pair = c(1, 0.90),
                           sizes = c(30L, 60L, 90L), reps = 8L,
                           n_grid = 40L, scope = c(0,100,0,100), seed = 4)
  expect_equal(dim(res$effects), c(3L, 8L, 2L))
  expect_equal(dim(res$summary), c(6L, 4L))
  expect_true(all(is.finite(res$summary$mean_d)))
  expect_true(all(is.finite(res$summary$sd_d)))
  expect_true(all(is.finite(res$beta_drift)))
  expect_true(is.finite(res$divergence))
})
