# Changelog

The Julia, R, and Python implementations share one version number. Each
release keeps the three implementations in agreement on shared inputs.

## 0.4.0

- The raw quantile frontier is made non-decreasing by its monotone envelope,
  the running maximum over outcome targets, instead of least-squares isotonic
  projection. Least-squares pooling averaged a low requirement at a sparse
  high target, where the tolerance cannot resolve a rank above the minimum,
  into the resolved targets below it, so one efficient case could lower the
  whole tolerant frontier. The envelope carries each requirement upward and
  leaves every fitted value on an observed order statistic. Results at
  `pi = 1` are unchanged. The previous behaviour remains available with
  `monotone = "isotonic"` (Julia `:isotonic`).
- New `qnca_resolution`: conditioning-set size, selected rank, raw and fitted
  requirement, inherited flag, realised exceptions, and resistance per target.
- New `qnca_outliers`: single and joint deletion influence on `d_pi` across a
  tolerance grid.
- New `qnca_sensitivity`: `d_pi` after adding cases at declared positions.
- New `monotone_envelope`, and a `monotone` argument on the frontier, the fit,
  the permutation test, and the reference band.

## 0.3.2

- The selected order statistic `max(1, ceil(k * (1 - pi)))` is now computed
  from the exact decimal value of the tolerance. In 0.3.1 the product was
  formed in binary floating point, so for some decimal tolerances it landed
  just above a whole number and the next rank was selected: with `pi = 0.95`
  and a conditioning set of `k = 100` cases, 0.3.1 returned the sixth order
  statistic instead of the fifth. Results at `pi = 1` are unchanged, and a
  frontier changes only at outcome targets where `k * (1 - pi)` is a whole
  number.
- Observations, tolerances, outcome grids, and fixed scopes are validated
  before estimation, and invalid input raises an error instead of returning a
  silent result.
- Julia: a seeded permutation test draws all permutations from one stream
  before the threaded evaluation, so the same seed gives the same p-value at
  any thread count. The reference-band replicates are evaluated in parallel
  with the same deterministic draws.
