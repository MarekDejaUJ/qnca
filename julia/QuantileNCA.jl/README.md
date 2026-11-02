# QNCA — Quantile Necessary Condition Analysis

Quantile Necessary Condition Analysis (QNCA) is a tolerance-parameterised
extension of Necessary Condition Analysis. It keeps the NCA bottleneck reading,
empty-zone effect size, and permutation logic. At each outcome target, the
frontier is a low type-1 quantile of the condition among the cases that attain
the target, indexed by a tolerance `pi` in (0, 1]; `1 - pi` is the share of
attaining cases allowed below the reported requirement.

At `pi = 1.00`, QNCA reproduces the CE-FDH ceiling of NCA. At lower tolerances,
a requirement no longer rests on the single most efficient attaining case. The
result is a set of tolerance-indexed frontiers reported beside the CE-FDH
anchor.

## Monotone step

The raw quantile frontier is made non-decreasing by its monotone envelope, the
running maximum over outcome targets: a requirement established at a lower
target is carried to every higher target. This is the default from 0.4.0. The
least-squares isotonic projection of 0.3 remains available with
`monotone = "isotonic"`. The two coincide at `pi = 1`. Least-squares pooling
averages an unresolved dip at a sparse high target into the resolved targets
below it; the envelope does not.

## Robustness tools

- `qnca_resolution` reports, target by target, the conditioning-set size, the
  selected rank, the raw and fitted requirement, whether the fitted value is
  inherited from a lower target, the realised share of attaining cases below
  it, and the resistance `rank - 1`.
- `qnca_outliers` deletes one case at a time, or `k` cases jointly, and reports
  the change in `d_pi` at every tolerance, in the form of the NCA outlier
  screen.
- `qnca_sensitivity` adds cases at declared positions and returns the
  empirical breakdown curve of `d_pi`.

## Data

The `data/` directory contains the synthetic datasets and result tables used by
the package examples and validation checks.

## Basic Use

### Python

```python
from qnca import generate_reverse_l, qnca, qnca_outliers, qnca_resolution

X, Y = generate_reverse_l(500, rng=42)
res = qnca(X, Y, pi=0.95, B=1999, scope=(0, 100, 0, 100), seed=1)
print(res.d_pi, res.p_pi)
table = qnca_resolution(X, Y, pi=0.95, scope=(0, 100, 0, 100))
screen = qnca_outliers(X, Y, pi_grid=(1.0, 0.95, 0.90), scope=(0, 100, 0, 100))
```

### R

```r
library(qnca)

d <- generate_reverse_L(500)
res <- qnca(d$X, d$Y, pi = 0.95, B = 1999, scope = c(0, 100, 0, 100))
res$d_pi
res$p_pi
table <- qnca_resolution(d$X, d$Y, pi = 0.95, scope = c(0, 100, 0, 100))
screen <- qnca_outliers(d$X, d$Y, pi_grid = c(1, 0.95, 0.90),
                        scope = c(0, 100, 0, 100))
```

### Julia

```julia
using QNCA
using Random

rng = MersenneTwister(42)
X, Y = generate_reverse_L(500; rng=rng)
res = qnca(X, Y; pi=0.95, B=1999, scope=(0.0, 100.0, 0.0, 100.0), seed=1)
res.d_pi, res.p_pi
table = qnca_resolution(X, Y; pi=0.95, scope=(0.0, 100.0, 0.0, 100.0))
screen = qnca_outliers(X, Y; pi_grid=[1.0, 0.95, 0.90],
                       scope=(0.0, 100.0, 0.0, 100.0))
```

## Versions

Changes between releases are listed in `CHANGELOG.md`.

## License

GPL-3.0-only.
