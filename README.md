# QNCA — Quantile Necessary Condition Analysis

Quantile Necessary Condition Analysis (QNCA) is a tolerance-parameterised
extension of Necessary Condition Analysis. It keeps the NCA bottleneck reading,
empty-zone effect size, and permutation logic, while replacing the deterministic
CE-FDH frontier with a family of conditional-quantile frontiers indexed by a
tolerance parameter `pi`.

At `pi = 1.00`, QNCA reproduces the CE-FDH member. At lower tolerances, the
frontier is less dependent on unusually efficient boundary observations. The
reported result is therefore a tolerance-indexed family of bottlenecks rather
than a single deterministic floor.

## Data

The `data/` directory contains the synthetic datasets and result tables used by
the package examples and validation checks.

## Basic Use

### Python

```python
from qnca import generate_reverse_l, qnca

X, Y = generate_reverse_l(500, rng=42)
res = qnca(X, Y, pi=0.95, B=1999, scope=(0, 100, 0, 100), seed=1)
print(res.d_pi, res.p_pi)
print(res.bottleneck[:5])
```

### R

```r
library(qnca)

d <- generate_reverse_L(500)
res <- qnca(d$X, d$Y, pi = 0.95, B = 1999, scope = c(0, 100, 0, 100))
res$d_pi
res$p_pi
head(res$bottleneck)
```

### Julia

```julia
using QNCA
using Random

rng = MersenneTwister(42)
X, Y = generate_reverse_L(500; rng=rng)
res = qnca(X, Y; pi=0.95, B=1999, scope=(0.0, 100.0, 0.0, 100.0), seed=1)
res.d_pi, res.p_pi
```

## Versions

Changes between releases are listed in `CHANGELOG.md`.

## License

GPL-3.0-only.
