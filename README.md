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

## Repository Layout

The `main` branch contains the shared datasets and validation results.

Language implementations are published on separate branches:

| Language | Branch | Install |
|---|---|---|
| R | `qnca@R` | `remotes::install_github("MarekDejaUJ/qnca", ref = "qnca@R")` |
| Python | `qnca@py` | `python -m pip install "git+https://github.com/MarekDejaUJ/qnca.git@qnca%40py"` |
| Julia | `qnca@Julia` | `Pkg.add(url="https://github.com/MarekDejaUJ/qnca.git", rev="qnca@Julia")` |

Quote branch names in shell commands because they contain `@`.

## Data

The `data/` directory contains the synthetic datasets and result tables used for
strict-limit validation, cross-language agreement, bottleneck
comparisons, outlier-fragility diagnostics, calibration checks, power checks,
and the pi-resolved spuriousness-band illustration.

## Method

Let `X` be the candidate necessary condition and `Y` the desired outcome, both
on a comparable scale. The scope is the rectangle
`[x_min, x_max] x [y_min, y_max]`.

For each target outcome level `y`, QNCA looks only at cases that reached that
level and computes a low quantile of their condition values:

```text
phi_pi(y) = Q_(1 - pi)({ X_i : Y_i >= y }),    pi in (0, 1].
```

At `pi = 1`, `Q_0` is the minimum, so `phi_1(y)` is the CE-FDH running-minimum
frontier. At `pi = 0.95`, the most efficient 5 percent of each conditioning set
no longer fix the frontier. The raw frontier is projected onto the
non-decreasing cone with pool-adjacent-violators isotonic regression:

```text
phi_tilde_pi = isoreg(y_grid, phi_pi(y_grid)), non-decreasing.
```

Under a fixed scope, outcome levels above the highest observed `Y` have no
conditioning set. Dul's CE-FDH counts that top band as fully empty, so the
shipped QNCA implementations carry the trailing frontier to `x_max`. This
fixed-scope convention is required for the `pi = 1` validation against Dul's
`NCA` package.

The effect size keeps Dul's empty-zone-over-scope form:

```text
d_pi = integral_y (phi_tilde_pi(y) - x_min) dy
       / ((x_max - x_min) * (y_max - y_min)).
```

The package computes this integral by the trapezoid rule over the outcome grid
and clamps the result to `[0, 1]`. The formula is inherited from NCA, but the
estimator is `pi`-relative: lowering `pi` raises the frontier and mechanically
enlarges the empty zone. Dul's conventional magnitude landmarks should therefore
be used only for the `pi = 1` CE-FDH comparison unless new calibration is done
for lower tolerances.

The permutation test also keeps Dul's finite-sample tail probability:

```text
p_pi = (1 + #{ b : d_pi^(b) >= d_pi^obs }) / (B + 1).
```

Every permutation recomputes the whole estimator: the quantile frontier,
isotonic projection, area, and `d_pi`. The p-value supports a necessity reading
when paired with a meaningful effect size; it does not prove the substantive
theory. With `B = 1999`, the smallest reportable p-value is exactly `0.0005`,
which means that no permutation matched or exceeded the observed statistic.

The practitioner-facing bottleneck table is the frontier read at decision
levels:

```text
T_x(t_y; pi) = phi_tilde_pi(t_y).
```

It is the minimum condition level required for a realistic chance of reaching
`Y >= t_y` at the selected tolerance. Necessity gives a floor, not sufficiency:
being above the threshold does not guarantee the outcome.

## Current Results

All numbers below come from the tracked result files under `data/`. The full
audit file is `data/qnca_comparison.log`.

### Gate Summary

| Check | Result |
|---|---:|
| Cross-language max `d_pi` difference | `4.540e-13` |
| Cross-language p-values equal | `true` |
| Cross-language bottleneck max difference | `1.000e-12` |
| Native CE-FDH vs Dul CE-FDH max `d` difference | `1.166e-15` |

Interpretation: Python, R, and Julia agree to floating point on the same data,
scope, `pi` grid, and permutation stream. The native CE-FDH implementation also
matches Dul's `NCA` CE-FDH effect size under the fixed `[0,100]^2` scope.

### Validation Invariant

`d_pi(pi = 1)` is compared to the native exact CE-FDH effect size using
`n_grid = 400`.

| Dataset | n | QNCA `d_pi(1)` | Native CE-FDH `d` | Absolute difference |
|---|---:|---:|---:|---:|
| n100 | 100 | 0.495553 | 0.495353 | 2.00e-04 |
| n200 | 200 | 0.492750 | 0.492739 | 1.10e-05 |
| n500 | 500 | 0.450484 | 0.450813 | 3.29e-04 |
| n1000 | 1000 | 0.451198 | 0.451110 | 8.78e-05 |

The remaining difference is grid discretization in the QNCA area integral. The
native CE-FDH value is exact on the observed step frontier.

### QNCA d and p Across pi

These values use `B = 1999`, scope `[0,100]^2`, and `n_grid = 50`. The Python,
R, and Julia values are identical to floating point; one column is shown for
readability.

| Dataset | n | pi | `d_pi` | `p_pi` |
|---|---:|---:|---:|---:|
| n100 | 100 | 1.00 | 0.496231 | 0.000500 |
| n100 | 100 | 0.95 | 0.522017 | 0.000500 |
| n100 | 100 | 0.90 | 0.541757 | 0.000500 |
| n200 | 200 | 1.00 | 0.494148 | 0.000500 |
| n200 | 200 | 0.95 | 0.542434 | 0.000500 |
| n200 | 200 | 0.90 | 0.584832 | 0.000500 |
| n500 | 500 | 1.00 | 0.451974 | 0.000500 |
| n500 | 500 | 0.95 | 0.510040 | 0.000500 |
| n500 | 500 | 0.90 | 0.556598 | 0.000500 |
| n1000 | 1000 | 1.00 | 0.451163 | 0.000500 |
| n1000 | 1000 | 0.95 | 0.514367 | 0.000500 |
| n1000 | 1000 | 0.90 | 0.556499 | 0.000500 |

Interpretation: the empty-zone signal is strong on the shipped reverse-L
datasets. Lower `pi` increases `d_pi` because the frontier rises after efficient
boundary cases are tolerated. The p-values sit at the `1/(1999+1)` floor, so the
correct reading is "no permutation matched or exceeded the observed effect",
not "p = 0".

### Dul NCA Comparison

The CE-FDH rows are the validation target. CR-FDH is included because Dul's
package reports it as a standard alternative ceiling.

| Dataset | n | Method | QNCA `d_pi(1)` | Native CE-FDH `d` | Dul `d` | Dul `p` | Native vs Dul diff |
|---|---:|---|---:|---:|---:|---:|---:|
| n100 | 100 | ce_fdh | 0.495553 | 0.495353 | 0.495353 | 0.0005002501 | 0 |
| n100 | 100 | cr_fdh | 0.495553 | 0.495353 | 0.468644 | 0.0005002501 | 0.0267089700 |
| n200 | 200 | ce_fdh | 0.492750 | 0.492739 | 0.492739 | 0.0005002501 | 1.665e-16 |
| n200 | 200 | cr_fdh | 0.492750 | 0.492739 | 0.467374 | 0.0005002501 | 0.0253647516 |
| n500 | 500 | ce_fdh | 0.450484 | 0.450813 | 0.450813 | 0.0005002501 | 2.776e-16 |
| n500 | 500 | cr_fdh | 0.450484 | 0.450813 | 0.434528 | 0.0005002501 | 0.0162848199 |
| n1000 | 1000 | ce_fdh | 0.451198 | 0.451110 | 0.451110 | 0.0005002501 | 1.166e-15 |
| n1000 | 1000 | cr_fdh | 0.451198 | 0.451110 | 0.441200 | 0.0005002501 | 0.0099101730 |

Interpretation: the CE-FDH equality holds against Dul's package. QNCA
`d_pi(1)` differs slightly from Dul CE-FDH because QNCA integrates over a finite
grid; the native CE-FDH function integrates the same CE-FDH geometry exactly and
matches Dul to floating point. CR-FDH is not the reduction target, so its
different `d` values are expected.

### Bottlenecks

The full bottleneck comparison is in `data/qnca_comparison.log`. For the n1000
dataset, the QNCA thresholds below show how the tolerance dial changes the
required condition floor.

| Outcome level | QNCA pi=1.00 | QNCA pi=0.95 | QNCA pi=0.90 | Dul CE-FDH |
|---:|---:|---:|---:|---:|
| 0 | 0.082 | 4.907 | 9.423 | 0.082 |
| 10 | 0.271 | 9.675 | 14.923 | 0.271 |
| 20 | 7.023 | 17.293 | 25.519 | 7.023 |
| 30 | 15.703 | 28.059 | 39.408 | 15.703 |
| 40 | 28.059 | 39.940 | 47.176 | 28.059 |
| 50 | 43.935 | 52.673 | 58.618 | 43.935 |
| 60 | 56.060 | 58.618 | 60.292 | 56.060 |
| 70 | 67.140 | 69.108 | 74.373 | 67.140 |
| 80 | 88.656 | 88.656 | 88.656 | 88.656 |
| 90 | 100.000 | 100.000 | 100.000 | NA |
| 100 | 100.000 | 100.000 | 100.000 | NA |

Interpretation: at `pi = 1`, QNCA reproduces Dul's CE-FDH bottleneck through the
observed outcome range. At lower `pi`, thresholds rise at lower and middle
outcome levels because unusually efficient cases no longer set the floor. The
`100.000` values at the highest fixed-scope levels are the required fixed-scope
top-band convention: above the highest observed outcome, the method counts the
whole band as outside the observed attainable set.

### Robustness and Calibration

The outlier diagnostic plants one efficient case with high `Y` and low `X` into
an otherwise clean reverse-L sample.

| Case | `d_pi(1.00)` | `d_pi(0.95)` | Gap |
|---|---:|---:|---:|
| clean | 0.457386 | 0.517198 | 0.059812 |
| with efficient outlier | 0.064953 | 0.250286 | 0.185333 |

Interpretation: the strict CE-FDH-like verdict collapses when one unusually
efficient case is added, while the tolerant frontier degrades less severely. The
widening gap is the diagnostic: a wide gap means the strict floor is being held
down by unusually efficient boundary cases.

The null calibration run uses 300 independent null samples and `B = 999`.

| pi | reps | rejections at 0.05 | false-positive rate |
|---:|---:|---:|---:|
| 1.00 | 300 | 22 | 0.0733 |
| 0.95 | 300 | 8 | 0.0267 |
| 0.90 | 300 | 10 | 0.0333 |

Interpretation: under independence, the tolerant frontiers sit below the 0.05
nominal level in this fixed-seed run, while the strict `pi = 1` frontier is
somewhat liberal. That is consistent with the method's limitation that the strict
CE-FDH limit inherits boundary fragility.

The power study varies ceiling strength directly. At strength 0, there is no
necessity ceiling; at strength 1, the attainable ceiling is tightest. Each row
uses 150 samples and `B = 499`.

| Ceiling strength | Rejections | Power | Mean `d_pi(0.95)` |
|---:|---:|---:|---:|
| 0.0 | 7 | 0.0467 | 0.063818 |
| 0.2 | 69 | 0.4600 | 0.135803 |
| 0.4 | 137 | 0.9133 | 0.235061 |
| 0.6 | 149 | 0.9933 | 0.345156 |
| 0.8 | 150 | 1.0000 | 0.445371 |
| 1.0 | 150 | 1.0000 | 0.532609 |

Interpretation: both the rejection rate and mean tolerant effect size rise as
the empty corner becomes a stronger structural feature of the data.

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

## License

GPL-3.0-only.
