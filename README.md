# QNCA — Quantile Necessary Condition Analysis

Quantile Necessary Condition Analysis (QNCA) is a tolerance-parameterised
extension of Necessary Condition Analysis. It keeps the NCA bottleneck reading,
empty-zone effect size, and permutation logic. At each outcome target, the
frontier is a low type-1 quantile of the condition among the cases that attain
the target, indexed by a tolerance `pi` in (0, 1]; `1 - pi` is the share of
attaining cases allowed below the reported requirement.

At `pi = 1.00`, QNCA reproduces the CE-FDH ceiling of NCA. At lower tolerances,
a requirement no longer rests on the single most efficient attaining case, and
each bottleneck absorbs a stated number of added efficient cases. The result is
a set of tolerance-indexed frontiers reported beside the CE-FDH anchor.

## Repository Layout

The `main` branch contains the shared datasets and validation results.

Language implementations are published on separate branches. The current
release of all three is 0.4.0, and each branch lists its changes in
`CHANGELOG.md`.

| Language | Branch | Install |
|---|---|---|
| R | `qnca@R` | `remotes::install_github("MarekDejaUJ/qnca", ref = "qnca@R")` |
| Python | `qnca@py` | `python -m pip install "git+https://github.com/MarekDejaUJ/qnca.git@qnca%40py"` |
| Julia | `qnca@Julia` | `Pkg.add(url="https://github.com/MarekDejaUJ/qnca.git", rev="qnca@Julia")` |

Quote branch names in shell commands because they contain `@`. The three
implementations agree to floating point on shared data, scopes, tolerances, and
permutation streams.

## Data

The `data/` directory contains the synthetic datasets and result tables used for
strict-limit validation, cross-language agreement, bottleneck comparisons,
outlier-fragility diagnostics, calibration checks, power checks, and the
pi-resolved spuriousness-band illustration. The same files are shipped on the
three package branches.

## Method

Let `X` be the candidate necessary condition and `Y` the desired outcome, both
on a comparable scale. The scope is the rectangle
`[x_min, x_max] x [y_min, y_max]`.

### Raw frontier

For each target outcome level `y`, QNCA looks only at the `k` cases that reached
that level and selects a low order statistic of their condition values:

```text
phi_pi(y) = X_(h; y),    h = max(1, ceil(k * (1 - pi))),    pi in (0, 1].
```

This is the type-1 empirical quantile. The rank `h` is evaluated in exact
arithmetic from the decimal value of `pi`, so `pi = 0.95` with `k = 100` selects
the fifth order statistic. At `pi = 1`, `h = 1` and `phi_1(y)` is the CE-FDH
running minimum. When `k < 1 / (1 - pi)`, the selected rank is still the
minimum, so tolerant frontiers coincide with CE-FDH where few cases attain the
target.

### Monotone envelope

A necessity boundary is non-decreasing. Every case attaining a higher target
also attains every lower one, so a requirement established at a lower target
also holds at a higher one. The fitted frontier is the monotone envelope of the
raw frontier, its running maximum over the outcome grid:

```text
phi_tilde_pi(y_j) = max_{l <= j} phi_pi(y_l).
```

The envelope carries each requirement upward and never lowers it, and every
fitted value is an observed order statistic. Least-squares isotonic projection
(pool-adjacent-violators) is available as an alternative
(`monotone = "isotonic"`; Julia `:isotonic`). It replaces each decreasing block
by its mean, so a low requirement at a sparse high target, where the tolerance
cannot resolve a rank above the minimum, is averaged into the resolved targets
below it. The two fits coincide at `pi = 1`.

Under a fixed scope, outcome levels above the highest observed `Y` have no
conditioning set. Dul's CE-FDH counts that top band as fully empty, so the
implementations carry the trailing frontier to `x_max`. This fixed-scope
convention is required for the `pi = 1` validation against Dul's `NCA` package.

### Resistance

Add `c` cases at arbitrary positions in the scope. At every target where
`h > c`, the raw requirement cannot fall below the `(h - c)`th smallest
condition value of the original attaining cases, nor rise above the `(h + c)`th
when that many cases attain the target. The envelope carries the lower bound to
every higher target. At `pi = 0.95`, one added case
is absorbed wherever more than 20 cases attain the target, and three wherever
more than 60 do. At `pi = 1`, `h = 1` everywhere, and one added case can lower
every bottleneck it attains to its own condition value. The fitted frontier and
`d_pi` cannot decrease as `pi` falls, so the tolerance moves the result in one
known direction.

### Effect size, permutation test, and bottlenecks

The effect size keeps Dul's empty-zone-over-scope form:

```text
d_pi = integral_y (phi_tilde_pi(y) - x_min) dy
       / ((x_max - x_min) * (y_max - y_min)).
```

The implementations compute this integral by the trapezoid rule over the
outcome grid and clamp the result to `[0, 1]`. The formula is inherited from
NCA, but the estimator is `pi`-relative: lowering `pi` raises the frontier and
enlarges the empty zone. Dul's conventional magnitude landmarks should therefore
be used only for the `pi = 1` CE-FDH comparison unless new calibration is done
for lower tolerances.

The permutation test keeps Dul's finite-sample tail probability:

```text
p_pi = (1 + #{ b : d_pi^(b) >= d_pi^obs }) / (B + 1).
```

Every permutation recomputes the whole estimator: the quantile frontier, the
monotone envelope, the area, and `d_pi`. The p-value supports a necessity
reading when paired with a meaningful effect size; it does not prove the
substantive theory. With `B = 1999`, the smallest reportable p-value is exactly
`0.0005`, which means that no permutation matched or exceeded the observed
statistic. The Julia implementation evaluates permutations in parallel and
draws them from one seeded stream, so the same seed gives the same p-value at
any thread count.

The practitioner-facing bottleneck table is the frontier read at decision
levels:

```text
T_x(t_y; pi) = phi_tilde_pi(t_y).
```

It is the minimum condition level required for a realistic chance of reaching
`Y >= t_y` at the selected tolerance. Necessity gives a floor, not sufficiency:
being above the threshold does not guarantee the outcome.

## Robustness Tools

- `qnca_resolution` reports, target by target, the conditioning-set size, the
  selected rank, the raw and fitted requirement, whether the fitted value is
  inherited from a lower target, the realised share of attaining cases below
  it, and the resistance `rank - 1`.
- `qnca_outliers` deletes one case at a time, or `k` cases jointly, and reports
  the change in `d_pi` at every tolerance, in the form of the NCA outlier
  screen.
- `qnca_sensitivity` adds cases at declared positions and returns the
  empirical breakdown curve of `d_pi`.

The proposed `spuriousness_band` compares the observed `d_pi` with a
Gaussian-copula benchmark that preserves the empirical marginals and rank
dependence but imposes no attainability floor. It is a diagnostic with
unestablished operating characteristics, reported beside the permutation test.

## Current Results

All numbers below come from the tracked result files under `data/`.

### Gate Summary

| Check | Result |
|---|---:|
| Cross-language max `d_pi` difference | `4.390e-13` |
| Cross-language p-values equal | `true` |
| Native CE-FDH vs Dul CE-FDH max `d` difference | `1.166e-15` |

Interpretation: Python, R, and Julia agree to floating point on the same data,
scope, `pi` grid, and permutation stream (`data/cross_language_agreement.csv`).
The native CE-FDH implementation matches the CE-FDH effect size of Dul's `NCA`
package (version 5.0.2) under the fixed `[0,100]^2` scope
(`data/dul_nca_results.csv`).

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

These values use `B = 1999`, scope `[0,100]^2`, and `n_grid = 50`.

| Dataset | n | pi | `d_pi` | `p_pi` |
|---|---:|---:|---:|---:|
| n100 | 100 | 1.00 | 0.496231 | 0.000500 |
| n100 | 100 | 0.95 | 0.521828 | 0.000500 |
| n100 | 100 | 0.90 | 0.546651 | 0.000500 |
| n200 | 200 | 1.00 | 0.494148 | 0.000500 |
| n200 | 200 | 0.95 | 0.540643 | 0.000500 |
| n200 | 200 | 0.90 | 0.591639 | 0.000500 |
| n500 | 500 | 1.00 | 0.451974 | 0.000500 |
| n500 | 500 | 0.95 | 0.509731 | 0.000500 |
| n500 | 500 | 0.90 | 0.558381 | 0.000500 |
| n1000 | 1000 | 1.00 | 0.451163 | 0.000500 |
| n1000 | 1000 | 0.95 | 0.515166 | 0.000500 |
| n1000 | 1000 | 0.90 | 0.558699 | 0.000500 |

Interpretation: the empty-zone signal is strong on the shipped reverse-L
datasets. Lower `pi` increases `d_pi` because the frontier rises once efficient
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

For the n1000 dataset, the thresholds below show how the tolerance changes the
required condition floor (`data/bottleneck_comparison_n1000.csv`).

| Outcome level | Dul CE-FDH | QNCA pi=1.00 | QNCA pi=0.95 | QNCA pi=0.90 |
|---:|---:|---:|---:|---:|
| 0 | 0.082 | 0.082 | 4.888 | 9.423 |
| 20 | 7.023 | 7.023 | 17.293 | 25.519 |
| 50 | 43.935 | 43.935 | 54.840 | 58.618 |
| 70 | 67.140 | 67.140 | 71.076 | 77.669 |
| 80 | 88.656 | 88.656 | 88.656 | 88.656 |

Interpretation: at `pi = 1`, QNCA reproduces Dul's CE-FDH bottleneck. At lower
`pi`, thresholds rise at lower and middle outcome levels because unusually
efficient cases no longer set the floor. At the top of the observed outcome
range, few cases attain the target, the selected rank is the minimum, and all
tolerances return the CE-FDH value.

### Outlier-Fragility Diagnostic

The outlier diagnostic adds one efficient case at `(X, Y) = (2, 95)` to an
otherwise clean reverse-L sample of 300 cases (`data/outlier_diagnostic.csv`).

| Case | `d_pi(1.00)` | `d_pi(0.95)` | Gap |
|---|---:|---:|---:|
| clean | 0.457386 | 0.515699 | 0.058312 |
| with efficient outlier | 0.064953 | 0.451162 | 0.386209 |

Interpretation: the strict CE-FDH effect collapses from 0.457 to 0.065 when one
unusually efficient case is added. The `pi = 0.95` effect moves from 0.516 to
0.451, because the frontier absorbs the case wherever more than 20 cases attain
the target and the envelope holds the requirement reached there above that
level. The widening gap is the diagnostic: a wide gap means the strict floor is
held down by unusually efficient boundary cases.

### Calibration and Power

The null calibration run uses 300 independent uniform samples with `n = 150` and
`B = 999` (`data/calibration_null.csv`).

| pi | reps | rejections at 0.05 | false-positive rate |
|---:|---:|---:|---:|
| 1.00 | 300 | 22 | 0.0733 |
| 0.95 | 300 | 12 | 0.0400 |
| 0.90 | 300 | 8 | 0.0267 |

Interpretation: under independence, the tolerant frontiers sit below the 0.05
nominal level in this fixed-seed run, while the strict `pi = 1` frontier is
somewhat liberal. Other null designs give different rates
(`data/calibration_reconciliation.csv`), so these values are screening checks,
not fixed decision thresholds.

The power study varies ceiling strength directly. At strength 0, there is no
necessity ceiling; at strength 1, the attainable ceiling is tightest. Each row
uses 150 samples with `n = 120` and `B = 499` (`data/power_curve.csv`).

| Ceiling strength | Rejections | Power | Mean `d_pi(0.95)` |
|---:|---:|---:|---:|
| 0.0 | 8 | 0.0533 | 0.074720 |
| 0.2 | 67 | 0.4467 | 0.144981 |
| 0.4 | 135 | 0.9000 | 0.241819 |
| 0.6 | 149 | 0.9933 | 0.350038 |
| 0.8 | 150 | 1.0000 | 0.448798 |
| 1.0 | 150 | 1.0000 | 0.534241 |

Interpretation: both the rejection rate and the mean tolerant effect size rise
as the empty corner becomes a stronger structural feature of the data.

## Basic Use

### Python

```python
from qnca import generate_reverse_l, qnca, qnca_outliers, qnca_resolution

X, Y = generate_reverse_l(500, rng=42)
res = qnca(X, Y, pi=0.95, B=1999, scope=(0, 100, 0, 100), seed=1)
print(res.d_pi, res.p_pi)
print(res.bottleneck[:5])
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
head(res$bottleneck)
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

## License

GPL-3.0-only.
