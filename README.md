# qnca: Quantile Necessary Condition Analysis

Quantile Necessary Condition Analysis (QNCA) extends Necessary Condition Analysis (NCA) with a tolerance. For each outcome target, the QNCA frontier is a low type-1 quantile of the condition among the cases that attain the target. The tolerance `pi` in (0, 1] sets the share `1 - pi` of attaining cases allowed below the reported requirement, and the analyst fixes `pi` before the analysis.

At `pi = 1`, the QNCA frontier equals the CE-FDH ceiling of NCA. At `pi < 1`, each requirement rests on a stated order statistic of the attaining cases, so each requirement absorbs a known number of unusually efficient cases. QNCA keeps the bottleneck table, the empty-zone effect size and the permutation test of NCA, and the frontiers at several tolerances are reported beside the CE-FDH ceiling.

The guide to QNCA is at <https://marekdejauj.github.io/qnca/>.

## Packages

The repository holds three implementations of release 0.4.0, one per folder, and each folder lists its changes in its changelog.

| Language | Package | Folder | Release |
|---|---|---|---|
| R | `qnca` | [`r/`](https://github.com/MarekDejaUJ/qnca/tree/main/r) | `install.packages("qnca")` |
| Python | `qnca` | [`python/`](https://github.com/MarekDejaUJ/qnca/tree/main/python) | `python -m pip install qnca` |
| Julia | `QuantileNecessity` | [`julia/QuantileNecessity.jl/`](https://github.com/MarekDejaUJ/qnca/tree/main/julia/QuantileNecessity.jl) | `Pkg.add("QuantileNecessity")` |

The development versions install from the folders:

```r
remotes::install_github("MarekDejaUJ/qnca", subdir = "r")
```

```bash
python -m pip install "git+https://github.com/MarekDejaUJ/qnca.git#subdirectory=python"
```

```julia
using Pkg
Pkg.add(url="https://github.com/MarekDejaUJ/qnca.git", subdir="julia/QuantileNecessity.jl")
```

The branches `qnca@R`, `qnca@py` and `qnca@Julia` keep the 0.4.0 code as first published. The three implementations agree to floating point on shared data, scopes, tolerances and permutation streams.

## Method

### Raw frontier

Let `X` be the candidate necessary condition and `Y` the outcome, both on a comparable scale, and let the scope be the rectangle `[x_min, x_max] x [y_min, y_max]`. For each outcome target `y`, QNCA takes the `k` cases with `Y >= y` and selects a low order statistic of their condition values:

```text
phi_pi(y) = X_(h; y),    h = max(1, ceil(k * (1 - pi))),    pi in (0, 1].
```

The selection is the type-1 empirical quantile. The rank `h` comes from the exact decimal value of `pi`, so `pi = 0.95` with `k = 100` selects the fifth order statistic. At `pi = 1`, the rank `h` equals 1 and `phi_1(y)` is the CE-FDH running minimum. Where fewer than `1 / (1 - pi)` cases attain a target, the selected rank is the minimum, and the tolerant frontier coincides with the CE-FDH ceiling at the target.

### Monotone envelope

A necessity frontier is non-decreasing: every case that attains a higher target attains every lower target, so a requirement established at a lower target holds at every higher target. The fitted frontier is the monotone envelope of the raw frontier, the running maximum over the outcome grid:

```text
phi_tilde_pi(y_j) = max_{l <= j} phi_pi(y_l).
```

The envelope carries each requirement upward, and every fitted value is an observed order statistic. The least-squares isotonic projection (pool-adjacent-violators) is an option, with `monotone = "isotonic"` in R and Python and `monotone = :isotonic` in Julia. The projection replaces each decreasing block by its mean, which averages a low requirement at a sparse high target into the resolved targets below the block. The envelope and the projection coincide at `pi = 1`.

Under a fixed scope, outcome levels above the highest observed `Y` have an empty conditioning set. The CE-FDH ceiling counts the top band as fully empty, and the implementations carry the last frontier value to `x_max`, which the comparison with the 'NCA' package at `pi = 1` requires.

### Resistance

Suppose `c` cases are added at arbitrary positions in the scope. At every target where `h > c`, the raw requirement stays at or above the `(h - c)`th smallest condition value of the original attaining cases, and at or below the `(h + c)`th value when at least `h + c` cases attain the target. The envelope carries the lower bound to every higher target. At `pi = 0.95`, the frontier absorbs one added case wherever more than 20 cases attain the target, and three added cases wherever more than 60 cases attain the target. At `pi = 1`, the rank `h` equals 1 everywhere, and one added case can lower every bottleneck the case attains to the condition value of the case. The fitted frontier and `d_pi` rise or stay level as `pi` falls, so the tolerance moves the result in one known direction.

### Effect size, permutation test and bottlenecks

The effect size keeps the empty-zone-over-scope form of NCA:

```text
d_pi = integral_y (phi_tilde_pi(y) - x_min) dy
       / ((x_max - x_min) * (y_max - y_min)).
```

The implementations integrate by the trapezoid rule over the outcome grid and clamp the result to `[0, 1]`. The formula comes from NCA, and the estimator is relative to `pi`: a lower `pi` raises the frontier and enlarges the empty zone. The magnitude landmarks of NCA therefore apply to the comparison at `pi = 1`, and lower tolerances call for a calibration of their own.

The permutation test keeps the finite-sample tail probability of NCA:

```text
p_pi = (1 + #{ b : d_pi^(b) >= d_pi^obs }) / (B + 1).
```

Every permutation recomputes the whole estimator: the quantile frontier, the monotone envelope, the area and `d_pi`. A small p-value with a meaningful effect size supports a necessity interpretation. With `B = 1999`, the smallest reportable p-value is exactly 0.0005, the value at which every permuted effect size stays below the observed effect size. The Julia implementation evaluates permutations in parallel from one seeded stream, so one seed gives one p-value at any number of threads.

The bottleneck table reads the frontier at decision levels:

```text
T_x(t_y; pi) = phi_tilde_pi(t_y).
```

Each entry is the condition level required for a realistic chance of reaching `Y >= t_y` at the selected tolerance. A bottleneck states a necessary floor, and a case at or above the floor may still miss the outcome.

## Diagnostics

`qnca_resolution` reports, target by target, the conditioning-set size, the selected rank, the raw and fitted requirement, a flag for a fitted value inherited from a lower target, the realised share of attaining cases below the requirement, and the resistance `rank - 1`.

`qnca_outliers` deletes one case at a time, or `k` cases jointly, and reports the change in `d_pi` at every tolerance, in the form of the outlier screen of NCA.

`qnca_sensitivity` adds cases at declared positions and returns the empirical breakdown curve of `d_pi`.

The proposed `spuriousness_band` compares the observed `d_pi` with a Gaussian-copula benchmark that preserves the empirical marginals and the rank dependence and imposes no attainability floor. The band is a diagnostic with operating characteristics still to be established, reported beside the permutation test.

## Results

Every number in the tables comes from the result files in [`data/`](https://github.com/MarekDejaUJ/qnca/tree/main/data).

### Agreement

| Check | Result |
|---|---:|
| Cross-language maximum `d_pi` difference | `4.390e-13` |
| Cross-language p-values equal | `true` |
| Native CE-FDH and 'NCA' CE-FDH maximum `d` difference | `1.166e-15` |

Python, R and Julia agree to floating point on the same data, scope, `pi` grid and permutation stream (`data/cross_language_agreement.csv`). The native CE-FDH implementation matches the CE-FDH effect size of the 'NCA' package, version 5.0.2, under the fixed scope `[0,100]^2` (`data/dul_nca_results.csv`).

### Strict limit

The table compares `d_pi` at `pi = 1` with the native exact CE-FDH effect size, with `n_grid = 400` (`data/validation_invariant.csv`).

| Dataset | n | QNCA `d_pi(1)` | Native CE-FDH `d` | Absolute difference |
|---|---:|---:|---:|---:|
| n100 | 100 | 0.495553 | 0.495353 | 2.00e-04 |
| n200 | 200 | 0.492750 | 0.492739 | 1.10e-05 |
| n500 | 500 | 0.450484 | 0.450813 | 3.29e-04 |
| n1000 | 1000 | 0.451198 | 0.451110 | 8.78e-05 |

The difference is the discretisation of the QNCA area integral on the outcome grid. The native CE-FDH value is exact on the observed step frontier.

### Effect size and p-value across tolerances

The values use `B = 1999`, the scope `[0,100]^2` and `n_grid = 50` (`data/pi_grid_results.csv`).

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

The empty zone is large on the four reverse-L datasets. A lower `pi` gives a larger `d_pi`, because the frontier rises once the tolerance admits efficient boundary cases below the requirement. Every p-value sits at the floor of `1/(1999 + 1)`, which means that every permuted effect size stayed below the observed effect size.

### Comparison with the 'NCA' package

The CE-FDH rows are the target of the comparison. The CR-FDH rows show the second standard ceiling of the 'NCA' package.

| Dataset | n | Method | QNCA `d_pi(1)` | Native CE-FDH `d` | 'NCA' `d` | 'NCA' `p` | Native and 'NCA' difference |
|---|---:|---|---:|---:|---:|---:|---:|
| n100 | 100 | ce_fdh | 0.495553 | 0.495353 | 0.495353 | 0.0005002501 | 0 |
| n100 | 100 | cr_fdh | 0.495553 | 0.495353 | 0.468644 | 0.0005002501 | 0.0267089700 |
| n200 | 200 | ce_fdh | 0.492750 | 0.492739 | 0.492739 | 0.0005002501 | 1.665e-16 |
| n200 | 200 | cr_fdh | 0.492750 | 0.492739 | 0.467374 | 0.0005002501 | 0.0253647516 |
| n500 | 500 | ce_fdh | 0.450484 | 0.450813 | 0.450813 | 0.0005002501 | 2.776e-16 |
| n500 | 500 | cr_fdh | 0.450484 | 0.450813 | 0.434528 | 0.0005002501 | 0.0162848199 |
| n1000 | 1000 | ce_fdh | 0.451198 | 0.451110 | 0.451110 | 0.0005002501 | 1.166e-15 |
| n1000 | 1000 | cr_fdh | 0.451198 | 0.451110 | 0.441200 | 0.0005002501 | 0.0099101730 |

The native CE-FDH effect size equals the CE-FDH effect size of the 'NCA' package to floating point. QNCA `d_pi(1)` differs slightly from the CE-FDH value because QNCA integrates over a finite grid. The CR-FDH ceiling is a straight line fitted to the corner points of the CE-FDH ceiling and gives a smaller empty zone on every dataset.

### Bottlenecks

The thresholds for the n1000 dataset show how the tolerance changes the required condition level (`data/bottleneck_comparison_n1000.csv`).

| Outcome level | 'NCA' CE-FDH | QNCA pi=1.00 | QNCA pi=0.95 | QNCA pi=0.90 |
|---:|---:|---:|---:|---:|
| 0 | 0.082 | 0.082 | 4.888 | 9.423 |
| 20 | 7.023 | 7.023 | 17.293 | 25.519 |
| 50 | 43.935 | 43.935 | 54.840 | 58.618 |
| 70 | 67.140 | 67.140 | 71.076 | 77.669 |
| 80 | 88.656 | 88.656 | 88.656 | 88.656 |

At `pi = 1`, QNCA reproduces the CE-FDH bottlenecks of the 'NCA' package. At lower `pi`, the thresholds rise at low and middle outcome levels, where unusually efficient cases set the CE-FDH floor. At the top of the observed outcome range, few cases attain the target, the selected rank is the minimum, and all tolerances return the CE-FDH value.

### Outlier-fragility

The diagnostic adds one efficient case at `(X, Y) = (2, 95)` to a reverse-L sample of 300 cases (`data/outlier_diagnostic.csv`).

| Case | `d_pi(1.00)` | `d_pi(0.95)` | Gap |
|---|---:|---:|---:|
| clean | 0.457386 | 0.515699 | 0.058312 |
| with efficient outlier | 0.064953 | 0.451162 | 0.386209 |

The strict CE-FDH effect falls from 0.457 to 0.065 when the single efficient case joins the sample. The effect at `pi = 0.95` moves from 0.516 to 0.451, because the frontier absorbs the case wherever more than 20 cases attain the target and the envelope carries the requirements of the resolved targets to every higher level. A wide gap between the two effects marks a strict floor held down by unusually efficient boundary cases.

### Calibration and power

The null calibration uses 300 independent uniform samples with `n = 150` and `B = 999` (`data/calibration_null.csv`).

| pi | reps | rejections at 0.05 | false-positive rate |
|---:|---:|---:|---:|
| 1.00 | 300 | 22 | 0.0733 |
| 0.95 | 300 | 12 | 0.0400 |
| 0.90 | 300 | 8 | 0.0267 |

Under independence, the tolerant frontiers reject below the nominal level of 0.05 in the fixed-seed run, and the strict frontier at `pi = 1` rejects above the nominal level. Other null designs give other rates (`data/calibration_reconciliation.csv`), so the rates serve as screening checks and leave the decision thresholds open.

The power study varies the strength of the ceiling, from 0 for a scatter with no necessity floor to 1 for the tightest attainable floor. Each row uses 150 samples with `n = 120` and `B = 499` (`data/power_curve.csv`).

| Ceiling strength | Rejections | Power | Mean `d_pi(0.95)` |
|---:|---:|---:|---:|
| 0.0 | 8 | 0.0533 | 0.074720 |
| 0.2 | 67 | 0.4467 | 0.144981 |
| 0.4 | 135 | 0.9000 | 0.241819 |
| 0.6 | 149 | 0.9933 | 0.350038 |
| 0.8 | 150 | 1.0000 | 0.448798 |
| 1.0 | 150 | 1.0000 | 0.534241 |

The rejection rate and the mean tolerant effect size both rise as the empty corner becomes a stronger feature of the data.

## Use

### R

```r
library(qnca)

set.seed(1)
d <- generate_reverse_L(500)
res <- qnca(d$X, d$Y, pi = 0.95, B = 1999, scope = c(0, 100, 0, 100))
res$d_pi
res$p_pi
head(res$bottleneck)
table <- qnca_resolution(d$X, d$Y, pi = 0.95, scope = c(0, 100, 0, 100))
screen <- qnca_outliers(d$X, d$Y, pi_grid = c(1, 0.95, 0.90),
                        scope = c(0, 100, 0, 100))
```

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

### Julia

```julia
using QuantileNecessity
using Random

X, Y = generate_reverse_L(500; rng=MersenneTwister(42))
res = qnca(X, Y; pi=0.95, B=1999, scope=(0.0, 100.0, 0.0, 100.0), seed=1)
res.d_pi, res.p_pi
table = qnca_resolution(X, Y; pi=0.95, scope=(0.0, 100.0, 0.0, 100.0))
screen = qnca_outliers(X, Y; pi_grid=[1.0, 0.95, 0.90],
                       scope=(0.0, 100.0, 0.0, 100.0))
```

## Data

The `data/` folder holds the synthetic datasets and the result tables of the strict-limit validation, the cross-language agreement, the bottleneck comparison, the outlier-fragility diagnostic, the calibration and power checks, and the illustration of the reference band.

## Citation

Deja, M. (2026). Quantile Necessary Condition Analysis: A Tolerance-Parameterized Generalization of Necessary Condition Analysis. *Quality & Quantity* (accepted for publication).

## Licence

GPL-3.0-only.
