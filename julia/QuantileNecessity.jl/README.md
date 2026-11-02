# QuantileNecessity.jl

Quantile Necessary Condition Analysis (QNCA) for Julia.

Quantile Necessary Condition Analysis extends Necessary Condition Analysis (NCA) with a tolerance. For each outcome target, the QNCA frontier is a low type-1 quantile of the condition among the cases that attain the target. The tolerance `pi` in (0, 1] sets the share `1 - pi` of attaining cases allowed below the reported requirement, and the analyst fixes `pi` before the analysis.

At `pi = 1`, the QNCA frontier equals the CE-FDH ceiling of NCA. At `pi < 1`, each requirement rests on a stated order statistic of the attaining cases, so each requirement absorbs a known number of unusually efficient cases. QNCA keeps the bottleneck table, the empty-zone effect size and the permutation test of NCA.

## Installation

The release installs from the General registry:

```julia
using Pkg
Pkg.add("QuantileNecessity")
```

The development version installs from the `julia/QuantileNecessity.jl` folder of the repository:

```julia
using Pkg
Pkg.add(url="https://github.com/MarekDejaUJ/qnca.git", subdir="julia/QuantileNecessity.jl")
```

The package needs Julia 1.6 or later and depends on standard libraries only.

## Example

`generate_reverse_L()` draws a reverse-L scatter, the pattern of a necessary condition. `qnca()` estimates the frontier at `pi = 0.95`, the effect size `d_pi` and the permutation p-value `p_pi` from 1999 permutations.

```julia
using QuantileNecessity, Random

X, Y = generate_reverse_L(500; rng=MersenneTwister(42))
res = qnca(X, Y; pi=0.95, B=1999, scope=(0.0, 100.0, 0.0, 100.0), seed=1)
res.d_pi, res.p_pi
# (0.5097469505372707, 0.0005)
```

A p-value of 0.0005 is the smallest value that 1999 permutations can report: every permuted effect size stayed below the observed effect size. The fields `res.y_grid` and `res.x_required` hold the bottleneck table: each outcome level with the condition level required for the outcome level at the chosen tolerance.

The permutation test runs on all Julia threads (start Julia with `julia -t auto`) and draws every permutation from one seeded stream, so one seed gives one p-value at any number of threads.

## Functions

| Function | Result |
|---|---|
| `qnca` | effect size, permutation p-value and bottleneck table at one tolerance |
| `qnca_frontier` | fitted frontier on a grid of outcome levels |
| `qnca_resolution` | per outcome target: conditioning-set size, selected rank, raw and fitted requirement, inherited flag, share of attaining cases below the requirement, resistance |
| `qnca_outliers` | change of `d_pi` across tolerances after the deletion of one case or of `k` cases |
| `qnca_sensitivity` | `d_pi` after the addition of cases at declared positions |
| `spuriousness_band` | proposed dependence-matched reference band for `d_pi` |
| `consistency_probe` | `d_pi` across subsample sizes at two tolerances |
| `nca_ce_fdh_d` | exact CE-FDH effect size of NCA |

## Agreement with NCA and across languages

On four synthetic datasets with 100 to 1000 cases, the native CE-FDH effect size equals the effect size of the R package 'NCA' (version 5.0.2) to within 1.2e-15, and QNCA `d_pi` at `pi = 1` differs from the exact CE-FDH value by at most 3.3e-04, the error of the outcome grid. The R, Python and Julia implementations give the same `d_pi` to within 4.4e-13 and the same p-values on shared data, tolerances and permutation streams. The result tables are in the [data folder](https://github.com/MarekDejaUJ/qnca/tree/main/data) of the repository.

## Documentation and citation

The guide to QNCA is at <https://marekdejauj.github.io/qnca/>, and the source is at <https://github.com/MarekDejaUJ/qnca>. The reference to cite is:

Deja, M. (2026). Quantile Necessary Condition Analysis: A Tolerance-Parameterized Generalization of Necessary Condition Analysis. *Quality & Quantity* (accepted for publication).

The software of release 0.4.0 is archived at Zenodo: Deja, M. (2026). *qnca: Quantile Necessary Condition Analysis*, version 0.4.0. Zenodo. <https://doi.org/10.5281/zenodo.23889752>.

## Licence

GPL-3.0-only.
