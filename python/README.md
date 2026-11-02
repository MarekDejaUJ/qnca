# qnca

Quantile Necessary Condition Analysis (QNCA) for Python.

Quantile Necessary Condition Analysis extends Necessary Condition Analysis (NCA) with a tolerance. For each outcome target, the QNCA frontier is a low type-1 quantile of the condition among the cases that attain the target. The tolerance `pi` in (0, 1] sets the share `1 - pi` of attaining cases allowed below the reported requirement, and the analyst fixes `pi` before the analysis.

At `pi = 1`, the QNCA frontier equals the CE-FDH ceiling of NCA. At `pi < 1`, each requirement rests on a stated order statistic of the attaining cases, so each requirement absorbs a known number of unusually efficient cases. QNCA keeps the bottleneck table, the empty-zone effect size and the permutation test of NCA.

## Installation

The release installs from PyPI:

```bash
python -m pip install qnca
```

The development version installs from the `python` folder of the repository:

```bash
python -m pip install "git+https://github.com/MarekDejaUJ/qnca.git#subdirectory=python"
```

The package needs Python 3.9 or later and NumPy.

## Example

`generate_reverse_l()` draws a reverse-L scatter, the pattern of a necessary condition. `qnca()` estimates the frontier at `pi = 0.95`, the effect size `d_pi`, the permutation p-value `p_pi` from 1999 permutations, and the bottleneck table.

```python
from qnca import generate_reverse_l, qnca

X, Y = generate_reverse_l(500, rng=42)
res = qnca(X, Y, pi=0.95, B=1999, scope=(0, 100, 0, 100), seed=1)
print(round(res.d_pi, 6), res.p_pi)
# 0.503332 0.0005
print(res.bottleneck[:3].round(6))
# [[0.       5.338193]
#  [2.040816 5.539541]
#  [4.081633 8.444432]]
```

A p-value of 0.0005 is the smallest value that 1999 permutations can report: every permuted effect size stayed below the observed effect size. Each row of the bottleneck table gives an outcome level and the condition level required for the outcome level at the chosen tolerance.

## Functions

| Function | Result |
|---|---|
| `qnca()` | effect size, permutation p-value and bottleneck table at one tolerance |
| `qnca_frontier()` | fitted frontier on a grid of outcome levels |
| `qnca_resolution()` | per outcome target: conditioning-set size, selected rank, raw and fitted requirement, inherited flag, share of attaining cases below the requirement, resistance |
| `qnca_outliers()` | change of `d_pi` across tolerances after the deletion of one case or of `k` cases |
| `qnca_sensitivity()` | `d_pi` after the addition of cases at declared positions |
| `spuriousness_band()` | proposed dependence-matched reference band for `d_pi` |
| `consistency_probe()` | `d_pi` across subsample sizes at two tolerances |
| `nca_ce_fdh_d()` | exact CE-FDH effect size of NCA |
| `qnca.nca_bridge.run_r_nca()` | QNCA at `pi = 1` beside the CE-FDH results of the R package 'NCA', where R and 'NCA' are installed |

## Agreement with NCA and across languages

On four synthetic datasets with 100 to 1000 cases, the native CE-FDH effect size equals the effect size of the R package 'NCA' (version 5.0.2) to within 1.2e-15, and QNCA `d_pi` at `pi = 1` differs from the exact CE-FDH value by at most 3.3e-04, the error of the outcome grid. The R, Python and Julia implementations give the same `d_pi` to within 4.4e-13 and the same p-values on shared data, tolerances and permutation streams. The result tables are in the [data folder](https://github.com/MarekDejaUJ/qnca/tree/main/data) of the repository.

## Documentation and citation

The guide to QNCA is at <https://marekdejauj.github.io/qnca/>, and the source is at <https://github.com/MarekDejaUJ/qnca>. The reference to cite is:

Deja, M. (2026). Quantile Necessary Condition Analysis: A Tolerance-Parameterized Generalization of Necessary Condition Analysis. *Quality & Quantity* (accepted for publication).

## Licence

GPL-3.0-only.
