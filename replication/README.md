# Replication

The folder holds the scripts and the results behind every number and figure of the article
"Quantile Necessary Condition Analysis: A Tolerance-Parameterized Generalization of Necessary
Condition Analysis" (Quality & Quantity). The scripts load the R, Python and Julia packages of
release 0.4.0 from the folders `r/`, `python/` and `julia/QuantileNecessity.jl/` of the
repository, and the stored tables of the validation studies are in `data/` at the repository root.

## Contents

- `run_all.sh`: regenerates the stored studies in `data/` (strict-limit validation, tolerance
  grid, bottleneck comparison, cross-language agreement, reference-band checks, calibration and
  power screens) and the four figures in `figures/`.
- `scripts/`: the Python drivers of `run_all.sh`; the R and Julia parts run through `Rscript` and
  `julia` on the same inputs and permutation streams.
- `R/figures/`: the figure scripts of the validation studies.
- `figures/`: the four figures and the bottleneck table of the validation studies.
- `studies/robustness.jl`: generators and the CE-FDH, QNCA, order-m and screen-and-delete
  results for the added-case, noise, variability and tolerance-path studies.
- `studies/robustness_nca.R`: CE-FDH, CR-FDH and quantile-regression ceilings and the outlier
  screen of the NCA package (5.0.2) on the same samples.
- `studies/summarise.py`, `studies/tables.py`, `studies/figures.R`, `studies/shape_check.R`:
  summaries, the retention table, the figures of the planted case and of the comparison, and the
  required condition of the quantile-regression ceiling at three outcome targets.
- `studies/application.R`: technology-acceptance application (cIPMA scores, NCA ceilings and
  outlier screen, QNCA effects, permutation tests, reference band, resolution and deletion
  tables).
- `studies/example1_data.py`: the sample of Example 1, stored in `data/example1.csv`.
- `results/`: results, summaries, figures and version records of the studies.

## Evaluation design

All samples lie in the fixed scope `[0, 100]^2`, which is also the effect-size
scope; `U` is uniform on `(0, 1)` and values are clamped to the scope.

Generators. G1, linear ceiling: `X ~ U(0, 100)`, `Y = U^2 (15 + 0.85 X) + e`,
`e ~ N(0, 4^2)`; the noise-free boundary is `x = (y - 15) / 0.85` for `y > 15`
and its empty-zone effect is `d* = 0.425`. G2, concave ceiling: `X ~ U(0, 100)`,
`Y = U^2 10 sqrt(X) + e`, `e ~ N(0, 4^2)`; the boundary is `x = y^2 / 100` and
`d* = 1/3`.

S1, added efficient cases: `n` in {100, 300, 1000}, 500 replicates per
generator and size. Each replicate draws one clean sample and five planted
positions per band, with `x ~ U(0, 10)` and `y ~ U(40, 70)` (middle band) or
`y ~ U(80, 100)` (top band). Contaminated samples add the first `m` positions,
`m` in {1, 2, 3, 5}, so every contaminated sample is paired with its clean
sample.

S2, measurement noise: noise-free G1 and G2, `n` in {100, 300}, 500
replicates; observed values `X + u`, `Y + v` with independent
`u, v ~ N(0, sigma^2)`, `sigma` in {0, 2, 5, 10}, the same base sample and the
same standard-normal draws scaled by `sigma` at every level.

S3, sampling variability: the clean samples of S1 give the Monte Carlo mean
and standard deviation of the effect and of the bottleneck at `y = 50`.

S4, tolerance path: QNCA with the envelope at `pi` in {1, 0.99, 0.975, 0.95,
0.925, 0.90, 0.85, 0.80} on the clean and the `m` in {1, 3, 5} samples of S1.

Methods. CE-FDH (exact step area); CR-FDH (NCA ceiling regression through the
CE-FDH corner points); QR (NCA quantile-regression ceiling, `qr.tau = 0.95`);
screen and delete (remove the case with the largest single-deletion CE-FDH
gain, recompute, and repeat while the gain is at least 1% of the current
effect, for at most ten deletions; a second variant deletes exactly as many
cases as were added); order-m (at each target, the expected minimum of `m`
draws with replacement from the attaining condition values, `m` in {19, 9},
made monotone by the envelope); QNCA with the envelope and with least-squares
pooling at `pi` in {0.95, 0.90}. Grid-based methods use 200 equally spaced
outcome targets and the fixed-scope carry-forward. Line ceilings use NCA's
area of the scope above the line; a quantile-regression line without a
positive slope describes no upper-left empty zone, so its effect is recorded as
0 and its bottleneck at `y = 50` as the scope minimum when the intercept is at
least 50. The bottleneck at `y = 50` is the frontier at that target, or
`(50 - intercept) / slope` for a line, clamped to the scope.

Outcomes. Retention (contaminated effect over the clean effect of the same
replicate: mean, 10th percentile, share at or above 0.90); bottleneck shift at
`y = 50`; noise shift (effect at `sigma` minus effect at `sigma = 0` for the
same base sample); error against `d*` (mean and root mean squared error);
Monte Carlo SD of the effect and of the bottleneck across clean replicates;
screen detection (share of planted cases among the deleted cases and number of
genuine cases deleted).

Application. The data are the survey of French adopters of e-book readers (n = 174) of Richter, N.F., Hauff, S., Kolev, A.E., Schubring, S. (2023). Dataset on an extended technology acceptance model: a combined application of PLS-SEM and NCA. Data in Brief, 48, 109190, https://doi.org/10.1016/j.dib.2023.109190. The file is "Extended TAM.csv" of Mendeley Data, https://doi.org/10.17632/pd5dp3phx2.4, under the CC BY 4.0 licence. `studies/application.R` downloads the file to `data/tam.csv` and checks its SHA-256 checksum (48930f478d415c2724135493c77b101a4aa1954795d8dcc035ec835bd600c7de); the repository keeps neither the file nor the scores computed from the file.

The rescaled latent scores are rebuilt from the data with the published cIPMA workflow (`seminr` 2.4.2, `cIPMA` 0.1.6). For the five direct antecedents of use, each pair keeps the empirical scope of its rescaled scores. The NCA package (5.0.2) gives the CE-FDH, CR-FDH and QR effects and the outlier screen with `k` in {1, 2}. QNCA uses 400 outcome targets, `pi` in {1, 0.95, 0.90}, permutation p-values with `B = 1999` (seed 20260930 plus the pair index), the reference band with `M = 999` at all three tolerances (seed 20261030 plus the pair index), and `qnca_outliers` with `k` in {1, 2}.

Seeds. The base seed of the simulation is 20260930, and the hash in the
scripts derives each cell seed from the base seed, the generator, the sample
size, the band and the replicate. Every run writes CSV results and a version record to
`results/`.

## Reproduction

The commands run from the root of the repository.

```bash
julia -t auto replication/studies/robustness.jl 500
Rscript replication/studies/robustness_nca.R 8
python replication/studies/summarise.py
python replication/studies/tables.py
Rscript replication/studies/shape_check.R
Rscript replication/studies/application.R
python replication/studies/example1_data.py
Rscript replication/studies/figures.R
PYTHON=python bash replication/run_all.sh
```

Requirements: Julia 1.9 or later; R 4.5 or later with NCA 5.0.2, seminr 2.4.2, cIPMA 0.1.6 and
ggplot2; Python 3.9 or later with numpy and pandas. The application needs seminr 2.4.2 and cIPMA
0.1.6:

```r
remotes::install_version("seminr", version = "2.4.2")
remotes::install_github("MarekDejaUJ/cIPMA@2e240a3a183db4216742e241ca7d68db7d6356de")
```

Seeds are fixed in the scripts; `results/*_provenance.txt` records versions and run times.

On Windows, `robustness_nca.R` runs with one core (`Rscript replication/studies/robustness_nca.R 1`),
because forked parallel processing is unavailable there. Runs under R 4.6.1 agree with the stored
results to floating point. Where cases tie to the fifteenth significant digit, the order of the
tied rows in the deletion screen can change, and a tail probability of the reference band can move
by one or two of the 999 null replicates.
