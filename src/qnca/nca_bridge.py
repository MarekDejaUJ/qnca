"""Bridge between QNCA and Dul's NCA, for the publication comparison.

Two comparison paths are provided, and they answer the same question at two
levels of trust:

1. :func:`compare_native` -- pure Python. It pairs ``qnca(pi=1)`` with the native
   CE-FDH effect size in :func:`qnca.core.nca_ce_fdh_d`. This always runs and is
   what the validation invariant test relies on, because it needs no R toolchain.

2. :func:`run_r_nca` -- shells out to the real ``NCA`` R package via ``Rscript``.
   This is the authoritative external check used for the published tables. It runs
   only where R and the NCA package are installed; otherwise it reports cleanly
   that it was skipped, so a CI run without R still succeeds.

The design separation is deliberate: the native path checks the *mathematics*
of the reduction; the R path checks *agreement with the reference software a
reviewer already trusts*. Both should be reported.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np

from .core import nca_ce_fdh_d, qnca

__all__ = ["compare_native", "R_BRIDGE_TEMPLATE", "write_r_bridge", "run_r_nca"]


def compare_native(
    X: np.ndarray,
    Y: np.ndarray,
    scope: tuple[float, float, float, float] | None = None,
    n_grid: int = 200,
    tol: float = 1e-2,
) -> dict[str, float | bool]:
    """Compare QNCA(pi=1) against the native CE-FDH effect size.

    A fine grid (default 200) is used because the QNCA effect size is a grid
    approximation of the continuous empty area, while the native CE-FDH value is
    exact on the data; refining the grid drives the two together. The returned
    ``within_tol`` is the validation invariant the test suite asserts.
    """
    res = qnca(X, Y, pi=1.0, n_grid=n_grid, B=None, scope=scope)
    d_ce_fdh = nca_ce_fdh_d(X, Y, scope=scope)
    diff = abs(res.d_pi - d_ce_fdh)
    return {
        "d_pi_pi1": res.d_pi,
        "d_ce_fdh": d_ce_fdh,
        "abs_diff": diff,
        "within_tol": bool(diff <= tol),
        "tol": tol,
        "n_grid": n_grid,
    }


# R script template. Uses the documented NCA 5.x interface: nca_analysis() with
# ceilings, scope, bottleneck.x/y, steps and test.rep, then nca_extract() to pull
# the effect size and p-value per ceiling. Datasets are read from CSV so the same
# simulated data feeds QNCA and NCA.
R_BRIDGE_TEMPLATE = r"""
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript nca_dul_bridge.R <data_dir> [test_rep]")
data_dir <- args[[1L]]
test_rep <- if (length(args) >= 2L) as.integer(args[[2L]]) else {TEST_REP}

if (!requireNamespace("NCA", quietly = TRUE)) {
  stop("R package 'NCA' is not installed. Run install.packages('NCA') first.")
}
library(NCA)
set.seed(42L)

repo_root <- normalizePath(file.path(data_dir, ".."), mustWork = TRUE)
source(file.path(repo_root, "r", "R", "qnca.R"))

sizes <- c({SIZES})
ceilings_use <- c("ce_fdh", "cr_fdh")
rows <- list(); idx <- 1L
r_version <- paste(R.version$major, R.version$minor, sep = ".")
nca_version <- as.character(utils::packageVersion("NCA"))

for (n in sizes) {
  f <- file.path(data_dir, paste0("qnca_dataset_n", n, ".csv"))
  d <- read.csv(f)
  scope_use <- c(0, 100, 0, 100)
  q1 <- qnca(d$X, d$Y, pi = 1, n_grid = 400L, B = NA, scope = scope_use)
  native <- nca_ce_fdh_d(d$X, d$Y, scope = scope_use)
  model <- nca_analysis(
    data = d, x = "X", y = "Y",
    ceilings = ceilings_use,
    scope = scope_use,
    bottleneck.x = "actual", bottleneck.y = "actual",
    steps = seq(0, 100, 10),
    test.rep = test_rep
  )
  for (ceiling in ceilings_use) {
    eff  <- tryCatch(nca_extract(model, x = "X", ceiling = ceiling, param = "Effect size"),
                     error = function(e) NA_real_)
    pval <- tryCatch(nca_extract(model, x = "X", ceiling = ceiling, param = "p-value"),
                     error = function(e) NA_real_)
    rows[[idx]] <- data.frame(
      dataset = paste0("n", n),
      n = n,
      method = ceiling,
      qnca_d_pi1 = q1$d_pi,
      native_ce_fdh_d = native,
      dul_d = as.numeric(eff),
      dul_p = as.numeric(pval),
      abs_diff_qnca_vs_dul = abs(q1$d_pi - as.numeric(eff)),
      abs_diff_native_vs_dul = abs(native - as.numeric(eff)),
      test_rep = test_rep,
      qnca_n_grid = 400L,
      r_version = r_version,
      nca_version = nca_version
    )
    idx <- idx + 1L
  }
}
out <- do.call(rbind, rows)
write.csv(out, file.path(data_dir, "dul_nca_results.csv"), row.names = FALSE)
print(out)
"""


def write_r_bridge(path: Path, sizes: list[int], test_rep: int) -> Path:
    """Materialise the R bridge script with concrete sizes and test.rep."""
    code = (
        R_BRIDGE_TEMPLATE.replace("{TEST_REP}", str(int(test_rep)))
        .replace("{SIZES}", ", ".join(str(int(s)) for s in sizes))
        .strip()
        + "\n"
    )
    path.write_text(code)
    return path


def run_r_nca(data_dir: Path, sizes: list[int], test_rep: int = 1999) -> tuple[bool, str]:
    """Run Dul's NCA via Rscript if available; report cleanly if not.

    Returns (ok, message). ``ok`` is False (not an exception) when R or the NCA
    package is missing, so callers and CI can treat 'no R here' as an expected,
    non-fatal outcome.
    """
    rscript = shutil.which("Rscript")
    if rscript is None:
        return False, "Rscript not found in PATH; R NCA comparison skipped."
    data_dir.mkdir(parents=True, exist_ok=True)
    script = write_r_bridge(data_dir / "nca_dul_bridge.R", sizes, test_rep)
    try:
        proc = subprocess.run(
            [rscript, str(script), str(data_dir), str(int(test_rep))],
            capture_output=True,
            text=True,
            timeout=3600,
            check=False,
        )
    except Exception as exc:  # pragma: no cover - environment dependent
        return False, f"Failed to launch Rscript: {exc}"
    (data_dir / "dul_nca_r_bridge.log").write_text(
        "STDOUT\n======\n" + proc.stdout + "\n\nSTDERR\n======\n" + proc.stderr
    )
    if proc.returncode != 0:
        return False, f"Rscript exited {proc.returncode}; see dul_nca_r_bridge.log."
    return True, f"R NCA comparison written to {data_dir / 'dul_nca_results.csv'}."
