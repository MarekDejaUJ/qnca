"""Quantile Necessary Condition Analysis (QNCA).

A tolerance-parameterised generalisation of Dul's Necessary Condition Analysis.
The deterministic NCA ceiling is the pi == 1 member of the tolerance-indexed
QNCA frontiers.
"""

from .core import (
    ConsistencyProbeResult,
    OutlierScreenResult,
    QNCAResult,
    ResolutionResult,
    SensitivityResult,
    SpuriousnessBandResult,
    consistency_probe,
    isotonic_increasing,
    monotone_envelope,
    nca_ce_fdh_d,
    permutation_test,
    qnca,
    qnca_d,
    qnca_frontier,
    qnca_outliers,
    qnca_rank,
    qnca_resolution,
    qnca_sensitivity,
    quantile_type1,
    quantile_type1_pi,
    spuriousness_band,
)
from .datasets import generate_reverse_l

__version__ = "0.4.0"

__all__ = [
    "QNCAResult",
    "SpuriousnessBandResult",
    "ConsistencyProbeResult",
    "ResolutionResult",
    "OutlierScreenResult",
    "SensitivityResult",
    "monotone_envelope",
    "qnca_resolution",
    "qnca_outliers",
    "qnca_sensitivity",
    "qnca",
    "qnca_frontier",
    "qnca_d",
    "permutation_test",
    "spuriousness_band",
    "consistency_probe",
    "quantile_type1",
    "quantile_type1_pi",
    "qnca_rank",
    "isotonic_increasing",
    "nca_ce_fdh_d",
    "generate_reverse_l",
    "__version__",
]
