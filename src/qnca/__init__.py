"""Quantile Necessary Condition Analysis (QNCA).

A tolerance-parameterised generalisation of Dul's Necessary Condition Analysis.
The deterministic NCA ceiling is the pi == 1 limit of the QNCA frontier family.
"""

from .core import (
    ConsistencyProbeResult,
    QNCAResult,
    SpuriousnessBandResult,
    consistency_probe,
    isotonic_increasing,
    nca_ce_fdh_d,
    permutation_test,
    qnca,
    qnca_d,
    qnca_frontier,
    qnca_rank,
    quantile_type1,
    quantile_type1_pi,
    spuriousness_band,
)
from .datasets import generate_reverse_l

__version__ = "0.3.2"

__all__ = [
    "QNCAResult",
    "SpuriousnessBandResult",
    "ConsistencyProbeResult",
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
