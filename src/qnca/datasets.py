"""Synthetic data generators for QNCA demonstrations and tests.

The reverse-L (or 'ceiling') design is the canonical necessity scatter: points
fill the lower-right triangle and the upper-left corner stays empty, because a
high outcome requires a high condition but a high condition does not guarantee a
high outcome. This is the shape NCA and QNCA are built to detect.
"""

from __future__ import annotations

import numpy as np

__all__ = ["generate_reverse_l"]


def generate_reverse_l(
    n: int,
    ci: float = 15.0,
    cs: float = 0.85,
    k: float = 2.0,
    noise: float = 4.0,
    rng: np.random.Generator | int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a reverse-L necessity scatter.

    The condition X is uniform on [0, 100]. The ceiling is a line ``ci + cs * X``;
    the realised outcome is a fraction ``U**k`` of that ceiling (U uniform), plus
    Gaussian noise. Raising ``k`` pushes points toward the floor and sharpens the
    empty upper-left corner; ``noise`` softens the boundary. Both variables are
    clamped to [0, 100] so the scope rectangle is fixed and comparable across n.

    Parameters
    ----------
    n : int
        Number of cases.
    ci, cs : float
        Intercept and slope of the deterministic ceiling.
    k : float
        Efficiency exponent; larger k means fewer cases sit near the ceiling.
    noise : float
        Standard deviation of additive Gaussian noise.
    rng : Generator, int, or None
        Random source. An int is used as a seed.
    """
    if not isinstance(rng, np.random.Generator):
        rng = np.random.default_rng(rng)
    X = rng.uniform(0.0, 100.0, n)
    ceiling = ci + cs * X
    eff = rng.uniform(0.0, 1.0, n) ** k
    Y = eff * ceiling + rng.normal(0.0, noise, n)
    return np.clip(X, 0.0, 100.0), np.clip(Y, 0.0, 100.0)
