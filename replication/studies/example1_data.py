"""Write the Example 1 sample: the manuscript's floor generator, n = 300, seed 47.

The seed is the study seed 42 plus 5, as in the stored planted-case run; the
planted case is (2, 95).

    python scripts/example1_data.py
"""

from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT.parent / "python" / "src"))
from qnca import generate_reverse_l  # noqa: E402

X, Y = generate_reverse_l(300, rng=np.random.default_rng(47))
out = ROOT / "data" / "example1.csv"
with out.open("w") as f:
    f.write("X,Y,planted\n")
    for x, y in zip(X, Y):
        f.write(f"{float(x)!r},{float(y)!r},0\n")
    f.write("2.0,95.0,1\n")
print(f"wrote {out}")
