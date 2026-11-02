"""LaTeX rows for the robustness table, built from results/summary_retention.csv.

    python scripts/tables.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"

METHODS = [
    ("CE-FDH", "CE-FDH"),
    ("CR-FDH", "CR-FDH"),
    ("QR", "QR ceiling, \\(\\tau=0.95\\)"),
    ("Screen 1% rule", "Screen and delete"),
    ("Order-m 19", "Order-\\(m\\), \\(m=19\\)"),
    ("QNCA least squares 0.95", "QNCA \\(0.95\\), least squares"),
    ("QNCA envelope 0.95", "QNCA \\(0.95\\), envelope"),
    ("QNCA envelope 0.9", "QNCA \\(0.90\\), envelope"),
]
CELLS = [(g, b, m) for g in ("G1", "G2") for b in ("mid", "top") for m in (1, 3)]


def main() -> int:
    r = pd.read_csv(RES / "summary_retention.csv")
    r = r[r.n == 300]
    lines = []
    for key, name in METHODS:
        vals = []
        for g, b, m in CELLS:
            v = r[(r.gen == g) & (r.band == b) & (r.level == m) & (r.label == key)].mean_retention
            if len(v) != 1:
                raise SystemExit(f"missing cell {key} {g} {b} {m}")
            vals.append(f"{float(v.iloc[0]):.2f}")
        lines.append(name + " & " + " & ".join(vals) + " \\\\")
    out = RES / "table_retention_n300.tex"
    out.write_text("\n".join(lines) + "\n")
    print(out.read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
