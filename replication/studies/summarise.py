"""Summaries of the robustness evaluation (outcomes defined in README.md).

Merges the Julia and NCA result files and writes one CSV per outcome to
results/summary_*.csv.

    python scripts/summarise.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
D_STAR = {"G1": 0.425, "G2": 1.0 / 3.0}
KEY = ["study", "gen", "n", "rep", "band", "level"]


def label(row) -> str:
    if row.method.startswith("QNCA") and np.isfinite(row.pi):
        return f"{row.method} {row.pi:.3g}"
    return row.method


def load() -> pd.DataFrame:
    j = pd.read_csv(RES / "robustness_julia.csv")
    r = pd.read_csv(RES / "robustness_nca.csv")
    ce = j[j.method == "CE-FDH"].merge(
        r[r.method == "NCA CE-FDH"], on=KEY, suffixes=("", "_nca"))
    diff = (ce.d - ce.d_nca).abs()
    (RES / "summary_cefdh_agreement.txt").write_text(
        f"samples = {len(ce)}\nmax_abs_diff = {diff.max():.3e}\n"
        f"samples_above_1e-12 = {int((diff > 1e-12).sum())}\n"
        f"max_abs_bottleneck_diff = {(ce.b50 - ce.b50_nca).abs().max():.3e}\n")
    if not diff.max() < 1e-6:
        raise SystemExit(f"native and NCA CE-FDH differ by {diff.max()}")
    lines = r[r.method.isin(["CR-FDH", "QR"])].copy()
    lines["slope_lost"] = lines.d.isna()
    lines["d"] = lines.d.fillna(0.0)
    lines["pi"] = np.nan
    lines["deleted"] = -1
    lines["planted_deleted"] = -1
    j["slope_lost"] = False
    out = pd.concat([j, lines[j.columns]], ignore_index=True)
    out["label"] = out.apply(label, axis=1)
    return out


def retention(df: pd.DataFrame) -> pd.DataFrame:
    s1 = df[df.study == "S1"]
    clean = s1[s1.band == "clean"][["gen", "n", "rep", "label", "d", "b50"]]
    clean = clean.rename(columns={"d": "d0", "b50": "b0"})
    # Screening is compared with the unscreened strict effect of the clean
    # sample, so deletions of genuine cases in the clean sample count.
    strict = clean[clean.label == "CE-FDH"].drop(columns="label")
    screens = sorted(l for l in s1.label.unique() if l.startswith("Screen"))
    clean = pd.concat([clean[~clean.label.isin(screens)]] +
                      [strict.assign(label=l) for l in screens], ignore_index=True)
    m = s1[s1.band != "clean"].merge(clean, on=["gen", "n", "rep", "label"])
    screened_clean = s1[(s1.band == "clean") & s1.label.isin(screens)].merge(
        clean, on=["gen", "n", "rep", "label"])
    m = pd.concat([m, screened_clean], ignore_index=True)
    m = m[m.d0 > 0]
    m["ret"] = m.d / m.d0
    m["shift50"] = m.b50 - m.b0
    g = m.groupby(["gen", "band", "n", "level", "label"])
    return pd.DataFrame({
        "mean_retention": g.ret.mean(),
        "p10_retention": g.ret.quantile(0.10),
        "share_ge_090": g.ret.apply(lambda v: float(np.mean(v >= 0.90))),
        "mean_shift_b50": g.shift50.mean(),
        "reps": g.ret.size(),
    }).reset_index()


def noise(df: pd.DataFrame) -> pd.DataFrame:
    s2 = df[df.study == "S2"]
    base = s2[s2.level == 0][["gen", "n", "rep", "label", "d"]].rename(columns={"d": "d0"})
    m = s2.merge(base, on=["gen", "n", "rep", "label"])
    m["shift"] = m.d - m.d0
    m["err"] = m.d - m.gen.map(D_STAR)
    g = m.groupby(["gen", "n", "level", "label"])
    return pd.DataFrame({
        "mean_shift": g["shift"].mean(),
        "sd_shift": g["shift"].std(),
        "mean_error": g.err.mean(),
        "rmse": g.err.apply(lambda e: float(np.sqrt(np.mean(e ** 2)))),
        "reps": g.err.size(),
    }).reset_index()


def variability(df: pd.DataFrame) -> pd.DataFrame:
    c = df[(df.study == "S1") & (df.band == "clean")]
    c = c.assign(err=c.d - c.gen.map(D_STAR))
    g = c.groupby(["gen", "n", "label"])
    return pd.DataFrame({
        "mean_d": g.d.mean(),
        "sd_d": g.d.std(),
        "bias_vs_dstar": g.err.mean(),
        "mean_b50": g.b50.mean(),
        "sd_b50": g.b50.std(),
        "reps": g.d.size(),
    }).reset_index()


def screen(df: pd.DataFrame) -> pd.DataFrame:
    s = df[(df.study == "S1") & df.label.str.startswith("Screen")]
    s = s.assign(genuine=s.deleted - s.planted_deleted,
                 planted_found=np.where(s.level > 0, s.planted_deleted / s.level.replace(0, np.nan), np.nan))
    g = s.groupby(["gen", "band", "n", "level", "label"])
    return pd.DataFrame({
        "mean_deleted": g.deleted.mean(),
        "mean_genuine_deleted": g.genuine.mean(),
        "mean_share_planted_found": g.planted_found.mean(),
    }).reset_index()


def slope_loss(df: pd.DataFrame) -> pd.DataFrame:
    s = df[df.label.isin(["CR-FDH", "QR"])]
    return s.groupby(["study", "gen", "n", "band", "level", "label"]).slope_lost.sum().reset_index()


def main() -> int:
    df = load()
    retention(df).to_csv(RES / "summary_retention.csv", index=False, float_format="%.4f")
    noise(df).to_csv(RES / "summary_noise.csv", index=False, float_format="%.4f")
    variability(df).to_csv(RES / "summary_variability.csv", index=False, float_format="%.4f")
    screen(df).to_csv(RES / "summary_screen.csv", index=False, float_format="%.4f")
    slope_loss(df).to_csv(RES / "summary_slope_loss.csv", index=False)
    print("summaries written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
