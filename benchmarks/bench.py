"""Timing of twowayfeweights on synthetic panels of increasing size.

    python benchmarks/bench.py
"""

from __future__ import annotations

import time
import warnings

import numpy as np
import pandas as pd

from twowayfeweights import twowayfeweights


def panel(n_g: int, n_t: int, n_ctrl: int = 5, obs_per_cell: int = 1, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    g = np.repeat(np.arange(n_g), n_t * obs_per_cell)
    t = np.tile(np.repeat(np.arange(n_t), obs_per_cell), n_g)
    adopt = rng.integers(1, n_t + 5, n_g)[g]
    D = (t >= adopt).astype(float)
    df = pd.DataFrame({"g": g, "t": t, "D": D, "Y": D * rng.normal(1, 1, n_g)[g] + 0.1 * t + rng.normal(size=g.size)})
    for j in range(n_ctrl):
        df[f"x{j}"] = rng.normal(size=g.size)
    df["w"] = rng.uniform(0.5, 1.5, g.size)
    df = df.sort_values(["g", "t"])
    df["dY"] = df.groupby("g")["Y"].diff()
    df["dD"] = df.groupby("g")["D"].diff()
    return df


def timeit(f, repeat=3):
    best = np.inf
    for _ in range(repeat):
        t0 = time.perf_counter()
        f()
        best = min(best, time.perf_counter() - t0)
    return best


if __name__ == "__main__":
    warnings.simplefilter("ignore")
    print(f"{'rows':>10} {'groups':>7} {'periods':>7} {'feTR':>8} {'feS':>8} {'fdTR':>8} {'fdS':>8}  (seconds, 5 controls + weights)")
    for n_g, n_t in [(1_000, 10), (10_000, 20), (50_000, 20), (100_000, 30)]:
        df = panel(n_g, n_t)
        ctrl = [f"x{j}" for j in range(5)]
        kw = dict(controls=ctrl, weights="w")
        times = [
            timeit(lambda: twowayfeweights(df, "Y", "g", "t", "D", type="feTR", **kw)),
            timeit(lambda: twowayfeweights(df, "Y", "g", "t", "D", type="feS", **kw)),
            timeit(lambda: twowayfeweights(df, "dY", "g", "t", "dD", type="fdTR", D0="D", **kw)),
            timeit(lambda: twowayfeweights(df, "dY", "g", "t", "dD", type="fdS", **kw)),
        ]
        print(f"{len(df):>10,} {n_g:>7,} {n_t:>7} " + " ".join(f"{x:8.3f}" for x in times))
