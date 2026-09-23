"""Test cases shared by the reference generators (Stata / R) and the parity tests."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).parent / "data"
FIXTURES = Path(__file__).parent / "fixtures"

N_STYR = 683  # state-year dummies in the Gentzkow et al. (2011) data


def _case(name, data, Y, G, T, D, type, D0=None, controls=(), weights=None, ots=(), trw=(), beta_only=False):
    """``beta_only=True``: Stata's own output is not deterministic (it depends on Stata's random sort
    tie-breaking), so only beta is compared.
    """
    return dict(
        name=name, data=data, Y=Y, G=G, T=T, D=D, type=type, D0=D0,
        controls=list(controls), weights=weights, other_treatments=list(ots), test_random_weights=list(trw),
        beta_only=beta_only,
    )


STYR = [f"styr{i}" for i in range(1, N_STYR + 1)]

CASES = [
    # --- wagepan (Stata/R help-file examples + extensions)
    _case("wp_feTR", "wagepan", "lwage", "nr", "year", "union", "feTR", trw=["educ"]),
    _case("wp_feTR_ctrl_w", "wagepan", "lwage", "nr", "year", "union", "feTR", controls=["hours", "married"], weights="wt", trw=["educ", "exper"]),
    _case("wp_feTR_ot", "wagepan", "lwage", "nr", "year", "union", "feTR", ots=["married"], trw=["educ"]),
    _case("wp_feTR_ot_ctrl_w", "wagepan", "lwage", "nr", "year", "union", "feTR", controls=["hours"], weights="wt", ots=["married", "south"], trw=["educ"]),
    _case("wp_fdTR", "wagepan", "diff_lwage", "nr", "year", "diff_union", "fdTR", D0="union", trw=["educ"]),
    _case("wp_fdTR_ctrl_w", "wagepan", "diff_lwage", "nr", "year", "diff_union", "fdTR", D0="union", controls=["hours"], weights="wt", trw=["educ"]),
    _case("wp_feS", "wagepan", "lwage", "nr", "year", "union", "feS", trw=["educ"]),
    _case("wp_feS_ctrl_w", "wagepan", "lwage", "nr", "year", "union", "feS", controls=["hours", "married"], weights="wt", trw=["educ"]),
    _case("wp_fdS", "wagepan", "diff_lwage", "nr", "year", "diff_union", "fdS", trw=["educ"]),
    _case("wp_fdS_ctrl_w", "wagepan", "diff_lwage", "nr", "year", "diff_union", "fdS", controls=["hours"], weights="wt", trw=["educ"]),
    # --- Gentzkow, Shapiro & Sinkinson (2011): did_book chapter 5
    _case("gz_feTR", "gentzkow", "prestout", "cnty90", "year", "numdailies", "feTR"),
    _case("gz_feTR_styr", "gentzkow", "prestout", "cnty90", "year", "numdailies", "feTR", controls=STYR),
    _case("gz_fdTR", "gentzkow", "changeprestout", "cnty90", "year", "changedailies", "fdTR", D0="numdailies", trw=["year"]),
    _case("gz_fdTR_styr", "gentzkow", "changeprestout", "cnty90", "year", "changedailies", "fdTR", D0="numdailies", controls=STYR),
    _case("gz_fdTR_styr_rw", "gentzkow", "changeprestout", "cnty90", "year", "changedailies", "fdTR", D0="numdailies", controls=STYR, trw=["year"]),
    _case("gz_feS", "gentzkow", "prestout", "cnty90", "year", "numdailies", "feS", trw=["year"]),
    _case("gz_fdS", "gentzkow", "changeprestout", "cnty90", "year", "changedailies", "fdS", trw=["year"]),
    _case("gz_feTR_w", "gentzkow", "prestout", "cnty90", "year", "numdailies", "feTR", weights="share_urb_baseline", trw=["year"]),
    _case("gz_feS_styr", "gentzkow", "prestout", "cnty90", "year", "numdailies", "feS", controls=STYR, trw=["year"]),
    _case("gz_fdS_styr", "gentzkow", "changeprestout", "cnty90", "year", "changedailies", "fdS", controls=STYR, trw=["year"]),
    _case("gz_feS_w", "gentzkow", "prestout", "cnty90", "year", "numdailies", "feS", weights="share_urb_baseline", trw=["year"]),
    _case("gz_fdS_w", "gentzkow", "changeprestout", "cnty90", "year", "changedailies", "fdS", weights="share_urb_baseline", trw=["year"]),
    # --- simulated micro data: several obs per cell, within-cell variation, holes, missing values
    _case("sim_feTR", "sim", "Y", "g", "t", "D", "feTR", controls=["x1", "x2"], weights="w", trw=["z"]),
    _case("sim_feTR_ot", "sim", "Y", "g", "t", "D", "feTR", controls=["x1"], weights="w", ots=["D2"], trw=["z"]),
    _case("sim_fdTR", "sim", "dY", "g", "t", "dD", "fdTR", D0="D", controls=["x1"], weights="w", trw=["z"]),
    # Stata's feS on several obs per cell depends on its (random) sort order within cells.
    _case("sim_feS", "sim", "Y", "g", "t", "D", "feS", controls=["x1"], weights="w", trw=["z"], beta_only=True),
    _case("sim_fdS", "sim", "dY", "g", "t", "dD", "fdS", controls=["x1"], weights="w", trw=["z"]),
    _case("sim_feTR_now", "sim", "Y", "g", "t", "D", "feTR"),
    _case("sim_fdTR_now", "sim", "dY", "g", "t", "dD", "fdTR", D0="D"),
    _case("sim_fdS_now", "sim", "dY", "g", "t", "dD", "fdS"),
    # --- same panel collapsed to one observation per (g,t) cell (unbalanced, with holes)
    _case("simcell_feTR", "simcell", "Y", "g", "t", "D", "feTR", controls=["x1", "x2"], weights="w", trw=["z"]),
    _case("simcell_feS", "simcell", "Y", "g", "t", "D", "feS", controls=["x1"], weights="w", trw=["z"]),
    _case("simcell_feS_now", "simcell", "Y", "g", "t", "D", "feS"),
    _case("simcell_fdTR", "simcell", "dY", "g", "t", "dD", "fdTR", D0="D", weights="w", trw=["z"]),
    _case("simcell_fdS", "simcell", "dY", "g", "t", "dD", "fdS", weights="w", trw=["z"]),
    _case("simcell_fdS_now", "simcell", "dY", "g", "t", "dD", "fdS", controls=["x1", "x2"]),
]


def make_sim(seed: int = 12345) -> pd.DataFrame:
    """Unbalanced individual-level panel; values are float32-representable so Stata stores them exactly."""
    rng = np.random.default_rng(seed)
    rows = []
    n_g, n_t = 60, 8
    for g in range(1, n_g + 1):
        adopt = rng.integers(3, n_t + 3)  # some groups never treated
        periods = [t for t in range(1, n_t + 1) if not (g % 7 == 0 and t == 4)]  # holes
        for t in periods:
            for _ in range(rng.integers(1, 4)):
                rows.append((g, t, adopt))
    df = pd.DataFrame(rows, columns=["g", "t", "adopt"])
    n = len(df)
    base = (df["t"] >= df["adopt"]).astype(float)
    # within-cell variation of the treatment in some groups
    flip = (df["g"] % 5 == 0) & (df["t"] >= df["adopt"]) & (rng.random(n) < 0.4)
    df["D"] = np.where(flip, 0.0, base) * np.round(1 + rng.random(n), 2)
    df["D2"] = (rng.random(n) < 0.3).astype(float)
    df["x1"] = np.round(rng.normal(size=n), 3)
    df["x2"] = np.round(df["g"] % 3 + 0.5 * df["t"], 3)
    df["z"] = np.round(rng.normal(size=n) + 0.1 * df["t"], 3)
    df["w"] = np.round(0.5 + rng.random(n), 3)
    df["Y"] = np.round(0.3 * df["D"] + 0.1 * df["t"] + rng.normal(size=n), 4)
    # first differences at the (g,t) level
    cell = df.groupby(["g", "t"])[["Y", "D"]].transform("mean")
    df["_Ym"], df["_Dm"] = cell["Y"], cell["D"]
    cm = df.groupby(["g", "t"], as_index=False)[["_Ym", "_Dm"]].first().sort_values(["g", "t"])
    prev_ok = (cm["g"] == cm["g"].shift()) & (cm["t"] == cm["t"].shift() + 1)
    cm["dY"] = np.where(prev_ok, cm["_Ym"] - cm["_Ym"].shift(), np.nan)
    cm["dD"] = np.where(prev_ok, cm["_Dm"] - cm["_Dm"].shift(), np.nan)
    df = df.merge(cm[["g", "t", "dY", "dD"]], on=["g", "t"], how="left").drop(columns=["_Ym", "_Dm", "adopt"])
    # scattered missing values
    df.loc[rng.random(n) < 0.02, "Y"] = np.nan
    df.loc[rng.random(n) < 0.02, "x1"] = np.nan
    df.loc[rng.random(n) < 0.01, "w"] = np.nan
    for c in ["D", "D2", "x1", "x2", "z", "w", "Y", "dY", "dD"]:
        df[c] = df[c].astype(np.float32).astype(np.float64)
    return df


@lru_cache(maxsize=None)
def load(name: str) -> pd.DataFrame:
    if name == "wagepan":
        from twowayfeweights import load_wagepan

        df = load_wagepan()
        rng = np.random.default_rng(2020)
        df["wt"] = np.round(0.5 + rng.random(len(df)), 3).astype(np.float32).astype(np.float64)
        for c in ["lwage", "diff_lwage", "hours"]:
            df[c] = df[c].astype(np.float32).astype(np.float64)
        return df
    if name == "gentzkow":
        df = pd.read_csv(DATA_DIR / "gentzkow.csv.gz")
        levels = np.sort(df["styr"].dropna().unique())
        assert len(levels) == N_STYR
        dummies = (df["styr"].to_numpy()[:, None] == levels[None, :]).astype(np.float64)
        return pd.concat([df, pd.DataFrame(dummies, columns=STYR, index=df.index)], axis=1)
    if name == "sim":
        return pd.read_csv(DATA_DIR / "sim.csv.gz")
    if name == "simcell":
        df = load("sim").dropna(subset=["Y", "x1", "w"])
        out = df.groupby(["g", "t"], as_index=False).agg(
            Y=("Y", "mean"), D=("D", "mean"), x1=("x1", "mean"), x2=("x2", "mean"), z=("z", "mean"),
            w=("w", "sum"), dY=("dY", "first"), dD=("dD", "first"),
        )
        for c in ["Y", "D", "x1", "x2", "z", "w"]:
            out[c] = out[c].astype(np.float32).astype(np.float64)
        return out
    raise KeyError(name)


def run_python(case: dict):
    from twowayfeweights import twowayfeweights

    return twowayfeweights(
        load(case["data"]), case["Y"], case["G"], case["T"], case["D"], type=case["type"], D0=case["D0"],
        controls=case["controls"] or None, weights=case["weights"],
        other_treatments=case["other_treatments"] or None,
        test_random_weights=case["test_random_weights"] or None, summary_measures=True,
    )
