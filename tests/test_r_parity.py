"""Parity with the R package TwoWayFEWeights (CRAN 2.1.0).

R is compared only where it agrees with Stata: R errors on missing weights and on micro data with
missing values, and deviates from Stata for fdTR/feTR micro data with missing outcomes (Python follows
Stata there, see test_stata_parity.py). Where R agrees with Stata, Python and R (both double
precision) must agree to 1e-8.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
from cases import CASES, FIXTURES, run_python

R = pd.read_csv(FIXTURES / "r" / "results.csv").set_index("name")
RW = pd.read_csv(FIXTURES / "r" / "random_weights.csv")
STATA = pd.read_csv(FIXTURES / "stata" / "results.csv").set_index("name")
COLS = ["beta", "nr_plus", "nr_minus", "sum_plus", "sum_minus", "tot_cells", "sensibility", "sensibility2"]


def _rel(a, b):
    a = np.nan if a is None else float(a)
    b = float(b)
    if np.isnan(a) and np.isnan(b):
        return 0.0
    if np.isnan(a) or np.isnan(b):
        return np.inf
    return abs(a - b) / max(abs(b), 1e-12)


def _r_matches_stata(name):
    r, s = R.loc[name], STATA.loc[name]
    return max(_rel(r[k], s[k]) for k in COLS) < 1e-6


CMP = [c for c in CASES if c["name"] in R.index and not R.loc[c["name"], "error"] and _r_matches_stata(c["name"])]


def test_enough_r_cases():
    assert len(CMP) >= 15


@pytest.mark.parametrize("case", CMP, ids=[c["name"] for c in CMP])
def test_against_r(case):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = run_python(case)
    d = res.to_dict()
    r = R.loc[case["name"]]
    for k in COLS:
        assert _rel(d[k], r[k]) < 1e-8, k
    ref = RW[RW["name"] == case["name"]]
    for _, row in ref.iterrows():
        got = res.mat.loc[row["var"]]
        for k, kk in [("Coef", "Coef"), ("SE", "SE"), ("t-stat", "tstat"), ("Correlation", "Correlation")]:
            assert _rel(got[k], row[kk]) < 1e-8, (row["var"], k)
