"""Parity with the Stata command (fixtures produced by tests/reference/build_reference.py).

Stata stores intermediate variables (P_gt, residuals, W, weights) as 4-byte floats, so agreement is
checked to ~1e-6 relative; counts must match exactly.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
from cases import CASES, FIXTURES, run_python

RTOL = 1e-6
# Regressions on Stata's float-stored W amplify its rounding slightly (observed max 1.3e-6).
RTOL_RW = 5e-6
STATA = FIXTURES / "stata"
RESULTS = pd.read_csv(STATA / "results.csv").set_index("name")
RW = pd.read_csv(STATA / "random_weights.csv")
OT = pd.read_csv(STATA / "other_treatments.csv")

params = [pytest.param(c, id=c["name"]) for c in CASES]


@pytest.fixture(scope="module")
def results():
    cache = {}

    def get(case):
        if case["name"] not in cache:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                cache[case["name"]] = run_python(case)
        return cache[case["name"]]

    return get


def _close(a, b, rtol=RTOL):
    if b is None or (isinstance(b, float) and np.isnan(b)):
        return a is None or np.isnan(a)
    return a is not None and abs(a - b) <= rtol * max(abs(b), 1e-12)


@pytest.mark.parametrize("case", params)
def test_scalars(case, results):
    r = results(case)
    s = RESULTS.loc[case["name"]]
    assert _close(r.beta, s["beta"], 1e-7)
    if case["beta_only"]:
        return
    assert r.nr_plus == s["nr_plus"]
    assert r.nr_minus == s["nr_minus"]
    assert r.tot_cells == s["tot_cells"]
    assert _close(r.sum_plus, s["sum_plus"])
    assert _close(r.sum_minus, s["sum_minus"]) if s["sum_minus"] != 0 else r.sum_minus == 0
    assert _close(r.sensibility, s["sensibility"])
    assert _close(r.sensibility2, s["sensibility2"])


@pytest.mark.parametrize("case", [p for p in params if p.values[0]["test_random_weights"]])
def test_random_weights(case, results):
    if case["beta_only"]:
        pytest.skip("Stata output depends on its sort order")
    r = results(case)
    ref = RW[RW["name"] == case["name"]].set_index("var")
    for v, row in ref.iterrows():
        got = r.mat.loc[v]
        assert _close(got["Coef"], row["Coef"], RTOL_RW)
        assert _close(got["SE"], row["SE"], RTOL_RW)
        assert _close(got["t-stat"], row["tstat"], RTOL_RW)
        assert _close(got["Correlation"], row["Correlation"], RTOL_RW)


@pytest.mark.parametrize("case", [p for p in params if p.values[0]["other_treatments"]])
def test_other_treatments(case, results):
    r = results(case)
    ref = OT[OT["name"] == case["name"]].sort_values("j")
    assert len(ref) == len(r.other_treatments)
    for o, (_, row) in zip(r.other_treatments, ref.iterrows()):
        assert o.nr_plus == row["nr_plus"]
        assert o.nr_minus == row["nr_minus"]
        assert o.tot_cells == row["tot_cells"]
        assert _close(o.sum_plus, row["sum_plus"])
        assert _close(o.sum_minus, row["sum_minus"])


@pytest.mark.parametrize("case", params)
def test_cell_weights(case, results):
    """Every (g,t) weight saved by Stata's path() option is reproduced."""
    if case["beta_only"]:
        pytest.skip("Stata output depends on its sort order")
    r = results(case)
    ref = pd.read_csv(STATA / "weights" / f"{case['name']}.csv.gz")
    keys = [k for k in ref.columns if k in ("Group_TWFE", "Time_TWFE", "Group", "Time")]
    m = ref.merge(r.weights, on=keys, how="left", suffixes=("_st", "_py"), validate="one_to_one")
    wcols = [c[:-3] for c in m.columns if c.endswith("_st")]
    for c in wcols:
        st, py = m[f"{c}_st"].to_numpy(), m[f"{c}_py"].to_numpy()
        assert (np.isnan(st) == np.isnan(py)).all()
        scale = np.nanmax(np.abs(st)) or 1.0
        np.testing.assert_allclose(py, st, rtol=0, atol=2e-6 * scale)
    # cells Stata leaves out of its file (it drops zero weights when some weights are negative) are zero
    extra = r.weights.merge(ref[keys], on=keys, how="left", indicator=True)
    assert (extra.loc[extra["_merge"] == "left_only", "weight"].fillna(0) == 0).all()
