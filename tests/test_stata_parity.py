"""Parity with the Stata command, on every case in ``cases.py``.

The Stata results are stored in ``fixtures/stata/results.json`` (one entry per case: the Stata command,
every number it returns and the text it prints) and ``fixtures/stata/weights/`` (the (g,t) weights saved
by ``path()``). They are produced by ``tests/reference/build_reference.py``.

Stata stores intermediate variables (P_gt, residuals, W, weights) as 4-byte floats, so agreement is
checked to ~1e-6 relative; counts and all printed text must match exactly.
"""

from __future__ import annotations

import json
import re
import warnings

import numpy as np
import pandas as pd
import pytest
from cases import CASES, FIXTURES, run_python

RTOL = 1e-6
# Regressions on Stata's float-stored W amplify its rounding slightly (observed max 1.3e-6).
RTOL_RW = 5e-6
STATA = FIXTURES / "stata"
REF = json.loads((STATA / "results.json").read_text(encoding="utf-8"))["cases"]

params = [pytest.param(c, id=c["name"]) for c in CASES]
exact = [pytest.param(c, id=c["name"]) for c in CASES if not c["beta_only"]]


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
    if b is None:
        return a is None or np.isnan(a)
    return a is not None and abs(a - b) <= rtol * max(abs(b), 1e-12)


def test_every_case_has_stata_results():
    assert sorted(REF) == sorted(c["name"] for c in CASES)


@pytest.mark.parametrize("case", params)
def test_scalars(case, results):
    r = results(case)
    s = REF[case["name"]]
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


@pytest.mark.parametrize("case", [p for p in exact if p.values[0]["test_random_weights"]])
def test_random_weights(case, results):
    r = results(case)
    ref = REF[case["name"]]["random_weights"]
    assert list(r.mat.index) == list(ref)
    for v, row in ref.items():
        got = r.mat.loc[v]
        assert _close(got["Coef"], row["Coef"], RTOL_RW)
        assert _close(got["SE"], row["SE"], RTOL_RW)
        assert _close(got["t-stat"], row["tstat"], RTOL_RW)
        assert _close(got["Correlation"], row["Correlation"], RTOL_RW)


@pytest.mark.parametrize("case", [p for p in params if p.values[0]["other_treatments"]])
def test_other_treatments(case, results):
    r = results(case)
    ref = REF[case["name"]]["other_treatments"]
    assert len(ref) == len(r.other_treatments)
    for o, row in zip(r.other_treatments, ref):
        assert o.nr_plus == row["nr_plus"]
        assert o.nr_minus == row["nr_minus"]
        assert o.tot_cells == row["tot_cells"]
        assert _close(o.sum_plus, row["sum_plus"])
        assert _close(o.sum_minus, row["sum_minus"])


_NUMBER = re.compile(r"(-?\d*\.\d+(?:e[-+]\d+)?|-?\d+)")


def _lines(text: str) -> list[str]:
    return [" ".join(line.split()) for line in text.splitlines() if line.strip()]


def _same_line(py: str, st: str) -> bool:
    """Same text; integers identical; decimals equal up to the digits Stata shows (+ its float rounding)."""
    a, b = _NUMBER.split(py), _NUMBER.split(st)
    if len(a) != len(b):
        return False
    for i, (x, y) in enumerate(zip(a, b)):
        if i % 2 == 0 or "." not in y:  # text, or an integer
            if x != y:
                return False
        else:
            mantissa, _, exponent = y.partition("e")
            last_digit = 10.0 ** (int(exponent or 0) - len(mantissa.split(".")[1]))
            if abs(float(x) - float(y)) > 0.5 * last_digit + RTOL_RW * abs(float(y)):
                return False
    return True


@pytest.mark.parametrize("case", exact)
def test_printed_output(case, results):
    """``print(result)`` shows the same lines as Stata (spacing and Stata's line wrapping ignored)."""
    py = _lines(results(case).summary())
    st = _lines(REF[case["name"]]["output"])
    assert len(py) == len(st), "\n".join(py)
    for p, s in zip(py, st):
        assert _same_line(p, s), f"\npython: {p}\nstata:  {s}"


@pytest.mark.parametrize("case", exact)
def test_cell_weights(case, results):
    """Every (g,t) weight saved by Stata's path() option is reproduced."""
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
