"""Behavioural tests: inputs, errors, invariances, printing, output files."""

from __future__ import annotations

import io
import sys
import warnings

import numpy as np
import pandas as pd
import pytest
from cases import load

import twowayfeweights as tw
from twowayfeweights import print_twowayfeweights, twowayfeweights


@pytest.fixture(scope="module")
def wagepan():
    return tw.load_wagepan()


def test_wagepan_r_package_values(wagepan):
    """Known values from the R package's tinytest suite (tests/test_wagepan.R)."""
    r = twowayfeweights(wagepan, "lwage", "nr", "year", "union", type="feTR", summary_measures=True,
                        test_random_weights="educ")
    assert (r.nr_plus, r.nr_minus, r.nr_weights, r.tot_cells) == (820, 147, 967, 1016)
    np.testing.assert_allclose([r.beta, r.sum_plus, r.sum_minus, r.sensibility, r.sensibility2],
                               [0.10662746641838491, 1.010529, -0.01052899, 0.0968691, 3.175859], rtol=1e-4)
    np.testing.assert_allclose(r.mat.loc["educ"].to_numpy(),
                               [-0.13445527172928173, 0.07136021078287243, -1.88417705405031, -0.11825873818614765],
                               rtol=1e-5)
    r = twowayfeweights(wagepan, "diff_lwage", "nr", "year", "diff_union", type="fdTR", D0="union",
                        summary_measures=True, test_random_weights="educ")
    assert (r.nr_plus, r.nr_minus, r.tot_cells) == (611, 405, 1016)
    np.testing.assert_allclose([r.beta, r.sum_plus, r.sensibility, r.sensibility2],
                               [0.06009597, 1.047636, 0.03209515, 0.5799135], rtol=1e-5)


def test_did_book_chapter5():
    """Numbers printed in https://anzonyquispe.github.io/did_book/chapters/ch05.html (Stata output)."""
    df = load("gentzkow")
    styr = [c for c in df.columns if c.startswith("styr") and c != "styr"]
    r = twowayfeweights(df, "prestout", "cnty90", "year", "numdailies", type="feTR")
    assert (r.nr_weights, r.nr_plus, r.nr_minus) == (10378, 6180, 4198)
    assert (round(r.sum_plus, 4), round(r.sum_minus, 4), round(r.beta, 4)) == (1.4740, -0.4740, 0.0029)
    r = twowayfeweights(df, "prestout", "cnty90", "year", "numdailies", type="feTR", controls=styr)
    assert (r.nr_weights, r.nr_plus, r.nr_minus, r.tot_cells) == (10342, 6195, 4147, 10378)
    assert (round(r.sum_plus, 4), round(r.beta, 4)) == (1.5331, -0.0012)
    r = twowayfeweights(df, "changeprestout", "cnty90", "year", "changedailies", type="fdTR", D0="numdailies",
                        controls=styr, test_random_weights="year")
    assert (r.nr_weights, r.nr_plus, r.nr_minus, r.tot_cells) == (9876, 5371, 4505, 10378)
    assert (round(r.sum_plus, 4), round(r.sum_minus, 4), round(r.beta, 4)) == (2.4271, -1.4271, 0.0026)
    np.testing.assert_allclose(r.mat.loc["year"].to_numpy(), [-0.1674271, 0.05101171, -3.2821307, -0.0631614],
                               rtol=1e-6)


@pytest.mark.parametrize("type_", ["feTR", "feS", "fdTR", "fdS"])
def test_row_order_invariance(type_):
    df = load("simcell")
    kw = dict(type=type_, D0="D" if type_ == "fdTR" else None, controls=["x1"], weights="w", test_random_weights="z")
    y, d = ("Y", "D") if type_.startswith("fe") else ("dY", "dD")
    a = twowayfeweights(df, y, "g", "t", d, **kw)
    b = twowayfeweights(df.sample(frac=1, random_state=3), y, "g", "t", d, **kw)
    assert a.to_dict() == pytest.approx(b.to_dict(), rel=1e-10, nan_ok=True)
    pd.testing.assert_frame_equal(a.mat, b.mat, rtol=1e-9)


def test_string_and_categorical_ids_match_numeric():
    df = load("sim").copy()
    base = twowayfeweights(df, "Y", "g", "t", "D", controls=["x1"], weights="w", test_random_weights="z")
    df["gs"] = "grp_" + df["g"].astype(str).str.zfill(3)
    df["gc"] = df["gs"].astype("category")
    for gcol in ("gs", "gc"):
        r = twowayfeweights(df, "Y", gcol, "t", "D", controls=["x1"], weights="w", test_random_weights="z")
        assert r.to_dict() == pytest.approx(base.to_dict(), rel=1e-10, nan_ok=True)


def test_feS_micro_equals_collapsed():
    """feS on individual data equals feS on (g,t) cell means with summed weights (paper's formula)."""
    micro = load("sim").dropna(subset=["Y", "x1", "w"])
    cells = micro.groupby(["g", "t"], as_index=False).agg(Y=("Y", "mean"), D=("D", "mean"), x1=("x1", "mean"),
                                                          w=("w", "sum"))
    micro = micro.copy()
    wy = micro["Y"] * micro["w"]
    cells["Y"] = (wy.groupby([micro["g"], micro["t"]]).sum() / micro.groupby(["g", "t"])["w"].sum()).to_numpy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = twowayfeweights(micro, "Y", "g", "t", "D", type="feS", controls=["x1"], weights="w")
    b = twowayfeweights(cells, "Y", "g", "t", "D", type="feS", controls=["x1"], weights="w")
    assert a.to_dict() == pytest.approx(b.to_dict(), rel=1e-8, nan_ok=True)


def test_polars_input():
    pl = pytest.importorskip("polars")
    df = tw.load_wagepan()
    a = twowayfeweights(df, "lwage", "nr", "year", "union")
    b = twowayfeweights(pl.from_pandas(df), "lwage", "nr", "year", "union")
    assert a.to_dict() == b.to_dict()


def test_normalization_warning():
    with pytest.warns(UserWarning, match="varies within some group"):
        twowayfeweights(load("sim"), "Y", "g", "t", "D")


@pytest.mark.parametrize(
    "kwargs, err",
    [
        (dict(type="xx"), ValueError),
        (dict(type="fdTR"), ValueError),
        (dict(type="fdS", other_treatments="married"), ValueError),
        (dict(controls="nope"), KeyError),
    ],
)
def test_errors(wagepan, kwargs, err):
    with pytest.raises(err):
        twowayfeweights(wagepan, "lwage", "nr", "year", "union", **kwargs)


def test_negative_weights_rejected(wagepan):
    df = wagepan.assign(wneg=-1.0)
    with pytest.raises(ValueError, match="Negative weights"):
        twowayfeweights(df, "lwage", "nr", "year", "union", weights="wneg")


@pytest.mark.parametrize("ext", [".csv", ".dta"])
def test_path_output(tmp_path, wagepan, ext):
    p = tmp_path / f"w{ext}"
    r = twowayfeweights(wagepan, "lwage", "nr", "year", "union", path=p)
    back = pd.read_csv(p) if ext == ".csv" else pd.read_stata(p)
    assert list(back.columns) == ["Group_TWFE", "Time_TWFE", "weight"]
    np.testing.assert_allclose(back["weight"].sum(), r.sum_plus + r.sum_minus)


def test_printing_and_legacy_access(wagepan, capsys):
    r = twowayfeweights(wagepan, "lwage", "nr", "year", "union", summary_measures=True, test_random_weights="educ",
                        other_treatments=None)
    text = r.summary()
    assert "estimates a weighted sum of 967 ATTs" in text
    assert "1016 (g,t) cells receive the treatment" in text
    assert "Summary Measures:" in text and "3.1759" in text
    print_twowayfeweights(r, D_name="union", type="feTR")
    assert "Positive weights        820         1.0105" in capsys.readouterr().out
    assert r["n_pos"] == 820 and r["n_atts"] == 967 and r["n_zero"] == 49
    assert r["test_random_weights"]["educ"]["coef"] == pytest.approx(-0.1344553, rel=1e-6)
    assert list(r.M.index) == ["Pos_Weights", "Neg_Weights", "Tot"]


def test_printing_other_treatments(wagepan):
    r = twowayfeweights(wagepan, "lwage", "nr", "year", "union", other_treatments=["married", "south"])
    text = r.summary()
    assert "estimates the sum of several terms" in text
    assert "Other treat.: married" in text and "Other treat.: south" in text
    assert r["other_treatments_results"]["married"]["nr_plus"] == r.other_treatments[0].nr_plus


def test_print_on_cp1252_console(wagepan, monkeypatch):
    r = twowayfeweights(wagepan, "lwage", "nr", "year", "union", summary_measures=True)
    buf = io.TextIOWrapper(io.BytesIO(), encoding="cp1252")
    monkeypatch.setattr(sys, "stdout", buf)
    print(r)  # must not raise UnicodeEncodeError
    buf.flush()
    assert b"Sum weights" in buf.buffer.getvalue()
