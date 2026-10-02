"""Estimation of the weights attached to two-way fixed effects regressions.

This is a line-by-line port of the Stata command ``twowayfeweights`` (de Chaisemartin and
D'Haultfoeuille). Comments of the form ``[ado: ...]`` point to the Stata statement being replicated.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from ._linalg import (
    aw_mean,
    aw_sd,
    cluster_wls,
    factorize_sorted,
    fwl_beta,
    group_mean_replace,
    partial_out,
)
from ._result import OtherTreatmentResult, TwoWayFEWeightsResult

TYPES = ("feTR", "feS", "fdTR", "fdS")
LIMIT_SENSITIVITY = 1e-10


def twowayfeweights(
    data,
    Y: str,
    G: str,
    T: str,
    D: str,
    type: str = "feTR",
    D0: str | None = None,
    summary_measures: bool = False,
    controls: str | Sequence[str] | None = None,
    weights: str | None = None,
    other_treatments: str | Sequence[str] | None = None,
    test_random_weights: str | Sequence[str] | None = None,
    path: str | Path | None = None,
) -> TwoWayFEWeightsResult:
    """Estimate the weights attached to a two-way fixed effects (or first-difference) regression.

    Parameters
    ----------
    data : pandas.DataFrame (or anything convertible to one, e.g. a polars DataFrame)
    Y : str
        Outcome variable (first-differenced outcome for ``fdTR``/``fdS``).
    G : str
        Group identifier.
    T : str
        Time period identifier (must be sortable).
    D : str
        Treatment variable (first-differenced treatment for ``fdTR``/``fdS``).
    type : {"feTR", "feS", "fdTR", "fdS"}
        * ``feTR`` - fixed-effects regression under common trends.
        * ``feS``  - fixed-effects regression, common trends and treatment effects stable over time.
        * ``fdTR`` - first-difference regression under common trends (requires ``D0``).
        * ``fdS``  - first-difference regression, common trends and stable treatment effects.
    D0 : str, optional
        Treatment level (non-differenced) variable. Required for ``type="fdTR"``.
    summary_measures : bool
        Display the sensitivity measures when printing the result (always computed).
    controls : str or list of str, optional
        Control variables included in the regression.
    weights : str, optional
        Variable of (analytic) weights.
    other_treatments : str or list of str, optional
        Other treatment variables included in the regression (``type="feTR"`` only).
    test_random_weights : str or list of str, optional
        Variables regressed on the weights, to assess whether the weights are correlated with
        variables that may be correlated with the treatment effect.
    path : str or Path, optional
        Save the (g,t)-level weights to this file (``.dta``, ``.parquet`` or ``.csv``).

    Returns
    -------
    TwoWayFEWeightsResult
    """
    if type not in TYPES:
        raise ValueError(f"type must be one of {TYPES}, got {type!r}")
    if type == "fdTR" and D0 is None:
        raise ValueError("The D0 argument must be provided when type='fdTR'.")
    controls = _as_list(controls)
    ots = _as_list(other_treatments)
    rw_vars = _as_list(test_random_weights)
    if ots and type != "feTR":
        raise ValueError("When other_treatments is specified, type must be 'feTR'.")
    if D0 is not None and type != "fdTR":
        warnings.warn(f"D0 is only used when type='fdTR'; it is ignored for type={type!r}.", stacklevel=2)
        D0 = None

    df = _to_pandas(data)
    needed = [Y, G, T, D, *controls, *ots, *rw_vars] + ([D0] if D0 else []) + ([weights] if weights else [])
    missing_cols = [c for c in dict.fromkeys(needed) if c not in df.columns]
    if missing_cols:
        raise KeyError(f"Columns not found in data: {missing_cols}")

    # ------------------------------------------------------------------ sample selection
    y = _num(df[Y], Y)
    d = _num(df[D], D)
    d0 = _num(df[D0], D0) if D0 else None
    X = _num_matrix(df, controls)
    OT = _num_matrix(df, ots)
    w = _num(df[weights], weights) if weights else np.ones(len(df))
    g_raw = df[G].to_numpy()
    t_raw = df[T].to_numpy()
    g_na = pd.isna(df[G]).to_numpy()
    t_na = pd.isna(df[T]).to_numpy()

    if type == "fdTR":
        # [ado: keep if (time!=. & outcome!=. & meantreat!=.) | treatment!=.]
        core_ok = ~(t_na | np.isnan(y) | np.isnan(d))
        keep = core_ok | ~np.isnan(d0)
        if X.shape[1]:
            keep &= ~(core_ok & np.isnan(X).any(axis=1))
    else:
        keep = ~(np.isnan(y) | g_na | t_na | np.isnan(d))
        if X.shape[1]:
            keep &= ~np.isnan(X).any(axis=1)
        if OT.shape[1]:
            keep &= ~np.isnan(OT).any(axis=1)

    idx = np.flatnonzero(keep)
    y, d, w, X, OT = y[idx], d[idx], w[idx], X[idx], OT[idx]
    d0 = d0[idx] if d0 is not None else None
    g_raw, t_raw, g_na, t_na = g_raw[idx], t_raw[idx], g_na[idx], t_na[idx]
    rw_raw = {v: df[v].to_numpy()[idx] for v in rw_vars}

    # ------------------------------------------------------------------ (g,t)-level variables
    cell0, n0 = _cells(g_raw, t_raw)
    m, varies = group_mean_replace(d[:, None], cell0, n0)
    if varies[0]:
        warnings.warn(
            "The treatment variable in the regression varies within some group * period cells. "
            "The command replaces the treatment by its average value in each group * period cell.",
            stacklevel=2,
        )
        d = m[:, 0]
    for name, M in (("control", X), ("other treatment", OT)):
        if M.shape[1] == 0:
            continue
        m, varies = group_mean_replace(M, cell0, n0)
        for j in np.flatnonzero(varies):
            var = (controls if name == "control" else ots)[j]
            warnings.warn(
                f"The {name} variable {var!r} varies within some group * period cells; "
                "it is replaced by its average value in each group * period cell.",
                stacklevel=2,
            )
            M[:, j] = m[:, j]
    # Variables used in test_random_weights see the (g,t)-averaged controls, as in Stata.
    for v in rw_vars:
        if v in controls:
            rw_raw[v] = X[:, controls.index(v)]
        elif v in ots:
            rw_raw[v] = OT[:, ots.index(v)]

    # [ado: keep if weight_XX!=.]
    ok_w = ~np.isnan(w)
    if not ok_w.all():
        sel = np.flatnonzero(ok_w)
        y, d, w, X, OT = y[sel], d[sel], w[sel], X[sel], OT[sel]
        d0 = d0[sel] if d0 is not None else None
        g_raw, t_raw, g_na, t_na = g_raw[sel], t_raw[sel], g_na[sel], t_na[sel]
        rw_raw = {k: v[sel] for k, v in rw_raw.items()}
    if (w < 0).any():
        raise ValueError("Negative weights encountered.")
    if len(y) == 0:
        raise ValueError("No observations left after removing missing values.")

    g_codes, n_g = factorize_sorted(g_raw)
    t_codes, n_t = factorize_sorted(t_raw)
    key = g_codes * (n_t + 1) + t_codes
    _, first, cell = np.unique(key, return_index=True, return_inverse=True)
    n_cells = len(first)
    cell_w = np.bincount(cell, weights=w, minlength=n_cells)
    nobs = w.sum()
    P_row = cell_w[cell] / nobs
    XO = np.hstack([X, OT]) if OT.shape[1] else X

    ctx = dict(
        y=y, d=d, d0=d0, w=w, X=XO, g=g_codes, t=t_codes, n_t=n_t, n_g=n_g,
        cell=cell, first=first, n_cells=n_cells, cell_w=cell_w, P_row=P_row, t_na=t_na,
    )
    if type == "feTR":
        cells = _fe_tr(ctx, OT)
    elif type == "fdTR":
        cells = _fd_tr(ctx)
    elif type == "feS":
        cells = _fe_s(ctx)
    else:
        cells = _fd_s(ctx)

    # ------------------------------------------------------------------ results
    weight = cells["weight"]
    weight[np.abs(weight) < LIMIT_SENSITIVITY] = 0.0
    summary = _summarize(weight)
    beta = cells["beta"]
    W, nat, gc, tc = cells["W"], cells["nat"], cells["g"], cells["t"]
    rows_c = cells["rows"]  # representative row for each cell

    if ots:
        tot_cells = cells["tot_cells"]
        sensibility = sensibility2 = None
    else:
        tot_cells = int(np.count_nonzero(nat != 0))  # NaN counts as nonzero, as in Stata
        sd = aw_sd(W, nat)
        sensibility = _div(abs(beta), sd)
        sensibility2 = _sensibility2(beta, W, weight, nat, gc, tc, n_g, n_t) if summary["sum_minus"] < 0 else None

    mat = None
    if rw_vars:
        cl = np.where(gc < n_g, gc, np.nan).astype(float)
        res = []
        for v in rw_vars:
            vals = _num(pd.Series(rw_raw[v][rows_c]), v)
            b, se, r2 = cluster_wls(vals, W, nat, cl)
            corr = (1.0 if b >= 0 else -1.0) * np.sqrt(r2) if np.isfinite(r2) else np.nan
            res.append((b, se, _div(b, se), corr))
        mat = pd.DataFrame(res, index=rw_vars, columns=["Coef", "SE", "t-stat", "Correlation"])

    wdf = pd.DataFrame(
        {
            "Group_TWFE" if not ots else "Group": g_raw[rows_c],
            "Time_TWFE" if not ots else "Time": t_raw[rows_c],
            "weight": weight,
        }
    )
    other_results = []
    for j, name in enumerate(ots):
        wo = cells["weight_others"][:, j]
        wdf[f"weight_others{j + 1}"] = wo
        s = _summarize(wo)
        other_results.append(OtherTreatmentResult(name=name, tot_cells=cells["tot_cells_ot"][j], **s))

    result = TwoWayFEWeightsResult(
        type=type,
        treatment=D,
        beta=beta,
        nr_plus=summary["nr_plus"],
        nr_minus=summary["nr_minus"],
        nr_weights=summary["nr_plus"] + summary["nr_minus"],
        sum_plus=summary["sum_plus"],
        sum_minus=summary["sum_minus"],
        tot_cells=tot_cells,
        sensibility=sensibility,
        sensibility2=sensibility2,
        mat=mat,
        other_treatments=other_results,
        weights=wdf,
        summary_measures=summary_measures,
    )
    if path is not None:
        _save(wdf, Path(path))
    return result


# ---------------------------------------------------------------------- per-type computations


def _collapse(ctx, *arrays):
    """Value of each row-level array at the first observation of every (g,t) cell."""
    f = ctx["first"]
    return [a[f] for a in arrays]


def _fe_tr(ctx, OT):
    y, d, w, P_row = ctx["y"], ctx["d"], ctx["w"], ctx["P_row"]
    mean_D = aw_mean(d, w)  # [ado: sum meantreat [aw=weight_XX]]
    nat = P_row * d / mean_D
    fit = np.ones(len(y), dtype=bool)
    e = partial_out(d, ctx["X"], ctx["g"], ctx["t"], w, ctx["cell"], fit)  # [ado: areg D i.T X, absorb(G)]
    beta = fwl_beta(e, y, d, w)
    denom = aw_mean(e * d, w)
    W = e * mean_D / denom
    weight = W * nat
    out = {"beta": beta}
    if OT.shape[1]:
        wo = W[:, None] * P_row[:, None] * OT / mean_D
        wo[np.abs(wo) < LIMIT_SENSITIVITY] = 0.0
        # Stata counts these before collapsing to (g,t) cells, i.e. at the observation level.
        out["tot_cells"] = int(np.count_nonzero(nat != 0))
        out["tot_cells_ot"] = [int(np.count_nonzero(OT[:, j] != 0)) for j in range(OT.shape[1])]
        (out["weight_others"],) = _collapse(ctx, wo)
    W, weight, nat, g, t = _collapse(ctx, W, weight, nat, ctx["g"], ctx["t"])
    out.update(W=W, weight=weight, nat=nat, g=g, t=t, rows=ctx["first"])
    return out


def _fd_tr(ctx):
    y, d, d0, w, P_row, X = ctx["y"], ctx["d"], ctx["d0"], ctx["w"], ctx["P_row"], ctx["X"]
    mean_D = aw_mean(d0, w)
    nat = P_row * d0 / mean_D
    fit = ~(np.isnan(d) | ctx["t_na"])
    if X.shape[1]:
        fit &= ~np.isnan(X).any(axis=1)
    e = partial_out(d, X, ctx["t"], None, w, ctx["cell"], fit)  # [ado: reg D i.T X]
    fit_b = fit & ~np.isnan(y)
    e_b = e if (fit_b == fit).all() else partial_out(d, X, ctx["t"], None, w, ctx["cell"], fit_b)
    beta = fwl_beta(e_b, y, d, w)  # [ado: areg Y D X, absorb(T)]
    e = np.where(np.isnan(e), 0.0, e)  # [ado: replace eps_2=0 if eps_2==.]

    e, P, d0c, natc, g, t = _collapse(ctx, e, P_row, d0, nat, ctx["g"], ctx["t"])
    n_t = ctx["n_t"]
    wt = e.copy()
    if len(e) > 1:
        same = (g[1:] == g[:-1]) & (t[:-1] < n_t) & (t[1:] == t[:-1] + 1)
        with np.errstate(divide="ignore", invalid="ignore"):
            cand = e[:-1] - e[1:] * P[1:] / P[:-1]
        ok = same & np.isfinite(cand)
        wt[:-1][ok] = cand[ok]
    denom = aw_mean(wt * d0c, P)
    W = wt * mean_D / denom
    return dict(beta=beta, W=W, weight=W * natc, nat=natc, g=g, t=t, rows=ctx["first"])


def _fe_s(ctx):
    y, d, w, P_row = ctx["y"], ctx["d"], ctx["w"], ctx["P_row"]
    fit = np.ones(len(y), dtype=bool)
    e = partial_out(d, ctx["X"], ctx["g"], ctx["t"], w, ctx["cell"], fit)
    beta = fwl_beta(e, y, d, w)

    n_cells, cell = ctx["n_cells"], ctx["cell"]
    cw = ctx["cell_w"]
    cwe = np.bincount(cell, weights=np.nan_to_num(e * w), minlength=n_cells)
    dc, P, g, t = _collapse(ctx, d, P_row, ctx["g"], ctx["t"])
    # E(eps_1 | g, t' >= t), weighted: reverse cumulative sums within each group.
    with np.errstate(divide="ignore", invalid="ignore"):
        E = _rev_cumsum_by_group(cwe, g) / _rev_cumsum_by_group(cw, g)

    delta = np.full(n_cells, np.nan)
    if n_cells > 1:
        prev_ok = (g[1:] == g[:-1]) & (t[1:] == t[:-1] + 1)
        delta[1:][prev_ok] = dc[1:][prev_ok] - dc[:-1][prev_ok]
    k = ~np.isnan(delta)  # [ado: drop if Delta_D==.]
    delta, P, E, g, t = delta[k], P[k], E[k], g[k], t[k]
    s = np.sign(delta)
    nat = P * np.abs(delta)
    nat = nat / nat.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        om = s * E / P
    denom = aw_mean(om, nat)
    W = om / denom
    return dict(beta=beta, W=W, weight=W * nat, nat=nat, g=g, t=t, rows=ctx["first"][k])


def _fd_s(ctx):
    y, d, w, P_row = ctx["y"], ctx["d"], ctx["w"], ctx["P_row"]
    fit = np.ones(len(y), dtype=bool)
    e = partial_out(d, ctx["X"], ctx["t"], None, w, ctx["cell"], fit)
    beta = fwl_beta(e, y, d, w)
    e, dc, P, g, t = _collapse(ctx, e, d, P_row, ctx["g"], ctx["t"])
    s = np.sign(dc)
    nat = P * np.abs(dc)
    nat = nat / np.nansum(nat)
    W = s * e
    W = W / aw_mean(W, nat)
    return dict(beta=beta, W=W, weight=W * nat, nat=nat, g=g, t=t, rows=ctx["first"])


# ---------------------------------------------------------------------- helpers


def _sensibility2(beta, W, weight, nat, g, t, n_g, n_t):
    """Second sensitivity measure (Stata's cumulative-sum algorithm, including its sort order)."""
    k = weight != 0  # [ado: keep if weight!=0] (missing weights are kept)
    W, weight, nat, g, t = W[k], weight[k], nat[k], g[k], t[k]
    n = len(W)
    if n == 0:
        return None
    Wk = np.where(np.isnan(W), np.inf, W)
    Wmiss = np.isnan(W).astype(np.int8)
    # gsort W -G -T  (missing values last for both ascending and descending keys)
    g_desc = np.where(g >= n_g, np.iinfo(np.int64).max, -g)
    t_desc = np.where(t >= n_t, np.iinfo(np.int64).max, -t)
    asc = np.lexsort((t_desc, g_desc, Wk, Wmiss))
    P_k = np.empty(n)
    S_k = np.empty(n)
    T_k = np.empty(n)
    P_k[asc] = np.cumsum(np.nan_to_num(nat[asc]))
    S_k[asc] = np.cumsum(np.nan_to_num(weight[asc]))
    T_k[asc] = np.cumsum(np.nan_to_num((nat * W**2)[asc]))
    # gsort -W G T
    W_desc = np.where(np.isnan(W), np.inf, -W)
    desc = np.lexsort((t, g, W_desc, Wmiss))
    W, P_k, S_k, T_k = W[desc], P_k[desc], S_k[desc], T_k[desc]

    denom = 1.0 - P_k
    with np.errstate(divide="ignore", invalid="ignore"):
        q = np.where(denom == 0, np.nan, S_k / denom)
        inner = np.where(denom == 0, np.nan, T_k + S_k**2 / denom)
        sens = np.where(inner > 0, abs(beta) / np.sqrt(np.where(inner > 0, inner, 1.0)), np.nan)
    # Stata: nonmissing < missing is true; missing < anything is false.
    ind = ~np.isnan(W) & (np.isnan(q) | (W < -q))
    ind[0] = False
    ind = np.maximum.accumulate(ind)
    tot = int(ind.sum())
    if tot == 0:
        return None
    val = sens[n - tot]
    return None if np.isnan(val) else float(val)


def _summarize(weight):
    pos = weight > 0
    neg = weight < 0
    return dict(
        nr_plus=int(pos.sum()),
        nr_minus=int(neg.sum()),
        sum_plus=float(weight[pos].sum()),
        sum_minus=float(weight[neg].sum()),
    )


def _rev_cumsum_by_group(x, g):
    """Sum of ``x`` over rows at or after each row within its (contiguous, sorted) group."""
    cs = np.cumsum(x[::-1])[::-1]  # suffix sums over the whole array
    ends = np.flatnonzero(np.r_[g[1:] != g[:-1], True])  # last row of each group
    nxt = np.r_[cs[ends[:-1] + 1], 0.0] if len(ends) else np.array([])
    grp = np.cumsum(np.r_[0, g[1:] != g[:-1]])
    return cs - nxt[grp]


def _cells(g_raw, t_raw):
    gc, ng = factorize_sorted(g_raw)
    tc, nt = factorize_sorted(t_raw)
    _, cell = np.unique(gc * (nt + 1) + tc, return_inverse=True)
    return cell, int(cell.max()) + 1 if len(cell) else 0


def _div(a, b):
    if b is None or not np.isfinite(b) or b == 0 or not np.isfinite(a):
        return np.nan
    return float(a / b)


def _as_list(x) -> list[str]:
    if x is None:
        return []
    if isinstance(x, str):
        return [x]
    return list(x)


def _to_pandas(data) -> pd.DataFrame:
    if isinstance(data, pd.DataFrame):
        return data
    if hasattr(data, "to_pandas"):
        try:
            return data.to_pandas()
        except ImportError:  # e.g. polars without pyarrow installed
            if hasattr(data, "to_dict"):
                return pd.DataFrame(data.to_dict(as_series=False))
            raise
    return pd.DataFrame(data)


def _num(s: pd.Series, name: str) -> np.ndarray:
    if s.dtype == object or isinstance(s.dtype, pd.CategoricalDtype) or pd.api.types.is_string_dtype(s.dtype):
        try:
            return pd.to_numeric(s).to_numpy(dtype=float, na_value=np.nan)
        except (ValueError, TypeError) as err:
            raise TypeError(f"Variable {name!r} must be numeric.") from err
    return s.to_numpy(dtype=float, na_value=np.nan)


def _num_matrix(df: pd.DataFrame, cols: Iterable[str]) -> np.ndarray:
    cols = list(cols)
    if not cols:
        return np.zeros((len(df), 0))
    sub = df[cols]
    if all(pd.api.types.is_numeric_dtype(t) or pd.api.types.is_bool_dtype(t) for t in sub.dtypes):
        return sub.to_numpy(dtype=np.float64, na_value=np.nan)
    return np.column_stack([_num(df[c], c) for c in cols])


def _save(wdf: pd.DataFrame, path: Path) -> None:
    suffix = path.suffix.lower()
    if suffix == ".dta":
        wdf.to_stata(path, write_index=False)
    elif suffix == ".parquet":
        wdf.to_parquet(path, index=False)
    else:
        wdf.to_csv(path, index=False)
