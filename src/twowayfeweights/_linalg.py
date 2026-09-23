"""Numerical kernels: grouping helpers, fixed-effect residualization, clustered WLS.

Everything here works on plain numpy arrays. Missing values are ``NaN``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

# Regressors whose weighted sum of squares (after absorbing fixed effects) is below this fraction
# of the largest one are treated as collinear with the fixed effects.
_COLLINEAR_TOL = 1e-10
# Eigenvalues of the unit-diagonal Gram matrix below this fraction of the largest are dropped
# (exact collinearity among regressors), i.e. singular values below ~3e-6 of the largest.
_EIG_TOL = 1e-11


def factorize_sorted(values) -> tuple[np.ndarray, int]:
    """Dense ranks of ``values`` (0-based, in sorted order); missing values get code ``n_levels``.

    Placing missing values after every observed level mirrors Stata, where missing sorts last.
    """
    codes, uniques = pd.factorize(pd.Series(values), sort=True, use_na_sentinel=True)
    codes = codes.astype(np.int64)
    n_levels = len(uniques)
    codes[codes < 0] = n_levels
    return codes, n_levels


def indicator(codes: np.ndarray, n: int) -> sp.csr_matrix:
    """Sparse ``len(codes) x n`` 0/1 matrix mapping rows to groups."""
    m = len(codes)
    return sp.csr_matrix((np.ones(m), (np.arange(m), codes)), shape=(m, n))


def group_mean_replace(M: np.ndarray, cell: np.ndarray, n_cells: int) -> tuple[np.ndarray, np.ndarray]:
    """Cell means of each column of ``M`` (ignoring NaN) and a per-column flag for within-cell variation.

    Returns ``(means_by_row, varies)`` where ``means_by_row`` has the shape of ``M`` (NaN for cells
    without any non-missing value) and ``varies[j]`` is True when column ``j`` takes more than one
    non-missing value within at least one cell (Stata's ``sd > 0`` check).
    """
    k = M.shape[1]
    if n_cells == len(cell):  # one observation per cell: nothing can vary
        return M, np.zeros(k, dtype=bool)
    # Compare every row with the first row of its cell; columns with missing values need the
    # slower NaN-aware min/max check (the first row of a cell may be the missing one).
    first = np.full(n_cells, -1, dtype=np.int64)
    first[cell[::-1]] = np.arange(len(cell))[::-1]
    ref = M[first[cell]]
    has_nan = np.isnan(M).any(axis=0)
    varies = np.zeros(k, dtype=bool)
    clean = ~has_nan
    if clean.any():
        varies[clean] = (M[:, clean] != ref[:, clean]).any(axis=0)
    if has_nan.any():
        order = np.argsort(cell, kind="stable")
        sc = cell[order]
        starts = np.flatnonzero(np.r_[True, sc[1:] != sc[:-1]])
        Ms = M[order][:, has_nan]
        with np.errstate(invalid="ignore"):
            hi = np.fmax.reduceat(Ms, starts, axis=0)
            lo = np.fmin.reduceat(Ms, starts, axis=0)
        varies[has_nan] = np.any(hi > lo, axis=0)
    if not varies.any():
        return M, varies

    Mv = M[:, varies]
    S = indicator(cell, n_cells)
    ok = ~np.isnan(Mv)
    sums = np.asarray(S.T @ np.where(ok, Mv, 0.0))
    counts = np.asarray(S.T @ ok.astype(np.float64))
    with np.errstate(invalid="ignore", divide="ignore"):
        means = sums / counts
    out = M.copy()
    out[:, varies] = means[cell]
    return out, varies


def partial_out(
    v: np.ndarray,
    X: np.ndarray | None,
    absorb_codes: np.ndarray,
    dummy_codes: np.ndarray | None,
    w: np.ndarray,
    cell: np.ndarray,
    fit: np.ndarray,
) -> np.ndarray:
    """Residuals of a weighted regression of ``v`` on ``X``, dummies and an absorbed fixed effect.

    Mirrors ``areg v i.dummy X [aw=w], absorb(absorb)`` followed by ``predict, residuals`` (with
    ``dummy_codes=None`` it is ``areg v X, absorb(absorb)``, i.e. ``reg v i.absorb X``).
    ``v``, ``X`` and the codes must be constant within each (G,T) ``cell`` among the ``fit`` rows,
    which lets the regression run on collapsed cells with summed weights (identical coefficients,
    far fewer rows). Residuals are returned for every ``fit`` row, including zero-weight rows
    (as Stata's ``predict``); NaN elsewhere.

    The dummies are never materialized: after the within-transformation on ``absorb`` their Gram
    blocks follow from the sparse (absorb x dummy) weight table. The normal equations are solved with
    a scaled pseudo-inverse (dropping collinear directions) and refined twice, which recovers the
    accuracy of an orthogonal (QR) solver.
    """
    n = len(v)
    out = np.full(n, np.nan)
    rows = np.flatnonzero(fit)
    if rows.size == 0:
        return out
    _, first, inv = np.unique(cell[rows], return_index=True, return_inverse=True)
    Wc = np.bincount(inv, weights=w[rows])
    r0 = rows[first]
    vc = v[r0]
    Xc = X[r0] if X is not None and X.shape[1] else np.zeros((len(r0), 0))
    ac_all = absorb_codes[r0]
    est = np.flatnonzero(Wc > 0)
    if est.size == 0:
        return out

    We = Wc[est]
    a_lv, a = np.unique(ac_all[est], return_inverse=True)
    na = a_lv.size
    Wa = np.bincount(a, weights=We, minlength=na)

    def demean(M):  # weighted within-transformation over the absorbed groups
        if M.ndim == 1:
            return M - (np.bincount(a, weights=We * M, minlength=na) / Wa)[a]
        A = indicator(a, na)
        return M - (np.asarray(A.T @ (We[:, None] * M)) / Wa[:, None])[a]

    vt = demean(vc[est])
    Xt = demean(Xc[est]) if Xc.shape[1] else Xc[est]
    k = Xt.shape[1]

    # Dummies: base level = first level present in the estimation sample (dropped).
    nb = 0
    if dummy_codes is not None:
        b_lv, b = np.unique(dummy_codes[r0][est], return_inverse=True)
        nb = b_lv.size - 1
    p = nb + k
    coef_b = np.zeros(nb)
    coef_x = np.zeros(k)

    if p > 0:
        WX = We[:, None] * Xt
        if nb:
            Wb = np.bincount(b, weights=We, minlength=nb + 1)
            C = sp.csr_matrix((We, (a, b)), shape=(na, nb + 1))
            Gbb = np.diag(Wb) - np.asarray((C.T @ sp.diags(1.0 / Wa) @ C).todense())
            Gbb = Gbb[1:, 1:]
            B = indicator(b, nb + 1)
            Gbx = np.asarray(B.T @ WX)[1:] if k else np.zeros((nb, 0))
        Gxx = Xt.T @ WX if k else np.zeros((0, 0))
        if nb:
            G = np.block([[Gbb, Gbx], [Gbx.T, Gxx]])
        else:
            G = Gxx

        def cross(r):  # Z' W r  for the (dummy-demeaned, X) design
            parts = []
            if nb:
                parts.append(np.bincount(b, weights=We * r, minlength=nb + 1)[1:])
            if k:
                parts.append(WX.T @ r)
            return np.concatenate(parts)

        def fitted(cb, cx):  # Z @ coef on the estimation cells
            f = np.zeros(len(We))
            if nb:
                f += demean(np.r_[0.0, cb][b])
            if k:
                f += Xt @ cx
            return f

        d = np.diag(G).copy()
        ok = d > _COLLINEAR_TOL * max(d.max(), 1.0)
        s = np.zeros(p)
        s[ok] = 1.0 / np.sqrt(d[ok])
        Gs = G * s[:, None] * s[None, :]
        lam, V = np.linalg.eigh(Gs)
        keep = lam > _EIG_TOL * max(lam.max(), 1.0)
        Vk = V[:, keep] / np.sqrt(lam[keep])

        def solve(g):
            return s * (Vk @ (Vk.T @ (s * g)))

        coef = solve(cross(vt))
        for _ in range(2):  # iterative refinement
            coef += solve(cross(vt - fitted(coef[:nb], coef[nb:])))
        coef_b, coef_x = coef[:nb], coef[nb:]

    # Linear prediction for every cell in the fit sample, then the absorbed effects.
    xb = Xc @ coef_x if k else np.zeros(len(vc))
    if nb:
        pos = np.full(int(dummy_codes.max()) + 1, -1, dtype=np.int64)
        pos[b_lv] = np.arange(nb + 1)
        bi = pos[dummy_codes[r0]]
        cb_full = np.r_[0.0, coef_b]
        xb = xb + np.where(bi >= 0, cb_full[np.maximum(bi, 0)], 0.0)  # unseen levels act as the base
    fe = np.bincount(a, weights=We * (vc[est] - xb[est]), minlength=na) / Wa
    fe_map = np.full(int(absorb_codes.max()) + 1, np.nan)
    fe_map[a_lv] = fe
    resid_c = vc - xb - fe_map[ac_all]
    out[rows] = resid_c[inv]
    return out


def fwl_beta(e: np.ndarray, y: np.ndarray, d: np.ndarray, w: np.ndarray) -> float:
    """Coefficient on ``d`` via Frisch-Waugh-Lovell, given ``e`` = residualized ``d`` on the same sample."""
    ok = ~(np.isnan(e) | np.isnan(y) | np.isnan(d)) & (w > 0)
    we = w[ok] * e[ok]
    den = np.dot(we, d[ok])
    return float(np.dot(we, y[ok]) / den) if den != 0 else np.nan


def aw_mean(x: np.ndarray, w: np.ndarray) -> float:
    """``summarize x [aweight=w]``: r(mean). Missing and zero-weight observations are excluded."""
    ok = ~(np.isnan(x) | np.isnan(w)) & (w != 0)
    sw = w[ok].sum()
    return float(np.dot(w[ok], x[ok]) / sw) if ok.any() and sw != 0 else np.nan


def aw_sd(x: np.ndarray, w: np.ndarray) -> float:
    """``summarize x [aweight=w]``: r(sd) (weights rescaled to sum to N, divisor N - 1)."""
    ok = ~(np.isnan(x) | np.isnan(w)) & (w != 0)
    n = int(ok.sum())
    if n < 2:
        return np.nan
    xs, ws = x[ok], w[ok]
    sw = ws.sum()
    m = np.dot(ws, xs) / sw
    var = np.dot(ws, (xs - m) ** 2) / sw * n / (n - 1)
    return float(np.sqrt(var)) if var >= 0 else np.nan


def cluster_wls(y: np.ndarray, x: np.ndarray, w: np.ndarray, cluster: np.ndarray) -> tuple[float, float, float]:
    """``reg y x [pweight=w], cluster(cluster)``: returns (b_x, se_x, r2).

    Uses Stata's small-sample factor (N-1)/(N-k) * C/(C-1) with k = 2 (slope + constant).
    """
    cluster = np.asarray(cluster, dtype=float)
    ok = ~(np.isnan(y) | np.isnan(x) | np.isnan(w) | np.isnan(cluster)) & (w != 0)
    y, x, w, cl = y[ok], x[ok], w[ok], cluster[ok]
    n = y.size
    if n < 3:
        return np.nan, np.nan, np.nan
    Xm = np.column_stack([x, np.ones(n)])
    XtWX = Xm.T @ (w[:, None] * Xm)
    try:
        bread = np.linalg.inv(XtWX)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan
    b = bread @ (Xm.T @ (w * y))
    e = y - Xm @ b
    _, cl_codes = np.unique(cl, return_inverse=True)
    n_cl = int(cl_codes.max()) + 1
    scores = Xm * (w * e)[:, None]
    U = np.zeros((n_cl, 2))
    np.add.at(U, cl_codes, scores)
    meat = U.T @ U
    V = bread @ meat @ bread
    if n_cl > 1:
        V *= (n - 1) / (n - 2) * n_cl / (n_cl - 1)
    ybar = np.dot(w, y) / w.sum()
    tss = np.dot(w, (y - ybar) ** 2)
    r2 = 1.0 - np.dot(w, e**2) / tss if tss > 0 else np.nan
    return float(b[0]), float(np.sqrt(V[0, 0])), float(r2)
