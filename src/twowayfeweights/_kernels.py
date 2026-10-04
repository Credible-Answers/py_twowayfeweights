from __future__ import annotations

import numpy as np
import pandas as pd


def weighted_mean(x: pd.Series, w: pd.Series) -> float:
    """Moyenne pondérée, en ignorant les paires où x ou w est NA.
    Traduction de cpp_weighted_mean."""
    mask = x.notna() & w.notna()
    if not mask.any():
        return np.nan
    num = (x[mask] * w[mask]).sum()
    den = w[mask].sum()
    if den == 0:
        return np.nan
    return num / den


def rev_cumsum_by_group(g: pd.Series, x: pd.Series) -> pd.Series:
    """Cumul inversé (de la fin vers le début) à l'intérieur de chaque groupe.
    Traduction de cpp_rev_cumsum_by_group.
    IMPORTANT : suppose que les lignes sont déjà triées par (g, t)."""
    x_filled = x.fillna(0.0)
    reversed_cumsum = (
        x_filled[::-1]
        .groupby(g[::-1], sort=False)
        .cumsum()[::-1]
    )
    return reversed_cumsum

def fdtr_wtilde2(
    g: pd.Series, t: pd.Series, eps_2: pd.Series, P_gt: pd.Series
) -> pd.Series:
    """Pour chaque ligne, compare à la ligne suivante du même groupe si sa
    période suit bien (t+1). Sinon retombe sur eps_2.
    Traduction de cpp_fdtr_wtilde2.
    IMPORTANT : suppose que les lignes sont déjà triées par (g, t)."""
    g_next = g.shift(-1)
    t_next = t.shift(-1)
    eps_next = eps_2.shift(-1)
    P_next = P_gt.shift(-1)

    same_group = g_next == g
    consecutive = t_next == (t + 1)
    valid = (
        same_group
        & consecutive
        & P_gt.notna() & (P_gt != 0)
        & P_next.notna()
        & eps_next.notna()
        & eps_2.notna()
    )

    val = eps_2.copy()
    val[valid] = eps_2[valid] - eps_next[valid] * (P_next[valid] / P_gt[valid])
    val = val.where(val.notna() & np.isfinite(val), eps_2)
    return val


def feS_delta(
    g: pd.Series, t: pd.Series, D: pd.Series, P_gt: pd.Series
) -> pd.DataFrame:
    """Pour chaque ligne, compare à la ligne précédente du même groupe si sa
    période précède bien (t-1) : calcule delta_D, son signe, sa valeur
    absolue, et un poids naturel. Sinon la ligne est marquée à exclure.
    Traduction de cpp_feS_delta.
    IMPORTANT : suppose que les lignes sont déjà triées par (g, t)."""
    g_prev = g.shift(1)
    t_prev = t.shift(1)
    D_prev = D.shift(1)

    ok = (
        (g_prev == g)
        & (t - 1 == t_prev)
        & D.notna() & D_prev.notna()
        & t.notna() & t_prev.notna()
    )

    delta_D = pd.Series(np.nan, index=g.index)
    delta_D[ok] = D[ok] - D_prev[ok]

    abs_delta_D = delta_D.abs()

    s_gt = pd.Series(0, index=g.index, dtype=int)
    s_gt[ok] = np.sign(delta_D[ok]).astype(int)

    nat_weight = pd.Series(np.nan, index=g.index)
    nat_weight[ok] = P_gt[ok] * abs_delta_D[ok]

    return pd.DataFrame({
        "delta_D": delta_D,
        "s_gt": s_gt,
        "abs_delta_D": abs_delta_D,
        "nat_weight": nat_weight,
        "keep": ok,
    })
