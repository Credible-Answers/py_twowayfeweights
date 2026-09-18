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
