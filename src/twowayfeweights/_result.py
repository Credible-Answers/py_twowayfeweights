from __future__ import annotations

import numpy as np
import pandas as pd


def summarize_weights(weight: pd.Series) -> dict:
    """Compte les poids positifs/négatifs et leurs sommes.
    Traduction de twowayfeweights_summarize_weights (utils.R)."""
    ok = weight.notna()
    weight_plus = weight[ok & (weight > 0)]
    weight_minus = weight[ok & (weight < 0)]

    return {
        "nr_plus": len(weight_plus),
        "nr_minus": len(weight_minus),
        "nr_weights": len(weight_plus) + len(weight_minus),
        "sum_plus": weight_plus.sum(),
        "sum_minus": weight_minus.sum(),
    }

def compute_sensibility(W: pd.Series, nat_weight: pd.Series, beta: float) -> float:
    """Première mesure de robustesse (Corollaire 1, point (i)) : le minimum
    de l'écart-type du traitement en dessous duquel beta et l'ATT
    pourraient être de signes opposés.
    Traduction de la section 'W_mean, W_sd, sensibility' de
    twowayfeweights_result.R."""
    mask = W.notna() & nat_weight.notna()
    W_ok = W[mask]
    nat_weight_ok = nat_weight[mask]

    W_mean = np.average(W_ok, weights=nat_weight_ok)
    M = (nat_weight_ok != 0).sum()
    W_sd = np.sqrt((nat_weight_ok * (W_ok - W_mean) ** 2).sum()) * np.sqrt(M / (M - 1))

    return abs(beta) / W_sd

