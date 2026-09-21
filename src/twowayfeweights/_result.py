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

def compute_sensibility2(
    G: pd.Series, T: pd.Series, W: pd.Series, nat_weight: pd.Series,
    weight_result: pd.Series, beta: float
) -> float | None:
    """Deuxième mesure de robustesse (Corollaire 1, point (ii)), calculée
    uniquement s'il existe des poids négatifs : le minimum de l'écart-type
    du traitement en dessous duquel beta pourrait être de signe différent
    de TOUS les effets traités.
    Traduction de la section 'sensibility2' de twowayfeweights_result.R.

    Returns None si aucun poids négatif n'existe (mesure non applicable)."""
    if (weight_result[weight_result.notna() & (weight_result < 0)]).sum() >= 0:
        return None

    mask = weight_result.notna() & (weight_result != 0)
    dat = pd.DataFrame({
        "G": G[mask], "T": T[mask], "W": W[mask],
        "nat_weight": nat_weight[mask], "weight_result": weight_result[mask],
    }).reset_index(drop=True)

    # tri croissant sur W, décroissant sur G puis T (pour les cumuls)
    dat = dat.sort_values(["W", "G", "T"], ascending=[True, False, False]).reset_index(drop=True)
    dat["Wsq"] = dat["nat_weight"] * dat["W"] ** 2
    dat["P_k"] = dat["nat_weight"].cumsum()
    dat["S_k"] = dat["weight_result"].cumsum()
    dat["T_k"] = dat["Wsq"].cumsum()

    # re-tri décroissant sur W, croissant sur G puis T (pour la lecture finale)
    dat = dat.sort_values(["W", "G", "T"], ascending=[False, True, True]).reset_index(drop=True)

    N = len(dat)
    dat["sens_measure2"] = abs(beta) / np.sqrt(dat["T_k"] + dat["S_k"] ** 2 / (1 - dat["P_k"]))
    dat["indicator"] = (dat["W"] < -dat["S_k"] / (1 - dat["P_k"])).astype(float)
    dat.loc[0, "indicator"] = 0
    dat["indicator_l"] = dat["indicator"].shift(1, fill_value=-1)
    dat["indicator"] = np.maximum(dat["indicator"], dat["indicator_l"])

    total_indicator = dat["indicator"].sum()
    idx = int(N - total_indicator - 1)  # -1 pour l'indexation 0-based de Python
    return dat["sens_measure2"].iloc[idx]
