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
