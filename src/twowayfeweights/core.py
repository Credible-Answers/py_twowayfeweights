from __future__ import annotations

import pandas as pd

from twowayfeweights._prepare import rename_var, transform, filter_data
from twowayfeweights._calculate import (
    calculate_feTR,
    calculate_feS,
    calculate_fdTR,
    calculate_fdS,
)

_CALCULATORS = {
    "feTR": calculate_feTR,
    "feS": calculate_feS,
    "fdTR": calculate_fdTR,
    "fdS": calculate_fdS,
}

# Types qui ont besoin d'être triés par (G, TFactorNum) avant le calcul,
# car ils comparent chaque ligne à sa voisine du même groupe.
_NEEDS_SORT = {"feS", "fdTR"}


def twowayfeweights(
    data: pd.DataFrame,
    Y: str,
    G: str,
    T: str,
    D: str,
    type: str = "feTR",
    D0: str | None = None,
    summary_measures: bool = False,
    controls: list[str] | None = None,
    weights: str | None = None,
    other_treatments: list[str] | None = None,
    test_random_weights: list[str] | None = None,
    path: str | None = None,
):
    """
    Estime les poids attachés aux régressions à effets fixes bidirectionnels
    (two-way fixed effects), suivant de Chaisemartin & D'Haultfoeuille (2020a).

    Portage Python du package R/Stata twowayfeweights.

    Parameters
    ----------
    data : pd.DataFrame
        Le jeu de données.
    Y : str
        Nom de la colonne de la variable dépendante.
    G : str
        Nom de la colonne identifiant les groupes.
    T : str
        Nom de la colonne identifiant les périodes.
    D : str
        Nom de la colonne du traitement.
    type : str
        Un parmi {"feTR", "feS", "fdTR", "fdS"}.
    D0 : str, optional
        Requis si type="fdTR". Niveau initial du traitement.

    Returns
    -------
    dict
        Un dictionnaire contenant :
        - "weights" : pd.DataFrame avec une ligne par cellule (G, T) et son
          poids (colonne "weight").
        - "beta" : le coefficient TWFE estimé.
    """
    if type not in _CALCULATORS:
        raise ValueError(f"type doit être un de {list(_CALCULATORS)}, reçu '{type}'.")

    if type == "fdTR" and D0 is None:
        raise ValueError("Le paramètre D0 est requis quand type='fdTR'.")

    if other_treatments and type != "feTR":
        raise ValueError("other_treatments ne peut être utilisé qu'avec type='feTR'.")

    controls = controls or []
    other_treatments = other_treatments or []
    test_random_weights = test_random_weights or []

    weights_col = data[weights] if weights is not None else None

    renamed = rename_var(
        data, Y=Y, G=G, T=T, D=D, D0=D0,
        controls=controls, treatments=other_treatments,
        random_weights=test_random_weights,
    )
    transformed = transform(
        renamed, controls=[f"ctrl_{c}" for c in controls],
        weights=weights_col, treatments=[f"OT_{t}" for t in other_treatments],
    )
    filtered = filter_data(
        transformed, cmd_type=type,
        controls=[f"ctrl_{c}" for c in controls],
        treatments=[f"OT_{t}" for t in other_treatments],
    )

    if type in _NEEDS_SORT:
        filtered = filtered.sort_values(["G", "TFactorNum"]).reset_index(drop=True)

    calculator = _CALCULATORS[type]
    result, beta = calculator(filtered, controls=[f"ctrl_{c}" for c in controls])

    weights_df = result[["G", "T", "weight_result"]].rename(columns={"weight_result": "weight"})

    return {"weights": weights_df, "beta": beta}
