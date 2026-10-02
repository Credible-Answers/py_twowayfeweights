from __future__ import annotations

import pandas as pd

from twowayfeweights._prepare import rename_var, transform, filter_data
from twowayfeweights._calculate import (
    calculate_feTR,
    calculate_feS,
    calculate_fdTR,
    calculate_fdS,
)
from twowayfeweights._result import (
    summarize_weights,
    compute_sensibility,
    compute_sensibility2,
    test_random_weights as run_test_random_weights,
)

_CALCULATORS = {
    "feTR": calculate_feTR,
    "feS": calculate_feS,
    "fdTR": calculate_fdTR,
    "fdS": calculate_fdS,
}

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
    Y, G, T, D : str
        Noms des colonnes de la variable dépendante, du groupe, de la
        période, et du traitement.
    type : str
        Un parmi {"feTR", "feS", "fdTR", "fdS"}.
    D0 : str, optional
        Requis si type="fdTR". Niveau initial du traitement.
    summary_measures : bool
        Si True, ajoute au résultat les mesures de robustesse (Corollaire 1).
    controls : list[str], optional
        Variables de contrôle.
    weights : str, optional
        Nom d'une colonne de poids de régression.
    other_treatments : list[str], optional
        Autres traitements (feTR uniquement). Non encore implémenté.
    test_random_weights : list[str], optional
        Variables dont on teste la corrélation avec les poids.
    path : str, optional
        Chemin où sauvegarder les poids en CSV.

    Returns
    -------
    dict
        - "weights" : pd.DataFrame avec une ligne par cellule (G, T) et son
          poids (colonne "weight").
        - "beta" : le coefficient TWFE estimé.
        - "summary" : dict (nr_plus, nr_minus, sum_plus, sum_minus,
          sensibility, sensibility2), présent si summary_measures=True.
        - "random_weights_test" : pd.DataFrame, présent si
          test_random_weights est fourni.
    """
    if type not in _CALCULATORS:
        raise ValueError(f"type doit être un de {list(_CALCULATORS)}, reçu '{type}'.")

    if type == "fdTR" and D0 is None:
        raise ValueError("Le paramètre D0 est requis quand type='fdTR'.")

    if other_treatments:
        raise NotImplementedError("other_treatments n'est pas encore implémenté.")

    controls = controls or []
    test_random_weights = test_random_weights or []

    weights_col = data[weights] if weights is not None else None

    renamed = rename_var(
        data, Y=Y, G=G, T=T, D=D, D0=D0,
        controls=controls, treatments=[],
        random_weights=test_random_weights,
    )
    transformed = transform(
        renamed, controls=[f"ctrl_{c}" for c in controls],
        weights=weights_col, treatments=[],
    )
    filtered = filter_data(
        transformed, cmd_type=type,
        controls=[f"ctrl_{c}" for c in controls],
        treatments=[],
    )

    if type in _NEEDS_SORT:
        filtered = filtered.sort_values(["G", "TFactorNum"]).reset_index(drop=True)

    calculator = _CALCULATORS[type]
    result, beta = calculator(filtered, controls=[f"ctrl_{c}" for c in controls])

    weights_df = result[["G", "T", "weight_result"]].rename(columns={"weight_result": "weight"})

    output = {"weights": weights_df, "beta": beta}

    if summary_measures:
        summary = summarize_weights(result["weight_result"])
        summary["sensibility"] = compute_sensibility(result["W"], result["nat_weight"], beta)
        summary["sensibility2"] = compute_sensibility2(
            result["G"], result["T"], result["W"], result["nat_weight"],
            result["weight_result"], beta,
        )
        output["summary"] = summary

    if test_random_weights:
        rw_cols = [f"RW_{v}" for v in test_random_weights]
        for c in rw_cols:
            result[c] = renamed.loc[result.index, c]
        output["random_weights_test"] = run_test_random_weights(result, rw_cols)

    if path is not None:
        weights_df.to_csv(path, index=False)

    return output
