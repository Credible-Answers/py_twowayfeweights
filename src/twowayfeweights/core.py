import pandas as pd


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

    Returns
    -------
    dict
        Résultats de l'estimation (à définir plus précisément).
    """
    # TODO: implémenter la logique de calcul des poids
    raise NotImplementedError("La fonction n'est pas encore implémentée.")
