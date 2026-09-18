from __future__ import annotations

import pandas as pd


def normalize_var(df: pd.DataFrame, varname: str) -> tuple[bool, pd.DataFrame]:
    """Remplace les valeurs de `varname` par leur moyenne dans la cellule
    (G, T) si la variable varie à l'intérieur d'au moins une cellule.
    Traduction de twowayfeweights_normalize_var.R.

    Returns
    -------
    changed : bool
        True si une normalisation a été effectuée.
    df : pd.DataFrame
        Le DataFrame (modifié si changed=True).
    """
    cell_stats = df.groupby(["G", "T"])[varname].agg(["mean", "std"])
    any_within_var = cell_stats["std"].sum(skipna=True) > 0

    if any_within_var:
        cell_means = cell_stats["mean"]
        df = df.copy()
        df[varname] = df.set_index(["G", "T"]).index.map(cell_means)

    return any_within_var, df

def rename_var(
    df: pd.DataFrame,
    Y: str,
    G: str,
    T: str,
    D: str,
    D0: str | None = None,
    controls: list[str] | None = None,
    treatments: list[str] | None = None,
    random_weights: list[str] | None = None,
) -> pd.DataFrame:
    """Renomme les colonnes choisies par l'utilisateur vers les noms
    internes standardisés (Y, G, T, D, D0, ctrl_*, OT_*, RW_*).
    Traduction de twowayfeweights_rename_var.R."""
    controls = controls or []
    treatments = treatments or []
    random_weights = random_weights or []

    controls_renamed = [f"ctrl_{c}" for c in controls]
    treatments_renamed = [f"OT_{t}" for t in treatments]

    original_names = [Y, G, T, D] + controls + treatments
    new_names = ["Y", "G", "T", "D"] + controls_renamed + treatments_renamed

    if D0 is not None:
        original_names = original_names + [D0]
        new_names = new_names + ["D0"]

    out = df[original_names].copy()
    out.columns = new_names

    if random_weights:
        rw_renamed = [f"RW_{w}" for w in random_weights]
        rw_df = df[random_weights].copy()
        rw_df.columns = rw_renamed
        out = pd.concat([out, rw_df], axis=1)

    return out
