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

def transform(
    df: pd.DataFrame,
    controls: list[str] | None = None,
    weights: pd.Series | None = None,
    treatments: list[str] | None = None,
) -> pd.DataFrame:
    """Normalise D, les contrôles et les autres traitements ; ajoute les
    poids et les colonnes de période factorisées.
    Traduction de twowayfeweights_transform.R."""
    controls = controls or []
    treatments = treatments or []

    _, df = normalize_var(df, "D")

    for control in controls:
        _, df = normalize_var(df, control)

    for treatment in treatments:
        _, df = normalize_var(df, treatment)

    df = df.copy()
    if weights is None:
        df["weights"] = 1
    else:
        df["weights"] = weights

    df["Tfactor"] = df["T"].astype("category")
    df["TFactorNum"] = df["Tfactor"].cat.codes + 1

    return df

def filter_data(
    df: pd.DataFrame,
    cmd_type: str,
    controls: list[str] | None = None,
    treatments: list[str] | None = None,
) -> pd.DataFrame:
    """Supprime les lignes avec des valeurs manquantes, selon le type
    d'estimation choisi.
    Traduction de twowayfeweights_filter.R."""
    controls = controls or []
    treatments = treatments or []

    if cmd_type != "fdTR":
        cols = ["Y", "G", "T", "D"] + controls + treatments
        df = df[df[cols].notna().all(axis=1)]
    else:
        tag1_ok = df[["D", "T", "Y"]].notna().all(axis=1)
        tag2_ok = df["D0"].notna()
        keep = tag1_ok | tag2_ok
        df = df[keep]
        tag1_ok = tag1_ok[keep]

        if controls:
            tag3_ok = df[controls].notna().all(axis=1)
            df = df[(~tag1_ok) | tag3_ok]

    return df

