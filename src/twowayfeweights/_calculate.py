from __future__ import annotations

import pandas as pd


def compute_P_gt(df: pd.DataFrame) -> pd.Series:
    """Calcule P_gt : la part que représente chaque cellule (G,T) dans la
    somme totale des poids de l'échantillon.
    Traduction des lignes 'obs <- sum(dt$weights) ; P_gt := sum(weights)...'
    de twowayfeweights_calculate.R."""
    obs = df["weights"].sum()
    P_gt = df.groupby(["G", "T"])["weights"].transform("sum")
    return P_gt / obs

from twowayfeweights._kernels import weighted_mean


def compute_nat_weight(df: pd.DataFrame, D_col: str, P_gt: pd.Series) -> tuple[pd.Series, float]:
    """Calcule nat_weight = P_gt * D / mean_D, où mean_D est la moyenne
    pondérée de D sur tout l'échantillon.
    Traduction de : mean_D <- cpp_weighted_mean(...) ; nat_weight := P_gt * D / mean_D
    de twowayfeweights_calculate.R (branche type_TR).

    D_col : nom de la colonne à utiliser comme D (ex: "D" pour feTR, "D0" pour fdTR).

    Returns
    -------
    nat_weight : pd.Series
    mean_D : float
    """
    mean_D = weighted_mean(df[D_col], df["weights"])
    nat_weight = P_gt * df[D_col] / mean_D
    return nat_weight, mean_D

import pyfixest as pf


def fit_denom_regression(df: pd.DataFrame, controls: list[str] | None = None) -> pd.Series:
    """Régresse D sur les effets fixes (G, Tfactor) et d'éventuels contrôles ;
    renvoie les résidus (eps_1), alignés sur l'index de df.
    Traduction de la régression 'denom.lm' pour le cas type_fe (feTR/feS)
    dans twowayfeweights_calculate.R."""
    controls = controls or []
    rhs = " + ".join(controls) if controls else "1"
    formula = f"D ~ {rhs} | G + Tfactor"
    fit = pf.feols(formula, data=df, weights="weights")
    return pd.Series(fit.resid(), index=df.index)

def fit_beta_regression(df: pd.DataFrame, controls: list[str] | None = None) -> float:
    """Régresse Y sur D, les effets fixes (G, Tfactor) et d'éventuels
    contrôles ; renvoie le coefficient sur D (beta).
    Traduction de la régression 'beta.lm' dans twowayfeweights_calculate.R."""
    controls = controls or []
    rhs = " + ".join(["D"] + controls)
    formula = f"Y ~ {rhs} | G + Tfactor"
    fit = pf.feols(formula, data=df, weights="weights")
    return fit.coef()["D"]

def compute_W_feTR(
    df: pd.DataFrame, eps_1: pd.Series, D_col: str, mean_D: float
) -> pd.Series:
    """Calcule W = eps_1 * mean_D / denom_W, où denom_W est la moyenne
    pondérée de (eps_1 * D).
    Traduction de la branche feTR de twowayfeweights_calculate.R
    (denom_W et W)."""
    denom_W = weighted_mean(eps_1 * df[D_col], df["weights"])
    W = eps_1 * mean_D / denom_W
    return W

def finalize_feTR(df: pd.DataFrame, W: pd.Series, nat_weight: pd.Series) -> pd.DataFrame:
    """Calcule weight_result = W * nat_weight, l'attache au DataFrame, et
    ne garde qu'une ligne par cellule (G, Tfactor).
    Traduction des dernières lignes de la branche feTR dans
    twowayfeweights_calculate.R."""
    out = df.copy()
    out["W"] = W
    out["nat_weight"] = nat_weight
    out["weight_result"] = W * nat_weight
    out = out.groupby(["G", "Tfactor"], as_index=False, observed=True).first()
    return out

def calculate_feTR(df: pd.DataFrame, controls: list[str] | None = None) -> tuple[pd.DataFrame, float]:
    """Orchestre le calcul complet des poids pour le cas type='feTR'.
    Assemble compute_P_gt, compute_nat_weight, fit_denom_regression,
    fit_beta_regression, compute_W_feTR et finalize_feTR, dans le même
    ordre que la branche feTR de twowayfeweights_calculate.R.

    Returns
    -------
    result : pd.DataFrame
        Une ligne par cellule (G, Tfactor), avec les colonnes W, nat_weight,
        weight_result.
    beta : float
        Le coefficient sur D dans la régression Y ~ D + controls | G + Tfactor.
    """
    P_gt = compute_P_gt(df)
    nat_weight, mean_D = compute_nat_weight(df, "D", P_gt)

    eps_1 = fit_denom_regression(df, controls)
    beta = fit_beta_regression(df, controls)

    W = compute_W_feTR(df, eps_1, "D", mean_D)

    result = finalize_feTR(df, W, nat_weight)

    return result, beta

from twowayfeweights._kernels import rev_cumsum_by_group, feS_delta


def compute_E_eps_1_g_ge(df: pd.DataFrame, eps_1: pd.Series) -> pd.Series:
    """Calcule, pour chaque ligne, la moyenne pondérée des résidus eps_1
    à partir de cette période et pour toutes celles qui suivent, dans le
    même groupe (d'où 'g_ge' = 'groupe, greater or equal').
    Traduction de la section feS de twowayfeweights_calculate.R
    (E_eps_1_g_ge_aux, weights_aux, E_eps_1_g_ge).

    IMPORTANT : df doit déjà être trié par (G, Tfactor)."""
    eps_w = eps_1 * df["weights"]
    E_eps_1_g_ge_aux = rev_cumsum_by_group(df["G"], eps_w)
    weights_aux = rev_cumsum_by_group(df["G"], df["weights"])
    return E_eps_1_g_ge_aux / weights_aux

def apply_feS_delta(df: pd.DataFrame, P_gt: pd.Series) -> pd.DataFrame:
    """Identifie les cellules 'switchers' (où D change par rapport à la
    période précédente dans le même groupe), ne garde que ces lignes, et
    normalise nat_weight pour qu'il somme à 1 sur les switchers gardés.
    Traduction de la section feS de twowayfeweights_calculate.R
    (cpp_feS_delta, filtre keep, P_S).

    IMPORTANT : df doit déjà être trié par (G, TFactorNum)."""
    delta_res = feS_delta(df["G"], df["TFactorNum"], df["D"], P_gt)

    out = df.copy()
    out["delta_D"] = delta_res["delta_D"]
    out["s_gt"] = delta_res["s_gt"]
    out["abs_delta_D"] = delta_res["abs_delta_D"]
    out["nat_weight"] = delta_res["nat_weight"]

    out = out[delta_res["keep"].values]

    P_S = out["nat_weight"].sum()
    out["nat_weight"] = out["nat_weight"] / P_S

    return out

def finalize_feS(
    df: pd.DataFrame, E_eps_1_g_ge: pd.Series, P_gt: pd.Series
) -> pd.DataFrame:
    """Calcule om_tilde_1, le normalise en W via une moyenne pondérée par
    nat_weight, puis calcule weight_result = W * nat_weight.
    Traduction de la fin de la branche feS de twowayfeweights_calculate.R.

    IMPORTANT : df doit déjà être passé par apply_feS_delta (donc contenir
    s_gt et nat_weight, déjà filtré sur les switchers), et E_eps_1_g_ge /
    P_gt doivent être alignés sur le même index que df (déjà filtrés pareil)."""
    out = df.copy()
    out["om_tilde_1"] = out["s_gt"] * E_eps_1_g_ge / P_gt

    denom_W = weighted_mean(out["om_tilde_1"], out["nat_weight"])
    out["W"] = out["om_tilde_1"] / denom_W
    out["weight_result"] = out["W"] * out["nat_weight"]

    return out


