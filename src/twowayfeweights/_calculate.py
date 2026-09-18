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


