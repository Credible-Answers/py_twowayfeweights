"""Result container and Stata-style printing."""

from __future__ import annotations

import sys
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

_FUNDING = (
    "The development of this package was funded by the European Union "
    "(ERC, REALLYCREDIBLE,GA N°101043899)."
)


@dataclass
class OtherTreatmentResult:
    """Weights attached to one of the ``other_treatments``."""

    name: str
    nr_plus: int
    nr_minus: int
    sum_plus: float
    sum_minus: float
    tot_cells: int

    @property
    def nr_weights(self) -> int:
        return self.nr_plus + self.nr_minus


@dataclass
class TwoWayFEWeightsResult:
    """Output of :func:`twowayfeweights`.

    Attributes mirror the R package (``nr_plus``, ``sum_minus``, ``sensibility``, ``mat`` ...);
    Stata's ``e()`` results are available through :attr:`M`, :attr:`lb_se_te` and :attr:`lb_se_te2`.
    """

    type: str
    treatment: str
    beta: float
    nr_plus: int
    nr_minus: int
    nr_weights: int
    sum_plus: float
    sum_minus: float
    tot_cells: int
    sensibility: float | None = None
    sensibility2: float | None = None
    mat: pd.DataFrame | None = None
    other_treatments: list[OtherTreatmentResult] = field(default_factory=list)
    weights: pd.DataFrame | None = field(default=None, repr=False)
    summary_measures: bool = False

    # ------------------------------------------------------------------ Stata-style accessors
    @property
    def M(self) -> pd.DataFrame:
        """Stata's ``e(M)`` / ``e(M1)``: number and sum of positive/negative weights."""
        return _mret(self.nr_plus, self.nr_minus, self.sum_plus, self.sum_minus)

    @property
    def lb_se_te(self) -> float | None:
        return self.sensibility

    @property
    def lb_se_te2(self) -> float | None:
        return self.sensibility2

    @property
    def n_zero(self) -> int:
        return max(0, self.tot_cells - self.nr_weights)

    # ------------------------------------------------------------------ legacy dict-style access
    def __getitem__(self, key: str):
        legacy = {
            "beta": self.beta,
            "n_pos": self.nr_plus,
            "n_neg": self.nr_minus,
            "n_atts": self.nr_weights,
            "sum_pos": self.sum_plus,
            "sum_neg": self.sum_minus,
            "tot_cells": self.tot_cells,
            "n_zero": self.n_zero,
            "sensibility": self.sensibility,
            "sensibility2": self.sensibility2,
            "dat_result": self.weights,
            "test_random_weights": {} if self.mat is None else {
                k: {"coef": r["Coef"], "se": r["SE"], "tstat": r["t-stat"], "corr": r["Correlation"]}
                for k, r in self.mat.iterrows()
            },
            "other_treatments_results": {
                o.name: {
                    "nr_plus": o.nr_plus, "nr_minus": o.nr_minus, "nr_weights": o.nr_weights,
                    "sum_plus": o.sum_plus, "sum_minus": o.sum_minus, "tot_cells": o.tot_cells,
                }
                for o in self.other_treatments
            },
        }
        if key in legacy:
            return legacy[key]
        return getattr(self, key)

    def to_dict(self) -> dict:
        return {
            "type": self.type,
            "beta": self.beta,
            "nr_plus": self.nr_plus,
            "nr_minus": self.nr_minus,
            "nr_weights": self.nr_weights,
            "sum_plus": self.sum_plus,
            "sum_minus": self.sum_minus,
            "tot_cells": self.tot_cells,
            "sensibility": self.sensibility,
            "sensibility2": self.sensibility2,
        }

    # ------------------------------------------------------------------ printing
    def summary(self) -> str:
        lines: list[str] = [""]
        n_att = self.nr_weights
        beta_s = f"{self.beta:.4f}"
        if self.other_treatments:
            lines += [
                "Under the common trends assumption,",
                f"the TWFE coefficient beta, equal to {beta_s}, estimates the sum of several terms.",
                "",
                f"The first term is a weighted sum of {n_att} ATTs of the treatment.",
                f"{self.nr_plus} ATTs receive a positive weight, and {self.nr_minus} receive a negative weight.",
            ]
        else:
            if self.type in ("feTR", "fdTR"):
                lines.append("Under the common trends assumption,")
            else:
                lines += [
                    "Under the common trends assumption and",
                    "the assumption that groups' treatment effects do not change over time,",
                ]
            lines += [
                f"the TWFE coefficient beta, equal to {beta_s}, estimates a weighted sum of {n_att} ATTs.",
                f"{self.nr_plus} ATTs receive a positive weight, and {self.nr_minus} receive a negative weight.",
            ]
        lines += _zero_line(self.tot_cells, n_att)
        lines += _table(f"Treat. var: {self.treatment}", self.nr_plus, self.nr_minus, self.sum_plus, self.sum_minus)

        for j, o in enumerate(self.other_treatments, start=1):
            lines += [
                "",
                f"The next term is a weighted sum of {o.nr_weights} ATTs of treatment {j} "
                "included in the other_treatments option.",
                f"{o.nr_plus} ATTs receive a positive weight, and {o.nr_minus} receive a negative weight.",
            ]
            lines += _zero_line(o.tot_cells, o.nr_weights)
            lines += _table(f"Other treat.: {o.name}", o.nr_plus, o.nr_minus, o.sum_plus, o.sum_minus)

        if self.summary_measures and not self.other_treatments:
            sub = self.type[:2]
            lines += [
                "",
                "Summary Measures:",
                f"TWFE coefficient (β_{sub}) = {_f4(self.beta)}",
                f"min σ(Δ) compatible with β_{sub} and Δ_TR = 0: {_f4(self.sensibility)}",
            ]
            if self.sum_minus < 0:
                lines.append(
                    f"min σ(Δ) compatible with treatment effect of opposite sign than β_{sub} "
                    f"in all (g,t) cells: {_f4(self.sensibility2)}"
                )
            lines.append("Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)")

        if self.mat is not None:
            title = "Regression of variables possibly correlated with the treatment effect on the weights"
            if self.other_treatments:
                title += " attached to the treatment"
            lines += ["", title, _format_mat(self.mat)]

        lines += ["", "", _FUNDING]
        return "\n".join(lines)

    def __str__(self) -> str:
        return _encodable(self.summary())

    def _repr_pretty_(self, p, cycle):  # IPython
        p.text(self.summary())


def print_twowayfeweights(result: TwoWayFEWeightsResult, D_name: str | None = None, type: str | None = None) -> None:
    """Print a result (kept for backward compatibility; ``print(result)`` does the same)."""
    if D_name is not None:
        result.treatment = D_name
    print(_encodable(result.summary()))


# ---------------------------------------------------------------------- formatting helpers

_ASCII = {"Σ": "Sum", "β": "beta", "σ": "sd", "Δ": "Delta", "°": "o"}


def _encodable(text: str) -> str:
    """Replace Greek letters by ASCII names when stdout cannot encode them (e.g. cp1252 consoles)."""
    enc = getattr(sys.stdout, "encoding", None) or "utf-8"
    try:
        text.encode(enc)
        return text
    except (UnicodeEncodeError, LookupError):
        for k, v in _ASCII.items():
            text = text.replace(k, v)
        return text.encode(enc, errors="replace").decode(enc)


def _mret(npl, nmi, spl, smi) -> pd.DataFrame:
    return pd.DataFrame(
        {"N_ATTs": [npl, nmi, npl + nmi], "Sum_weights": [spl, smi, spl + smi]},
        index=["Pos_Weights", "Neg_Weights", "Tot"],
    )


def _f4(x) -> str:
    return "." if x is None or not np.isfinite(x) else f"{x:.4f}"


def _zero_line(tot_cells: int, n_att: int) -> list[str]:
    if tot_cells > n_att:
        return [
            f"{tot_cells} (g,t) cells receive the treatment, but the ATTs of "
            f"{tot_cells - n_att} cells receive a weight equal to zero."
        ]
    return []


def _table(label: str, npl: int, nmi: int, spl: float, smi: float) -> list[str]:
    w1 = max(24, len(label) + 2)
    width = w1 + 24
    rule = "-" * max(48, width)

    def row(a, b, c):
        return f"{a:<{w1}}{b:<12}{c:<12}".rstrip()

    return [
        rule,
        row(label, "# ATTs", "Σ weights"),
        rule,
        row("Positive weights", str(npl), f"{spl:.4f}"),
        row("Negative weights", str(nmi), f"{smi:.4f}"),
        rule,
        row("Total", str(npl + nmi), f"{spl + smi:.4f}"),
        rule,
    ]


def _format_mat(mat: pd.DataFrame) -> str:
    cols = list(mat.columns)
    w0 = max(8, max(len(str(i)) for i in mat.index) + 2)
    head = " " * w0 + "".join(f"{c:>13}" for c in cols)
    body = [
        f"{str(i):<{w0}}" + "".join(f"{v:>13.7g}" if np.isfinite(v) else f"{'.':>13}" for v in r)
        for i, r in zip(mat.index, mat.to_numpy())
    ]
    return "\n".join([head, *body])
