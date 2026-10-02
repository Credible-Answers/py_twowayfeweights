"""Weights attached to two-way fixed effects regressions (de Chaisemartin & D'Haultfoeuille, 2020).

Python port of the Stata / R package ``twowayfeweights``.
"""

from ._core import twowayfeweights
from ._result import OtherTreatmentResult, TwoWayFEWeightsResult, print_twowayfeweights
from .datasets import load_wagepan

__all__ = [
    "twowayfeweights",
    "print_twowayfeweights",
    "TwoWayFEWeightsResult",
    "OtherTreatmentResult",
    "load_wagepan",
]
__version__ = "0.1.0"
