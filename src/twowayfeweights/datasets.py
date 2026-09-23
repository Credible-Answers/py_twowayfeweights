"""Example datasets shipped with the package."""

from __future__ import annotations

from importlib import resources

import pandas as pd


def load_wagepan() -> pd.DataFrame:
    """Wooldridge's ``wagepan`` panel (545 men, 1980-1987), as used in the Stata/R examples."""
    with resources.files(__package__).joinpath("data/wagepan.csv.gz").open("rb") as fh:
        return pd.read_csv(fh, compression="gzip")
