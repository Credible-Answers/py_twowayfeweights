# twowayfeweights (Python)

[![tests](https://github.com/anzonyquispe/py_twowayfeweights/actions/workflows/tests.yml/badge.svg)](https://github.com/anzonyquispe/py_twowayfeweights/actions/workflows/tests.yml)

Estimates the weights and the sensitivity measures attached to two-way fixed effects (TWFE) and
first-difference regressions, following de Chaisemartin & D'Haultfoeuille (2020). This is a Python
port of the Stata and R command `twowayfeweights`. It reproduces the Stata results: see
[Parity with Stata](#parity-with-stata).

A TWFE coefficient estimates a weighted sum of the treatment effects in each (g,t) cell. Some of those
weights can be negative. When they are, the coefficient can have the opposite sign of every
individual treatment effect. `twowayfeweights` computes these weights, reports how many are negative
and how much they sum to, and measures how robust the coefficient is to heterogeneous effects.

## Installation

```bash
pip install twowayfeweights
```

Requires Python ≥ 3.9. The only dependencies are `numpy`, `pandas` and `scipy`.

## Quick start

```python
from twowayfeweights import twowayfeweights, load_wagepan

df = load_wagepan()
res = twowayfeweights(df, Y="lwage", G="nr", T="year", D="union",
                      type="feTR", summary_measures=True, test_random_weights="educ")
print(res)
```

```
Under the common trends assumption,
the TWFE coefficient beta, equal to 0.1066, estimates a weighted sum of 967 ATTs.
820 ATTs receive a positive weight, and 147 receive a negative weight.
1016 (g,t) cells receive the treatment, but the ATTs of 49 cells receive a weight equal to zero.
------------------------------------------------
Treat. var: union       # ATTs      Σ weights
------------------------------------------------
Positive weights        820         1.0105
Negative weights        147         -0.0105
------------------------------------------------
Total                   967         1.0000
------------------------------------------------

Summary Measures:
TWFE coefficient (β_fe) = 0.1066
min σ(Δ) compatible with β_fe and Δ_TR = 0: 0.0969
min σ(Δ) compatible with treatment effect of opposite sign than β_fe in all (g,t) cells: 3.1759
Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

Regression of variables possibly correlated with the treatment effect on the weights
                 Coef           SE       t-stat  Correlation
educ       -0.1344553   0.07136021    -1.884177   -0.1182587
```

## Usage

```python
twowayfeweights(data, Y, G, T, D, type="feTR", D0=None, summary_measures=False,
                controls=None, weights=None, other_treatments=None,
                test_random_weights=None, path=None)
```

| argument | Stata equivalent | description |
|---|---|---|
| `Y, G, T, D` | `varlist` | outcome, group, time period, treatment |
| `type` | `type()` | `"feTR"`, `"feS"`, `"fdTR"` or `"fdS"` (see below) |
| `D0` | 5th variable | treatment level (not differenced); required with `fdTR` |
| `controls` | `controls()` | control variables |
| `weights` | `weight()` | analytic weights |
| `other_treatments` | `other_treatments()` | other treatments in the regression (`feTR` only) |
| `test_random_weights` | `test_random_weights()` | variables regressed on the weights |
| `summary_measures` | `summary_measures` | print the sensitivity measures (they are always computed) |
| `path` | `path()` | save the (g,t) weights to `.csv`, `.dta` or `.parquet` |

`data` can be a pandas DataFrame, or anything with a `.to_pandas()` method (e.g. a polars DataFrame).
Group and time identifiers can be numeric, strings or categoricals.

**Estimation types**

* `feTR`: regression with group and period fixed effects, under the common trends assumption.
* `feS`: the same regression, also assuming that each group's treatment effect does not change over time.
* `fdTR`: first-difference regression under common trends. `Y` and `D` are the first differences,
  and `D0` is the treatment level.
* `fdS`: first-difference regression, also assuming stable treatment effects.

**The result** is a `TwoWayFEWeightsResult`. Its attributes follow the R package:

| attribute | content |
|---|---|
| `beta` | the TWFE / FD coefficient |
| `nr_plus`, `nr_minus`, `nr_weights` | number of positive, negative, and nonzero weights |
| `sum_plus`, `sum_minus` | sum of positive and of negative weights |
| `tot_cells` | number of (g,t) cells that receive the treatment |
| `sensibility`, `sensibility2` | the two sensitivity measures (Stata `e(lb_se_te)`, `e(lb_se_te2)`) |
| `mat` | `test_random_weights` table (Coef, SE, t-stat, Correlation) |
| `other_treatments` | one `OtherTreatmentResult` per other treatment |
| `weights` | DataFrame of the weight in every (g,t) cell (what Stata's `path()` saves) |
| `M` | Stata's `e(M)` matrix |

`print(res)` shows the Stata-style table. `res.to_dict()` returns the scalars.
`print_twowayfeweights(res)` and dictionary access such as `res["beta"]` are kept for backward compatibility.

### Example: first-difference regression with many controls

This reproduces chapter 5 of the [DiD book](https://anzonyquispe.github.io/did_book/chapters/ch05.html),
using data from Gentzkow, Shapiro & Sinkinson (2011):

```python
styr = [c for c in df.columns if c.startswith("styr") and c != "styr"]   # 683 state-year dummies
res = twowayfeweights(df, "changeprestout", "cnty90", "year", "changedailies",
                      type="fdTR", D0="numdailies", controls=styr, test_random_weights="year")
# beta = 0.0026: 5371 positive weights (sum 2.4271), 4505 negative weights (sum -1.4271)
```

## Parity with Stata

The test suite (`tests/`) checks 36 specifications (17 of them `feS`/`fdS`) against output from the Stata command, run in
batch mode by `tests/reference/build_reference.py`. The specifications cover all four types, with
and without controls, weights, other treatments and random-weight tests. They use the wagepan data,
the Gentzkow et al. data with 683 controls, and a simulated unbalanced panel with gaps, several
observations per cell and missing values. The checks:

* the numbers of positive, negative and treated-cell weights match **exactly**;
* beta matches to ~1e-8 relative;
* the weight sums, the sensitivity measures and every (g,t) weight match to ~1e-7 relative.

On the 20 specifications where the R package runs and agrees with Stata, Python matches R to 1e-8 or better.

The small remaining differences come from Stata. It stores intermediate variables (`gen`, `predict`,
`egen`) as 4-byte floats, while this package computes everything in double precision.

One known difference: for `type="feS"` with **several observations per (g,t) cell**, Stata builds the
running average of residuals one observation at a time, and then keeps an arbitrary observation of
each cell. Its output then depends on the row order, and it changes from run to run because Stata
breaks sort ties randomly. This package uses the formula from the paper: all observations in the
periods t' ≥ t. With one observation per cell, which is the usual case, the results are identical.

## Performance

Regressions run on (g,t) cells with summed weights, not on individual observations. The time
dummies are never built as a dense matrix: after the group effects are absorbed, their normal
equations come from a sparse group × period weight table. The system is solved with a scaled
pseudo-inverse, which handles collinearity, plus iterative refinement. Everything else is vectorized
numpy code.

Seconds per call on the same data, on one Windows laptop (Stata 18 MP, R 4.5 with TwoWayFEWeights 2.1.0):

| data | type | Stata | R | Python |
|---|---|---:|---:|---:|
| 200k rows (10k groups × 20 periods, 5 controls, weights) | feTR | 3.6 | 16.4 | 0.17 |
| | feS | 1.8 | 13.4 | 0.20 |
| | fdTR | 3.1 | 19.5 | 0.14 |
| | fdS | 0.9 | 10.5 | 0.13 |
| 1M rows (50k groups × 20 periods, 5 controls, weights) | feTR | 12.9 | 77.8 | 1.05 |
| | feS | 9.6 | 56.6 | 1.16 |
| | fdTR | 12.6 | 82.9 | 1.17 |
| | fdS | 7.1 | 50.2 | 0.76 |
| Gentzkow et al. (16,872 rows, 683 controls) | feTR | 27.4 | 1245 | 0.58 |
| | feS | 29.5 | – | 0.59 |
| | fdTR | 24.0 | – | 0.54 |
| | fdS | 24.0 | – | 0.55 |

`python benchmarks/bench.py` times panels of up to 3M rows.

## References

* de Chaisemartin, C. and D'Haultfoeuille, X. (2020). Two-Way Fixed Effects Estimators with
  Heterogeneous Treatment Effects. *American Economic Review*, 110(9), 2964–2996.
* de Chaisemartin, C. and D'Haultfoeuille, X. (2023). Two-way fixed effects and differences-in-differences
  estimators with several treatments. *Journal of Econometrics*.

Stata: `ssc install twowayfeweights`. R: `install.packages("TwoWayFEWeights")`.

The development of the original package was funded by the European Union (ERC, REALLYCREDIBLE, GA N°101043899).
