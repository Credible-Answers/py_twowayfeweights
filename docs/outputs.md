# Outputs

`twowayfeweights()` returns a result object. Printing it shows the Stata-style report; every number
in that report is also stored in the object, so you can use it in your own code.

```python
# pip install twowayfeweights
from twowayfeweights import twowayfeweights, load_wagepan

df = load_wagepan()
res = twowayfeweights(df, Y="lwage", G="nr", T="year", D="union",
                      test_random_weights=["educ", "exper"])

print(res)        # the full report, as in Stata
res.beta          # 0.10662746555843541
res.nr_minus      # 147 cells with a negative weight
res.sum_minus     # -0.010528987131443892
res.sensibility   # 0.09686909502992004
```

## All the statistics

| attribute | what it is | Stata |
|---|---|---|
| `res.beta` | The TWFE (or first-difference) coefficient. | `e(beta)` |
| `res.nr_plus` | Number of treated cells with a **positive** weight. | `e(M)[1,1]` |
| `res.nr_minus` | Number of treated cells with a **negative** weight. | `e(M)[2,1]` |
| `res.nr_weights` | Number of nonzero weights (`nr_plus + nr_minus`). | `e(M)[3,1]` |
| `res.sum_plus` | Sum of the positive weights. | `e(M)[1,2]` |
| `res.sum_minus` | Sum of the negative weights. | `e(M)[2,2]` |
| `res.tot_cells` | Number of (group, period) cells that receive the treatment. Cells whose weight is exactly zero are counted here but not in `nr_weights`. | |
| `res.sensibility` | Summary measure 1: the smallest standard deviation of the treatment effects across cells under which beta is compatible with an average effect of zero. Small values mean beta is fragile. | `e(lb_se_te)` |
| `res.sensibility2` | Summary measure 2: the smallest standard deviation under which the treatment effect could have the opposite sign of beta in every cell. `None` when no weight is negative. | `e(lb_se_te2)` |
| `res.mat` | The `test_random_weights` table (see below). `None` if the option is not used. | `e(randomweightstest1)` |
| `res.other_treatments` | One entry per variable in `other_treatments` (see below). | `e(M2)`, `e(M3)`... |
| `res.weights` | The weight of every (group, period) cell (see below). | file saved by `path()` |
| `res.M` | The table of counts and sums, laid out like Stata's `e(M)`. | `e(M)` |
| `res.notes` | The notes printed before the results (e.g. a treatment that varies within cells). | |
| `res.type`, `res.treatment` | The `type` used and the name of the treatment variable. | |

The summary measures are always computed; `summary_measures=True` only adds them to the printed
report. With `other_treatments`, the summary measures are not defined and are `None`.

## All the numbers at once

```python
res.to_dict()
```

```text
{'type': 'feTR', 'beta': 0.10662746555843541, 'nr_plus': 820, 'nr_minus': 147, 'nr_weights': 967,
 'sum_plus': 1.0105289871314438, 'sum_minus': -0.010528987131443892, 'tot_cells': 1016,
 'sensibility': 0.09686909502992004, 'sensibility2': 3.175858897581741}
```

To put several results in one table, collect the dictionaries:

```python
import pandas as pd

rows = []
for t in ["feTR", "feS"]:
    r = twowayfeweights(df, Y="lwage", G="nr", T="year", D="union", type=t)
    rows.append(r.to_dict())
pd.DataFrame(rows)
```

## Counts and sums (Stata's `e(M)`)

```python
res.M
```

```text
             N_ATTs  Sum_weights
Pos_Weights     820     1.010529
Neg_Weights     147    -0.010529
Tot             967     1.000000
```

## Are the weights correlated with other variables? (`res.mat`)

With `test_random_weights`, each variable is regressed on the weights. `res.mat` is a DataFrame with
one row per variable:

| column | meaning |
|---|---|
| `Coef` | coefficient of the regression of the variable on the weights |
| `SE` | its standard error, clustered by group |
| `t-stat` | `Coef / SE` |
| `Correlation` | correlation between the variable and the weights |

```python
res.mat
res.mat.loc["educ", "SE"]     # one number
```

```text
           Coef        SE    t-stat  Correlation
educ  -0.134455  0.071360 -1.884177    -0.118259
exper -0.201326  0.105611 -1.906288    -0.083915
```

Like Stata, the command does not report p-values. A two-sided p-value from the t-stat (normal
approximation), if you want one:

```python
from scipy import stats
2 * stats.norm.sf(abs(res.mat["t-stat"]))
```

## The weight of every cell (`res.weights`)

A DataFrame with one row per (group, period) cell, the same content that Stata saves with `path()`:

| column | meaning |
|---|---|
| `Group_TWFE` | group |
| `Time_TWFE` | period |
| `weight` | the weight of the cell (0 for untreated cells) |

With `other_treatments`, the columns are `Group`, `Time`, `weight` (the main treatment) and
`weight_others1`, `weight_others2`... (one per other treatment).

```python
w = res.weights
w[w["weight"] < 0]            # the cells with a negative weight
w.groupby("Time_TWFE")["weight"].sum()   # total weight by period
```

To save it to a file, pass `path` (or use pandas directly):

```python
res = twowayfeweights(df, Y="lwage", G="nr", T="year", D="union", path="weights.csv")  # or .dta, .parquet
res.weights.to_excel("weights.xlsx")   # needs openpyxl
```

## Other treatments (`res.other_treatments`)

With `other_treatments`, `res.other_treatments` is a list with one entry per other treatment, in
the order given. Each entry has the same fields as the main treatment:

```python
res = twowayfeweights(df, Y="lwage", G="nr", T="year", D="union", other_treatments=["married"])
o = res.other_treatments[0]
o.name, o.nr_plus, o.nr_minus, o.sum_plus, o.sum_minus, o.tot_cells
```

```text
('married', 917, 997, 0.4827791014372501, -0.4827791014372501, 1914)
```

## Code written for older versions

Dictionary access with the names of the previous Python version still works, e.g. `res["beta"]`,
`res["n_pos"]`, `res["n_neg"]`, `res["sum_neg"]`, `res["dat_result"]` (the cell weights) and
`res["test_random_weights"]["educ"]["coef"]`. `print_twowayfeweights(res)` prints the report.
