# Arguments

```python
from twowayfeweights import twowayfeweights

res = twowayfeweights(data, Y, G, T, D, type="feTR", D0=None, summary_measures=False,
                      controls=None, weights=None, other_treatments=None,
                      test_random_weights=None, path=None)
```

Variables are always given by their **column name** in `data`. Where an argument accepts several
variables, you can pass one name (`"educ"`) or a list (`["educ", "exper"]`).

## Required

| argument | Stata | what it is |
|---|---|---|
| `data` | the dataset in memory | A **pandas DataFrame** with one row per observation. Other data frames must be converted first, e.g. `df.to_pandas()` for polars. |
| `Y` | 1st variable | The outcome. For `fdTR` and `fdS`, the **first difference** of the outcome. |
| `G` | 2nd variable | The group identifier (person, county, firm...). Numbers, strings or categories. |
| `T` | 3rd variable | The time period. Must be sortable (years, months, 1, 2, 3...). |
| `D` | 4th variable | The treatment. For `fdTR` and `fdS`, the **first difference** of the treatment. |

## Options

| argument | Stata | default | what it does |
|---|---|---|---|
| `type` | `type()` | `"feTR"` | Which regression and which assumptions. See below. |
| `D0` | 5th variable | `None` | The treatment **level** (not differenced). Required with `type="fdTR"`, ignored otherwise. |
| `summary_measures` | `summary_measures` | `False` | Print the two sensitivity measures. They are always computed and stored, whatever this option. |
| `controls` | `controls()` | `None` | Control variables included in the regression. |
| `weights` | `weight()` | `None` | A variable of analytic weights. Must not be negative. |
| `other_treatments` | `other_treatments()` | `None` | Other treatment variables included in the regression. Only with `type="feTR"`. The output then shows the weights attached to each treatment. |
| `test_random_weights` | `test_random_weights()` | `None` | Variables regressed on the weights, to see whether the weights are correlated with things that may drive the treatment effect (e.g. education, age). |
| `path` | `path()` | `None` | Save the weight of every (group, period) cell to a file: `.csv`, `.dta` (Stata) or `.parquet` (needs the `pyarrow` package). |

## The four types

| `type` | regression | assumptions |
|---|---|---|
| `"feTR"` | outcome on treatment, with group and period fixed effects | common trends |
| `"feS"` | same as `feTR` | common trends, and each group's treatment effect does not change over time |
| `"fdTR"` | first difference of the outcome on first difference of the treatment, with period fixed effects | common trends (needs `D0`) |
| `"fdS"` | same as `fdTR` | common trends, and treatment effects stable over time |

Use the `fe` types when you estimated a fixed effects regression, and the `fd` types when you
estimated a first-difference regression. Use an `S` type if you are willing to assume that each
group's treatment effect is the same in every period.

## What the data can look like

* **Several observations per (group, period) cell** are allowed, like individual data inside
  counties. If the treatment or a control varies within a cell, it is replaced by its cell average,
  and the output starts with a note saying so, as in Stata.
* **Unbalanced panels and gaps** are fine.
* **Missing values**: rows with a missing outcome, group, period, treatment or control are dropped,
  following the same rules as Stata.

## Errors

| error | when |
|---|---|
| `TypeError` | `data` is not a pandas DataFrame |
| `KeyError` | a column name is not in `data` |
| `ValueError` | unknown `type`; `fdTR` without `D0`; `other_treatments` with a type other than `feTR`; negative weights; no observations left after dropping missing values |
