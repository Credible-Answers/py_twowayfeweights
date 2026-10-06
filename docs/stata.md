# Comparison with Stata

This package is tested against the Stata command `twowayfeweights` from SSC
(`ssc install twowayfeweights`). Every change to the package must keep these tests passing.

## What is compared

The test suite runs **36 specifications** in both Stata and Python: all four types, with and
without controls, weights, other treatments and `test_random_weights`, on three datasets:

* the `wagepan` data (the example data shipped with the package);
* the Gentzkow, Shapiro and Sinkinson (2011) newspapers data, including a regression with 683 controls;
* a simulated unbalanced panel with gaps, several observations per cell and missing values.

For each specification, the tests check that:

| what | how close |
|---|---|
| numbers of positive, negative and treated-cell weights | identical |
| beta | within 1e-7 relative |
| sums of weights, both summary measures | within 1e-6 relative |
| `test_random_weights` table (Coef, SE, t-stat, Correlation) | within 5e-6 relative |
| weight of every (group, period) cell | within 2e-6 of the largest weight |

The [examples](examples.md) are 11 of these specifications, written as copy-paste code.

## Why the numbers are not identical to the last digit

Stata stores intermediate variables created with `gen`, `predict` and `egen` as 4-byte floats,
which keep about 7 significant digits. This package computes everything in double precision. The
differences that remain are this rounding in Stata; they never show in the rounded numbers of the
report, except sometimes in the 8th digit of the random-weights table.

## One known difference

For `type="feS"` with **several observations per (group, period) cell**, Stata builds a running
average one observation at a time and then keeps an arbitrary observation of each cell. Its result
depends on the order of the rows, and changes from run to run because Stata breaks sort ties
randomly. This package uses the formula of the paper instead. With one observation per cell, the
usual case, the results are identical. For that one specification, the tests compare beta only.

## Where the Stata results are stored

* `tests/fixtures/stata/results.json`: for each specification, the Stata command, every number it
  returns and the text it prints;
* `tests/fixtures/stata/weights/`: the cell weights saved by Stata's `path()` option.

`tests/reference/build_reference.py` reruns Stata and regenerates these files, e.g. after a new
release of the Stata command.

## R package

`benchmarks/compare_r.py` compares the results and run times with the R package `TwoWayFEWeights`.
This is for information only; the tests compare with Stata.
