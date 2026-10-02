# Changelog

## 0.1.0

First release on PyPI.

* Port of the Stata command `twowayfeweights`, supporting all four estimation types: `feTR`, `feS`,
  `fdTR` and `fdS`.
* Options: `controls`, `weights`, `other_treatments`, `test_random_weights`, `summary_measures` and `path`.
* Checked against Stata on 36 specifications, 17 of them `feS`/`fdS` (see README, "Parity with Stata").
* Depends only on `numpy`, `pandas` and `scipy`. Regressions run at the (g,t) cell level.
* Result object with Stata-style printing. Legacy dictionary access and `print_twowayfeweights` are kept.
