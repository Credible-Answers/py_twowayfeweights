# twowayfeweights

**twowayfeweights** tells you whether you can trust a two-way fixed effects (TWFE) regression when
treatment effects differ across groups or over time.

It is the Python version of the Stata and R command `twowayfeweights`, by de Chaisemartin and
D'Haultfoeuille (2020), and it prints the same results as Stata.

## The idea in plain words

You have a panel: groups (people, counties, firms...) observed over several periods, some of them
treated. You run the usual regression of the outcome on the treatment with group and period fixed
effects, and you get a coefficient, beta.

Beta is a **weighted average of the treatment effects** of every treated (group, period) cell. The
problem: some of those weights can be **negative**. When they are, beta can be positive even if the
treatment effect is negative in every cell, or the other way around.

`twowayfeweights` computes those weights and tells you:

* how many cells get a positive weight and how many get a negative one, and how much each side sums to;
* how much the treatment effects would need to vary across cells for beta to be misleading
  (the *summary measures*): a small number means beta is fragile;
* optionally, whether the weights are correlated with variables that may drive the treatment effect.

## Install

```bash
pip install twowayfeweights
```

It needs Python 3.9 or newer, and only installs `numpy`, `pandas` and `scipy`.

## A first example

```python
# pip install twowayfeweights
from twowayfeweights import twowayfeweights, load_wagepan

df = load_wagepan()          # example data: 545 workers observed from 1980 to 1987
res = twowayfeweights(df, Y="lwage", G="nr", T="year", D="union", summary_measures=True)
print(res)                   # Stata-style report
print(res.beta)              # every number is also an attribute: see "Outputs"
```

```text
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


The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
0.10662746555843541
```

**How to read it.** Beta (0.1066) averages the effects of 967 treated cells. 147 of them get a
negative weight, but those weights only sum to -0.0105, so they matter little. The first summary
measure says beta could be produced by treatment effects that are zero on average if their
standard deviation across cells were only 0.0969. That is small compared with beta: the result is
sensitive to heterogeneous effects.

## Where to go next

* [Arguments](arguments.md): every input, with its Stata equivalent.
* [Outputs](outputs.md): where each number is stored and how to get it.
* [Examples](examples.md): copy-paste examples, each checked against Stata.
* [Comparison with Stata](stata.md): how the package is tested.

## Reference

de Chaisemartin, C. and D'Haultfoeuille, X. (2020). Two-Way Fixed Effects Estimators with
Heterogeneous Treatment Effects. *American Economic Review*, 110(9), 2964-2996.
