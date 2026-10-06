# Examples

Every example on this page is one of the specifications the test suite compares with Stata.
Each code block is complete: copy it, paste it into Python, and run it. The output shown is
what Python prints; click **Stata output** to see what the Stata command prints. The tests
check that the two match: same text, same counts, same rounded numbers.

!!! note
    In the random-weights table, the last of the eight digits can differ from Stata's.
    Stata stores some intermediate variables as 4-byte floats; Python uses double precision.

## Fixed effects regression (feTR)

The basic case: a regression of log wages on union membership with worker and year fixed effects. `test_random_weights` checks whether the weights are correlated with education.

**Stata**

```stata
twowayfeweights lwage nr year union, type(feTR) summary_measures test_random_weights(educ)
```

**Python**

```python
# pip install twowayfeweights
from twowayfeweights import twowayfeweights, load_wagepan

df = load_wagepan()

res = twowayfeweights(df, Y="lwage", G="nr", T="year", D="union", type="feTR",
                      test_random_weights="educ", summary_measures=True)
print(res)
```

**Output**

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

Regression of variables possibly correlated with the treatment effect on the weights

B[1,4]
             Coef           SE       t-stat  Correlation
educ   -.13445527    .07136021   -1.8841771   -.11825874


The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
```

??? quote "Stata output"

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
    TWFE coefficient (β_fe) =    0.1066
    min σ(Δ) compatible with β_fe and Δ_TR = 0:    0.0969
    min σ(Δ) compatible with treatment effect of opposite sign than β_fe in all (g,t) cells:    3.1759
    Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

    Regression of variables possibly correlated with the treatment effect on the weights

    B[1,4]
                 Coef           SE       t-stat  Correlation
    educ   -.13445527    .07136021    -1.884177   -.11825874


    The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
    ```

## Fixed effects with controls and weights

Same regression, adding two control variables and a weight variable `wt`, and testing two variables.

**Stata**

```stata
twowayfeweights lwage nr year union, type(feTR) summary_measures controls(hours married) weight(wt) test_random_weights(educ exper)
```

**Python**

```python
# pip install twowayfeweights
import numpy as np
from twowayfeweights import twowayfeweights, load_wagepan

df = load_wagepan()
# a weight variable (the same one is used in Stata)
df["wt"] = np.round(0.5 + np.random.default_rng(2020).random(len(df)), 3)

res = twowayfeweights(df, Y="lwage", G="nr", T="year", D="union", type="feTR",
                      controls=["hours", "married"], weights="wt",
                      test_random_weights=["educ", "exper"], summary_measures=True)
print(res)
```

**Output**

```text
Under the common trends assumption,
the TWFE coefficient beta, equal to 0.0979, estimates a weighted sum of 1015 ATTs.
831 ATTs receive a positive weight, and 184 receive a negative weight.
1016 (g,t) cells receive the treatment, but the ATTs of 1 cells receive a weight equal to zero.
------------------------------------------------
Treat. var: union       # ATTs      Σ weights
------------------------------------------------
Positive weights        831         1.0117
Negative weights        184         -0.0117
------------------------------------------------
Total                   1015        1.0000
------------------------------------------------

Summary Measures:
TWFE coefficient (β_fe) = 0.0979
min σ(Δ) compatible with β_fe and Δ_TR = 0: 0.0871
min σ(Δ) compatible with treatment effect of opposite sign than β_fe in all (g,t) cells: 2.5914
Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

Regression of variables possibly correlated with the treatment effect on the weights

B[2,4]
              Coef           SE       t-stat  Correlation
 educ   -.13860738    .06888104    -2.012272   -.12589627
exper   -.18929723     .1028961    -1.839693   -.08125612


The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
```

??? quote "Stata output"

    ```text
    Under the common trends assumption,
    the TWFE coefficient beta, equal to 0.0979, estimates a weighted sum of 1015 ATTs.
    831 ATTs receive a positive weight, and 184 receive a negative weight.
    1016 (g,t) cells receive the treatment, but the ATTs of 1 cells receive a weight equal to zero.
    ------------------------------------------------
    Treat. var: union       # ATTs      Σ weights
    ------------------------------------------------
    Positive weights        831         1.0117
    Negative weights        184         -0.0117
    ------------------------------------------------
    Total                   1015        1.0000
    ------------------------------------------------

    Summary Measures:
    TWFE coefficient (β_fe) =    0.0979
    min σ(Δ) compatible with β_fe and Δ_TR = 0:    0.0871
    min σ(Δ) compatible with treatment effect of opposite sign than β_fe in all (g,t) cells:    2.5914
    Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

    Regression of variables possibly correlated with the treatment effect on the weights

    B[2,4]
                  Coef           SE       t-stat  Correlation
     educ   -.13860738    .06888104   -2.0122719   -.12589627
    exper   -.18929724     .1028961   -1.8396931   -.08125612


    The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
    ```

## Several treatments (other_treatments)

The regression also includes marriage as a second treatment. The output then shows the weights attached to each treatment.

**Stata**

```stata
twowayfeweights lwage nr year union, type(feTR) summary_measures other_treatments(married) test_random_weights(educ)
```

**Python**

```python
# pip install twowayfeweights
from twowayfeweights import twowayfeweights, load_wagepan

df = load_wagepan()

res = twowayfeweights(df, Y="lwage", G="nr", T="year", D="union", type="feTR",
                      other_treatments="married", test_random_weights="educ",
                      summary_measures=True)
print(res)
```

**Output**

```text
Under the common trends assumption,
the TWFE coefficient beta, equal to 0.1038, estimates the sum of several terms.

The first term is a weighted sum of 1016 ATTs of the treatment.
852 ATTs receive a positive weight, and 164 receive a negative weight.
------------------------------------------------
Treat. var: union       # ATTs      Σ weights
------------------------------------------------
Positive weights        852         1.0117
Negative weights        164         -0.0117
------------------------------------------------
Total                   1016        1.0000
------------------------------------------------

The next term is a weighted sum of 1914 ATTs of treatment 1 included in the other_treatments option.
917 ATTs receive a positive weight, and 997 receive a negative weight.
------------------------------------------------
Other treat.: married   # ATTs      Σ weights
------------------------------------------------
Positive weights        917         0.4828
Negative weights        997         -0.4828
------------------------------------------------
Total                   1914        0.0000
------------------------------------------------

Regression of variables possibly correlated with the treatment effect on the weights attached to the treatment

B[1,4]
             Coef           SE       t-stat  Correlation
educ   -.13467697    .07124164   -1.8904249    -.1185583


The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
```

??? quote "Stata output"

    ```text
    Under the common trends assumption,
    the TWFE coefficient beta, equal to 0.1038, estimates the sum of several terms.

    The first term is a weighted sum of 1016 ATTs of the treatment.
    852 ATTs receive a positive weight, and 164 receive a negative weight.
    ------------------------------------------------
    Treat. var: union       # ATTs      Σ weights
    ------------------------------------------------
    Positive weights        852         1.0117
    Negative weights        164         -0.0117
    ------------------------------------------------
    Total                   1016        1.0000
    ------------------------------------------------

    The next term is a weighted sum of 1914 ATTs of treatment 1 included in the other_treatments option.
    917 ATTs receive a positive weight, and 997 receive a negative weight.
    ------------------------------------------------
    Other treat.: married   # ATTs      Σ weights
    ------------------------------------------------
    Positive weights        917         0.4828
    Negative weights        997         -0.4828
    ------------------------------------------------
    Total                   1914        0.0000
    ------------------------------------------------

    Regression of variables possibly correlated with the treatment effect on the weights attached to the treatment

    B[1,4]
                 Coef           SE       t-stat  Correlation
    educ   -.13467696    .07124164   -1.8904248    -.1185583


    The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
    ```

## Several treatments with controls and weights

Two other treatments, one control and weights.

**Stata**

```stata
twowayfeweights lwage nr year union, type(feTR) summary_measures controls(hours) weight(wt) other_treatments(married south) test_random_weights(educ)
```

**Python**

```python
# pip install twowayfeweights
import numpy as np
from twowayfeweights import twowayfeweights, load_wagepan

df = load_wagepan()
# a weight variable (the same one is used in Stata)
df["wt"] = np.round(0.5 + np.random.default_rng(2020).random(len(df)), 3)

res = twowayfeweights(df, Y="lwage", G="nr", T="year", D="union", type="feTR", controls="hours",
                      weights="wt", other_treatments=["married", "south"],
                      test_random_weights="educ", summary_measures=True)
print(res)
```

**Output**

```text
Under the common trends assumption,
the TWFE coefficient beta, equal to 0.0990, estimates the sum of several terms.

The first term is a weighted sum of 1016 ATTs of the treatment.
830 ATTs receive a positive weight, and 186 receive a negative weight.
------------------------------------------------
Treat. var: union       # ATTs      Σ weights
------------------------------------------------
Positive weights        830         1.0117
Negative weights        186         -0.0117
------------------------------------------------
Total                   1016        1.0000
------------------------------------------------

The next term is a weighted sum of 1914 ATTs of treatment 1 included in the other_treatments option.
862 ATTs receive a positive weight, and 1052 receive a negative weight.
------------------------------------------------
Other treat.: married   # ATTs      Σ weights
------------------------------------------------
Positive weights        862         0.4873
Negative weights        1052        -0.4873
------------------------------------------------
Total                   1914        0.0000
------------------------------------------------

The next term is a weighted sum of 1529 ATTs of treatment 2 included in the other_treatments option.
768 ATTs receive a positive weight, and 761 receive a negative weight.
------------------------------------------------
Other treat.: south     # ATTs      Σ weights
------------------------------------------------
Positive weights        768         0.3698
Negative weights        761         -0.3698
------------------------------------------------
Total                   1529        0.0000
------------------------------------------------

Regression of variables possibly correlated with the treatment effect on the weights attached to the treatment

B[1,4]
             Coef           SE       t-stat  Correlation
educ   -.13866343    .06899447   -2.0097759   -.12594082


The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
```

??? quote "Stata output"

    ```text
    Under the common trends assumption,
    the TWFE coefficient beta, equal to 0.0990, estimates the sum of several terms.

    The first term is a weighted sum of 1016 ATTs of the treatment.
    830 ATTs receive a positive weight, and 186 receive a negative weight.
    ------------------------------------------------
    Treat. var: union       # ATTs      Σ weights
    ------------------------------------------------
    Positive weights        830         1.0117
    Negative weights        186         -0.0117
    ------------------------------------------------
    Total                   1016        1.0000
    ------------------------------------------------

    The next term is a weighted sum of 1914 ATTs of treatment 1 included in the other_treatments option.
    862 ATTs receive a positive weight, and 1052 receive a negative weight.
    ------------------------------------------------
    Other treat.: married   # ATTs      Σ weights
    ------------------------------------------------
    Positive weights        862         0.4873
    Negative weights        1052        -0.4873
    ------------------------------------------------
    Total                   1914        0.0000
    ------------------------------------------------

    The next term is a weighted sum of 1529 ATTs of treatment 2 included in the other_treatments option.
    768 ATTs receive a positive weight, and 761 receive a negative weight.
    ------------------------------------------------
    Other treat.: south     # ATTs      Σ weights
    ------------------------------------------------
    Positive weights        768         0.3698
    Negative weights        761         -0.3698
    ------------------------------------------------
    Total                   1529        0.0000
    ------------------------------------------------

    Regression of variables possibly correlated with the treatment effect on the weights attached to the treatment

    B[1,4]
                 Coef           SE       t-stat  Correlation
    educ   -.13866343    .06899447   -2.0097759   -.12594082


    The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
    ```

## Fixed effects, stable treatment effects (feS)

Same regression as the first example, under the extra assumption that each group's treatment effect does not change over time.

**Stata**

```stata
twowayfeweights lwage nr year union, type(feS) summary_measures test_random_weights(educ)
```

**Python**

```python
# pip install twowayfeweights
from twowayfeweights import twowayfeweights, load_wagepan

df = load_wagepan()

res = twowayfeweights(df, Y="lwage", G="nr", T="year", D="union", type="feS",
                      test_random_weights="educ", summary_measures=True)
print(res)
```

**Output**

```text
Under the common trends assumption and
the assumption that groups' treatment effects do not change over time,
the TWFE coefficient beta, equal to 0.1066, estimates a weighted sum of 228 ATTs.
228 ATTs receive a positive weight, and 0 receive a negative weight.
------------------------------------------------
Treat. var: union       # ATTs      Σ weights
------------------------------------------------
Positive weights        228         1.0000
Negative weights        0           0.0000
------------------------------------------------
Total                   228         1.0000
------------------------------------------------

Summary Measures:
TWFE coefficient (β_fe) = 0.1066
min σ(Δ) compatible with β_fe and Δ_TR = 0: 0.1555
Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

Regression of variables possibly correlated with the treatment effect on the weights

B[1,4]
             Coef           SE       t-stat  Correlation
educ    .11004249    .14939977    .73656399    .04569715


The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
```

??? quote "Stata output"

    ```text
    Under the common trends assumption and
    the assumption that groups' treatment effects do not change over time,
    the TWFE coefficient beta, equal to 0.1066, estimates a weighted sum of 228 ATTs.
    228 ATTs receive a positive weight, and 0 receive a negative weight.
    ------------------------------------------------
    Treat. var: union       # ATTs      Σ weights
    ------------------------------------------------
    Positive weights        228         1.0000
    Negative weights        0           0.0000
    ------------------------------------------------
    Total                   228         1.0000
    ------------------------------------------------

    Summary Measures:
    TWFE coefficient (β_fe) =    0.1066
    min σ(Δ) compatible with β_fe and Δ_TR = 0:    0.1555
    Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

    Regression of variables possibly correlated with the treatment effect on the weights

    B[1,4]
                 Coef           SE       t-stat  Correlation
    educ    .11004248    .14939977    .73656395    .04569714


    The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
    ```

## feS with controls and weights

The feS version of the second example.

**Stata**

```stata
twowayfeweights lwage nr year union, type(feS) summary_measures controls(hours married) weight(wt) test_random_weights(educ)
```

**Python**

```python
# pip install twowayfeweights
import numpy as np
from twowayfeweights import twowayfeweights, load_wagepan

df = load_wagepan()
# a weight variable (the same one is used in Stata)
df["wt"] = np.round(0.5 + np.random.default_rng(2020).random(len(df)), 3)

res = twowayfeweights(df, Y="lwage", G="nr", T="year", D="union", type="feS",
                      controls=["hours", "married"], weights="wt", test_random_weights="educ",
                      summary_measures=True)
print(res)
```

**Output**

```text
Under the common trends assumption and
the assumption that groups' treatment effects do not change over time,
the TWFE coefficient beta, equal to 0.0979, estimates a weighted sum of 228 ATTs.
228 ATTs receive a positive weight, and 0 receive a negative weight.
------------------------------------------------
Treat. var: union       # ATTs      Σ weights
------------------------------------------------
Positive weights        228         1.0000
Negative weights        0           0.0000
------------------------------------------------
Total                   228         1.0000
------------------------------------------------

Summary Measures:
TWFE coefficient (β_fe) = 0.0979
min σ(Δ) compatible with β_fe and Δ_TR = 0: 0.1143
Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

Regression of variables possibly correlated with the treatment effect on the weights

B[1,4]
             Coef           SE       t-stat  Correlation
educ    .12252658    .10384546    1.1798934    .06419958


The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
```

??? quote "Stata output"

    ```text
    Under the common trends assumption and
    the assumption that groups' treatment effects do not change over time,
    the TWFE coefficient beta, equal to 0.0979, estimates a weighted sum of 228 ATTs.
    228 ATTs receive a positive weight, and 0 receive a negative weight.
    ------------------------------------------------
    Treat. var: union       # ATTs      Σ weights
    ------------------------------------------------
    Positive weights        228         1.0000
    Negative weights        0           0.0000
    ------------------------------------------------
    Total                   228         1.0000
    ------------------------------------------------

    Summary Measures:
    TWFE coefficient (β_fe) =    0.0979
    min σ(Δ) compatible with β_fe and Δ_TR = 0:    0.1143
    Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

    Regression of variables possibly correlated with the treatment effect on the weights

    B[1,4]
                 Coef           SE       t-stat  Correlation
    educ    .12252659    .10384546    1.1798936    .06419959


    The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
    ```

## First-difference regression (fdTR)

A regression of the change in log wages on the change in union membership. `D0` is the treatment level (not differenced).

**Stata**

```stata
twowayfeweights diff_lwage nr year diff_union union, type(fdTR) summary_measures test_random_weights(educ)
```

**Python**

```python
# pip install twowayfeweights
from twowayfeweights import twowayfeweights, load_wagepan

df = load_wagepan()

res = twowayfeweights(df, Y="diff_lwage", G="nr", T="year", D="diff_union", type="fdTR",
                      D0="union", test_random_weights="educ", summary_measures=True)
print(res)
```

**Output**

```text
Under the common trends assumption,
the TWFE coefficient beta, equal to 0.0601, estimates a weighted sum of 1016 ATTs.
611 ATTs receive a positive weight, and 405 receive a negative weight.
------------------------------------------------
Treat. var: diff_union  # ATTs      Σ weights
------------------------------------------------
Positive weights        611         1.0476
Negative weights        405         -0.0476
------------------------------------------------
Total                   1016        1.0000
------------------------------------------------

Summary Measures:
TWFE coefficient (β_fd) = 0.0601
min σ(Δ) compatible with β_fd and Δ_TR = 0: 0.0321
min σ(Δ) compatible with treatment effect of opposite sign than β_fd in all (g,t) cells: 0.5799
Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

Regression of variables possibly correlated with the treatment effect on the weights

B[1,4]
             Coef           SE       t-stat  Correlation
educ   -.06649017    .02893837   -2.2976472    -.0994798


The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
```

??? quote "Stata output"

    ```text
    Under the common trends assumption,
    the TWFE coefficient beta, equal to 0.0601, estimates a weighted sum of 1016 ATTs.
    611 ATTs receive a positive weight, and 405 receive a negative weight.
    ------------------------------------------------
    Treat. var: diff_union  # ATTs      Σ weights
    ------------------------------------------------
    Positive weights        611         1.0476
    Negative weights        405         -0.0476
    ------------------------------------------------
    Total                   1016        1.0000
    ------------------------------------------------

    Summary Measures:
    TWFE coefficient (β_fd) =    0.0601
    min σ(Δ) compatible with β_fd and Δ_TR = 0:    0.0321
    min σ(Δ) compatible with treatment effect of opposite sign than β_fd in all (g,t) cells:    0.5799
    Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

    Regression of variables possibly correlated with the treatment effect on the weights

    B[1,4]
                 Coef           SE       t-stat  Correlation
    educ   -.06649017    .02893837   -2.2976473    -.0994798


    The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
    ```

## fdTR with a control and weights

The first-difference regression with a control and weights.

**Stata**

```stata
twowayfeweights diff_lwage nr year diff_union union, type(fdTR) summary_measures controls(hours) weight(wt) test_random_weights(educ)
```

**Python**

```python
# pip install twowayfeweights
import numpy as np
from twowayfeweights import twowayfeweights, load_wagepan

df = load_wagepan()
# a weight variable (the same one is used in Stata)
df["wt"] = np.round(0.5 + np.random.default_rng(2020).random(len(df)), 3)

res = twowayfeweights(df, Y="diff_lwage", G="nr", T="year", D="diff_union", type="fdTR",
                      D0="union", controls="hours", weights="wt", test_random_weights="educ",
                      summary_measures=True)
print(res)
```

**Output**

```text
Under the common trends assumption,
the TWFE coefficient beta, equal to 0.0776, estimates a weighted sum of 1016 ATTs.
612 ATTs receive a positive weight, and 404 receive a negative weight.
------------------------------------------------
Treat. var: diff_union  # ATTs      Σ weights
------------------------------------------------
Positive weights        612         1.0480
Negative weights        404         -0.0480
------------------------------------------------
Total                   1016        1.0000
------------------------------------------------

Summary Measures:
TWFE coefficient (β_fd) = 0.0776
min σ(Δ) compatible with β_fd and Δ_TR = 0: 0.0397
min σ(Δ) compatible with treatment effect of opposite sign than β_fd in all (g,t) cells: 0.7201
Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

Regression of variables possibly correlated with the treatment effect on the weights

B[1,4]
             Coef           SE       t-stat  Correlation
educ   -.06724739    .02768205   -2.4292779   -.10634003


The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
```

??? quote "Stata output"

    ```text
    Under the common trends assumption,
    the TWFE coefficient beta, equal to 0.0776, estimates a weighted sum of 1016 ATTs.
    612 ATTs receive a positive weight, and 404 receive a negative weight.
    ------------------------------------------------
    Treat. var: diff_union  # ATTs      Σ weights
    ------------------------------------------------
    Positive weights        612         1.0480
    Negative weights        404         -0.0480
    ------------------------------------------------
    Total                   1016        1.0000
    ------------------------------------------------

    Summary Measures:
    TWFE coefficient (β_fd) =    0.0776
    min σ(Δ) compatible with β_fd and Δ_TR = 0:    0.0397
    min σ(Δ) compatible with treatment effect of opposite sign than β_fd in all (g,t) cells:    0.7201
    Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

    Regression of variables possibly correlated with the treatment effect on the weights

    B[1,4]
                 Coef           SE       t-stat  Correlation
    educ   -.06724739    .02768205   -2.4292779   -.10634003


    The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
    ```

## First differences, stable treatment effects (fdS)

First-difference regression under the extra assumption that treatment effects are stable over time.

**Stata**

```stata
twowayfeweights diff_lwage nr year diff_union, type(fdS) summary_measures test_random_weights(educ)
```

**Python**

```python
# pip install twowayfeweights
from twowayfeweights import twowayfeweights, load_wagepan

df = load_wagepan()

res = twowayfeweights(df, Y="diff_lwage", G="nr", T="year", D="diff_union", type="fdS",
                      test_random_weights="educ", summary_measures=True)
print(res)
```

**Output**

```text
Under the common trends assumption and
the assumption that groups' treatment effects do not change over time,
the TWFE coefficient beta, equal to 0.0601, estimates a weighted sum of 228 ATTs.
228 ATTs receive a positive weight, and 0 receive a negative weight.
------------------------------------------------
Treat. var: diff_union  # ATTs      Σ weights
------------------------------------------------
Positive weights        228         1.0000
Negative weights        0           0.0000
------------------------------------------------
Total                   228         1.0000
------------------------------------------------

Summary Measures:
TWFE coefficient (β_fd) = 0.0601
min σ(Δ) compatible with β_fd and Δ_TR = 0: 2.2081
Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

Regression of variables possibly correlated with the treatment effect on the weights

B[1,4]
             Coef           SE       t-stat  Correlation
educ    7.0513352    4.4381072    1.5888159    .11625131


The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
```

??? quote "Stata output"

    ```text
    Under the common trends assumption and
    the assumption that groups' treatment effects do not change over time,
    the TWFE coefficient beta, equal to 0.0601, estimates a weighted sum of 228 ATTs.
    228 ATTs receive a positive weight, and 0 receive a negative weight.
    ------------------------------------------------
    Treat. var: diff_union  # ATTs      Σ weights
    ------------------------------------------------
    Positive weights        228         1.0000
    Negative weights        0           0.0000
    ------------------------------------------------
    Total                   228         1.0000
    ------------------------------------------------

    Summary Measures:
    TWFE coefficient (β_fd) =    0.0601
    min σ(Δ) compatible with β_fd and Δ_TR = 0:    2.2081
    Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

    Regression of variables possibly correlated with the treatment effect on the weights

    B[1,4]
                 Coef           SE       t-stat  Correlation
    educ    7.0513258    4.4381053    1.5888144    .11625119


    The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
    ```

## fdS with a control and weights

The fdS version with a control and weights.

**Stata**

```stata
twowayfeweights diff_lwage nr year diff_union, type(fdS) summary_measures controls(hours) weight(wt) test_random_weights(educ)
```

**Python**

```python
# pip install twowayfeweights
import numpy as np
from twowayfeweights import twowayfeweights, load_wagepan

df = load_wagepan()
# a weight variable (the same one is used in Stata)
df["wt"] = np.round(0.5 + np.random.default_rng(2020).random(len(df)), 3)

res = twowayfeweights(df, Y="diff_lwage", G="nr", T="year", D="diff_union", type="fdS",
                      controls="hours", weights="wt", test_random_weights="educ",
                      summary_measures=True)
print(res)
```

**Output**

```text
Under the common trends assumption and
the assumption that groups' treatment effects do not change over time,
the TWFE coefficient beta, equal to 0.0776, estimates a weighted sum of 228 ATTs.
228 ATTs receive a positive weight, and 0 receive a negative weight.
------------------------------------------------
Treat. var: diff_union  # ATTs      Σ weights
------------------------------------------------
Positive weights        228         1.0000
Negative weights        0           0.0000
------------------------------------------------
Total                   228         1.0000
------------------------------------------------

Summary Measures:
TWFE coefficient (β_fd) = 0.0776
min σ(Δ) compatible with β_fd and Δ_TR = 0: 2.7278
Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

Regression of variables possibly correlated with the treatment effect on the weights

B[1,4]
             Coef           SE       t-stat  Correlation
educ    7.5566584    4.4677004    1.6913978    .13146924


The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
```

??? quote "Stata output"

    ```text
    Under the common trends assumption and
    the assumption that groups' treatment effects do not change over time,
    the TWFE coefficient beta, equal to 0.0776, estimates a weighted sum of 228 ATTs.
    228 ATTs receive a positive weight, and 0 receive a negative weight.
    ------------------------------------------------
    Treat. var: diff_union  # ATTs      Σ weights
    ------------------------------------------------
    Positive weights        228         1.0000
    Negative weights        0           0.0000
    ------------------------------------------------
    Total                   228         1.0000
    ------------------------------------------------

    Summary Measures:
    TWFE coefficient (β_fd) =    0.0776
    min σ(Δ) compatible with β_fd and Δ_TR = 0:    2.7278
    Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

    Regression of variables possibly correlated with the treatment effect on the weights

    B[1,4]
                 Coef           SE       t-stat  Correlation
    educ    7.5566602    4.4677008     1.691398    .13146928


    The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
    ```

## Newspapers and turnout, 683 controls

Chapter 5 of the DiD book, with data from Gentzkow, Shapiro and Sinkinson (2011): the effect of the number of daily newspapers in a county on presidential turnout, with 683 state-year dummies as controls.

**Stata**

```stata
twowayfeweights changeprestout cnty90 year changedailies numdailies, type(fdTR) summary_measures controls(styr1-styr683) test_random_weights(year)
```

**Python**

```python
# pip install twowayfeweights
import pandas as pd
from twowayfeweights import twowayfeweights

df = pd.read_csv("https://raw.githubusercontent.com/Credible-Answers/py_twowayfeweights/main/tests/data/gentzkow.csv.gz")
# one dummy per state-year (683 dummies), used as controls
styr = pd.get_dummies(df["styr"], prefix="styr", dtype=float)
df = pd.concat([df, styr], axis=1)

res = twowayfeweights(df, Y="changeprestout", G="cnty90", T="year", D="changedailies", type="fdTR",
                      D0="numdailies", controls=list(styr.columns), test_random_weights="year",
                      summary_measures=True)
print(res)
```

**Output**

```text
Under the common trends assumption,
the TWFE coefficient beta, equal to 0.0026, estimates a weighted sum of 9876 ATTs.
5371 ATTs receive a positive weight, and 4505 receive a negative weight.
10378 (g,t) cells receive the treatment, but the ATTs of 502 cells receive a weight equal to zero.
------------------------------------------------
Treat. var: changedail~s# ATTs      Σ weights
------------------------------------------------
Positive weights        5371        2.4271
Negative weights        4505        -1.4271
------------------------------------------------
Total                   9876        1.0000
------------------------------------------------

Summary Measures:
TWFE coefficient (β_fd) = 0.0026
min σ(Δ) compatible with β_fd and Δ_TR = 0: 0.0004
min σ(Δ) compatible with treatment effect of opposite sign than β_fd in all (g,t) cells: 0.0006
Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

Regression of variables possibly correlated with the treatment effect on the weights

B[1,4]
             Coef           SE       t-stat  Correlation
year    -.1674271    .05101171   -3.2821306    -.0631614


The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
```

??? quote "Stata output"

    ```text
    Under the common trends assumption,
    the TWFE coefficient beta, equal to 0.0026, estimates a weighted sum of 9876 ATTs.
    5371 ATTs receive a positive weight, and 4505 receive a negative weight.
    10378 (g,t) cells receive the treatment, but the ATTs of 502 cells receive a weight equal to zero.
    ------------------------------------------------
    Treat. var: changedail~s# ATTs      Σ weights
    ------------------------------------------------
    Positive weights        5371        2.4271
    Negative weights        4505        -1.4271
    ------------------------------------------------
    Total                   9876        1.0000
    ------------------------------------------------

    Summary Measures:
    TWFE coefficient (β_fd) =    0.0026
    min σ(Δ) compatible with β_fd and Δ_TR = 0:    0.0004
    min σ(Δ) compatible with treatment effect of opposite sign than β_fd in all (g,t) cells:    0.0006
    Reference: Corollary 1, de Chaisemartin, C and D'Haultfoeuille, X (2020a)

    Regression of variables possibly correlated with the treatment effect on the weights

    B[1,4]
                 Coef           SE       t-stat  Correlation
    year    -.1674271    .05101171   -3.2821307    -.0631614


    The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).
    ```
