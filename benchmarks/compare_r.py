"""Compare this package with the R package TwoWayFEWeights (CRAN): run time and results.

This is a comparison only; the test suite checks Python against Stata (tests/test_stata_parity.py).

    python benchmarks/compare_r.py                  # use the stored R results in benchmarks/r_results/
    python benchmarks/compare_r.py --r Rscript      # rerun R first (needs R with TwoWayFEWeights, haven, jsonlite)

Rerunning R needs the Stata data files built by ``tests/reference/build_reference.py``.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))

from cases import CASES, run_python  # noqa: E402

OUT = Path(__file__).resolve().parent / "r_results"
WORK = ROOT / "tests" / "reference" / "work"
COLS = ["beta", "nr_plus", "nr_minus", "sum_plus", "sum_minus", "tot_cells", "sensibility", "sensibility2"]


def r_script() -> str:
    q = lambda p: str(p).replace("\\", "/")  # noqa: E731
    specs = []
    for c in CASES:
        if len(c["controls"]) > 100:  # the R package takes >30 min with 683 controls; skipped
            continue
        specs.append({k: c[k] for k in ("name", "data", "Y", "G", "T", "D", "D0", "type", "controls", "weights",
                                         "other_treatments", "test_random_weights")})
    (WORK / "cases.json").write_text(json.dumps(specs), encoding="utf-8")
    return f"""
suppressPackageStartupMessages({{library(TwoWayFEWeights); library(haven); library(jsonlite)}})
setwd("{q(WORK)}")
cases <- fromJSON("cases.json", simplifyVector = FALSE)
data <- list(wagepan = as.data.frame(read_dta("wagepan.dta")), sim = as.data.frame(read_dta("sim.dta")), simcell = as.data.frame(read_dta("simcell.dta")),
             gentzkow = as.data.frame(read_dta("gentzkow.dta")))
res <- list(); rw <- list()
for (cs in cases) {{
  nz <- function(x) if (length(x) == 0) NULL else unlist(x)
  t0 <- Sys.time()
  r <- tryCatch(twowayfeweights(data[[cs$data]], cs$Y, cs$G, cs$T, cs$D, type = cs$type,
         D0 = if (is.null(cs$D0)) NULL else cs$D0, controls = nz(cs$controls),
         weights = if (is.null(cs$weights)) NULL else data[[cs$data]][[cs$weights]],
         other_treatments = nz(cs$other_treatments), test_random_weights = nz(cs$test_random_weights),
         summary_measures = TRUE), error = function(e) {{ message(cs$name, ": ", conditionMessage(e)); NULL }})
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(r)) {{ res[[length(res) + 1]] <- data.frame(name = cs$name, error = TRUE, seconds = el); next }}
  g <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)
  res[[length(res) + 1]] <- data.frame(name = cs$name, error = FALSE, seconds = el, beta = g(r$beta),
    nr_plus = g(r$nr_plus), nr_minus = g(r$nr_minus), sum_plus = g(r$sum_plus), sum_minus = g(r$sum_minus),
    tot_cells = g(r$tot_cells), sensibility = g(r$sensibility), sensibility2 = g(r$sensibility2))
  if (!is.null(r$mat)) {{
    m <- as.data.frame(r$mat); m$var <- unlist(cs$test_random_weights); m$name <- cs$name
    names(m)[1:4] <- c("Coef", "SE", "tstat", "Correlation"); rw[[length(rw) + 1]] <- m
  }}
}}
dir.create("r_out", showWarnings = FALSE)
write.csv(do.call(rbind, lapply(res, function(d) {{ for (k in c("beta","nr_plus","nr_minus","sum_plus","sum_minus","tot_cells","sensibility","sensibility2")) if (is.null(d[[k]])) d[[k]] <- NA; d }})),
          "r_out/results.csv", row.names = FALSE)
if (length(rw)) write.csv(do.call(rbind, rw), "r_out/random_weights.csv", row.names = FALSE)
"""


def run_r(exe: str) -> None:
    if not (WORK / "gentzkow.dta").exists():
        sys.exit("run tests/reference/build_reference.py --stata ... first (it builds the data files)")
    script = WORK / "run_r.R"
    script.write_text(r_script(), encoding="utf-8")
    subprocess.run([exe, str(script)], cwd=WORK, check=True)
    OUT.mkdir(exist_ok=True)
    for f in ("results.csv", "random_weights.csv"):
        p = WORK / "r_out" / f
        if p.exists():
            pd.read_csv(p).to_csv(OUT / f, index=False, float_format="%.17g")


def _rel(a, b) -> float:
    a = np.nan if a is None else float(a)
    b = float(b)
    if np.isnan(a) and np.isnan(b):
        return 0.0
    if np.isnan(a) or np.isnan(b):
        return np.inf
    return abs(a - b) / max(abs(b), 1e-12)


def compare() -> pd.DataFrame:
    r_res = pd.read_csv(OUT / "results.csv").set_index("name")
    rows = []
    for c in CASES:
        t0 = time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            py = run_python(c).to_dict()
        py_s = time.perf_counter() - t0
        row = {"case": c["name"], "python_s": py_s, "r_s": np.nan, "max_rel_diff": np.nan, "r_status": "not run"}
        if c["name"] in r_res.index:
            r = r_res.loc[c["name"]]
            row["r_s"] = r["seconds"]
            if r["error"]:
                row["r_status"] = "error"
            else:
                row["r_status"] = "ok"
                row["max_rel_diff"] = max(_rel(py[k], r[k]) for k in COLS)
        rows.append(row)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--r", help="path to Rscript (reruns the R package)")
    a = ap.parse_args()
    if a.r:
        run_r(a.r)
    with pd.option_context("display.width", 120, "display.max_rows", None):
        print(compare().to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print("\nmax_rel_diff: largest relative difference between Python and R over beta, the weight counts and sums,")
    print("and the two sensitivity measures. R deviates from Stata on some micro-data cases (missing values).")
