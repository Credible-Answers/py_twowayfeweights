"""Regenerate the Stata (and R) reference fixtures used by the parity tests.

Usage (from the repository root)::

    python tests/reference/build_reference.py --stata "C:/Program Files/Stata18/StataMP-64.exe"
    python tests/reference/build_reference.py --r Rscript

Requires Stata with ``gtools`` installed, and/or R with the CRAN package ``TwoWayFEWeights``.
The Stata program that is run is ``tests/reference/twowayfeweights.ado``.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))

from cases import CASES, DATA_DIR, FIXTURES, N_STYR, load, make_sim  # noqa: E402

WORK = HERE / "work"
ADO = HERE / "twowayfeweights.ado"  # Stata reference implementation (chaisemartinPackages/twowayfeweights)


def prepare_data(gentzkow_source: Path | None) -> None:
    DATA_DIR.mkdir(exist_ok=True)
    WORK.mkdir(exist_ok=True)
    if not (DATA_DIR / "gentzkow.csv.gz").exists():
        if gentzkow_source is None:
            raise SystemExit("tests/data/gentzkow.csv.gz missing: pass --gentzkow path/to/gentzkowetal_didtextbook.dta")
        g = pd.read_stata(gentzkow_source)
        cols = ["cnty90", "year", "prestout", "numdailies", "styr", "changedailies", "changeprestout", "share_urb_baseline"]
        g[cols].to_csv(DATA_DIR / "gentzkow.csv.gz", index=False, float_format="%.9g", compression="gzip")
    if not (DATA_DIR / "sim.csv.gz").exists():
        make_sim().to_csv(DATA_DIR / "sim.csv.gz", index=False, float_format="%.9g", compression="gzip")

    for name in ("wagepan", "sim", "simcell"):
        load(name).to_stata(WORK / f"{name}.dta", write_index=False, version=118)
    g = pd.read_csv(DATA_DIR / "gentzkow.csv.gz")
    g.to_stata(WORK / "gentzkow_raw.dta", write_index=False, version=118)


# ---------------------------------------------------------------------- Stata


def stata_do() -> str:
    q = lambda p: str(p).replace("\\", "/")  # noqa: E731
    out = WORK / "stata_out"
    out.mkdir(exist_ok=True)
    lines = [
        "clear all",
        "set more off",
        "set maxvar 5000",
        f'cd "{q(WORK)}"',
        f'qui do "{q(ADO)}"',
        # Gentzkow data with state-year dummies, built exactly as in the did_book chapter.
        'use "gentzkow_raw.dta", clear',
        "qui tab styr, gen(styr)",
        f"assert r(r) == {N_STYR}",
        'save "gentzkow.dta", replace',
        'log using "stata_out/stata.log", text replace',
        'file open fh using "stata_out/results.csv", write replace',
        'file write fh "name,beta,nr_plus,nr_minus,sum_plus,sum_minus,tot_cells,sensibility,sensibility2" _n',
        'file open fo using "stata_out/other_treatments.csv", write replace',
        'file write fo "name,j,nr_plus,nr_minus,sum_plus,sum_minus,tot_cells" _n',
        'file open fr using "stata_out/random_weights.csv", write replace',
        'file write fr "name,var,Coef,SE,tstat,Correlation" _n',
    ]
    g = "%24.17g"
    for c in CASES:
        n = c["name"]
        vl = " ".join(v for v in [c["Y"], c["G"], c["T"], c["D"], c["D0"]] if v)
        opts = [f"type({c['type']})", "summary_measures", f'path("stata_out/{n}_w.dta")']
        if c["controls"]:
            ctrl = c["controls"]
            opts.append(f"controls({'styr1-styr%d' % N_STYR if ctrl[0] == 'styr1' and len(ctrl) == N_STYR else ' '.join(ctrl)})")
        if c["weights"]:
            opts.append(f"weight({c['weights']})")
        if c["other_treatments"]:
            opts.append(f"other_treatments({' '.join(c['other_treatments'])})")
        if c["test_random_weights"]:
            opts.append(f"test_random_weights({' '.join(c['test_random_weights'])})")
        M = "e(M1)" if c["other_treatments"] else "e(M)"
        lines += [
            "",
            f"* ---- {n}",
            f'use "{c["data"]}.dta", clear',
            "scalar drop _all",
            f'di "CASE: {n}"',
            f"cap noi twowayfeweights {vl}, {' '.join(opts)}",
            "if _rc == 0 {",
            f'  file write fh "{n}," {g} (e(beta)) "," {g} (el({M},1,1)) "," {g} (el({M},2,1)) "," '
            f'{g} (el({M},1,2)) "," {g} (el({M},2,2)) "," {g} (tot_cells) "," {g} (e(lb_se_te)) "," {g} (e(lb_se_te2)) _n',
        ]
        for j, _ in enumerate(c["other_treatments"], start=1):
            lines.append(
                f'  file write fo "{n},{j}," {g} (el(e(M{j + 1}),1,1)) "," {g} (el(e(M{j + 1}),2,1)) "," '
                f'{g} (el(e(M{j + 1}),1,2)) "," {g} (el(e(M{j + 1}),2,2)) "," {g} (tot_cells_{j}) _n'
            )
        for i, v in enumerate(c["test_random_weights"], start=1):
            B = "e(randomweightstest1)"
            lines.append(
                f'  file write fr "{n},{v}," {g} (el({B},{i},1)) "," {g} (el({B},{i},2)) "," '
                f'{g} (el({B},{i},3)) "," {g} (el({B},{i},4)) _n'
            )
        lines += ["}", "else {", f'  file write fh "{n},ERROR" _n', "}"]
    lines += ["", "file close fh", "file close fo", "file close fr", "log close", "exit, clear STATA"]
    return "\n".join(lines) + "\n"


def run_stata(exe: str) -> None:
    do = WORK / "run_stata.do"
    do.write_text(stata_do(), encoding="utf-8")
    subprocess.run([exe, "/e", "do", str(do)], cwd=WORK, check=True)
    out = WORK / "stata_out"
    dest = FIXTURES / "stata"
    (dest / "weights").mkdir(parents=True, exist_ok=True)
    for f in ("results.csv", "other_treatments.csv", "random_weights.csv"):
        pd.read_csv(out / f, na_values=["."], skipinitialspace=True).to_csv(dest / f, index=False, float_format="%.17g")
    (dest / "stata.log").write_text((out / "stata.log").read_text(encoding="utf-8", errors="replace"), encoding="utf-8")
    for c in CASES:
        p = out / f"{c['name']}_w.dta"
        if p.exists():
            pd.read_stata(p).to_csv(dest / "weights" / f"{c['name']}.csv.gz", index=False, float_format="%.9g")


# ---------------------------------------------------------------------- R


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
res <- list(); rw <- list(); ot <- list()
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
    script = WORK / "run_r.R"
    script.write_text(r_script(), encoding="utf-8")
    subprocess.run([exe, str(script)], cwd=WORK, check=True)
    dest = FIXTURES / "r"
    dest.mkdir(parents=True, exist_ok=True)
    for f in ("results.csv", "random_weights.csv"):
        p = WORK / "r_out" / f
        if p.exists():
            pd.read_csv(p).to_csv(dest / f, index=False, float_format="%.17g")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stata", help="path to the Stata executable")
    ap.add_argument("--r", help="path to Rscript")
    ap.add_argument("--gentzkow", type=Path, help="gentzkowetal_didtextbook.dta (first run only)")
    a = ap.parse_args()
    prepare_data(a.gentzkow)
    if a.stata:
        run_stata(a.stata)
    if a.r:
        (WORK / "gentzkow.dta").exists() or sys.exit("run --stata first (it builds work/gentzkow.dta)")
        run_r(a.r)
