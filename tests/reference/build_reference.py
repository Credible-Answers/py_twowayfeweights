"""Regenerate the Stata results that the parity tests compare against.

Usage (from the repository root)::

    python tests/reference/build_reference.py --stata "C:/Program Files/Stata18/StataMP-64.exe"

Requires Stata. The script runs the official command from SSC (``ssc install twowayfeweights``,
installed automatically if missing, together with ``gtools``) on every case in ``tests/cases.py`` and
writes:

* ``tests/fixtures/stata/results.json``: for each case, the Stata command, every number it returns
  and the text it prints;
* ``tests/fixtures/stata/weights/<case>.csv.gz``: the (g,t) weights saved by Stata's ``path()`` option.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from cases import CASES, DATA_DIR, FIXTURES, N_STYR, load, make_sim  # noqa: E402

WORK = HERE / "work"
STATA_FIXTURES = FIXTURES / "stata"


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


def stata_do() -> str:
    q = lambda p: str(p).replace("\\", "/")  # noqa: E731
    out = WORK / "stata_out"
    out.mkdir(exist_ok=True)
    lines = [
        "clear all",
        "set more off",
        "set maxvar 5000",
        f'cd "{q(WORK)}"',
        "cap which gtools",
        "if _rc ssc install gtools",
        "cap which twowayfeweights",
        "if _rc ssc install twowayfeweights",
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


# ---------------------------------------------------------------------- collecting Stata's output


def _num(x):
    """JSON value of a Stata number (missing -> None, whole numbers -> int)."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return None
    x = float(x)
    return int(x) if x.is_integer() and abs(x) < 2**53 else x


def parse_log(text: str) -> dict[str, dict[str, str]]:
    """Command and printed output of each case in Stata's log."""
    text = re.sub(r"\n> ", "", text.replace("\r\n", "\n"))  # undo Stata's 80-column line wrapping
    out = {}
    for block in re.split(r"\nCASE: ", text)[1:]:
        name = block.split("\n", 1)[0].strip()
        m = re.search(r"\n\. cap noi (twowayfeweights[^\n]*)\n(.*?)\n\. if _rc == 0", block, re.S)
        if m is None:
            continue
        command = re.sub(r'\s*path\("[^"]*"\)', "", m.group(1)).strip()
        output = "\n".join(line.rstrip() for line in m.group(2).strip("\n").split("\n"))
        out[name] = {"command": command, "output": output}
    return out


def build_json(out_dir: Path, run_date: str) -> dict:
    """Assemble ``results.json`` from the files the Stata do-file writes."""
    res = pd.read_csv(out_dir / "results.csv", na_values=["."], skipinitialspace=True, float_precision="round_trip").set_index("name")
    ot = pd.read_csv(out_dir / "other_treatments.csv", na_values=["."], skipinitialspace=True, float_precision="round_trip")
    rw = pd.read_csv(out_dir / "random_weights.csv", na_values=["."], skipinitialspace=True, float_precision="round_trip")
    log = parse_log((out_dir / "stata.log").read_text(encoding="utf-8", errors="replace"))
    cases = {}
    for c in CASES:
        n = c["name"]
        r = res.loc[n]
        entry = {**log[n], **{k: _num(r[k]) for k in res.columns}}
        entry["random_weights"] = {
            row["var"]: {k: _num(row[k]) for k in ("Coef", "SE", "tstat", "Correlation")}
            for _, row in rw[rw["name"] == n].iterrows()
        }
        entry["other_treatments"] = [
            {k: _num(row[k]) for k in ("nr_plus", "nr_minus", "sum_plus", "sum_minus", "tot_cells")}
            for _, row in ot[ot["name"] == n].sort_values("j").iterrows()
        ]
        cases[n] = entry
    return {
        "source": "Stata command twowayfeweights from SSC (ssc install twowayfeweights)",
        "generated": run_date,
        "generator": "tests/reference/build_reference.py",
        "cases": cases,
    }


def write_json(data: dict) -> None:
    path = STATA_FIXTURES / "results.json"
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def run_stata(exe: str) -> None:
    do = WORK / "run_stata.do"
    do.write_text(stata_do(), encoding="utf-8")
    subprocess.run([exe, "/e", "do", str(do)], cwd=WORK, check=True)
    out = WORK / "stata_out"
    write_json(build_json(out, date.today().isoformat()))
    (STATA_FIXTURES / "weights").mkdir(parents=True, exist_ok=True)
    for c in CASES:
        p = out / f"{c['name']}_w.dta"
        if p.exists():
            pd.read_stata(p).to_csv(STATA_FIXTURES / "weights" / f"{c['name']}.csv.gz", index=False, float_format="%.9g")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stata", required=True, help="path to the Stata executable")
    ap.add_argument("--gentzkow", type=Path, help="gentzkowetal_didtextbook.dta (first run only)")
    a = ap.parse_args()
    prepare_data(a.gentzkow)
    run_stata(a.stata)
