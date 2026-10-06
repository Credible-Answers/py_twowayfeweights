"""Write ``docs/examples.md`` from the Stata comparison cases.

    python tests/reference/build_docs_examples.py

Each example shows the Stata command, the equivalent copy-paste Python code, the Python output and
(collapsed) the Stata output stored in ``tests/fixtures/stata/results.json``.
"""

from __future__ import annotations

import contextlib
import io
import json
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))

from cases import CASES, DATA_DIR  # noqa: E402

DOCS = ROOT / "docs" / "examples.md"
REF = json.loads((ROOT / "tests" / "fixtures" / "stata" / "results.json").read_text(encoding="utf-8"))["cases"]
GENTZKOW_URL = "https://raw.githubusercontent.com/Credible-Answers/py_twowayfeweights/main/tests/data/gentzkow.csv.gz"

# (case name, title, plain-language description)
EXAMPLES = [
    ("wp_feTR", "Fixed effects regression (feTR)",
     "The basic case: a regression of log wages on union membership with worker and year fixed effects. "
     "`test_random_weights` checks whether the weights are correlated with education."),
    ("wp_feTR_ctrl_w", "Fixed effects with controls and weights",
     "Same regression, adding two control variables and a weight variable `wt`, and testing two variables."),
    ("wp_feTR_ot", "Several treatments (other_treatments)",
     "The regression also includes marriage as a second treatment. The output then shows the weights "
     "attached to each treatment."),
    ("wp_feTR_ot_ctrl_w", "Several treatments with controls and weights",
     "Two other treatments, one control and weights."),
    ("wp_feS", "Fixed effects, stable treatment effects (feS)",
     "Same regression as the first example, under the extra assumption that each group's treatment "
     "effect does not change over time."),
    ("wp_feS_ctrl_w", "feS with controls and weights", "The feS version of the second example."),
    ("wp_fdTR", "First-difference regression (fdTR)",
     "A regression of the change in log wages on the change in union membership. `D0` is the treatment "
     "level (not differenced)."),
    ("wp_fdTR_ctrl_w", "fdTR with a control and weights", "The first-difference regression with a control and weights."),
    ("wp_fdS", "First differences, stable treatment effects (fdS)",
     "First-difference regression under the extra assumption that treatment effects are stable over time."),
    ("wp_fdS_ctrl_w", "fdS with a control and weights", "The fdS version with a control and weights."),
    ("gz_fdTR_styr_rw", "Newspapers and turnout, 683 controls",
     "Chapter 5 of the DiD book, with data from Gentzkow, Shapiro and Sinkinson (2011): the effect of the "
     "number of daily newspapers in a county on presidential turnout, with 683 state-year dummies as controls."),
]


def _arg(v) -> str:
    """``"educ"`` for one variable, ``["educ", "exper"]`` for several."""
    v = [v] if isinstance(v, str) else list(v)
    return json.dumps(v[0]) if len(v) == 1 else json.dumps(v)


def python_code(case: dict) -> str:
    c = case
    lines = ["# pip install twowayfeweights"]
    if c["data"] == "wagepan":
        if c["weights"]:
            lines.append("import numpy as np")
        lines += ["from twowayfeweights import twowayfeweights, load_wagepan", "", "df = load_wagepan()"]
        if c["weights"]:
            lines += ["# a weight variable (the same one is used in Stata)",
                      'df["wt"] = np.round(0.5 + np.random.default_rng(2020).random(len(df)), 3)']
        controls = c["controls"]
    else:
        lines += [
            "import pandas as pd",
            "from twowayfeweights import twowayfeweights",
            "",
            f'df = pd.read_csv("{GENTZKOW_URL}")',
            "# one dummy per state-year (683 dummies), used as controls",
            'styr = pd.get_dummies(df["styr"], prefix="styr", dtype=float)',
            "df = pd.concat([df, styr], axis=1)",
        ]
        controls = None
    args = [f'Y="{c["Y"]}"', f'G="{c["G"]}"', f'T="{c["T"]}"', f'D="{c["D"]}"', f'type="{c["type"]}"']
    if c["D0"]:
        args.append(f'D0="{c["D0"]}"')
    if c["data"] != "wagepan" and c["controls"]:
        args.append("controls=list(styr.columns)")
    elif controls:
        args.append(f"controls={_arg(controls)}")
    if c["weights"]:
        args.append(f'weights="{c["weights"]}"')
    if c["other_treatments"]:
        args.append(f"other_treatments={_arg(c['other_treatments'])}")
    if c["test_random_weights"]:
        args.append(f"test_random_weights={_arg(c['test_random_weights'])}")
    args.append("summary_measures=True")
    rows, cur = [], "res = twowayfeweights(df, "
    for a in args:
        piece = a + ", "
        if len(cur) + len(piece) > 100:
            rows.append(cur.rstrip())
            cur = " " * len("res = twowayfeweights(") + piece
        else:
            cur += piece
    rows.append(cur.rstrip().rstrip(",") + ")")
    lines += ["", *rows, "print(res)"]
    return "\n".join(lines)


def run(code: str, local_data: bool = True) -> str:
    """Run an example and return what it prints (the Gentzkow data is read from the repository copy)."""
    if local_data:
        code = code.replace(GENTZKOW_URL, str(DATA_DIR / "gentzkow.csv.gz").replace("\\", "/"))
    buf = io.StringIO()
    with warnings.catch_warnings(), contextlib.redirect_stdout(buf):
        warnings.simplefilter("ignore")
        exec(compile(code, "<example>", "exec"), {})
    return buf.getvalue().strip("\n")


def build() -> str:
    cases = {c["name"]: c for c in CASES}
    out = [
        "# Examples",
        "",
        "Every example on this page is one of the specifications the test suite compares with Stata.",
        "Each code block is complete: copy it, paste it into Python, and run it. The output shown is",
        "what Python prints; click **Stata output** to see what the Stata command prints.",
        "",
        "!!! note",
        "    In the random-weights table, the last of the eight digits can differ from Stata's.",
        "    Stata stores some intermediate variables as 4-byte floats; Python uses double precision.",
        "",
    ]
    for name, title, text in EXAMPLES:
        c = cases[name]
        ref = REF[name]
        out += [
            f"## {title}",
            "",
            text,
            "",
            "**Stata**",
            "",
            "```stata",
            ref["command"],
            "```",
            "",
            "**Python**",
            "",
            "```python",
            python_code(c),
            "```",
            "",
            "**Output**",
            "",
            "```text",
            run(python_code(c)),
            "```",
            "",
            '??? quote "Stata output"',
            "",
            "    ```text",
            *[("    " + line).rstrip() for line in ref["output"].split("\n")],
            "    ```",
            "",
        ]
    return "\n".join(out).rstrip() + "\n"


if __name__ == "__main__":
    DOCS.write_text(build(), encoding="utf-8")
    print(f"wrote {DOCS.relative_to(ROOT)}")
