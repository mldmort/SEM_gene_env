#!/usr/bin/env python3
"""
Per-cohort and meta-analysed SEM results side by side.

Reads results/sem/sem_estimates.tsv and sem_meta.tsv (02_sem_lavaan.R) and writes:
  results/sem/sem_report.tsv   long: parent, sample, label, cohort (or "meta"), est, se, CI, p
  results/sem/sem_report.md    one table per parent x sample: rows = paths, columns = cohorts + meta
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import pandas as pd

from gxe import log
from gxe.analysis import results_path

COHORTS = ["spark_wgs", "ssc", "spark_wes"]
LABELS = {
    "g": "γ  parent age → DNM",
    "d": "δ  DNM → ASD",
    "c": "c  parent age → ASD (direct)",
    "asd~bo_linear": "birth order → ASD (per rank)",
    "asd~bo_2": "2nd-born vs 1st → ASD",
    "asd~bo_3": "3rd-born vs 1st → ASD",
    "asd~bo_4plus": "4th+ vs 1st → ASD",
    "rp": "r_p  RV → parent age",
    "ra": "r_a  RV → ASD",
    "re": "r_e  RV → EDU",
    "ep": "e_p  EDU → parent age",
    "ea": "e_a  EDU → ASD",
    "pp": "p_p  PGS → parent age",
    "pa": "p_a  PGS → ASD",
    "pe": "p_e  PGS → EDU",
    "ind_dnm": "indirect causal γ·δ",
    "conf_rv": "confounding via RV",
    "conf_edu": "confounding via EDU",
    "conf_pgs": "confounding via PGS",
    "conf_total": "confounding, total",
    "total_causal": "c + γ·δ",
}


def cell(r):
    star = "**" if r["pvalue"] < 0.05 else ""
    return f"{star}{r['est']:+.4f}{star} [{r['ci.lower']:+.4f}, {r['ci.upper']:+.4f}] p={r['pvalue']:.2g}"


def main():
    est = pd.read_csv(results_path("sem", "sem_estimates.tsv"), sep="\t")
    meta = pd.read_csv(results_path("sem", "sem_meta.tsv"), sep="\t")
    fit = pd.read_csv(results_path("sem", "sem_fit.tsv"), sep="\t")
    keep = ["rv_def", "parent", "sample", "label", "est", "se", "ci.lower", "ci.upper", "pvalue"]
    long = pd.concat([est[["cohort"] + keep],
                      meta.assign(cohort="meta")[["cohort"] + keep + ["het_p"]]], ignore_index=True)
    long = long[long["label"].isin(LABELS)]
    long.to_csv(results_path("sem", "sem_report.tsv"), sep="\t", index=False)

    lines = ["# Parental-age SEM: per cohort and fixed-effect meta",
             "", "Parental age in decades; PGS z-scored within cohort; ASD paths on the probit scale. "
             "Cells: estimate [95% CI] p; bold = p < 0.05. het p = between-cohort heterogeneity.", ""]
    # Side-by-side comparison of the RV definitions (meta estimates, all children).
    rv_defs = list(dict.fromkeys(long["rv_def"]))
    lines += ["## RV definitions compared (meta, all children)", ""]
    for parent in dict.fromkeys(long["parent"]):
        lines.append(f"**{parent.replace('_', ' ')}**")
        lines.append("")
        lines.append("| | " + " | ".join(rv_defs) + " |")
        lines.append("|---" * (len(rv_defs) + 1) + "|")
        for lab, name in LABELS.items():
            cells = []
            for rd in rv_defs:
                r = long[(long["rv_def"] == rd) & (long["parent"] == parent) & (long["sample"] == "all")
                         & (long["label"] == lab) & (long["cohort"] == "meta")]
                cells.append(cell(r.iloc[0]) if len(r) else "")
            lines.append(f"| {name} | " + " | ".join(cells) + " |")
        lines.append("")
    for (rv_def, parent, sample), g in long.groupby(["rv_def", "parent", "sample"], sort=False):
        f = fit[(fit["rv_def"] == rv_def) & (fit["parent"] == parent) & (fit["sample"] == sample)].set_index("cohort")
        lines += [f"## RV = {rv_def} | {parent.replace('_', ' ')} — {'all children' if sample == 'all' else 'controls only'}", ""]
        lines.append("| | " + " | ".join(f"{c} (n={int(f.loc[c, 'n']):,})" if c in f.index else c for c in COHORTS)
                     + " | meta | het p |")
        lines.append("|---" * (len(COHORTS) + 3) + "|")
        for lab, name in LABELS.items():
            rows = g[g["label"] == lab].set_index("cohort")
            if rows.empty:
                continue
            cells = [cell(rows.loc[c]) if c in rows.index else "" for c in COHORTS + ["meta"]]
            het = rows.loc["meta", "het_p"] if "meta" in rows.index else float("nan")
            lines.append(f"| {name} | " + " | ".join(cells) + f" | {het:.2g} |")
        lines.append("")
    md = results_path("sem", "sem_report.md")
    md.write_text("\n".join(lines))
    log(f"wrote {md} and sem_report.tsv")


if __name__ == "__main__":
    main()
