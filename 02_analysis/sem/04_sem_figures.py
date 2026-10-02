#!/usr/bin/env python3
"""
Figures for the parental-age SEM.

  sem_path_<parent>.png   path diagram (layout of the conceptual model) with the
                          fixed-effect meta estimate on every arrow; blue = positive,
                          orange = negative, grey dashed = p >= 0.05; width ~ |z|
  sem_forest.png          per-cohort and meta estimates (95% CI) for the key paths,
                          one row per parent

Inputs : results/sem/sem_estimates.tsv, sem_meta.tsv, sem_fit.tsv
Outputs: results/sem/figures/
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from gxe import load_config, log
from gxe.analysis import results_path

INK, INK_2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
POS, NEG, NS = "#2a78d6", "#eb6834", "#a3a29c"
COHORT_COLORS = {"spark_wgs": "#2a78d6", "ssc": "#eb6834", "spark_wes": "#1baf7a"}
COHORT_NAMES = {"spark_wgs": "SPARK WGS", "ssc": "SSC", "spark_wes": "SPARK WES", "meta": "Meta"}
PARENT_NAMES = {"father_age": "Paternal age", "mother_age": "Maternal age"}

# Node centres (x, y) and labels, following the conceptual figure.
NODES = {
    "age": ((0.15, 0.60), "{parent}\n(decades)"),
    "dnm": ((0.50, 0.88), "DNM\ndamaging de novo"),
    "asd": ((0.85, 0.60), "ASD\n(probit)"),
    "rv": ((0.50, 0.46), "RV\ninherited damaging"),
    "edu": ((0.50, 0.27), "EDU\nEA PGS"),
    "pgs": ((0.50, 0.07), "PGS\nASD PGS (loSSC)"),
}
NODE_W, NODE_H = 0.17, 0.10
# label -> (polyline waypoints from source border to target border, label position)
EDGES = {
    "g": ([(0.20, 0.65), (0.415, 0.88)], (0.27, 0.80)),
    "d": ([(0.585, 0.88), (0.80, 0.65)], (0.73, 0.80)),
    "c": ([(0.235, 0.60), (0.765, 0.60)], (0.50, 0.635)),
    "rp": ([(0.415, 0.46), (0.21, 0.55)], (0.33, 0.535)),
    "ra": ([(0.585, 0.46), (0.79, 0.55)], (0.67, 0.535)),
    "re": ([(0.50, 0.41), (0.50, 0.32)], (0.565, 0.365)),
    "ep": ([(0.415, 0.27), (0.33, 0.27), (0.15, 0.55)], (0.235, 0.36)),
    "ea": ([(0.585, 0.27), (0.67, 0.27), (0.85, 0.55)], (0.765, 0.36)),
    "pe": ([(0.50, 0.12), (0.50, 0.22)], (0.565, 0.17)),
    "pp": ([(0.415, 0.07), (0.25, 0.07), (0.075, 0.55)], (0.12, 0.22)),
    "pa": ([(0.585, 0.07), (0.75, 0.07), (0.925, 0.55)], (0.88, 0.22)),
}
GREEK = {"g": "γ", "d": "δ", "c": "c", "rp": "r_p", "ra": "r_a", "re": "r_e", "ep": "e_p",
         "ea": "e_a", "pe": "p_e", "pp": "p_p", "pa": "p_a"}


def stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else ""


def path_diagram(meta, parent, fit, out):
    m = meta[(meta["parent"] == parent) & (meta["sample"] == "all")].set_index("label")
    fig, ax = plt.subplots(figsize=(11, 7.2))
    ax.set_xlim(0, 1); ax.set_ylim(-0.07, 1.04); ax.axis("off")

    for name, ((x, y), text) in NODES.items():
        ax.add_patch(FancyBboxPatch((x - NODE_W / 2, y - NODE_H / 2), NODE_W, NODE_H,
                                    boxstyle="round,pad=0.008,rounding_size=0.02",
                                    fc="white", ec=INK, lw=1.4, zorder=3))
        ax.text(x, y, text.format(parent=PARENT_NAMES[parent]), ha="center", va="center",
                fontsize=11, color=INK, zorder=4, linespacing=1.3)

    zmax = m.loc[list(EDGES), "z"].abs().max()
    for lab, (pts, (lx, ly)) in EDGES.items():
        r = m.loc[lab]
        sig = r["pvalue"] < 0.05
        color = (POS if r["est"] > 0 else NEG) if sig else NS
        lw = 1.2 + 4.0 * min(abs(r["z"]) / zmax, 1) if sig else 1.2
        ls = "-" if sig else (0, (4, 3))
        for a, b in zip(pts[:-2], pts[1:-1]):
            ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=lw, ls=ls, zorder=2,
                    solid_capstyle="round")
        ax.add_patch(FancyArrowPatch(pts[-2], pts[-1], arrowstyle="-|>,head_length=9,head_width=5",
                                     color=color, lw=lw, ls=ls, zorder=2, shrinkA=0, shrinkB=1))
        ax.text(lx, ly, f"{GREEK[lab]} {r['est']:+.3f}{stars(r['pvalue'])}", ha="center", va="center",
                fontsize=10.5, color=INK, zorder=5,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=GRID, lw=0.6))

    f = fit[(fit["parent"] == parent) & (fit["sample"] == "all")]
    ind, conf, c = (m.loc[k] for k in ("ind_dnm", "conf_total", "c"))
    ax.text(0.0, 1.035, f"{PARENT_NAMES[parent]} → ASD: fixed-effect meta of SPARK WGS, SSC, SPARK WES "
            f"(n = {int(f['n'].sum()):,} EUR children)", fontsize=12, color=INK, va="top")
    ax.text(0.0, -0.065,
            f"Indirect causal γ·δ = {ind['est']:+.4f}{stars(ind['pvalue'])}    "
            f"Genetic confounding, total = {conf['est']:+.4f}{stars(conf['pvalue'])}    "
            f"Direct c = {c['est']:+.4f}{stars(c['pvalue'])}        "
            "\nArrows: blue = positive, orange = negative (p < 0.05), width ∝ |z|; grey dashed = p ≥ 0.05.  * p<0.05  ** p<0.01  *** p<0.001",
            fontsize=8.5, color=INK_2, va="bottom")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log(f"wrote {out}")


def forest(est, meta, out):
    paths = [("g", "γ  age → DNM"), ("d", "δ  DNM → ASD"), ("ind_dnm", "γ·δ  indirect"),
             ("c", "c  direct"), ("conf_total", "genetic confounding"), ("ep", "e_p  EDU → age")]
    rows = list(COHORT_NAMES)
    fig, axes = plt.subplots(2, len(paths), figsize=(17, 6.2), sharey=True, constrained_layout=True)
    for i, parent in enumerate(PARENT_NAMES):
        for j, (lab, title) in enumerate(paths):
            ax = axes[i, j]
            e = est[(est["parent"] == parent) & (est["sample"] == "all") & (est["label"] == lab)].set_index("cohort")
            mm = meta[(meta["parent"] == parent) & (meta["sample"] == "all") & (meta["label"] == lab)]
            ax.axvline(0, color=INK_2, lw=0.8, ls=":")
            for k, cohort in enumerate(rows):
                r = mm.iloc[0] if cohort == "meta" else e.loc[cohort]
                y = len(rows) - 1 - k
                color = INK if cohort == "meta" else COHORT_COLORS[cohort]
                ax.plot([r["ci.lower"], r["ci.upper"]], [y, y], color=color, lw=2)
                ax.plot(r["est"], y, marker="D" if cohort == "meta" else "o", ms=8 if cohort == "meta" else 7,
                        color=color, mec="white", mew=1.2, zorder=3)
            het = mm.iloc[0]["het_p"]
            ax.set_title(f"{title}\nhet p = {het:.2g}", fontsize=10, color=INK)
            ax.set_yticks(range(len(rows)), [COHORT_NAMES[c] for c in reversed(rows)], fontsize=9)
            ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            ax.tick_params(colors=INK_2, labelsize=8)
        axes[i, 0].set_ylabel(PARENT_NAMES[parent], fontsize=11, color=INK)
    fig.suptitle("SEM paths per cohort and fixed-effect meta (95% CI). Age in decades; ASD paths probit",
                 x=0.01, ha="left", fontsize=12, color=INK)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    log(f"wrote {out}")


def main():
    rv = load_config()["analysis"]["sem"]["rv_primary"]
    est = pd.read_csv(results_path("sem", "sem_estimates.tsv"), sep="\t").query("rv_def == @rv")
    meta = pd.read_csv(results_path("sem", "sem_meta.tsv"), sep="\t").query("rv_def == @rv")
    fit = pd.read_csv(results_path("sem", "sem_fit.tsv"), sep="\t").query("rv_def == @rv")
    log(f"figures for RV definition: {rv}")
    for parent in PARENT_NAMES:
        path_diagram(meta, parent, fit, results_path("sem", "figures", f"sem_path_{parent}.png"))
    forest(est, meta, results_path("sem", "figures", "sem_forest.png"))


if __name__ == "__main__":
    main()
