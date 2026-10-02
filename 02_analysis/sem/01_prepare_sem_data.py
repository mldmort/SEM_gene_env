#!/usr/bin/env python3
"""
Analysis-ready input tables for the parental-age SEM (one per cohort).

Rows: DNM children (gxe.analysis.dnm_children: in the callset, zero-count and
outlier filters), restricted to config `analysis.sem.ancestry` and to
case_status in {case, control}.

Columns written (variable names used by 02_sem_lavaan.R):
  FID              family, for clustered standard errors
  asd              1 = case, 0 = control
  sex_male         1 / 0 (manifest SEX is PLINK-coded: 1 = male, 2 = female)
  father_age, mother_age   age at the child's birth, in DECADES (years / 10)
  dnm              damaging de novo count (sum of config `sem.dnm`)
  rv_<name>        damaging inherited rare-variant count for each definition in config
                   `sem.rv_definitions` (sum of its columns)
  edu, pgs         child's EA and ASD PGS (already z-scored within cohort)
  seq covariates   coverage (child), mother_coverage, father_coverage, mother_mapped_pct,
                   father_mapped_pct: mean sequencing depth and CRAM mapped-read %, z-scored
                   within this table (config `sem.seq_covariates`)
  PC1..PCk         proj_PC1..k, z-scored within this table
  birth_order      among the family's children, capped at 4 (descriptive)
  bo_*             birth-order terms for the parental-age and ASD equations: bo_linear, or
                   bo_2 / bo_3 / bo_4plus indicators vs first-born (config sem.birth_order_coding)
  src_<source>     0/1 genotype-source indicators, reference = most common source
                   (only when the cohort has more than one source, i.e. SPARK WES)
Missing values stay missing; each model drops incomplete rows for its own variables.

Output: results/sem/input_<cohort>.tsv, results/sem/input_summary.tsv
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import pandas as pd

from gxe import load_config, log
from gxe.analysis import dnm_children, results_path


def zscore(s):
    s = pd.to_numeric(s, errors="coerce")
    return (s - s.mean()) / s.std()


def main():
    cfg = load_config()
    sem = cfg["analysis"]["sem"]
    pcs = [f"proj_PC{i}" for i in range(1, sem["n_pcs"] + 1)]
    rv_defs = sem["rv_definitions"]
    rv_cols = sorted({c for cols in rv_defs.values() for c in cols})
    seq = sem["seq_covariates"]
    need = ["FID", "source", "case_status", "SEX", "predicted_ancestry", "birth_order", *seq.values(),
            *sem["parent_ages"], *sem["dnm"], *rv_cols, sem["pgs_asd"], sem["pgs_edu"], *pcs]

    summary = []
    for cohort in sem["cohorts"]:
        df, flow = dnm_children(cohort, columns=need, cfg=cfg)
        df = df[df["predicted_ancestry"] == sem["ancestry"]]
        df = df[df["case_status"].isin(["case", "control"])]
        log(f"  {sem['ancestry']} cases + controls: {len(df):,}")

        out = pd.DataFrame({
            "FID": df["FID"].astype(str),
            "asd": (df["case_status"] == "case").astype(int),
            "sex_male": pd.to_numeric(df["SEX"], errors="coerce").map({1: 1, 2: 0}),
            "father_age": df["father_age_at_child_birth"] / 10,
            "mother_age": df["mother_age_at_child_birth"] / 10,
            "dnm": df[sem["dnm"]].sum(axis=1),
            **{f"rv_{name}": df[cols].sum(axis=1) for name, cols in rv_defs.items()},
            "edu": df[sem["pgs_edu"]],
            "pgs": df[sem["pgs_asd"]],
        })
        for name, col in seq.items():
            out[name] = zscore(df[col])
        # Birth order among the family's children, capped at 4 (4 = fourth or later), coded per
        # config `sem.birth_order_coding`; the R script uses every bo_* column it finds.
        bo = pd.to_numeric(df["birth_order"], errors="coerce").clip(upper=4)
        out["birth_order"] = bo
        if sem["birth_order_coding"] == "linear":
            out["bo_linear"] = bo
        else:
            for k, name in ((2, "bo_2"), (3, "bo_3"), (4, "bo_4plus")):
                out[name] = (bo == k).astype(float).where(bo.notna())
        # Exome / genotyping source: indicator per source except the most common (WES only).
        sources = df["source"].value_counts()
        for src in sources.index[1:]:
            out[f"src_{src}"] = (df["source"] == src).astype(int)
        if len(sources) > 1:
            log(f"  source indicators vs {sources.index[0]}: {', '.join(sources.index[1:])}")
        for i, c in enumerate(pcs, 1):
            out[f"PC{i}"] = zscore(df[c])

        path = results_path("sem", f"input_{cohort}.tsv")
        out.to_csv(path, sep="\t", index=False, na_rep="NA")
        log(f"  wrote {path}")
        summary.append({
            "cohort": cohort, **{f"flow_{k}": v for k, v in flow.items() if k in ("children", "in_dnm_callset")},
            "n": len(out), "cases": int(out["asd"].sum()), "controls": int((out["asd"] == 0).sum()),
            "families": out["FID"].nunique(),
            "father_age_n": int(out["father_age"].notna().sum()),
            "mother_age_n": int(out["mother_age"].notna().sum()),
            "birth_order_n": int(out["birth_order"].notna().sum()),
            "dnm_mean": round(out["dnm"].mean(), 4), "dnm_carriers": int((out["dnm"] > 0).sum()),
            **{f"rv_{name}_mean": round(out[f"rv_{name}"].mean(), 4) for name in rv_defs},
            **{f"rv_{name}_carriers": int((out[f"rv_{name}"] > 0).sum()) for name in rv_defs},
        })
    s = pd.DataFrame(summary)
    s.to_csv(results_path("sem", "input_summary.tsv"), sep="\t", index=False)
    log("\n" + s.to_string(index=False))


if __name__ == "__main__":
    main()
