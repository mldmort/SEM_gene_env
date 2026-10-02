#!/usr/bin/env python3
"""
SSC phenotypes, one row per sequenced SSC individual (key: sample_id = SS_RUID, e.g. SSC02219).

Sources (SSC Version 15 phenotype set, config `inputs.ssc.pheno_dir`):
  SSC IDs Sanders.csv            SS_RUID -> individual ("11000.p1"), family, relation, sex
  SSC_Age_of_Subjects_&_...csv   birth month (YYYY-MM), parental ages at birth, and
                                 birth order among the family's children
  */ssc_core_descriptive.csv     proband + sibling core measures (IQ, ADOS, SRS, Vineland, ...)
  Proband Data/ssc_background_hx highest_edu_mother / highest_edu_father, one row per family;
                                 applied to every child in the family
  Air Pollution/...csv           prenatal and enrolment PM2.5 / O3 / NO2
  SFARI smoking file             parental smoking, keyed by the child

SSC has no ADI and no household income, so those columns are absent here.
IQ columns are renamed to the SPARK names: fsiq, viq, nviq.

Output: data/intermediate/phenotypes/ssc_phenotypes.parquet
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from gxe import intermediate_path, left_merge, load_config, log, write_table
from gxe.phenotypes import AQ_COLS, birth_order, edu_code, load_sfari_smoking, numeric, parent_age

IQ_RENAME = {
    "ssc_diagnosis_full_scale_iq": "fsiq",
    "ssc_diagnosis_verbal_iq": "viq",
    "ssc_diagnosis_nonverbal_iq": "nviq",
    "ssc_diagnosis_full_scale_iq_type": "fsiq_test_type",
    "ssc_diagnosis_verbal_iq_type": "viq_test_type",
    "ssc_diagnosis_nonverbal_iq_type": "nviq_test_type",
}


def year_month_to_dec(s):
    """'1999-03' -> 1999.208 (mid-month)."""
    parts = s.astype("string").str.extract(r"^(\d{4})-(\d{1,2})$")
    return (pd.to_numeric(parts[0], errors="coerce")
            + (pd.to_numeric(parts[1], errors="coerce") - 0.5) / 12).round(3)


def load_sanders(d):
    df = pd.read_csv(d / "SSC IDs Sanders.csv", dtype=str, encoding="utf-8-sig",
                     usecols=["SS_RUID", "individual", "family", "relation", "sex"])
    df = df.rename(columns={"SS_RUID": "sample_id", "individual": "ssc_individual",
                            "family": "ssc_family", "relation": "ssc_relation",
                            "sex": "ssc_sex"})
    df["ssc_role"] = df["ssc_relation"].str[0].map({"p": "proband", "s": "sibling"}).fillna(
        df["ssc_relation"].map({"fa": "father", "mo": "mother"}))
    log(f"  Sanders: {len(df):,} sequenced individuals")
    return df.drop_duplicates("sample_id")


def load_ages(d, min_age):
    df = pd.read_csv(d / "SSC_Age_of_Subjects_&_ParentalAgeAtBirth.csv", dtype=str,
                     usecols=["id()", "birth", "father_birth", "mother_birth",
                              "father_age_birth", "mother_age_birth"])
    df = df.rename(columns={"id()": "ssc_individual"}).drop_duplicates("ssc_individual")
    birth = year_month_to_dec(df["birth"])
    df["birth_year_dec"] = birth
    for side in ("mother", "father"):
        age = parent_age(birth, year_month_to_dec(df[f"{side}_birth"]), min_age)
        df[f"{side}_age_at_child_birth"] = age
        # Cross-check against the age the SSC file ships (units not documented).
        shipped = pd.to_numeric(df[f"{side}_age_birth"], errors="coerce")
        both = age.notna() & (shipped > 0)
        log(f"  {side} age from birth months: n={age.notna().sum():,}; "
            f"median shipped/derived ratio {(shipped[both] / age[both]).median():.2f}")
    # Birth order among the family's children (probands p*, siblings s*) with a birth month.
    fam, rel = df["ssc_individual"].str.split(".", n=1).str[0], df["ssc_individual"].str.split(".", n=1).str[1]
    is_child = rel.fillna("").str.match(r"^[ps]\d")
    df["birth_order"], df["n_children_known_birth"] = birth_order(fam.where(is_child), birth.where(is_child))
    log(f"  birth_order: {df['birth_order'].notna().sum():,}")
    return df[["ssc_individual", "birth_year_dec", "birth_order", "n_children_known_birth",
               "mother_age_at_child_birth", "father_age_at_child_birth"]]


def load_core_descriptive(d):
    parts = [pd.read_csv(p, dtype=str, low_memory=False) for p in sorted(d.glob("*/ssc_core_descriptive.csv"))]
    df = pd.concat(parts, ignore_index=True).drop(columns=["sex", "measure"], errors="ignore")
    df = df.rename(columns={"individual": "ssc_individual", **IQ_RENAME})
    df = df.drop_duplicates("ssc_individual")
    log(f"  core descriptive: {len(df):,} individuals from {len(parts)} files")
    return df


def load_parent_education(d):
    df = pd.read_csv(d / "Proband Data" / "ssc_background_hx.csv", dtype=str,
                     usecols=["individual", "highest_edu_mother", "highest_edu_father"])
    df["ssc_family"] = df["individual"].str.split(".").str[0]
    df = df.rename(columns={"highest_edu_mother": "mother_highest_education",
                            "highest_edu_father": "father_highest_education"})
    df = df.drop(columns="individual").drop_duplicates("ssc_family")
    df["mother_edu_code"] = edu_code(df["mother_highest_education"])
    df["father_edu_code"] = edu_code(df["father_highest_education"])
    return df


def load_air_pollution(d):
    df = pd.read_csv(d / "Air Pollution" / "ssc_subjects_with_aq_exposure.csv", dtype=str,
                     usecols=["subject_id"] + AQ_COLS)
    df["ssc_individual"] = df.pop("subject_id").str.replace("-", ".", regex=False)
    df = numeric(df, AQ_COLS).drop_duplicates("ssc_individual")
    log(f"  air pollution: {len(df):,} individuals")
    return df


def main():
    cfg = load_config()
    d = Path(cfg["inputs"]["ssc"]["pheno_dir"])
    log("SSC phenotypes")
    df = load_sanders(d)
    for part in (load_ages(d, cfg["phenotypes"]["min_parent_age"]),
                 load_core_descriptive(d), load_air_pollution(d)):
        df = left_merge(df, part, on="ssc_individual")

    # Family-level parental education goes on children only.
    edu = load_parent_education(d)
    is_child = df["ssc_role"].isin(["proband", "sibling"])
    df = left_merge(df, edu, on="ssc_family")
    df.loc[~is_child, edu.columns.drop("ssc_family")] = pd.NA

    df = left_merge(df, load_sfari_smoking(cfg["inputs"]["smoking"], "SSC"))
    for c in ("birth_year_dec", "mother_age_at_child_birth", "father_age_at_child_birth",
              "mother_edu_code", "fsiq", "prenatal_pm25_mean", "father_ever_smoke"):
        log(f"    {c}: {df[c].notna().sum():,}")
    write_table(df, intermediate_path("phenotypes", "ssc_phenotypes.parquet"))


if __name__ == "__main__":
    main()
