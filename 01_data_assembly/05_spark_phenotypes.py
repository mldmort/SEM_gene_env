#!/usr/bin/env python3
"""
SPARK phenotypes, one row per registered individual (children and adults).

Sources (release set in config `inputs.spark`):
  individuals_registration   IDs, sex, ASD status, race/ethnicity, birth year
                             birth_year_dec = registration_year - age_at_registration_months / 12
                             (unusable when age_validity_flag == 1 or either field is missing)
  core_descriptive_variables IQ, cognitive impairment, SCQ, RBS-R, Vineland, DCDQ, milestones
  background_history_*       parental education and household income
                             (child form, then sibling form, then the adult's own form)
  area_deprivation_index     ADI 2015 / 2019 / 2020 / 2022 (latest evaluation per subject)
  air pollution              prenatal and enrolment PM2.5 / O3 / NO2
  asr                        Adult Self-Report substance use, on the adult's own row and
                             copied to their children as mother_asr_* / father_asr_*
  SFARI smoking file         parental smoking, keyed by the child
  withdrawn / paused lists   flags spark_withdrawn / spark_paused

Derived per child: birth_order (among all registered children of the same biological
mother, or father if no mother is listed), mother_/father_age_at_child_birth, mother_/father_edu_code (1..10),
annual_household_income_ordinal (1..9).

Output: data/intermediate/phenotypes/spark_phenotypes.parquet  (key: sample_id = subject_sp_id)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from gxe import intermediate_path, left_merge, load_config, log, write_table
from gxe.phenotypes import (AQ_COLS, birth_order, edu_code, income_ordinal, load_sfari_smoking,
                            numeric, parent_age)

REG_COLS = [
    "subject_sp_id", "family_sf_id", "biomother_sp_id", "biofather_sp_id", "sex", "asd",
    "age_at_registration_months", "registration_year", "age_validity_flag",
    "individual_data_validity", "asd_validity_flag", "current_depend_adult",
    "diagnosis_age", "cognitive_impairment_initial", "language_level_initial",
    "num_asd_parents", "num_asd_siblings", "multiple_birth",
    "race_asian", "race_african_amer", "race_native_amer", "race_native_hawaiian",
    "race_white", "race_other", "race_more_than_one_calc", "hispanic",
]
# Core descriptive columns already taken from the registration file.
CDV_SKIP = {"family_sf_id", "biomother_id", "biofather_id", "sex", "asd", "current_depend_adult",
            "age_at_registration_months", "age_at_registration_years", "registration_year",
            "diagnosis_age", "enrolled_multiple", "asd_validity_flag"}
ADI_SPECIAL = {"GQ", "PH", "GQ-PH", "QDI"}
ASR_ITEMS = ["q006_uses_drugs", "q090_drinks_too_much", "q124_used_tobacco", "q126_used_drugs"]
ASR_SCALES = ["tobacco_raw_score", "tobacco_t_score", "alcohol_raw_score", "alcohol_t_score",
              "drugs_raw_score", "drugs_t_score", "substance_mean", "substance_t_score"]
ASR_TO_PARENTS = ["tobacco_user", "alcohol_user", "q124_used_tobacco",
                  "tobacco_raw_score", "alcohol_raw_score", "drugs_raw_score", "substance_t_score"]


def release_file(cfg, name):
    s = cfg["inputs"]["spark"]
    return Path(s["release_dir"]) / f"{name}-{s['release_tag']}.csv"


def supplement_file(cfg, name):
    s = cfg["inputs"]["spark"]
    return Path(s["release_dir"]) / "Supplemental Information and Data Dictionary" / f"{name}-{s['release_tag']}.csv"


def load_registration(cfg):
    df = pd.read_csv(release_file(cfg, "individuals_registration"), dtype=str,
                     usecols=REG_COLS, low_memory=False)
    df = df.rename(columns={"subject_sp_id": "sample_id"}).drop_duplicates("sample_id")
    age_mo = pd.to_numeric(df["age_at_registration_months"], errors="coerce")
    reg_yr = pd.to_numeric(df["registration_year"], errors="coerce")
    usable = (df["age_validity_flag"].fillna("").str.strip() != "1") & age_mo.notna() & reg_yr.notna()
    df["birth_year_dec"] = (reg_yr - age_mo / 12.0).where(usable).round(3)
    log(f"  registration: {len(df):,} individuals, {usable.sum():,} with a usable birth year")

    birth = df.set_index("sample_id")["birth_year_dec"]
    min_age = cfg["phenotypes"]["min_parent_age"]
    for side, col in (("mother", "biomother_sp_id"), ("father", "biofather_sp_id")):
        df[f"{side}_age_at_child_birth"] = parent_age(df["birth_year_dec"], df[col].map(birth), min_age)
        log(f"    {side}_age_at_child_birth: {df[f'{side}_age_at_child_birth'].notna().sum():,}")

    # Birth order among all registered children of the same biological mother
    # (biological father when no mother is listed), sequenced or not.
    mom = df["biomother_sp_id"].where(df["biomother_sp_id"].fillna("").str.strip() != "")
    dad = df["biofather_sp_id"].where(df["biofather_sp_id"].fillna("").str.strip() != "")
    group = ("M:" + mom).fillna("F:" + dad)
    df["birth_order"], df["n_children_known_birth"] = birth_order(group, df["birth_year_dec"])
    log(f"    birth_order: {df['birth_order'].notna().sum():,}")
    return df


def load_core_descriptive(cfg):
    df = pd.read_csv(release_file(cfg, "core_descriptive_variables"), dtype=str, low_memory=False)
    df = df.drop(columns=[c for c in CDV_SKIP if c in df.columns])
    df = df.rename(columns={"subject_sp_id": "sample_id"}).drop_duplicates("sample_id")
    log(f"  core descriptive: {len(df):,} individuals, {df.shape[1] - 1} columns")
    return df


def load_background(cfg):
    """Parental education and household income: child form > sibling form > adult form."""
    cols = ["subject_sp_id", "mother_highest_education", "father_highest_education",
            "annual_household_income"]
    parts = []
    for form in ("child", "sibling", "adult"):
        extra = ["highest_level_education"] if form == "adult" else []
        f = pd.read_csv(release_file(cfg, f"background_history_{form}"), dtype=str,
                        usecols=cols + extra, low_memory=False)
        f["bghx_form"] = form
        parts.append(f)
        log(f"  background_history_{form}: {len(f):,} rows")
    df = pd.concat(parts, ignore_index=True).rename(columns={"subject_sp_id": "sample_id"})
    df["annual_household_income_ordinal"] = income_ordinal(df["annual_household_income"])
    df["mother_edu_code"] = edu_code(df["mother_highest_education"])
    df["father_edu_code"] = edu_code(df["father_highest_education"])
    # Take each field from the first form that has a usable value.
    out = df.groupby("sample_id", sort=False).first()
    out["bghx_form"] = df.groupby("sample_id", sort=False)["bghx_form"].agg(",".join)
    log(f"    income known for {out['annual_household_income_ordinal'].notna().sum():,}")
    return out.reset_index()


def load_adi(cfg):
    df = pd.read_csv(release_file(cfg, "area_deprivation_index"), dtype=str, low_memory=False)
    df = df.rename(columns={"subject_sp_id": "sample_id"})
    df = df.sort_values("evaluation_date").drop_duplicates("sample_id", keep="last")
    adi_cols = [c for c in df.columns if c.startswith("X20")]
    for c in adi_cols:
        df[c] = pd.to_numeric(df[c].where(~df[c].isin(ADI_SPECIAL)), errors="coerce")
    df = df.rename(columns={c: "adi_" + c[1:].replace("_adi_", "_").replace("_rank_", "_")
                                .replace("percentile", "pct") for c in adi_cols})
    df = df.rename(columns={"evaluation_date": "adi_evaluation_date"})
    log(f"  ADI: {len(df):,} individuals")
    return df[["sample_id", "adi_evaluation_date"] + [c for c in df.columns if c.startswith("adi_20")]]


def load_air_pollution(cfg):
    df = pd.read_csv(cfg["inputs"]["spark"]["air_pollution"], dtype=str,
                     usecols=["subject_id"] + AQ_COLS, low_memory=False)
    df = numeric(df.rename(columns={"subject_id": "sample_id"}), AQ_COLS).drop_duplicates("sample_id")
    log(f"  air pollution: {len(df):,} individuals")
    return df


def load_asr(cfg):
    df = pd.read_csv(release_file(cfg, "asr"), dtype=str, low_memory=False,
                     usecols=["subject_sp_id", "age_at_eval_years", "highest_ed"] + ASR_ITEMS + ASR_SCALES)
    df = numeric(df, ["age_at_eval_years"] + ASR_ITEMS + ASR_SCALES)
    df["tobacco_user"] = (df["q124_used_tobacco"] > 0).astype("Int64").where(df["q124_used_tobacco"].notna())
    alc = (df["q090_drinks_too_much"].fillna(0) > 0) | (df["alcohol_raw_score"].fillna(0) > 0)
    df["alcohol_user"] = alc.astype("Int64").where(
        df["q090_drinks_too_much"].notna() | df["alcohol_raw_score"].notna())
    df = df.rename(columns={"subject_sp_id": "sample_id"}).drop_duplicates("sample_id")
    df = df.rename(columns={c: f"asr_{c}" for c in df.columns if c != "sample_id"})
    log(f"  ASR: {len(df):,} adults")
    return df


def id_list(path):
    df = pd.read_csv(path, dtype=str)
    return set(df.iloc[:, 0].str.strip())


def main():
    cfg = load_config()
    log("SPARK phenotypes")
    df = load_registration(cfg)
    asr = load_asr(cfg)
    for part in (load_core_descriptive(cfg), load_background(cfg), load_adi(cfg),
                 load_air_pollution(cfg), asr):
        df = left_merge(df, part)
    asr_keyed = asr.set_index("sample_id")
    for side, col in (("mother", "biomother_sp_id"), ("father", "biofather_sp_id")):
        for c in ASR_TO_PARENTS:
            df[f"{side}_asr_{c}"] = df[col].map(asr_keyed[f"asr_{c}"])

    df = left_merge(df, load_sfari_smoking(cfg["inputs"]["smoking"], "SPARK"))

    # Withdrawn individuals are usually absent from the release files, so add them as
    # ID-only rows; step 08 uses the flag to drop them from the assembled tables.
    withdrawn = id_list(supplement_file(cfg, "withdrawn_individuals"))
    paused = id_list(supplement_file(cfg, "paused_individuals"))
    missing = sorted((withdrawn | paused) - set(df["sample_id"]))
    df = pd.concat([df, pd.DataFrame({"sample_id": missing})], ignore_index=True)
    df["spark_withdrawn"] = df["sample_id"].isin(withdrawn)
    df["spark_paused"] = df["sample_id"].isin(paused)
    log(f"  withdrawn: {len(withdrawn):,}; paused: {len(paused):,}; "
        f"{len(missing):,} of these are not in the registration file")

    write_table(df, intermediate_path("phenotypes", "spark_phenotypes.parquet"))


if __name__ == "__main__":
    main()
