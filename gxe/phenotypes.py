"""Phenotype codings shared by the SPARK and SSC loaders."""
import pandas as pd

from gxe.io import log

# Harmonised parental education, 1 (least) .. 10 (most). See docs/education_coding.md.
EDU_CODE = {
    # SPARK
    "did_not_attend_high_school": 1,
    "some_high_school": 3,
    "ged_diploma": 4,
    "high_school_graduate": 5,
    "trade_school": 6,
    "some_college": 7,
    "associate_degree": 8,
    "baccalaureate_degree": 9,
    "graduate_or_professional_degree": 10,
    # SSC
    "less-ninth": 1,
    "up-ninth": 2,
    "some-hs": 3,
    "ged": 4,
    "high-school": 5,
    "some-college": 7,
    "associate": 8,
    "baccalaureate": 9,
    "graduate": 10,
}

# SPARK annual_household_income brackets -> ordinal 1..9. Anything else
# (e.g. "na_survey_logic", a survey-skip marker) is missing.
INCOME_BRACKETS = {
    "less_than_20000": 1, "21000_35000": 2, "36000_50000": 3,
    "51000_65000": 4, "66000_80000": 5, "81000_100000": 6,
    "101000_130000": 7, "131000_160000": 8, "over_161000": 9,
}

AQ_COLS = ["enroll_pm25_mean", "enroll_o3_mean", "enroll_no2_mean",
           "prenatal_pm25_mean", "prenatal_o3_mean", "prenatal_no2_mean"]

SMOKING_COLS = {
    "Father_Ever_Smoke": "father_ever_smoke",
    "Mother_Ever_Smoke": "mother_ever_smoke",
    "Father_Pack_Years": "father_pack_years",
    "Mother_Pack_Years": "mother_pack_years",
    "Parental_Smoking_Status": "parental_smoking_status",
}


def edu_code(s):
    return s.astype("string").str.strip().map(EDU_CODE).astype("Int64")


def income_ordinal(s):
    return s.astype("string").str.strip().str.lower().map(INCOME_BRACKETS).astype("Int64")


def numeric(df, cols):
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def load_sfari_smoking(path, cohort_label):
    """SFARI parental-smoking file rows for one study ("SPARK" or "SSC"), keyed by sample_id."""
    df = pd.read_csv(path, dtype=str, low_memory=False,
                     usecols=["Cohort", "Sample_ID", *SMOKING_COLS])
    df = df[df["Cohort"].str.upper() == cohort_label.upper()]
    df = df.rename(columns={"Sample_ID": "sample_id", **SMOKING_COLS}).drop(columns="Cohort")
    df = numeric(df, ["father_ever_smoke", "mother_ever_smoke",
                      "father_pack_years", "mother_pack_years"])
    df = df.drop_duplicates("sample_id")
    log(f"  SFARI smoking ({cohort_label}): {len(df):,} children")
    return df


def parent_age(child_birth, parent_birth, min_age):
    """Parent age at the child's birth; values below `min_age` are set to missing."""
    age = child_birth - parent_birth
    return age.where(age >= min_age).round(2)


def birth_order(group, birth):
    """1-based birth order within `group` (ties, e.g. twins, share the lower rank) and the
    number of children in the group with a known birth date. Missing where either is missing."""
    known = group.notna() & birth.notna()
    order = pd.Series(pd.NA, index=group.index, dtype="Float64")
    n = pd.Series(pd.NA, index=group.index, dtype="Float64")
    g = group[known]
    order[known] = birth[known].groupby(g).rank(method="min")
    n[known] = g.map(g.value_counts())
    return order.astype("Int64"), n.astype("Int64")
