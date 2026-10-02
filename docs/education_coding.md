# Parental education coding

SPARK and SSC record parental education with different labels. `gxe/phenotypes.py`
maps both to one ordinal code, 1 (least) to 10 (most), in `mother_edu_code` /
`father_edu_code`. The raw labels stay in `mother_highest_education` /
`father_highest_education`.

| Code | SPARK label | SSC label |
|---|---|---|
| 1 | did_not_attend_high_school | less-ninth |
| 2 | — | up-ninth |
| 3 | some_high_school | some-hs |
| 4 | ged_diploma | ged |
| 5 | high_school_graduate | high-school |
| 6 | trade_school | — |
| 7 | some_college | some-college |
| 8 | associate_degree | associate |
| 9 | baccalaureate_degree | baccalaureate |
| 10 | graduate_or_professional_degree | graduate |

Blank cells, the SSC values `true` / `false`, and any other label are missing.

Sources:
- SPARK: `background_history_child`, then `background_history_sibling`, then `background_history_adult`.
- SSC: `Proband Data/ssc_background_hx.csv`. This is one row per family, applied to every child in the family.

Judgement calls:
- `trade_school` (SPARK only) sits at 6, between high school (5) and some college (7).
- `up-ninth` (SSC only) sits at 2, between less than 9th grade (1) and some high school (3).
